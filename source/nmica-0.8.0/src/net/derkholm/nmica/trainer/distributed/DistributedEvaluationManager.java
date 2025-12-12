/*
 * NestedMICA Motif Inference Toolkit
 *
 * Copyright (c) 2004-2007: Genome Research Ltd.
 *
 * NestedMICA is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * or see the on-line version at http://www.gnu.org/copyleft/lgpl.txt
 *
 */

package net.derkholm.nmica.trainer.distributed;

import java.io.IOException;
import java.net.InetSocketAddress;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import net.derkholm.nmica.maths.DoubleProcedure;
import net.derkholm.nmica.matrix.ObjectMatrix2D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix2D;
import net.derkholm.nmica.model.ContributionGroup;
import net.derkholm.nmica.model.Datum;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.FacetteMap;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.trainer.EvaluationManager;
import net.derkholm.nmica.trainer.TrainableState;
import net.derkholm.nmica.trainer.TrainableStateContext;
import net.derkholm.nmica.trainer.distributed.messages.ContributionRequest;
import net.derkholm.nmica.trainer.distributed.messages.ContributionResponse;
import net.derkholm.nmica.trainer.distributed.messages.DatumRequest;
import net.derkholm.nmica.trainer.distributed.messages.DatumResponse;
import net.derkholm.nmica.trainer.distributed.messages.LikelihoodRequest;
import net.derkholm.nmica.trainer.distributed.messages.LikelihoodResponse;
import net.derkholm.nmica.trainer.distributed.messages.NotReady;
import net.derkholm.nmica.trainer.distributed.messages.Ready;
import net.derkholm.nmica.trainer.distributed.messages.Shutdown;
import net.derkholm.nmica.trainer.distributed.messages.TrainerConfigRequest;
import net.derkholm.nmica.trainer.distributed.messages.TrainerConfigResponse;
import net.derkholm.nmica.utils.WorkQueue;
import net.derkholm.nmica.utils.mq.MessageHandler;
import net.derkholm.nmica.utils.mq.MessageQueue;
import net.derkholm.nmica.utils.mq.Packable;
import net.derkholm.nmica.utils.mq.QueueDeadException;
import net.derkholm.nmica.utils.mq.MessageQueue.Message;

/**
 * EvaluationManager which farms out jobs to remote nodes
 * 
 * @author thomas
 * @version branch-jdk15
 */

public class DistributedEvaluationManager implements EvaluationManager {
    private List<WorkUnit> workList = new ArrayList<WorkUnit>();
    private TrainableStateContext trainer;
    private TrainableState currentState;
    private short currentSid = 0;
    private MessageQueue<Packable> messageQ;
    private volatile boolean inShutdown = false;
    private Set<WorkerRecord> readySet = new HashSet<WorkerRecord>();
    private boolean seedContributions = true;
    private boolean debug = false;
    private Thread ticker;
    private int hoods = 0;
    private double crabSignifier = 4.0;
    private double crabRate = 0.001;
    private WorkQueue workQueue;
    private int crabSteps = 0;
    private int minOffload = 0;
    private boolean disableBalancer = false;
    
    private int localThreads = 0;
    private WorkerRecord localWorker = new WorkerRecord(null);
    private ObjectMatrix2D hoodCalcs;

    private volatile FinishLine finishLine = null;

    private Thread shutdownHook;

    /**
     * Damping parameter controlling the rate at which node performance weights
     * are adjusted.
     */
    public void setCrabRate(double d) {
        this.crabRate = d;
    }
    public void setCrabSignifier(double d) {
    	this.crabSignifier = d;
    }
    public void setMinOffload(int i) {
    	this.minOffload = i;
    }
    public void setDisableBalancer(boolean b) {
    	this.disableBalancer = b;
    }
    
    /**
     * Number of local threads to use for evaluating 'hoods, in addition to remote workers.
     * 
     * @param i
     */
    public void setLocalThreads(int i) {
    	this.localThreads = i;
    	this.localWorker.advertisedPower = i;
    	if (localThreads == 0) {
    		readySet.remove(localWorker);
    	} else {
    		readySet.add(localWorker);
    	}
    }
    
    private void makeHoodCalcs() {
        FacetteMap fm = trainer.getFacetteMap();
        Facette[] facettes = fm.getFacettes();
        Datum[] data = trainer.getDataSet();
        
        hoodCalcs = new SimpleObjectMatrix2D(facettes.length, trainer.getDataSet().length);
        for (int d = 0; d < data.length; ++d) {
            Object[] fd = data[d].getFacettedData();
            for (int f = 0; f < facettes.length; ++f) {
                if (fd[f] != null) {
                    hoodCalcs.set(f, d, facettes[f].getLikelihoodCalculator(fd[f]));
                }
            }
        }
    }
    
    private class BandwidthThread extends Thread {
        private long oldTime = -1;
        private int oldTxPackets = 0;
        private int oldRxPackets = 0;
        private int oldTxBytes = 0;
        private int oldRxBytes = 0;
        private int oldHoods = 0;

        private final NumberFormat FORMAT = new DecimalFormat("######0.0");

        public void run() {
            while (!inShutdown) {
                long time = System.currentTimeMillis();
                int newTxPackets = messageQ.getTxPackets();
                int newTxBytes = messageQ.getTxBytes();
                int newRxPackets = messageQ.getRxPackets();
                int newRxBytes = messageQ.getRxBytes();
                int newHoods = hoods;
                if (oldTime > 0) {
                    int deltaTime = (int) (time - oldTime);

                    StringBuffer sb = new StringBuffer();
                    sb.append("Out ");
                    sb.append(FORMAT.format((1.0 * (newTxBytes - oldTxBytes))
                            / deltaTime));
                    sb.append("kb/s ");
                    sb.append(FORMAT
                            .format((1000.0 * (newTxPackets - oldTxPackets))
                                    / deltaTime));
                    sb.append("packs/s       In ");
                    sb.append(FORMAT.format((1.0 * (newRxBytes - oldRxBytes))
                            / deltaTime));
                    sb.append("kb/s ");
                    sb.append(FORMAT
                            .format((1000.0 * (newRxPackets - oldRxPackets))
                                    / deltaTime));
                    sb.append("packs/s       Throughput ");
                    sb
                            .append(FORMAT
                                    .format((1000.0 * (newHoods - oldHoods) / deltaTime)));
                    sb.append("lc/s");
                    System.err.println(sb.toString());
                    
                    double pps = (1000.0 * (newTxPackets - oldTxPackets)) / deltaTime;
                    
                    for (WorkerRecord w : readySet) {
                    	System.out.println(w.toString() + ": " + w.weight);
                    }
                    
                    // messageQ.setSlug(pps < 20);
                }
                oldTxPackets = newTxPackets;
                oldTxBytes = newTxBytes;
                oldRxPackets = newRxPackets;
                oldRxBytes = newRxBytes;
                oldHoods = newHoods;
                oldTime = time;

                try {
                    Thread.sleep(1000L);
                } catch (InterruptedException ex) {
                }
            }
        }

    }

    public void shutdown() {
        Runtime.getRuntime().removeShutdownHook(shutdownHook);
        doShutdown();
    }

    private void doShutdown() {
        for (WorkerRecord w : readySet) {
            try {
            	if (w.node != null) {
            		messageQ.sendMessage(w.node, new Shutdown());
            	}
            } catch (QueueDeadException ex) {
            }
        }
        messageQ.shutdown();
        this.inShutdown = true;
    }

    public DistributedEvaluationManager(int port, boolean bandwidthTicker, boolean crabMonitor) throws IOException {
        messageQ = new MessageQueue<Packable>(port);
        messageQ.setCodec(Protocol.CODEC);
        messageQ.setFlushThreshold(40);
        messageQ.setSynchronousHandler(new TMonHandler());
        messageQ.start();

        if (bandwidthTicker) {
            ticker = new BandwidthThread();
            ticker.start();
        }
        
        if (crabMonitor) {
        	Thread cm = new Thread() {
        		public void run() {
        			while (true) {
	        			System.err.println("---------- CRAB weights -----------");
	        			System.err.printf("Total steps     : %d%n", crabSteps);
	        			for (WorkerRecord r : readySet) {
	        				System.err.printf("%-20s: %g (%d)%n", r.node == null ? "<local>" : r.node.getEndpoint().toString(), r.weight, r.wuCount);
	        				r.wuCount = 0;
	        			}
	        			System.err.println("-----------------------------------");
	        			try {
	        				Thread.sleep(30000L);
	        			} catch (InterruptedException ex) {}
        			}
        		}
        	};
        	cm.setDaemon(true);
        	cm.start();
        }

        shutdownHook = new Thread() {
            public void run() {
                doShutdown();
            }
        };
        Runtime.getRuntime().addShutdownHook(shutdownHook);
    }

    public InetSocketAddress getDatagramEndpoint() {
        return messageQ.getEndpoint();
    }

    private WorkQueue getWorkQueue() {
        if (workQueue == null) {
            workQueue = WorkQueue.create(localThreads);
        }
        return workQueue;
    }
    
    public void startLikelihoodCalculations(TrainableState state) {
        if (currentState != null) {
            throw new IllegalStateException("Can't start likelihood calculations while queues are busy");
        }
        if (trainer == null) {
            trainer = state.getContext();
            makeHoodCalcs();
        } else if (trainer != state.getContext()) {
            throw new RuntimeException("Already bound to a trainer");
        }

        currentState = state;
        ++currentSid;
        if (currentSid > 1 << 14) {
            currentSid = 0; // sid wrap-around for the sake of sanity.
        }
        workList.clear();
    }

    
    
    public void enqueueLikelihoodCalculation(TrainableState state, int d, int f, DoubleProcedure writeback) {
        if (state != currentState) {
            throw new IllegalStateException(
                    "Requested state is not open for likelihood evaluations");
        }

        // synchronized (workList) {
            // this needs to be monitored, to ensure unique WID generation.
            int wid = workList.size();
            WorkUnit wu = new WorkUnit(
            		wid, 
            		currentSid, 
            		f, 
            		trainer.facetteIndexToContributionIndex(f), 
            		d, 
            		/* new SimpleMatrix1D(state.getMixture(d)), // defensive copy to de-view this */
            		state.getMixtureReadOnly(d),                        // No longer copy this.  Should be okay since the net protocol doesn't use serialization any more
                    writeback
            );

            workList.add(wu);
        // }
    }

    /**
     * Helper class which serves as a latch for monitoring the completion of a
     * job-set, and also logs completion ranks.
     */

    private static final class FinishLine {
        private CountDownLatch readyLatch = new CountDownLatch(1);
        private CountDownLatch finishedLatch = new CountDownLatch(1);
        private int[] counts;
        private int total = 0;
        private int completed = 0;
        private int winner = -1;
        private int loser = -1;

        FinishLine() {
        }

        public void setCounts(int[] counts) {
            this.counts = new int[counts.length];
            for (int i = 0; i < counts.length; ++i) {
            	this.counts[i] = counts[i];
                total += counts[i];
            }
            readyLatch.countDown();
            if (total == 0) {
                finishedLatch.countDown();
            }
        }

        public synchronized void commitAndDecrement(WorkUnit wu, double likelihood)
                throws InterruptedException 
        {
            readyLatch.await();
            if (wu.writeback != null) {
	            wu.writeback.run(likelihood);
	            wu.writeback = null;
	            int i = wu.assignedWorkerID;
	            --counts[i];
	            if (counts[i] == 0 && winner < 0) {
	                winner = i;
	            }
	            ++completed;
	            if (completed == total) {
	                loser = i;
	                finishedLatch.countDown();
	            }
            }
        }

        public int getTotal() {
            return total;
        }

        public int getCompleted() {
            return completed;
        }

        public int getWinner() {
            return winner;
        }

        public int getLoser() {
            return loser;
        }

        public void await() 
            throws InterruptedException 
        {
            finishedLatch.await();
        }

        public boolean await(long timeout, TimeUnit unit)
            throws InterruptedException 
        {
            return finishedLatch.await(timeout, unit);
        }
    }

    private long lastRemoteWarning = 0;
    
    public void endLikelihoodCalculations(TrainableState state) {
        if (state != currentState) {
            throw new IllegalStateException("Requested state is not open for likelihood evaluations");
        }

        int[] wuCounts;
        
        if (workList.size() > 0) {
        
	        RELIABLE_HOOD_LOOP: while (true) {
	            WorkerRecord[] workers = null;
	            do {
	                if (workers != null) {
	                    try {
	                        Thread.sleep(1000L);
	                    } catch (Exception ex) {
	                    }
	                }
	                workers = readySet.toArray(new WorkerRecord[readySet.size()]);
	                if (workers.length < 2 && System.currentTimeMillis() - lastRemoteWarning >= 30000L) {
	                	boolean hasRemote = false;
	                	for (WorkerRecord w : workers) {
	                		hasRemote = w.node != null;
	                	}
	                	if (!hasRemote) {
	                		System.err.printf("Warning: running in distributed mode but no dlepnodes are currently connected%n");
	                		lastRemoteWarning = System.currentTimeMillis();
	                	}
	                }
	            } while (workers.length == 0);
	
	            if (seedContributions) {
	                ContributionGroup[] cgs = trainer.getFacetteMap().getContributionGroups();
	                for (int cg = 0; cg < cgs.length; ++cg) {
	                    for (int c = 0; c < trainer.getComponents(); ++c) {
	                        Object seedItem = currentState.getContribution(cg, c).getItem();
	                        ContributionResponse resp = new ContributionResponse();
	                        resp.component = c;
	                        resp.contributionGroup = cg;
	                        resp.sid = currentSid;
	                        resp.contribution = seedItem;
	                        for (int n = 0; n < workers.length; ++n) {
	                            try {
	                            	if (workers[n].node != null) {
	                            		messageQ.sendMessage(workers[n].node, resp);
	                            	}
	                            } catch (QueueDeadException ex) {
	
	                            }
	                        }
	                    }
	                }
	            }
	            // messageQ.flush();
	
	            boolean countsForCrabUpdate = !disableBalancer;
	            finishLine = new FinishLine();
	            WorkQueue localQueue = null;
	            // once a new finishline is installed, wu's can't be checked in
	            // until we
	            // release the latch by setting a new count array.
	            synchronized (workList) { // is this necessary
	            	if (disableBalancer) {
	            		double maxPower = 0;
	            		for (WorkerRecord w : workers) {
	            			w.weight = w.advertisedPower;
	            			maxPower = Math.max(maxPower, w.weight);
	            		}
	            		for (WorkerRecord w : workers) {
	            			w.weight /= maxPower;
	            		}
	            	}
	                wuCounts = new int[workers.length];
	                double[] crabWeights = new double[workers.length];
	                double crabTotal = 0.0;
                    for (int w = 0; w < workers.length; ++w) {
                        crabWeights[w] = workers[w].weight * 1.3142;  // random-ish non-integer to help with error diffusion.
                        crabTotal += crabWeights[w];
                    }
	                int totalSent = 0;
	                double crabCounter = Math.random() * crabTotal; // easy way to
	                                                                // start the
	                                                                // round-robin
	                                                                // on a sensible
	                                                                // random place
	                int wls = workList.size();
	                for (int wui = 0; wui < wls; ++wui) {
	                	final WorkUnit wu = workList.get(wui);
	                    if (wu.writeback != null) {
	                    	int targetWorker = -1;
                            while (crabCounter >= crabTotal) {
                                crabCounter -= crabTotal;
                            }
                            double tt = 0.0;
                            for (int w = 0; w < workers.length; ++w) {
                                tt += crabWeights[w];
                                if (tt >= crabCounter) {
                                    targetWorker = w;
                                    break;
                                }
                            }
                            crabCounter += 1.0;
                            wu.assignedWorkerID = targetWorker;

                            WorkerRecord twr = workers[targetWorker];
                            if (twr.node == null) {
                            	if (localQueue == null) {
                            		localQueue = getWorkQueue();
                            	}
                            	localQueue.add(new Runnable() {
									public void run() {
						                LikelihoodCalculator lc = (LikelihoodCalculator) hoodCalcs.get(wu.facette, wu.datum);
						                double term = lc.likelihood(currentState.getContributions(trainer.facetteIndexToContributionIndex(wu.facette)), wu.weights);
						                try {
											finishLine.commitAndDecrement(wu, term);
										} catch (InterruptedException e) {
											// TODO Auto-generated catch block
											e.printStackTrace();
										}
									}
                            	});
                            } else {
		                        LikelihoodRequest req = new LikelihoodRequest();
		                        req.sid = currentSid;
		                        req.wid = wu.wid;
		                        req.contributionGroup = wu.contributionGroup;
		                        req.datum = wu.datum;
		                        req.facette = wu.facette;
		                        req.weights = wu.weights;
		                        try {
		                            messageQ.sendMessage(workers[targetWorker].node, req);
		                        } catch (QueueDeadException ex) {
		
		                        }
                            }
                            ++twr.wuCount;
                            ++wuCounts[targetWorker];
                            ++totalSent;
	                    }
	                }
	                // used to unCrab a cycle if any node had less that 4 wus assigned,
	                // but this meant that if a node got a pathologically low weight, the
	                // whole CRAB system could get stuck. Hopefully this will have a 
	                // similar effect.
	                if (totalSent < (crabSignifier * workers.length)) {
	                    countsForCrabUpdate = false;
	                }
	                for (int wuCount : wuCounts) {
	                	if (wuCount == 0) {
	                		countsForCrabUpdate = false;
	                	}
	                }
	                // this implicitly releases the comms thread if it's already
	                // waiting to cross the finish line.
	                finishLine.setCounts(wuCounts);
	            }
	
	            // long calcPhaseStart = System.currentTimeMillis();
	            messageQ.flush();
	            if (localQueue != null) {
	            	localQueue.flush();
	            }
	            messageQ.flush();
	            int completedBeforeAwait = 0;
	            int misTicks = 0;
	
	          EVAL_WATCHDOG_LOOP: 
	        	while (true) {
	                try {
	                    if (finishLine.await(500L, TimeUnit.MILLISECONDS)) {
	                        if (countsForCrabUpdate) {
	                        	/*
	                        	
	                            // Currently ignoring the loser, since it's hard (impossible?)
	                            // to distinguish beween slow nodes and one-off glitches.
	                            int winner = finishLine.getWinner();
	                            if (winner >= 0) {
	                            	boolean counts2 = true;
	                            	if (workers[winner].weight >= 1.0) {
	                            		int m = 0;
	                            		for (int wc : wuCounts) {
	                            			m = Math.max(m, wc);
	                            		}
	                            		if (m > wuCounts[winner]) {
	                            			counts2 = false;
	                            		}
	                            	}
	                            	if (counts2) {
		                                workers[winner].weight += crabRate;
		                                double maxWeight = 0.0;
		                                for (int w = 0; w < workers.length; ++w) {
		                                    maxWeight = Math.max(maxWeight, workers[w].weight);
		                                }
		                                for (int w = 0; w < workers.length; ++w) {
		                                    workers[w].weight /= maxWeight;
		                                }
		                                ++crabSteps;
	                            	}
	                            }
	                            
	                            */
	                        	
	                            int loser = finishLine.getLoser();
	                            if (loser >= 0) {
	                            	boolean counts2 = true;
	                            	if (counts2) {
		                                workers[loser].weight -= crabRate;
		                                double maxWeight = 0.0;
		                                for (int w = 0; w < workers.length; ++w) {
		                                    maxWeight = Math.max(maxWeight, workers[w].weight);
		                                }
		                                for (int w = 0; w < workers.length; ++w) {
		                                    workers[w].weight /= maxWeight;
		                                }
		                                ++crabSteps;
	                            	}
	                            }
	                        }
	                        break RELIABLE_HOOD_LOOP;
	                    } else {
	                        int completedAfterAwait = finishLine.getCompleted();
	                        if (completedAfterAwait == completedBeforeAwait) {
	                            // System.err.println("Tick has passed without any work checked in");
	                        	countsForCrabUpdate = false;
	                            ++misTicks;
	                            if (misTicks > 1) {
	                                break EVAL_WATCHDOG_LOOP;
	                            }
	                        }
	                        completedBeforeAwait = completedAfterAwait;
	                    }
	                } catch (Exception ex) {
	                    throw new RuntimeException(ex);
	                }
	            }
	            int total = finishLine.getTotal();
	            int remaining = finishLine.getTotal() - finishLine.getCompleted();
	            System.err.println("Some likelihood calculations (" + remaining + "/" + total + ") timed out, restarting");
	            for (WorkUnit wu : workList) {
	                if (wu.writeback != null) {
	                    // System.err.printf("Unit assigned to %s%n", workers[wu.assignedWorkerID].node);
	                	++workers[wu.assignedWorkerID].wuLost;
	                }
	            }
	
	            /*
	             * for (Iterator ri = readySet.iterator(); ri.hasNext(); ) {
	             * WorkerRecord r = (WorkerRecord) ri.next(); try {
	             * messageQ.sendMessage(r.node, new Flush()); } catch
	             * (QueueDeadException ex) { } }
	             */
	
	            long now = System.currentTimeMillis();
	            for (Iterator<WorkerRecord> ri = readySet.iterator(); ri.hasNext();) {
	                WorkerRecord r = ri.next();
	                if (r.node != null && now - r.lastPing > 3000) {
	                    System.err.printf("Lost contact with %s, removing from ready-set%n", r.node.getEndpoint().toString());
	                    ri.remove();
	                }
	            }
	        }
        }

        currentState = null;
    }

    private static class WorkerRecord {
        public final MessageQueue.Peer node;
        public long lastPing = 0;
        public double weight = 1.0;
        public double advertisedPower = 1.0;
        
        public int wuCount = 0;
        public int wuLost = 0;
        
        public WorkerRecord(MessageQueue.Peer node) {
            this.node = node;
        }
        
        public String toString() {
        	return node == null ? "local-threads" : node.getEndpoint().toString();
        }
    }

    private class TMonHandler implements MessageHandler<Packable> {
		public void handleMessage(MessageQueue<Packable> q, Message<Packable> msg) {
			try {
	            MessageQueue.Peer sender = (MessageQueue.Peer) msg.getSender(); // needed to keep ecj happy...
	            WorkerRecord worker = (WorkerRecord) sender.getUserData();
	            if (worker == null) {
	                System.err.printf("New node joined: %s%n", sender.getEndpoint().toString());
	                worker = new WorkerRecord(sender);
	                sender.setUserData(worker);
	            }
	
	            Packable body = msg.getBody();
	            if (debug) {
	                System.err.println("Got " + body.getClass().getName());
	            }
	            if (body instanceof TrainerConfigRequest) {
	                TrainerConfigResponse resp = new TrainerConfigResponse();
	                resp.components = trainer.getComponents();
	                resp.dataSetSize = trainer.getDataSet().length;
	                resp.facetteMap = trainer.getFacetteMap();
	
	                messageQ.sendMessage(sender, resp);
	                messageQ.flush();
	            } else if (body instanceof Ready) {
	                readySet.add(worker);
	                worker.lastPing = System.currentTimeMillis();
	                worker.advertisedPower = ((Ready) body).power;
	            } else if (body instanceof NotReady) {
	                // if (!readySet.contains(worker)) {
	                //     System.err.println("Strange, node isn't in ready set");
	                // }
	                readySet.remove(worker);
	            } else if (body instanceof DatumRequest) {
	                DatumRequest req = (DatumRequest) body;
	                DatumResponse resp = new DatumResponse();
	                resp.datumIndex = req.datumIndex;
	                resp.facette = req.facette;
	                resp.datum = trainer.getDataSet()[req.datumIndex].getFacettedData()[req.facette];
	                messageQ.sendMessage(sender, resp);
	                messageQ.flush();
	            } else if (body instanceof ContributionRequest) {
	            	try {
	                    ContributionRequest req = (ContributionRequest) body;
	                    ContributionResponse resp = new ContributionResponse();
	                    resp.sid = currentSid;
	                    resp.component = req.component;
	                    resp.contributionGroup = req.contributionGroup;
	                    resp.contribution = currentState.getContribution(req.contributionGroup, req.component).getItem();
	                    messageQ.sendMessage(sender, resp);
	            	} catch (Exception ex) {
	            		ex.printStackTrace();
	            	}
	            } else if (body instanceof LikelihoodResponse) {
	                LikelihoodResponse resp = (LikelihoodResponse) body;
	                if (resp.sid != currentSid) {
	                    // System.err.println("Yuck, got back an out-of-date reponse");
	                } else {
		                WorkUnit wu = workList.get(resp.wid);
		                if (wu.writeback != null) {
		                    // now checks in the result and diddles the latches
		                    // atomically. ooops.
		                    finishLine.commitAndDecrement(wu, resp.likelihood);
		                    ++hoods;
		                }
	                }
	            } else {
	                System.err.println("Unrecognized message class: " + body.getClass().getName());
	            }
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
    	
    }
}
