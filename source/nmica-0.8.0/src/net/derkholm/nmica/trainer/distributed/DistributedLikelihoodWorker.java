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

import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.ObjectMatrix2D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix2D;
import net.derkholm.nmica.model.*;
import net.derkholm.nmica.trainer.distributed.messages.*;
import net.derkholm.nmica.trainer.distributed.messages.Shutdown;
import net.derkholm.nmica.utils.mq.MessageHandler;
import net.derkholm.nmica.utils.mq.MessageQueue;
import net.derkholm.nmica.utils.mq.QueueDeadException;
import net.derkholm.nmica.utils.mq.Packable;
import net.derkholm.nmica.utils.mq.MessageQueue.Message;
import net.derkholm.nmica.utils.tracker.SimpleTracker;
import net.derkholm.nmica.utils.tracker.Task;
import net.derkholm.nmica.utils.tracker.Tracker;
import org.biojava.bio.symbol.SymbolList;

import java.io.IOException;
import java.lang.reflect.Method;
import java.net.InetSocketAddress;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * Object which connects to a remote trainer and performs likelihood calculations.
 * 
 * @author thomas
 */
public class DistributedLikelihoodWorker implements Runnable, Worker {
    private InetSocketAddress endpoint;
    private int numThreads;
    private int lruSize = 1000;
    
    private MessageQueue<Packable> messageQ;
    private MessageQueue.Peer master;
    private EvaluationQueue evalQ;
    
    private boolean running = true;
    
    // long-term state
    private int components;
    private Facette[] facettes;
    private ContributionGroup[] contributionGroups;
    private ObjectMatrix2D hoodCalcs;
    
    // transient state
    private int currentSid = -1;
    private ObjectMatrix2D contributions = null;
    private Map<Object,ContributionItem> contributionCache = new LinkedHashMap<Object,ContributionItem>(1000, 0.75F, true) {
        protected boolean removeEldestEntry(Map.Entry<Object,ContributionItem> e) {
            return size() > lruSize;
        }
    };
    
    // Monitoring
    private Task currentDatumFetchTask;
    
    // Throughput testing
    private boolean throughputMonitor = false;
    private boolean throughputWatchdog = false;
    private int tasksRun = 0;
    private int bases = 0;
    private int baseCards = 0;
    
    // Debugging
    private int sleepy = 0;
    
    public void setThroughputMonitor(boolean b) {
        this.throughputMonitor = b;
    }
    
    public void setThroughputWatchdog(boolean b) {
        this.throughputWatchdog = b;
    }
    
    public void setSleepy(int i) {
        this.sleepy = i;
    }
    
    public void setLruSize(int i) {
        this.lruSize = i;
    }
    
    private class ThroughputMonitor extends Thread {
        private int oldTasksRun;
        private int oldBases;
        private long oldTime = -1;
        
        public void run() {
            long lastWorkTime = -1;
            while (running) {
                long newTime = System.currentTimeMillis();
                int newTasksRun = tasksRun;
                int newBases = bases;
                if (throughputMonitor && oldTime > 0) {
                    int deltaTime = (int) (newTime - oldTime);
                    System.err.print("" + ((1000 * (newTasksRun - oldTasksRun)) / deltaTime) + "\t" + ((newBases - oldBases) / deltaTime) + "\t" + evalQ.size());
                    for (EvaluationQueue.EQThread thread : evalQ.threads) {
                        System.err.print("\t" + thread.state.toString());
                    }
                    System.err.println();
                }
                if (oldTasksRun != newTasksRun) {
                    lastWorkTime = newTime;
                }
                oldTasksRun = newTasksRun;
                oldBases = newBases;
                oldTime = newTime;
                
                if (throughputWatchdog && lastWorkTime > 0 && (newTime - lastWorkTime) > 20000L) {
                    System.err.println("Throughtput watchdog triggered, killing old workers!");
                    evalQ.repopulate();
                }
                
                try {
                    Thread.sleep(1000L);
                } catch (InterruptedException ex) {}
            }
        }
    }
    
    private LikelihoodCalculator getLikelihoodCalculator(int f, int d) {
        return (LikelihoodCalculator) hoodCalcs.get(f, d);
    }
    
    private void initLikelihoodCalculators() 
    		throws Exception
    {
        for (int f = 0; f < hoodCalcs.rows(); ++f) {
            for (int d = 0; d < hoodCalcs.columns(); ++d) {
                LikelihoodCalculator calc = facettes[f].getLikelihoodCalculator(getDatum(d, f));
                hoodCalcs.set(f, d, calc);
            }
        }
            
        try {
            Method shortCircuitCall = org.biojava.utils.ChangeSupport.class.getMethod(
                    "setGlobalChangeBypass", 
                    new Class[] {Boolean.TYPE}
            );
            shortCircuitCall.invoke(null, new Object[] {Boolean.TRUE});
        } catch (NoSuchMethodException ex) {
            System.err.println("Short-circuiting isn't available");
        }
    }
    
    private Object getDatum(int d, int f) 
    	    throws Exception
    {
    	long start = System.currentTimeMillis();
        // System.err.println("Getting datum " + d);
        
        Tracker tracker = new SimpleTracker();
        currentDatumFetchTask = tracker.newTask();
        
        DatumRequest request = new DatumRequest();
        request.datumIndex = d;
        request.facette = f;
        currentDatumFetchTask.setData(request);
        messageQ.sendMessage(master, request);
        messageQ.flush();
        
        tracker.waitForTasks(2000L);
        
        DatumResponse resp = (DatumResponse) currentDatumFetchTask.getResult();
        if (resp == null) {
            System.err.println("Didn't get datum");
            return getDatum(d, f);
        } else {
        	// System.err.printf("Got datum (%dms)%n", System.currentTimeMillis() - start);
        }
        return resp.datum;
    }
    
    DistributedLikelihoodWorker(InetSocketAddress endpoint, int numThreads) {
        this.endpoint = endpoint;
        this.numThreads = numThreads;
    }
    
    public static DistributedLikelihoodWorker connect(String regServer, int regPort, int numThreads)
    		throws Exception
    {
        InetSocketAddress endpoint = new InetSocketAddress(regServer, regPort);
        return new DistributedLikelihoodWorker(endpoint, numThreads);
    }
    
    public void start() 
    {
        try {
            messageQ = new MessageQueue<Packable>();
            messageQ.setCodec(Protocol.CODEC);
            messageQ.setFlushThreshold(200);
            messageQ.start();
            
            master = messageQ.getPeer(endpoint);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        
        evalQ = new EvaluationQueue(numThreads);
        evalQ.start();
        
        if (throughputMonitor || throughputWatchdog) {
            new ThroughputMonitor().start();
        }
        
        this.run();
    }
    
    public void stop() {
        running = false;
    }
    
    private boolean updateSid(int sid) 
    		throws Exception
    {
        if (sid != currentSid) {
            if (sid < currentSid && currentSid - sid < 1000) {
                System.err.println("Hit a request for an out-of-date SID (currentSid=" + currentSid + " requestedSid=" + sid + "), bailing out...");
                return false;
            }
            
            contributions = new SimpleObjectMatrix2D(contributionGroups.length, components);
            currentSid = sid;
        }
        return true;
    }
    
    public void run() {
        try {
            messageQ.sendMessage(master, new TrainerConfigRequest());
            messageQ.flush();
            Object o = messageQ.next().getBody();
            if (o instanceof TrainerConfigResponse) {
                TrainerConfigResponse trc = (TrainerConfigResponse) o;
                components = trc.components;
                facettes = trc.facetteMap.getFacettes();
                contributionGroups = trc.facetteMap.getContributionGroups();
                hoodCalcs = new SimpleObjectMatrix2D(facettes.length, trc.dataSetSize);
                
                System.err.println("We're configured!");
            } else {
                throw new Exception("Unexpected response: " + o);
            }
            
            Thread readyPinger = new Thread("Readiness pinger") {
                public void run() {
                    try {
                        initLikelihoodCalculators();
                    } catch (Exception ex) {
                        throw new RuntimeException(ex);
                    }
                    while (true) {
                        try {
                        	Ready ready = new Ready();
                        	ready.power = (float) numThreads;
                            messageQ.sendMessage(master, ready);
                            Thread.sleep(1000L);
                        } catch (InterruptedException ex) {
                            
                        } catch (QueueDeadException ex) {
                            return;
                        }
                    }
                }
            };
            readyPinger.setDaemon(true);
            readyPinger.start();
            
            messageQ.setSynchronousHandler(new MessageHandler<Packable>() {
				public void handleMessage(MessageQueue<Packable> q, Message<Packable> msg) {
	               	Packable body = msg.getBody();
	                
	               	try {
	                if (body instanceof LikelihoodRequest) {
	                    LikelihoodRequest req = (LikelihoodRequest) body;
	                    if (updateSid(req.sid)) {
	                        evalQ.enqueueWork(msg);
	                    }
	                } else if (body instanceof DatumResponse) {
	                    currentDatumFetchTask.completed(body);
	                } else if (body instanceof ContributionResponse) {
	                    ContributionResponse resp = (ContributionResponse) body;
	                    if (updateSid(resp.sid)) {
		                    Object obj = contributions.get(resp.contributionGroup, resp.component);
		                    Task task = null;
		                    if (obj instanceof Task) {
		                        task = (Task) obj;
		                    }
	                        // we're the only thread that touches the LRU contrib cache.
	                        ContributionItem ci = contributionCache.get(resp.contribution);
	                        if (ci == null) {
	                            ci = new SimpleContributionItem(resp.contribution);
	                            contributionCache.put(resp.contribution, ci);
	                        } 
		                    contributions.set(resp.contributionGroup, resp.component, ci);
		                    if (task != null) {
		                        task.completed();
		                    }
	                    }
	                } else if (body instanceof Flush) {
	                    evalQ.flush();
	                } else if (body instanceof Shutdown) {
	                    running = false;
	                } else {
	                    System.err.println("Unexpected message type " + body.getClass().getName());
	                }	
	               	} catch (Exception ex) {
	               		throw new RuntimeException("Error communicating with trainer", ex);
	               	}
				}
            	
            });
            
        } catch (Exception ex) {
        	throw new RuntimeException("Error communicating with trainer", ex);
        }
    }
    
    private static enum EQThreadState {
        WAIT, PREP, CALC, RESP;
    }
    
    private class EvaluationQueue {
        private BlockingQueue<MessageQueue.Message<Packable>> work
                = new LinkedBlockingQueue<MessageQueue.Message<Packable>>();
        private EQThread[] threads;
	private boolean run = true;
        
        public EvaluationQueue(int nthreads) {
            this.threads = new EQThread[nthreads];
        }
        
        public void start() {
            for (int t = 0; t < threads.length; ++t) {
                threads[t] = new EQThread();
                threads[t].start();
            }
        }
        
        public void repopulate() {
            work = new LinkedBlockingQueue<MessageQueue.Message<Packable>>();
            for (int t = 0; t < threads.length; ++t) {
                threads[t].localRun = false;
                threads[t].interrupt();
                threads[t] = new EQThread();
                threads[t].start();
            }
        }
        
        public void flush() {
            synchronized (work) {
                work.clear();
            }
        }
        
        public int size() {
            return work.size();
        }
        
        public void enqueueWork(MessageQueue.Message<Packable> wu) {
	    while (true) {
		try {
		    work.put(wu);
		    return;
		} catch (InterruptedException ex) {}
	    }
        }
        
        public void shutdown() {
            run = false;
        }
        

        
        private class EQThread extends Thread {
            public boolean localRun = true;
            public boolean running = true;
            public EQThreadState state = EQThreadState.WAIT;
            
	        public void run() {
                int suspicious = 0;
	          HOOD_LOOP:
	            while (run && localRun) {
                    try {
                        state = EQThreadState.WAIT;
                        MessageQueue.Message<Packable> mreq = null;
                        try {
                            mreq = work.poll(200L, TimeUnit.MICROSECONDS);
                        } catch (InterruptedException ex) {
                            System.err.println("EQThread was interrupted");
                        }
                        if (mreq != null) {
                            state = EQThreadState.PREP;
                            LikelihoodRequest req = (LikelihoodRequest) mreq.getBody();

                            LikelihoodCalculator calc = getLikelihoodCalculator(req.facette, req.datum);
                            if (throughputMonitor) {
                                int reqBases = ((SymbolList) calc.getData()).length();
                                bases += reqBases;
                            }
                            state = EQThreadState.CALC;
                            double hood = calc.likelihood(
                                    new ContributionGopher(
                                            req.contributionGroup),
                                    req.weights);
                            state = EQThreadState.RESP;

                            LikelihoodResponse resp = new LikelihoodResponse();
                            resp.sid = req.sid;
                            resp.wid = req.wid;
                            resp.likelihood = hood;
                            if (sleepy > 0 && Math.random() < (1.0 / sleepy)) {
                                Thread.sleep(1L);
                            }
                            messageQ.sendMessage(master, resp);
                            if (work.isEmpty()) {
                                messageQ.flush();
                            }
                            ++tasksRun;
                            suspicious = 0;
                        } else if (work.size() > 0) {
                            ++suspicious;
                            System.err.printf("Poll timed out with queue non-empty, suspicious (%d)%n", suspicious);
                        }
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
	        
	            // System.err.println("EQThread is shutting down");
	        }
        }
    }
    
    private class ContributionGopher implements ObjectMatrix1D {
        private final int cg;
        
        public ContributionGopher(int cg) {
            this.cg = cg;
        }

        /* (non-Javadoc)
         * @see net.derkholm.nmica.matrix.ObjectMatrix1D#size()
         */
        public int size() {
            return components;
        }

        /* (non-Javadoc)
         * @see net.derkholm.nmica.matrix.ObjectMatrix1D#get(int)
         */
        public Object get(int pos) {
            Object o = contributions.get(cg, pos);
            if (o == null) {
                // System.err.println("Fetching contribs");
	            try {
	                Tracker tracker = new SimpleTracker();
	                
	                int cnt = 0;
	                for (int p = 0; p < components; ++p) {
	                    if (contributions.get(cg, p) == null) {
			                Task t = tracker.newTask();
			                
			                ContributionRequest req = new ContributionRequest();
			                req.sid = (short) currentSid;
			                req.component = p;
			                req.contributionGroup = cg;
			                
			                t.setData(req);
			                contributions.set(cg, p, t);
			                messageQ.sendMessage(master, req);
			                messageQ.flush();
			                ++cnt;
	                    }
	                }
	                
	                System.err.println("Sent out " + cnt + " datum requests with sid=" + currentSid);
	                //long then = System.currentTimeMillis();
	                
	                messageQ.flush();
	                tracker.waitForTasks(2000L);
	                
	                //long now = System.currentTimeMillis();
	                // System.err.println("Completed in " + (now - then) + "ms");
	                
		            o = contributions.get(cg, pos);
	            } catch (Exception ex) {
	                throw new RuntimeException("Error fetching contribution", ex);
	            }
            }
            
            return o;
        }
            

        /* (non-Javadoc)
         * @see net.derkholm.nmica.matrix.ObjectMatrix1D#set(int, java.lang.Object)
         */
        public void set(int pos, Object v) {
            throw new RuntimeException("ContributionGophers are read-only!");
        }
        
    }
}
