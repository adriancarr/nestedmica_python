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

package net.derkholm.nmica.apps;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.trainer.distributed.DistributedLikelihoodWorker;
import net.derkholm.nmica.trainer.distributed.Worker;

/**
 * @author thomas
 */

@App(overview="A worker process for distributed NestedMICA operation", generateStub=true)
@NMApp(launchName="nmworker", vm=VirtualMachine.SERVER)
public class MotifFinderNode {
    private String server;
    private int port = 1099;
    private int timeLimit = 0;
    private boolean nag = false;
    private int threads = 1;
    private int sleepy = 0;
    private int lruSize = -1;
    private boolean throughputMonitor = false;
    private boolean throughputWatchdog = true;
    
    @Option(help="Kill process if throughtput falls to zero", optional=true, userLevel=UserLevel.DEBUG)
    public void setThroughputWatchdog(boolean b) {
        this.throughputWatchdog = b;
    }
    
    @Option(help="Periodically log DLEP throughput", optional=true, userLevel=UserLevel.DEBUG)
    public void setThroughputMonitor(boolean b) {
        this.throughputMonitor = b;
    }
    
    @Option(help="The size of the worker nodes' contribution cache", optional=true, userLevel=UserLevel.EXPERT)
    public void setLruSize(int i) {
        this.lruSize = i;
    }

    @Option(help="Performance debugging", optional=true, userLevel=UserLevel.DEBUG)
    public void setSleepy(int i) {
        this.sleepy = i;
    }
    
    @Option(help="The number of threads to use to use.  For optimal performance, this should usually equal the number of CPU cores in your computer", optional=true)
    public void setThreads(int i) {
        this.threads = i;
    }
    
    public void setNag(boolean b) {
        this.nag = b;
    }
    
    @Option(help="The time (in seconds) for this node to run before automatically exiting", optional=true)
    public void setTimeLimit(int i) {
        this.timeLimit = i;
    }
    
    @Option(help="The name of the server where the main motiffinder process is running", optional=false)
    public void setServer(String s) {
        this.server = s;
    }
    
    @Option(help="The network port where the motiffinder process is listening.  Should match the '-port' option used when starting motiffinder", optional=false)
    public void setPort(int p) {
        this.port = p;
    }
    
    public void main(String[] args)
    		throws Exception
    {
        Worker w = null;
        do {
            w = DistributedLikelihoodWorker.connect(server, port, threads);
            if (sleepy > 0) {
                ((DistributedLikelihoodWorker) w).setSleepy(sleepy);
            }
            if (lruSize >= 0) {
                ((DistributedLikelihoodWorker) w).setLruSize(lruSize);
            }
            if (throughputMonitor) {
                ((DistributedLikelihoodWorker) w).setThroughputMonitor(throughputMonitor);
            }
            if (throughputWatchdog) {
                ((DistributedLikelihoodWorker) w).setThroughputWatchdog(throughputWatchdog);
            }
        } while (nag && w == null);
        if (w == null) {
            System.err.println("Couldn't connect to " + server + ":" + port);
            return;
        }
        w.start();

        if (timeLimit > 0) {
            long now = System.currentTimeMillis();
            long then = now + (1000L * timeLimit);
            while (now < then) {
                try {
                    Thread.sleep(then - now);
                } catch (Exception ex) {
                    // ignored.
                }
                now = System.currentTimeMillis();
            }
            w.stop();
        }
    }
}
