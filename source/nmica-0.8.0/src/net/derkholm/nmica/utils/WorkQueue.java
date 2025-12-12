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

package net.derkholm.nmica.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.locks.AbstractQueuedSynchronizer;

/**
 * Simple queue for executing a group of <code>Runnable</code> instances in a multithreaded environment.
 *
 * @author Thomas Down
 */

public abstract class WorkQueue {    
    public static WorkQueue create(int cpus) {
        if (cpus < 2 && !Boolean.getBoolean("nmica.workqueue_always_thread")) {
            return new SimpleWorkQueue();
        } else {
            return new ThreadedWorkQueue(cpus);
        }
    }
    
    private static class SimpleWorkQueue extends WorkQueue {
        private List<Runnable> work = new ArrayList<Runnable>();
        
        public void add(Runnable r) {
            work.add(r);
        }
        
        public void flush() {
            for (Runnable r : work) {
                r.run();
            }
            work.clear();
        }
    }
    
    private static class ThreadedWorkQueue extends WorkQueue {
        private class Sync extends AbstractQueuedSynchronizer {

            private static final long serialVersionUID = 3761121626695480368L;
            
            public void increase() {
                while (true) {
                    int cur = getState();
                    if (compareAndSetState(cur, cur + 1)) {
                        return;
                    }
                }
            }
            
            public boolean tryReleaseShared(int i) {
                while (true) {
                    int cur = getState();
                    if (compareAndSetState(cur, cur - 1)) {
                        return cur == 1;
                    }
                }
            }
            
            public int tryAcquireShared(int i) {
                if (getState() == 0) {
                    return 1;
                } else {
                    return -1;
                }
            }
        }
        
        private BlockingQueue<Runnable> work = new LinkedBlockingQueue<Runnable>();
        private Worker[] workers;
        private Sync sync = new Sync();
        
        public ThreadedWorkQueue(int cpus) {
            if (cpus < 1) {
                cpus = 1;
            }
                workers = new Worker[cpus];
                for (int w = 0; w < workers.length; ++w) {
                    workers[w] = new Worker("" + w);
                    workers[w].start();
                }
        }
        
        public void add(Runnable r) {
            sync.increase();
            work.add(r);
        }
        
        public void flush() {  
            // if (watchdogTimer > 0) {
            //    try {
            //        if (!sync.tryAcquireSharedNanos(1, 1000000000L * watchdogTimer)) {
            //            System.err.println("Watchdog timeout in WorkQueue");
            //            for (Worker w : workers) {
            //                System.err.println(w.getName() + ": " + w.currentJob);
            //            }
            //        }
            //    } catch (InterruptedException e) {
            //        e.printStackTrace();
            //    }
            // } else {
                sync.acquireShared(1);
            // }
        }
        
        private class Worker extends Thread {
            private final String name;
            private Runnable currentJob;
            
            public Worker(String name) {
                super();
                setDaemon(true);
                this.name = name;
            }

            public void run() {
                while (true) {
                    try {
                        currentJob = work.take();
                        currentJob.run();
                        currentJob = null;
                        sync.releaseShared(1);
                    } catch (InterruptedException ex) {}
                }
            }
        }
    }

    public abstract void add(Runnable r);
    public abstract void flush();
}
