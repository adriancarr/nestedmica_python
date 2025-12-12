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

package net.derkholm.nmica.utils.tracker;

/**
 * General-purpose Tracker implementation.
 * 
 * @author thomas
 */
public class SimpleTracker implements Tracker {
    private int tasks = 0;
    private int completed = 0;
    private TaskFailedException failure = null;
    private Object lock = new Object();
    private boolean active = false;
    private long wakeupInstant;
    
    private class StTask implements Task {
        private boolean isComplete = false;
        private Object data;
        private Object result;
        
        public void setData(Object o) {
            this.data = o;
        }
        
        public Object getData() {
            return data;
        }
        
        /* (non-Javadoc)
         * @see net.derkholm.nmica.utils.tracker.Task#getTracker()
         */
        public Tracker getTracker() {
            return SimpleTracker.this;
        }

        /* (non-Javadoc)
         * @see net.derkholm.nmica.utils.tracker.Task#completed()
         */
        public void completed() {
            synchronized (lock) {
                if (!isComplete) {
	                isComplete = true;
	                ++completed;
	                if (completed >= tasks) {
	                    // wakeupInstant = System.currentTimeMillis();
	                    lock.notify();
	                }
                }
            }
        }
        
        /* (non-Javadoc)
         * @see net.derkholm.nmica.utils.tracker.Task#completed()
         */
        public void completed(Object result) {
            synchronized (lock) {
                if (!isComplete) {
	                isComplete = true;
	                ++completed;
	                this.result = result;
	                if (completed >= tasks) {
	                    // wakeupInstant = System.currentTimeMillis();
	                    lock.notify();
	                }
                }
            }
        }
        
        public boolean isComplete() {
            return isComplete;
        }

        public Object getResult() {
            return result;
        }
        
        /* (non-Javadoc)
         * @see net.derkholm.nmica.utils.tracker.Task#failed()
         */
        public void failed() {
            synchronized (lock) {
                failure = new TaskFailedException();
                lock.notifyAll();
            }
        }

        /* (non-Javadoc)
         * @see net.derkholm.nmica.utils.tracker.Task#failed(java.lang.Throwable)
         */
        public void failed(Throwable cause) {
            synchronized (lock) {
                failure = new TaskFailedException(cause);
                lock.notifyAll();
            }
        }
        
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.tracker.Tracker#newTask()
     */
    public Task newTask() {
        synchronized (lock) {
            if (active) {
                throw new IllegalStateException("Can't add extra tasks once tracking has started");
            }
            
            ++tasks;
            return new StTask();
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.tracker.Tracker#totalTasks()
     */
    public int totalTasks() {
        synchronized (lock) {
            return tasks;
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.tracker.Tracker#pendingTasks()
     */
    public int pendingTasks() {
        synchronized (lock) {
            return tasks - completed;
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.tracker.Tracker#waitForTasks()
     */
    public void waitForTasks() throws TaskFailedException {
        synchronized (lock) {
            active = true;
            
            while (tasks != completed && failure == null) {
                try {
                    lock.wait();
                } catch (InterruptedException ex) {}
            }
            // System.err.println("Tracker woken up, pingtime=" + (System.currentTimeMillis() - wakeupInstant));
            
            if (failure != null) {
                throw failure;
            }
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.tracker.Tracker#waitForTasks(long)
     */
    public boolean waitForTasks(long timeout) throws TaskFailedException {
        synchronized (lock) {
            active = true;
            
            if (tasks != completed && failure == null) {
                try {
                    lock.wait(timeout);
                } catch (InterruptedException ex) {}
            }
            
            if (failure != null) {
                throw failure;
            }
            
            return (completed == tasks);
        }
    }

}
