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

package net.derkholm.nmica.utils.mq;

import java.io.IOException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.SocketAddress;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.CancelledKeyException;
import java.nio.channels.ClosedSelectorException;
import java.nio.channels.DatagramChannel;
import java.nio.channels.SelectionKey;
import java.nio.channels.Selector;
import java.util.*;
import java.util.concurrent.*;

/**
 * Object for asynchronously sending and receiving messages via a UPD socket.
 * Supports boxcaring of multiple messages into one packet, and tracking of
 * peers.
 * 
 * @author thomas
 * @version branch-jdk15 (massive re-writes).
 */
public class MessageQueue<T>
{
    private final static int STATE_NEW = 0;
    private final static int STATE_RUNNING = 1;
    private final static int STATE_DRAINING = 2;
    private final static int STATE_DEAD = 3;
    
    private int state = STATE_NEW;
    
    private DatagramChannel channel;
    private InetSocketAddress endpoint;
    private Selector selector;
    private SelectionKey key;
    
    private MessageCodec<T> codec;
    
    private ConcurrentMap<SocketAddress,Peer> peers = new ConcurrentHashMap<SocketAddress,Peer>();
    private Peer[] peerArray = new Peer[0];
    private int peerRoundRobinCounter = 0;
    
    private Queue<Message<T>> inQ = new LinkedList<Message<T>>();
    private Semaphore inQnonEmpty = new Semaphore(0);
    private MessageHandler<T> handler;
    private boolean handling = false;
    
    private int packetSize = 7900;
    private int flushQueueLength = 500;

    private int txBytes = 0;
    private int txPackets = 0;
    private int rxBytes = 0;
    private int rxPackets = 0;
    private boolean slug = false;
    
    private QMon monitor;
    
    // Held true while waiting on a write operation, therefore
    // generally not worth trying to wake the selector even if
    // a flush has been suggested.
    private boolean currentlyInWriteMode = false;
    
    public void setSynchronousHandler(MessageHandler<T> handler) {
    	this.handler = handler;
    }
    
    public void setSlug(boolean b) {
    	this.slug = b;
    }
    
    // access to throughput counters.

    public int getTxBytes() {
        return txBytes;
    }
    public int getRxBytes() {
        return rxBytes;
    }
    public int getTxPackets() {
        return txPackets;
    }
    public int getRxPackets() {
        return rxPackets;
    }
    
    /**
     * Memento for a remote system which is accessible via
     * this MessageQueue.
     */
    
    public static interface Peer {
    	public SocketAddress getEndpoint();
    	public void setUserData(Object o);
    	public Object getUserData();
    }
    
    private static final class PeerImpl<T> implements Peer {
        private final SocketAddress endpoint;
        private Object userData;
        Queue<T> queue;
        
        private PeerImpl(SocketAddress endpoint) {
            this.endpoint = endpoint;
            queue = new LinkedList<T>();
        }
        
        /**
         * Get the network endpoint associated with this Peer.
         */
        
        public SocketAddress getEndpoint() {
            return endpoint;
        }

        /**
         * Attach an arbitrary piece of end-user data to this peer object.
         */
        
        public void setUserData(Object o) {
            this.userData = o;
        }
        
        /**
         * Fetch the end-user data attached to this object
         * @return the userData, or <code>null</code> if it has not been set.
         */
        
        public Object getUserData() {
            return userData;
        }
    }
    
    /**
     * Message that has arrived via the queue.
     */
    
    public static final class Message<T> {
        private final Peer sender;
        private final T body;
        
        private Message(Peer sender, T body) {
            this.sender = sender;
            this.body = body;
        }
        
        public Peer getSender() {
            return sender;
        }
        
        public T getBody() {
            return body;
        }
    }
    
    /**
     * Set the <code>MessageCodec</code> object which is responsible for encoding an decoding messages sent via this queue.
     * 
     * @param codec
     */
    
    public void setCodec(MessageCodec<T> codec) {
        this.codec = codec;
    }
    
    /**
     * Return the current codec for this queue.
     * @return
     */
    
    public MessageCodec<T> getCodec() {
        return codec;
    }
    
    /**
     * Set a threshold beyond which a flush will be automatically initiated.  If the number
     * of outgoing messages for any peer exceeds this threshold, the queue manager will make
     * an attempt to start transmitting outgoing messages.
     */

    public void setFlushThreshold(int i) {
	this.flushQueueLength = i;
    }

    /**
     * Create a new MessageQueue which listens on the specified UDP port
     * @param port
     * @throws IOException
     */
    
    public synchronized void reset() 
    	throws Exception
    {
    	System.err.printf("Resetting mq%n");
    	selector.close();
    	selector = Selector.open();
        key = channel.register(selector, SelectionKey.OP_READ);
    }
    
    public MessageQueue(int port) 
    		throws IOException
    {
        super();
        channel = DatagramChannel.open();
        channel.socket().bind(new InetSocketAddress(port));
        channel.configureBlocking(false);
        int realPort = channel.socket().getLocalPort();
        endpoint = new InetSocketAddress(InetAddress.getLocalHost(), realPort);
        
        selector = Selector.open();
        key = channel.register(selector, SelectionKey.OP_READ);
        
        /*
        new Thread("Buzzer!") {
        	public void run() {
        		while (true) {
        			try {
        				Thread.sleep(13000L);
        				System.err.println("Running buzzer");
        				for (int i = 0; i < 50; ++i) {
        					selector.wakeup();
        					Thread.sleep(20L);
        				}
        			} catch (InterruptedException ex) {}
        		}
        	}
        }.start();
        
        */
    }
    
    /**
     * Create a new MessageQueue using some free port number chosen arbitrarily by the system.
     * @throws IOException
     */
    
    public MessageQueue() 
    		throws IOException
    {
        this(0);
    }
    
    /**
     * Get the peer object associated with a particular communication endpoint.
     * If there has been no previous communication with this endpoint, a new
     * peer is created.
     * 
     * @param peerAddr
     * @return
     */
    
    public Peer getPeer(SocketAddress peerAddr) 
    {
	    // This is threadsafe without synchronization.
        Peer p = peers.get(peerAddr);
		if (p == null) {
			p = new PeerImpl<T>(peerAddr);
			Peer oldPeer = peers.putIfAbsent(peerAddr, p);
			if (oldPeer != null) {
				p = oldPeer;
			} else {
				synchronized (this) {
					peerArray = peers.values().toArray(new Peer[0]);
				}
			}
		}
		return p;
    }
    
    /**
	 * Enqueue a message for transmission the the specified peer.
	 * 
	 * @param p
	 * @param body
	 */
    
    public void sendMessage(Peer p, T body)
            throws QueueDeadException 
    {
        if (state == STATE_DEAD) {
            throw new QueueDeadException();
        }

        boolean needsFlush;
        synchronized (this) {
	        Queue<T> outQ = ((PeerImpl<T>) p).queue;
	        outQ.add(body);
	        needsFlush = outQ.size() >= flushQueueLength;
        }
        if (needsFlush) {
        	flush();
        }
    }
        
    /**
     * Get the next inbound message from the queue, waiting until one is
     * ready if necessary.
     * 
     * @return a Message
     * @throws QueueDeadException if the queue is in the dead state.  
     *                            This includes situations where the queue dies 
     *                            after the <code>next</code> method is called
     *                            but before a message becomes available.
     */
    
    public Message<T> next() 
    	throws QueueDeadException 
    {
		while (true) {
			try {
				inQnonEmpty.acquire();
				synchronized (inQ) {
					if (!inQ.isEmpty()) {
						Message<T> m = inQ.remove();
						if (!inQ.isEmpty()) {
							inQnonEmpty.release();
						} else {
							flush();
						}
						return m;
					}
					if (state == STATE_DEAD) {
						inQnonEmpty.release(); // wake up another waiter, if
												// necessary
						throw new QueueDeadException();
					}
				}
			} catch (InterruptedException ex) {
			}
		}
	}

    public Message<T> next(long timeout)
	throws QueueDeadException, InterruptedException
    {
	if (inQnonEmpty.tryAcquire(timeout, TimeUnit.MICROSECONDS)) {
	    synchronized (inQ) {
		Message<T> m = inQ.remove();
		if (inQ.poll() != null) {
		    inQnonEmpty.release();
		} else if (m != null && state == STATE_DEAD) {
		    inQnonEmpty.release();
		    throw new QueueDeadException();
		}
		return m;
	    }
	}
	return null;
    }

    public Message<T> nextNonBlocking()
	throws QueueDeadException
    {
	if (inQnonEmpty.tryAcquire()) {
	    synchronized (inQ) {
		Message<T> m = inQ.remove();
		if (inQ.poll() != null) {
		    inQnonEmpty.release();
		} else {
		    flush();
                }

		if (m == null && state == STATE_DEAD) {
		    inQnonEmpty.release();
		    throw new QueueDeadException();
		}
		return m;
	    }
	}
	return null;
    }

    public void next(Collection<MessageQueue.Message<T>> sink)
		throws QueueDeadException
    {
		while (true) {
		    try {
			    long before = System.currentTimeMillis();
				inQnonEmpty.tryAcquire(100, TimeUnit.MILLISECONDS);
				long after = System.currentTimeMillis();
				// if ((after - before) > 50) {
				// 	System.err.printf("inQnonEmpty aquisition took %dms%n", (after - before));
				// }
				synchronized (inQ) {
				    boolean sunk = false;
				    while (!inQ.isEmpty()) {
						sink.add(inQ.remove());
						sunk = true;
				    }
				    if (!sunk && state == STATE_DEAD) {
						inQnonEmpty.release();
						throw new QueueDeadException();
				    }
				    flush();
				}
				return;
		    } catch (InterruptedException ex) {}
		}
    }
    
    /**
     * Attempt to inititate the transmission of outbound messages which
     * are currently in the queue.  This method is just a hint to the queue,
     * and may return without doing anything.
     */
    
    private long wakeupTime = -1;
    
    
    public void flush() {
    	try {
    		if (monitor.doOutput()) {
    			if (slug) {
    				System.err.println("Output still pending, going for a wakeup");
    			}
    			wakeupTime = System.currentTimeMillis();
    			selector.wakeup();
    		}
    	} catch (Exception ex) {
    		ex.printStackTrace();
    	}
    	/*
    	int maxLength = 0, totLength = 0;
    	synchronized (this) {
			for (Queue<T> q : outQbyPeer.values()) {
				maxLength = Math.max(maxLength, q.size());
				totLength += q.size();
			}
    	}
    	*/
    	/*
    	if (Math.random() < 0.001) {
    		int maxLength = 0, totLength = 0;
    		for (BlockingQueue<T> q : outQbyPeer.values()) {
    			maxLength = Math.max(maxLength, q.size());
    			totLength += q.size();
    		}
    		// System.err.printf("Flushing outq, maxLength=%d, totLength=%d%n", maxLength, totLength);
    		if (maxLength == 0) {
    			new Exception("Flushing empty queue").printStackTrace();
    		}
    	}
    	*/
		// System.err.printf("Flushing outq, maxLength=%d, totLength=%d, cwm=%s%n", maxLength, totLength, currentlyInWriteMode);
    	/*
        if (!currentlyInWriteMode && totLength >= 1 ) {
        	try {
				if (monitor.doOutput()) {
					wakeupTime = System.currentTimeMillis();
				    selector.wakeup();
				} else {
					// System.err.println("Barging send was successful!");
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }*/
    }
    
    /**
     * Start sending and receiving messages via this queue.
     *
     */
    
    public synchronized void start() {
        if (state != STATE_NEW) {
            throw new IllegalStateException("Only newly created MessageQueues can be started");
        }
        monitor = new QMon();
        monitor.start();
    }
    
    /**
     * Wait until all current outgoing messages have been sent, then shut down
     * the thread which handles messaging for this queue.
     *
     */
    
    public synchronized void shutdown() {
        if (state != STATE_RUNNING) {
            throw new IllegalStateException("MessageQueue is not currently running, can't shutdown");
        }
        state = STATE_DRAINING;
    }
    
    /**
     * Get the network endpoint on which this message queue is listening.
     * @return
     */
    
    public InetSocketAddress getEndpoint() {
        return endpoint;
    }
    
    private class QMon extends Thread {
    	ByteBuffer inBuffer;
    	ByteBuffer outBuffer;
    	Peer currentDestination = null;
    	private long packetPrepTime = -1;
    	private int sendAttempts = 0;
    	
    	public synchronized boolean doOutput()
    		throws Exception
    	{
    		if (state != STATE_RUNNING) {
    			return false;
    		}
    		
    		int cnt = 0;
    		boolean isClear = true; // we hope!
    		while (isClear && prepareOutgoing() && cnt < 3) {
    			isClear = sendOutgoing();
    			if (isClear) {
    				++cnt;
    			}
    		}
    		if (cnt > 3) {
    			System.err.printf("Send %d packets%n", cnt);
    		}
    		return currentDestination != null;
    	}
    	
    	private boolean prepareOutgoing() 
    		throws CodingException 
    	{
    		synchronized (MessageQueue.this) {
	    		if (currentDestination != null) {
	    			return true;
	    		}
	    		
	    		Queue<T> outQ = null;
				int pal = peerArray.length;
		      FIND_SENDER_LOOP: 
		    	for (int prr = 1; prr <= pal; ++prr) {
					int t = peerRoundRobinCounter + prr;
					while (t >= pal) {
						t -= pal;
					}
					outQ = ((PeerImpl<T>) peerArray[t]).queue;
					if (outQ != null && !outQ.isEmpty()) {
						currentDestination = peerArray[t];
						peerRoundRobinCounter = t;
						break FIND_SENDER_LOOP;
					}
				}
				if (outQ != null && !outQ.isEmpty()) {
					outBuffer.clear();
					int numMessages = 0;
					T o;
					PACK_LOOP: while ((o = outQ.peek()) != null) {
						int size = codec.sizeMessage(o);
						int totalPacket = outBuffer.position() + size;
						if (numMessages == 0 || (totalPacket < packetSize)) {
							codec.writeMessage(outBuffer, o);
							outQ.remove();
							++numMessages;
						} else {
							break PACK_LOOP;
						}
					}
					outBuffer.flip();
					packetPrepTime = System.currentTimeMillis();
					sendAttempts = 0;
					return true;
				} else {
					return false;
				}
    		}
    	}
    	
    	private boolean sendOutgoing() 
    		throws IOException
    	{
    		int bytesSent = channel.send(
    				outBuffer,
					currentDestination.getEndpoint()
			);
    		++sendAttempts;
			if (bytesSent > 0) {
				currentDestination = null;
				++txPackets;
				txBytes += bytesSent;
				
				/*
				long delay = System.currentTimeMillis() - packetPrepTime;
				if (delay > 5) {
					System.err.printf("Packet delayed by %dms%n", delay);
				}
				if (sendAttempts > 1) {
					System.err.printf("Packet send at %d attempt%n", sendAttempts);
				}
				*/
				return true;
			} else {
				// System.err.printf("Packet postponed%n");
				return false;
			}
    	}
    	
		public void run() {
			try {
				synchronized (MessageQueue.this) {
					if (state != STATE_NEW) {
						throw new IllegalStateException();
					}
					
					inBuffer = ByteBuffer.allocateDirect(1 << 16);
					outBuffer = ByteBuffer.allocateDirect(1 << 16);
					inBuffer.order(ByteOrder.LITTLE_ENDIAN);
					outBuffer.order(ByteOrder.LITTLE_ENDIAN);
					
					state = STATE_RUNNING;
				}

		      ETERNITY_LOOP: 
				while (true) {
					try {
						boolean shutdownMode;
						boolean wantToWrite = prepareOutgoing();
						shutdownMode = (state == STATE_DRAINING);
	
						if (!wantToWrite && currentDestination == null && shutdownMode) {
							break ETERNITY_LOOP;
						}
	
						if (wantToWrite) {
							key.interestOps(SelectionKey.OP_READ | SelectionKey.OP_WRITE);
							currentlyInWriteMode = true;
						} else {
							key.interestOps(SelectionKey.OP_READ);
							currentlyInWriteMode = false;
						}
	
						long preTime = 0;
						wakeupTime = -1;
						if (slug) {
							System.err.printf("SLUG: going to sleep (write=%s)%n",
									currentlyInWriteMode);
							preTime = System.currentTimeMillis();
						}
						int numKeys = selector.select(200L);
						if (slug) {
							if (wakeupTime > 0) {
								System.err.printf("SLUG: tried to awake %d ago%n",
										System.currentTimeMillis() - wakeupTime);
							}
							System.err.printf("SLUG: sleepTime=%d%n", System
									.currentTimeMillis()
									- preTime);
						}
						if (wakeupTime > 0 && (System.currentTimeMillis() - wakeupTime > 80L)) {
							reset();
						}
						wakeupTime = -1;
	
						for (Iterator<SelectionKey> ki = selector.selectedKeys().iterator(); ki.hasNext();) {
							SelectionKey sk = ki.next();
							if (sk != key) {
								System.err.printf("Bad key%n");
								// Something's wrong here, but let's ignore the
								// erroneous key and carry on.
								continue;
							}
							if (slug) {
								System.err.printf("SLUG: got key, writable=%s readable=%s%n", key.isWritable(), key.isReadable());
							}
	
							if (key.isWritable()) {
								doOutput();
							}
	
							if (key.isReadable()) {
								List<Message<T>> received = new ArrayList<Message<T>>();
								RX_PACKET_LOOP: while (true) {
									inBuffer.clear();
									SocketAddress address = channel
											.receive(inBuffer);
									inBuffer.flip();
	
									if (inBuffer.limit() == 0) {
										break RX_PACKET_LOOP;
									}
	
									++rxPackets;
									rxBytes += inBuffer.limit();
	
									Peer p = getPeer(address);
									while (inBuffer.position() < inBuffer.limit()) {
										T body = codec.readMessage(inBuffer);
										Message<T> msg = new Message<T>(p, body);
										received.add(msg);
									}
								}
	
								if (handler != null) {
									handling = true;
									int rc = received.size();
									for (int ri = 0; ri < rc; ++ri) {
										Message<T> msg = received.get(ri);
										handler.handleMessage(MessageQueue.this, msg);
									}
									handling = false;
								} else {
									// anything that changes the size of inQ needs to hold this at the
									// moment.
									synchronized (inQ) {
										int oldSize = inQ.size();
										inQ.addAll(received);
										if (oldSize == 0 && !received.isEmpty()) {
											// ensure that a single permit is available when we release the inQ mutexQ.
											inQnonEmpty.drainPermits();
											inQnonEmpty.release();
										}
									}
								}
							}
	
							ki.remove();
						}
					} catch (ClosedSelectorException e) {
						// System.err.printf("Selector was replaced!%n");
					} catch (CancelledKeyException e) {
						// System.err.printf("Cancelled key (avoid this with better syncing?)%n");
					}
				}
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
	}
}
