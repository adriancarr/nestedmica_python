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
/*
 * Created on Aug 10, 2004
 */
package net.derkholm.nmica.demos;

import java.net.InetSocketAddress;
import java.nio.ByteBuffer;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

/**
 * @author thomas
 */
public class StreamServer {
    public static void main(String[] args) 
    		throws Exception
    {
        ServerSocketChannel ssc = ServerSocketChannel.open();
        ssc.configureBlocking(true);
        ssc.socket().bind(new InetSocketAddress(13335));
        while (true) {
            final SocketChannel sc = ssc.accept();
            sc.configureBlocking(true);
            Thread socketServiceThread = new Thread() {
                public void run() {
                    try {
                        ByteBuffer buffer = ByteBuffer.allocate(1024);
                        while (true) {
                            buffer.clear();
                            buffer.limit(4);
                            while (buffer.position() < 4) {
                                sc.read(buffer);
                            }
                            buffer.rewind();
                            int i = buffer.getInt();
                            // System.err.println("Read: " + i);
                            buffer.rewind();
                            buffer.putInt(i + 1);
                            buffer.rewind();
                            sc.write(buffer);
                        }
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            };
            socketServiceThread.start();
        }
    }
}
