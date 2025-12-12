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
import java.net.SocketAddress;
import java.nio.ByteBuffer;
import java.nio.channels.DatagramChannel;

/**
 * Simple ping server that listens to a datagram socket.
 * 
 * @author thomas
 */
public class DatagramServer {
    public static void main(String[] args) 
    		throws Exception
    {
        DatagramChannel channel = DatagramChannel.open();
        channel.socket().bind(new InetSocketAddress(13335));
        int port = channel.socket().getLocalPort();
        System.out.println("Port is: " + port);
        
        ByteBuffer buffer = ByteBuffer.allocateDirect(1 << 16);
        while (true) {
            buffer.clear();
            SocketAddress source = channel.receive(buffer);
            buffer.flip();
            int i = buffer.getInt();
            buffer.clear();
            buffer.putInt(i + 1);
            buffer.flip();
            channel.send(buffer, source);
        }
    }
}
