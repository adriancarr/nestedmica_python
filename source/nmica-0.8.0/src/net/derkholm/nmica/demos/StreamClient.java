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
import java.nio.channels.SocketChannel;

/**
 * @author thomas
 */
public class StreamClient {
    public static void main(String[] args) 
    		throws Exception
    {
        String serverHost = args[0];
        int pings = Integer.parseInt(args[1]);
        SocketChannel channel = SocketChannel.open();
        channel.configureBlocking(true);
        channel.connect(new InetSocketAddress(serverHost, 13335));
        
        ByteBuffer buffer = ByteBuffer.allocate(1024);
        
        long before = System.currentTimeMillis();
        for (int p = 0; p < pings; ++p) {
	        buffer.putInt(p);
	        buffer.rewind();
	        buffer.limit(4);
	        channel.write(buffer);
	        buffer.rewind();
	        while (buffer.position() < 4) {
	            channel.read(buffer);
	        }
	        buffer.rewind();
        }
        long after = System.currentTimeMillis();
        
        System.out.println("Stream roundtrip takes " + ((1.0 * (after - before)) / pings) + "ms");
    }
}
