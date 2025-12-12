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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.ByteBuffer;

/**
 * @author thomas
 */
public class PackingTools {
    private PackingTools() {
        
    }
    
    public static byte[] stream(Object o)
    		throws CodingException
    {
        try {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ObjectOutputStream oos = new ObjectOutputStream(baos);
            oos.writeObject(o);
            oos.close();
            return baos.toByteArray();
        } catch (Exception ex) {
            throw new CodingException("Couldn't stream object", ex);
        }
    }
    
    public static void stream(Object o, ByteBuffer b)
    		throws CodingException
    {
        byte[] s = stream(o);
        b.putInt(s.length);
        b.put(s);
    }
    
    public static Object unstream(byte[] b)
    		throws CodingException
    {
        try {
            ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(b));
            Object o = ois.readObject();
            return o;
        } catch (Exception ex) {
            throw new CodingException("Couldn't unpack", ex);
        }
    }
    
    public static Object unstream(ByteBuffer buffer)
    		throws CodingException
    {
        int length = buffer.getInt();
        byte[] b = new byte[length];
        buffer.get(b);
        return unstream(b);
    }
}
