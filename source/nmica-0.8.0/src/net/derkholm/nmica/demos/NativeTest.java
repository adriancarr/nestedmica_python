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
 * Created on May 20, 2005
 */
package net.derkholm.nmica.demos;

import java.nio.ByteBuffer;
import java.nio.IntBuffer;

import net.derkholm.nmica.utils.CliTools;

public class NativeTest {
    static {
        System.loadLibrary("nativetest");
    }
    
    private String mode = "array";
    private int repeats = 10000;
    private int elements = 100;
    
    public void setElements(int i) {
        this.elements = i;
    }
    
    public void setMode(String buffer) {
        this.mode = buffer;
    }

    public void setRepeats(int repeats) {
        this.repeats = repeats;
    }

    /**
     * @param args
     */
    public static void main(String[] args) 
        throws Exception
    {
        NativeTest app = new NativeTest();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }

    public void run(String[] args)
        throws Exception
    {
        if ("buffer".equals(mode)) {
            testBuffer();
        } else if ("bounce".equals(mode)) {
            testBounce();
        } else {
            testArray();
        }
    }
    
    public void testArray() {
        int[] array = new int[elements];
        for (int i = 0; i < array.length; ++i) {
            array[i] = i;
        }
        
        System.out.println(sumArray(array));
        for (int c = 1; c < repeats; ++c) {
            sumArray(array);
        }
    }
    
    public void testBuffer() {
        IntBuffer buffer = ByteBuffer.allocateDirect(elements << 2).asIntBuffer();
        for (int i = 0; i < elements; ++i) {
            buffer.put(i);
        }
        
        System.out.println(sumBuffer(buffer));
        for (int c = 1; c < repeats; ++c) {
            sumBuffer(buffer);
        }
    }
    
    public void testBounce() {
        int[] array = new int[elements];
        for (int i = 0; i < array.length; ++i) {
            array[i] = i;
        }
        IntBuffer buffer = ByteBuffer.allocateDirect(elements << 2).asIntBuffer();
        
        for (int c = 0; c < repeats; ++c) {
            buffer.clear();
            // for (int i = 0; i < array.length; ++i) {
            //    buffer.put(array[i]);
            // }
            buffer.put(array);
            sumBuffer(buffer);
        }
    }
    
    private native int sumArray(int[] array);
    private native int sumBuffer(IntBuffer buffer);
}
