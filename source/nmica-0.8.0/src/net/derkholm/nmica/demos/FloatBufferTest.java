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
 * Created on Jan 4, 2005
 */
package net.derkholm.nmica.demos;

import java.nio.ByteBuffer;

import net.derkholm.nmica.utils.CliTools;

/**
 * @author thomas
 */
public class FloatBufferTest {
    private int size = 30;
    private int cycles = 10000;
    private boolean direct = false;
    
    public void setDirect(boolean b) {
        this.direct = b;
    }
    public void setCycles(int cycles) {
        this.cycles = cycles;
    }
    public void setSize(int size) {
        this.size = size;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        FloatBufferTest app = new FloatBufferTest();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        ByteBuffer buffer;
        if (direct) {
            buffer = ByteBuffer.allocateDirect(size * 4);
        } else {
            buffer = ByteBuffer.allocate(size * 4);
        }
        for (int c = 0; c < cycles; ++c) {
            buffer.clear();
            for (int i = 0; i < size; ++i) {
                buffer.putFloat(3.142F);
            }
        }
    }
}
