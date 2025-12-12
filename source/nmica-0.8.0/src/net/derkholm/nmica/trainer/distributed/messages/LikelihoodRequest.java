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

package net.derkholm.nmica.trainer.distributed.messages;

import java.nio.ByteBuffer;

import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix1D;
import net.derkholm.nmica.utils.mq.CodingException;
import net.derkholm.nmica.utils.mq.Packable;

/**
 * @author thomas
 */
public class LikelihoodRequest implements Packable {
    private static final boolean TX_PACKED = true;
    private static final int PACKED_FLAG = 0x100000;
    
    public int wid;
    public short sid;
    public int facette;
    public int contributionGroup;
    public int datum;
    public Matrix1D weights;
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#getPackedSize()
     */
    public int getPackedSize() throws CodingException {
        if (TX_PACKED) {
            int words = (weights.size() / 30) + 1;
            return 4 + 2 + 2 + 2 + 2 + 4 + (4 * words);
        } else {
            return 4 + 2 + 2 + 2 + 2 + 4 + (4 * weights.size());
        }
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#pack(java.nio.ByteBuffer)
     */
    public void pack(ByteBuffer buffer) throws CodingException {
        buffer.putInt(wid);
        buffer.putShort(sid);
        buffer.putShort((short) facette);
        buffer.putShort((short) contributionGroup);
        buffer.putShort((short) datum);
        
        int wc = weights.size();
        if (TX_PACKED) {
            buffer.putInt(wc | PACKED_FLAG);
            /* if (weights instanceof BinaryMatrixSlice) {
            	BinaryMatrixSlice bms = (BinaryMatrixSlice) weights;
            	byte[] data = bms.getByteData();
            	int offset = bms.getByteOffset();
	            int word = 0;
	            int bits = 0;
	            for (int w = 0; w < wc; ++w) {
	                word = (word << 1) | data[offset + w];
	                ++bits;
	                if (bits == 30) {
	                    buffer.putInt(word);
	                    word = 0;
	                    bits = 0;
	                }
	            }
	            if (bits != 0) {
	                buffer.putInt(word);
	            }
            } else { */
	            int word = 0;
	            int bits = 0;
	            for (int w = 0; w < wc; ++w) {
	                word = (word << 1);
	                if (weights.get(w) != 0.0) {
	                    word |= 1;
	                }
	                
	                ++bits;
	                if (bits == 30) {
	                    buffer.putInt(word);
	                    word = 0;
	                    bits = 0;
	                }
	            }
	            if (bits != 0) {
	                buffer.putInt(word);
	            }
            // }
        } else {
	        buffer.putInt(wc);
	        for (int w = 0; w < wc; ++w) {
	            buffer.putFloat((float) weights.get(w));
	        }
        }
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#unpack(java.nio.ByteBuffer)
     */
    public void unpack(ByteBuffer buffer) throws CodingException {
        wid = buffer.getInt();
        sid = buffer.getShort();
        facette = buffer.getShort();
        contributionGroup = buffer.getShort();
        datum = buffer.getShort();
        int wl = buffer.getInt();
        if ((wl & PACKED_FLAG) != 0) {
        	wl = wl & (PACKED_FLAG - 1);
            weights = new SimpleMatrix1D(wl);
            for (int wb = 0; wb < wl; wb += 30) {
                int word = buffer.getInt();
                int bm = Math.min(29, wl - wb - 1);
                while (bm >= 0) {
                    if ((word & 1) != 0) {
                        weights.set(wb + bm, 1.0);
                    }
                    word = word >> 1;
                    --bm;
                }
            }
            // System.err.println("done!");
        } else {
            weights = new SimpleMatrix1D(wl);
            for (int w = 0; w < wl; ++w) {
                weights.set(w, buffer.getFloat());
            }
        }
    }
}
