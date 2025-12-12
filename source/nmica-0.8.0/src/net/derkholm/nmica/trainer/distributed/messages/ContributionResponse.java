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

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.utils.mq.CodingException;
import net.derkholm.nmica.utils.mq.Packable;
import net.derkholm.nmica.utils.mq.PackingTools;

/**
 * @author thomas
 */
public class ContributionResponse implements Packable {
    private static final byte PAYLOAD_SERIALIZED = 0;
    private static final byte PAYLOAD_WM = 0x10;
    
    private static final Symbol[] SYMBOL_INDICES;
    
    static {
        try {
            SYMBOL_INDICES = (Symbol[]) DNATools.createDNA("acgt").toList().toArray(new Symbol[0]);
        } catch (IllegalSymbolException e) {
            throw new Error(e);
        }
    }
    
    public short sid;
    public int component;
    public int contributionGroup;
    public Object contribution;
    
    private byte[] contributionStream;
    private float[] wmContributionCache;

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#getPackedSize()
     */
    public int getPackedSize() throws CodingException {
        if (contribution instanceof WeightMatrix) {
            return 2 + 2 + 2 + 1 + 2 + (12 * ((WeightMatrix) contribution).columns()); 
        } else {
            streamObjects();
            return 2 + 2 + 2 + 1 + 4 + contributionStream.length;
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#pack(java.nio.ByteBuffer)
     */
    public void pack(ByteBuffer buffer) throws CodingException {
        buffer.putShort(sid);
        buffer.putShort((short) component);
        buffer.putShort((short) contributionGroup);
        if (contribution instanceof WeightMatrix) {
            WeightMatrix wm = (WeightMatrix) contribution;
            if (wmContributionCache == null) {
                int cols = wm.columns();
                float[] _wmContributionCache = new float[cols * 3];
                int z = 0;
                try {
	                for (int c = 0; c < wm.columns(); ++c) {
	                    Distribution col = wm.getColumn(c);
		                for (int i = 0; i < SYMBOL_INDICES.length - 1; ++i) {
		                    _wmContributionCache[z++] = (float) col.getWeight(SYMBOL_INDICES[i]);
		                }
	                }
                } catch (IllegalSymbolException ex) {
                    throw new CodingException(ex);
                }
                wmContributionCache = _wmContributionCache;
            }
            
            buffer.put(PAYLOAD_WM);
            buffer.putShort((short) wm.columns());
            for (int z = 0; z < wmContributionCache.length; ++z) {
	            buffer.putFloat(wmContributionCache[z]);
            }
        } else {
            streamObjects();
            buffer.put(PAYLOAD_SERIALIZED);
            buffer.putInt(contributionStream.length);
            buffer.put(contributionStream);
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#unpack(java.nio.ByteBuffer)
     */
    public void unpack(ByteBuffer buffer) throws CodingException {
        sid = buffer.getShort();
        component = buffer.getShort();
        contributionGroup = buffer.getShort();
        byte type = buffer.get();
        if (type == PAYLOAD_SERIALIZED) {
            contribution = PackingTools.unstream(buffer);
        } else if (type == PAYLOAD_WM) {
            int columns = buffer.getShort();
            try {
                Distribution[] cols = new Distribution[columns];
                for (int c = 0; c < columns; ++c) {
                    cols[c] = new NMSimpleDistribution(DNATools.getDNA());
                    double tot = 0;
                    for (int i = 0; i < SYMBOL_INDICES.length - 1; ++i) {
                        double x = buffer.getFloat();
                        cols[c].setWeight(SYMBOL_INDICES[i], x);
                        tot += x;
                    }
                    cols[c].setWeight(SYMBOL_INDICES[SYMBOL_INDICES.length - 1], 1.0 - tot);
                }
                contribution = new SimpleWeightMatrix(cols);
            } catch (Exception ex) {
                throw new CodingException(ex);
            }
        } else {
            throw new CodingException("Unknown payload type " + type);
        }
    }

    private void streamObjects() 		
    		throws CodingException
	{
        if (contributionStream == null) {
            contributionStream = PackingTools.stream(contribution);
        }
    }
    
}
