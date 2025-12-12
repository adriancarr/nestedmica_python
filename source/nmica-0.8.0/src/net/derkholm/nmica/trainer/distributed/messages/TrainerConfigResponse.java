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

import net.derkholm.nmica.model.FacetteMap;
import net.derkholm.nmica.utils.mq.CodingException;
import net.derkholm.nmica.utils.mq.Packable;
import net.derkholm.nmica.utils.mq.PackingTools;

/**
 * @author thomas
 */
public class TrainerConfigResponse implements Packable {
    public int components;
    public int dataSetSize;
    public FacetteMap facetteMap;
    private byte[] facetteMapStream;
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#getPackedSize()
     */
    public int getPackedSize() 
    		throws CodingException
    {
        streamObjects();
        return 4 + 4 + 4 + facetteMapStream.length;
    }
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#pack(java.nio.ByteBuffer)
     */
    public void pack(ByteBuffer buffer) throws CodingException {
        streamObjects();
        buffer.putInt(components);
        buffer.putInt(dataSetSize);
        buffer.putInt(facetteMapStream.length);
        buffer.put(facetteMapStream);
    }
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#unpack(java.nio.ByteBuffer)
     */
    public void unpack(ByteBuffer buffer) throws CodingException {
        components = buffer.getInt();
        dataSetSize = buffer.getInt();
        facetteMap = (FacetteMap) PackingTools.unstream(buffer);
    }
    
    private void streamObjects() 
    		throws CodingException
    {
        if (facetteMapStream == null) {
            facetteMapStream = PackingTools.stream(facetteMap);
        }
    }
}
