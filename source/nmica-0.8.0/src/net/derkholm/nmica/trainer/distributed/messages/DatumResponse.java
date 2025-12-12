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

import net.derkholm.nmica.utils.mq.CodingException;
import net.derkholm.nmica.utils.mq.Packable;
import net.derkholm.nmica.utils.mq.PackingTools;

/**
 * @author thomas
 */
public class DatumResponse implements Packable {
    public int datumIndex;
    public int facette;
    public Object datum;
    private byte[] datumStream;
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#getPackedSize()
     */
    public int getPackedSize() throws CodingException {
        streamObjects();
        return 4 + 4 + 4 + datumStream.length;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#pack(java.nio.ByteBuffer)
     */
    public void pack(ByteBuffer buffer) throws CodingException {
        streamObjects();
        buffer.putInt(datumIndex);
        buffer.putInt(facette);
        buffer.putInt(datumStream.length);
        buffer.put(datumStream);
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#unpack(java.nio.ByteBuffer)
     */
    public void unpack(ByteBuffer buffer) throws CodingException {
        datumIndex = buffer.getInt();
        facette = buffer.getInt();
        datum = PackingTools.unstream(buffer);
    }

    private void streamObjects() 
		throws CodingException
	{
        if (datumStream == null) {
            datumStream = PackingTools.stream(datum);
        }
    }
}
