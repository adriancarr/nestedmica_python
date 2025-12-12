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

import java.nio.ByteBuffer;

/**
 * Stub implementation of Packable, useful for implementing
 * messages with no payload.
 * 
 * @author thomas
 */
public class AbstractPackable implements Packable {
    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#getPackedSize()
     */
    public int getPackedSize() throws CodingException {
        return 0;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#pack(java.nio.ByteBuffer)
     */
    public void pack(ByteBuffer buffer) throws CodingException {
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.Packable#unpack(java.nio.ByteBuffer)
     */
    public void unpack(ByteBuffer buffer) throws CodingException {
    }
}
