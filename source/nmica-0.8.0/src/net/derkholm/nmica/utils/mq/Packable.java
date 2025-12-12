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
 * API for objects that know how to encode and decode themselves. Typically,
 * instances of Packable will have a no-args constructor used prior to
 * unpacking.
 *
 * @author thomas
 * @author Matthew Pocock
 */
public interface Packable
{
	/**
	 * Calculate the number of bytes needed to pack this instance.
	 *
	 * @return  the number of bytes needed to back this instance
	 * @throws CodingException    if there was a problem calculating the size
	 *    needed
	 */
	public int getPackedSize()
	        throws CodingException;

	/**
	 * Pack this instance into the buffer. This should write exactly
	 * getPackedSize() bytes.
	 *
	 * @param buffer  the ByteBuffer to write to
	 * @throws CodingException    if there was a problem writing this object to
	 *    buffer
	 */
	public void pack(ByteBuffer buffer)
	        throws CodingException;

	/**
	 * Unpack this instance from the buffer. This is expected to consume as many
	 * bytes as returned by getPackedSize(). The previous state of the object will
	 * be replaced by the unpacked state.
	 *
	 * @param buffer  the ByteBuffer to read from
	 * @throws CodingException    if there was a problem reding this object from
	 *    the buffer
	 */
	public void unpack(ByteBuffer buffer)
	        throws CodingException;
}
