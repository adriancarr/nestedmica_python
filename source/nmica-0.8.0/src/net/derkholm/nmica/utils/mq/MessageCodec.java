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
 * Object which encapsulates some standard for transmitting messages in
 * datagrams. A MessageCodec is often restricted to some type of objects. For
 * example, it may only operate on serializable instances or on objects
 * inheriting from a particular base class. One codec is typically used for
 * each domain of communication, so typically a Codec will handle many different
 * types of object.
 * 
 * @author thomas
 * @author Matthew Pocock
 */
public interface MessageCodec<T>
{
	/**
	 * Calculate the size (in bytes) that the message would be for encoding obj.
	 *
	 * @param obj   the instance to serialize
	 * @return      the number of bytes needed to serialize it
	 * @throws CodingException  if there was an error calculating the size
	 *    required
	 */
	public int sizeMessage(T obj)
	        throws CodingException;

	/**
	 * Encode obj to the ByteBuffer. This should use exactly the number of bytes
	 * indicated by sizeMessage().
	 *
	 * @param buffer  the ByteBuffer to write obj to
	 * @param obj     the instance to encode
	 * @throws CodingException  if there was an error while writing obj to the
	 *    ByteBuffer
	 */
	public void writeMessage(ByteBuffer buffer, T obj)
	        throws CodingException;

  /**
   * Decode an object from the ByteBuffer. This is expected to consume as
   * many bytes as sizeMessage() indicates.
   *
   * @param buffer  the ByteBuffer to read an object from
   * @return        the instance read
   * @throws CodingException  if there was an error while reading the instance
   *    from the ByteBuffer
   */
	public T readMessage(ByteBuffer buffer)
          throws CodingException;
}
