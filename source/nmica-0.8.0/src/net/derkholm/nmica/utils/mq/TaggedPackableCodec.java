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

import java.lang.reflect.Constructor;
import java.nio.ByteBuffer;
import java.util.HashMap;
import java.util.Map;

/**
 * General purpose codec for transmitting <code>Packable</code>
 * messages. This codec places a single-byte tag at the beginning of each
 * object and follows that by the
 * 
 * @author thomas
 */
public class TaggedPackableCodec
        implements MessageCodec<Packable>
{
    private final Map<Class<? extends Packable>,Byte> classToTag;
    private final Map<Byte,Class<? extends Packable>> tagToClass;
    private final Map<Byte,Constructor<? extends Packable>> tagToCt;
    
    public TaggedPackableCodec() {
        classToTag = new HashMap<Class<? extends Packable>,Byte>();
        tagToClass = new HashMap<Byte,Class<? extends Packable>>();
        tagToCt = new HashMap<Byte,Constructor<? extends Packable>>();
    }
    
    public void addMessageType(byte tag, Class<? extends Packable> messageClass) 
    		throws IllegalArgumentException
    {
        Byte tagKey = new Byte(tag);
        if (classToTag.containsKey(messageClass)) {
            throw new IllegalArgumentException("Message class " + messageClass.getName() + " is already registered");
        }
        if (tagToClass.containsKey(tagKey)) {
            throw new IllegalArgumentException("Tag " + tag + " already has an associated message type");
        }
        if (!Packable.class.isAssignableFrom(messageClass)) {
            throw new IllegalArgumentException("Message class " + messageClass.getName() + " is not Packable");
        }
        
        try {
            Constructor<? extends Packable> con = messageClass.getConstructor();
            tagToCt.put(tagKey, con);
            tagToClass.put(tagKey, messageClass);
            classToTag.put(messageClass, tagKey);
        } catch (Exception ex) {
            throw new IllegalArgumentException("Message class " + messageClass.getName() + " does not have an accessible no-args constructor");
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.MessageCodec#sizeMessage(java.lang.Object)
     */
    public int sizeMessage(Packable pmsg)
            throws CodingException
    {
        return 1 + pmsg.getPackedSize();
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.MessageCodec#writeMessage(java.nio.ByteBuffer, java.lang.Object)
     */
    public void writeMessage(ByteBuffer buffer, Packable pmsg) throws CodingException {
        Class msgClass = pmsg.getClass();
        Byte tagKey = classToTag.get(msgClass);
        if (tagKey == null) {
            throw new CodingException("Message class " + msgClass + " is not supported");
        }
        buffer.put(tagKey.byteValue());
        pmsg.pack(buffer);
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.mq.MessageCodec#readMessage(java.nio.ByteBuffer)
     */
    public Packable readMessage(ByteBuffer buffer)
    		throws CodingException
    {
        byte tag = buffer.get();
        Byte tagKey = new Byte(tag);
        Constructor<? extends Packable> ct = tagToCt.get(tagKey);
        if (ct == null) {
            throw new CodingException("Unrecognized tag " + tag + " in stream");
        }
        Packable message;
        try {
            message = ct.newInstance();
        } catch (Exception ex) {
            throw new CodingException("Unable to instantiate message");
        }
        message.unpack(buffer);
        return message;
    }

}
