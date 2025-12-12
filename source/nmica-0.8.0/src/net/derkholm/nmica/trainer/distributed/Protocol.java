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

package net.derkholm.nmica.trainer.distributed;

import net.derkholm.nmica.utils.mq.MessageCodec;
import net.derkholm.nmica.utils.mq.TaggedPackableCodec;
import net.derkholm.nmica.utils.mq.Packable;

import net.derkholm.nmica.trainer.distributed.messages.*;

/**
 * Constants for implementing the standard distributed training
 * wire protocol.
 * 
 * @author thomas
 */
public class Protocol {
    /**
     * The current protocol version, should be bumped whenever
     * the messages or tags change.
     */
        
   public final static int VERSION = 4;
    
    /**
     * Codec for sending packets over the wire.
     */
    
    public final static MessageCodec<Packable> CODEC;
    
    static {
        TaggedPackableCodec codec = new TaggedPackableCodec();
        codec.addMessageType((byte) 0, Ready.class);
        codec.addMessageType((byte) 0x1, NotReady.class);
        codec.addMessageType((byte) 0x2, Shutdown.class);
        codec.addMessageType((byte) 0x3, Flush.class);
        
        codec.addMessageType((byte) 0x10, TrainerConfigRequest.class);
        codec.addMessageType((byte) 0x11, TrainerConfigResponse.class);
        
        codec.addMessageType((byte) 0x20, LikelihoodRequest.class);
        codec.addMessageType((byte) 0x21, LikelihoodResponse.class);
        
        codec.addMessageType((byte) 0x30, DatumRequest.class);
        codec.addMessageType((byte) 0x31, DatumResponse.class);
        
        codec.addMessageType((byte) 0x40, ContributionRequest.class);
        codec.addMessageType((byte) 0x41, ContributionResponse.class);
        
        // codec.addMessageType((byte) 0x7f, LikelihoodRequestExptl.class);
        
        CODEC = codec;
    }
    
    private Protocol() {
    }
}
