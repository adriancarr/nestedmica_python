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

import java.io.Serializable;

import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.maths.*;

/**
 * Memento which defines a likelihood calculation operation in a distributed
 * nstrainer run.  WorkUnits should only be interpreted in the context of
 * the <code>TrainingServices</code> instance which issued them.  They are
 * invalidated by the issue of WorkUnits with higher SIDs.
 * 
 * @author Thomas Down
 */
final class WorkUnit implements Serializable {
    int assignedWorkerID = -1;
    
    private static final long serialVersionUID = 7409264031892618335L;
    
    /**
     * Work unit ID.
     */
    
    public final int wid;
    
    /**
     * Work session ID.
     */
    
    public final int sid;
    
    /**
     * Index of requested facette in the trainer's FacetteMap
     */
    
    public final int facette;
    
    /**
     * Index of the requested ContributionGroup in the trainer's FacetteMap.  This
     * must be linked to the requested Facette.
     */
    
    public final int contributionGroup;
    
    /**
     * Index of the requested Datum in the trainer's dataset.
     */
    
    public final int datum;
    
    /**
     * Weight vector for this likelihood cell.
     */
    
    public final Matrix1D weights;

    /**
     * Completion procedure
     */

    public DoubleProcedure writeback;    


    public WorkUnit(int wid, int sid, int facette, int contributionGroup, int datum, Matrix1D weights, DoubleProcedure writeback) {
        this.wid = wid;
        this.sid = sid;
        this.facette = facette;
        this.contributionGroup = contributionGroup;
        this.datum = datum;
        this.weights = weights;
        this.writeback = writeback;
    }
}
