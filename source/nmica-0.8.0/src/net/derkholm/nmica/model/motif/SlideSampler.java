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

package net.derkholm.nmica.model.motif;

import java.io.Serializable;

import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.PenalizedVariate;

import org.biojava.bio.BioError;

/**
 * Motif sampler which slides the active window.
 *
 * @author Thomas Down
 */

public class SlideSampler implements ContributionSampler, Serializable {
    private final static long serialVersionUID = 1000000L;
    
    public SlideSampler() {
    }
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles)
    {
        try {
            NMWeightMatrix wmp = (NMWeightMatrix) parent;
            int offset = Math.random() < 0.5 ? 1 : -1;
            return new PenalizedVariate(parent, new NMWeightMatrix(wmp.rawColumnArray(), wmp.columns(), wmp.offset() + offset), 0.0, false, this);
        } catch (Exception ex) {
            throw new BioError("Assertion failed: error sampling preference", ex);
        }
    }
    
    public String toString() {
        return "SlideSampler()";
    }
}
