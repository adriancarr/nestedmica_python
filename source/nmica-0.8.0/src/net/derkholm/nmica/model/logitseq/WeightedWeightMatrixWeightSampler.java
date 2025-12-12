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

package net.derkholm.nmica.model.logitseq;

import java.io.Serializable;

import net.derkholm.nmica.maths.Gaussian;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.PenalizedVariate;

/**
 * Sample the weight component of a WeightedWeightMatrix. 
 * 
 * @author thomas
 */
public class WeightedWeightMatrixWeightSampler implements ContributionSampler, Serializable {
    private static final long serialVersionUID = 3066772733918120817L;
    
    public WeightedWeightMatrixWeightSampler() {
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionSampler#sample(java.lang.Object, net.derkholm.nmica.matrix.ObjectMatrix1D)
     */
    public PenalizedVariate sample(final Object parent, final ObjectMatrix1D uncles)
            throws IllegalArgumentException 
    {
        WeightedWeightMatrix wwm = (WeightedWeightMatrix) parent;
        return new PenalizedVariate(
                parent,
                new WeightedWeightMatrix(
                        wwm.getWeightMatrix(),
                        wwm.getWeight() + Gaussian.standardVariate()
                ),
                0,
                false,
                this
        );
    }
}
