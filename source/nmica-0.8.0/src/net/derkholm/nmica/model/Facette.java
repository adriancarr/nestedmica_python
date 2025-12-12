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

package net.derkholm.nmica.model;

/**
 * Definition for one facette of a multi-ICA dataset.
 *
 * @author Thomas Down
 */

public interface Facette {
    /**
     * Get the type of objects which have defined likelihoods under this
     * facette's model.
     */
    
    public Class getDataType();
    
    /**
     * Get the type of objects which can contribute to this facette.
     */
    
    public Class getContributionType();
    
    /**
     * Determine if a given value in the mixing matrix implies that a contribution is
     * decoupled under this facette.
     */
    
    public boolean isContributionDecoupled(double weight);
    
    /**
     * Return an object for calculating the likelihood of some piece of data under
     * this model.
     *
     * @throws IllegalArgumentException if the data is not of an appropriate type.
     */
    
    public LikelihoodCalculator getLikelihoodCalculator(Object data);
}
