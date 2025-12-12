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
 * Memento for a proposal in an MCMC sampling process.
 * The penalty term is used to ensure detailed balance. 
 *
 * @author Thomas Down
 */

public final class PenalizedVariate {
    private final Object parent;
    private final Object variate;
    private final double balancePenalty;
    private final ContributionSampler sampler;
    private final boolean isSilent;
    
    public PenalizedVariate(Object parent, Object variate, double balancePenalty, boolean isSilent, ContributionSampler sampler) {
        this.parent = parent;
        this.variate = variate;
        this.balancePenalty = balancePenalty;
        this.sampler = sampler;
        this.isSilent = isSilent;
    }
    
    /**
     * Flag to indicate whether or not a change could possibly effect likelihoods
     */
    
    public boolean isSilent() {
    	return isSilent;
    }
    
    /**
     * Return the previous state in the chain
     */
    
    public Object getParent() {
        return parent;
    }
    
    /**
     * Return the proposed next state in the chain.
     */
    
    public Object getVariate() {
        return variate;
    }
    
    /**
     * Return a penalty term to preserve detailed balance in the metropolis-hastings process.
     * This should be the log Q ratio of the metropolis-hastings acceptance criterion.
     * A value of zero implies that the process is perfectly reversable.
     */
    
    public double getBalancePenalty() {
        return balancePenalty;
    }
    
    /**
     * Return the sampler that created this variate.
     */
    
    public ContributionSampler getSampler() {
        return sampler;
    }
}
