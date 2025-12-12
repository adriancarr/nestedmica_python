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

import java.io.Serializable;

import net.derkholm.nmica.matrix.ObjectMatrix1D;

import org.biojava.bio.BioError;

/**
 * Implementation of ContributionSampler which randomly dispatches sample requests
 * to a set of underlying samplers.
 *
 * @author Thomas Down
 */

public class MultiplexContributionSampler implements ContributionSampler, Serializable {
	private static final long serialVersionUID = -581309554566836936L;
	
	private ContributionSampler[] samplers = new ContributionSampler[0];
    private double[] weights = new double[0];
    double totalWeight = 0.0;
    
    public MultiplexContributionSampler() {
    }
    
    /**
     * Add a sampler to the multiplex
     */
    
     public void addSampler(ContributionSampler cs) {
         addSampler(cs, 1.0);
     }
     
    public void addSampler(ContributionSampler cs, double weight) {
        int len = samplers.length;
        ContributionSampler[] _samplers = new ContributionSampler[len + 1];
        System.arraycopy(samplers, 0, _samplers, 0, len);
        double[] _weights = new double[len + 1];
        System.arraycopy(weights, 0, _weights, 0, len);
        _samplers[len] = cs;
        _weights[len] = weight;
        totalWeight += weight;
        samplers = _samplers;
        weights = _weights;
    }
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles) 
        throws IllegalArgumentException
    {
        double slice = Math.random() * totalWeight;
        for (int i = 0; i < samplers.length; ++i) {
            slice -= weights[i];
            if (slice <= 0) {
                return samplers[i].sample(parent, uncles);
            }
        }
        
        throw new BioError("Assertion failed: didn't pick a ContributionSampler");
    }
}
