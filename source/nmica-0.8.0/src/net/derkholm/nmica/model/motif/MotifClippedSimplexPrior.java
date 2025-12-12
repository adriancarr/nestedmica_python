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


import net.derkholm.nmica.seq.NMSimpleDistribution;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 * @author thomas
 */
public class MotifClippedSimplexPrior extends MotifUniformSimplexPrior {
	private static final long serialVersionUID = 6591858390403070238L;
	
	private double minWeight = 0;
    private double maxWeight = 1.0;
    
    /**
     * @param alpha
     * @param length
     */
    public MotifClippedSimplexPrior(FiniteAlphabet alpha, int minLength, int maxLength, int unmaskedLength, double minWeight, double maxWeight) 
    {
        super(alpha, minLength, maxLength, unmaskedLength);
        this.minWeight = minWeight;
        this.maxWeight = maxWeight;
    }

    public double probabilityColumn(Distribution d) 
    {
    	// Disabled support for anything other than NMSimpleDistributions to see if the
    	// instanceof was expensive.
    	
    	/* if (d instanceof NMSimpleDistribution) { */
    		double[] wa = ((NMSimpleDistribution) d)._getWeightsArray();
    		for (double w : wa) {
    			if (w < minWeight || w > maxWeight) {
                    return Double.NEGATIVE_INFINITY;
                }
    		}
    		return 0;
    	/* } else {
	        try {
	            for (Iterator<?> i = ((FiniteAlphabet) d.getAlphabet()).iterator(); i.hasNext(); ) {
	                Symbol s = (Symbol) i.next();
	                double w = d.getWeight(s);
	                if (w < minWeight || w > maxWeight) {
	                    return Double.NEGATIVE_INFINITY;
	                }
	            }
	            return 0;
	        } catch (IllegalSymbolException ex) {
	            throw new BioError("Assertion failed: couldn't iterator distribution", ex);
	        }
    	} */
    }
}
