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

import net.derkholm.nmica.maths.MathsTools;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.PenalizedVariate;

import org.biojava.bio.BioError;

/**
 * Motif sampler which adjusts the length of a motif
 *
 * @author Thomas Down
 */

public class RetrimSampler implements ContributionSampler, Serializable {
    private final static long serialVersionUID = -2562505375552600554L;
    
    private final int minLength;
    private final int maxLength;
    
    public RetrimSampler(WeightMatrixPrior wmp) {
    	this.minLength = wmp.getMinLength();
    	this.maxLength = wmp.getMaxLength();
    }
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles)
    {
        try {
            NMWeightMatrix wmp = (NMWeightMatrix) parent;
            int oldLength = wmp.columns();
            while (true) {
            	switch (MathsTools.randomInt(4)) {
            	case 0:
            		if (oldLength < maxLength) {
            			return new PenalizedVariate(
            					parent,
            					new NMWeightMatrix(wmp.rawColumnArray(), oldLength + 1, wmp.offset() - 1),
            					0.0,
            					false,
            					this
            			);
            		}
            		break;
            	case 1:
            		if (oldLength < maxLength) {
            			return new PenalizedVariate(
            					parent,
            					new NMWeightMatrix(wmp.rawColumnArray(), oldLength + 1, wmp.offset()),
            					0.0,
            					false,
            					this
            			);
            		}
            		break;
            	case 2:
            		if (oldLength > minLength) {
            			return new PenalizedVariate(
            					parent,
            					new NMWeightMatrix(wmp.rawColumnArray(), oldLength - 1, wmp.offset() + 1),
            					0.0,
            					false,
            					this
            			);
            		}
            		break;
            	case 3:
            		if (oldLength > minLength) {
            			return new PenalizedVariate(
            					parent,
            					new NMWeightMatrix(wmp.rawColumnArray(), oldLength - 1, wmp.offset()),
            					0.0,
            					false,
            					this
            			);
            		}
            		break;
            	default:
            		throw new BioError("Out-of-bounds in RetrimSampler");
            	}
            }
        } catch (Exception ex) {
            throw new BioError("Assertion failed: error sampling preference");
        }
    }
    
    public String toString() {
        return "RetrimSampler()";
    }
}
