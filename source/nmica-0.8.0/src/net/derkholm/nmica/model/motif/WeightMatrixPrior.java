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
import net.derkholm.nmica.model.ContributionPrior;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;

/**
 * Simple non-informative prior over motifs. 
 *
 * @author Thomas Down
 */

public abstract class WeightMatrixPrior implements ContributionPrior, Serializable {
	private static final long serialVersionUID = 1000002L;
	
    private final int minLength;
    private final int maxLength;
    private final int unmaskedLength;
    private final Alphabet alpha;
    
    protected WeightMatrixPrior(Alphabet alpha, int minLength, int maxLength, int unmaskedLength) {
        this.alpha = alpha;
        this.minLength = minLength;
        this.maxLength = maxLength;
        this.unmaskedLength = unmaskedLength;
    }
    
    public final int getMinLength() {
    	return minLength;
    }
    
    public final int getMaxLength() {
    	return maxLength;
    }
    
    public final int getUnmaskedLength() {
    	return unmaskedLength;
    }
    
    public final Alphabet getAlphabet() {
        return alpha;
    }
    
    public final Class getContributionType() {
        return WeightMatrix.class;
    }
    
    public final double probability(Object o) 
        throws IllegalArgumentException
    {
        NMWeightMatrix wm = (NMWeightMatrix) o;
        int cc = wm.columns(), allcc = wm.unmaskedColumns();
        if (cc < minLength || cc > maxLength) {
            return Double.NEGATIVE_INFINITY;
        }
        double P = 0;
        for (int c = 0; c < allcc; ++c) {
            try {
                P += probabilityColumn(wm.getUnmaskedColumn(c));
            } catch (IllegalAlphabetException ex) {
                throw new IllegalArgumentException("Couldn't handle column");
            }
        }
        return P;
    }
    
    public final Object variate() {
        Distribution[] dists = new Distribution[unmaskedLength];
        int length = minLength + MathsTools.randomInt(maxLength - minLength + 1);
        int offset = MathsTools.randomInt(unmaskedLength);
        for (int c = 0; c < dists.length; ++c) {
            dists[c] = variateColumn();
        }
        try {
            return new NMWeightMatrix(dists, length, offset);
        } catch (IllegalAlphabetException ex) {
            throw new BioError("Assertion failed: couldn't build motif", ex);
        }
    }
    
    public abstract double probabilityColumn(Distribution dist) throws IllegalAlphabetException;
    
    public abstract Distribution variateColumn();
    
}
