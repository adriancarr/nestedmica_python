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
import java.util.Arrays;
import java.util.Iterator;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.model.ContributionPrior;
import net.derkholm.nmica.seq.NMSimpleDistribution;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * Simple non-informative prior over variable length motifs in a particular
 * length range.  
 * (Currently only does DNA motifs, but could probabably be generalized)
 *
 * @author Thomas Down
 */

public class MotifDropoffPrior implements ContributionPrior, Serializable {
    private static final long serialVersionUID = 7211783954582017859L;
    
    private final AtomicSymbol[] symArray;
    private final FiniteAlphabet alpha;
    private final int length;
    private final int maxDropoff;
    private final double dropoffThreshold = 0.3;
    
    public MotifDropoffPrior(FiniteAlphabet alpha, int length, int maxDropoff) {
        super();
        this.alpha = alpha;
        this.length = length;
        this.maxDropoff = maxDropoff;
        
        symArray = new AtomicSymbol[alpha.size()];
        {
            int p = 0;
            for (Iterator si = alpha.iterator(); si.hasNext(); ) {
                symArray[p++] = (AtomicSymbol) si.next();
            }
        }
    }
    
    public final int getLength() {
        return length;
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
        WeightMatrix wm = (WeightMatrix) o;
        if (wm.columns() != length) {
            return Double.NEGATIVE_INFINITY;
        }
        double P = 0;
        for (int c = 0; c < wm.columns(); ++c) {
            try {
                int drops = Math.max(0, Math.max(maxDropoff - c - 1, maxDropoff - (length - c - 1) - 1));
                if (drops > 0) {
                    Distribution dist = wm.getColumn(c);
                    double maxWeight = 0;
                    for (int x = 0; x < symArray.length; ++x) {
                        maxWeight = Math.max(maxWeight, dist.getWeight(symArray[x]));
                    }
                    if (maxWeight > dropoffThreshold) {
                        P += NativeMath.log2((1.0 * (maxDropoff - drops)) / maxDropoff);
                    }
                }
            } catch (IllegalSymbolException ex) {
                throw new IllegalArgumentException("Couldn't handle column");
            }
        }
        return P;
    }
    
    public final Object variate() {
        int trimLeft = (int) Math.floor(Math.random() * maxDropoff);
        int trimRight = (int) Math.floor(Math.random() * maxDropoff);
        Distribution[] dists = new Distribution[length];
        for (int c = 0; c < dists.length; ++c) {
            if (c < trimLeft || (length - c - 1) < trimRight) {
                dists[c] = droppedVariateColumn();
            } else {
                dists[c] = variateColumn();
            }
        }
        try {
            return new SimpleWeightMatrix(dists);
        } catch (IllegalAlphabetException ex) {
            throw new BioError("Assertion failed: couldn't build motif", ex);
        }
    }
    
    private Distribution droppedVariateColumn() {
        double uniform = 1.0 / symArray.length;
        double variableContribution = dropoffThreshold - uniform;
        double uniformContribution = 1.0 - variableContribution;
        try {
            Distribution d = variateColumn();
            for (int x = 0; x < symArray.length; ++x) {
                d.setWeight(
                        symArray[x],
                        uniformContribution * uniform +
                        variableContribution * d.getWeight(symArray[x])
                );
            }
            return d;
        } catch (Exception ex) {
            throw new BioError("Assertion failed: coudln't modify column", ex);
        }
    }
    
    private Distribution variateColumn() {
        double[] cutArray = new double[symArray.length];
        cutArray[0] = 1.0;
        for (int x = 1; x < cutArray.length; ++x) {
            cutArray[x] = Math.random();
        }
        Arrays.sort(cutArray);
        
        try {
            Distribution d = new NMSimpleDistribution((FiniteAlphabet) getAlphabet());
            double last = 0;
            for (int x = 0; x < cutArray.length; ++x) {
                d.setWeight(symArray[x], cutArray[x] - last);
                last = cutArray[x];
            }
            return d;
        } catch (Exception ex) {
            throw new BioError(ex);
        }
    }
}
