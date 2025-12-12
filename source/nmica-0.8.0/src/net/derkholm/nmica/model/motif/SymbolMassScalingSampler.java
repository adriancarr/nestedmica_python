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
import java.util.Iterator;

import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.PenalizedVariate;
import net.derkholm.nmica.seq.NMSimpleDistribution;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 * Motif sampler which adjusts the symbol preference vector by sampling from a
 * gaussian on the surface of the preference hypersphere.
 *
 * @author Thomas Down
 */

public class SymbolMassScalingSampler implements ContributionSampler, Serializable {
    private final static long serialVersionUID = -3756360057487333269L;
    
    private final double maxStep;
    
    public SymbolMassScalingSampler(double maxStep) {
        this.maxStep = maxStep;
    }
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles)
    {
        double step = Math.random() * maxStep;
        
        try {
            NMWeightMatrix wmp = (NMWeightMatrix) parent;
            Distribution[] dists = new Distribution[wmp.unmaskedColumns()];
            int vCol = (int) (Math.random() * wmp.unmaskedColumns());
            for (int c = 0; c < wmp.unmaskedColumns(); ++c) {
                if (c == vCol) {
                    Distribution old = wmp.getUnmaskedColumn(c);
                    try {
                        Distribution newDist;
                        if (Math.random() < 0.5) {
                            newDist = sampleDist(old, step);
                        } else {
                            newDist = sampleDistBack(old, step);
                        }
                        dists[c] = newDist;
                    } catch (Exception ex) {
                        throw new BioError(ex);
                    }
                } else {
                    dists[c] = wmp.getUnmaskedColumn(c);
                }
            }
            return new PenalizedVariate(
            		parent, 
            		new NMWeightMatrix(dists, wmp.columns(), wmp.offset()),
            		0.0, 
            		wmp.isHidden(vCol),
            		this
            );
        } catch (Exception ex) {
            throw new BioError("Assertion failed: error sampling preference");
        }
    }
    
    private Distribution sampleDist(Distribution dist, double step) 
        throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) dist.getAlphabet();
        Distribution nudist = new NMSimpleDistribution(alpha);
        AtomicSymbol beneficiary = (AtomicSymbol) new UniformDistribution(alpha).sampleSymbol();
        double inc = dist.getWeight(beneficiary) * step;
        for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) si.next();
            double p = dist.getWeight(s);
            double w;
            if (s == beneficiary) {
                w = (p + inc) / (1.0 + inc);
            } else {
                w = p / (1.0 + inc);
            }
            nudist.setWeight(s, w);
        }
        return nudist;
    }
    
    private Distribution sampleDistBack(Distribution dist, double step) 
        throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) dist.getAlphabet();
        Distribution nudist = new NMSimpleDistribution(alpha);
        AtomicSymbol beneficiary = (AtomicSymbol) new UniformDistribution(alpha).sampleSymbol();
        
        double bFwd = dist.getWeight(beneficiary);
        double bBack = bFwd / (1.0 + step * (1 - bFwd));
        for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) si.next();
            double p = dist.getWeight(s);
            double w;
            if (s == beneficiary) {
                w = bBack;
            } else {
                w = p * (1.0 - bBack) / (1.0 - bFwd);
            }
            nudist.setWeight(s, w);
        }
        return nudist;
    }
    
    public String toString() {
        return String.format("SymbolMassScalingSampler(maxStep=%g)", maxStep);
    }
}
