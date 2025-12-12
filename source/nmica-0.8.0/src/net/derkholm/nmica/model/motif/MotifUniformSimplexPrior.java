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

import java.util.Arrays;
import java.util.Iterator;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.seq.NMSimpleDistribution;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 * Simple non-informative prior over motifs.  (Currently only does DNA
 * motifs, but could probabably be generalized)
 *
 * @author Thomas Down
 */

public class MotifUniformSimplexPrior extends WeightMatrixPrior {
	private static final long serialVersionUID = 1L;
	
	private AtomicSymbol[] symArray;
    
    public MotifUniformSimplexPrior(FiniteAlphabet alpha, int minLength, int maxLength, int unmaskedLength) {
        super(alpha, minLength, maxLength, unmaskedLength);
        
        symArray = new AtomicSymbol[alpha.size()];
        {
            int p = 0;
            for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
                symArray[p++] = (AtomicSymbol) si.next();
            }
        }
        
    }
    
    public double probabilityColumn(Distribution d) 
    {
        return 0.0;
    }
    
    public Distribution variateColumn() {
        while (true) {
            Distribution d = rawVariateColumn();
            double p = probabilityColumn(d);
            if (p == Double.NEGATIVE_INFINITY || (p < 0.0 && Math.random() > NativeMath.exp2(p))) {
                // resample
            } else {
                return d;
            }
        }
    }
    
    public Distribution rawVariateColumn() {
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
