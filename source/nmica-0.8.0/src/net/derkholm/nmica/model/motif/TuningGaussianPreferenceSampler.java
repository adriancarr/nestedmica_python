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

import net.derkholm.nmica.maths.Gaussian;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.PenalizedVariate;
import net.derkholm.nmica.seq.consensus.ConsensusDistribution;
import org.biojava.bio.BioError;
import org.biojava.bio.dist.Count;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Sampler which adjusts base preferences while taking into account other
 * members of the ensemble
 *
 * @author Thomas Down
 */

public class TuningGaussianPreferenceSampler implements ContributionSampler, Serializable {
    public TuningGaussianPreferenceSampler() {
    	throw new RuntimeException("Not updated for NMWeightMatrix");
    }
    
    private double dotProduct(Count c1, Count c2)
        throws Exception
    {
        double sp = 0;
        for (Iterator i = ((FiniteAlphabet) c1.getAlphabet()).iterator(); i.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) i.next();
            sp += c1.getCount(s) * c2.getCount(s);
        }
        return sp;
    }
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles)
    {   
        try {
            WeightMatrix wmp = (WeightMatrix) parent;
            Distribution[] dists = new Distribution[wmp.columns()];
            int vCol = (int) Math.floor(Math.random() * wmp.columns());
            
            ConsensusDistribution oldCD = (ConsensusDistribution) wmp.getColumn(vCol);
            List<Double> angles = new ArrayList<Double>();
            {
                boolean[] dirt = new boolean[uncles.size()];
                int cnt = 0;
                while (cnt < 20) {
                    int index = (int) Math.floor(Math.random() * dirt.length);
                    if (dirt[index]) {
                        continue;
                    }
                    dirt[index] = true;
                    double angle = Math.acos(
                        dotProduct(
                            oldCD.getDirection(),
                            ((ConsensusDistribution) ((WeightMatrix) uncles.get(index)).getColumn(vCol)).getDirection()
                        )
                    );
                    if (Math.abs(angle) < 0.001) {
                        continue;
                    }
                    
                    angles.add(new Double(angle));
                    ++cnt;
                }
            }
            Collections.sort(angles);
            angles = angles.subList(0, 10);
            
            double precision;
            {
                double var = 0;
                for (Double a: angles) {
                    //var += Math.pow(a, 2.0);
	                var += a*a; // mrp: performance fix
                }
                precision = 4.0 / (var / angles.size());
            }
            
            if (Math.random() < 0.01) {
                System.err.println("Tuned to: " + precision);
            }
            
            for (int c = 0; c < wmp.columns(); ++c) {
                if (c == vCol) {
                    ConsensusDistribution cd = (ConsensusDistribution) wmp.getColumn(c);
                    dists[c] = new ConsensusDistribution(sampleDirection(cd.getDirection(), precision), cd.getHardness());
                } else {
                    dists[c] = wmp.getColumn(c);
                }
            }
            return new PenalizedVariate(parent, new SimpleWeightMatrix(dists), 0.0, false, this);
        } catch (Exception ex) {
            throw new BioError("Assertion failed: error sampling preference");
        }
    }

    private Count sampleDirection(Count c, double precision) 
        throws Exception
    {
        double xA = c.getCount(DNATools.a());
        double xC = c.getCount(DNATools.c());
        double xG = c.getCount(DNATools.g());
        double xT = c.getCount(DNATools.t());
        
        double a1 = Math.atan2(xC, xA);
        double a3 = Math.acos(xT);
        double a2 = Math.acos(xG / Math.sin(a3));
        
        double a1p = a1 + Gaussian.standardVariate() / precision;
        double a2p = a2 + Gaussian.standardVariate() / precision;
        double a3p = a3 + Gaussian.standardVariate() / precision;
        
        Count nc = new IndexedCount(DNATools.getDNA());
        nc.setCount(DNATools.a(), Math.sin(a3p) * Math.sin(a2p) * Math.cos(a1p));
        nc.setCount(DNATools.c(), Math.sin(a3p) * Math.sin(a2p) * Math.sin(a1p));
        nc.setCount(DNATools.g(), Math.sin(a3p) * Math.cos(a2p));
        nc.setCount(DNATools.t(), Math.cos(a3p));
        return nc;
    }
}
