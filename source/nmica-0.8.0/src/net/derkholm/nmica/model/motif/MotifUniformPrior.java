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

import net.derkholm.nmica.seq.consensus.ConsensusDistribution;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Count;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.seq.DNATools;

/**
 * Simple non-informative prior over motifs.  (Currently only does DNA
 * motifs, but could probabably be generalized)
 *
 * @author Thomas Down
 */

public class MotifUniformPrior extends WeightMatrixPrior {
    public MotifUniformPrior(int minLength, int maxLength, int unmaskedLength) {
        super(DNATools.getDNA(), minLength, maxLength, unmaskedLength);
    }
    
    public double probabilityColumn(Distribution d) 
    {
        ConsensusDistribution cd = (ConsensusDistribution) d;
        double hard = cd.getHardness();
        if (hard < 0 || hard > 10.0) {
            return Double.NEGATIVE_INFINITY;
        }
        return 0.0;
    }
    
    public Distribution variateColumn() {
        try {
            return new ConsensusDistribution(randomUnitVector(), 10.0 * Math.random());
        } catch (Exception ex) {
            throw new BioError("Assertion failed: couldn't build consensus", ex);
        }
    }
    
    private Count randomUnitVector() {
        try {
            double t0 = 2.0 * Math.PI * Math.random();
            double t1 = 2.0 * Math.PI * Math.random();
            double t2 = 2.0 * Math.PI * Math.random();
            
            Count dir = new IndexedCount(DNATools.getDNA());
            dir.setCount(DNATools.a(), Math.sin(t0) * Math.sin(t1) * Math.cos(t2));
            dir.setCount(DNATools.c(), Math.sin(t0) * Math.sin(t1) * Math.sin(t2));
            dir.setCount(DNATools.g(), Math.sin(t0) * Math.cos(t1));
            dir.setCount(DNATools.t(), Math.cos(t0));
            return dir;
        } catch (Exception ex) {
            throw new BioError("Assertion failed: couldn't build consensus", ex);
        }
    }
}
