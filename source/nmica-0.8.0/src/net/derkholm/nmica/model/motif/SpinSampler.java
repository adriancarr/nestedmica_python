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

import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.PenalizedVariate;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;

/**
 * Motif sampler which 'spins' the motif by one column left or right, replacing
 * the lost column with a random sample.
 *
 * @author Thomas Down
 */

public class SpinSampler implements ContributionSampler, Serializable {
    private final static long serialVersionUID = -7946579690595955650L;
    
    private final WeightMatrixPrior prior;
    
    public SpinSampler(WeightMatrixPrior prior) {
        this.prior = prior;
        throw new RuntimeException("Not updated for NMWeightMatrix");
    }
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles)
    {
        try {
            WeightMatrix wmp = (WeightMatrix) parent;
            Distribution[] columns = new Distribution[wmp.columns()];
            
            boolean newStart = (Math.random() < 0.5);
            Distribution newDist = prior.variateColumn();
            if (newStart) {
                columns[0] = newDist;
                for (int i = 0; i < (wmp.columns() - 1); ++i) {
                    columns[i + 1] = wmp.getColumn(i);
                }
            } else {
                columns[wmp.columns() - 1] = newDist;
                for (int i = 1; i < wmp.columns(); ++i) {
                    columns[i - 1] = wmp.getColumn(i);
                }
            }
            
            return new PenalizedVariate(parent, new SimpleWeightMatrix(columns), 0.0, false, this);
        } catch (Exception ex) {
            throw new BioError("Assertion failed: error sampling preference");
        }
    }
    
    public String toString() {
        return "SpinSampler()";
    }
}
