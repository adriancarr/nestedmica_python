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

/**
 * Motif sampler which randomly replaces columns with alternatives borrowed from an
 * uncle.
 *
 * @author Thomas Down
 */

public class CrossColumnSampler implements ContributionSampler, Serializable {
	private static final long serialVersionUID = 1L;

	public CrossColumnSampler() {}
    
    public PenalizedVariate sample(Object parent, ObjectMatrix1D uncles)
    {
        try {
            NMWeightMatrix wmp = (NMWeightMatrix) parent;
            NMWeightMatrix wmu = (NMWeightMatrix) uncles.get((int) Math.floor(Math.random() * uncles.size()));
            Distribution[] columns = new Distribution[wmp.unmaskedColumns()];
            
            boolean vc = false;   // track visible changes.
            for (int c = 0; c < columns.length; ++c) {
                if (Math.random() < 0.8) {
                    columns[c] = wmp.getUnmaskedColumn(c);
                } else {
                    columns[c] = wmu.getUnmaskedColumn(c);
                    if (!wmp.isHidden(c)) {
                    	vc = true;
                    }
                }
            }
            
            return new PenalizedVariate(parent, new NMWeightMatrix(columns, wmp.columns(), wmp.offset()), 0.0, !vc, this);
        } catch (Exception ex) {
            throw new BioError("Assertion failed: error sampling preference");
        }
    }
}
