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

package net.derkholm.nmica.model.logitseq;

import java.io.Serializable;

import org.biojava.bio.dp.WeightMatrix;

import net.derkholm.nmica.maths.Gaussian;
import net.derkholm.nmica.model.ContributionPrior;
import net.derkholm.nmica.model.motif.WeightMatrixPrior;

/**
 * @author thomas
 */
public class WeightedWeightMatrixPrior implements ContributionPrior, Serializable {
    private static final long serialVersionUID = 2788635461542022614L;
    
    private final WeightMatrixPrior wmPrior;
    private final double precision = 0.1;
    
    public WeightedWeightMatrixPrior(WeightMatrixPrior wmPrior) {
        super();
        this.wmPrior = wmPrior;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionPrior#getContributionType()
     */
    public Class getContributionType() {
        return WeightedWeightMatrix.class;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionPrior#probability(java.lang.Object)
     */
    public double probability(Object o) throws IllegalArgumentException {
        WeightedWeightMatrix wwm = (WeightedWeightMatrix) o;
        return wmPrior.probability(wwm.getWeightMatrix()) + Gaussian.probability(wwm.getWeight(), 0, precision);
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionPrior#variate()
     */
    public Object variate() {
        return new WeightedWeightMatrix(
               (WeightMatrix) wmPrior.variate(),
               Gaussian.standardVariate() / precision
        );
    }

}
