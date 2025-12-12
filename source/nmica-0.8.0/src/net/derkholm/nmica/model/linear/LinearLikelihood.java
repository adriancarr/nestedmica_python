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

package net.derkholm.nmica.model.linear;

import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;

/**
 * Linear model of numerical (vector) data.
 *
 * @author Thomas Down
 */

public class LinearLikelihood implements LikelihoodCalculator {
    private final LinearFacette facette;
    private final Matrix1D data;
    
    LinearLikelihood(LinearFacette facette, Matrix1D data) {
        this.facette = facette;
        this.data = data;
    }
    
    public Facette getFacette() {
        return facette;
    }
    
    public Object getData() {
        return data;
    }
    
    public double likelihood(ObjectMatrix1D contributions, Matrix1D mix) {
        double L = 0;
        double noise = facette.getNoisePrecision();
        
        for (int t = 0; t < data.size(); ++t) {
            double x = 0;
            for (int c = 0; c < contributions.size(); ++c) {
                x += ((Matrix1D) contributions.get(c)).get(t) * mix.get(c);
            }
            L += gauss(data.get(t), x, noise);
        }
        
        return L;
    }
    
 	public static double gauss(double x, double mean, double precision) {
	    return Math.log(Math.sqrt(precision / 2 / Math.PI) * Math.exp(-precision / 2 * Math.pow(x - mean, 2)));
	}
}

