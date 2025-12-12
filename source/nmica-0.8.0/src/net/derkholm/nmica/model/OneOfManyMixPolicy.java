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

package net.derkholm.nmica.model;

import java.io.Serializable;
import net.derkholm.nmica.matrix.*;

/**
 * Prior and sampler for mixing matrices
 *
 * @author Thomas Down
 */
 
public class OneOfManyMixPolicy implements MixPolicy, Serializable {
	private static final long serialVersionUID = 6265414984872659046L;

	public OneOfManyMixPolicy() {
    }
    
    public void variate(Matrix1D target) {
        int tv = (int) Math.floor(Math.random() * target.size());
        for (int i = 0; i < target.size(); ++i) {
            target.set(i, i == tv ? 1.0 : 0.0);
        }
    }
    
    public double prior(Matrix1D m) {
        int numNonZero = 0;
        for (int i = 0; i < m.size(); ++i) {
            if (m.get(i) != 0.0) {
                ++numNonZero;
            }
        }
        return numNonZero == 1 ? 0.0 : Double.NEGATIVE_INFINITY;
    }
    
    public void sample(Matrix1D target) {
        if (target.size() == 1) {
            return;
        }
        int cnz = -1;
        for (int i = 0; i < target.size(); ++i) {
            if (target.get(i) != 0.0) {
                cnz = i;
            }
        }
        if (cnz < 0) {
            cnz = 0;
        }
        int nnz = cnz;
        while (nnz == cnz) {
            nnz = (int) Math.floor(Math.random() * target.size());
        }
        target.set(cnz, 0.0);
        target.set(nnz, 1.0);
    }
    
    public void sampleComponent(Matrix1D target, int c) {
        // not good.
    }
    
    public boolean isBinary() {
        return true;
    }
    
    public Matrix2D createCompatibleMatrix(int rows, int columns) {
        return new SimpleMatrix2D(rows, columns);
    }
    
    public CommitableMatrix2D createCommitableMatrix(int rows, int columns) {
        return new SimpleCommitableMatrix2D(rows, columns);
    }
}
