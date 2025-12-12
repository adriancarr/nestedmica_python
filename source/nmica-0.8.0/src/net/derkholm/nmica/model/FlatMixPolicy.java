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
 * Mixing policy which leaves the mixing matrix fixed with all elements equal
 * to 1.0
 *
 * @author Thomas Down
 */
 
public class FlatMixPolicy implements MixPolicy, Serializable {
	private static final long serialVersionUID = 1L;

	public FlatMixPolicy() {
    }
    
    public void variate(Matrix1D target) {
        for (int i = 0; i < target.size(); ++i) {
            target.set(i, 1.0);
        }
    }
    
    public double prior(Matrix1D m) {
        for (int i = 0; i < m.size(); ++i) {
            if (m.get(i) != 0.0) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0;
    }
    
    public void sample(Matrix1D target) {
    }
    
    public void sampleComponent(Matrix1D target, int c) {
    }
    
    public Matrix2D createCompatibleMatrix(int rows, int columns) {
        return new FlatMixingMatrix(rows, columns);
    }
    
    public CommitableMatrix2D createCommitableMatrix(int rows, int columns) {
        return new FlatMixingMatrix(rows, columns);
    }
    
    public boolean isBinary() {
        return true;
    }
}
