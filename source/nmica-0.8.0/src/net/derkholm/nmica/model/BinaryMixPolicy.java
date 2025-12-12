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

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.*;

/**
 * Prior and sampler for mixing matrices
 *
 * @author Thomas Down
 */
 
public class BinaryMixPolicy implements MixPolicy, Serializable {
    private static final long serialVersionUID = -3957640888327427114L;
    
    private final double positiveWeight;
    private final double posl;
    private final double negl;
    
    public BinaryMixPolicy(double positiveWeight) {
        this.positiveWeight = positiveWeight;
        posl = NativeMath.log2(positiveWeight);
        negl = NativeMath.log2(1.0 - positiveWeight);
    }
    
    public void variate(Matrix1D target) {
        for (int i = 0; i < target.size(); ++i) {
            target.set(i, Math.random() < positiveWeight ? 1.0 : 0.0);
        }
    }
    
    public double prior(Matrix1D m) {
        double L = 0;
	int mSize = m.size();
        for (int i = 0; i < mSize; ++i) {
            if (m.get(i) == 0.0) {
                L += negl;
            } else {
                L += posl;
            }
        }
        return L;
    }
    
    public void sample(Matrix1D target) {
        int flipPos = (int) Math.floor(Math.random() * target.size());
        if (target.get(flipPos) == 0.0) {
            target.set(flipPos, 1.0);
        } else {
            target.set(flipPos, 0.0);
        }
    }
    
    public void sampleComponent(Matrix1D target, int c) {
        target.set(c, Math.random() < positiveWeight ? 1.0 : 0.0);
    }
    
    public boolean isBinary() {
        return true;
    }
    
    public Matrix2D createCompatibleMatrix(int rows, int columns) {
        // return new SimpleMatrix2D(rows, columns);
        return new BinaryMatrix2D(rows, columns);
    }
    
    public CommitableMatrix2D createCommitableMatrix(int rows, int columns) {
        // return new SimpleCommitableMatrix2D(rows, columns);
        return new BinaryCommitableMatrix2D(rows, columns);
    }
}
