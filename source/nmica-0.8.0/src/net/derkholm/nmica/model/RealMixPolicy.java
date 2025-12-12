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

import net.derkholm.nmica.maths.Cauchy;
import net.derkholm.nmica.matrix.CommitableMatrix2D;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleCommitableMatrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;

/**
 * Prior and sampler for mixing matrices
 *
 * @author Thomas Down
 */
 
public class RealMixPolicy implements MixPolicy, Serializable {
    private static final long serialVersionUID = 301925405194405767L;
    
    final double minValue;
    final double maxValue;
    double stepSize = 0.1;
    
    public RealMixPolicy(double min, double max) {
        this.minValue = min;
        this.maxValue = max;
    }
    
    public RealMixPolicy(double min, double max, double stepSize) {
        this.minValue = min;
        this.maxValue = max;
        this.stepSize = stepSize;
    }
    
    public void variate(Matrix1D target) {
        for (int i = 0; i < target.size(); ++i) {
            target.set(i, minValue + Math.random() * (maxValue - minValue));
        }
    }
    
    public double prior(Matrix1D m) {
        for (int i = 0; i < m.size(); ++i) {
            double v = m.get(i);
            if (v < minValue || v > maxValue) {
                return Double.NEGATIVE_INFINITY;
            } 
        }
        return 0;
    }
    
    public void sample(Matrix1D target) {
        int flipPos = (int) Math.floor(Math.random() * target.size());
        target.set(flipPos, target.get(flipPos) + 0.1 * Cauchy.standardVariate());
    }
    
    public void sampleComponent(Matrix1D target, int c) {
        target.set(c, target.get(c) + 0.1 * Cauchy.standardVariate());
    }
    
    public boolean isBinary() {
        return false;
    }
    
    public Matrix2D createCompatibleMatrix(int rows, int columns) {
        return new SimpleMatrix2D(rows, columns);
    }
    
    public CommitableMatrix2D createCommitableMatrix(int rows, int columns) {
        return new SimpleCommitableMatrix2D(rows, columns);
    }
}
