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

import net.derkholm.nmica.matrix.CommitableMatrix2D;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;


/**
 * Mix policy which imposes a constant mixing matrix
 * 
 * @author thomas
 */
public class ConstantMixPolicy implements MixPolicy {

	private final CommitableMatrix2D matrix;
	
	public ConstantMixPolicy(Matrix2D m)
	{
		final Matrix2D cm = new SimpleMatrix2D(m);
		matrix = new CommitableMatrix2D() {
			public double getCommitted(int row, int col) {
				return cm.get(row, col);
			}

			public int columns() {
				return cm.columns();
			}

			public double get(int row, int col) {
				return cm.get(row, col);
			}

			public double[] getRaw() {
				throw new UnsupportedOperationException();
			}

			public int rows() {
				return cm.rows();
			}

			public void set(int row, int col, double v) {
				if (v != get(row, col)) {
					throw new UnsupportedOperationException();
				}
				// otherwise succeed silently.
			}

			public void commit() {
			}

			public boolean isDirty() {
				return false;
			}

			public void rollback() {
			}
		};
	}
	
	public CommitableMatrix2D createCommitableMatrix(int rows, int columns) {
		if (rows != matrix.rows() || columns != matrix.columns()) {
			throw new IllegalArgumentException("Matrix sizes don't match");
		}
		return matrix;
	}

	public Matrix2D createCompatibleMatrix(int rows, int columns) {
		return createCommitableMatrix(rows, columns);
	}

	public boolean isBinary() {
		return true;
	}

	public double prior(Matrix1D m) {
		return 0;
	}

	public void sample(Matrix1D target) {
	}

	public void sampleComponent(Matrix1D target, int c) {
	}

	public void variate(Matrix1D target) {
	}

}
