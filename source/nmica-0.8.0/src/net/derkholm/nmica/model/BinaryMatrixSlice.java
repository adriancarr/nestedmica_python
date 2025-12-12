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

import net.derkholm.nmica.matrix.Matrix1D;

public class BinaryMatrixSlice implements Matrix1D {
	private final byte[] data;
	private final int offset;
	private final int size;
	
	BinaryMatrixSlice(byte[] data, int offset, int size) {
		this.data = data;
		this.offset = offset;
		this.size = size;
	}

	public double get(int pos) {
		return data[offset + pos];
	}

	public double[] getRaw() {
		throw new UnsupportedOperationException();
	}

	public void set(int pos, double v) {
		throw new UnsupportedOperationException();
	}

	public int size() {
		return size;
	}
	
	public byte[] getByteData() {
		return data;
	}
	
	public int getByteOffset() {
		return offset;
	}
}
