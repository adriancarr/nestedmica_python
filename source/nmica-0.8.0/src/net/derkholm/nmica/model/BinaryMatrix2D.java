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

import net.derkholm.nmica.matrix.Matrix2D;

/**
 * @author thomas
 */
class BinaryMatrix2D implements Matrix2D, Serializable {
	private static final long serialVersionUID = 29187494607829404L;
	
	private final int rows;
    private final int columns;
    private byte[] data;
    
    public BinaryMatrix2D(int rows, int columns) {
        this.rows = rows;
        this.columns = columns;
        this.data = new byte[rows * columns];
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#rows()
     */
    public int rows() {
        return rows;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#columns()
     */
    public int columns() {
        return columns;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#get(int, int)
     */
    public double get(int row, int col) {
        return data[row * columns + col];
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#set(int, int, double)
     */
    public void set(int row, int col, double v) {
        data[row * columns + col] = (v == 0.0 ? (byte) 0 : (byte) 1);
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#getRaw()
     */
    public double[] getRaw() {
        throw new UnsupportedOperationException();
    }

}
