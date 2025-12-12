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

/**
 * @author thomas
 */
public class BinaryCommitableMatrix2D implements CommitableMatrix2D {
    private static final int EDIT_LIST_SIZE = 5;
    
    private final int rows;
    private final int columns;
    private byte[] foreground;
    private byte[] background;
    
    private int editCount;
    private int[] editIndex;
    
    public BinaryCommitableMatrix2D(int rows, int columns) {
        this.rows = rows;
        this.columns = columns;
        int cells = rows * columns;
        this.foreground = new byte[cells];
        this.background = new byte[cells];
        editIndex = new int[EDIT_LIST_SIZE];
        editCount = 0;
    }
    
    public BinaryMatrixSlice sliceRow(int r) {
    	return new BinaryMatrixSlice(foreground, r * columns, columns);
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
        return foreground[row * columns + col];
    }

    public double getCommitted(int row, int col) {
        return background[row * columns + col];
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#set(int, int, double)
     */
    public void set(int row, int col, double d) {
        int cell = row * columns + col;
        foreground[cell] = (d == 0.0 ? (byte) 0 : (byte) 1);
        if (editCount < EDIT_LIST_SIZE) {
            editIndex[editCount] = cell;
        }
        ++editCount;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#getRaw()
     */
    public double[] getRaw() {
        throw new UnsupportedOperationException();
    }
    
    public void commit() {
        sync(background, foreground);
    }
    
    public void rollback() {
        sync(foreground, background);
    }
    
    public boolean isDirty() {
        return editCount > 0;
    }
    
    private void sync(byte[] to, byte[] from) {
        if (editCount <= EDIT_LIST_SIZE) {
            for (int e = 0; e < editCount; ++e) {
                int i = editIndex[e];
                to[i] = from[i];
            }
        } else {
            System.arraycopy(from, 0, to, 0, to.length);
        }
        editCount = 0;
    }

}
