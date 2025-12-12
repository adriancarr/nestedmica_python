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

import net.derkholm.nmica.matrix.CommitableMatrix2D;

/**
 * Dummy Commitable implementation of Matrix2D where all elmenents are fixed at 1.0
 * 
 * @author thomas
 */
class FlatMixingMatrix implements CommitableMatrix2D, Serializable {
	private static final long serialVersionUID = -719377647744578596L;
	
	private final int rows;
    private final int columns;
    
    public FlatMixingMatrix(int rows, int columns) {
        this.rows = rows;
        this.columns = columns;
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
        return 1.0;
    }
    
    public double getCommitted(int row, int col) {
        return 1.0;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#set(int, int, double)
     */
    public void set(int row, int col, double v) {
        if (v != 1.0) {
            throw new UnsupportedOperationException();
        }
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.matrix.Matrix2D#getRaw()
     */
    public double[] getRaw() {
        throw new UnsupportedOperationException();
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.Commitable#commit()
     */
    public void commit() {
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.Commitable#rollback()
     */
    public void rollback() {
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.utils.Commitable#isDirty()
     */
    public boolean isDirty() {
        return false;
    }

}
