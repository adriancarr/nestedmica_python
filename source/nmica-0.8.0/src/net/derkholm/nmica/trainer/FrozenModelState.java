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

package net.derkholm.nmica.trainer;

import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.MatrixTools;
import net.derkholm.nmica.matrix.ObjectMatrix2D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix2D;
import net.derkholm.nmica.utils.ArrayTools;

/**
 * Opaque representation of a Frozen NSTrainer state.  Intended for checkpointing
 * during a training process
 *
 * @author Thomas Down
 */

class FrozenModelState implements java.io.Serializable {
    FrozenModelState(TrainableState ts) {
        contributions = new SimpleObjectMatrix2D(ts.getContributions());
        // mixingMatrix  = new SimpleMatrix2D(ts.getMixingMatrix());
        {
            Matrix2D tMatrix = ts.getMixingMatrix();
            mixingMatrix = ts.getContext().getMixPolicy().createCompatibleMatrix(tMatrix.rows(), tMatrix.columns());
            MatrixTools.copy(mixingMatrix, tMatrix);
        }
        permute = (int []) ArrayTools.copy(ts.getPermutation());
    }
    
    ObjectMatrix2D contributions;
    Matrix2D mixingMatrix;
    int[] permute;
}
