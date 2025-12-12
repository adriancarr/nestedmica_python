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

package net.derkholm.nmica.model.motif;

import net.derkholm.nmica.matrix.Matrix2D;
import org.biojava.bio.symbol.SymbolList;

/**
 * Version of MotifUncountedLikelihood which uses native code to implement
 * the dynamic programming step.
 * 
 * @author thomas
 */

class MotifUncountedLikelihoodNative extends MotifUncountedLikelihood {
    /**
     * @param mmf
     * @param seq
     * @throws Exception if there is an error with the precomputed parts of
     *         the likelihood calculation
     */
    MotifUncountedLikelihoodNative(MotifFacette mmf, SymbolList seq) throws Exception {
        super(mmf, seq);
    }

    static {
        System.loadLibrary("nmica");
    }
    
    protected double sumMotifs(
            double[] bgScores, 
            int numMotifs,
            Matrix2D wmScores, 
            int[] advances, 
            double[] motifTrans
    ) 
    {
        return nativeSumMotifs(bgScores.length, bgScores, numMotifs, wmScores.getRaw(), wmScores.rows(), wmScores.columns(), advances, motifTrans);
    }
    
    private native double nativeSumMotifs(
            int length,
            double[] bgScores, 
            int numMotifs,
            double[] wmScores, 
            int wmScores_rows, 
            int wmScores_columns, 
            int[] advances, 
            double[] motifTrans
    );
}
