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

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;

import org.biojava.bio.symbol.SymbolList;

/**
 * Calculate likelihood functions for short sequence motifs embedded in a background,
 * assuming some degree of clustering of motif instances.
 * 
 * <p>
 * <strong>WARNING:</strong> MotifUncountedLikelihood instances are not threadsafe.
 * </p>
 *
 * @author Thomas Down
 */

class MotifUncountedClusterLikelihood extends MotifUncountedLikelihood {
    MotifUncountedClusterLikelihood(MotifFacette mmf, SymbolList seq) 
        throws Exception
    {
        super(mmf, seq);
    }
    
    protected double sumMotifs(double[] bgScores, int numMotifs, Matrix2D wmScores, int[] advances, double[] motifTrans)
        throws Exception
    {
        int length = bgScores.length;
        Matrix2D matrix = new SimpleMatrix2D(2, length + 1);
        matrix.set(0, 0, NativeMath.log2(0.5));
        matrix.set(1, 0, NativeMath.log2(0.5));
        
        MotifFacette facette = getMotifFacette();
        
        double sumTrans = 0;
        double[] motifPenalty = new double[numMotifs];
        for (int m = 0; m < numMotifs; ++m) {
            motifPenalty[m] = NativeMath.log2(motifTrans[m]);
            sumTrans += motifTrans[m];
        }
        double clusterInPenalty = NativeMath.log2(facette.getClusterIn());
        double clusterOutPenalty = NativeMath.log2(facette.getClusterOut());
        double outsideBasePenalty = NativeMath.log2(1.0 - facette.getClusterIn());
        double insideBasePenalty = NativeMath.log2(1.0 - sumTrans - facette.getClusterOut());
        
        for (int i = 1; i <= length; ++i) {
        	double bgEmit = bgScores[i - 1];
        	{
        		// inside
	            double score = matrix.get(0, i - 1) + bgEmit + insideBasePenalty;
	            score = NativeMath.addLog2(score, matrix.get(1, i - 1) + bgEmit + clusterInPenalty);
	            for (int m = 0; m < numMotifs; ++m) {
	                int wml = advances[m];
	                if (i >= wml) {
	                    double fromScore = matrix.get(0, i - wml);
	                    double emitScore =  wmScores.get(i - 1, m);
	                    score = NativeMath.addLog2(score, fromScore + emitScore + motifPenalty[m]);
	                }
	            }
	            matrix.set(0, i, score);
        	}
        	{
        		// outside
        		double score = matrix.get(1, i - 1) + bgEmit + outsideBasePenalty;
        		score = NativeMath.addLog2(score, matrix.get(0, i - 1) + bgEmit + clusterOutPenalty);
        		matrix.set(1, i, score);
        	}
        }
        
        return NativeMath.addLog2(matrix.get(0, length), matrix.get(1, length));
    }
}
