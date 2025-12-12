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

package net.derkholm.nmica.model.logitseq;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionItem;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.utils.CollectTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.util.ArrayList;
import java.util.List;

/**
 * @author thomas
 */

class LogisticSequenceLikelihoodCalculator implements LikelihoodCalculator {
    private final LogisticSequenceFacette facette;
    private final Sequence sequence;
    private final int[] trimmedSeqIndices;
    private final int label;
    
    public LogisticSequenceLikelihoodCalculator(LogisticSequenceFacette facette, Sequence seq) 
    		throws IllegalSymbolException
    {
        this.facette = facette;
        this.sequence = seq;
        
        AlphabetIndex index = facette.getAlphabetIndex();
        
        List<Integer> indexList = new ArrayList<Integer>();
        // mrp: not used? List hoodList = new ArrayList();
        for (int pos = 1; pos <= seq.length(); ++pos) {
            Symbol s = seq.symbolAt(pos);
            if (s instanceof AtomicSymbol) {
                indexList.add(new Integer(index.indexForSymbol(s)));
            } else {
                indexList.add(new Integer(-1));
            }
        }
        
        trimmedSeqIndices = CollectTools.toIntArray(indexList);
        
        if (seq.getAnnotation().containsProperty("mocca.label")) {
            label = ((Integer) seq.getAnnotation().getProperty("mocca.label")).intValue();
        } else {
            throw new IllegalArgumentException("Unlabelled sequence");
        }
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.LikelihoodCalculator#getFacette()
     */
    public Facette getFacette() {
        return facette;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.LikelihoodCalculator#getData()
     */
    public Object getData() {
        return sequence;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.LikelihoodCalculator#likelihood(net.derkholm.nmica.matrix.ObjectMatrix1D, net.derkholm.nmica.matrix.Matrix1D)
     */
    public double likelihood(ObjectMatrix1D contributions, Matrix1D weights) {
        double eta = 0;
        for (int i = 0; i < contributions.size(); ++i) {
            ContributionItem ci = (ContributionItem) contributions.get(i);
            WeightedWeightMatrix wwm = (WeightedWeightMatrix) ci.getItem();
            Matrix2D bm = (Matrix2D) ci.getItemView(facette.forwardBmView());
            
            int wml = bm.columns();
            double odds = wml * NativeMath.fastlog2(0.25);
            int maxPos = trimmedSeqIndices.length - wml + 1;
            double[] xTerms = new double[maxPos]; 
            for (int p = 0; p < maxPos; ++p) {
                try {
	                xTerms[p] = scoreWM(bm, p);
                } catch (Exception ex) {
                    throw new IllegalArgumentException();
                }
            }
            double x = NativeMath.addLog2(xTerms);
            
            eta += wwm.getWeight() * (x - odds);
        }
        
        double pi = logit(eta);
        if (label >= 0) {
            return NativeMath.log2(pi);
        } else {
            return NativeMath.log2(1.0 - pi);
        }
    }

    private double logit(double x) {
        return 1.0 / (1 + Math.exp(-x));
    }

    
    private double scoreWM(Matrix2D iwm, int pos)
		throws Exception
	{
        double score = 0.0;
		int maxCol = iwm.columns();
		for (int col = 0; col < maxCol; ++col) {
		    int i = trimmedSeqIndices[pos + col];
		    if (i < 0) {
		        return Double.NEGATIVE_INFINITY;
		    } else {
		        score += iwm.get(i, col);
		    }
		}
		return score;
	}
    
}
