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

import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.maths.MathsTools;
import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.model.ContributionItem;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.utils.CollectTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Calculate likelihood functions for short sequence motifs against a mosaic-
 * type background.
 *
 * @author Thomas Down
 */

class MotifDiscriminitiveLikelihood implements LikelihoodCalculator {
    private final MotifFacette facette;
    private final SymbolList data;
    private final int[] trimmedSeqIndices;
    private final AlphabetIndex index;
    private final int maxIndex;
    private double bgHood;
    private double[] symbolHood;
    private double[] symbolHoodLog;
    
    MotifDiscriminitiveLikelihood(MotifFacette mmf, SymbolList seq) 
        throws Exception
    {
        this.facette = mmf;
        this.data = seq;
        
        FiniteAlphabet fa = (FiniteAlphabet) seq.getAlphabet();
        this.index = AlphabetManager.getAlphabetIndex(fa);
        this.maxIndex = fa.size();
        
        double[] tempHood = new double[seq.length()];
        Location seqMask = facette.getBackground().backgroundSymbolLikelihood(seq, tempHood);
        
        List<Integer> indexList = new ArrayList<Integer>();
        List<Double> hoodList = new ArrayList<Double>();
        bgHood = 0;
        for (int pos = 1; pos <= seq.length(); ++pos) {
            if (seqMask.contains(pos)) {
                hoodList.add(new Double(tempHood[pos - 1]));
                bgHood += tempHood[pos - 1];
                Symbol s = seq.symbolAt(pos);
                if (s instanceof AtomicSymbol) {
                    indexList.add(new Integer(index.indexForSymbol(s)));
                } else {
                    indexList.add(new Integer(-1));
                }
            }
        }
        
        trimmedSeqIndices = CollectTools.toIntArray(indexList);
        symbolHoodLog = CollectTools.toDoubleArray(hoodList);
        
        symbolHood = new double[symbolHoodLog.length];
        for (int i = 0; i < symbolHood.length; ++i) {
            symbolHood[i] = NativeMath.exp2(symbolHoodLog[i]);
        }
    }
    
    public Facette getFacette() {
        return facette;
    }
    
    public Object getData() {
        return data;
    }
    
    public double likelihood(ObjectMatrix1D contributions, Matrix1D mix) {
        int card = contributions.size(); 
        int target = 0;
        int offTarget = 0;
        for (int i = contributions.size() - 1; i >= 0; --i) {
        	target = target << 1;
        	offTarget = offTarget << 1;
            if (mix.get(i) != 0.0) {
                target |= 1;
            } else {
            	offTarget |= 1;
            }
        }
        // System.err.println("Sequence = " + ((Sequence) data).getName() + " Target = " + target);
        
        try {
        	/*
        	
            Matrix2D[] motifs = new Matrix2D[card];
            for (int c = 0; c < card; ++c) {
                motifs[c] = motifToMatrix((WeightMatrix) ((ContributionItem) contributions.get(c)).getItem());
            }
            
            Matrix2D wmScores = new SimpleMatrix2D(trimmedSeqIndices.length, motifs.length);
            int[] advances = new int[motifs.length];
            for (int m = 0; m < motifs.length; ++m) {
                Matrix2D iwm = motifs[m];
                int wml = iwm.columns();
                advances[m] = wml;
                int maxPos = trimmedSeqIndices.length - wml + 1;
                for (int i = 0; i < maxPos; ++i) {
                    double rscore = scoreWM(iwm, i);
                    double score = NativeMath.log2(rscore);
                    wmScores.set(i + wml - 1, m, score);
                }
            }
            
            */
            
        	Matrix2D[] motifs;
            if (facette.getRevComp()) {
            	motifs = new Matrix2D[card << 1];
	            for (int m = 0, c = 0; m < card; ++m) {
	            	ContributionItem ci = (ContributionItem) contributions.get(m);
	                motifs[c++] = (Matrix2D) ci.getItemView(facette.getForwardBitMatrixView());
	                motifs[c++] = (Matrix2D) ci.getItemView(facette.getReverseBitMatrixView());
	            }
            } else {
	            motifs = new Matrix2D[card];
	            for (int m = 0, c = 0; c < card; ++m) {
	            	ContributionItem ci = (ContributionItem) contributions.get(m);
	                motifs[c++] = (Matrix2D) ci.getItemView(facette.getForwardBitMatrixView());
	            }
            }
            
            Matrix2D wmScores = new SimpleMatrix2D(trimmedSeqIndices.length, motifs.length);
            int[] advances = new int[motifs.length];
            for (int m = 0; m < motifs.length; ++m) {
                Matrix2D iwm = motifs[m];
                int wml = iwm.columns();
                advances[m] = wml;
                int maxPos = trimmedSeqIndices.length - wml + 1;
                for (int i = 0; i < maxPos; ++i) {
                    double score = scoreWM(iwm, i);
                    wmScores.set(i + wml - 1, m, score);
                }
            }
            
            
            double score =  weaveMotifs(symbolHoodLog, wmScores, advances, 1.0 * motifs.length / symbolHood.length, facette.getSoftness(), target, offTarget);
            // System.err.print("." + score);
            return score;
        } catch (Exception ex) {
            ex.printStackTrace();
            throw new IllegalArgumentException("Error calculating likelihood");
        }
    }
    
    private double weaveMotifs(double[] bgScores, Matrix2D wmScores, int[] advances, double motifTransition, double softness, int target, int offTarget) {
        int length = bgScores.length;
        int numMotifs = wmScores.columns();
        boolean revComp = facette.getRevComp();
        if (revComp) {
        	numMotifs = numMotifs >> 1;
        }
        int stateSpace = 1 << numMotifs;
        
        double basePenalty = NativeMath.log2(1.0 - motifTransition);
        double[] motifPenalty = new double[numMotifs];
        for (int n = 0; n < numMotifs; ++n) {
            motifPenalty[n] = NativeMath.log2(motifTransition / (numMotifs - n));
        }
        
        Matrix2D matrix = new SimpleMatrix2D(length + 1, stateSpace);
        matrix.set(0, 0, 0.0); // at the start of the sequence, no motifs used.
        for (int s = 1; s < stateSpace; ++s) {
            matrix.set(0, s, Double.NEGATIVE_INFINITY);
        }
        
        if (revComp) {
        	javaWeaveBidiDP(matrix, advances, basePenalty, motifPenalty, bgScores, wmScores);
        } else {
        	javaWeaveDP(matrix, advances, basePenalty, motifPenalty, bgScores, wmScores);
        }
        // System.out.printf("Target=%d, on=%g, off=%g%n", target, matrix.get(length, target), matrix.get(length, target ^ 1));
        
        double onTargetScore = matrix.get(length, target);
        double offTargetScore = matrix.get(length, offTarget);
        if (offTargetScore == Double.NEGATIVE_INFINITY) {
        	return Double.NEGATIVE_INFINITY;
        } else {
        	return (2 * onTargetScore) - offTargetScore;
        }
    }
    
    protected void weaveDP(Matrix2D matrix, int[] advances, double basePenalty, double[] motifPenalty, double[] bgScores, Matrix2D wmScores) {
        javaWeaveDP(matrix, advances, basePenalty, motifPenalty, bgScores, wmScores);
    }
    
    private void javaWeaveDP(Matrix2D matrix, int[] advances, double basePenalty, double[] motifPenalty, double[] bgScores, Matrix2D wmScores) {
        int length = bgScores.length;
        int numMotifs = wmScores.columns();
        int stateSpace = 1 << numMotifs;
        
        for (int i = 1; i <= length; ++i) {
            for (int s = 0; s < stateSpace; ++s) {
                int iminus = i - 1;
                double score = matrix.get(iminus, s) + bgScores[iminus] + basePenalty;
                for (int m = 0; m < numMotifs; ++m) {
                    int wml = advances[m];
                    if (i >= wml) {
                        int mState = 1 << m;
                        if ((s & mState) != 0) {
                            int from = s & ~mState;
                            double fromScore = matrix.get(i - wml, from);
                            double emitScore =  wmScores.get(iminus, m);
                            score = NativeMath.addLog2(score, fromScore + emitScore + motifPenalty[MathsTools.popcnt(from)]);
                        }
                    }
                }
                matrix.set(i, s, score);
            }
        }
    }
    
    private void javaWeaveBidiDP(Matrix2D matrix, int[] advances, double basePenalty, double[] motifPenalty, double[] bgScores, Matrix2D wmScores) {
        int length = bgScores.length;
        int numMotifs = wmScores.columns() >> 1;
        int stateSpace = 1 << numMotifs;
        
        for (int i = 1; i <= length; ++i) {
            for (int s = 0; s < stateSpace; ++s) {
                int iminus = i - 1;
                double score = matrix.get(iminus, s) + bgScores[iminus] + basePenalty;
                for (int m = 0; m < numMotifs; ++m) {
                    int wml = advances[m];
                    if (i >= wml) {
                        int mState = 1 << m;
                        if ((s & mState) != 0) {
                            int from = s & ~mState;
                            int mdash = m << 1;
                            double fromScore = matrix.get(i - wml, from);
                            double emitScore =  NativeMath.addLog2(wmScores.get(iminus, mdash), wmScores.get(iminus, mdash + 1)) - 1;
                            score = NativeMath.addLog2(score, fromScore + emitScore + motifPenalty[MathsTools.popcnt(from)]);
                        }
                    }
                }
                matrix.set(i, s, score);
            }
        }
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
