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

import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Calculate likelihood functions for short sequence motifs against a mosaic-
 * type background.
 *
 * @author Thomas Down
 */

class MotifLikelihood implements LikelihoodCalculator {
    private final MotifFacette facette;
    private final SymbolList data;
    private final int[] trimmedSeqIndices;
    private final AlphabetIndex index;
    private final int maxIndex;
    private double bgHood;
    private double[] symbolHood;
    private double[] symbolHoodLog;
    
    MotifLikelihood(MotifFacette mmf, SymbolList seq) 
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
        // OPTIMIZATION: work out the cardinality first (is this sane?)
        
        int card = 0;
        for (int i = 0; i < contributions.size(); ++i) {
            if (mix.get(i) != 0.0) {
                ++card;
            }
        }
        
        if (card == -1) {
            return bgHood;
  /*      } else if (card == 1) {
            try {
                double motifTransition = 1.0 / symbolHood.length;
                Matrix2D iwm = null;
                for (int i = 0; i < contributions.size(); ++i) {
                    if (mix.get(i) != 0.0) {
                    	ContributionItem ci = (ContributionItem) contributions.get(i);
                        iwm = (Matrix2D) ci.getItemView(facette.getForwardBitMatrixView());
                        break;
                    }
                }
                
                double[] sh = symbolHoodLog;
                double bg = bgHood;
                
                double tot = 0;
                int wml = iwm.columns();
                int maxPos = trimmedSeqIndices.length - wml + 1;
                
                for (int pos = 0; pos < maxPos; ++pos) {
                    double score = scoreWM(iwm, pos);
                    double bs = 1.0;
                    for (int z = 0; z < wml; ++z) {
                        bs += sh[pos + z];
                    }
                    
                    tot = NativeMath.addLog2(tot, score - bs);
                }
                
                double hardL = bg + tot + (maxPos - 1) * NativeMath.fastlog2(1.0 - motifTransition) + NativeMath.fastlog2(motifTransition);
                double softness = facette.getSoftness();
                if (softness == 0) {
                    return hardL;
                } else {
                    return NativeMath.addLog2(hardL + NativeMath.log2(1.0 - softness), bg + NativeMath.log2(softness));
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                throw new IllegalArgumentException("Error calculating likelihood");
            }
            */
            
        } else {
            try {
                Matrix2D[] motifs = new Matrix2D[card];
                for (int m = 0, c = 0; c < card; ++m) {
                    if (mix.get(m) != 0.0) {
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
                
                return weaveMotifs(symbolHoodLog, wmScores, advances, 1.0 * motifs.length / symbolHood.length, facette.getSoftness());
            } catch (Exception ex) {
                ex.printStackTrace();
                throw new IllegalArgumentException("Error calculating likelihood");
            }
        }
    }
    
    private double weaveMotifs(double[] bgScores, Matrix2D wmScores, int[] advances, double motifTransition, double softness) {
        int length = bgScores.length;
        int numMotifs = wmScores.columns();
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
        
        
        /*
        
        javaWeaveDP(matrix, advances, basePenalty, motifPenalty, bgScores, wmScores);
        {
            double marginaL;
            if (softness == 0.0) {
                marginaL = matrix.get(length, stateSpace - 1);
            } else {
                double tot = Double.NEGATIVE_INFINITY;
                for (int s = 0; s < stateSpace; ++s) {
                    int pop = MathsTools.popcnt(s);
                    tot = MathsTools.addLog(tot, matrix.get(length, s) + pop * Math.log(1.0 - softness) + (numMotifs - pop) * Math.log(softness));
                }
                marginaL = tot;
            }
            System.err.println("Pure java: " + marginaL);
        }
        
        */
        
        weaveDP(matrix, advances, basePenalty, motifPenalty, bgScores, wmScores);
        if (softness == 0.0) {
            return matrix.get(length, stateSpace - 1);
        } else {
            double tot = Double.NEGATIVE_INFINITY;
            for (int s = 0; s < stateSpace; ++s) {
                int pop = MathsTools.popcnt(s);
                tot = NativeMath.addLog2(tot, matrix.get(length, s) + pop * NativeMath.log2(1.0 - softness) + (numMotifs - pop) * NativeMath.log2(softness));
            }
            return tot;
        }
        // System.err.println("Native: " + marginaL);
        // return marginaL;
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
                            score = MathsTools.addLog(score, fromScore + emitScore + motifPenalty[MathsTools.popcnt(from)]);
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
