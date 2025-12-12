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

import net.derkholm.nmica.maths.DoubleFunction;
import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.model.ContributionItem;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.utils.ArrayTools;
import net.derkholm.nmica.utils.CollectTools;

import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Calculate likelihood functions for short sequence motifs against a mosaic-
 * type background.
 * 
 * <p>
 * <strong>WARNING:</strong> MotifUncountedLikelihood instances are not threadsafe.
 * </p>
 *
 * @author Thomas Down
 */

class MotifUncountedLikelihood implements LikelihoodCalculator {
    private final MotifFacette facette;
    private final SymbolList data;
    private final int[] trimmedSeqIndices;
    private double bgHood;
    private double[] symbolHood;
    private double[] symbolHoodLog;
    
    // storage for the weaving step is pre-allocated so we're nice to the GC
    
    private double[] matrix;
    private Matrix2D wmScores;
    private Matrix2D[] lastSeenCache;
    
    // Only used in cache-instrumentation mode.
    
    private static boolean instrumentCache;
    private static int cacheHits = 0;
    private static int cacheMisses = 0;
    private static int cacheNearMissesHigher = 0;
    private static int cacheNearMissesLower = 0;
    private static int cacheNearMissCycles = 0;
    private static int cacheNearMissesAmbig = 0;
    
    static {
        instrumentCache = Boolean.getBoolean("nmica.instrument_cache");
        if (instrumentCache) {
            System.err.println("Motif facette cache instrumentation is enabled");
            System.err.println("Warning: this probably isn't threadsafe...");
            Runtime.getRuntime().addShutdownHook(new Thread() {
                public void run() {
                    System.err.println("Cache instrumentation");
                    System.err.println("---------------------");
                    System.err.println();
                    System.err.printf("Hits = %d\n", cacheHits);
                    System.err.printf("Misses = %d\n", cacheMisses);
                    System.err.printf("Near-misses(higher) = %d\n", cacheNearMissesHigher);
                    System.err.printf("Near-misses(lower) = %d\n", cacheNearMissesLower);
                    System.err.printf("Cycles with at least one near-miss = %d\n", cacheNearMissCycles);
                    System.err.printf("Cycles with both higher and lower misses = %d\n", cacheNearMissesAmbig);
                }
            }  );
        }
    }
    
    MotifUncountedLikelihood(MotifFacette mmf, SymbolList seq) 
        throws Exception
    {
        // System.err.println("Initializing a " + getClass().getName());
        
        this.facette = mmf;
        this.data = seq;
        
        AlphabetIndex index = facette.getAlphabetIndex();
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
        
        matrix = new double[symbolHood.length + 1];
        
        // System.err.println("Background likelihood is " + bgHood);
    }
    
    public MotifFacette getMotifFacette() {
    	return facette;
    }
    
    public Facette getFacette() {
        return facette;
    }
    
    public Object getData() {
        return data;
    }
    
    public double likelihood(ObjectMatrix1D contributions, Matrix1D mix) {
        double mixThreshold = 0.001;
        
        DoubleFunction mtf = facette.getMixTransferFunction();
        double[] mixVals = new double[contributions.size()];
        int card = 0;
        for (int i = 0; i < contributions.size(); ++i) {
            mixVals[i] = mtf.eval(mix.get(i));
            if (mixVals[i] > mixThreshold) {
                ++card;
            }
        }
        
        if (card == 0) {
            return bgHood;
        } else {
            try {
                double edgePruneThreshold = facette.getEdgePruneThreshold();
                Matrix2D[] motifs;
                double[] transitions;
                if (facette.getRevComp()) {
                    motifs = new Matrix2D[card << 1];
                    transitions = new double[card << 1];
                    for (int m = 0, c = 0; c < card; ++m) {
                        if (mixVals[m] > mixThreshold) {
                            int i = (c++) << 1;
                            ContributionItem ci = (ContributionItem) contributions.get(m);
                            motifs[i] = (Matrix2D) ci.getItemView(facette.getForwardBitMatrixView());
                            motifs[i + 1] = (Matrix2D) ci.getItemView(facette.getReverseBitMatrixView());
                            transitions[i] = transitions[i + 1] = mixVals[m] * facette.getUncountedExpectation() / symbolHood.length / 2;
                        }
                    }
                } else {
                    motifs = new Matrix2D[card];
                    transitions = new double[card];
                    for (int m = 0, c = 0; c < card; ++m) {
                        if (mixVals[m] > mixThreshold) {
                            int i = c++;

                            ContributionItem ci;
                            try {
                                ci = (ContributionItem) contributions.get(m);
                            } catch (ClassCastException ex) {
                                System.err.println("Not a ContributionItem, got a " + contributions.get(m).getClass().getName() + " instead!");
                                throw ex;
                            }
                            motifs[i] = (Matrix2D) ci.getItemView(facette.getForwardBitMatrixView());
                            transitions[i] = mixVals[m] * facette.getUncountedExpectation() / symbolHood.length;
                            // System.err.println("transitions[" + i + "] = " + transitions[i]);
                        }
                        
                    }
                }
                                
                if (wmScores == null || wmScores.columns() < motifs.length) {
                    /* Matrix2D */ wmScores = new SimpleMatrix2D(trimmedSeqIndices.length, motifs.length);
                    lastSeenCache = new Matrix2D[motifs.length];
                }
                
                if (instrumentCache) {
                    boolean hasHigher = false;
                    boolean hasLower = false;
                    for (int m = 0; m < motifs.length; ++m) {
                        int hitIndex = ArrayTools.indexOf(lastSeenCache, motifs[m]);
                        if (hitIndex < 0) {
                            ++cacheMisses;
                        } else if (hitIndex == m) {
                            ++cacheHits;
                        } else if (hitIndex > m){
                            ++cacheNearMissesHigher;
                            hasHigher = true;
                        } else {
                            ++cacheNearMissesLower;
                            hasLower = true;
                        }
                    }
                    if (hasHigher || hasLower) {
                        ++cacheNearMissCycles;
                    }
                    if (hasHigher && hasLower) {
                        ++cacheNearMissesAmbig;
                    }
                }
                
                int[] advances = new int[motifs.length];
                for (int m = 0; m < motifs.length; ++m) {
                    Matrix2D iwm = motifs[m];
                    int wml = iwm.columns();
                    advances[m] = wml;
                    
                    int cacheHit = -1;
                  CACHESEARCH_LOOP:
                    for (int cm = m; cm < lastSeenCache.length; ++cm) {
                        if (lastSeenCache[cm] == iwm) {
                            cacheHit = cm;
                            break CACHESEARCH_LOOP;
                        }
                    }
                    
                    if (cacheHit == m) {
                        // nothing to do :-)
                    } else if (cacheHit > m) {
                        // near-miss, copy a cache line across
                        for (int i = 0; i < trimmedSeqIndices.length; ++i) {
                            wmScores.set(i, m, wmScores.get(i, cacheHit));
                        }
                        lastSeenCache[m] = iwm;
                    } else {
                        // full-miss, need to recalculate
                        int maxPos = trimmedSeqIndices.length - wml + 1;
                        for (int i = 0; i < maxPos; ++i) {
                            double score = scoreWM(iwm, i);
                            wmScores.set(i + wml - 1, m, score);
                        }
                        lastSeenCache[m] = iwm;
                    }
                }

                double score = sumMotifs(symbolHoodLog, motifs.length, wmScores, advances, transitions);
                // System.err.println("Score = " + score);
                return score;
            } catch (Exception ex) {
                ex.printStackTrace();
                throw new IllegalArgumentException("Error calculating likelihood");
            }
        }
    }
    
    protected double sumMotifs(double[] bgScores, int numMotifs, Matrix2D wmScores, int[] advances, double[] motifTrans)
        throws Exception
    {
        int length = bgScores.length;
        matrix[0] = 0;
        
        double sumTrans = 0;
        double[] motifPenalty = new double[numMotifs];
        for (int m = 0; m < numMotifs; ++m) {
            motifPenalty[m] = NativeMath.log2(motifTrans[m]);
            sumTrans += motifTrans[m];
        }
        double basePenalty = NativeMath.log2(1.0 - sumTrans);
        
        for (int i = 1; i <= length; ++i) {
            double score = matrix[i - 1] + bgScores[i - 1] + basePenalty;
            for (int m = 0; m < numMotifs; ++m) {
                int wml = advances[m];
                if (i >= wml) {
                    double fromScore = matrix[i - wml];
                    double emitScore =  wmScores.get(i - 1, m);
                    score = NativeMath.addLog2(score, fromScore + emitScore + motifPenalty[m]);
                }
            }
            matrix[i] = score;
        }
        
        return matrix[length];
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
