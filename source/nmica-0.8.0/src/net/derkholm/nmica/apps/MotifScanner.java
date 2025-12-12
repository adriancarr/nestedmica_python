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

package net.derkholm.nmica.apps;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.utils.CollectTools;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author thomas
 */

@App(overview="Scan sequences with a set of weight matrix motifs", generateStub=true)
@NMApp(launchName="nmscan", vm=VirtualMachine.SERVER)
public class MotifScanner {
    public static enum TargetStrand {
        FORWARD,
        REVERSE,
        BOTH
    }
    
    public static enum Format {
    	GFF,
    	FASTA,
    	TABLE
    }
    
    private static final double LOG_2 = Math.log(2.0);
    
    private SequenceDB seqs;
    private Motif[] motifs;
    private double scoreThreshold = Double.NaN;
    private double softness = 0.0;
    private TargetStrand strand = TargetStrand.BOTH;
    private double[] thresholdList;
    private boolean gapMap = false;
    private boolean maxPerSeq = false;
    private AlphabetIndex index;
    private Format format = Format.GFF;
    private boolean maxAlwaysOutput = false;
    
    @Option(help="When in -maxPerSeq mode, include output for every sequence", optional=true, userLevel=UserLevel.EXPERT)
    public void setMaxAlwaysOutput(boolean b) {
    	this.maxAlwaysOutput = b;
    }
    
    @Option(help="Output format (default=gff)", optional=true)
    public void setFormat(Format b) {
    	this.format = b;
    }
    
    @Option(help="Report the highest motif score found anywhere on an input sequence", optional=true)
    public void setMaxPerSeq(boolean b) {
        this.maxPerSeq = b;
    }
    
    @Option(help="Handle gaps in the sequence", optional=true, userLevel=UserLevel.EXPERT)
    public void setGapMap(boolean b) {
        this.gapMap = b;
    }
    
    @Option(help="Individually set thresholds for each motif", optional=true, userLevel=UserLevel.EXPERT)
    public void setThresholdList(Reader r)
        throws IOException
    {
        BufferedReader br = new BufferedReader(r);
        List<Double> l = new ArrayList<Double>();
        for (String line = br.readLine(); line != null; line = br.readLine()) {
            l.add(new Double(line));
        }
        thresholdList = CollectTools.toDoubleArray(l);
    }
    
    @Option(help="Strand of a DNA sequence to analyse (default=both)", optional=true)
    public void setStrand(TargetStrand t) {
        this.strand = t;
    }
    
    @Option(help="The lowest score to report", optional=true)
    public void setScoreThreshold(double d) {
        this.scoreThreshold = d;
    }
    
    @Option(help="The sequences to analyse", optional=false)
    public void setSeqs(Reader r) throws Exception {
        seqs = new HashSequenceDB();
        SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(r));
        while (si.hasNext()) {
            seqs.addSequence(si.nextSequence());
        }
    }
    
    @Option(help="The motifs to scan", optional=false)
    public void setMotifs(File f)
		throws Exception
	{
		if (f.getName().endsWith(".jos")) {
            throw new Exception("JOS format motif-sets are no longer supported");
            
            /*
            
		    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(f));
		    List<WeightMatrix> l = new ArrayList<WeightMatrix>();
		    try {
		        while (true) {
		            l.add((WeightMatrix) ois.readObject());
		        }
		    } catch (Exception ex) {
		    }
		    motifs = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
            
            */
		} else {
		    motifs = MotifIOTools.loadMotifSetXML(new FileInputStream(f));
		}
	}
    
    private void recordHit(GFFWriter gffw, Sequence seq, boolean isReverseStrand, int start, int stop, String name, double score)
        throws Exception
    {
    	if (format == Format.FASTA) {
    		SymbolList sl = seq.subList(start, stop);
    		if (isReverseStrand) {
    			sl = DNATools.reverseComplement(sl);
    		}
    		System.out.printf(">%s_%s_%d_%d%n", name, seq.getName(), start, stop);
    		System.out.println(sl.seqString());
    	} else if (format == Format.TABLE) {
    		SymbolList sl = seq.subList(start, stop);
    		if (isReverseStrand) {
    			sl = DNATools.reverseComplement(sl);
    		}
    		System.out.printf("%s\t%s\t%d\t%d\t%s\t%g%n", name, seq.getName(), start, stop, sl.seqString(), score);
    	} else {
	        SimpleGFFRecord record = new SimpleGFFRecord();
	        record.setSeqName(seq.getName());
	        record.setStart(start);
	        record.setEnd(stop);
	        record.setStrand(isReverseStrand ? StrandedFeature.NEGATIVE : StrandedFeature.POSITIVE);
	        record.setFeature(name);
	        record.setSource("MotifScanner");
	        record.setScore(score);
	        gffw.recordLine(record);
    	}
    }
    
    private interface HitRecorder {
        public void recordHit(boolean isReverseStrand, int min, int max, String name, double score)
            throws Exception;
    }
    
    public void main(String[] args)
    		throws Exception
    	{
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(System.out));
        final GFFWriter gffw = new GFFWriter(pw);
        
        boolean forward = (strand == TargetStrand.FORWARD || strand == TargetStrand.BOTH);
        boolean reverse = (strand == TargetStrand.REVERSE || strand == TargetStrand.BOTH);
        
        BM[] fwdMotifs = new BM[motifs.length];
        BM[] backMotifs = new BM[motifs.length];
        for (int w = 0; w < motifs.length; ++w) {
            if (index == null) {
                index = AlphabetManager.getAlphabetIndex((FiniteAlphabet) motifs[w].getWeightMatrix().getAlphabet());
            }
            fwdMotifs[w] = wmToBm(index, motifs[w].getWeightMatrix());
            backMotifs[w] = wmToBm(index, WmTools.reverseComplement(motifs[w].getWeightMatrix()));
        }
        
        for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
            final Sequence rawSl = si.nextSequence();
            byte[] sl;
            HitRecorder hr;
            sl = indexSeq(index, rawSl);
            hr = new HitRecorder() {
                public void recordHit(boolean isReverseStrand, int min, int max, String name, double score) 
                    throws Exception
                {
                    MotifScanner.this.recordHit(gffw, rawSl, isReverseStrand, min, max, name, score);
                }
            };
                
            for (int w = 0; w < motifs.length; ++w) {
                int maxPos = sl.length - motifs[w].getWeightMatrix().columns();
                BM fwdMatrix = fwdMotifs[w];
                BM backMatrix = backMotifs[w];
                double maxScore = fwdMatrix.getMaxScore();
                double threshold = (thresholdList == null ? scoreThreshold : thresholdList[w]);
                if (Double.isNaN(threshold)) {
                	threshold = motifs[w].getThreshold();
                }
                
                if (maxPerSeq) {
                    double max = maxAlwaysOutput ? threshold : Double.NEGATIVE_INFINITY;
                    for (int pos = 1; pos <= maxPos; ++pos) {
                         if (forward) {
                            double bitsSubOptimal = maxScore - fwdMatrix.score(sl, pos, gapMap);
                            max = Math.max(max, -bitsSubOptimal);
                        }
                        if (reverse) {
                            double bitsSubOptimal = maxScore - backMatrix.score(sl, pos, gapMap);
                            max = Math.max(max, -bitsSubOptimal);
                        }
                    }
                    
                    if (max >= threshold) {
                        System.out.println(rawSl.getName() + "\t" + motifs[w].getName() + "\t" + max);
                    }
                } else {
                    for (int pos = 1; pos <= maxPos; ++pos) {
                         if (forward) {
                            double bitsSubOptimal = maxScore - fwdMatrix.score(sl, pos, gapMap);
                            
                            if (-bitsSubOptimal >= threshold) {
                                hr.recordHit(false, pos, fwdMatrix.endPos(), motifs[w].getName(), -bitsSubOptimal);
        	                	   }
                        }
                        if (reverse) {
                            double bitsSubOptimal = maxScore - backMatrix.score(sl, pos, gapMap);
                        
                            if (-bitsSubOptimal >= threshold) {
                                hr.recordHit(true, pos, backMatrix.endPos(), motifs[w].getName(), -bitsSubOptimal);
                            }
                        }
                    }
                }
            }
        }
        pw.flush();
    	}

    private class BM {
        private Matrix2D bm;
        private AlphabetIndex index;
        private String name;
        private double maxScore;
        private int endPos;
        
        public int endPos() {
            return endPos;
        }
        
        BM(Matrix2D bm, AlphabetIndex index, String name) 
            throws Exception
        {
            this.bm = bm;
            this.index = index;
            this.name = name;
            this.maxScore = bmMaxScore(bm);
        }
        
        public double getMaxScore() {
            return maxScore;
        }
        
        public double score(byte[] sl, int pos, boolean gm)
            throws Exception
        {
            return scoreBM(sl, bm, pos, gm);
        }
        
        private  double scoreBM(byte[] sl, Matrix2D bm, int pos, boolean gm)
                throws Exception
        {
            double score = 0.0;
            int scol = 0;
            try {
                for (int col = 0; col < bm.columns(); ++col, ++scol) {
                    byte s = sl[pos + scol];
                    while (gm && scol > 0 && s < 0) {
                        s = sl[pos + ++scol];
                    }
                    if (s >= 0) {
                        score += bm.get(s, col);
                    } else {
                        return Double.NEGATIVE_INFINITY;
                    }
                }
            } catch (IndexOutOfBoundsException ex) {
                return Double.NEGATIVE_INFINITY;
            }
            endPos = pos + scol - 1;
            return score;
        }
        
        
        private double bmMaxScore(Matrix2D bm)
                throws Exception
        {
            double wmScore = 0.0;
            for (int c = 0; c < bm.columns(); ++c) {
                double colScore = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < bm.rows(); ++i) {
                    colScore = Math.max(colScore, bm.get(i, c));
                }
                wmScore += colScore;
            }
            return wmScore;
        }
    }
    
    private BM wmToBm(AlphabetIndex index, WeightMatrix wm)
    		throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
        Matrix2D bm = new SimpleMatrix2D(alpha.size(), wm.columns());
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution wmCol = wm.getColumn(c);
            
            for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
                Symbol s = (Symbol) si.next();
                double baseWeight = wmCol.getWeight(s);
                double moderatedWeight = (softness / 4) + baseWeight * (1.0 - softness);
                bm.set(index.indexForSymbol(s), c, Math.log(moderatedWeight) / LOG_2);
            }
        }
        return new BM(bm, index, null);
    }
    

    
    private SymbolList degappify(SymbolList sl)
        throws Exception
    {
        Symbol gap = sl.getAlphabet().getGapSymbol();
        List<Symbol> l = new ArrayList<Symbol>();
        for (Iterator<?> si = sl.iterator(); si.hasNext(); ) {
            Symbol s = (Symbol) si.next();
            if (s != gap) {
                l.add(s);
            }
        }
        return new SimpleSymbolList(sl.getAlphabet(), l);
    }
    
    private byte[] indexSeq(AlphabetIndex index, SymbolList sl)
        throws Exception
    {
        Symbol gap = sl.getAlphabet().getGapSymbol();
        
        byte[] bsl = new byte[sl.length() + 1];
        for (int i = 1; i <= sl.length(); ++i) {
            Symbol s = sl.symbolAt(i);
            if (s == gap) {
                bsl[i] = -2;
            } else if (s instanceof AtomicSymbol) {
                bsl[i] = (byte) index.indexForSymbol(s);
            } else {
                bsl[i] = -1;
            }
        }
        return bsl;
    }
}
