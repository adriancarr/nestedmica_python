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

import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.utils.CliTools;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * @author thomas
 */
public class MotifPickThresholds {
    private static final double LOG_2 = Math.log(2);
    
    private SequenceDB seqs;
    private WeightMatrix[] motifs;
    private double prop = 0.001;
    
    public void setProp(double d) {
        this.prop = d;
    }
    
    public void setSeqs(Reader r) throws Exception {
        seqs = new HashSequenceDB();
        SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(r));
        while (si.hasNext()) {
            seqs.addSequence(si.nextSequence());
        }
    }
    
    public void setMotifs(InputStream is)
    		throws Exception
    {
        ObjectInputStream ois = new ObjectInputStream(is);
        List l = new ArrayList();
        try {
            while (true) {
                l.add(ois.readObject());
            }
        } catch (Exception ex) {
        }
        motifs = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
    }
    
    public static void main(String[] args) throws Exception {
        MotifPickThresholds app = new MotifPickThresholds();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        int[][] b = new int[motifs.length][];
        for (int m = 0; m < motifs.length; ++m) {
            WeightMatrix fwdMatrix = wmToBm(motifs[m]);
            double maxScore = bmMaxScore(fwdMatrix);
            List<Double> scoreList = new ArrayList<Double>();
            for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
                Sequence seq = si.nextSequence();
                int maxPos = seq.length() - motifs[m].columns() + 1;
                for (int pos = 1; pos <= maxPos; ++pos) {
                    double bitsSubOptimal = maxScore - scoreBM(seq, fwdMatrix, pos);
                    scoreList.add(new Double(bitsSubOptimal));
                }
            }
            Collections.sort(scoreList);
            System.out.println(scoreList.get((int) Math.floor(prop * scoreList.size())));
        }
        
    }
    
    private WeightMatrix wmToBm(WeightMatrix wm)
    		throws Exception
    {
        double softness = 0;
        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
        
        WeightMatrix bm = new SimpleWeightMatrix(wm.getAlphabet(), wm.columns(), NMSimpleDistribution.FACTORY);
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution wmCol = wm.getColumn(c);
            Distribution bmCol = bm.getColumn(c);
            
            for (Iterator si = alpha.iterator(); si.hasNext(); ) {
                Symbol s = (Symbol) si.next();
                double baseWeight = wmCol.getWeight(s);
                double moderatedWeight = (softness / 4) + baseWeight * (1.0 - softness);
                bmCol.setWeight(s, Math.log(moderatedWeight) / LOG_2);
            }
        }
        return bm;
    }
    
    private static double scoreBM(SymbolList sl, WeightMatrix wm, int pos)
    		throws Exception
    {
        double score = 0.0;
        Symbol gap = sl.getAlphabet().getGapSymbol();
        try {
            for (int col = 0, scol = 0; col < wm.columns(); ++col, ++scol) {
                Symbol s = sl.symbolAt(pos + scol);
                while (s == gap) {
                    s = sl.symbolAt(pos + ++scol);
                }
                if (s instanceof AtomicSymbol) {
                    score += wm.getColumn(col).getWeight(s);
                } else {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        } catch (IndexOutOfBoundsException ex) {
            return Double.NEGATIVE_INFINITY;
        }
        return score;
    }
    
    
    private static double bmMaxScore(WeightMatrix wm)
    		throws Exception
    {
        double wmScore = 0.0;
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution col = wm.getColumn(c);
            double colScore = Double.NEGATIVE_INFINITY;
            for (Iterator i = ((FiniteAlphabet) col.getAlphabet()).iterator(); i.hasNext(); ) {
                colScore = Math.max(colScore, col.getWeight((Symbol) i.next()));
            }
            wmScore += colScore;
        }
        return wmScore;
    }

}
