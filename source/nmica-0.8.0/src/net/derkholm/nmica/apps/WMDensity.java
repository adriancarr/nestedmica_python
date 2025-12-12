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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

public class WMDensity{
    
    public static void printUsage(){
        System.err.println("Usage: java net.derkholm.nmica.apps.WMDensity <options>");
        System.err.println("Where the options are:");
        System.err.println(" -seq			a file of sequences to be scanned");
        System.err.println(" -motifs 			a weightmatrix file in jos format");
        System.err.println(" -threshold 	the weightmatrix threshold required");
        System.err.println(" -5 prime		count from start of sequence (default");
        System.err.println(" -3 prime		count from end of sequence");
        System.err.println(" -bucketSize    default 15");
        System.err.println(" -numBuckets    default 10");
        System.err.println(" -minSeqLength  default 150");
        System.err.println(" -maxSeqLength  default 300");
        System.err.println(" -motifNumber   which motif in the list of N motifs (starts from 0)");
    }
    
    private static final double LOG_2 = Math.log(2);
    
    public static void main(String[] args)
    throws Exception{

        File seqFile = null;
        InputStream motifsFile = null;
        double threshold = 1.0;
        boolean threePrime = false;
        int bucketSize = 15;
        int numBuckets = 10;
        int minSeqLength = 150;
        int maxSeqLength = 300;
        int motifNumber = 0;
        try {
            for (int i = 0; i < args.length; ++i) {
                if (args[i].equals("-seq") || args[i].equals("-fasta")) {
                    seqFile = new File(args[++i]);
                } else if(args[i].equals("-motifs")){ 
                    motifsFile = new FileInputStream(args[++i]);
                } else if (args[i].equals("-threshold")) {
                    threshold = Double.parseDouble(args[++i]);
                } else if (args[i].equals("-3prime")) {
                    threePrime = true;
                } else if (args[i].equals("-5prime")) {
                    threePrime = false;
                } else if (args[i].equals("-bucketSize")) {
                    bucketSize = Integer.parseInt(args[++i]);
                } else if (args[i].equals("-numBuckets")) {
                    numBuckets = Integer.parseInt(args[++i]);
                } else if (args[i].equals("-minSeqLength")) {
                    minSeqLength = Integer.parseInt(args[++i]);
                } else if (args[i].equals("-maxSeqLength")) {
                    maxSeqLength = Integer.parseInt(args[++i]);
                } else if (args[i].equals("-motifNumber")) {
                	   motifNumber = Integer.parseInt(args[++i]);
                } else {
                    System.err.println("Unrecognized option " + args[i]);
                    printUsage();
                    return;
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            printUsage();
            return;
        }
        
        if (seqFile == null || motifsFile == null) {
            printUsage();
            return;
        }

        SequenceDB seqs = loadDB(seqFile, minSeqLength, maxSeqLength);
        WeightMatrix[] wms;
        {
            ObjectInputStream ois = new ObjectInputStream(motifsFile);
            List l = new ArrayList();
            try {
                while (true) {
                    l.add(ois.readObject());
                }
            } catch (Exception ex) {
            }
            wms = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
        }

        double[] realProfile = wmProfile(seqs, wms[motifNumber], threshold, bucketSize, numBuckets, threePrime);
        for (int i = 0; i < numBuckets; ++i) {
            System.out.println(
                "" +
                i + '\t' +
                realProfile[i] + '\t' 
            );
        }

	
    }
 
    
    private static double[] wmProfile(SequenceDB seqs, WeightMatrix wm, double threshold, 
            							int bucketSize, int buckets, boolean threePrime)
    throws Exception
    {
        int maxPos = bucketSize * buckets;
        double[] counts = new double[buckets];
        int num = 0;
       
        for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
            
            Sequence seq = (Sequence) si.nextSequence();
            //double bitsSubOptimal = -Math.log(score / wmMaxScore(motifs[w])) / LOG_2;
            int seqMaxPos = Math.min(maxPos, seq.length() - wm.columns() + 1);
            
            for (int pos = 1; pos <= seqMaxPos; ++pos) {
                if (threePrime) {
                  if (-Math.log(scoreWMprobs(wm, seq, seq.length() - pos - wm.columns() + 2)/wmMaxScore(wm))/LOG_2 <= threshold) {
                        counts[(pos - 1) / bucketSize] += 1.0;                        	
                  } 
                } else {
                    double bitsSubOptimal = -Math.log(scoreWMprobs(wm, seq, pos)/wmMaxScore(wm))/LOG_2;
                    // System.err.println("bso = " + bitsSubOptimal);
                    if (bitsSubOptimal <= threshold) {
                        counts[(pos - 1) / bucketSize] += 1.0;
                    }
                  
               }
               ++num;
            }
        }
        
        for (int i = 0; i < counts.length; ++i) {
            counts[i] /= num;
        }
        return counts;
    }

    private static double wmMaxScore(WeightMatrix wm)
	throws Exception
	{
        double wmScore = 1.0;
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution col = wm.getColumn(c);
            double colScore = 0;
            for (Iterator i = ((FiniteAlphabet) col.getAlphabet()).iterator(); i.hasNext(); ) {
                colScore = Math.max(colScore, col.getWeight((Symbol) i.next()));
            }
            wmScore *= colScore;
        }
        return wmScore;
	}

    public static double scoreWMprobs(	WeightMatrix twm, Sequence seq, int pos)
    throws Exception
    {
        double score = 1.0;
        int seqLength = seq.length();
        
        for (int i = 0; i < twm.columns(); ++i) {
            int ipos = pos + i;
            if (ipos < 1 || ipos > seqLength) {
                return 0;
            } else {
                Symbol s = seq.symbolAt(pos + i);
                if (! (s instanceof AtomicSymbol)) {
                    return 0;
                }
                score = score * twm.getColumn(i).getWeight(s);
            }
        }
        return score;
    }
    
    
    
    private static SequenceDB loadDB(File f, int minLength, int maxLength)
    throws Exception{
        int cnt = 0;
        SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(f)));
        SequenceDB seqDB = new HashSequenceDB();
        while(si.hasNext()){
            Sequence seq = (Sequence) si.nextSequence();
           
            if(seq.length() >= minLength && seq.length() <= maxLength){
                seqDB.addSequence(seq);
                ++cnt;
            }
        }
        
        System.err.println("Loaded " + cnt + " sequences");
        
        return seqDB;
    }
    	
}