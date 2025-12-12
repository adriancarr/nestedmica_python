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
import net.derkholm.nmica.model.motif.MosaicTools;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import java.io.BufferedReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

@App(overview="Perform a parameter sweep over possible mosaic background models", generateStub=true)
@NMApp(launchName="nmevaluatebg", vm=VirtualMachine.SERVER)
public class EvaluateMosaicClasses {
    private Reader seqs;
    private Reader trainSeqs;
    private Reader testSeqs;
    private int minClasses = 1;
    private int maxClasses = 8;
    private int mosaicOrder = 1;
    private double mosaicTransition = 0.005;
    private int evalSeqs = 50;
    private int trim = 0;
    private FiniteAlphabet alpha = DNATools.getDNA();
    
    @Option(help="Alphabet of the input sequences (default=DNA)", optional=true)
    public void setAlphabet(MicaAlphabet alpha) 
    {
    	this.alpha = alpha.biojavaAlphabet();
    }
    
    @Option(help="Fasta file of sequences to use for the training phases", optional=true)
    public void setTrainSeqs(Reader r) {
        this.trainSeqs = r;
    }

    @Option(help="Fasta file of sequences to use for the evalutation phases", optional=true)
    public void setTestSeqs(Reader r) {
        this.testSeqs = r;
    }
    
    @Option(help="Number of sequences to use for evalutation", optional=true)
    public void setEvalSeqs(int i) {
        this.evalSeqs = i;
    }
    
    @Option(help="Fasta file of sequences to split between testing and training", optional=true)
    public void setSeqs(Reader r) {
        this.seqs = r;
    }
    
    @Option(help="The minimum number of classes to evaluate", optional=true)
    public void setMinClasses(int i) {
        this.minClasses = i;
    }
    
    @Option(help="The maximum number of classes to evaluate", optional=true)
    public void setMaxClasses(int i) {
        this.maxClasses = i;
    }
    
    @Option(help="Markov chain order (zero-based)", optional=true)
    public void setOrder(int i) {
        this.mosaicOrder = i + 1;
    }
    
    @Option(help="Mosaic transition probability", optional=true)
    public void setMosaicTransition(double d) {
        this.mosaicTransition = d;
    }

    @Option(help="Number of bases to trim from the start of sequences (useful when comparing multiple -mosaicOrder values)", optional=true)
    public void setTrim(int i) {
        this.trim = i;
    }
    
    public SequenceDB makeDB(Collection<Sequence> c)
    		throws Exception
    {
        SequenceDB db = new HashSequenceDB();
        for (Sequence s: c) {
            db.addSequence(s);
        }
        return db;
    }
    
    private SequenceDB loadDB(Reader r)
    		throws Exception
    {
        SequenceDB db = new HashSequenceDB();
        SequenceIterator si = SeqIOTools.readFasta(new BufferedReader(r), alpha.getTokenization("token"));
        while (si.hasNext()) {
            db.addSequence(si.nextSequence());
        }
        return db;
    }
    
    public void main(String[] args) 
        throws Exception
    {
        SequenceDB testDB, trainDB;
        if (testSeqs != null) {
            testDB = loadDB(testSeqs);
            trainDB = loadDB(trainSeqs);
        } else if (seqs != null) {
        	SequenceDB allDB = loadDB(seqs);
	        List<Sequence> seqList = new ArrayList<Sequence>();
	        {
	            SequenceIterator si = allDB.sequenceIterator();
	            while (si.hasNext()) {
	            	seqList.add(si.nextSequence());
	            }
	        }
	        Collections.shuffle(seqList);
	        testDB = makeDB(seqList.subList(0, evalSeqs));
	        trainDB = makeDB(seqList.subList(evalSeqs, seqList.size()));
        } else {
        	System.err.println("You must specify either the -seqs option or both -trainSeqs and -testSeqs");
        	return;
        }
        
        for (int mosaicClasses = minClasses; mosaicClasses <= maxClasses; ++mosaicClasses) {
	        Distribution[] patches = MosaicTools.optimizePatches(
	            trainDB,
	            mosaicClasses,
	            mosaicOrder, 
	            mosaicTransition,
	            false
	        );
	        
	        MarkovModel mm = MosaicTools.makePatchModel(patches, mosaicTransition);
            DP dp = new org.biojava.bio.dp.onehead.SingleDP(mm);
            double score = 0;
            for (SequenceIterator si = testDB.sequenceIterator(); si.hasNext(); ) {
                Sequence seq = si.nextSequence();
                SymbolList sl = seq.subList(trim + 1, seq.length());
                if (mosaicOrder > 1) {
                    sl = SymbolListViews.orderNSymbolList(sl, mosaicOrder);
                }
                score += dp.forward(new SymbolList[] {sl}, ScoreType.PROBABILITY);
            }
            System.out.println("" + mosaicClasses + "\t" + score);
        }
    }
}


