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

import java.util.Iterator;

import net.derkholm.nmica.seq.NMSimpleDistribution;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Count;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.dist.OrderNDistributionFactory;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.dist.UntrainableDistribution;
import org.biojava.bio.dp.BaumWelchTrainer;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.StoppingCriteria;
import org.biojava.bio.dp.TrainingAlgorithm;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolListViews;

/**
 * Support class for the mosaic family of background models.  Probably
 * needs refactoring
 *
 * @author Thomas Down
 * @since mica2
 */

public class MosaicTools {
    private final MarkovModel mm;
    private final Distribution[] readouts;
    
    private MosaicTools(MarkovModel mm, Distribution[] readouts) {
        this.mm = mm;
        this.readouts = readouts;
    }
   
    
    public static Distribution[] optimizePatches(SequenceDB seqs, int patches, int word, double transToAlt, boolean match)
        throws Exception
    {
        if (word > 1) {
            SequenceDB newSeqs = new HashSequenceDB();
            for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
                Sequence s = si.nextSequence();
                Sequence s2 = new SimpleSequence(
                    SymbolListViews.orderNSymbolList(s, word),
                    null,
                    s.getName() + "-o" + word,
                    Annotation.EMPTY_ANNOTATION
                );
                newSeqs.addSequence(s2);
            }
            seqs = newSeqs;
        }
        
        MosaicTools patchy = makePatchModel((FiniteAlphabet) seqs.sequenceIterator().nextSequence().getAlphabet(), patches, transToAlt, seqs, match);
        DP dp = new org.biojava.bio.dp.onehead.SingleDP(patchy.mm);
        TrainingAlgorithm ta = new BaumWelchTrainer(dp);
        
        ta.train(
            seqs,
            0.01,
            new StoppingCriteria() {
                public boolean isTrainingComplete(TrainingAlgorithm ta) {
                    // System.err.println(ta.getCycle() + "\t" + ta.getCurrentScore());
                    System.err.print('.');
                    return Math.abs(ta.getLastScore() - ta.getCurrentScore()) < 0.001;
                }
            }
        );
        System.err.println();
        return patchy.readouts;
    }
    
    private final static int[] ADVANCE_ONE = new int[] {1};
    
    public static MarkovModel makePatchModel(Distribution[] dists, double transToAlt)
        throws Exception
    {
        int patches = dists.length;
        FiniteAlphabet alpha = (FiniteAlphabet) dists[0].getAlphabet();
        
        // double transToAlt = 0.005;
        double transToEnd = 0.001;
        double transToSelf = 1.0 - transToEnd - (patches - 1) * transToAlt;
        
        // System.err.println("transToSelf = " + transToSelf);
        
        MarkovModel mm = new SimpleMarkovModel(
            1,
            alpha,
            "Patchy composition model"
        );
        State magic = mm.magicalState();
        
        EmissionState[] patchStates = new EmissionState[patches];
        for (int p = 0; p < patches; ++p) {
            patchStates[p] = new SimpleEmissionState(
                "patch" + p,
                Annotation.EMPTY_ANNOTATION,
                ADVANCE_ONE,
                dists[p]
            );
            mm.addState(patchStates[p]);
        }
        
        for (int p = 0; p < patches; ++p) {
            mm.createTransition(magic, patchStates[p]);
            mm.createTransition(patchStates[p], magic);
            for (int q = 0; q < patches; ++q) {
                mm.createTransition(patchStates[p], patchStates[q]);
            }
        }
        
        {
            Distribution dist = new UntrainableDistribution(mm.transitionsFrom(magic));
            for (int p = 0; p < patches; ++p) {
                dist.setWeight(patchStates[p], 1.0 / patches);
            }
            mm.setWeights(magic, dist);
        }
        
        for (int p = 0; p < patches; ++p) {
            Distribution dist = new UntrainableDistribution(mm.transitionsFrom(patchStates[p]));
            for (int q = 0; q < patches; ++q) {
                double pr;
                if (p == q) {
                    pr = transToSelf;
                } else {
                    pr = transToAlt;
                }
                dist.setWeight(patchStates[q], pr);
            }
            dist.setWeight(magic, transToEnd);
            mm.setWeights(patchStates[p], dist);
        }
        
        return mm;
    }
    
    private static MosaicTools makePatchModel(FiniteAlphabet alpha, int patches, double transToAlt, SequenceDB seqs, boolean matchDists)
        throws Exception
    {
        Distribution[] dists = new Distribution[patches];
        for (int p = 0; p < patches; ++p) {
        	if (matchDists) {
        		dists[p] = matchedRandomDist(alpha, seqs);
        	} else {
        		dists[p] = randomDist(alpha);
        	}
        }

        return new MosaicTools(makePatchModel(dists, transToAlt), dists);
    }
    
    private static Distribution randomDist(FiniteAlphabet fa)
	    throws Exception
	{
	    int samples = 10;
	    int wordSize = fa.getAlphabets().size();
	    Distribution dist;
	    if (wordSize == 1) {
	        dist = new NMSimpleDistribution(fa);
	    } else {
	        dist = new OrderNDistributionFactory(DistributionFactory.DEFAULT).createDistribution(fa);
	    }
	    DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
	    dtc.registerDistribution(dist);
	    for (Iterator i = fa.iterator(); i.hasNext(); ) {
	        Symbol s = (Symbol) i.next();
	        for (int j = 0; j < samples; ++j) {
	            dtc.addCount(dist, s, Math.random());
	        }
	    }
	    dtc.train();
	    return dist;
	}

    
    private static Distribution matchedRandomDist(FiniteAlphabet fa,SequenceDB seqs)
        throws Exception
    {
        int samples = 5;
        int wordSize = fa.getAlphabets().size();
        Distribution dist;
        if (wordSize == 1) {
            dist = new NMSimpleDistribution(fa);
        } else {
            dist = new OrderNDistributionFactory(DistributionFactory.DEFAULT).createDistribution(fa);
        }
        DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
        dtc.registerDistribution(dist);
// Mutlu -->

	Count count = new IndexedCount(fa);
	int allCounts = 0;
	for (SequenceIterator seqit = seqs.sequenceIterator();seqit.hasNext();) {
		Sequence seq = seqit.nextSequence();
		for (int i =1;i<seq.length();i++) {
			Symbol s = seq.symbolAt(i);
			if (s instanceof AtomicSymbol) {
				count.increaseCount((AtomicSymbol)s,1.0);
				allCounts++;
			}
		}
	}	
        for (Iterator i = fa.iterator(); i.hasNext(); ) {
            Symbol s = (Symbol) i.next();
            for (int j = 0; j < samples; ++j) {
                dtc.addCount(dist, s, Math.random());
            }
            dtc.addCount(dist, s, 3.0*count.getCount((AtomicSymbol) s) / allCounts * 1.0);
        }
// <-- Mutlu
        dtc.train();
        return dist;
    }

}
