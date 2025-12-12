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

package net.derkholm.nmica.seq.align;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.GapDistribution;
import org.biojava.bio.dist.PairDistribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.DotState;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleDotState;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.dp.twohead.CellCalculatorFactoryMaker;
import org.biojava.bio.dp.twohead.DPInterpreter;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;

/**
 * @author thomas
 */
public class OverhangAligner implements Aligner {
    private final FiniteAlphabet alpha = DNATools.getDNA();
    private final DP aligner;
    
    public OverhangAligner() 
    		throws Exception
    {
        // CellCalculatorFactoryMaker cfFactM = new DPCompiler(false);
        CellCalculatorFactoryMaker cfFactM = new DPInterpreter.Maker();
        DPFactory fact = new DPFactory.DefaultFactory(cfFactM);
        
        MarkovModel model = generateAligner(
          alpha,
          0.79, 0.98,
          0.1, 0.95
        );
        
        aligner = fact.createDP(model);
    }
    
    public SymbolList align(SymbolList seq0, SymbolList seq1)
    		throws Exception
    	{
        StatePath result = aligner.viterbi(new SymbolList[] {seq0, seq1}, ScoreType.PROBABILITY);
        SymbolList rawAlign = result.symbolListForLabel(StatePath.SEQUENCE);
        return rawAlign;
    	}
    
    private static MarkovModel generateAligner(
            FiniteAlphabet alpha,
            double pMatch, double pExtendMatch,
            double pGap, double pExtendGap
    ) throws Exception {
        double pEndMatch = 1.0 - pExtendMatch;
        double pEndGap = 1.0 - pExtendGap;
        double pEnd = 1.0 - pMatch - 2.0*pGap; 
        
        FiniteAlphabet dna = alpha;
        FiniteAlphabet dna2 =
            (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(
                    Collections.nCopies(2, dna)
            );
        
        MarkovModel model = new SimpleMarkovModel(2, dna2, "pair-wise aligner with overhangs");
        
        Distribution nullModel = new UniformDistribution(dna);
        Distribution gap = new GapDistribution(dna);
        Distribution matchDist = generateMatchDist((FiniteAlphabet) dna2);
        Distribution nullModel2 = new PairDistribution(nullModel, nullModel);
        Distribution insert1Dist = new PairDistribution(nullModel, gap);
        Distribution insert2Dist = new PairDistribution(gap, nullModel);
        
        DotState hub = new SimpleDotState("hub");
        EmissionState match = new SimpleEmissionState(
                "match",
                Annotation.EMPTY_ANNOTATION,
                new int [] { 1, 1 },
                matchDist
        );
        
        EmissionState insert1 = new SimpleEmissionState(
                "insert1",
                Annotation.EMPTY_ANNOTATION,
                new int [] { 1, 0 },
                insert1Dist
        );
        
        EmissionState insert2 = new SimpleEmissionState(
                "insert2",
                Annotation.EMPTY_ANNOTATION,
                new int [] { 0, 1 },
                insert2Dist
        );
        
        EmissionState preOverhang1 = new SimpleEmissionState(
                "preoverhang1",
                Annotation.EMPTY_ANNOTATION,
                new int[] {1, 0},
                insert1Dist
        );
        EmissionState preOverhang2 = new SimpleEmissionState(
                "preoverhang2",
                Annotation.EMPTY_ANNOTATION,
                new int[] {0, 1},
                insert2Dist
        );
        EmissionState postOverhang1 = new SimpleEmissionState(
                "postoverhang1",
                Annotation.EMPTY_ANNOTATION,
                new int[] {1, 0},
                insert1Dist
        );
        EmissionState postOverhang2 = new SimpleEmissionState(
                "postoverhang2",
                Annotation.EMPTY_ANNOTATION,
                new int[] {0, 1},
                insert2Dist
        );
        
        model.addState(hub);
        model.addState(match);
        model.addState(insert1);
        model.addState(insert2);
        model.addState(preOverhang1);
        model.addState(preOverhang2);
        model.addState(postOverhang1);
        model.addState(postOverhang2);
        
        Distribution dist;
        
        model.createTransition(model.magicalState(), hub);
        model.createTransition(model.magicalState(), preOverhang1);
        model.createTransition(model.magicalState(), preOverhang2);
        
        model.createTransition(preOverhang1, preOverhang1);
        model.createTransition(preOverhang2, preOverhang2);
        model.createTransition(preOverhang1, hub);
        model.createTransition(preOverhang2, hub);
        
        model.createTransition(hub, match);
        model.createTransition(hub, insert1);
        model.createTransition(hub, insert2);
        
        model.createTransition(hub, postOverhang1);
        model.createTransition(hub, postOverhang2);
        model.createTransition(hub, model.magicalState());
        
        model.createTransition(postOverhang1, postOverhang1);
        model.createTransition(postOverhang2, postOverhang2);
        model.createTransition(postOverhang1, model.magicalState());
        model.createTransition(postOverhang2, model.magicalState());
        
        model.createTransition(match, match);
        model.createTransition(match, hub);
        
        model.createTransition(insert1, insert1);
        model.createTransition(insert1, hub);
        
        model.createTransition(insert2, insert2);
        model.createTransition(insert2, hub);
        
        model.getWeights(model.magicalState()).setWeight(hub, 1.0 / 3);
        model.getWeights(model.magicalState()).setWeight(preOverhang1, 1.0 / 3);
        model.getWeights(model.magicalState()).setWeight(preOverhang2, 1.0 / 3);
        
        {
            dist = model.getWeights(preOverhang1);
            dist.setWeight(preOverhang1, 0.99);
            dist.setWeight(hub, 0.01);
        }
        {
            dist = model.getWeights(preOverhang2);
            dist.setWeight(preOverhang2, 0.99);
            dist.setWeight(hub, 0.01);
        }
        
        dist = model.getWeights(hub);
        dist.setWeight(match, pMatch);
        dist.setWeight(insert1, pGap);
        dist.setWeight(insert2, pGap);
        dist.setWeight(model.magicalState(), pEnd / 3);
        dist.setWeight(postOverhang1, pEnd / 3);
        dist.setWeight(postOverhang2, pEnd / 3);
        
        {
            dist = model.getWeights(postOverhang1);
            dist.setWeight(postOverhang1, 0.99);
            dist.setWeight(model.magicalState(), 0.01);
        }
        {
            dist = model.getWeights(postOverhang2);
            dist.setWeight(postOverhang2, 0.99);
            dist.setWeight(model.magicalState(), 0.01);
        }
        
        dist = model.getWeights(match);
        dist.setWeight(match, pExtendMatch);
        dist.setWeight(hub, pEndMatch);
        
        dist = model.getWeights(insert1);    
        dist.setWeight(insert1, pExtendGap);
        dist.setWeight(hub, pEndGap);
        
        dist = model.getWeights(insert2);    
        dist.setWeight(insert2, pExtendGap);
        dist.setWeight(hub, pEndGap);
        
        return model;
    }
    
    private static Distribution generateMatchDist(FiniteAlphabet dna2)
    throws Exception {
        Distribution dist = DistributionFactory.DEFAULT.createDistribution(dna2);
        int size = dna2.size();
        int matches = (int) Math.sqrt(size);
        
        double pMatch = 0.45;
        
        double matchWeight = pMatch / matches;
        double missWeigth = (1.0 - pMatch) / (size - matches);
        
        for(Iterator i = dna2.iterator(); i.hasNext(); ) {
            BasisSymbol cps = (BasisSymbol) i.next();
            List sl = cps.getSymbols();
            if(sl.get(0) == sl.get(1)) {
                dist.setWeight(cps, matchWeight);
            } else {
                dist.setWeight(cps, missWeigth);
            }
        }
        
        return dist;
    }
}
