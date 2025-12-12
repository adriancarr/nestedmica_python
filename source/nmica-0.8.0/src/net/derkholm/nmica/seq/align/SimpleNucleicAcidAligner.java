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
/*
 *                      NestedMICA
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 * 
 * Created on Jul 12, 2004
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
import org.biojava.bio.dp.twohead.DPCompiler;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;

/**
 * @author thomas
 */
public class SimpleNucleicAcidAligner implements Aligner {
    private final FiniteAlphabet alpha = DNATools.getDNA();
    private final DP aligner;
    
    public SimpleNucleicAcidAligner() 
    		throws Exception
    {
        CellCalculatorFactoryMaker cfFactM = new DPCompiler(false);
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
        
        MarkovModel model = new SimpleMarkovModel(2, dna2, "pair-wise aligner");
        
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
        
        model.addState(hub);
        model.addState(match);
        model.addState(insert1);
        model.addState(insert2);
        
        Distribution dist;
        
        model.createTransition(model.magicalState(), hub);
        
        model.createTransition(hub, match);
        model.createTransition(hub, insert1);
        model.createTransition(hub, insert2);
        model.createTransition(hub, model.magicalState());
        
        model.createTransition(match, match);
        model.createTransition(match, hub);
        
        model.createTransition(insert1, insert1);
        model.createTransition(insert1, hub);
        
        model.createTransition(insert2, insert2);
        model.createTransition(insert2, hub);
        
        model.getWeights(model.magicalState()).setWeight(hub, 1.0);
        
        dist = model.getWeights(hub);
        dist.setWeight(match, pMatch);
        dist.setWeight(insert1, pGap);
        dist.setWeight(insert2, pGap);
        dist.setWeight(model.magicalState(), pEnd);
        
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
        
        double pMatch = 0.6;
        
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
