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

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.dp.twohead.CellCalculatorFactoryMaker;
import org.biojava.bio.dp.twohead.DPCompiler;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ListTools;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
/**
 * @author mutlu
 */
public class DiNucleotideAligner implements Aligner {
    private final FiniteAlphabet alpha = DNATools.getDNA();
    private final FiniteAlphabet diNucleotide = (FiniteAlphabet)AlphabetManager.getCrossProductAlphabet(Collections.nCopies(2, DNATools.getDNA()));
    private final Alphabet quartetNucleotide = AlphabetManager.getCrossProductAlphabet(Collections.nCopies(2, diNucleotide));
    private final DP aligner;

    public DiNucleotideAligner() 
    		throws Exception
    {
        CellCalculatorFactoryMaker cfFactM = new DPCompiler(false);
        DPFactory fact = new DPFactory.DefaultFactory(cfFactM);

        MarkovModel model = generateAligner(
          diNucleotide,
          0.79, 0.98,
          0.1, 0.95
        );

        aligner = fact.createDP(model);
    }

    public SymbolList align(SymbolList seq0, SymbolList seq1)
    		throws Exception
    	{
        SymbolList viewSourceSeq = SymbolListViews.orderNSymbolList(seq0, 2);
        SymbolList viewTargetSeq = SymbolListViews.orderNSymbolList(seq1, 2);
        StatePath result = aligner.viterbi(new SymbolList[] {viewSourceSeq, viewTargetSeq}, ScoreType.ODDS);
        SymbolList rawAlign = result.symbolListForLabel(StatePath.SEQUENCE);
        List<Symbol> pruneDiGapsList = new ArrayList<Symbol>();

        for (int i = 1; i < rawAlign.length(); i++) {
            Symbol x =  rawAlign.symbolAt (i);
             if (x instanceof BasisSymbol) {
                BasisSymbol xBasis= (BasisSymbol) x;
                Symbol topComponent    = (Symbol) xBasis.getSymbols().get(0);
                Symbol bottomComponent = (Symbol) xBasis.getSymbols().get(1);

                if ((topComponent instanceof BasisSymbol)&&(bottomComponent instanceof BasisSymbol)) {
                   Symbol alignTopSym    = (Symbol) ((BasisSymbol) topComponent   ).getSymbols().get(1);
                   Symbol alignBottomSym = (Symbol) ((BasisSymbol) bottomComponent).getSymbols().get(1);

                   if ((alignTopSym.getName().equals("[]")) && ((alignBottomSym.getName().equals("[]")))) {
                      continue;
                   }
                   pruneDiGapsList.add(x);
                   if (i==1) {pruneDiGapsList.add(x);}
                }
            }
        }
        SymbolList niceAlign = new SimpleSymbolList(quartetNucleotide, pruneDiGapsList); 
        List<Symbol> mapToMonoAlignment = new ArrayList<Symbol>();

        for (int i = 1; i < niceAlign.length(); i++) {
            Symbol x =  niceAlign.symbolAt (i);
            if (x instanceof BasisSymbol) {
                BasisSymbol xBasis= (BasisSymbol) x;
                Symbol topComponent    = (Symbol) xBasis.getSymbols().get(0);
                Symbol bottomComponent = (Symbol) xBasis.getSymbols().get(1);

                if (topComponent instanceof BasisSymbol) {
                     Symbol alignTopSym    = (Symbol) ((BasisSymbol) topComponent   ).getSymbols().get(1);
                     Symbol alignBottomSym = (Symbol) ((BasisSymbol) bottomComponent).getSymbols().get(1);
                     if ((alignTopSym.getName().equals("[]")) && ((alignBottomSym.getName().equals("[]")))) continue;
                     Symbol aligned = diNucleotide.getSymbol (new ListTools.Doublet (alignTopSym, alignBottomSym));
                     mapToMonoAlignment.add(aligned);
                }
                else {
                     continue;
                }
            }
            else {
                continue;
            }
        }

        Alphabet alignAlpha = rawAlign.getAlphabet();
        SymbolList monoAlign = new SimpleSymbolList(alignAlpha, mapToMonoAlignment); 

        return monoAlign;
    }

    private static MarkovModel generateAligner(
            FiniteAlphabet diNucleotide,
            double pMatch, double pExtendMatch,
            double pGap, double pExtendGap
    ) throws Exception {
        double pEndMatch = 1.0 - pExtendMatch;
        double pEndGap = 1.0 - pExtendGap;
        double pEnd = 1.0 - pMatch - 2.0*pGap; 


    FiniteAlphabet dna2 = diNucleotide;
    FiniteAlphabet dna2x2 =
      (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(
        Collections.nCopies(2, dna2)
      );


        MarkovModel model = new SimpleMarkovModel(2, dna2x2, "pair-wise aligner");

        Distribution nullModel = new UniformDistribution(dna2);
        Distribution gap = new GapDistribution(dna2);
        Distribution matchDist = generateMatchDist((FiniteAlphabet) dna2x2);
        Distribution nullModel2 = new PairDistribution(nullModel, nullModel);
        Distribution gap2 = new PairDistribution(gap, gap);
        Distribution insert1Dist = new PairDistribution(nullModel2, gap2);
        Distribution insert2Dist = new PairDistribution(gap2, nullModel2);

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
    FileInputStream disk  = new FileInputStream ("SavedDistributionObject.dat");
    ObjectInputStream obj = new ObjectInputStream (disk);
    dist=(OrderNDistribution) obj.readObject();

    return dist;
  }
}
