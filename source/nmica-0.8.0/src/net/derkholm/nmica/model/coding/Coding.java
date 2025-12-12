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

package net.derkholm.nmica.model.coding;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ListTools;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;
/**
 * 
 * Support class for the coding family of background models.  Probably
 * needs refactoring. A background model to model internal coding exons 
 * based on distribution of amino acids and (DNA * (DNA)^3) codon usage.
 * nmica.model.coding
 * @author Bernard Leong
 * @since mica2
 */

public class Coding {

    /**
     * Calculate Codon distributions (conditioned by previous amino acid and last base of the
     * previous amino acid for an input SequenceDB
     * @param seqs
     * @return OrderNDistribution
     * @throws Exception
     */
    public static OrderNDistribution makeCodonDistributions(SequenceDB seqs)
		throws Exception	
    {
        FiniteAlphabet alpha = ProteinTools.getTAlphabet();
        FiniteAlphabet alphanew = createProteinAlphabet(alpha);
		List<Symbol> am = new ArrayList<Symbol>();

		for (Iterator<Symbol> si = alphanew.iterator(); si.hasNext(); ){
		    am.add(si.next());
		}
		int patches = am.size();
	
		OrderNDistribution codons = createCodonDist(seqs, alphanew);
		
		for(Iterator<Symbol> gi = ((FiniteAlphabet) codons.getAlphabet()).iterator(); gi.hasNext(); ){
		    Symbol a = gi.next();
		    
		    if(Double.isNaN(codons.getWeight(a))){
			codons.setWeight(a, 0.0);
		    }
	
		    //System.out.println(a.getName() + "\t" + codons.getWeight(a)); //check weights
		}
	
		return codons;
    }

    
    public static OrderNDistribution makePlainCodonDistributions(SequenceDB seqs)
		throws Exception	
	{
        //System.err.println("Plain Codon");
	    FiniteAlphabet alpha = ProteinTools.getTAlphabet();
	    FiniteAlphabet alphanew = createProteinAlphabet(alpha);
		List<Symbol> am = new ArrayList<Symbol>();
	
		for (Iterator<Symbol> si = alphanew.iterator(); si.hasNext(); ){
		    am.add(si.next());
		}
		int patches = am.size();
	
		OrderNDistribution codons = createPlainCodonDist(seqs, alphanew);
		
		for(Iterator<Symbol> gi = ((FiniteAlphabet) codons.getAlphabet()).iterator(); gi.hasNext(); ){
		    Symbol a = gi.next();
		    
		    if(Double.isNaN(codons.getWeight(a))){
			codons.setWeight(a, 0.0);
		    }
	
		    //System.out.println(a.getName() + "\t" + codons.getWeight(a)); //check weights
		}
	
		return codons;
	}
    
    /**
     * Calculates an OrderNDistribution distribution of the probability of finding a codon given 
     * the previous amino acid and the last base of the previous codon.
     * 
     * @param seqs the input coding SequenceDB
     * @param alphanew the proteinT alphabet
     * @return OrderNDistribution codons
     * @throws Exception
     */
    
    private static OrderNDistribution createCodonDist(SequenceDB seqs, FiniteAlphabet alphanew)
		throws Exception
    {
		DistributionFactory conditionalFactory 
		    = new OrderNDistributionFactory(DistributionFactory.DEFAULT);
		
		Alphabet codonAlphabet = AlphabetManager.getCrossProductAlphabet(
	            Collections.nCopies(3, DNATools.getDNA())
		    ); //conditioning alphabet
	
		Alphabet conditioningAlphabet =  AlphabetManager.getCrossProductAlphabet(
		    new ListTools.Doublet(alphanew,
					  DNATools.getDNA()
	        ));
	
		Alphabet conditionalAlphabet = AlphabetManager.getCrossProductAlphabet(
	            new ListTools.Doublet(
	                conditioningAlphabet,
	                codonAlphabet
			)
		    ); //conditional alphabet
		
	        OrderNDistribution conditionalCodons 
		    = (OrderNDistribution) conditionalFactory.createDistribution(
	            conditionalAlphabet
		    );
	        
	        DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
	        dtc.registerDistribution(conditionalCodons);
	        for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
	            Sequence seq = si.nextSequence();
	            
		    int startPhase = getStartPhase(seq);
	            int firstPos = 4 - startPhase;
	            int lastPos = seq.length() - 1;
	            while (((lastPos - firstPos + 1) % 3) != 0) {
	                --lastPos;
	            }
	        
	         
		    //System.err.println("Start phase =" + startPhase + "\t" + "firstPos=" + firstPos + "\t lastPos=" + lastPos);
		    //System.err.println(seq.getName() + "\t" + seq.getAnnotation());
		    
	        SymbolList codonSL = SymbolListViews.windowedSymbolList(
					 seq.subList(firstPos, lastPos), 3);
		    
		    SymbolList aminocodonSL = SymbolListViews.windowedSymbolList(
		            						DNATools.toRNA(
		            						        seq.subList(firstPos, lastPos)), 3);
	
		    SymbolList aminoacid = RNATools.translate(aminocodonSL);
		    
		    for (int pos = 2; pos <= codonSL.length(); ++pos) {
			
			Symbol codon = codonSL.symbolAt(pos); //codon at position i
			BasisSymbol oldcodon = (BasisSymbol) codonSL.symbolAt(pos-1); //codon at position (i-1)
			Symbol lastbase = extractLastBase(oldcodon); //last base in old codon
			Symbol amino = aminoacid.symbolAt(pos-1); // translated old codon
			Symbol combinebases = (Symbol) conditioningAlphabet.getSymbol(new ListTools.Doublet(amino, lastbase));
			/*
			if(amino == ProteinTools.ter()){
			    throw new Exception("Hit nonsense codon"); //Paranoia check
			}
			*/
			/*System.err.println(codon.getName() + "\t" +
					   lastbase.getName() + "\t" +
					   oldcodon.getName() + "\t" +
					   amino.getName());*/
		     
			dtc.addCount(conditionalCodons,
				     conditionalAlphabet.getSymbol(
					   new ListTools.Doublet(
				       		 combinebases, 
						 codon
						 )
					   ),
				     1.0
				     );
	
			
			
		    }
		    dtc.train();
		}
		
	        return conditionalCodons;
		
    }
    
    /**
     * This method creates a OrderNDistribution with the conditional probability
     * P(current codon|last base from previous codon). 
     * @param seqs
     * @param alphanew
     * @return OrderNDistribution
     * @throws Exception
     */

    private static OrderNDistribution createPlainCodonDist(SequenceDB seqs, FiniteAlphabet alphanew)
    throws Exception
    {
    		DistributionFactory conditionalFactory = new OrderNDistributionFactory(DistributionFactory.DEFAULT);
	
    		Alphabet codonAlphabet = AlphabetManager.getCrossProductAlphabet(
    		        Collections.nCopies(3, DNATools.getDNA())
    		); //conditioning alphabet

    		Alphabet conditionalAlphabet = AlphabetManager.getCrossProductAlphabet(
            new ListTools.Doublet(
                DNATools.getDNA(),
                codonAlphabet
		)
	    ); //conditional alphabet
	
        OrderNDistribution conditionalCodons = (OrderNDistribution) conditionalFactory.createDistribution(conditionalAlphabet);
        
        DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
        dtc.registerDistribution(conditionalCodons);
        for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
            Sequence seq = si.nextSequence();
            
	    int startPhase = getStartPhase(seq);
            int firstPos = 4 - startPhase;
            int lastPos = seq.length() - 1;
            while (((lastPos - firstPos + 1) % 3) != 0) {
                --lastPos;
            }
        
         
	    //System.err.println("Start phase =" + startPhase + "\t" + "firstPos=" + firstPos + "\t lastPos=" + lastPos);
	    //System.err.println(seq.getName() + "\t" + seq.getAnnotation());
	    
        SymbolList codonSL = SymbolListViews.windowedSymbolList(
				 seq.subList(firstPos, lastPos), 3);
	    
	    SymbolList aminocodonSL = SymbolListViews.windowedSymbolList(
	            						DNATools.toRNA(
	            						        seq.subList(firstPos, lastPos)), 3);

	    SymbolList aminoacid = RNATools.translate(aminocodonSL);
	    
	    for (int pos = 2; pos <= codonSL.length(); ++pos) {
		
		Symbol codon = codonSL.symbolAt(pos); //codon at position i
		BasisSymbol oldcodon = (BasisSymbol) codonSL.symbolAt(pos-1); //codon at position (i-1)
		Symbol lastbase = extractLastBase(oldcodon); //last base in old codon
		
		
		/*
		if(amino == ProteinTools.ter()){
		    throw new Exception("Hit nonsense codon"); //Paranoia check
		}
		*/
		
		dtc.addCount(conditionalCodons,
			     conditionalAlphabet.getSymbol(
				   new ListTools.Doublet(
				     lastbase, 
					 codon
					 )
				   ),
			     1.0
			     );

		
		
	    }
	    dtc.train();
        }
	
        return conditionalCodons;
        
    }
    
    /**
     * Extract the last base from the previous codon
     * @param codon
     * @return Symbol
     */
    
    
    public static Symbol extractLastBase(BasisSymbol codon){
	
	List<Symbol> l = new ArrayList<Symbol>();
	for(Symbol s: (List<Symbol>) codon.getSymbols()) {
	    l.add(s);
	    
	}

        return l.get(2);

    }

    /**
     * Extracts the last base from the codon BasisSymbol object 
     * 
     */
    
    public static final Pattern PHASE_PATTERN;
    static {
        try {
            PHASE_PATTERN = Pattern.compile("start_phase=([012])");
        } catch (PatternSyntaxException pse) {
            throw new RuntimeException(pse);
        }
    }




    public static int getStartPhase(Sequence seq)
        throws Exception
    {
        Annotation ann = seq.getAnnotation();
        if (ann.containsProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE)) {
            String dline = (String) ann.getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE);
            Matcher pmatcher = PHASE_PATTERN.matcher(dline);
            if (pmatcher.find()) {
                // System.err.println("Found some phase information");
                return Integer.parseInt(pmatcher.group(1));
            }
        }
        return 0;
    }
    
    public static FiniteAlphabet createProteinAlphabet(FiniteAlphabet alpha)
	throws Exception
    {
 
      Set<Symbol> symbols = new HashSet<Symbol>();
      for (Iterator<Symbol> si = alpha.iterator(); si.hasNext(); ){
	    Symbol s = si.next();
	    if (s != ProteinTools.sec()){
		//Remove Selenocysteine from the alphabet list
		symbols.add(s);
		//System.err.println(s.getName());
	    }
	}

      //System.err.println();
      FiniteAlphabet alphanew = new SimpleAlphabet(symbols, "ProteinT20");
      AlphabetManager.registerAlphabet(alphanew.getName(), alphanew);
      return alphanew;
    }


}


