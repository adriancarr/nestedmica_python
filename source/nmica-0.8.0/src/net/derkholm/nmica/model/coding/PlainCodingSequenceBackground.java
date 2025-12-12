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

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.model.motif.SequenceBackground;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.dist.OrderNDistributionFactory;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ListTools;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Interface for basic coding sequence models.
 *
 * @author Bernard Leong
 */

public class PlainCodingSequenceBackground implements SequenceBackground, Serializable {

    static final long serialVersionUID = -742638802183258198L;
    
    private final static Pattern namePattern = Pattern.compile("patch([0-9]+)");
    private transient OrderNDistribution plainBackgroundCodingDistribution;

    public PlainCodingSequenceBackground(OrderNDistribution plainBackgroundCodingDistribution) 
    {
        this.plainBackgroundCodingDistribution = plainBackgroundCodingDistribution;
    }
    
    public OrderNDistribution getPlainBackgroundCodingDistribution() {
        return plainBackgroundCodingDistribution;
    }
    
    /**
     * returns the background coding distribution
     */
    
    public Location backgroundSymbolLikelihood(SymbolList sl, double[] likelihoods)
	throws IllegalAlphabetException
    {
	try{
	    	   
	    FiniteAlphabet alpha = ProteinTools.getTAlphabet();
	    FiniteAlphabet alphanew = Coding.createProteinAlphabet(alpha);

	    Alphabet codonAlphabet = AlphabetManager.getCrossProductAlphabet(
				     Collections.nCopies(3, DNATools.getDNA())
	                             ); //conditioning alphabet



	    Alphabet conditionalAlphabet = AlphabetManager.getCrossProductAlphabet(
					   new ListTools.Doublet(
                                           DNATools.getDNA(),
                                           codonAlphabet
		                           )
	                                   );
 
	    SymbolList seq = sl;	  
	    Sequence cloneSeq = (Sequence) seq;
	    
	    int startPhase = Coding.getStartPhase(cloneSeq);
        int firstPos = 4 - startPhase;
        int lastPos = seq.length() - 1;
            
	    while (((lastPos - firstPos + 1) % 3) != 0) {
                --lastPos;
            }
	    
	    SymbolList codonSL = SymbolListViews.windowedSymbolList(
				 seq.subList(firstPos, lastPos), 3);
	    
	    SymbolList aminocodonSL = SymbolListViews.windowedSymbolList(DNATools.toRNA(seq.subList(firstPos, lastPos)), 3);

	    SymbolList aminoacid = RNATools.translate(aminocodonSL);

	    OrderNDistribution codonBackground = getPlainBackgroundCodingDistribution();

	    List<Location> locList = new ArrayList<Location>();

	    for (int pos = 2; pos < codonSL.length(); ++pos) 
	    	{
	    		BasisSymbol codon = (BasisSymbol) codonSL.symbolAt(pos);
	    		BasisSymbol oldcodon = (BasisSymbol) codonSL.symbolAt(pos-1);
	    		List<Symbol> symList = codon.getSymbols();
	    		Symbol lastbase = Coding.extractLastBase(oldcodon);
	    		
	    		
	    		
			/*if(amino == ProteinTools.ter()){
			    throw new Exception("Hit nonsense codon"); //Paranoia check
			}*/

	    		int basePos = ((pos - 1) * 3) + 1;

	    		for (int c = 0; c < 3; ++c) {
	    			List<Symbol> newSymList = new ArrayList<Symbol>(Collections.nCopies(3, DNATools.n()));
	    			newSymList.set(c, symList.get(c));
	    			Symbol maskedCodon = codonAlphabet.getSymbol(newSymList);
	    			Symbol conditioner = conditionalAlphabet.getSymbol(new ListTools.Doublet(lastbase, maskedCodon));
	    			double sh = 0;
	    			sh += codonBackground.getWeight(conditioner);
	    			likelihoods[basePos + c] = NativeMath.log2(sh);
	    			locList.add(new PointLocation(basePos + c));
	    			
	        }
	    
	    }

	    
	    return LocationTools.union(locList);

	} catch (Exception ex) {
            throw new IllegalAlphabetException(ex, "Couldn't evaluated background");
         }
        
	
	
    }
    
        
    
    private void writeObject(ObjectOutputStream oos)
    throws IOException {

        oos.defaultWriteObject();
        OrderNDistribution plainCodonBackground = getPlainBackgroundCodingDistribution();
        
        try{
	        Symbol[] orderDNA = new Symbol[] {DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t()};
	        Alphabet codonAlphabet = AlphabetManager.getCrossProductAlphabet(
			        Collections.nCopies(3, DNATools.getDNA())
			); //conditioning alphabet
	
	        Alphabet conditionalAlphabet = AlphabetManager.getCrossProductAlphabet(
	                new ListTools.Doublet(
	                        DNATools.getDNA(),
	                        codonAlphabet
	                )
			); //conditional alphabet
	        //List dummyList = new ArrayList();
	        
	        for(int i = 0; i < orderDNA.length; ++i){
	            for(int j = 0; j < orderDNA.length; ++j) {
	                for(int k = 0; k < orderDNA.length; k++){
	                    for(int l = 0; l < orderDNA.length; l++){
	
	                        Symbol simpleCodon = codonAlphabet.getSymbol(Arrays.
	                                asList(new Object[] {(Symbol) orderDNA[i], 
	                                        			  (Symbol) orderDNA[j], 
	                                        			  (Symbol) orderDNA[k]})); 
	                        
	                        Symbol combinedCodon = conditionalAlphabet.getSymbol(Arrays.asList(
	                                new Object[]{(Symbol) orderDNA[l], (Symbol) simpleCodon}));
	        
	                        //dummyList.add((Symbol) combinedCodon);
	                        double codonBackgroundScore = (double) plainCodonBackground.getWeight(combinedCodon);
	                        System.out.println(combinedCodon.getName() + "\t" + plainCodonBackground.getWeight(combinedCodon));
	                        oos.writeDouble(codonBackgroundScore);
	                        
	                        
	                        
	                    }
	                    
	                }
	            }
	        }
	                
        } catch(Exception e){
          throw new IOException("Couldn't serialize biojava crap");
        }
     
    }
        
    private void readObject(ObjectInputStream ois)
		throws IOException, ClassNotFoundException, IllegalAlphabetException, IllegalSymbolException, ChangeVetoException
	{       
        ois.defaultReadObject();
        DistributionFactory conditionalFactory  = new OrderNDistributionFactory(DistributionFactory.DEFAULT);
        Symbol[] orderDNA = new Symbol[] {DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t()};
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
                
        for(int i = 0; i < orderDNA.length; ++i){
            for(int j = 0; j < orderDNA.length; ++j) {
                for(int k = 0; k < orderDNA.length; k++){
                    for(int l = 0; l < orderDNA.length; l++){

                        Symbol simpleCodon = codonAlphabet.getSymbol(Arrays.
                                asList(new Object[] {(Symbol) orderDNA[i], 
                                        			  (Symbol) orderDNA[j], 
                                        			  (Symbol) orderDNA[k]})); 
                        
                        Symbol combinedCodon = conditionalAlphabet.getSymbol(Arrays.asList(
                                new Object[]{(Symbol) orderDNA[l], (Symbol) simpleCodon}));
                        
                        double plainCodonScore = (double) ois.readDouble();
                        conditionalCodons.setWeight(combinedCodon, plainCodonScore);
                        
                        //System.err.println(combinedCodon.getName() + (double) ois.readDouble());
                        //modelToLikelihood.put(combinedCodon, new Double(codonBackgroundscore[])) 
                    }
                }
            }
        }
       
        plainBackgroundCodingDistribution = conditionalCodons;   
        
    }
    

}

