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

import net.derkholm.nmica.utils.CliTools;
import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.*;

import java.io.BufferedReader;
import java.io.Reader;
import java.util.*;

/**
 * @author thomas
 */
public class SpikeSequences {
    private String[] motifs;
    private int[] counts;
    private SequenceDB seqs;
    
    public void setMotif(String[] motifs) {
        this.motifs = motifs;
    }
    
    public void setCount(int[] counts) {
        this.counts = counts;
    }
    
    public void setSeqs(Reader r)
    		throws Exception
    {
        seqs = new HashSequenceDB();
        SequenceIterator si = SeqIOTools.readFastaProtein(new BufferedReader(r));
        while (si.hasNext()) {
            seqs.addSequence(si.nextSequence());
        }
    }
    
    public static void main(String[] args)
	    throws Exception
	{
        SpikeSequences app = new SpikeSequences();
        CliTools.configureBean(app, args);
        app.run();
	}
    
    public Set randomSubset(Set s, int count) {
        List l = new ArrayList(s);
        Collections.shuffle(l);
        return new HashSet(l.subList(0, count));
    }
	
    private int findSafePosition(int featureLength, int seqLength, Location mask) {
        while (true) {
            int pos = 1 + (int) Math.floor(Math.random() * (seqLength - featureLength));
            Location feature = new RangeLocation(pos, pos + featureLength - 1);
            if (!LocationTools.overlaps(feature, mask)) {
                return pos;
            }
        }
    }
    
    public void run()
    		throws Exception
    {
        SymbolList[] sMotifs = new SymbolList[motifs.length];
        Set[] idSets = new Set[motifs.length];
        for (int i = 0; i < motifs.length; ++i) {
            sMotifs[i] = ProteinTools.createProtein(motifs[i]);
            idSets[i] = randomSubset(seqs.ids(), counts[i]);
        }
        
	    for (SequenceIterator si = seqs.sequenceIterator(); si.hasNext(); ) {
	        Sequence seq = si.nextSequence();
	        Location mask = LocationTools.union(
	                new RangeLocation(1, 20),
	                new RangeLocation(seq.length() - 19, seq.length())
	        );
	        
	        SimpleSymbolList sl = new SimpleSymbolList(seq);
	        String descLine = seq.getName();
	        
	        for (int m = 0; m < sMotifs.length; ++m) {
	            if (idSets[m].contains(seq.getName())) {
	                int pos = findSafePosition(sMotifs[m].length(), seq.length(), mask);
	                
		            SymbolList atomicMotif = atomize(sMotifs[m]);
		            sl.edit(new Edit(pos, atomicMotif.length(), atomicMotif));
		            mask = LocationTools.union(
		                    mask,
		                    new RangeLocation(pos - 5, pos + atomicMotif.length() + 4)
		            );
		            
		            descLine += " motif" + m + "=" + pos;
	            }
	        }
	        
	        Annotation anno = new SmallAnnotation();
	        anno.setProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE, descLine);
	        new FastaFormat().writeSequence(
	            new SimpleSequence(
	                sl,
	                null,
	                seq.getName(),
	                anno
	            ), 
	            System.out
	        );
	    }
	}

	private static SymbolList atomize(SymbolList sl)
	    throws Exception
	{
	    List<Symbol> atoms = new ArrayList<Symbol>();
	    for (Iterator<Symbol> i = sl.iterator(); i.hasNext(); ) {
	        Symbol s = i.next();
	        Symbol[] options = alphaToArray((FiniteAlphabet) s.getMatches());
	        atoms.add(options[(int) Math.floor(Math.random() * options.length)]);
	    }
	    return new SimpleSymbolList(sl.getAlphabet(), atoms);
	}
	
	private static Symbol[] alphaToArray(FiniteAlphabet a)
	    throws Exception
	{
	    List l = new ArrayList();
	    for (Iterator ai = a.iterator(); ai.hasNext(); ) {
	        l.add(ai.next());
	    }
	    return (Symbol[]) l.toArray(new Symbol[0]);
	}
}
