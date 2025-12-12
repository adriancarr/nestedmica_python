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

import net.derkholm.nmica.seq.align.SimpleNucleicAcidAligner;
import net.derkholm.nmica.utils.CliTools;
import org.biojava.bio.Annotation;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.*;

import java.io.*;
import java.util.*;

/**
 * Align orthologous pairs of sequences, then shred the aligning regions
 * into convenient fragments.
 * 
 * @author Thomas Down
 */

public class ShredAlign {
    private int shredThreshold = 500;
    private int tileSize = 400;
    
    public void setShredThreshold(int i) {
        this.shredThreshold = i;
    }
    
    public void setTileSize(int i) {
        this.tileSize = i;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        ShredAlign app = new ShredAlign();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        if (args.length != 2) {
            System.err.println("Needs two sequence files");
        }
        
        SequenceDB seqs0 = loadDB(new File(args[0]));
        SequenceDB seqs1 = loadDB(new File(args[1]));
        PrintStream out0 = new PrintStream(new FileOutputStream("shred-" + args[0]));
        PrintStream out1 = new PrintStream(new FileOutputStream("shred-" + args[1]));
        
        SimpleNucleicAcidAligner ali = new SimpleNucleicAcidAligner();
        
        Set allIds = new TreeSet();
        allIds.addAll(seqs0.ids());
        allIds.addAll(seqs1.ids());
        for (Iterator idi = allIds.iterator(); idi.hasNext(); ) {
            String id = (String) idi.next();
            if (seqs0.ids().contains(id) && seqs1.ids().contains(id)) {
                Sequence seq0 = seqs0.getSequence(id);
                Sequence seq1 = seqs1.getSequence(id);
                if (seq0.length() > 3000 || seq1.length() > 3000) {
                    System.err.println("Skipping " + id + " because it's ridiculously long.  Sorry.");
                    continue;
                }
                
                System.err.println("Aligning " + id);
                SymbolList alignment = processAlignment(ali.align(seq0, seq1));
                
                if (alignment.length() < shredThreshold) {
                    dumpAlign(alignment, id, out0, out1);
                } else {
                    int numTiles = (int) Math.ceil((1.0 * alignment.length()) / tileSize);
                    int realTileSize = (int) Math.floor((1.0 * alignment.length()) / numTiles);
                    int pos = 1;
                    while (pos < alignment.length()) {
                        int endPos = pos + realTileSize - 1;
                        if (alignment.length() - endPos < 10) {
                            endPos = alignment.length();
                        }
                        dumpAlign(alignment.subList(pos, endPos), id + "__" + pos, out0, out1);
                        pos = endPos + 1;
                    }
                }
            } else {
                System.err.println("Skipping " + id + " because there isn't an ortholog");
            }
        }
    }
    
    private boolean isMatch(BasisSymbol sym)
		throws Exception
	{
		Symbol gap = ProteinTools.getAlphabet().getGapSymbol();
		for (Iterator si = sym.getSymbols().iterator(); si.hasNext(); ) {
		    if (!isRealProtein((Symbol) si.next())) {
		        return false;
		    }
		}
		return true;
	}
    
    private boolean isRealProtein(Symbol s) {
        if (s == ProteinTools.a()) {
            return true;
        } else if (s == ProteinTools.c()) {
            return true;
        } else if (s == ProteinTools.g()) {
            return true;
        } else if (s == ProteinTools.t()) {
            return true;
        }
        return false;
    }
    
    private void dumpAlign(SymbolList align, String name, PrintStream out0, PrintStream out1)
    		throws Exception
    	{
        int realBases = 0;
        List l0 = new ArrayList();
        List l1 = new ArrayList();
        boolean needGapRun = false;
        
        for (Iterator i = align.iterator(); i.hasNext(); ) {
            BasisSymbol sym = (BasisSymbol) i.next();
            if (isMatch(sym)) {
                if (needGapRun) {
                    List nList = Collections.nCopies(3, ProteinTools.n());
                    l0.addAll(nList);
                    l1.addAll(nList);
                    needGapRun = false;
                }
                List symList = sym.getSymbols();
                l0.add(symList.get(0));
                l1.add(symList.get(1));
                ++realBases;
            } else {
                needGapRun = true;
            }
        }
        
        if (realBases > 50) {
            new FastaFormat().writeSequence(
                    new SimpleSequence(
                            new SimpleSymbolList(ProteinTools.getAlphabet(), l0),
                            null,
                            name,
                            Annotation.EMPTY_ANNOTATION
                    ),
                    out0
           );
           new FastaFormat().writeSequence(
                    new SimpleSequence(
                            new SimpleSymbolList(ProteinTools.getAlphabet(), l1),
                            null,
                            name,
                            Annotation.EMPTY_ANNOTATION
                    ),
                    out1
           );
        }
    	}
    
    private SymbolList processAlignment(SymbolList rawAlign)
		throws Exception
	{
		Alphabet alignAlpha = rawAlign.getAlphabet();
		Symbol alignGap = alignAlpha.getGapSymbol();
		List<Symbol> sl = new ArrayList<Symbol>();
		for (Iterator<Symbol> i = rawAlign.iterator(); i.hasNext(); ) {
		    Symbol s = i.next();
		    if (! (s instanceof BasisSymbol)) {
		        continue;
		    }
		    if (s == alignGap) {
		        continue;
		    }
		    sl.add(s);
		}
		SymbolList niceAlign = new SimpleSymbolList(alignAlpha, sl);
		return niceAlign;
	}
    
    private static SequenceDB loadDB(File f)
	    throws Exception
	{
	    SequenceIterator si = SeqIOTools.readFastaProtein(new BufferedReader(new FileReader(f)));
	    SequenceDB seqDB = new HashSequenceDB();
	    while (si.hasNext()) {
	        Sequence seq = si.nextSequence();
	        seqDB.addSequence(seq);
	    }
	    return seqDB;
	}
}
