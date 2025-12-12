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

import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.utils.CliTools;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author thomas
 */
public class ConsensusToWM {
    private File out;
    
    public void setOut(File f) {
        this.out = f;
    }
    
    public static void main(String[] args)
	    throws Exception
	{
        ConsensusToWM app = new ConsensusToWM();
        args = CliTools.configureBean(app, args);
        app.run(args);
	}
    
    public void run(String[] motifs)
    		throws Exception
    {
        ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(out));
        
        for (int i = 0; i < motifs.length; ++i) {
            SymbolList sm = DNATools.createDNA(motifs[i]);
            WeightMatrix wm = new SimpleWeightMatrix(sm.getAlphabet(), sm.length(), NMSimpleDistribution.FACTORY);
            DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
            for (int c = 0; c < wm.columns(); ++c) {
                dtc.registerDistribution(wm.getColumn(c));
            }
            for (int c = 0; c < wm.columns(); ++c) {
                Symbol s =  sm.symbolAt(c + 1);
                Distribution dist = wm.getColumn(c);
                for (Iterator mi = ((FiniteAlphabet) s.getMatches()).iterator(); mi.hasNext(); ) {
                    dtc.addCount(dist, (Symbol) mi.next(), 1.0);
                }
            }
            dtc.train();
            
            oos.writeObject(wm);
        }
        
        oos.close();
	}

	private static SymbolList atomize(SymbolList sl)
	    throws Exception
	{
	    List<Symbol> atoms = new ArrayList<Symbol>();
	    for (Iterator i = sl.iterator(); i.hasNext(); ) {
	        Symbol s = (Symbol) i.next();
	        Symbol[] options = alphaToArray((FiniteAlphabet) s.getMatches());
	        atoms.add(options[(int) Math.floor(Math.random() * options.length)]);
	    }
	    return new SimpleSymbolList(sl.getAlphabet(), atoms);
	}
	
	private static Symbol[] alphaToArray(FiniteAlphabet a)
	    throws Exception
	{
	    List<Symbol> l = new ArrayList<Symbol>();
	    for (Iterator<Symbol> ai = a.iterator(); ai.hasNext(); ) {
	        l.add(ai.next());
	    }
	    return l.toArray(new Symbol[0]);
	}
}
