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

import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.utils.CliTools;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author td2
 */
public class MotifSetComparisonLGL {
    private boolean revComp = true;
    private double baseWeight = 100;

    public void setBaseWeight(double d) {
        this.baseWeight = d;
    }
    
    public void setRevComp(boolean b) {
        this.revComp = b;
    }
    
	public static void main(String[] args) 
		throws Exception
	{
		MotifSetComparisonLGL app = new MotifSetComparisonLGL();
		args = CliTools.configureBean(app, args);
		app.run(args);
	}
	
	public void run(String[] args)
		throws Exception
	{
		WeightMatrix[] set0 = loadWms(new File(args[0]));
		// WeightMatrix[] set1 = loadWms(new File(args[1]));
		String set0prefix = args[1];
		// String set1prefix = args[3];
		
		WeightMatrix[] set;
		String[] names;
		{
		    List<WeightMatrix> setL = new ArrayList<WeightMatrix>();
		    List<String> namesL = new ArrayList<String>();
		    for (int w = 0; w < set0.length; ++w) {
		        setL.add(set0[w]);
		        namesL.add(set0prefix + w);
		    }
		    /*
		    
		    for (int w = 0; w < set1.length; ++w) {
		        setL.add(set1[w]);
		        namesL.add(set1prefix + w);
		    }
		    */
		    
		    set = setL.toArray(new WeightMatrix[0]);
		    names = namesL.toArray(new String[0]);
		}
		
		Distribution elsewhere = new UniformDistribution((FiniteAlphabet) set0[0].getAlphabet());
		
		for (int i = 1; i < set.length; ++i) {
		    System.out.println("# " + names[i]);
		    for (int j = 0; j < i; ++j) {
				double score = cfWms(set[i], elsewhere, set[j], elsewhere);
				double hScore = 10000.0 / Math.pow(score, 4);
				if (hScore < 10.0) {
					System.out.println(
						names[j] + "\t" +
						hScore
					);
				}
			}
		}
	}
	
	public double cfWms(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
		throws Exception
	{
	    double score = doCfWms(wm0, pad0, wm1, pad1);
	    if (revComp) {
	        score = Math.min(score, doCfWms(wm0, pad0, WmTools.reverseComplement(wm1), pad1));
	    }
	    return score;
	}
	
	public double doCfWms(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
		throws Exception
	{
		double bestScore = Double.POSITIVE_INFINITY;
		int minPos = -wm1.columns();
		int maxPos = wm0.columns() + wm1.columns();
		for (int offset = -wm1.columns(); offset <= wm0.columns(); ++offset) {
			double score = 0.0;
			for (int pos = minPos; pos <= maxPos; ++pos) {
				Distribution col0 = pad0, col1 = pad1;
				if (pos >= 0 && pos < wm0.columns()) {
					col0 = wm0.getColumn(pos);
				}
				int opos = pos - offset;
				if (opos >= 0 && opos < wm1.columns()) {
					col1 = wm1.getColumn(opos);
				}
				double cScore = 0;
				for (Iterator i = ((FiniteAlphabet) col0.getAlphabet()).iterator(); i.hasNext(); ) {
					Symbol s = (Symbol) i.next();
					cScore += Math.pow(col0.getWeight(s) - col1.getWeight(s), 2);
				}
				score += Math.sqrt(cScore);
			}
			bestScore = Math.min(score, bestScore);
		}
		return bestScore;
	}
	
	public static WeightMatrix[] loadWms(File f)
		throws Exception
	{
		if (f.getName().endsWith(".jos")) {
		    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(f));
		    List l = new ArrayList();
		    try {
		        while (true) {
		            l.add(ois.readObject());
		        }
		    } catch (Exception ex) {
		    }
		    return (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
		} else {
		    Motif[] d = MotifIOTools.loadMotifSetXML(new FileInputStream(f));
		    WeightMatrix[] motifs = new WeightMatrix[d.length];
		    for (int m = 0; m < d.length; ++m) {
		        motifs[m] = d[m].getWeightMatrix();
		    }
		    return motifs;
		}
}
}
