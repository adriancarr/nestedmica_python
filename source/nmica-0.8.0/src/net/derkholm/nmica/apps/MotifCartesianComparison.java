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

import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

/**
 * @author thomas
 */
public class MotifCartesianComparison {
    private double threshold = Double.POSITIVE_INFINITY;
    private boolean best = false;
    private boolean revComp = false;
    
    public void setThreshold(double d) {
        this.threshold = d;
    }
    
    public void setBest(boolean b) {
        this.best = b;
    }
    
    public void setRevComp(boolean b) {
        this.revComp = b;
    }
    
	public static void main(String[] args) 
		throws Exception
	{
		MotifCartesianComparison app = new MotifCartesianComparison();
		args = CliTools.configureBean(app, args);
		app.run(args);
	}
	
	public void run(String[] args)
		throws Exception
	{
		WeightMatrix[] set0 = loadWms(new File(args[0]));
		WeightMatrix[] set1 = loadWms(new File(args[1]));
		
		Distribution elsewhere = new UniformDistribution((FiniteAlphabet) set0[0].getAlphabet());
		
		for (int i = 0; i < set0.length; ++i) {
		    if (best) {
		        int bestJ = -1;
		        double bestScore = Double.NEGATIVE_INFINITY;
		        for (int j = 0; j < set1.length; ++j) {
				    double score = cfWms(set0[i], elsewhere, set1[j], elsewhere);
				    if (score > bestScore) {
				        bestScore = score;
				        bestJ = j;
				    }
		        }
		        if (bestScore < threshold) {
		            System.out.println("" + i + "->" + bestJ + ": " + bestScore);
		        }
		    } else {
				for (int j = 0; j < set1.length; ++j) {
				    double score = cfWms(set0[i], elsewhere, set1[j], elsewhere);
				    if (score < threshold) {
						System.out.println(
								"" + i + "->" + j + ": " +
								score
					    );
				    }
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
