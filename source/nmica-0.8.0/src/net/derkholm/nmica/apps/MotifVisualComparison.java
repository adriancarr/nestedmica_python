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

import net.derkholm.nmica.gui.WMPanel;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.utils.CliTools;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

import javax.swing.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author td2
 */
public class MotifVisualComparison {
    private double threshold = Double.NEGATIVE_INFINITY;
    private boolean revComp = false;
    private boolean brh = false;
    
    public void setBrh(boolean b) {
        this.brh = b;
    }
    
    public void setThreshold(double d) {
        this.threshold = d;
    }
    
    public void setRevComp(boolean b) {
        this.revComp = b;
    }
    
	public static void main(String[] args) 
		throws Exception
	{
		MotifVisualComparison app = new MotifVisualComparison();
		args = CliTools.configureBean(app, args);
		app.run(args);
	}
	
	private static class Pair implements Comparable<Pair> {
	    public final WeightMatrix wm0;
	    public final WeightMatrix wm1;
	    public final boolean flip;
	    public final double score;
	    
	    public Pair(WeightMatrix wm0, WeightMatrix wm1, double score, boolean flip) {
	        this.wm0 = wm0;
	        this.wm1 = wm1;
	        this.score = score;
	        this.flip = flip;
	    }


        public int compareTo(Pair po) {
            double dif = score - po.score;
            if (dif < 0) {
                return -1;
            } else if (dif > 0) {
                return 1;
            } else {
                return 0;
            }
        }
        
        public boolean equals(Object o) {
            Pair po = (Pair) o;
            return wm0 == po.wm0 && wm1 == po.wm1;
        }
        
        public int hashCode() {
            return wm0.hashCode() * 37 + wm1.hashCode();
        }
	}
	
	public List<Pair> bestHits(WeightMatrix[] set0, WeightMatrix[] set1, boolean flip)
		throws Exception
	{
	    List<Pair> pairList = new ArrayList<Pair>();
	    Distribution elsewhere = new UniformDistribution((FiniteAlphabet) set0[0].getAlphabet());
	    for (int i = 0; i < set0.length; ++i) {
	        int bestJ = -1;
	        double bestScore = Double.NEGATIVE_INFINITY;
	        boolean bestIsFlipped = false;
	        for (int j = 0; j < set1.length; ++j) {
	            {
	                double score = doCfWms(set0[i], elsewhere, set1[j], elsewhere) /* - nullCfWms(set0[i], elsewhere, set1[j], elsewhere) */;
	                if (score > bestScore) {
	                    bestScore = score;
	                    bestJ = j;
	                    bestIsFlipped = false;
	                }
	            }
	            {
	                double score = doCfWms(set0[i], elsewhere, WmTools.reverseComplement(set1[j]), elsewhere) /* - nullCfWms(set0[i], elsewhere, WmTools.reverseComplement(set1[j]), elsewhere) */;
	                if (score > bestScore) {
	                    bestScore = score;
	                    bestJ = j;
	                    bestIsFlipped = true;
	                }
	            }
	            
	        }
	        if (bestScore > threshold) {
	            if (!flip) {
	                pairList.add(new Pair(set0[i], set1[bestJ], bestScore, bestIsFlipped));
	            } else {
	                pairList.add(new Pair(set1[bestJ], set0[i], bestScore, bestIsFlipped));
	            }
	        }
	    }
	    return pairList;
	}
	
	public void run(String[] args)
		throws Exception
	{
	    File f0 = new File(args[0]);
	    File f1 = new File(args[1]);
		WeightMatrix[] set0 = loadWms(f0);
		WeightMatrix[] set1 = loadWms(f1);
		
		DecimalFormat fmt = new DecimalFormat("###,###.00;-###,###.00");
		
		JFrame viewers = new JFrame("MotifVisualComparison: " + f0.getName() + " vs " + f1.getName());
		Box b = Box.createVerticalBox();
		List<Pair> pairList;
		
		if (brh) {
		    Set<Pair> h0 = new HashSet<Pair>(bestHits(set0, set1, false));
		    Set<Pair> h1 = new HashSet<Pair>(bestHits(set1, set0, true));
		    h0.retainAll(h1);
		    pairList = new ArrayList<Pair>(h0);
		} else {
		    pairList = bestHits(set0, set1, false);
		}
		
	    Collections.sort(pairList);
		Collections.reverse(pairList);
		for (Pair p: pairList) {
		    Box h = Box.createHorizontalBox();
		    h.add(new WMPanel(p.wm0));
		    h.add(Box.createHorizontalStrut(20));
	        h.add(new JLabel(fmt.format(p.score)));
	        h.add(Box.createHorizontalStrut(20));
	        h.add(new WMPanel(p.flip ? WmTools.reverseComplement(p.wm1) : p.wm1));
	        b.add(h);
		}
		
		viewers.getContentPane().add(b);
        viewers.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent wev) {
                System.exit(0);
            }
        } );
        viewers.pack();
		viewers.setVisible(true);
	}
	
	public double nullCfWms(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
		throws Exception
	{
	    double cPadPad = 0;
	    for (Iterator i = ((FiniteAlphabet) pad0.getAlphabet()).iterator(); i.hasNext(); ) {
	        Symbol s = (Symbol) i.next();
	        cPadPad += pad0.getWeight(s) * pad1.getWeight(s);
	    }
	    double padPad = Math.log(cPadPad) / Math.log(2.0);
	    double tot = 0;
	    int minPos = -wm1.columns();
		int maxPos = wm0.columns() + wm1.columns();
		for (int offset = -wm1.columns(); offset <= wm0.columns(); ++offset) {
		    tot += padPad;
		}
		return tot;
	}
	
	public double doCfWms(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
		throws Exception
	{
		double bestScore = Double.NEGATIVE_INFINITY;
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
					cScore += col0.getWeight(s) * col1.getWeight(s);
				}
				score += Math.log(cScore) / Math.log(2.0);
			}
			bestScore = Math.max(score, bestScore);
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
