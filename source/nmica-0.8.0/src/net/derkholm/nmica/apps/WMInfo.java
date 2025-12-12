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

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

/**
 * @author td2
 *
 */
public class WMInfo {
	private static final double LOG_2 = Math.log(2.0);
	
   public static void main(String[] args)
       throws Exception
   {
       for (int f = 0; f < args.length; ++f) {
	       File model = new File(args[f]);
	       
	       WeightMatrix[] wms;
	       {
	           ObjectInputStream ois = new ObjectInputStream(new FileInputStream(model));
	           List l = new ArrayList();
	           try {
	               while (true) {
	                   l.add(ois.readObject());
	               }
	           } catch (Exception ex) {
	           }
	           wms = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
	       }
	       
	       double i = 0;
	       for (int w = 0; w < wms.length; ++w) {
	       		double wi = wmInfo(wms[w]);
	       		if (args.length == 1) {
	       		    System.out.println("" + w + '\t' + wi);
	       		}
	       		i += wi;
	       }
	       if (args.length == 1) {
	           System.out.println("Total: " + i);
	       } else {
	           System.out.println(model.getName() + '\t' + i);
	       }
       }
   }
   
    private static double wmInfo(WeightMatrix wm)
   	    throws Exception
    {
    	double i = 0;
    	for (int c = 0; c < wm.columns(); ++c) {
    		i += info(wm.getColumn(c));
    	}
    	return i;
	}

	/**
	 * @param column
	 * @return
	 */
	private static double info(Distribution column) 
		throws Exception
	{
		double i = 0;
		FiniteAlphabet alpha = (FiniteAlphabet) column.getAlphabet();
		for (Iterator si = alpha.iterator(); si.hasNext(); ) {
			double w = column.getWeight((Symbol) si.next());
			if (w > 0) {
				i += w * Math.log(1/w) / LOG_2;
			}
		}
		return 2.0 - i;
	}
}
