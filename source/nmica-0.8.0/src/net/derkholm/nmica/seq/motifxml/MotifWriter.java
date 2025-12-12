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

package net.derkholm.nmica.seq.motifxml;

import java.util.Iterator;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.xml.XMLWriter;

/**
 * Write out a weightmatrix in XML format.
 * 
 * @author thomas
 */
public class MotifWriter {
    public void writeMatrix(WeightMatrix wm, XMLWriter xw)
	    throws Exception
	{
	    xw.openTag("weightmatrix");
	    xw.attribute("alphabet", wm.getAlphabet().getName());
	    xw.attribute("columns", "" + wm.columns());
	    for (int col = 0; col < wm.columns(); ++col) {
	        xw.openTag("column");
	        xw.attribute("pos", "" + col);
	        Distribution dist = wm.getColumn(col);
	        for (Iterator si = ((FiniteAlphabet) dist.getAlphabet()).iterator(); si.hasNext(); ) {
	            Symbol s = (Symbol) si.next();
	            xw.openTag("weight");
	            xw.attribute("symbol", s.getName());
	            xw.print("" + dist.getWeight(s));
	            xw.closeTag("weight");
	        }
	        xw.closeTag("column");
	    }
	    xw.closeTag("weightmatrix");
	}
}
