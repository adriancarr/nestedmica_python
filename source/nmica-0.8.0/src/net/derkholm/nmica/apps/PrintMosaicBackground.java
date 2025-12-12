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
import java.util.Iterator;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.bjv2.util.cli.App;

import net.derkholm.nmica.model.motif.MosaicSequenceBackground;

/**
 * @author thomas
 */

@App(overview="Print a human-readable version of a NestedMICA background model", generateStub=true)
public class PrintMosaicBackground {
    public static void main(String[] args)
    		throws Exception
    {
        File bgFile = new File(args[0]);
        MosaicSequenceBackground bg = (MosaicSequenceBackground) new ObjectInputStream(new FileInputStream(bgFile)).readObject();
        Distribution[] dists = bg.getBackgroundDistributions();
        FiniteAlphabet alpha = (FiniteAlphabet) dists[0].getAlphabet();
        
        int maxSymName = 0;
        for (Iterator<?> i = alpha.iterator(); i.hasNext(); ) {
            maxSymName = Math.max(maxSymName, ((Symbol) i.next()).getName().length());
        }
        
        for (int d = 0; d < dists.length; ++d) {
            System.out.println("Mosaic class " + d);
            for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
                Symbol s = (Symbol) si.next();
                System.out.println("    " + pad(s.getName(), maxSymName + 4) + dists[d].getWeight(s));
            }
            System.out.println();
        }
    }
    
    private static String pad(String s, int len)
    		throws Exception
    	{
        StringBuffer sb = new StringBuffer();
        sb.append(s);
        while (sb.length() < len) {
            sb.append(' ');
        }
        return sb.toString();
    	}
}
