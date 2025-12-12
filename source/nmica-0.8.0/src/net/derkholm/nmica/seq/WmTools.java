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

package net.derkholm.nmica.seq;

import java.util.*;

import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.*;
import org.biojava.bio.dist.*;
import org.biojava.bio.dp.*;
import org.biojava.utils.*;

/**
 * Utility methods for working with weight matrices.
 *
 * @author Thomas Down
 */

public class WmTools {
    private WmTools() {
    }
    
    public static WeightMatrix reverseComplement(WeightMatrix wm)
        throws IllegalAlphabetException
    {
        try {
            int cols = wm.columns();
            Distribution[] newDists = new Distribution[cols];
            for (int c = 0; c < cols; ++c) {
                Distribution oldDist = wm.getColumn(cols - c - 1);
                FiniteAlphabet alpha = (FiniteAlphabet) oldDist.getAlphabet();
                Distribution newDist = newDists[c] = new NMSimpleDistribution(alpha);
                for (Iterator si = alpha.iterator(); si.hasNext(); ) {
                    Symbol s = (Symbol) si.next();
                    newDist.setWeight(s, oldDist.getWeight(DNATools.complement(s)));
                }
            }
            return new SimpleWeightMatrix(newDists);
        } catch (IllegalSymbolException ex) {
            throw new IllegalAlphabetException(ex, "Couldn't reverse complement");
        } catch (ChangeVetoException ex) {
            throw new BioError("Assertion failed: couldn't modify newly created distribution", ex);
        }
    }

    /**
     * @param matrix
     * @return
     */
    public static String consensus(WeightMatrix wm) {
        StringBuffer sb = new StringBuffer();
        
        try {
            SymbolTokenization toke = null;
            try {
                toke = wm.getAlphabet().getTokenization("token");
            } catch (Exception ex) {
            }
            for (int pos = 0; pos < wm.columns(); ++pos) {
                double max = 0.0;
                Symbol maxSym = null;
                Distribution col = wm.getColumn(pos);
                for (Iterator si = ((FiniteAlphabet) wm.getAlphabet()).iterator(); si.hasNext(); ) {
                    Symbol sym = (Symbol) si.next();
                    double w = col.getWeight(sym);
                    if (w > max) {
                        max = w;
                        maxSym = sym;
                    }
                }
                String token = "x";
                if (toke != null) {
                    token = toke.tokenizeSymbol(maxSym);
                }
                sb.append(token);
            }
        } catch (Exception ex) {
            throw new BioError(ex);
        }
        
        return sb.toString();
    }
    
    
    /**
     * @param matrix
     * @return
     */
    public static String fluffyConsensus(WeightMatrix wm) {
        StringBuffer sb = new StringBuffer();
        
        try {
            SymbolTokenization toke = null;
            try {
                toke = wm.getAlphabet().getTokenization("token");
            } catch (Exception ex) {
            }
            for (int pos = 0; pos < wm.columns(); ++pos) {
                double max = 0.0;
                Symbol maxSym = null;
                Distribution col = wm.getColumn(pos);
                for (Iterator si = ((FiniteAlphabet) wm.getAlphabet()).iterator(); si.hasNext(); ) {
                    Symbol sym = (Symbol) si.next();
                    double w = col.getWeight(sym);
                    if (w > max) {
                        max = w;
                        maxSym = sym;
                    }
                }
                String token = "x";
                if (toke != null) {
                    token = toke.tokenizeSymbol(maxSym);
                }
                if (max > 0.8) {
                    token = token.toUpperCase();
                } else if (max > 0.3) {
                    token = token.toLowerCase();
                } else {
                    token = "n";
                }
                sb.append(token);
            }
        } catch (Exception ex) {
            throw new BioError(ex);
        }
        
        return sb.toString();
    }
}
