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

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.symbol.FiniteAlphabet;

import java.util.HashMap;
import java.util.Map;

/**
 * Utility methods for working with weight matrices.
 *
 * @author Thomas Down
 */

public class DistTools {
    private final static Map<FiniteAlphabet, Distribution> uniDists = new HashMap<FiniteAlphabet, Distribution>();
    
    private DistTools() {
    }
    
    public static Distribution getUniformDistribution(FiniteAlphabet fa) {
        Distribution uni = uniDists.get(fa);
        if (uni == null) {
            uni = new UniformDistribution(fa);
            uniDists.put(fa, uni);
        }
        return uni;
    }
}
