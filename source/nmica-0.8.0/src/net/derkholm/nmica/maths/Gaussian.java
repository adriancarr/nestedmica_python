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

package net.derkholm.nmica.maths;

/**
 * Gaussian (normal) distribution.
 *
 * @author Thomas Down
 */

public class Gaussian {
    private Gaussian() {
    }
    
    public static double standardVariate() {
        return Math.cos(Math.random() * 6.28) * 
               Math.sqrt(2 * Math.log(1 / Math.random()));
    }
    
    public static double probability(double x, double mean, double precision) {
        return NativeMath.log2(Math.sqrt(precision / 2 / Math.PI) * Math.exp(-precision / 2 * Math.pow(x - mean, 2)));
    }
}
