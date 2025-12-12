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

import java.io.Serializable;

/**
 * @author thomas
 */
public class LogisticFunction implements DoubleFunction, Serializable {
    public static final DoubleFunction INSTANCE;
    
    static {
        INSTANCE = new LogisticFunction();
    }
    
    private LogisticFunction() {
    }
    
    private Object readResolve() {
        return INSTANCE;
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.maths.DoubleFunction#eval(double)
     */
    public double eval(double d) {
        return 1.0 / (1 + Math.exp(-d));
    }
}
