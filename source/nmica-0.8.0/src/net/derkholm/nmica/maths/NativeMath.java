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
 * General-purpose mathematical utility routines, including base-2 logarithms.
 * This class includes some fast approximation routines.
 *
 * @author Thomas Down
 */

public class NativeMath {
    static {
        System.loadLibrary("nmica");
    }

    
    private NativeMath() {
    }
    
    /**
     * Logarithm to base e using standard system libraries.
     * 
     * @param d
     * @return
     */
    
    public static native double log(double d);
    
    /**
     * e**d using standard system libraries.
     * @param d
     * @return
     */
    
    public static native double exp(double d);
    
    public static native double loopLog(double d, int iter);
    
    public static native double loopExp(double d, int iter);
    
    /**
     * Calculate <code>log(e**x + e**y)</code>.  Implemented in
     * native code so it only requires one JNI call.
     */
    
    public static native double addLog(double x, double y);
    
    /**
     * Logarithm to base 2.
     * 
     * @param d
     * @return
     */
    
    public static native double log2(double d);
    
    /**
     * Calculate <code>2**d</code>.
     * 
     * @param d
     * @return
     */
    
    public static native double exp2(double d);
    
    /**
     * Calculate <code>log_2(2**x + 2**y)</code>.
     * 
     * @param x
     * @param y
     * @return
     */
    
    public static native double addLog2(double x, double y);
    
    public static native double addLog2(double[] x);
    
    /**
     * Calculate a fast approximation of <code>log_2(d)</code>.
     * 
     * @param d
     * @return
     */
    
    public static native double fastlog2(double d);
}
