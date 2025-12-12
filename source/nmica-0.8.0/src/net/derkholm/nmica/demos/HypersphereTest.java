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
package net.derkholm.nmica.demos;

public class HypersphereTest {
    public static void main(String[] args)
        throws Exception
    {
        HypersphereTest app = new HypersphereTest();
        app.run();
    }
    
    public void run() {
        double a1 = randomAngle(), a2 = randomAngle() + Math.PI, a3 = randomAngle() + Math.PI;
        System.out.println("a1 = " + a1);
        System.out.println("a2 = " + a2);
        System.out.println("a3 = " + a3);
        
        double f1 = Math.sin(a3) * Math.sin(a2) * Math.cos(a1);
        double f2 = Math.sin(a3) * Math.sin(a2) * Math.sin(a1);
        double f3 = Math.sin(a3) * Math.cos(a2);
        double f4 = Math.cos(a3);
        
        System.out.println();
        
        System.out.println("f1 = " + f1);
        System.out.println("f2 = " + f2);
        System.out.println("f3 = " + f3);
        System.out.println("f4 = " + f4);
        
        System.out.println();
        System.out.println("Sanity check: " + Math.sqrt(f1*f1 + f2*f2 + f3*f3 + f4*f4));
        
        double a1p = Math.atan2(f2, f1);
        double a3p = Math.acos(f4);
        double a2p = Math.acos(f3 / Math.sin(a3p));
        
        System.out.println();
        System.out.println("a1p = " + a1p);
        System.out.println("a2p = " + a2p);
        System.out.println("a3p = " + a3p);
        
        double f1p = Math.sin(a3p) * Math.sin(a2p) * Math.cos(a1p);
        double f2p = Math.sin(a3p) * Math.sin(a2p) * Math.sin(a1p);
        double f3p = Math.sin(a3p) * Math.cos(a2p);
        double f4p = Math.cos(a3p);
        
                System.out.println();
        
        System.out.println("f1p = " + f1p);
        System.out.println("f2p = " + f2p);
        System.out.println("f3p = " + f3p);
        System.out.println("f4p = " + f4p);
        
        System.err.println("Am I a deviant? " + Math.sqrt(Math.pow(f1p - f1, 2) + Math.pow(f2p - f2, 2) + Math.pow(f3p - f3, 2) + Math.pow(f4p - f4, 2)));
    }
    
    private double randomAngle() {
        return Math.PI * (1.0 - 2 * Math.random());
    }
}
