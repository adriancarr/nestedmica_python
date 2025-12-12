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
/*
 * Created on Sep 13, 2004
 */
package net.derkholm.nmica.demos;

import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix1D;
import net.derkholm.nmica.model.BinaryMixPolicy;
import net.derkholm.nmica.model.MixPolicy;
import net.derkholm.nmica.utils.CliTools;

/**
 * @author thomas
 */
public class MixPolicyTest {
    private int samples = 5000;
    private int components = 10;
    private double occupancy = 0.5;
    
    public void setSamples(int d) {
        this.samples = d;
    }
    public void setComponents(int d) {
        this.components = d;
    }
    public void setOccupancy(double d) {
        this.occupancy = d;
    }
     
    public static void main(String[] args) 
    		throws Exception
    {
        MixPolicyTest app = new MixPolicyTest();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        	MixPolicy policy = new BinaryMixPolicy(occupancy);
        	Matrix1D trix = new SimpleMatrix1D(components);
        	int[] counts = new int[components + 1];
        	
        	for (int c = 0; c < samples; ++c) {
        	    policy.variate(trix);
        	    int cnt = 0;
        	    for (int i = 0; i < components; ++i) {
        	        if (trix.get(i) != 0.0) {
        	            ++cnt;
        	        }
        	    }
        	    ++counts[cnt];
        	}
        	
        	for (int i = 0; i < components; ++i) {
        	    System.out.println("" + i + '\t' + ((1.0 * counts[i]) / samples));
        	}
    }
}
