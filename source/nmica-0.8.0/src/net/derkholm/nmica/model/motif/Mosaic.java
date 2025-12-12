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

package net.derkholm.nmica.model.motif;

import java.io.Serializable;

import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.Alphabet;

public class Mosaic implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private final Distribution[] distributions;
    private final double transition;
    private final Annotation annotation;
    
    public Mosaic(Distribution[] dists, double trans) {
    	this.distributions = dists;
    	this.transition = trans;
    	this.annotation = new SmallAnnotation();
    }
    
    public Mosaic(Distribution[] dists, double trans, Annotation props) {
    	this.distributions = dists;
    	this.transition = trans;
    	this.annotation = new SmallAnnotation(props);
    }
    
    public Distribution[] getDistributions() {
    	return distributions;
    }
    
    public double getTransition() {
    	return transition;
    }
    
    public Annotation getAnnotation() {
    	return annotation;
    }
    
    public Alphabet getAlphabet() {
    	return distributions[0].getAlphabet();
    }
}
