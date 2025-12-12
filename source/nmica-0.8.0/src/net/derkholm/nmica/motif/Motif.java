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

package net.derkholm.nmica.motif;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.utils.AbstractChangeable;

/**
 * Sequence motif, maybe with some annotations attached.
 * 
 * @author thomas
 */

public class Motif extends AbstractChangeable implements Annotatable, Cloneable {
	public static final int NO_VERSION = -1;
	
    private Annotation annotation;
    private WeightMatrix weightMatrix;
    private String name;
    private double threshold;
    private int version = NO_VERSION;
    
    public Motif() {
        annotation = new SmallAnnotation();
    }
    
    public Motif(Motif m) {
    	this.annotation = new SmallAnnotation(m.getAnnotation());
    	this.name = m.getName();
    	this.threshold = m.getThreshold();
    	this.weightMatrix = m.getWeightMatrix();
    	this.version = m.getVersion();
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String name) {
        this.name = name;
    }
    public double getThreshold() {
        return threshold;
    }
    public void setThreshold(double threshold) {
        this.threshold = threshold;
    }
    public WeightMatrix getWeightMatrix() {
        return weightMatrix;
    }
    public void setWeightMatrix(WeightMatrix weightMatrix) {
        this.weightMatrix = weightMatrix;
    }
    public void setVersion(int i) {
    	this.version = i;
    }
    public int getVersion() {
    	return version;
    }

    /* (non-Javadoc)
     * @see org.biojava.bio.Annotatable#getAnnotation()
     */
    public Annotation getAnnotation() {
        return annotation;
    }
}
