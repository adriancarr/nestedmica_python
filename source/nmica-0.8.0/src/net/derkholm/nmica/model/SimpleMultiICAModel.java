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

package net.derkholm.nmica.model;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.MatrixTools;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix1D;

/**
 * Stupid implementation of the MultiICAModel interface, which doesn't do
 * any caching.  This implementation is not <code>Commitable</code>, so it
 * probably can't be used for training.  The main use is for model persistance.
 *
 * @author Thomas Down
 */

public class SimpleMultiICAModel implements MultiICAModel, Serializable {
    static final long serialVersionUID = -7179481801889622800L;
    
    private final FacetteMap facettes;
    private final Datum[] data;
    private final int components;
    private Map<ContributionGroup,ObjectMatrix1D> contributions;
    private Matrix2D mixingMatrix;
    
    /**
     * Construct a new model with no contributions or mixing coefficients.
     *
     * <p>
     * (is this actually useful?)
     * </p>
     */
    
    public SimpleMultiICAModel(
        FacetteMap facettes,
        Datum[] data, 
        int components
    ) 
    {
        this.facettes = facettes;
        this.data = data;
        this.components = components;
        mixingMatrix = new SimpleMatrix2D(data.length, components);
        contributions = new HashMap<ContributionGroup, ObjectMatrix1D>();
        for (ContributionGroup cg : facettes.getContributionGroups()) {
            contributions.put(cg, new SimpleObjectMatrix1D(components));
        }
    }
    
    /**
     * Construct a new model object as a copy of the specified object
     */
    
    public SimpleMultiICAModel(MultiICAModel m) {
        this.facettes = m.getFacetteMap();
        this.data = m.getDataSet();
        this.components = m.getComponents();
        
        this.mixingMatrix = new SimpleMatrix2D(data.length, components);
        for (int d = 0; d < data.length; ++d) {
            MatrixTools.copy(MatrixTools.viewRow(mixingMatrix, d), m.getMixture(d));
        }
        this.contributions = new HashMap<ContributionGroup,ObjectMatrix1D>();
        for (ContributionGroup cg : facettes.getContributionGroups()) {
            contributions.put(cg, new SimpleObjectMatrix1D(m.getContributions(cg)));
        }
    }
     
    public FacetteMap getFacetteMap() {
        return facettes;
    }
    
    public int getComponents() {
        return components;
    }
    
    public Datum[] getDataSet() {
        return data;
    }
    
    public ContributionItem getContribution(ContributionGroup cg, int c) {
        return (ContributionItem) getContributions(cg).get(c);
    }
    
    public ObjectMatrix1D getContributions(ContributionGroup cg) {
        ObjectMatrix1D c = (ObjectMatrix1D) contributions.get(cg);
        if (c == null) {
            throw new IllegalArgumentException("This model doesn't include contribution group " + cg.getName());
        }
        return c;
    }
    
    public double likelihood() {
        Facette[] faces = facettes.getFacettes();
        
        double L = 0;
        for (int d = 0; d < data.length; ++d) {
            for (int f = 0; f < faces.length; ++f) {
                L += faces[f].getLikelihoodCalculator(data[d].getFacettedData()[f]).
                     likelihood(getContributions(facettes.getContributionForFacette(faces[f])), getMixture(d));
            }
        }
        return L;
    }
    
    public Matrix2D getMixingMatrix() {
        return mixingMatrix;
    }
    
    public Matrix1D getMixture(int datumIndex) {
        return MatrixTools.viewRow(mixingMatrix, datumIndex);
    }
}
