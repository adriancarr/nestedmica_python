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

import net.derkholm.nmica.matrix.*;

/**
 * General purpose mixture model for some dataset.
 * The data involved can have multiple "facettes", each of which is
 * modelled separately, based on the same mixing process.  The precise
 * details of how multiple contributions can be mixed together to model
 * a particular facette of the data can be defined in implementations of
 * the <code>Facette</code> interface.
 *
 * @author Thomas Down
 */

public interface MultiICAModel {
    /**
     * Return an object which specifies the layout of this model.
     */
    
    public FacetteMap getFacetteMap();
    
    /**
     * Return the number of components in the model
     */
     
    public int getComponents();
    
    /**
     * Return the contributions for the requested facette
     */
     
    public ObjectMatrix1D getContributions(ContributionGroup f);
    
    /**
     * Return the contribution of a given component to the requested
     * model facette.
     */
    
    public ContributionItem getContribution(ContributionGroup cg, int component);
    
    /**
     * Return this model's dataset.
     */
    
    public Datum[] getDataSet();
    
    /**
     * Return the (log-2) likelihood of the dataset given the model.
     */
     
    public double likelihood();
    
    /**
     * Return the mixing coefficients for this model.
     */
     
    public Matrix1D getMixture(int datumIndex);
}
