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

import java.util.*;
import net.derkholm.nmica.utils.*;

/**
 * Default implementation of FacetteMap.
 *
 * @author Thomas Down
 */

public class SimpleFacetteMap implements FacetteMap, java.io.Serializable {
	private static final long serialVersionUID = 9080147573813775721L;
	
	private final ContributionGroup[] contributions;
    private final Facette[] facettes;
    private final boolean[] matrix;
    
    /**
     * Construct a new FacetteMap as a copy of the specificed map.
     */
    
    public SimpleFacetteMap(FacetteMap fm) {
        this.contributions = (ContributionGroup[]) ArrayTools.copy(fm.getContributionGroups());
        this.facettes = (Facette[]) ArrayTools.copy(fm.getFacettes());
        matrix = new boolean[contributions.length * facettes.length];
        for (int c = 0; c < contributions.length; ++c) {
            for (int f = 0; f < facettes.length; ++f) {
                matrix[c * contributions.length + f] = fm.contributesToFacette(contributions[c], facettes[f]);
            }
        }
    }
    
    public SimpleFacetteMap(ContributionGroup[] contributions, Facette[] facettes) {
        this.contributions = contributions;
        this.facettes = facettes;
        matrix = new boolean[contributions.length * facettes.length];
    }
    
    public ContributionGroup[] getContributionGroups() {
        return contributions;
    }
    
    public Facette[] getFacettes() {
        return facettes;
    }
    
    public ContributionGroup getContributionForFacette(Facette f) 
        throws IllegalArgumentException
    {
        for (int c = 0; c < contributions.length; ++c) {
            if (contributesToFacette(contributions[c], f)) {
                return contributions[c];
            }
        }
        throw new NoSuchElementException("No group contributing to " + f.toString());
    }
    
    public Facette[] getFacettesForContribution(ContributionGroup cg)
        throws IllegalArgumentException
    {
        List<Facette> l = new ArrayList<Facette>();
        for (int f = 0; f < facettes.length; ++f) {
            if (contributesToFacette(cg, facettes[f])) {
                l.add(facettes[f]);
            }
        }
        return l.toArray(new Facette[l.size()]);
    }
    
    public boolean contributesToFacette(ContributionGroup cg, Facette f)
        throws IllegalArgumentException
    {
        int cgi = indexOf(contributions, cg);
        int fi = indexOf(facettes, f);
        return matrix[cgi * contributions.length + fi];
    }
    
    public void setContributesToFacette(ContributionGroup cg, Facette f, boolean b) {
        int cgi = indexOf(contributions, cg);
        int fi = indexOf(facettes, f);
        
        matrix[cgi * contributions.length + fi] = b;
    }
    
    private int indexOf(Object[] array, Object o) 
        throws IllegalArgumentException
    {
        for (int i = 0; i < array.length; ++i) {
            if (array[i] == o) {
                return i;
            }
        }
        throw new IllegalArgumentException();
    }
}
