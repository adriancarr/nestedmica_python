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

package net.derkholm.nmica.seq.consensus;

import java.io.Serializable;
import java.util.*;

import org.biojava.utils.*;
import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.dist.*;
import net.derkholm.nmica.seq.*;

/**
 * Distribution object which separates the base preference from the information
 * content.
 *
 * @author Thomas Down
 */

public class ConsensusDistribution extends AbstractDistribution implements Serializable {
    static final long serialVersionUID = 7736616807552328793L;
    
    private transient Distribution scores;
    private transient Distribution nullModel;
    private final Count direction;
    private final double hardness;
    
    public ConsensusDistribution(Count direction, double hardness)
    {
        this.direction = direction;
        this.hardness = hardness;
    }
    
    /**
     * Return the direction-vector for this distribution.
     */
    
    public Count getDirection() {
        return direction;
    }
    
    /**
     * Return the hardness parameter for this distribution.
     */
     
    public double getHardness() {
        return hardness;
    }
    
    public Alphabet getAlphabet() {
        return direction.getAlphabet();
    }
    
    protected void setWeightImpl(AtomicSymbol s, double d) 
        throws IllegalSymbolException, ChangeVetoException
    {
        throw new ChangeVetoException("Can't modify consensus distributions directly");
    }
    
    protected double getWeightImpl(AtomicSymbol s)
        throws IllegalSymbolException
    {
        return getScores().getWeight(s);
    }
    
    public Distribution getNullModel() {
        if (nullModel == null) {
            nullModel = new UniformDistribution((FiniteAlphabet) getAlphabet());
        }
        return nullModel;
    }
    
    public void setNullModelImpl(Distribution dist)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("That's soooo last year");
    }        
    
    protected Distribution getScores() {
        if (scores == null) {
            try {
                FiniteAlphabet fa = (FiniteAlphabet) direction.getAlphabet();
                Distribution _scores = new NMSimpleDistribution(fa);
                Distribution _nullo = getNullModel();
                double norm = 0;
                for (Iterator si = fa.iterator(); si.hasNext(); ) {
                    AtomicSymbol s = (AtomicSymbol) si.next();
                    norm += Math.exp(
                        Math.log(_nullo.getWeight(s)) + 
                        hardness * direction.getCount(s)
                    );
                }
                for (Iterator si = fa.iterator(); si.hasNext(); ) {
                    AtomicSymbol s = (AtomicSymbol) si.next();
                    _scores.setWeight(
                        s,
                        Math.exp(
                            Math.log(_nullo.getWeight(s)) +
                            hardness * direction.getCount(s)
                        ) / norm
                    );
                }
                scores = _scores;
            } catch (Exception ex) {
                throw new BioError("Assertion failed: madness in score calculation", ex);
            }
        }
        return scores;
    }
}
