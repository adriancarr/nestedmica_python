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

import net.derkholm.nmica.maths.NativeMath;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.onehead.SingleDPMatrix;
import org.biojava.bio.symbol.*;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Interface for basic sequence models.
 *
 * @author Thomas Down
 */

public class MosaicSequenceBackground implements SequenceBackground, Serializable {
    static final long serialVersionUID = -742638802183258198L;
    
    private final static Pattern namePattern = Pattern.compile("patch([0-9]+)");
    private Distribution[] backgroundDistributions;
    private double backgroundTransition;
    private transient DP backgroundDP;
    
    public MosaicSequenceBackground(
        Distribution[] backgroundDistributions,
        double backgroundTransition
    ) {
        this.backgroundDistributions = backgroundDistributions;
        this.backgroundTransition = backgroundTransition;
    }
    
    public Distribution[] getBackgroundDistributions() {
        return backgroundDistributions;
    }
    
    public double getBackgroundTransition() {
        return backgroundTransition;
    }
    
    public DP getBackgroundDP() 
        throws Exception
    {
        if (backgroundDP == null) {
            MarkovModel mm = MosaicTools.makePatchModel(backgroundDistributions, backgroundTransition);
            backgroundDP = new org.biojava.bio.dp.onehead.SingleDP(mm);
        }
        return backgroundDP;
    }
    
    public int getMosaicClasses() {
        return backgroundDistributions.length;
    }
    
    public int getMosaicOrder() {
		int order = 1;
		{
			List al = backgroundDistributions[0].getAlphabet().getAlphabets();
			if (al.size() > 1) {
				order = al.size() + ((Alphabet) al.get(0)).getAlphabets().size() - 1;
			}
		}
		return order;
    }
    
    public Location backgroundSymbolLikelihood(SymbolList sl, double[] likelihoods)
        throws IllegalAlphabetException
    {
        try {
            DP dp = getBackgroundDP();
            Distribution[] background = getBackgroundDistributions();
            int backgroundOrder = getMosaicOrder();
           
            SymbolList seq = sl;
            if (backgroundOrder > 1) {
                seq = SymbolListViews.orderNSymbolList(seq, backgroundOrder);
            } 
            
            SingleDPMatrix forwardMatrix = (SingleDPMatrix) dp.forwardMatrix(
                new SymbolList[] {seq},
                ScoreType.PROBABILITY
            );
            double score = forwardMatrix.getScore();
            SingleDPMatrix backwardMatrix = (SingleDPMatrix) dp.backwardMatrix(
                new SymbolList[] {seq},
                ScoreType.PROBABILITY
            );
            
            List<Location> locList = new ArrayList<Location>();
            for (int pos = 1; pos <= seq.length(); ++pos) {
                Symbol sym = seq.symbolAt(pos);
                if (sym instanceof AtomicSymbol) {
                    State[] states = forwardMatrix.states();
                    double[] fcol = forwardMatrix.scores[pos];
                    double[] bcol = backwardMatrix.scores[pos];
                    double sh = 0;
                    for (int s = 0; s < states.length; ++s) {
                        Matcher nameMatcher = namePattern.matcher(states[s].getName());
                        if (nameMatcher.matches()) {
                            int p = Integer.parseInt(nameMatcher.group(1));
                            double weight = Math.exp(fcol[s] + bcol[s] - score);
                            sh += weight * background[p].getWeight(sym);
                        }
                    }
                    likelihoods[pos + backgroundOrder - 2] = NativeMath.log2(sh);
                    locList.add(new PointLocation(pos + backgroundOrder - 1));
                }
            }
            return LocationTools.union(locList);
        } catch (Exception ex) {
            throw new IllegalAlphabetException(ex, "Couldn't evaluated background");
        }
        
    }
}

