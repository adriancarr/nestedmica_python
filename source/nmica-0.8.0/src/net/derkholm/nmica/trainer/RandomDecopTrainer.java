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

package net.derkholm.nmica.trainer;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.ContributionGroup;
import net.derkholm.nmica.model.ContributionItem;
import net.derkholm.nmica.model.ContributionPrior;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.Datum;
import net.derkholm.nmica.model.FacetteMap;
import net.derkholm.nmica.model.MixPolicy;
import net.derkholm.nmica.model.PenalizedVariate;
import net.derkholm.nmica.model.SimpleContributionItem;

/**
 * Trainer which decorrelates states by applying a random sequence
 * of Monte-carlo operations.
 * 
 * @author thomas
 */
public class RandomDecopTrainer extends Trainer {
    private static final long serialVersionUID = 8805382875360428627L;
    
    private int minMixtureMoves = 100;
    private int minContributionMoves = 100;
    private int minMixtureProposals = 100;
    private int minContributionProposals = 100;
    
    /**
     * @param facetteMap
     * @param data
     * @param components
     * @param priors
     * @param samplers
     * @param mixPolicy
     * @param ensembleSize
     */
    public RandomDecopTrainer(FacetteMap facetteMap, Datum[] data,
            int components, ContributionPrior[] priors,
            ContributionSampler[] samplers, MixPolicy mixPolicy,
            int ensembleSize) 
    {
        super(facetteMap, data, components, priors, samplers, mixPolicy,
                ensembleSize);
    }
    
    public void setMinMixtureMoves(int i) {
        minMixtureMoves = i;
    }
    
    public int getMinMixtureMoves() {
        return minMixtureMoves;
    }
    
    public void setMinContributionMoves(int i) {
        minContributionMoves = i;
    }
    
    public int getMinContributionMoves() {
        return minContributionMoves;
    }
    
    /**
     * @return Returns the minContributionProposals.
     */
    public int getMinContributionProposals() {
        return minContributionProposals;
    }
    /**
     * @param minContributionProposals The minContributionProposals to set.
     */
    public void setMinContributionProposals(int minContributionProposals) {
        this.minContributionProposals = minContributionProposals;
    }
    /**
     * @return Returns the minMixtureProposals.
     */
    public int getMinMixtureProposals() {
        return minMixtureProposals;
    }
    /**
     * @param minMixtureProposals The minMixtureProposals to set.
     */
    public void setMinMixtureProposals(int minMixtureProposals) {
        this.minMixtureProposals = minMixtureProposals;
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.trainer.Trainer#decorrelateState(net.derkholm.nmica.trainer.TrainableState, double)
     */
    protected double decorrelateState(TrainableState newModel, double minLikelihood, int seeded, int unseeded) 
    {
    	// System.err.printf("dcs, ml=%g l=%g%n", minLikelihood, newModel.likelihood());
    	
        Datum[] data = getDataSet();
        int components = getComponents();
        ContributionGroup[] cgs = getFacetteMap().getContributionGroups();
        ContributionSampler[] samplers = getSamplers();
        MixPolicy mixPolicy = getMixPolicy();
        
        double oldPrior = prior(newModel);
        double newHood = Double.NEGATIVE_INFINITY;
        int mixMoves = 0;
        int conMoves = 0;
        int mixProposals = 0;
        int conProposals = 0;
        
        while (true) {
        	int mixNeeded = Math.max(0, minMixtureMoves - mixMoves) + Math.max(0, minMixtureProposals - mixProposals);
            int conNeeded = Math.max(0, minContributionMoves - conMoves) + Math.max(0, minContributionProposals - conProposals);
            if (mixNeeded <= 0 || conNeeded <= 0) {
            	break;
            }
        	// System.err.printf("mixProposals=%d/%d mixMoves=%d/%d conProposals=%d/%d conMoves=%d/%d%n", mixProposals, minMixtureProposals, mixMoves, minMixtureMoves, conProposals, minContributionProposals, conMoves, minContributionProposals);
        	
            boolean mixSampled = false;
            boolean wasSilent = false;
            double balancePenalty = 0;
            if (Math.random() < ((1.0 * mixNeeded) / (mixNeeded + conNeeded))) {
            	// System.err.println("MixtureSample");
                do {
                    int d = (int) Math.floor(Math.random() * data.length);
                    mixPolicy.sample(newModel.getMixture(d));
                } while (Math.random() < 0.5);
                mixSampled = true;
                ++mixProposals;
            } else {
            	// System.err.println("ContributionSample");
                if (unseeded > 1 && Math.random() < 0.1) {
                    int a, b;
                    a = b = seeded + (int) Math.floor(Math.random() * unseeded);
                    while (a == b) {
                        b = seeded + (int) Math.floor(Math.random() * unseeded);
                    }
                    newModel.permuteContributions(a, b);
                    double newPrior = prior(newModel);
                    if (newPrior > oldPrior || Math.random() < NativeMath.exp2(newPrior - oldPrior)) {
                        newModel.commit();
                        oldPrior = newPrior;
                    } else {
                        newModel.rollback();
                    }
                }
                
                int g = (int) Math.floor(Math.random() * cgs.length);
                int c = seeded + (int) Math.floor(Math.random() * unseeded);
                ObjectMatrix1D contrib = newModel.getContributions(g);
                PenalizedVariate pv = samplers[g].sample(
                        ((ContributionItem) contrib.get(c)).getItem(), 
                        getUncleVector(g, c)
                );
                // contrib.set(c, new SimpleContributionItem(pv.getVariate()));
                newModel.setContributionMaybeClean(g, c, new SimpleContributionItem(pv.getVariate()), pv.isSilent());
                balancePenalty = pv.getBalancePenalty();
                ++conProposals;
                wasSilent = pv.isSilent();
            }
            
            double newPrior = prior(newModel);
            
            if (Math.random() > NativeMath.exp2(newPrior - oldPrior + balancePenalty)) {
            	// System.err.println("Balance rollback");
                newModel.rollback();
                continue;
            }
            
            double propHood = newModel.likelihood();
            if (propHood < minLikelihood) {
            	// System.err.println("Contour rollback");
                newModel.rollback();
                continue;
            }
            newModel.commit();
            newHood = propHood;
            oldPrior = newPrior;
            
            if (mixSampled) {
                ++mixMoves;
            } else {
            	if (!wasSilent) {
            		++conMoves;
            	}
            }
        }
        
        return newHood;
    }
}
