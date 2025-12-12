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
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Trainer which attempts to efficiently utilize large-scale likelihood calculation facilities
 * by maintaining a queue of pre-evaluated decorrelation ops which can be applied in the future.
 * 
 * @author thomas
 */
public class QueuedDecopTrainer extends Trainer {
    private static final long serialVersionUID = -409176609194689042L;
    
    private int minContributionMoves = 10;
    private int minContributionProposals = 20;
    
    private int mixtureDecopSessions = 2;
    private double mixtureFractionPerSession = 0.5;
    
    /**
     * @param facetteMap
     * @param data
     * @param components
     * @param priors
     * @param samplers
     * @param mixPolicy
     * @param ensembleSize
     */
    public QueuedDecopTrainer(FacetteMap facetteMap, Datum[] data,
            int components, ContributionPrior[] priors,
            ContributionSampler[] samplers, MixPolicy mixPolicy,
            int ensembleSize) 
    {
        super(facetteMap, data, components, priors, samplers, mixPolicy,
                ensembleSize);
    }
    
    public void setMixtureDecopSessions(int s) {
        this.mixtureDecopSessions = s;
    }
    
    public void setMixtureFractionPerSession(double f) {
        this.mixtureFractionPerSession = f;
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
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.trainer.Trainer#decorrelateState(net.derkholm.nmica.trainer.TrainableState, double)
     */
    protected double decorrelateState(TrainableState newModel, double minLikelihood, int seeded, int unseeded) 
    {
        if (seeded > 0) {
            throw new RuntimeException("FIXME");
        }
        
        Datum[] data = getDataSet();
        int components = getComponents();
        ContributionGroup[] cgs = getFacetteMap().getContributionGroups();
        ContributionSampler[] samplers = getSamplers();
        MixPolicy mixPolicy = getMixPolicy();
        
        double newHood = newModel.likelihood();
        int conMoves = 0;
        int conProposals = 0;

        // System.err.println("Starting likelihood is " + newHood);
        
        newModel.commit();
        for (int mixSession = 0; mixSession < mixtureDecopSessions; ++mixSession) {
            List<Integer> proposedOffsets = new ArrayList<Integer>();
            // propose a queue on new mixture decops
            for (int d = 0; d < data.length; ++d) {
                if (Math.random() < mixtureFractionPerSession) {
                    proposedOffsets.add(new Integer(d));
                    
                    Matrix1D mixture = newModel.getMixture(d);
                    double oldDatumPrior = mixPolicy.prior(mixture);
                  SAMPLE_MIX_LOOP:
                    while (true) {
                        mixPolicy.sample(mixture);
                        double newDatumPrior = mixPolicy.prior(mixture);
                        if (newDatumPrior > oldDatumPrior || Math.random() < NativeMath.exp2(newDatumPrior - oldDatumPrior)) {
                            break SAMPLE_MIX_LOOP;
                        }
                        newModel.rollbackDatum(d);
                    }
                }
            }
            
            // force all the dirty likelihoods to be recalculated
            double dummyHood = newModel.likelihood();
            
            Collections.shuffle(proposedOffsets);
            // debug: int accepted = 0;
            for (Integer d: proposedOffsets) {
                double oldCellHood = newModel.getOldDatumLikelihood(d);
                double newCellHood = newModel.getNewDatumLikelihood(d);
                double hoodChange = newCellHood - oldCellHood;
                // System.err.println("HoodChange=" + hoodChange);
                
                double propHood = newHood + hoodChange;
                if (propHood > minLikelihood) {
                    newHood = propHood;
                    // debug: ++accepted;
                } else {
                    newModel.rollbackDatum(d);
                }
            }
            // System.err.println("Accepted " + accepted + "/" + proposedOffsets.size() + " queued decops");
            newModel.commit();
        } 
        
        double oldPrior = prior(newModel);
        int conNeeded = Math.max(0, minContributionMoves - conMoves) + Math.max(0, minContributionProposals - conProposals);
        while (conNeeded > 0) {
            double balancePenalty;
            
            if (components > 1 && Math.random() < 0.1) {
                int a, b;
                a = b = (int) Math.floor(Math.random() * components);
                while (a == b) {
                    b = (int) Math.floor(Math.random() * components);
                }
                newModel.permuteContributions(a, b);
                double newPrior = prior(newModel);
                
                if (newPrior > oldPrior || Math.random() < NativeMath.exp2(newPrior - oldPrior)) {
                    newModel.commit();
                    oldPrior = newPrior;
                    // System.err.println("XO success");
                } else {
                    newModel.rollback();
                    // System.err.println("XO failure");
                }
            }
            
            int g = (int) Math.floor(Math.random() * cgs.length);
            int c = (int) Math.floor(Math.random() * components);
            ObjectMatrix1D contrib = newModel.getContributions(g);
            PenalizedVariate pv = samplers[g].sample(
                    ((ContributionItem) contrib.get(c)).getItem(), 
                    getUncleVector(g, c)
            );
            // contrib.set(c, new SimpleContributionItem(pv.getVariate()));
            newModel.setContributionMaybeClean(g, c, new SimpleContributionItem(pv.getVariate()), pv.isSilent());
            balancePenalty = pv.getBalancePenalty();
            ++conProposals;
            
            double newPrior = prior(newModel);
            
            if (Math.random() > NativeMath.exp2(newPrior - oldPrior + balancePenalty)) {
                newModel.rollback();
                continue;
            }
            
            double propHood = newModel.likelihood();
            if (propHood < minLikelihood) {
                newModel.rollback();
                continue;
            }
            newModel.commit();
            newHood = propHood;
            oldPrior = newPrior;
            
            if (!pv.isSilent()) {
            	++conMoves;
            }
            
            conNeeded = Math.max(0, minContributionMoves - conMoves) + Math.max(0, minContributionProposals - conProposals);
        }
        // System.err.println("Total contribution proposals: " + conProposals);
        
        return newHood;
    }
}
