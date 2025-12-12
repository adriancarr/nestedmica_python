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

import java.util.ConcurrentModificationException;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix1D;
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
 * Trainer which efficiently resamples columns of the mixing matrix by setting
 * them to their maxima then applying bounded monte-carlo.
 * 
 * @author thomas
 */
public class MixtureResamplingTrainer extends Trainer {
    private static final long serialVersionUID = -409176601994689042L;
    
    private int minContributionMoves = 10;
    private int minContributionProposals = 20;
    
    private int mixtureDecopSessions = 5;
    
    private double moveFraction = 0.5;
    private double proposalFraction = 2.0;
    private double permuteProbability = 0.1;
    
    /**
     * @param facetteMap
     * @param data
     * @param components
     * @param priors
     * @param samplers
     * @param mixPolicy
     * @param ensembleSize
     */
    public MixtureResamplingTrainer(FacetteMap facetteMap, Datum[] data,
            int components, ContributionPrior[] priors,
            ContributionSampler[] samplers, MixPolicy mixPolicy,
            int ensembleSize) 
    {
        super(facetteMap, data, components, priors, samplers, mixPolicy, ensembleSize);
    }
    
    public void setPermuteProbability(double d) {
    	this.permuteProbability = d;
    }
    
    public void setMixtureDecopSessions(int s) {
        this.mixtureDecopSessions = s;
    }
    
    public void setMinContributionMoves(int i) {
        minContributionMoves = i;
    }
    
    public int getMinContributionMoves() {
        return minContributionMoves;
    }
    
    public void setMoveFraction(double d) {
        this.moveFraction = d;
    }
    
    public void setProposalFraction(double d) {
        this.proposalFraction = d;
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
        Datum[] data = getDataSet();
        int components = getComponents();
        ContributionGroup[] cgs = getFacetteMap().getContributionGroups();
        ContributionSampler[] samplers = getSamplers();
        MixPolicy mixPolicy = getMixPolicy();
        
        double newHood = newModel.likelihood();
        int conMoves = 0;
        int conProposals = 0;

        for (int mixSession = 0; mixSession < mixtureDecopSessions; ++mixSession) {
            // System.err.printf("MixtureSession %d, O=%g...", mixSession, newHood);
            newHood = resampleMixtureComponent(newModel, (int) (Math.random() * components), minLikelihood);
            // System.err.printf("L=%g...Done!%n", newHood);
        }
        
        // if (Math.abs(newHood - newModel.likelihood()) > 0.01) {
        // 	System.err.printf("resampleMixtureComponent oddity: newHood=%g newModel.likelihood()=%g%n", newHood, newModel.likelihood());
        // }
        
        double oldPrior = prior(newModel);
        // int conNeeded = Math.max(0, minContributionMoves - conMoves) + Math.max(0, minContributionProposals - conProposals);
        while (true) {
            int conNeeded = Math.max(0, minContributionMoves - conMoves) + Math.max(0, minContributionProposals - conProposals);
            if (conNeeded <= 0) {
                break;
            }
            
            // System.err.printf("ContributionSession(%d, %d) O=%g Or=%g...", conMoves, conProposals, newHood, newModel.likelihood());
            double balancePenalty = 0;
            
            if (unseeded > 1 && Math.random() < permuteProbability) {
                int a, b;
                a = b = seeded + (int) (Math.random() * unseeded);
                while (a == b) {
                    b = seeded + (int) (Math.random() * unseeded);
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
            
            int g = 0;
            if (cgs.length > 1) {
            	g = (int) Math.random() * cgs.length;
            }
            int c = seeded + (int) (Math.random() * unseeded);
            ObjectMatrix1D contrib = newModel.getContributions(g);
            ContributionItem oldci = (ContributionItem) contrib.get(c);
            double oldElementPrior = contributionPrior(newModel, g, c);
            PenalizedVariate pv = samplers[g].sample(
                    oldci.getItem(), 
                    getUncleVector(g, c)
            );
            // contrib.set(c, new SimpleContributionItem(pv.getVariate()));
            newModel.setContributionMaybeClean(g, c, new SimpleContributionItem(pv.getVariate(), oldci.getHistoryThread()), pv.isSilent());
            double newElementPrior = contributionPrior(newModel, g, c);
            balancePenalty = pv.getBalancePenalty();
            ++conProposals;
            
            double acceptChance = (newElementPrior - oldElementPrior + balancePenalty);
            if (acceptChance < 0 && Math.random() > NativeMath.exp2(acceptChance)) {
            // if (Math.random() > NativeMath.exp2(newPrior - oldPrior + balancePenalty)) {
            // if (NativeMath.fastlog2(Math.random()) > ) {
                // System.err.printf("Rejecting due to prior%n");
                newModel.rollback();
                continue;
            }
            
            double propHood = newModel.likelihood();
            
            // if (pv.isSilent() && propHood != newHood) {
            //   	System.err.printf("Apparently-silent change wasn't: old=%g new=%g dif=%g sampler=%s%n", newHood, propHood, propHood-newHood, pv.getSampler());
            // }
            
            if (propHood < minLikelihood) {
                // System.err.printf("Rejecting due to likelihood %g%n", propHood);
                newModel.rollback();
                continue;
            }
            newModel.commit();
            newHood = propHood;
            oldPrior = oldPrior - oldElementPrior + newElementPrior;
            
            if (!pv.isSilent()) {
            	++conMoves;
            }
            
            
            // System.err.printf("Success!%n");
        }
        
        
        return newHood;
    }
    
    private double resampleMixtureComponent(TrainableState newModel, int component, double minLikelihood) {
        boolean ignorePrior = getIgnoreMixturePrior();
        
        Datum[] data = getDataSet();
        MixPolicy mixPolicy = getMixPolicy();
        
        newModel.commit();
        
        boolean[] oldState = new boolean[data.length];
        double[] priorRatios = new double[data.length];
        for (int d = 0; d < data.length; ++d) {
            Matrix1D mm = newModel.getMixture(d);
            oldState[d] = mm.get(component) != 0.0;
            mm.set(component, oldState[d] ? 0.0 : 1.0);
            if (!ignorePrior) {
            	double oldPrior = mixPolicy.prior(mm);
                double newPrior = mixPolicy.prior(mm);
                priorRatios[d] = oldState[d] ? oldPrior - newPrior : newPrior - oldPrior;
            }
        }
        newModel.likelihood();
        
        boolean[] flagVector = new boolean[data.length];
        double[] difs = new double[data.length];
        
        double newHood = 0;                           
        for (int d = 0; d < flagVector.length; ++d) {
            double oldCellHood = newModel.getOldDatumLikelihood(d);
            double newCellHood = newModel.getNewDatumLikelihood(d);
            double hoodChange = newCellHood - oldCellHood;
            difs[d] = oldState[d] ? -hoodChange : hoodChange;
            if (difs[d] > 0) {
                flagVector[d] = true;
            }
            newHood += Math.max(oldCellHood, newCellHood);
        }
        
        int stepTarget = (int) Math.ceil(moveFraction * data.length);
        int propTarget = (int) Math.ceil(proposalFraction * data.length);
        int steps = 0;
        int props = 0;
        
      STEP_LOOP:
        while (steps < stepTarget && props < propTarget) {
                  if (Math.random() < 0.2) {
                      try {
	                      int dPlus = pick(flagVector, true);
	                      int dMinus = pick(flagVector, false);
	                      
	                      // this ought to take prior into account, just in case it's odd.
	                      
	                      double stepHood = newHood - difs[dPlus] + difs[dMinus];
	                      double pc = priorRatios[dMinus] - priorRatios[dPlus];
	                      if (stepHood > minLikelihood && (pc >= 0 || Math.random() < NativeMath.exp2(pc))) {
	                          flagVector[dPlus] = false;
	                          flagVector[dMinus] = true;
	                          newHood = stepHood;
	                          ++steps;
	                      }
                      } catch (ZeroCardinalityException ex) {
                    	  // System.err.println("ZCE");
                          continue STEP_LOOP;
                      }
                  } else {
  	                 int d = (int) (Math.random() * flagVector.length);
  	                 boolean current = flagVector[d];
  	                 
  	                 double stepHood;
  	                 double pc;
  	                 if (current) {
  	                     stepHood = newHood - difs[d];
  	                     pc = -priorRatios[d];
  	                 } else {
  	                     stepHood = newHood + difs[d];
  	                     pc = priorRatios[d];
  	                 }
  	                 
  	                 if (stepHood >= minLikelihood && (pc >= 0 || Math.random() < NativeMath.exp2(pc))) {
  	                     newHood = stepHood;
  	                     if (current) {
  	                         flagVector[d] = false;
  	                     } else {
  	                         flagVector[d] = true;
  	                     }
  	                     ++steps;
  	                 }
                  }
                  ++props;
        }
        
        for (int d = 0; d < flagVector.length; ++d) {
        	// System.err.printf("%g\t%s%n", difs[d], flagVector[d] ^ oldState[d]);
            if (!(flagVector[d] ^ oldState[d])) {
                newModel.rollbackDatum(d, component);
            }
        }
        newModel.commit();
        
        return newHood;
    }
        
    private static int pick(boolean[] flagVector, boolean val) 
    		throws ZeroCardinalityException
    {
        int card = 0;
        for (int d = 0; d < flagVector.length; ++d) {
            if (flagVector[d] == val) {
                ++card;
            }
        }
        if (card == 0) {
            throw new ZeroCardinalityException();
        } 
        int hit = (int) (Math.random() * card);
        for (int d = 0; d < flagVector.length; ++d) {
            if (flagVector[d] == val) {
                if (hit == 0) {
                    return d;
                }
                --hit;
            }
        }
        throw new ConcurrentModificationException(String.format("Strange.  card=%d hit=%d", card, hit));
    }
    
    private static class ZeroCardinalityException extends Exception {
		private static final long serialVersionUID = 1L;

		public ZeroCardinalityException() {
            super();
        }
    }
}
