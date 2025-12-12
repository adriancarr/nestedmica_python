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

import net.derkholm.nmica.maths.MathsTools;
import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.MatrixTools;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.model.*;
import net.derkholm.nmica.utils.CollectTools;
import org.biojava.bio.BioError;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.*;

/**
 * Generate a series of MultiICA models for the specified dataset using a
 * Nested Sampling strategy.
 * 
 * <p>
 * This is now an abstract base class, providing framework code and the basic
 * nested sampling algorithm.  Concrete subclasses provide specific strategies
 * for decorrelating the states within the NS ensemble.
 * </p>
 *
 * @author Thomas Down
 */

public abstract class Trainer implements Serializable, TrainableStateContext {
    private static final long serialVersionUID = 100000002L;
    
    // Configuration
    
    private final FacetteMap facetteMap;
    private final Datum[] data;
    private final int components;
    private final ContributionPrior[] priors;
    private final ContributionSampler[] samplers;
    private final MixPolicy mixPolicy;
    private final int ensembleSize;
    
    // Seeding
    
    private ContributionItem[][] seedValues = null;
    private transient int[] seedCounts = null;
    
    private transient Facette[] facettes;
    private transient ContributionGroup[] cgs;
    private transient Map<Facette, Integer> facetteIndices;
    private transient Map<ContributionGroup, Integer> contributionIndices;
    private transient int[] facetteIndicesToContributionIndex;
    private transient int[][] contributionIndicesToFacetteIndices;
    
    private double crossOverProb = 0.0;
    private double replaceComponentProb = 0.0;
    private boolean ignoreMixturePrior = false;
    private int samplesToRemove = 1;
    private boolean restrictCrossover = false;
    
    // State
    
    private transient Map<MultiICAModel, Double> modelToLikelihood;
    private int step = 0;
    private int stepsSinceSuccessfulDirectSample = 0;
    private double accumulatedEvidence = Double.NEGATIVE_INFINITY;

    // Evaluation support
    
    private transient EvaluationManager evalManager = new LocalEvaluationManager();
    
    // Monitoring
    
    private transient double[] sortedLikelihoods;
    
    /**
     * A model which has been discarded by the nested sampling training process.
     * The model can be used as a sample which is representative of
     * <code>weight / evidence</code> of the total posterior.
     */
    
    public static class WeightedModel {
        private final MultiICAModel model;
        private final double likelihood;
        private final double priorMass;
        
        WeightedModel(MultiICAModel model, double likelihood, double priorMass) {
            this.model = model;
            this.likelihood = likelihood;
            this.priorMass = priorMass;
        }
        
        public MultiICAModel getModel() {
            return model;
        }
        
        public double getLikelihood() {
            return likelihood;
        }
        
        public double getPriorMass() {
            return priorMass;
        }
        
        public double getWeight() {
            return priorMass + likelihood;
        }
    }
    
    private void makeIndices() {
        this.facettes = facetteMap.getFacettes();
        this.cgs = facetteMap.getContributionGroups();
        
        facetteIndices = new HashMap<Facette, Integer>();
        for (int f = 0; f < facettes.length; ++f) {
            facetteIndices.put(facettes[f], new Integer(f));
        }
        contributionIndices = new HashMap<ContributionGroup, Integer>();
        for (int c = 0; c < cgs.length; ++c) {
            contributionIndices.put(cgs[c], new Integer(c));
        }
        
        facetteIndicesToContributionIndex = new int[facettes.length];
        for (int f = 0; f < facettes.length; ++f) {
            facetteIndicesToContributionIndex[f] = contributionGroupToIndex(facetteMap.getContributionForFacette(facettes[f]));
        }
        
        contributionIndicesToFacetteIndices = new int[cgs.length][];
        for (int g = 0; g < cgs.length; ++g) {
            Facette[] fList = facetteMap.getFacettesForContribution(cgs[g]);
            int[] fIndList = new int[fList.length];
            for (int f = 0; f < fList.length; ++f) {
                fIndList[f] = facetteToIndex(fList[f]);
            }
            contributionIndicesToFacetteIndices[g] = fIndList;
        }
    }
    
    public Trainer(
        FacetteMap facetteMap,
        Datum[] data,
        int components,
        ContributionPrior[] priors,
        ContributionSampler[] samplers,
        MixPolicy mixPolicy,
        int ensembleSize
    ) 
    {
        this.facetteMap = facetteMap;
        this.data = data;
        this.components  = components;
        this.priors = priors;
        this.samplers = samplers;
        this.mixPolicy = mixPolicy;
        this.ensembleSize = ensembleSize;
        
        makeIndices();
        
        seedValues = new ContributionItem[cgs.length][];
    }
    
    public void setSeedContributions(ContributionGroup cg, ContributionItem[] seeds)
        throws IllegalArgumentException
    {
        int ci = contributionGroupToIndex(cg);
        if (ci < 0) {
            throw new IllegalArgumentException("Unknown contribution group " + cg.toString());
        }
        if (seeds.length > components) {
            throw new IllegalArgumentException(String.format("Specified %d seeds for a model with %d components", seeds.length, components));
        }
        seedValues[ci] = seeds;
    }
    
    private void writeObject(ObjectOutputStream oos)
    		throws IOException
    {
        oos.defaultWriteObject();
        TrainableState[] tsa = modelToLikelihood.keySet().toArray(new TrainableState[0]);
        FrozenModelState[] fma = new FrozenModelState[tsa.length];
        for (int m = 0; m < tsa.length; ++m) {
            fma[m] = new FrozenModelState(tsa[m]);
        }
        oos.writeObject(fma);
    }
    
    private void readObject(ObjectInputStream ois)
    		throws IOException, ClassNotFoundException
    {
        ois.defaultReadObject();
        makeIndices();
        evalManager = new LocalEvaluationManager();
        
        // backwards compatibility
        if (seedValues == null) {
            seedValues = new ContributionItem[cgs.length][];
            for (int c = 0; c < cgs.length; ++c) {
                seedValues[c] = new ContributionItem[0];
                seedCounts[c] = 0;
            }
        }
        seedCounts = new int[cgs.length];
        for (int c = 0; c < cgs.length; ++c) {
            seedCounts[c] = seedValues[c].length;
        }
        
        FrozenModelState[] fma = (FrozenModelState[]) ois.readObject();
        modelToLikelihood = new HashMap<MultiICAModel, Double>();
        for (int m = 0; m < fma.length; ++m) {
            TrainableState ts = new TrainableState(this, fma[m]);
            modelToLikelihood.put(ts, new Double(ts.likelihood()));
        }
    }

    public void setRestrictCrossover(boolean b) {
    	this.restrictCrossover = b;
    }
    
    public void setEvaluationManager(EvaluationManager em) {
        this.evalManager = em;
    }
    
    public EvaluationManager getEvaluationManager() {
        return evalManager;
    }
    
    public void setIgnoreMixturePrior(boolean b) {
        this.ignoreMixturePrior = b;
    }
    
    public boolean getIgnoreMixturePrior() {
        return ignoreMixturePrior;
    }
    
    public void setReplaceComponentProb(double d) {
        this.replaceComponentProb = d;
    }
    
    /**
     * @return Returns the samplesToRemove.
     */
    public int getSamplesToRemove() {
        return samplesToRemove;
    }
    
    /**
     * @param samplesToRemove The samplesToRemove to set.
     */
    public void setSamplesToRemove(int samplesToRemove) {
        this.samplesToRemove = samplesToRemove;
    }
    
    public void setCrossOverProb(double d) {
        this.crossOverProb = d;
    }
    
    int facetteToIndex(Facette f)
        throws IllegalArgumentException
    {
        Integer i = facetteIndices.get(f);
        if (i == null) {
            throw new IllegalArgumentException("This model doesn't know anything about " + f.toString());
        }
        return i.intValue();
    }
    
    public int contributionGroupToIndex(ContributionGroup cg)
        throws IllegalArgumentException
    {
        Integer i = contributionIndices.get(cg);
        if (i == null) {
            throw new IllegalArgumentException("This model doesn't know anything about " + cg.toString());
        }
        return i.intValue();
    }
    
    /**
     *  
     * @param fi
     * @return
     */
    
    public int facetteIndexToContributionIndex(int fi) {
        return facetteIndicesToContributionIndex[fi];
    }
    
    int[] contributionIndexToFacetteIndices(int gi) {
        return contributionIndicesToFacetteIndices[gi];
    }
    
    public FacetteMap getFacetteMap() {
        return facetteMap;
    }
    
    public Facette[] getFacettes() {
        return facettes;
    }
    
    public int getComponents() {
        return components;
    }
    
    public Datum[] getDataSet() {
        return data;
    }
    
    public ContributionPrior[] getPriors() {
        return priors;
    }
    
    public ContributionSampler[] getSamplers() {
        return samplers;
    }
    
    public MixPolicy getMixPolicy() {
        return mixPolicy;
    }
    
    protected double contributionPrior(TrainableState ts, int g, int c) {
    	return priors[g].probability(ts.getContribution(g, c).getItem());
    }
    
    protected double prior(TrainableState ts) {
        double L = 0;
        if (!ignoreMixturePrior) {
	        for (int d = 0; d < data.length; ++d) {
	            L += mixPolicy.prior(ts.getMixture(d));
	        }
        }
        for (int g = 0; g < cgs.length; ++g) {
            for (int c = seedCounts[g]; c < components; ++c) {
                L += contributionPrior(ts, g, c);
            }
        }
        if (Double.isNaN(L)) {
            L = Double.NEGATIVE_INFINITY;
        }
        return L;
    }
    
    private TrainableState directSampleModel() {
        TrainableState ts = new TrainableState(this);
        for (int d = 0; d < data.length; ++d) {
            mixPolicy.variate(ts.getMixture(d));
        }
        for (int g = 0; g < cgs.length; ++g) {
            ObjectMatrix1D contrib = ts.getContributions(g); 
            for (int c = 0; c < components; ++c) {
                if (c < seedCounts[g]) {
                    contrib.set(c, seedValues[g][c]);
                } else {
                    contrib.set(c, new SimpleContributionItem(priors[g].variate()));
                }
            }
        }
        ts.commit();
        return ts;
    }
    
    private TrainableState copyState(TrainableState parent) {
        TrainableState ts = new TrainableState(this);
        for (int d = 0; d < data.length; ++d) {
            MatrixTools.copy(ts.getMixture(d), parent.getMixture(d));
        }
        for (int f = 0; f < cgs.length; ++f) {
            MatrixTools.copy(ts.getContributions(f), parent.getContributions(f));
        }
        ts.commit();
        return ts;
    }
    
    private TrainableState crossover(TrainableState parent, TrainableState partner, int partnerTarget) {
    	int parentTarget = partnerTarget;
    	if (restrictCrossover) {
    		int crossInHistoryId = ((ContributionItem) partner.getContributions(0).get(partnerTarget)).getHistoryThread().getId();
    		ObjectMatrix1D src0 = parent.getContributions(0);
        	for (int i = 0; i < src0.size(); ++i) {
        		if (((ContributionItem) src0.get(i)).getHistoryThread().getId() == crossInHistoryId) {
        			parentTarget = i;
        			break;
        		}
        	}
    	}
    	
    	// if (parentTarget != partnerTarget) {
    	// 	System.err.println("Retargetted crossover!");
    	// }
    	
        TrainableState xo = new TrainableState(this);
        for (int d = 0; d < data.length; ++d) {
        	Matrix1D dest = xo.getMixture(d);
        	Matrix1D src0 = parent.getMixture(d);
        	Matrix1D src1 = partner.getMixture(d);
        	for (int i = 0; i < dest.size(); ++i) {
        		double v;
        		if (i == parentTarget) {
        			v = src1.get(partnerTarget);
        		} else {
        			v = src0.get(i);
        		}
                dest.set(i, v);
            }
        }
        for (int f = 0; f < cgs.length; ++f) {
        	ObjectMatrix1D dest = xo.getContributions(f);
        	ObjectMatrix1D src0 = parent.getContributions(f);
        	ObjectMatrix1D src1 = partner.getContributions(f);
        	int crossInHistoryId = ((ContributionItem) src1.get(partnerTarget)).getHistoryThread().getId();
        	boolean historyIdClash = false;
        	for (int i = 0; i < src0.size(); ++i) {
        		if (i != parentTarget && ((ContributionItem) src0.get(i)).getHistoryThread().getId() == crossInHistoryId) {
        			historyIdClash = true;
        		}
        	}
        	 for (int i = 0; i < dest.size(); ++i) {
        		 ContributionItem xci;
        		 if (i == parentTarget) {
        			 xci = (ContributionItem) src1.get(partnerTarget);
        			 if (historyIdClash) {
        				 xci = new SimpleContributionItem(xci.getItem(), new HistoryThread(xci.getHistoryThread()));
        			 }
        		 } else {
        			 xci = (ContributionItem) src0.get(i);
        			 if (xci.getHistoryThread().getId() == crossInHistoryId) {
        				 xci = new SimpleContributionItem(xci.getItem(), new HistoryThread(xci.getHistoryThread()));
        			 }
        		 }
                 dest.set(i, xci);
             }
        }
        xo.commit();
        return xo;
    }
    
    public void init() {
        seedCounts = new int[cgs.length];
        for (int c = 0; c < cgs.length; ++c) {
            if (seedValues[c] == null) {
                seedValues[c] = new ContributionItem[0];
            }
            seedCounts[c] = seedValues[c].length;
        }
        
        modelToLikelihood = new HashMap<MultiICAModel, Double>();
        for (int i = 0; i < ensembleSize; ++i) {
            double hood = Double.NEGATIVE_INFINITY;
            TrainableState ts = null;
            while (hood == Double.NEGATIVE_INFINITY) {
                ts = directSampleModel();
                hood = ts.likelihood();
                // if (hood == Double.NEGATIVE_INFINITY) {
                //    System.err.println("Rejecting improbable sample");
                //}
            }
            modelToLikelihood.put(ts, new Double(hood));
        }
    }
    
    private transient TrainableState[] models;
    
    public double getPriorResidue() {
        return step * NativeMath.log2((1.0 * ensembleSize) / (ensembleSize + 1));
    }
    
    public WeightedModel next() {
        ++step;
        
        // remove the least likely model
        
        double minHood = Double.POSITIVE_INFINITY;
        MultiICAModel minInstance = null;
        for (Map.Entry<MultiICAModel, Double> ensme: modelToLikelihood.entrySet()) {
            double hood = ensme.getValue();
            if (hood < minHood) {
                minHood = hood;
                minInstance = ensme.getKey();
            }
        }
        if (minInstance == null) {
            throw new BioError("Assertion failed: wasn't able to remove a model.  modelToLikelihood.size() = " + modelToLikelihood.size() + " minHood = " + minHood);
        }
        modelToLikelihood.remove(minInstance);
        
        double penalty = -NativeMath.log2(ensembleSize) + step * NativeMath.log2((1.0 * ensembleSize) / (ensembleSize + 1));
        WeightedModel wm = new WeightedModel(new SimpleMultiICAModel(minInstance), minHood, penalty);
        accumulatedEvidence = NativeMath.addLog2(accumulatedEvidence, wm.getWeight());
        
        // generate a new sample
        
        TrainableState newModel = null;
        double newHood = Double.NEGATIVE_INFINITY;
        if (stepsSinceSuccessfulDirectSample < 20) {
            newModel = directSampleModel();
            newHood = newModel.likelihood();
        } 
        
        if (newHood < minHood) {
            int seeded = MathsTools.max(seedCounts);
            int unseeded = components - seeded;
            
            {
                Set<MultiICAModel> s = modelToLikelihood.keySet();
                models = s.toArray(new TrainableState[s.size()]);
            }
            newModel = copyState(models[(int) Math.floor(Math.random() * models.length)]);
            if (!isStateHistoryLegal(newModel)) {
            	System.err.println("Ambigious history before decorrelation");
            }
            if (unseeded > 1) {
                if (Math.random() < replaceComponentProb) {
                    int c = seeded + (int) Math.floor(Math.random() * unseeded);
                    for (int g = 0; g < cgs.length; ++g) {
                        ObjectMatrix1D contrib = newModel.getContributions(g);
                        contrib.set(c, new SimpleContributionItem(priors[g].variate()));
                    }
                    for (int d = 0; d < data.length; ++d) {
                        mixPolicy.sampleComponent(newModel.getMixture(d), c);
                    }
                    newHood = newModel.likelihood();
                    
                    if (newHood >= minHood) {
                        newModel.commit();
                    } else {
                        newModel.rollback();
                    }
                } else if (Math.random() < crossOverProb) {
	                TrainableState partner = models[(int) (Math.random() * models.length)];
	                TrainableState xoNewModel = crossover(newModel, partner, (int) Math.floor(Math.random() * components));
	                newHood = xoNewModel.likelihood();
	                if (newHood >= minHood) {
	                    newModel = xoNewModel;
	                } else {
	                }
	                if (!isStateHistoryLegal(newModel)) {
	                	System.err.println("Ambigious history after XO");
	                }
                }
            }
            newHood = decorrelateState(newModel, minHood, seeded, unseeded);
            if (!isStateHistoryLegal(newModel)) {
            	System.err.println("Ambigious history after decops");
            }
            ++stepsSinceSuccessfulDirectSample;
        } else {
            stepsSinceSuccessfulDirectSample = 0;
        }
	    double newModelLikelihood = newModel.likelihood();
	    if (newModelLikelihood < minHood) {
	      System.err.println("Ugh, likelihood drop " + minHood + " -> " + newModelLikelihood);
	      System.err.println("But we *think* the hood was " + newHood);
	    }
        modelToLikelihood.put(newModel, new Double(newModel.likelihood()));
        
        sortedLikelihoods = null;
        return wm;
    }
    
    protected boolean isStateHistoryLegal(TrainableState state) {
    	
    	/* NO-OP except while debugging
    	
    	Set<Integer> is = new HashSet<Integer>();
    	for (int f = 0; f < cgs.length; ++f) {
    		is.clear();
        	ObjectMatrix1D dest = state.getContributions(f);
        	for (int c = 0; c < dest.size(); ++c) {
        		if(!is.add(((ContributionItem) dest.get(c)).getHistoryThread().getId())) {
        			return false;
        		}
        	}
    	}
    	
    	*/
    	
    	return true;
    }
    
    /**
     * 
     * @param newModel the state to decorrelate
     * @param minLikelihood the currently active likelihood contour.
     * @return the final likelihood of the decorrelated state (which should be >= minLikelihood)
     */
    
    protected abstract double decorrelateState(TrainableState newModel, double minLikelihood, int seeded, int unseeded);
    
    protected ObjectMatrix1D getUncleVector(int g, int c) {
        return new UncleVector(models, g, c);
    }
    
    public void terminate() {
        throw new RuntimeException("Termination conditions aren't implemented yet...");
    }
    
    public double getAccumulatedEvidence() {
        return accumulatedEvidence;
    }
    
    /**
     * Return the current pool of models.  Note that these should not be
     * modified at all if the trainer is still running, otherwise Bad Things
     * are likely to happen.
     */
    
    public MultiICAModel[] getCurrentEnsemble() {
        return modelToLikelihood.keySet().toArray(new MultiICAModel[0]);
    }
    
    /**
     * Return the number of cycles for which this trainer has been running.
     */
    
    public int getCycle() {
        return step;
    }
    
    // Monitoring
    
    private double[] getSortedLikelihoods() {
        if (sortedLikelihoods == null) {
            double[] l = CollectTools.toDoubleArray(modelToLikelihood.values());
            Arrays.sort(l);
            sortedLikelihoods = l;
        }
        return sortedLikelihoods;
    }
    
    /**
     * Get the likelihood of the worst model currently present in the nested ensemble
     */
    
    public double getMinimumLikelihood() {
        return getSortedLikelihoods()[0];
    }
    
    /**
     * Get the likelihood of the best model currently present in the nested ensemble
     */
    
    public double getMaximumLikelihood() {
        double[] l = getSortedLikelihoods();
        return l[l.length - 1];
    }
    
    /**
     * Get the interquartile range of (log-2) likelihoods for the ensemble.
     */
    
    public double getLikelihoodIQR() {
        double[] l = getSortedLikelihoods();
        int loQ = (int) (0.25 * l.length);
        int hiQ = (int) (0.75 * l.length);
        return l[hiQ] - l[loQ];
    }
    
    public MultiICAModel getBestModel() {
        double bestScore = Double.NEGATIVE_INFINITY;
        MultiICAModel bestModel = null;
        for (Map.Entry<MultiICAModel, Double> me: modelToLikelihood.entrySet()) {
            double score = me.getValue();
            if (score > bestScore) {
                bestScore = score;
                bestModel = me.getKey();
            }
        }
        return bestModel;
    }
    
    private class UncleVector implements ObjectMatrix1D {
        private final TrainableState[] models;
        private final int group;
        private final int component;
        
        UncleVector(TrainableState[] models, int group, int component) {
            this.models = models;
            this.group = group;
            this.component = component;
        }
        
        public int size() {
            return models.length;
        }
        
        public Object get(int pos) {
            return models[pos].getContribution(group, component);
        }
        
        public void set(int pos, Object o) {
            throw new UnsupportedOperationException("This is a read-only matrix for accessing other model states");
        }
    }
}
