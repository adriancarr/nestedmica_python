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

import net.derkholm.nmica.maths.DoubleProcedure;
import net.derkholm.nmica.matrix.CommitableMatrix2D;
import net.derkholm.nmica.matrix.CommitableObjectMatrix2D;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.MatrixTools;
import net.derkholm.nmica.matrix.MatrixWrapper1D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.ObjectMatrix2D;
import net.derkholm.nmica.matrix.ObjectMatrixWrapper1D;
import net.derkholm.nmica.matrix.SimpleCommitableMatrix2D;
import net.derkholm.nmica.matrix.SimpleCommitableObjectMatrix2D;
import net.derkholm.nmica.model.ContributionGroup;
import net.derkholm.nmica.model.ContributionItem;
import net.derkholm.nmica.model.Datum;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.FacetteMap;
import net.derkholm.nmica.model.MultiICAModel;
import net.derkholm.nmica.utils.ArrayTools;
import net.derkholm.nmica.utils.Commitable;

/**
 * Commitable implementation of MultiICAModel, optimized for use in sampling
 * processes.
 *
 * @author Thomas Down
 */

public class TrainableState implements MultiICAModel, Commitable {
    private final TrainableStateContext trainer;
    private CommitableObjectMatrix2D contributions;
    private CommitableMatrix2D mixingMatrix;
    private CommitableMatrix2D hoodCache;
    
    private int[] permute;
    private int[] permuteBack;
    boolean permutationDirty = false;
    
    private int nf;
    
    public TrainableState(TrainableStateContext t) {
        this.trainer = t;
        contributions = new SimpleCommitableObjectMatrix2D(getFacetteMap().getContributionGroups().length, getComponents());
        mixingMatrix = t.getMixPolicy().createCommitableMatrix(getDataSet().length, getComponents());
        hoodCache = new SimpleCommitableMatrix2D(getDataSet().length, getFacetteMap().getFacettes().length, Double.NaN);
        
        permute = new int[getComponents()];
        for (int p = 0; p < permute.length; ++p) {
            permute[p] = p;
        }
        permuteBack = (int[]) ArrayTools.copy(permute);
        nf = t.getFacetteMap().getFacettes().length;
    }
    
    public TrainableState(TrainableStateContext t, MultiICAModel template)
    {
        if (!template.getFacetteMap().equals(t.getFacetteMap())) {
            throw new IllegalArgumentException("Trying to copy-construct a TrainableState in an incompatible context");
        }
        this.trainer = t;
        
        contributions = new SimpleCommitableObjectMatrix2D(getFacetteMap().getContributionGroups().length, getComponents());
        for (int c = 0; c < contributions.rows(); ++c) {
            ObjectMatrix1D tContrib = template.getContributions(getFacetteMap().getContributionGroups()[c]);
            for (int i = 0; i < tContrib.size(); ++i) {
                contributions.set(c, i, tContrib.get(i));
            }
        }
        contributions.commit();
        
        mixingMatrix = t.getMixPolicy().createCommitableMatrix(getDataSet().length, getComponents());
        for (int d = 0; d < getDataSet().length; ++d) {
            Matrix1D tMix = template.getMixture(d);
            for (int i = 0; i < tMix.size(); ++i) {
                mixingMatrix.set(d, i, tMix.get(i));
            }
        }
        mixingMatrix.commit();
      
        hoodCache = new SimpleCommitableMatrix2D(getDataSet().length, getFacetteMap().getFacettes().length, Double.NaN);
        
        permute = new int[getComponents()];
        for (int p = 0; p < permute.length; ++p) {
            permute[p] = p;
        }
        permuteBack = (int[]) ArrayTools.copy(permute);
        nf = t.getFacetteMap().getFacettes().length;
    }
    
    TrainableState(TrainableStateContext t, FrozenModelState fms) {
        this.trainer = t;
        contributions = new SimpleCommitableObjectMatrix2D(fms.contributions);
        {
            Matrix2D fMatrix = fms.mixingMatrix;
            mixingMatrix = t.getMixPolicy().createCommitableMatrix(fMatrix.rows(), fMatrix.columns());
            MatrixTools.copy(mixingMatrix, fMatrix);
        }
        hoodCache = new SimpleCommitableMatrix2D(getDataSet().length, getFacetteMap().getFacettes().length, Double.NaN);
        
        permute = (int[]) ArrayTools.copy(fms.permute);
        permuteBack = (int[]) ArrayTools.copy(fms.permute);
        nf = t.getFacetteMap().getFacettes().length;
    }
    
    void setContributionMaybeClean(int g, int c, ContributionItem item, boolean isClean)
    {
    	c = permute[c];
        contributions.set(g, c, item);
        if (!isClean) {
	        int dl = getDataSet().length;
	        FacetteMap fm = trainer.getFacetteMap();
	        ContributionGroup cg = fm.getContributionGroups()[g];
	        
	        for (int d = 0; d < dl; ++d) {
	            /* if (mixingMatrix.get(d, c) != 0.0) {
	                for (int f = 0; f < hoodCache.columns(); ++f) {
	                    hoodCache.set(d, f, Double.NaN);
	                }
	            } */
	        	
	        	for (int fi = 0; fi < hoodCache.columns(); ++fi) {
	        		Facette f = fm.getFacettes()[fi];
	        		if (fm.contributesToFacette(cg, f) && !f.isContributionDecoupled(mixingMatrix.get(d, c))) {
	        			hoodCache.set(d, fi, Double.NaN);
	        		}
	        	}
	        }
        }
    }
    
    int[] getPermutation() {
        return permute;
    }
    
    public TrainableStateContext getContext() {
        return trainer;
    }
    
    public void permuteContributions(int a, int b) 
    {
        int t = permute[a];
        permute[a] = permute[b];
        permute[b] = t;
        permutationDirty = true;
    }
    
    public FacetteMap getFacetteMap() {
        return trainer.getFacetteMap();
    }
    
    public int getComponents() {
        return trainer.getComponents();
    }
    
    public Datum[] getDataSet() {
        return trainer.getDataSet();
    }
    
    public ContributionItem getContribution(ContributionGroup f, int c) {
        return (ContributionItem) getContribution(trainer.contributionGroupToIndex(f), c);
    }
    
    public ContributionItem getContribution(int f, int c) {
        return (ContributionItem) contributions.get(f, permute[c]);
    }
    
    public ObjectMatrix1D getContributions(int cg) {
        return new TaintingContributionView(cg);
    }
    
    public ObjectMatrix1D getContributions(ContributionGroup cg) {
        return getContributions(trainer.contributionGroupToIndex(cg));
    }
    
    ObjectMatrix2D getContributions() {
        return contributions;
    }
    
    public double likelihood() {
        Datum[] data = getDataSet();
        final Facette[] facettes = getFacetteMap().getFacettes();
        
        EvaluationManager eval = trainer.getEvaluationManager();
        eval.startLikelihoodCalculations(this);
        for (int d = 0; d < data.length; ++d) {
            final Object[] fd = data[d].getFacettedData();
            for (int f = 0; f < facettes.length; ++f) {
	            double term = hoodCache.get(d, f);
	            if (Double.isNaN(term)) {
	                if (fd[f] == null) {
	                    hoodCache.set(d, f, 0.0);
	                } else {
	                    final int upd = d;
	                    final int upf = f;
	                    
	                    eval.enqueueLikelihoodCalculation(
	                            this, 
	                            d,
	                            f,
	                            new DoubleProcedure() {
	                                public void run(double d) {
                                        if (Double.isNaN(d)) {
                                            // System.err.println("Got a NaN");
                                            d = Double.NEGATIVE_INFINITY;
                                        }
	                                    hoodCache.set(upd, upf, d);
	                                }
	                            }
	                    );
	                }
	            }
            }
        }
        eval.endLikelihoodCalculations(this);
        
        double L = 0;
        for (int d = 0; d < data.length; ++d) {
            for (int f = 0; f < facettes.length; ++f) {
                L += hoodCache.get(d, f);
            }
        }
        return L;
    }
    
    double getOldDatumLikelihood(int d) {
        final Facette[] facettes = getFacetteMap().getFacettes();
        double L = 0;
        for (int f = 0; f < facettes.length; ++f) {
            L += hoodCache.getCommitted(d, f);
        }
        return L;
    }
    
    double getNewDatumLikelihood(int d) {
        final Facette[] facettes = getFacetteMap().getFacettes();
        double L = 0;
        for (int f = 0; f < facettes.length; ++f) {
            L += hoodCache.get(d, f);
        }
        return L;
    }

    Matrix2D getMixingMatrix() {
        return mixingMatrix;
    }
    
    public Matrix1D getMixture(int datumIndex) {
        return new TaintingMixtureView(datumIndex);
    }
    
    public Matrix1D getMixtureReadOnly(int datumIndex) {
    	// Doesn't work, permutation is lost :-(
    	
    	// if (mixingMatrix instanceof BinaryCommitableMatrix2D) {
    	//	return ((BinaryCommitableMatrix2D) mixingMatrix).sliceRow(datumIndex);
    	//} else {
    		return getMixture(datumIndex); // TODO: make this RO?
    	//}
    }
    
    public double getMixture(int datumIndex, int contribIndex) {
        return mixingMatrix.get(datumIndex, permute[contribIndex]);
    }
    
    public void commit() {
        contributions.commit();
        ((Commitable) mixingMatrix).commit();
        hoodCache.commit();
        if (permutationDirty) {
            System.arraycopy(permute, 0, permuteBack, 0, permute.length);
            permutationDirty = false;
        }
    }
    
    public void rollback() {
        contributions.rollback();
        ((Commitable) mixingMatrix).rollback();
        hoodCache.rollback();
        if (permutationDirty) {
            System.arraycopy(permuteBack, 0, permute, 0, permute.length);
            permutationDirty = false;
        }
    }
    
    void rollbackDatum(int d) {
        int nc = getComponents();
        for (int c = 0; c < nc; ++c) {
            mixingMatrix.set(d, c, mixingMatrix.getCommitted(d, c));
        }
        for (int f = 0; f < nf; ++f) {
            hoodCache.set(d, f, hoodCache.getCommitted(d, f));
        }
    }
    
    void rollbackDatum(int d, int c) {
        mixingMatrix.set(d, c, mixingMatrix.getCommitted(d, c));
        for (int f = 0; f < nf; ++f) {
            hoodCache.set(d, f, hoodCache.getCommitted(d, f));
        }
    }
    
    public boolean isDirty() {
        return contributions.isDirty() || mixingMatrix.isDirty() || permutationDirty;
    }
    
    private class TaintingMixtureView extends MatrixWrapper1D {
        private final int taintD;
        
        public TaintingMixtureView(int d) {
            super(MatrixTools.viewRow(mixingMatrix, d));
            this.taintD = d;
        }
        
        public void set(int i, double v) {
            super.set(permute[i], v);
            for (int f = 0; f < hoodCache.columns(); ++f) {
                hoodCache.set(taintD, f, Double.NaN);
            }
        }
        
        public double get(int i) {
            return super.get(permute[i]);
        }
    }
    
    private class TaintingContributionView extends ObjectMatrixWrapper1D {
    	private final int g;
    	
        public TaintingContributionView(int g) {
            super(MatrixTools.viewRow(contributions, g));
            this.g = g;
        }
        
        public void set(int c, Object o) {
            setContributionMaybeClean(g, c, (ContributionItem) o, false);
        }
        
        public Object get(int c) {
            return super.get(permute[c]);
        }
    }
}
