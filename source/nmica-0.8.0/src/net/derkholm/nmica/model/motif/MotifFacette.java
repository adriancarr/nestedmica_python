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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.Iterator;

import net.derkholm.nmica.maths.DoubleFunction;
import net.derkholm.nmica.maths.IdentityDoubleFunction;
import net.derkholm.nmica.model.ContributionView;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Model of sequence motifs against an arbitrary background.
 *
 * @author Thomas Down
 */

public class MotifFacette implements Facette, Serializable {
    static final long serialVersionUID = -3139788280715954374L;
    
    private SequenceBackground background;
    private double softness;
    private boolean nativeDP;
    private boolean uncounted = false;
    private boolean revComp = false;
    @SuppressWarnings("unused")
	private boolean realMixingCoefficients = true;
    private double uncountedExpectation = 1.0;
    private double edgePruneThreshold = 0;
    private double clusterIn = 0.01;
    private double clusterOut = 0.02;
    private boolean cluster = false;
    private FiniteAlphabet alpha;
    private DoubleFunction mixTransferFunction = IdentityDoubleFunction.INSTANCE;
    private boolean discriminate = false;
    
    private int order = 1;
    
    private transient AlphabetIndex index;
    private transient int maxIndex;
    private transient ContributionView forwardBitMatrixView;
    private transient ContributionView reverseBitMatrixView;
 
    public MotifFacette(SequenceBackground background, double softness, boolean nativeDP, boolean uncounted, boolean revComp, double uncountedExpectation, boolean realMixingCoefficients, double edgePruneThreshold, int order, FiniteAlphabet alpha) {
        this.background = background;
        this.softness = softness;
        this.nativeDP = nativeDP;
        this.uncounted = uncounted;
        this.revComp = revComp;
        this.uncountedExpectation = uncountedExpectation;
        this.realMixingCoefficients = realMixingCoefficients;
        this.edgePruneThreshold = edgePruneThreshold;
        this.order = order;
        this.alpha= alpha;
        
        init();
    }
    
    public MotifFacette(SequenceBackground background, double softness, boolean nativeDP, boolean uncounted, boolean revComp, double uncountedExpectation, boolean realMixingCoefficients, double edgePruneThreshold, FiniteAlphabet alpha) {
        this.background = background;
        this.softness = softness;
        this.nativeDP = nativeDP;
        this.uncounted = uncounted;
        this.revComp = revComp;
        this.uncountedExpectation = uncountedExpectation;
        this.realMixingCoefficients = realMixingCoefficients;
        this.edgePruneThreshold = edgePruneThreshold;
        this.order = 1;
	this.alpha = alpha;
        // System.err.println("Initializing a MotifFacette, revComp=" + revComp);
        
        init();
    }
    
    private void init() {
        this.index = AlphabetManager.getAlphabetIndex(alpha);
        this.maxIndex = alpha.size();
        final int MAX_PRUNE = 3;
        
        BitMatrixContributionView.Environment env = new BitMatrixContributionView.Environment() {
            public AlphabetIndex getAlphabetIndex() {
                return index;
            }

            public int getMaxIndex() {
                return maxIndex;
            }

            public int getPruneLeft(WeightMatrix wm) 
		    		throws IllegalSymbolException
		    {
		        double threshold = getEdgePruneThreshold();
		        if (threshold == 0) {
		            return 0;
		        }
		        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
		        
		        for (int p = 0; p < MAX_PRUNE; ++p) {
		            Distribution dist = wm.getColumn(p);
		            for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
		                if (dist.getWeight((Symbol) si.next()) >= threshold) {
		                    return p;
		                }
		            }
		        }
		        
		        return MAX_PRUNE;
		    }
		    
		    public int getPruneRight(WeightMatrix wm) 
		    		throws IllegalSymbolException
		    {
		        double threshold = getEdgePruneThreshold();
		        if (threshold == 0) {
		            return 0;
		        }
		        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
		        
		        for (int p = 0; p < MAX_PRUNE; ++p) {
		            Distribution dist = wm.getColumn(wm.columns() - p - 1);
		            for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
		                if (dist.getWeight((Symbol) si.next()) >= threshold) {
		                    return p;
		                }
		            }
		        }
		        
		        return MAX_PRUNE;
		    }

            public Class getItemType() {
                return WeightMatrix.class;
            }

            public WeightMatrix getWeightMatrix(Object item) {
                return (WeightMatrix) item;
            }
            
        };
        
        forwardBitMatrixView = BitMatrixContributionView.forward(env);
        reverseBitMatrixView = BitMatrixContributionView.reverse(env);
    }
    
    private void readObject(ObjectInputStream stream)
    		throws IOException, ClassNotFoundException
    {
        stream.defaultReadObject();
        init();
    }
    
    public MotifFacette(SequenceBackground background, double softness, boolean nativeDP, boolean uncounted, boolean revComp, FiniteAlphabet alpha) {
        this(background, softness, nativeDP, uncounted, revComp, 1.0, false, 0.0, alpha);
    }
    
    public MotifFacette(SequenceBackground background, double softness, boolean nativeDP, boolean uncounted, FiniteAlphabet alpha) {
        this(background, softness, nativeDP, uncounted, false, alpha);
    }
    
    public boolean isContributionDecoupled(double weight) {
    	if (!discriminate && weight == 0.0) {
    		return true;
    	} else {
    		return false;
    	}
    }
    
    public void setDiscriminate(boolean b) {
    	this.discriminate = b;
    }
    
    public int getOrder() {
        return order;
    }
    
    public DoubleFunction getMixTransferFunction() {
        return mixTransferFunction;
    }
    
    public void setMixTransferFunction(DoubleFunction f) {
        this.mixTransferFunction = f;
    }
    
	/**
	 * @return Returns the cluster.
	 */
	public boolean isCluster() {
		return cluster;
	}
	/**
	 * @param cluster The cluster to set.
	 */
	public void setCluster(boolean cluster) {
		this.cluster = cluster;
	}
	/**
	 * @return Returns the clusterIn.
	 */
	public double getClusterIn() {
		return clusterIn;
	}
	/**
	 * @param clusterIn The clusterIn to set.
	 */
	public void setClusterIn(double clusterIn) {
		this.clusterIn = clusterIn;
	}
	/**
	 * @return Returns the clusterOut.
	 */
	public double getClusterOut() {
		return clusterOut;
	}
	/**
	 * @param clusterOut The clusterOut to set.
	 */
	public void setClusterOut(double clusterOut) {
		this.clusterOut = clusterOut;
	}
    ContributionView getForwardBitMatrixView() {
        return forwardBitMatrixView;
    }
    ContributionView getReverseBitMatrixView() {
        return reverseBitMatrixView;
    }
    AlphabetIndex getAlphabetIndex() {
        return index;
    }
    int getMaxIndex() {
        return maxIndex;
    }
    
    public SequenceBackground getBackground() {
        return background;
    }
    
    public double getEdgePruneThreshold() {
        return edgePruneThreshold;
    }
    
    public double getUncountedExpectation() {
        return uncountedExpectation;
    }
    
    public double getSoftness() {
        return softness;
    }
    
    public boolean getRevComp() {
        return revComp;
    }
    
    public Class getDataType() {
        return SymbolList.class;
    }
    
    public Class getContributionType() {
         return WeightMatrix.class;
    }
    
    public LikelihoodCalculator getLikelihoodCalculator(Object data) {
        try {
        	if (discriminate) {
        		return new MotifDiscriminitiveLikelihood(this, (SymbolList) data);
        	} else if (uncounted) {
            	if (cluster) {
            		if (nativeDP) {
            			return new MotifUncountedClusterLikelihoodNative(this, (SymbolList) data);
            		} else {
            			return new MotifUncountedClusterLikelihood(this, (SymbolList) data);
            		}
            	} else {
            		if (nativeDP) {
            			return new MotifUncountedLikelihoodNative(this, (SymbolList) data);
            		} else {
            			return new MotifUncountedLikelihood(this, (SymbolList) data);
            		}
            	}
            } else {
                // no native counted implementation yet.
                return new MotifLikelihood(this, (SymbolList) data);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            throw new IllegalArgumentException("Couldn't create likelihood calculator");
        }
    }
}

