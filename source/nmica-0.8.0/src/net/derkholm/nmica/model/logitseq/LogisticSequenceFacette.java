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

package net.derkholm.nmica.model.logitseq;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;

import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;

import net.derkholm.nmica.model.ContributionView;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.model.motif.BitMatrixContributionView;

/**
 * @author thomas
 */
public class LogisticSequenceFacette implements Facette, Serializable {
    private static final long serialVersionUID = -2815106384218772101L;
    
    private transient AlphabetIndex index;
    private transient int maxIndex;
    private transient ContributionView forward;
    private transient ContributionView reverse;
    
    public LogisticSequenceFacette() {
        init();
    }
    
    private void readObject(ObjectInputStream stream)
		throws IOException, ClassNotFoundException
	{
		stream.defaultReadObject();
		init();
	}

    AlphabetIndex getAlphabetIndex() {
        return index;
    }
    
    ContributionView forwardBmView() {
        return forward;
    }
    
    ContributionView reverseBmView() {
        return reverse;
    }
    
    private void init() {
        FiniteAlphabet fa = DNATools.getDNA();
        this.index = AlphabetManager.getAlphabetIndex(fa);
        this.maxIndex = fa.size();
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
                return 0;
		    }
		    
		    public int getPruneRight(WeightMatrix wm) 
		    		throws IllegalSymbolException
		    {
		        return 0;
		    }

            public Class getItemType() {
                return WeightedWeightMatrix.class;
            }

            public WeightMatrix getWeightMatrix(Object item) {
                return ((WeightedWeightMatrix) item).getWeightMatrix();
            }
            
        };
        
        forward = BitMatrixContributionView.forward(env);
        reverse = BitMatrixContributionView.reverse(env);
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.Facette#getDataType()
     */
    public Class getDataType() {
        return Sequence.class;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.Facette#getContributionType()
     */
    public Class getContributionType() {
        return WeightedWeightMatrix.class;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.Facette#getLikelihoodCalculator(java.lang.Object)
     */
    public LikelihoodCalculator getLikelihoodCalculator(Object data) {
        try {
            return new LogisticSequenceLikelihoodCalculator(this, (Sequence) data);
        } catch (Exception ex) {
            ex.printStackTrace();
            throw new IllegalArgumentException("Couldn't create likelihood calculator");
        }
    }

	public boolean isContributionDecoupled(double weight) {
		return false;
	}

}
