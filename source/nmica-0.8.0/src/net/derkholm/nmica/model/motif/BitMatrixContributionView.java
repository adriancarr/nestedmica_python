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

import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.model.ContributionView;

/**
 * @author thomas
 */

public abstract class BitMatrixContributionView implements ContributionView {
    public static interface Environment {
        public AlphabetIndex getAlphabetIndex();
        public int getMaxIndex();
        public int getPruneLeft(WeightMatrix wm) throws IllegalSymbolException;
        public int getPruneRight(WeightMatrix wm) throws IllegalSymbolException;
        public Class getItemType();
        public WeightMatrix getWeightMatrix(Object item);
    }
    
    public static ContributionView forward(Environment env) {
        return new Forward(env);
    }
    
    public static ContributionView reverse(Environment env) {
        return new Reverse(env);
    }
    
    @SuppressWarnings("unchecked")
	private static class Forward extends BitMatrixContributionView {
        public Forward(Environment env) {
            super(env);
        }
        
        protected double getWeight(WeightMatrix wm, int c, int pruneLeft, int pruneRight, Symbol s)
            throws Exception
        {
            Distribution col = wm.getColumn(c + pruneLeft);
            return col.getWeight(s);
        }
    }
    
    @SuppressWarnings("unchecked")
	private static class Reverse extends BitMatrixContributionView {
        public Reverse(Environment env) {
            super(env);
        }
        
        protected double getWeight(WeightMatrix wm, int c, int pruneLeft, int pruneRight, Symbol s)
            throws Exception
        {
            Distribution col = wm.getColumn(wm.columns() - c - 1 - pruneLeft);
            return col.getWeight(DNATools.complement(s));
        }
    }
    
    protected final Environment facette;
    
    private BitMatrixContributionView(Environment facette) {
        this.facette = facette;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionView#getItemType()
     */
    public Class<?> getItemType() {
        return facette.getItemType();
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionView#getViewType()
     */
    public Class<?> getViewType() {
        return Matrix2D.class;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.model.ContributionView#makeView(java.lang.Object)
     */
    public Object makeView(Object item) {
        try {
	        WeightMatrix wm = facette.getWeightMatrix(item);
	        AlphabetIndex index = facette.getAlphabetIndex();
	        int maxIndex = facette.getMaxIndex();
	        int pruneLeft = facette.getPruneLeft(wm);
	        int pruneRight = facette.getPruneRight(wm);
	        
	        Matrix2D iwm = new SimpleMatrix2D(maxIndex, wm.columns() - pruneLeft - pruneRight);
	        for (int c = 0; c < iwm.columns(); ++c) {
	            // Distribution col = wm.getColumn(c + pruneLeft);
	            for (int i = 0; i < maxIndex; ++i) {
	                iwm.set(i, c, NativeMath.fastlog2(getWeight(wm, c, pruneLeft, pruneRight, index.symbolForIndex(i))));
	            }
	        }
	        return iwm;
        } catch (Exception ex) {
            throw new BioRuntimeException(ex);
        }
    }
    
    protected abstract double getWeight(WeightMatrix wm, int c, int pruneLeft, int pruneRight, Symbol s) throws Exception;
    
 
}
