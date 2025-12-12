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
import net.derkholm.nmica.matrix.ObjectMatrix2D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix2D;
import net.derkholm.nmica.model.Datum;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.FacetteMap;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.utils.WorkQueue;

/**
 * Basic implementation of EvaluationManager which uses in-process work queues to perform likelihood evaluations.
 * 
 * @author thomas
 */
public class LocalEvaluationManager implements EvaluationManager {
    private int workerThreads;
    private ObjectMatrix2D hoodCalcs;
    private transient WorkQueue workQueue;
    
    private TrainableStateContext trainer;
    private TrainableState currentState;
    
    public LocalEvaluationManager() {
        this(1);
    }
    
    public LocalEvaluationManager(int workers) {
        this.workerThreads = workers;
    }
    
    private WorkQueue getWorkQueue() {
        if (workQueue == null) {
            workQueue = WorkQueue.create(workerThreads);
        }
        return workQueue;
    }
    
    private void makeHoodCalcs() {
        FacetteMap fm = trainer.getFacetteMap();
        Facette[] facettes = fm.getFacettes();
        Datum[] data = trainer.getDataSet();
        
        hoodCalcs = new SimpleObjectMatrix2D(facettes.length, trainer.getDataSet().length);
        for (int d = 0; d < data.length; ++d) {
            Object[] fd = data[d].getFacettedData();
            for (int f = 0; f < facettes.length; ++f) {
                if (fd[f] != null) {
                    hoodCalcs.set(f, d, facettes[f].getLikelihoodCalculator(fd[f]));
                }
            }
        }
    }
    
    /* (non-Javadoc)
     * @see net.derkholm.nmica.nstrainer.EvaluationManager#startLikelihoodCalculations(net.derkholm.nmica.nstrainer.TrainableState)
     */
    public void startLikelihoodCalculations(TrainableState state) {
        if (currentState != null) {
            throw new IllegalStateException("Can't start likelihood calculations while queues are busy");
        }
        if (trainer == null) {
            trainer = state.getContext();
            makeHoodCalcs();
        } else if (trainer != state.getContext()) {
            throw new RuntimeException("Already bound to a trainer");
        }
        currentState = state;
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.nstrainer.EvaluationManager#enqueueLikelihoodCalculation(net.derkholm.nmica.nstrainer.TrainableState, int, int, net.derkholm.nmica.maths.DoubleProcedure)
     */
    public void enqueueLikelihoodCalculation(TrainableState state, final int d, final int f, final DoubleProcedure writeback) 
    {
        getWorkQueue().add(new Runnable() {
            public void run() {
                LikelihoodCalculator lc = (LikelihoodCalculator) hoodCalcs.get(f, d);
                double term = lc.likelihood(currentState.getContributions(trainer.facetteIndexToContributionIndex(f)), currentState.getMixture(d));
                writeback.run(term);
            }
            
            public String toString() {
                 return String.format("LikelihoodJob(%s,%d)", trainer.getDataSet()[d].getName(), f);
            }
        } );
    }

    /* (non-Javadoc)
     * @see net.derkholm.nmica.nstrainer.EvaluationManager#endLikelihoodCalculations(net.derkholm.nmica.nstrainer.TrainableState)
     */
    public void endLikelihoodCalculations(TrainableState state) {
        if (workQueue != null) {
            workQueue.flush();
        }
        currentState = null;
    }

}
