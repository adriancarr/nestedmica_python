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

package net.derkholm.nmica.seq;

import java.io.*;
import java.util.Map;

import org.biojava.bio.*;
import org.biojava.bio.dist.*;
import org.biojava.bio.symbol.*;
import org.biojava.utils.*;
import org.bjv2.util.SmallMap;

/**
 * A simple implementation of a distribution, which works with any finite alphabet.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Mark Schreiber
 * @since 1.0
 * @serial WARNING serialized versions of this class may not be compatible with later versions of BioJava
 */
public class NMSimpleDistribution
extends AbstractDistribution implements Serializable{
	private static final Map<FiniteAlphabet,AlphabetIndex> LOCKED_INDICES = new SmallMap<FiniteAlphabet, AlphabetIndex>();
	private static final ChangeListener LOCK = new ChangeAdapter() {
        public void preChange(ChangeEvent ce) throws ChangeVetoException {
              throw new ChangeVetoException(
                ce,
                "Can't allow the index to change as we have probabilities."
              );
          }
        };
	
	private static synchronized AlphabetIndex getAndLockIndex(FiniteAlphabet alpha) {
		AlphabetIndex i = LOCKED_INDICES.get(alpha);
		if (i == null) {
			i = AlphabetManager.getAlphabetIndex(alpha);
			i.addChangeListener(LOCK);
			LOCKED_INDICES.put(alpha, i);
		}
		return i;
	}
	
    public static final DistributionFactory FACTORY = new DistributionFactory() {
        public Distribution createDistribution(Alphabet alpha) throws IllegalAlphabetException {
            return new NMSimpleDistribution((FiniteAlphabet) alpha);
        }
    };
    
    static final long serialVersionUID = 7252850540926095729L;
    
  private transient AlphabetIndex indexer;
  private transient double[] weights = null;//because indexer is transient.
  private Distribution nullModel;
  private FiniteAlphabet alpha;
  
  public double[] _getWeightsArray() {
	  return weights;
  }
  
  public int hashCode() {
      int hc = 0;
      double[] weights = getWeights();
      for (int i = 0; i < weights.length; ++i) {
          hc = hc ^ (int) (weights[i] * (1 << 30));
      }
      return hc;
  }
  
  public boolean equals(Object o) {
      if (o instanceof NMSimpleDistribution) {
          NMSimpleDistribution nmd = (NMSimpleDistribution) o;
          if (nmd.getAlphabet() != this.getAlphabet()) {
              return false;
          }
          
          double[] weights = getWeights();
          double[] nmdw = nmd._friendlyGetWeights();
          for (int i = 0; i < weights.length; ++i) {
              if (nmdw[i] != weights[i]) {
                  if (!Double.isNaN(weights[i])) {
                      return false;
                  }
              }
          }
          return true;
      } else if (o instanceof Distribution) {
          return super.equals(o);
      } else {
          return false;
      }
  }
  
  private static class SymbolWeightMemento implements Serializable {
      static final long serialVersionUID = 5223128163879670657L;
      
      public final Symbol symbol;
      public final double weight;
      
      public SymbolWeightMemento(Symbol s, double weight) {
          this.symbol = s;
          this.weight = weight;
      }
  }
  
  private void writeObject(ObjectOutputStream oos)
      throws IOException
  {
      oos.defaultWriteObject();
      
      SymbolWeightMemento[] swm = new SymbolWeightMemento[weights.length];
      for (int w = 0; w < swm.length; ++w) {
          swm[w] = new SymbolWeightMemento(indexer.symbolForIndex(w), weights[w]);
      }
      oos.writeObject(swm);
  }

  private void readObject(ObjectInputStream stream)
    throws IOException, ClassNotFoundException
  {
    stream.defaultReadObject();
    indexer = NMSimpleDistribution.getAndLockIndex(alpha);
    weights = new double[alpha.size()];
    
    SymbolWeightMemento[] swm = (SymbolWeightMemento[]) stream.readObject();
    for (int m = 0; m < swm.length; ++m) {
        try {
            weights[indexer.indexForSymbol(swm[m].symbol)] = swm[m].weight;
        } catch (IllegalSymbolException ex) {
            throw new IOException("Symbol in serialized stream can't be found in the alphabet");
        }
    }
  }

  public Alphabet getAlphabet() {
    return indexer.getAlphabet();
  }

  public Distribution getNullModel() {
    return this.nullModel;
  }



  protected void setNullModelImpl(Distribution nullModel)

  throws IllegalAlphabetException, ChangeVetoException {
    this.nullModel = nullModel;
  }


  /**
   * Indicate whether the weights array has been allocated yet.
   *
   * @return  true if the weights are allocated
   */
  protected boolean hasWeights() {
    return weights != null;
  }


  /**
   * Get the underlying array that stores the weights.
   *
   * <p>
   * Modifying this will modify the state of the distribution.
   * </p>
   *
   * @return  the weights array
   */
  protected double[] getWeights() {
    if(weights == null) {
      weights = new double[((FiniteAlphabet)getAlphabet()).size()];
      for(int i = 0; i < weights.length; i++) {
        weights[i] = Double.NaN;

      }
    }
     return weights;
  }

  double[] _friendlyGetWeights() {
      return getWeights();
  }

  public double getWeightImpl(AtomicSymbol s)

  throws IllegalSymbolException {
    if(!hasWeights()) {
      return Double.NaN;
    } else {
      int index = indexer.indexForSymbol(s);
      return weights[index];
    }
  }


  protected void setWeightImpl(AtomicSymbol s, double w)
  throws IllegalSymbolException, ChangeVetoException {
    double[] weights = getWeights();
    
    /*
    
    if(w < 0.0) {
      throw new IllegalArgumentException(
        "Can't set weight to negative score: " +
        s.getName() + " -> " + w
      );
    }
    
    */
    
    weights[indexer.indexForSymbol(s)] = w;
  }

  private void initialise(FiniteAlphabet alphabet) {
    this.alpha = alphabet;
    this.indexer = NMSimpleDistribution.getAndLockIndex(alpha);

    try {
      setNullModel(DistTools.getUniformDistribution(alphabet));
    } catch (Exception e) {
      throw new BioError("This should never fail. Something is screwed!", e);
    }
  }

  /**
   * make an instance of SimpleDistribution for the specified Alphabet.
   */
  public NMSimpleDistribution(FiniteAlphabet alphabet)
  {
    initialise(alphabet);
  }

  /**
   * make an instance of SimpleDistribution with weights identical
   * to the specified Distribution.
   *
   * @param dist Distribution to copy the weights from.
   */
  public NMSimpleDistribution(Distribution dist)
  {
    try {
    initialise((FiniteAlphabet) dist.getAlphabet());

    // now copy over weights
    int alfaSize = ((FiniteAlphabet)getAlphabet()).size();

    for (int i = 0; i < alfaSize; i++) {
      weights = new double[alfaSize];
      weights[i] = dist.getWeight(indexer.symbolForIndex(i));
    }
    }
    catch (IllegalSymbolException ise) {
      System.err.println("an impossible error surely! "); ise.printStackTrace();
    }
  }

  /**
   * Register an SimpleDistribution.Trainer instance as the trainer for this distribution.
   */
  public void registerWithTrainer(DistributionTrainerContext dtc) {
   dtc.registerTrainer(this, new Trainer());
  }


  /**
   * A simple implementation of a trainer for this class.
   *
   * @author Matthew Pocock
   * @since 1.0
   */
  protected class Trainer implements DistributionTrainer {
    private final Count counts;

    /**
     * Create a new trainer.
     */
    public Trainer() {
      counts = new IndexedCount(indexer);
    }

    public void addCount(DistributionTrainerContext dtc, AtomicSymbol sym, double times)
    throws IllegalSymbolException {
      try {
          counts.increaseCount(sym, times);
      } catch (ChangeVetoException cve) {
        throw new BioError(
          "Assertion Failure: Change to Count object vetoed", cve
        );
      }
    }

    public double getCount(DistributionTrainerContext dtc, AtomicSymbol sym)
    throws IllegalSymbolException {
      return counts.getCount(sym);
    }



    public void clearCounts(DistributionTrainerContext dtc) {
      try {
        int size = ((FiniteAlphabet) counts.getAlphabet()).size();
        for(int i = 0; i < size; i++) {
          counts.zeroCounts();
        }
      } catch (ChangeVetoException cve) {
        throw new BioError(
          "Assertion Failure: Change to Count object vetoed",cve
        );
      }
    }



    public void train(DistributionTrainerContext dtc, double weight)
    throws ChangeVetoException {
      if(!hasListeners(Distribution.WEIGHTS))  {
        trainImpl(dtc, weight);
      } else {
        ChangeSupport changeSupport = getChangeSupport(Distribution.WEIGHTS);
        synchronized(changeSupport) {
          ChangeEvent ce = new ChangeEvent(
            NMSimpleDistribution.this,
            Distribution.WEIGHTS
          );
          changeSupport.firePreChangeEvent(ce);
          trainImpl(dtc, weight);
          changeSupport.firePostChangeEvent(ce);
        }
      }
    }



    protected void trainImpl(DistributionTrainerContext dtc, double weight) {
      //System.out.println("Training");
      try {
        Distribution nullModel = getNullModel();
        double[] weights = getWeights();
        double[] total = new double[weights.length];
        double sum = 0.0;

        for(int i = 0; i < total.length; i++) {
          AtomicSymbol s = (AtomicSymbol) indexer.symbolForIndex(i);
          sum +=
            total[i] =
              getCount(dtc, s) +
              nullModel.getWeight(s) * weight;
        }
        double sum_inv = 1.0 / sum;
        for(int i = 0; i < total.length; i++) {
          //System.out.println("\t" + weights[i] + "\t" + total[i] * sum_inv);
          weights[i] = total[i] * sum_inv;
        }
      } catch (IllegalSymbolException ise) {
        throw new BioError(
          "Assertion Failure: Should be impossible to mess up the symbols.",ise
        );
      }
    }
  }
}



