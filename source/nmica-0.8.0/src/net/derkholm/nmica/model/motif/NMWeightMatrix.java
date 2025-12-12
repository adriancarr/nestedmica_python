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

import java.io.Serializable;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;

/**
 * Weight matrix with possible hidden flanking columns.
 * 
 * <p>
 * This is based on SimpleWeightMatrix from BioJava 1.5
 * </p>
 * 
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class NMWeightMatrix implements WeightMatrix, Serializable {
  private static final long serialVersionUID = 100000L;
    
  private final Distribution [] columns;
  private final Alphabet alpha;
  private final int offset;
  private final int length;

  Distribution[] rawColumnArray() {
	  // Evil but useful...
	  return columns;
  }
  
  public Alphabet getAlphabet() {
    return alpha;
  }

  public int columns() {
    return length;
  }
  
  public Distribution getColumn(int column) {
	column += offset;
	while (column >= columns.length) {
		column -= columns.length;
	}
    return columns[column];
  }
  
  public int unmaskedColumns() {
	  return this.columns.length;
  }
  
  public Distribution getUnmaskedColumn(int column) {
	    return columns[column];
	  }
  
  public int offset() {
	  return offset;
  }
  
  public boolean isHidden(int unmaskedColumn) {
	  int viscol = unmaskedColumn - offset;
	  while (viscol < 0) {
		  viscol += columns.length;
	  }
	  return viscol >= length;
  }
  
  public NMWeightMatrix(Distribution[] columns, int length, int offset)
  throws IllegalAlphabetException {
	if (length < 1 || length > columns.length) {
		throw new IllegalArgumentException("Bad length " + length);
	}
    this.alpha = columns[0].getAlphabet();
    for(int c = 0; c < columns.length; c++) {
      if(columns[c].getAlphabet() != alpha) {
        throw new IllegalAlphabetException(
          "All columns must emit the same alphabet. Expecting " +
          alpha.getName() + ", but found " + columns[c].getAlphabet().getName()
        );
      }
    }
    while (offset > columns.length) {
    	offset -= columns.length;
    }
    while (offset < 0) {
    	offset += columns.length;
    }
    this.columns = columns;
    this.offset = offset;
    this.length = length;
  }
  
  /*
  
  public int hashCode() {
      int hc = 0;
      for (int c = 0; c < columns.length; ++c) {
          hc = (23 * hc) + columns[c].hashCode();
      }
      return hc;
  }
  
  public boolean equals(Object o) {
      if (o instanceof WeightMatrix) {
          WeightMatrix wm = (WeightMatrix) o;
          if (wm.columns() != this.columns()) {
              return false;
          }
          if (wm.getAlphabet() != this.getAlphabet()) {
              return false;
          }
          
          for (int c = 0; c < columns(); ++c) {
              if (! this.getColumn(c).equals(wm.getColumn(c))) {
                  return false;
              }
          }
          return true;
      }
      return false;
  }
  
  */
}
