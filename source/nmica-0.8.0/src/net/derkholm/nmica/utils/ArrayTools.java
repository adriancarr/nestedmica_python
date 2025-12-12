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

package net.derkholm.nmica.utils;

import java.lang.reflect.*;

/**
 * High-level operations on arrays.
 *
 * @author Thomas Down
 */

public class ArrayTools {
    private ArrayTools() {
    }
    
    /**
     * Copy an array.  This method works on both primative and reference arrays
     *
     * @param a An array to copy
     * @return a copy of a
     * @throws IllegalArgumentException if a is not an array
     */
    
    public static Object copy(Object a) {
        Class<?> clazz = a.getClass();
        if (!clazz.isArray()) {
            throw new IllegalArgumentException("This method can only act on arrays.");
        }
        
        Class<?> comp = clazz.getComponentType();
        int length = Array.getLength(a);
        
        Object na = Array.newInstance(comp, length);
        System.arraycopy(a, 0, na, 0, length);
        return na;
    }

    /**
     * @param species
     * @param speciesName
     * @return
     */
    public static int indexOf(Object[] a, Object target) {
        for (int i = 0; i < a.length; ++i) {
            if (a[i] != null && a[i].equals(target)) {
                return i;
            }
        }
        return -1;
    }
}
