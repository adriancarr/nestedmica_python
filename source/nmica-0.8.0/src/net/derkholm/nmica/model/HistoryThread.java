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

package net.derkholm.nmica.model;

import java.io.Serializable;
import java.util.Iterator;

public final class HistoryThread implements Serializable {
	private static final long serialVersionUID = 1000000L;
	
	private final HistoryThread parent;
	private final int id;
	
	public HistoryThread() {
		this.parent = null;
		this.id = newId();
	}
	
	public HistoryThread(HistoryThread parent) {
		this.parent = parent;
		this.id = newId();
	}
	
	public int getId() {
		return id;
	}
	
	public Iterator<HistoryThread> lineage() {
		return new Iterator<HistoryThread>() {
			private HistoryThread state = HistoryThread.this;

			public boolean hasNext() {
				return state != null;
			}

			public HistoryThread next() {
				HistoryThread ht = state;
				state = state.parent;
				return ht;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
	
	private static int idSeed = 0;
	private static synchronized int newId() {
		return ++idSeed;
	}
}
