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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * General-purpose implementation of ContributionItem, which maintains a
 * thread-safe lazy cache of views.
 * 
 * @author thomas
 */
public class SimpleContributionItem implements ContributionItem, Serializable {
    private static final long serialVersionUID = 10000001L;
    
    private final Object item;
    private final HistoryThread historyThread;
    private transient Map<ContributionView,Object> views = new HashMap<ContributionView,Object>();
    
    private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
        stream.defaultReadObject();
        views = new HashMap<ContributionView,Object>();
    }
    
    public SimpleContributionItem(Object item) {
        this.item = item;
        this.historyThread = new HistoryThread();
    }
    
    public SimpleContributionItem(Object item, HistoryThread historyThread) {
        this.item = item;
        this.historyThread = historyThread;
    }
    
    public final Object getItem() {
        return item;
    }
    
    public final synchronized Object getItemView(ContributionView view) {
        Object viewItem = views.get(view);
        if (viewItem == null) {
            viewItem = view.makeView(item);
            views.put(view, viewItem);
        } 
        return viewItem;
    }
    
    public final HistoryThread getHistoryThread() {
    	return historyThread;
    }
}
