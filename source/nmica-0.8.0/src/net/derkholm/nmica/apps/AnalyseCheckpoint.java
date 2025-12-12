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

package net.derkholm.nmica.apps;

import net.derkholm.nmica.model.MultiICAModel;
import net.derkholm.nmica.trainer.Trainer;
import net.derkholm.nmica.utils.CliTools;

import java.io.*;
import java.util.*;

/**
 * @author td2
 *
 */
public class AnalyseCheckpoint {
	private File matrixOut;
	private boolean usageOut = false;
	private File checkpoint;
	private int rank = 1;
	
	public void setMatrixOut(File f) {
		this.matrixOut = f;
	}
	
	public void setUsageOut(boolean b) {
		this.usageOut = b;
	}
	
	public void setCheckpoint(File f) {
		this.checkpoint = f;
	}
	
	public void setRank(int i) {
		this.rank = i;
	}
	
	public void run()
		throws Exception
    {
        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(checkpoint));
        Trainer trainer = (Trainer) ois.readObject();
        ois.close();
        
        List<MultiICAModel> models = new ArrayList<MultiICAModel>(
                Arrays.asList(trainer.getCurrentEnsemble()));
        Collections.sort(models, new Comparator<MultiICAModel>() {
        	public int compare(MultiICAModel m1, MultiICAModel m2) {
        		double dif = m1.likelihood() - m2.likelihood();
        		if (dif < 0) {
        			return -1;
        		} else if (dif > 0) {
        			return 1;
        		} else {
        			return m1.hashCode() - m2.hashCode();
        		}
        	}
        });
        MultiICAModel model;
        int realRank = rank;
        if (rank > 0) {
        	model = (MultiICAModel) models.get(models.size() - rank);
        } else {
        	realRank = models.size() + rank;
        	model = (MultiICAModel) models.get(-rank - 1);
        }
        
        System.err.println("Rank: " + realRank);
        System.err.println("Likelihood: " + model.likelihood());
        
        if (matrixOut != null) {
            ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(matrixOut));
            for (int m = 0; m < model.getComponents(); ++m) {
                oos.writeObject(model.getContribution(model.getFacetteMap().getContributionGroups()[0], m).getItem());
            }
            oos.close();
        }
        
        if (usageOut) {
        	System.out.print("Usage");
        	for (int m = 0; m < model.getComponents(); ++m) {
                System.out.print("\t" + countMix(model, m));
            }
            System.out.println();
        }
    }
	
    private int countMix(MultiICAModel m, int comp) {
        int c = 0;
        int d = m.getDataSet().length;
        for (int i = 0; i < d; ++i) {
            if (m.getMixture(i).get(comp) != 0.0) {
                ++c;
            }
        }
        return c;
    }
	
	public static void main(String[] args) throws Exception {
		AnalyseCheckpoint app = new AnalyseCheckpoint();
		CliTools.configureBean(app, args);
		app.run();
    }
}
