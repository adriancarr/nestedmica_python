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
/*
 * Created on May 19, 2005
 */
package net.derkholm.nmica.demos;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix1D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix1D;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.model.SimpleContributionItem;
import net.derkholm.nmica.model.motif.MotifFacette;
import net.derkholm.nmica.model.motif.SequenceBackground;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

public class MotifScoringHarness {
    private int repeats = 10000;
    private WeightMatrix[] motifs;
    private SequenceBackground backgroundModel;
    
    public void setBackgroundModel(InputStream is)
        throws Exception
    {
        ObjectInputStream ois = new ObjectInputStream(is);
        backgroundModel = (SequenceBackground) ois.readObject();
        ois.close();
    }
    
    public void setRepeats(int i) {
        this.repeats = i;
    }
    
    public void setMotifs(File f)
        throws Exception
    {
        if (f.getName().endsWith(".jos")) {
            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(f));
            List<WeightMatrix> l = new ArrayList<WeightMatrix>();
            try {
                while (true) {
                    l.add((WeightMatrix) ois.readObject());
                }
            } catch (Exception ex) {
            }
            motifs = l.toArray(new WeightMatrix[0]);
        } else {
            Motif[] d = MotifIOTools.loadMotifSetXML(new FileInputStream(f));
            motifs = new WeightMatrix[d.length];
            for (int m = 0; m < d.length; ++m) {
                motifs[m] = d[m].getWeightMatrix();
            }
        }
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) 
        throws Exception
    {
        MotifScoringHarness app = new MotifScoringHarness();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }

    public void run(String[] args)
        throws Exception
    {
        File seqFile = new File(args[0]);
        
        MotifFacette facette = new MotifFacette(
                backgroundModel,
                0.0,
                true,
                true,
                true,
                1.0,
                false,
                0.0,
                DNATools.getDNA()
            );
        
        List<LikelihoodCalculator> calcs = new ArrayList<LikelihoodCalculator>();
        {
            SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(seqFile)));
            while (si.hasNext()) {
                calcs.add(facette.getLikelihoodCalculator(si.nextSequence()));
            }
        }
        
        for (int c = 0; c < repeats; ++c) {
            ObjectMatrix1D contributions = new SimpleObjectMatrix1D(motifs.length);
            for (int m = 0; m < motifs.length; ++m) {
                contributions.set(m, new SimpleContributionItem(motifs[m]));
            }
            Matrix1D weights = new SimpleMatrix1D(motifs.length, 1.0);
            for (LikelihoodCalculator calc : calcs) {
                calc.likelihood(contributions, weights);
            }
        }
    }
}
