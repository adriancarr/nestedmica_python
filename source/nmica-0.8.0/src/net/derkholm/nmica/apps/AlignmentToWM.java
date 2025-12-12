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

import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author thomas
 */
public class AlignmentToWM {
    public static void main(String[] args) 
    		throws Exception
    {
        File seqFile = new File(args[0]);
        File outFile = new File(args[1]);
        
        List<Sequence> seqList = new ArrayList<Sequence>();
        SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(seqFile)));
        int maxLength = -1;
        while (si.hasNext()) {
            Sequence seq = si.nextSequence();
            maxLength = Math.max(maxLength, seq.length());
            seqList.add(seq);
        }
        
        WeightMatrix wm = new SimpleWeightMatrix(DNATools.getDNA(), maxLength, DistributionFactory.DEFAULT);
        DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
        for (int c = 0; c < wm.columns(); ++c) {
            dtc.registerDistribution(wm.getColumn(c));
        }
        
        for (Sequence seq: seqList) {
            for (int c = 0; c < seq.length(); ++c) {
                dtc.addCount(wm.getColumn(c), seq.symbolAt(c + 1), 1.0);
            }
        }
        dtc.train();
        
        ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outFile));
        oos.writeObject(wm);
        oos.close();
    }
}
