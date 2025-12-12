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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.Reader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;

import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix1D;
import net.derkholm.nmica.matrix.ObjectMatrix1D;
import net.derkholm.nmica.matrix.SimpleMatrix1D;
import net.derkholm.nmica.matrix.SimpleObjectMatrix1D;
import net.derkholm.nmica.model.Facette;
import net.derkholm.nmica.model.LikelihoodCalculator;
import net.derkholm.nmica.model.SimpleContributionItem;
import net.derkholm.nmica.model.logitseq.LogisticSequenceFacette;
import net.derkholm.nmica.model.logitseq.WeightedWeightMatrix;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.utils.CliTools;

/**
 * @author thomas
 */
public class MoccaClassify {
    private InputStream motifs;
    private BufferedReader seqs;
    private double threshold = 0.5;
    
    public void setMotifs(InputStream motifs) {
        this.motifs = motifs;
    }
    public void setSeqs(Reader seqs) {
        this.seqs = new BufferedReader(seqs);
    }
    public void setThreshold(double threshold) {
        this.threshold = threshold;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        MoccaClassify app = new MoccaClassify();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        Motif[] mots = MotifIOTools.loadMotifSetXML(motifs);
        
        ObjectMatrix1D contributions = new SimpleObjectMatrix1D(mots.length);
        Matrix1D flatMix = new SimpleMatrix1D(mots.length);
        for (int m = 0; m < mots.length; ++m) {
            WeightMatrix wm = mots[m].getWeightMatrix();
            double weight = Double.parseDouble(mots[m].getAnnotation().getProperty("mocca.weight").toString());
            contributions.set(m, new SimpleContributionItem(new WeightedWeightMatrix(wm, weight)));
            flatMix.set(m, 1.0);
        }
        
        Facette logitSeq = new LogisticSequenceFacette();
        
        SequenceIterator si = SeqIOTools.readFastaDNA(seqs);
        while (si.hasNext()) {
            Sequence seq = si.nextSequence();
            
            int label = 0;
            Pattern LABEL_PATTERN = Pattern.compile("label=(1|-1)");
            String descLine = (String) seq.getAnnotation().getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE);
	        Matcher matcher = LABEL_PATTERN.matcher(descLine);
	        if (matcher.find()) {
	            label = Integer.parseInt(matcher.group(1));
	        }
            
            seq.getAnnotation().setProperty("mocca.label", new Integer(1));
            LikelihoodCalculator calc = logitSeq.getLikelihoodCalculator(seq);
            
            System.out.println(seq.getName() + "\t" + label + "\t" + NativeMath.exp2(calc.likelihood(contributions, flatMix)));
        }
    }
}
