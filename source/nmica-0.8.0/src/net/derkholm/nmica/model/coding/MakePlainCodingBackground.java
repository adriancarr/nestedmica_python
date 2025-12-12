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

package net.derkholm.nmica.model.coding;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.Reader;

import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.model.coding.PlainCodingSequenceBackground;

import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;


/**
 * 
 * Support class for creating mosaic backgrounds for coding sequences.  Probably
 * needs refactoring. A background model to model exons based on distribution 
 * of amino acids and (DNA * (DNA)^3) codon usage (DeathStar Model).
 * nmica.model.coding
 * @author Bernard Leong
 * @since mica2
 */


public class MakePlainCodingBackground {
    private Reader seqs;
    private OutputStream out = null;
    private InputStream ois = null;
    
    public void setSeqs(Reader r) {
        this.seqs = r;
    }
    
    public void setOut(OutputStream os) {
        this.out = os;
    }
    
    public void run() 
        throws Exception
    {
        SequenceDB seqDB = new HashSequenceDB();
        {
            SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(seqs));
            while (si.hasNext()) {
                seqDB.addSequence(si.nextSequence());
            }
        }
        
        OrderNDistribution plainbackgroundCodonDistribution = Coding.makePlainCodonDistributions(seqDB);
               
        ObjectOutputStream oos = new ObjectOutputStream(out);
        oos.writeObject(new PlainCodingSequenceBackground(plainbackgroundCodonDistribution));
        oos.close();
        
    }

    
    public static void main(String[] args)
        throws Exception
    {
        MakePlainCodingBackground app = new MakePlainCodingBackground();
        CliTools.configureBean(app, args);
        app.run();
    }
}

