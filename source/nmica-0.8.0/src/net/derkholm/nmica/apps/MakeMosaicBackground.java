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
import java.io.File;
import java.io.FileReader;
import java.io.OutputStream;
import java.text.SimpleDateFormat;
import java.util.Date;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Train mosaic background models for use with the NestedMICA motif finder", generateStub=true)
@NMApp(launchName="nmmakebg", vm=VirtualMachine.SERVER)
public class MakeMosaicBackground {

    private File seqs;
    private int mosaicClasses = 1;
    private int mosaicOrder = 1;
    private double mosaicTransition = 0.005;
    private OutputStream out = null;
    private FiniteAlphabet alpha = DNATools.getDNA();
    private boolean match = false;
    
    @Option(help="Initialize background distributions somewhere close to the average composition of the input sequence.  This option is recommended when working with protein sequences (default=false).", optional=true)
    public void setMatch(boolean b) {
    	this.match = b;
    }
    
    @Option(help="Alphabet of the input sequences (default=DNA)", optional=true)
    public void setAlphabet(MicaAlphabet alpha) 
    {
    	this.alpha = alpha.biojavaAlphabet();
    }
    
    @Option(help="A FASTA file of sequences to use when training the background", optional=false)
    public void setSeqs(File r) {
        this.seqs = r;
    }
    
    @Option(help="The number of distinct sequence classes to model", optional=true)
    public void setClasses(int i) {
        this.mosaicClasses = i;
    }
    
    @Option(help="The order of Markov chain to use when modeling background sequences (i.e. a value of one means that the background model will consider dinucleotide frequencies).", optional=true)
    public void setOrder(int i) {
        this.mosaicOrder = i + 1;
    }
    
    public void setMosaicTransition(double d) {
        this.mosaicTransition = d;
    }
    
    @Option(help="Filename for writing the trained background model", optional=false)
    public void setOut(OutputStream os) {
        this.out = os;
    }
    
    public void main(String[] args) 
        throws Exception
    {
        SequenceDB seqDB = new HashSequenceDB();
        {
            SequenceIterator si = SeqIOTools.readFasta(new BufferedReader(new FileReader(seqs)), alpha.getTokenization("token"));
            while (si.hasNext()) {
                seqDB.addSequence(si.nextSequence());
            }
        }
        
        Distribution[] patches = MosaicTools.optimizePatches(
            seqDB,
            mosaicClasses,
            mosaicOrder, 
            mosaicTransition,
            match
        );
 
        Mosaic mosaic = new Mosaic(patches, mosaicTransition);
        mosaic.getAnnotation().setProperty("creator.name", "nmmakebg");
        mosaic.getAnnotation().setProperty("creator.version", Version.VERSION);
        mosaic.getAnnotation().setProperty("input", seqs.getName());
        mosaic.getAnnotation().setProperty("date", new SimpleDateFormat().format(new Date()));
        
        XMLOutputFactory factory = XMLOutputFactory.newInstance();
        XMLStreamWriter xw = factory.createXMLStreamWriter(out);
        MosaicIO.writeMosaic(xw, mosaic);
        xw.close();
    }
}


