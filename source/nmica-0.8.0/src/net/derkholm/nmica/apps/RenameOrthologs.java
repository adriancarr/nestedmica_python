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
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.StringTokenizer;

import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;

/**
 * @author thomas
 */
public class RenameOrthologs {
    private String prefix = "ortho";
    private BufferedReader orthoReader = null;
    
    public void setOrthologList(Reader orthoReader) {
        this.orthoReader = new BufferedReader(orthoReader);
    }
    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        RenameOrthologs app = new RenameOrthologs();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        SequenceDB[] dbs = new SequenceDB[args.length];
        PrintStream[] out = new PrintStream[args.length];
        for (int d = 0; d < args.length; ++d) {
            dbs[d] = loadDB(new File(args[d]));
            out[d] = new PrintStream(new FileOutputStream(args[d] + ".ortho"));
        }
        
        int seed = 0;
        for (String line = orthoReader.readLine(); line != null; line = orthoReader.readLine()) {
            String name = prefix + (++seed);
            StringTokenizer toke = new StringTokenizer(line);
            int d = 0;
            while (toke.hasMoreTokens()) {
                String realName = toke.nextToken();
                Sequence seq = dbs[d].getSequence(realName);
                new FastaFormat().writeSequence(
                        new SimpleSequence(
                                seq,
                                null,
                                name,
                                Annotation.EMPTY_ANNOTATION
                       ),
                       out[d]
                );
                ++d;
            }
            
        }
        
        for (int d = 0; d < out.length; ++d) {
            out[d].close();
        }
    }
    
    private static SequenceDB loadDB(File f)
	    throws Exception
	{
	    SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(f)));
	    SequenceDB seqDB = new HashSequenceDB();
	    while (si.hasNext()) {
	        Sequence seq = si.nextSequence();
	        seqDB.addSequence(seq);
	    }
	    return seqDB;
	}
}
