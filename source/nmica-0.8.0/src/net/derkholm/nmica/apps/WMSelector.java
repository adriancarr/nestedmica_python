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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.dp.WeightMatrix;

/**
 * @author thomas
 */
public class WMSelector {
    private File out;
    
    public void setOut(File f) {
        this.out = f;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        WMSelector app = new WMSelector();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(out));
        
        for (int i = 0; i < args.length; i += 2) {
            File mFile = new File(args[i]);
            int index = Integer.parseInt(args[i + 1]);
            
            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(mFile));
            List l = new ArrayList();
            try {
                while (true) {
                    l.add(ois.readObject());
                }
            } catch (Exception ex) {
            }
            WeightMatrix[] motifs = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
            oos.writeObject(motifs[index]);
        }
        
        oos.close();
    }
}
