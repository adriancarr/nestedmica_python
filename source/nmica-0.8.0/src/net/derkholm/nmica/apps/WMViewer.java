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

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JFrame;

import net.derkholm.nmica.gui.WMPanel;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.dp.WeightMatrix;

/**
 * Simple weight-matrix viewer
 *
 * @author Thomas Down
 */

public class WMViewer {
    public boolean revComp = false;
    
    public void setRevComp(boolean b) {
        this.revComp = b;
    }
    
    public static void main(String[] args)
    		throws Exception
    	{
        WMViewer app = new WMViewer();
        args = CliTools.configureBean(app, args);
        app.run(args);
    	}
    
    public void run(String[] args)
        throws Exception
    {
        File model = new File(args[0]);
        
        WeightMatrix[] wms;
        {
            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(model));
            List l = new ArrayList();
            try {
                while (true) {
                    l.add(ois.readObject());
                }
            } catch (Exception ex) {
            }
            wms = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
        }
        
        JFrame viewers = new JFrame("WeightMatrixViewer - " + model.getName());
        Box b = Box.createVerticalBox();
        for (int i = 0; i < wms.length; ++i) {
            Box h = Box.createHorizontalBox();
            h.add(new WMPanel(wms[i]));
            if (revComp) {
                h.add(Box.createHorizontalStrut(50));
                h.add(new WMPanel(WmTools.reverseComplement(wms[i])));
            }
            b.add(h);
        }
        viewers.getContentPane().add(b);
        
        viewers.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent wev) {
                System.exit(0);
            }
        } );
        
        viewers.pack();
        viewers.setVisible(true);
    }
}
