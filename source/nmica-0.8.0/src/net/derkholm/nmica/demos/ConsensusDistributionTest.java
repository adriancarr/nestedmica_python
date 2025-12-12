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
package net.derkholm.nmica.demos;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;

import net.derkholm.nmica.maths.Gaussian;
import net.derkholm.nmica.seq.consensus.ConsensusDistribution;
import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.dist.Count;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.gui.DNAStyle;
import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.gui.TextLogoPainter;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;

/**
 * Demo program for ConsensusDistribution objects
 *
 * @author Thomas Down
 */

public class ConsensusDistributionTest {
    public static void main(String[] args)
        throws Exception
    {
        ConsensusDistributionTest cdt = new ConsensusDistributionTest();
        CliTools.configureBean(cdt, args);
        cdt.run();
    }
    
    private Count dir = new IndexedCount(DNATools.getDNA());
    private double hardness = 1.0;
    
    public void setA(double d) 
        throws IllegalSymbolException, ChangeVetoException
    {
        dir.setCount(DNATools.a(), d);
    }
    
    public void setC(double d) 
        throws IllegalSymbolException, ChangeVetoException
    {
        dir.setCount(DNATools.c(), d);
    }

    public void setG(double d) 
        throws IllegalSymbolException, ChangeVetoException
    {
        dir.setCount(DNATools.g(), d);
    }

    public void setT(double d) 
        throws IllegalSymbolException, ChangeVetoException
    {
        dir.setCount(DNATools.t(), d);
    }    
    
    public void setHardness(double d) {
        this.hardness = d;
    }
    
    public void run()
        throws Exception
    {
        FiniteAlphabet dna = DNATools.getDNA();
        
        double normSquared = 0;
        for (Iterator si = dna.iterator(); si.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) si.next();
            normSquared += Math.pow(dir.getCount(s), 2);
        }
        if (normSquared > 0.01) {
            for (Iterator si = dna.iterator(); si.hasNext(); ) {
                AtomicSymbol s = (AtomicSymbol) si.next();
                dir.setCount(s, dir.getCount(s) / Math.sqrt(normSquared));
            }
        } else {
            dir.setCount(DNATools.t(), 1.0);
        }
        
        Distribution cd = new ConsensusDistribution(dir, hardness);
        final DistributionLogo dl = new DistributionLogo();
        dl.setBackground(Color.white);
        dl.setOpaque(true);
        dl.setDistribution(cd);
        dl.setPreferredSize(new Dimension(200, 200));
        dl.setLogoPainter(new TextLogoPainter());
        dl.setStyle(new DNAStyle());
        
        JPanel controls = new JPanel();
        controls.setLayout(new FlowLayout());
        final JSlider slider = new JSlider(0, 100, 10);
        slider.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent cev) {
                    hardness = 0.1 * slider.getValue();
                    try {
                        dl.setDistribution(new ConsensusDistribution(dir, hardness));
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            }
        );
        controls.add(slider);
        JButton directionButton = new JButton("Dir");
        directionButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aev) {
                try {
                    dir = sampleDirection(dir);
                    dl.setDistribution(new ConsensusDistribution(dir, hardness));
                } catch (Exception ex) {
                    ex.printStackTrace();
                }   
            }
        } );
        controls.add(directionButton);
        
        JFrame f = new JFrame("ConsensusDistribution");
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(dl, "Center");
        f.getContentPane().add(controls, "North");
        f.pack();
        f.setVisible(true);
    }
    
    private Count sampleDirection(Count c) 
        throws Exception
    {
        List l = new ArrayList(DNATools.createDNA("acgt").toList());
        Collections.shuffle(l);
        AtomicSymbol planeX = (AtomicSymbol) l.get(0);
        AtomicSymbol planeY = (AtomicSymbol) l.get(1);
        double oldX = c.getCount(planeX);
        double oldY = c.getCount(planeY);
        double norm = Math.sqrt(Math.pow(oldX, 2) + Math.pow(oldY, 2));
        if (norm < 0.00001) {
            System.err.println("Sampling on a silly plane");
            return c;
        }
        
        double theta = Math.atan2(oldY, oldX);
        double thetaNew = theta + Gaussian.standardVariate() / Math.sqrt(10);

        double newX = norm * Math.cos(thetaNew);     
        double newY = norm * Math.sin(thetaNew);

        Count newC = new IndexedCount(DNATools.getDNA());
        for (Iterator i = DNATools.getDNA().iterator(); i.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) i.next();
            double val;
            if (s == planeX) {
                val = newX;
            } else if (s == planeY) {
                val = newY;
            } else {
                val = c.getCount(s);
            }
            newC.setCount(s, val);
        }
        return newC;
    }
}
