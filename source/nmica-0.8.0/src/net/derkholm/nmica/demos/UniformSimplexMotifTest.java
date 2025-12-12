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
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Iterator;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import net.derkholm.nmica.gui.WMPanel;
import net.derkholm.nmica.model.ContributionPrior;
import net.derkholm.nmica.model.motif.MotifDropoffPrior;
import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.utils.CliTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 * Demo program for the motif priors and samplers objects
 *
 * @author Thomas Down
 */

public class UniformSimplexMotifTest {
    public static void main(String[] args)
        throws Exception
    {
        UniformSimplexMotifTest cdt = new UniformSimplexMotifTest();
        CliTools.configureBean(cdt, args);
        cdt.run();
    }
    
    private int length = 12;
    
    public void setLength(int l) {
        this.length = l;
    }
    
    private WeightMatrix randomMatrix(FiniteAlphabet alpha, int length)
        throws Exception
    {
        ContributionPrior prior = new MotifDropoffPrior(alpha, length, 3);
        return (WeightMatrix) prior.variate();
    }
    
    private WeightMatrix sampleWmNew(WeightMatrix wm)
        throws Exception
    {
        Distribution[] dists = new Distribution[wm.columns()];
        int victim = (int) Math.floor(Math.random() * dists.length);
        for (int c = 0; c < dists.length; ++c) {
            //if (c == victim) {
                dists[c] = sampleDistNew(wm.getColumn(c));
            // } else {
            //    dists[c] = wm.getColumn(c);
            // }
        }
        return new SimpleWeightMatrix(dists);
    }
    
    private WeightMatrix sampleWmNewBack(WeightMatrix wm)
        throws Exception
    {
        Distribution[] dists = new Distribution[wm.columns()];
        int victim = (int) Math.floor(Math.random() * dists.length);
        for (int c = 0; c < dists.length; ++c) {
            //if (c == victim) {
                dists[c] = sampleDistNewBack(wm.getColumn(c));
            // } else {
            //    dists[c] = wm.getColumn(c);
            // }
        }
        return new SimpleWeightMatrix(dists);
    }
    
    private Distribution sampleDistNew(Distribution dist) 
        throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) dist.getAlphabet();
        Distribution nudist = new NMSimpleDistribution(alpha);
        AtomicSymbol beneficiary = (AtomicSymbol) new UniformDistribution(alpha).sampleSymbol();
        double inc = dist.getWeight(beneficiary) * 0.3;
        for (Iterator si = alpha.iterator(); si.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) si.next();
            double p = dist.getWeight(s);
            double w;
            if (s == beneficiary) {
                w = (p + inc) / (1.0 + inc);
            } else {
                w = p / (1.0 + inc);
            }
            nudist.setWeight(s, w);
        }
        return nudist;
    }
    
    private Distribution sampleDistNewBack(Distribution dist) 
        throws Exception
    {
        double amount = 0.3;
        
        FiniteAlphabet alpha = (FiniteAlphabet) dist.getAlphabet();
        Distribution nudist = new NMSimpleDistribution(alpha);
        AtomicSymbol beneficiary = (AtomicSymbol) new UniformDistribution(alpha).sampleSymbol();
        
        double bFwd = dist.getWeight(beneficiary);                        
        double bBack = bFwd / (1.0 + amount * (1 - bFwd));
        for (Iterator si = alpha.iterator(); si.hasNext(); ) {
            AtomicSymbol s = (AtomicSymbol) si.next();
            double p = dist.getWeight(s);
            double w;
            if (s == beneficiary) {
                w = bBack;
            } else {
                w = p * (1.0 - bBack) / (1.0 - bFwd);
            }
            nudist.setWeight(s, w);
        }
        return nudist;
    }
    
    public void run()
        throws Exception
    {
        final FiniteAlphabet dna = DNATools.getDNA();
        WeightMatrix wm = randomMatrix(dna, length);
        final WMPanel viewer = new WMPanel(wm);
        
        JPanel controls = new JPanel();
        controls.setLayout(new FlowLayout());
        JButton restartButton = new JButton("Restart");
        restartButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aev) {
                try {
                    viewer.setMatrix(randomMatrix(dna, length));
                } catch (Exception ex) {
                    ex.printStackTrace();
                }   
            }
        } );
        controls.add(restartButton);
        
        JButton prefButton = new JButton("Sample F");
        prefButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aev) {
                try {
                    viewer.setMatrix(sampleWmNew(viewer.getMatrix()));
                } catch (Exception ex) {
                    ex.printStackTrace();
                }   
            }
        } );
        controls.add(prefButton);
        
        JButton prefButtonB = new JButton("Sample B");
        prefButtonB.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aev) {
                try {
                    viewer.setMatrix(sampleWmNewBack(viewer.getMatrix()));
                } catch (Exception ex) {
                    ex.printStackTrace();
                }   
            }
        } );
        controls.add(prefButtonB);
        
        JButton runButton = new JButton("Run sampler");
        runButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aev) {
                new Thread() {
                    public void run() {
                while (true) {
                    try {
                        WeightMatrix oldMatrix = viewer.getMatrix();
                        WeightMatrix newMatrix;
                        if (Math.random() < 0.5) {
                            newMatrix = sampleWmNew(oldMatrix);
                        } else {
                            newMatrix = sampleWmNewBack(oldMatrix);
                        }
                        viewer.setMatrix(newMatrix);
                        Thread.sleep(200L);
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
                    }
            }.start();
        }
        });
        controls.add(runButton);
        
        /*
        
        JButton hardButton = new JButton("Sample hardness");
        hardButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aev) {
                try {
                    viewer.setMatrix((WeightMatrix) hardSampler.sample(viewer.getMatrix(), EMPTY_VECTOR));
                } catch (Exception ex) {
                    ex.printStackTrace();
                }   
            }
        } );
        controls.add(hardButton);
        
        */
        
        JFrame f = new JFrame("Motif priors");
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(viewer, "Center");
        f.getContentPane().add(controls, "North");
        f.pack();
        f.setVisible(true);
    }
}
