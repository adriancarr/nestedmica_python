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

import net.derkholm.nmica.gui.WMPanel;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.seq.align.Aligner;
import net.derkholm.nmica.utils.ArrayTools;
import net.derkholm.nmica.utils.CollectTools;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.*;
import org.biojava.bio.dp.onehead.SingleDPMatrix;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.NameTokenization;
import org.biojava.bio.symbol.*;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.MutableTreeNode;
import javax.swing.tree.TreePath;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.GeneralPath;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Browser for viewing motif instances
 * 
 * @author thomas
 */

@App(overview="Show motif hits on a sequence",generateStub=true)
public class MotifExplorer {
    /**
     * @author thomas
     */
    private class AlignMemento {
        public final int species0;
        public final int species1;
        
        /**
         * @param i
         * @param j
         */
        public AlignMemento(int i, int j) {
            this.species0 = i;
            this.species1 = j;
        }
        
        public String toString() {
            return "Align " + species[species0] + " - " + species[species1];
        }

    }
        
    private static final Pattern MOTIF_PATTERN;
    private static final double LOG_2 = Math.log(2);
    private static FiniteAlphabet alpha = DNATools.getDNA();
    
    static {
        try {
	        MOTIF_PATTERN = Pattern.compile("motif([0-9]+)");
        } catch (Exception ex) {
            throw new BioError(ex);
        }
    }
    
    private Color[] palette = new Color[] {
            Color.red,
            Color.green,
            Color.blue,
            new Color(0.0F, 0.5F, 0.0F),   // dark green
            new Color(0.7F, 0.0F, 0.5F),   // grape
            Color.orange,
            new Color(0.4F, 0.0F, 0.6F),   // purple
            new Color(0.5F, 0.7F, 0.25F),  // chartreuse
            Color.cyan,
            new Color(0.6F, 0.3F, 0.1F),   // brown
            new Color(1.0F, 0.6F, 0.5F),    // salmon
            Color.magenta,
            Color.yellow,
            Color.lightGray,
            Color.black
    } ;
    
    /**
     * @author thomas
     */
    private class MotifBeadPanel extends JPanel {
        private Sequence seq = null;
        private Matrix2D mosaic = null;
        
        public MotifBeadPanel() {
            super();
            setPreferredSize(new Dimension(1000, 220));
        }
        
        public void setSequence(Sequence seq) {
            this.seq = seq;
            if (seq.getAnnotation().containsProperty("mex.mosaic")) {
                mosaic = (Matrix2D) seq.getAnnotation().getProperty("mex.mosaic");
            } else {
                mosaic = null;
            }
            repaint();
        }
        
        public Color mixColor(Color from, Color to, double amount) {
      	  float x = (float) amount;
      	  float y = (float) (1.0 - amount);

      	  return new Color((int) (y * from.getRed() + x * to.getRed()),
      			   (int) (y * from.getGreen() + x * to.getGreen()),
      			   (int) (y * from.getBlue() + x * to.getBlue()));
        }

        private Color idColor(double id) {
            return mixColor(Color.yellow, Color.red, Math.max(0, id - 0.25) / 0.75);
        }
        
        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            
            int topBit = 80;
            int width = 960;
            
            if (seq != null) {
                int barH;
                if (mosaic == null) {
                    if (seq.filter(new FeatureFilter.ByType("ungap")).countFeatures() == 0) {
                        barH = 10;
                        g.fillRect(20, topBit, width, 10);
                    } else {
                        barH = 20;
                    }
                } else {
                    g.setColor(Color.black);
                    g.fillRect(20, topBit, width, 10);
                    int classes = mosaic.columns();
                    double scale = 1.0 / (width / seq.length());
                    for (int c = 0; c < classes; ++c) {
                        for (int pos = 0; pos < width; ++pos) {
                            int minSp = 1 + (int) Math.floor(scale * pos);
                            int maxSp = (int) Math.floor(scale * (pos + 1));
                            
                            int points = 0;
                            double tot = 0;
                            for (int spos = minSp; spos <= maxSp; ++spos) {
                                if (spos >= 1 && spos <= mosaic.rows()) {
                                    ++points;
                                    tot += mosaic.get(spos - 1, c);
                                }
                            }
                            g.setColor(mixColor(Color.white, palette[c], tot / points));
                            g.fillRect(20 + pos, topBit + 10 + (8 * c), 1, 8);
                        }
                    }
                    g.setColor(Color.black);
                    g.fillRect(20, topBit + 10 + (8 * classes), width, 10);
                    barH = 20 + (8 * classes);
                }
	            for (Iterator fi = seq.features(); fi.hasNext(); ) {
	                Feature f = (Feature) fi.next();
	                
	                Matcher m = MOTIF_PATTERN.matcher(f.getType());
	                if (m.matches()) {
	                    int w = Integer.parseInt(m.group(1));
	                    
		                int mid = (f.getLocation().getMax() + f.getLocation().getMin()) / 2;
		                int pos = 20 + (int) ((width * mid) / seq.length());
		                int minPos = 20 + (int) ((width * f.getLocation().getMin()) / seq.length());
		                int maxPos = 20 + (int) ((width * f.getLocation().getMax()) / seq.length());
		                if (maxPos - minPos < 6) {
		                    minPos = pos - 3;
		                    maxPos = pos + 3;
		                }
		                
		                // Color c = Color.gray;
		                while (w >= palette.length) {
		                    w -= palette.length;
		                }
		                Color c = palette[w];
		                g.setColor(c);
		                
		                double score = ((Double) f.getAnnotation().getProperty("motif.score")).doubleValue();
		                int height = 10 + (int) Math.floor(score * 5);
		                int barTopPos;
		                int barMidPos;
		                if (!f.getAnnotation().containsProperty("motif.draw_hint") || !"bottom".equals(f.getAnnotation().getProperty("motif.draw_hint"))) {
		                    g.fillRect(minPos, topBit - height, maxPos - minPos, height);
		                    barTopPos = topBit - height - 10;
		                    barMidPos = barTopPos + 3;
		                } else {
		                    g.fillRect(minPos, topBit + barH, maxPos - minPos, height);
		                    barTopPos = topBit + barH + height + 10;
		                    barMidPos = barTopPos - 3;
		                }
		                if (f instanceof StrandedFeature) {
		                    StrandedFeature.Strand strand = ((StrandedFeature) f).getStrand();
		                    GeneralPath arrow = new GeneralPath();
		                    int vector = barTopPos - barMidPos;
		                    if (strand == StrandedFeature.POSITIVE) {
		                        arrow.moveTo(pos - 8, barTopPos);
		                        arrow.lineTo(pos + 3, barTopPos);
		                        arrow.lineTo(pos + 3, barTopPos + vector);
		                        arrow.lineTo(pos + 8, (barTopPos + barMidPos) / 2);
		                        arrow.lineTo(pos + 3, barMidPos - vector);
		                        arrow.lineTo(pos + 3, barMidPos);
		                        arrow.lineTo(pos - 8, barMidPos);
		                        arrow.closePath();
		                    } else {
		                        arrow.moveTo(pos + 8, barTopPos);
		                        arrow.lineTo(pos - 3, barTopPos);
		                        arrow.lineTo(pos - 3, barTopPos + vector);
		                        arrow.lineTo(pos - 8, (barTopPos + barMidPos) / 2);
		                        arrow.lineTo(pos - 3, barMidPos - vector);
		                        arrow.lineTo(pos - 3, barMidPos);
		                        arrow.lineTo(pos + 8, barMidPos);
		                        arrow.closePath();
		                    }
		                    // g.setColor(Color.black);
		                    ((Graphics2D) g).fill(arrow);
		                }
	                } else if ("id".equals(f.getType())) {
	                    int start = 20 + (int) ((width * f.getLocation().getMin()) / seq.length());
	                    int end = 21 + (int) Math.ceil((width * f.getLocation().getMax()) / seq.length());
	                    g.setColor(idColor(((Double) f.getAnnotation().getProperty("id")).doubleValue()));
	                    g.fillRect(start, topBit + 3, Math.max(1, end-start), 14);
	                    // g.setColor(Color.black);
	                    // g.fillRect(start, topBit, Math.max(1, end-start), 3);
	                    // g.fillRect(start, topBit + 17, Math.max(1, end-start), 3);
	                } else if ("gap".equals(f.getType())) {
	                    
	                    /*
	                    
	                    StrandedFeature sf = (StrandedFeature) f;
	                    Color c = Color.black;
	                    int start = 20 + (int) ((460.0 * sf.getLocation().getMin()) / seq.length());
	                    int end = 20 + (int) ((460.0 * sf.getLocation().getMax()) / seq.length());
	                    int height = 20;
	                    g.setColor(c);
		                if (sf.getStrand() == StrandedFeature.POSITIVE) {
		                    g.fillRect(start, topBit - height, Math.max(1, end - start), height);
		                } else {
		                    g.fillRect(start, topBit + barH, Math.max(1, end - start), height);
		                }
		                
		                */
		                
	                } else if ("ungap".equals(f.getType())) {
	                    StrandedFeature.Strand strand = ((StrandedFeature) f).getStrand();
	                    
	                    int start = 20 + (int) ((width * f.getLocation().getMin()) / seq.length());
	                    int end = 20 + (int) Math.ceil((width * f.getLocation().getMax()) / seq.length());
	                    g.setColor(Color.black);
	                    if (strand == StrandedFeature.POSITIVE) {
	                        g.fillRect(start, topBit, Math.max(1, end-start), 3);
	                    }
	                    if (strand == StrandedFeature.NEGATIVE)  {
	                        g.fillRect(start, topBit + 17, Math.max(1, end-start), 3);
	                    }
	                }else if ("annotation".equals(f.getSource())) {
	                    String type = "" + f.getType();
	                    int typeWidth = g.getFontMetrics().stringWidth(type);
	                    // int mid = (f.getLocation().getMax() + f.getLocation().getMin()) / 2;
		                int minPos = 20 + (int) ((width * f.getLocation().getMin()) / seq.length());
		                int maxPos = 20 + (int) ((width * f.getLocation().getMax()) / seq.length());
		                int midPos = (minPos + maxPos) / 2;
		                g.setColor(Color.red);
		                g.fillRect(minPos - 3, topBit + barH, maxPos - minPos, 40);
		                g.setColor(Color.black);
		                g.drawString("" + f.getType(), midPos - typeWidth / 2, topBit + 55 + barH);
	                }
	            }
	            
	            g.setColor(Color.black);
	            int segLength = (int) (width / 2 - 100);
	            g.drawLine(20, topBit + 80 + barH, 20 + segLength, topBit + 80 + barH);
	            g.drawLine(20, topBit + 80 + barH, 40, topBit + 70 + barH);
	            g.drawLine(20, topBit + 80  + barH, 40, topBit + 90 + barH);
	            
	            g.drawLine(980 - segLength, topBit + 80 + barH, 980, topBit + 80 + barH);
	            g.drawLine(980, topBit + 80 + barH, 960, topBit + 70 + barH);
	            g.drawLine(980, topBit + 80 + barH, 960, topBit + 90 + barH);
	            
	            g.drawString("" + seq.length() + "bp", 470, topBit + 85 + barH);

	            if (seq.filter(new FeatureFilter.ByType("id")).countFeatures() > 0) {
		            g.drawString("%id", 5, topBit + 145);
		            int[] keyIds = new int[] {25, 40, 55, 70, 85, 100};
		            for (int k = 0; k < keyIds.length; ++k) {
		                double id = (1.0 * keyIds[k]) / 100;
		                Color c = idColor(id);
		                g.setColor(c);
		                g.fillRect(40 + k * 40, topBit + 120, 40, 40);
		                g.setColor(Color.black);
		                g.drawRect(40 + k * 40, topBit + 120, 40, 40);
		                g.drawString("" + keyIds[k] + "%", 45 + k * 40, topBit + 145);
		            }
	            }
            }
        }
    }
    
    private class FixedListCellRenderer extends DefaultListCellRenderer {
        public Component getListCellRendererComponent(
            JList list,
        		Object value,
                int index,
                boolean isSelected,
                boolean cellHasFocus)
            {


        	    setText((value == null) ? "" : value.toString());

        	setEnabled(list.isEnabled());
        	setFont(list.getFont());

        	return this;
            }
    
        public FixedListCellRenderer() {
            super();
        }
        
        public boolean isOpaque() {
            return true;
        }
    }
    
    private class MyCellRenderer extends JLabel implements ListCellRenderer {
        public MyCellRenderer() {
            setOpaque(true);
        }
        public Component getListCellRendererComponent(
            JList list,
            Object value,
            int index,
            boolean isSelected,
            boolean cellHasFocus)
        {
            setText(value.toString());
            Color c = palette[index % palette.length];
            setForeground(c);
            return this;
        }
    }
    
    private String[] species;
    private SequenceDB[] seqs;
    private WeightMatrix[] motifs;
    private double scoreThreshold = 2;
    private File cache = null;
    private double softness = 0.0;
    private String alignerClass = "net.derkholm.nmica.seq.align.SimpleNucleicAcidAligner";
    private MosaicSequenceBackground mosaic = null;
    private File[] annotations = new File[0];
    private File[] gffAnnotations = new File[0];
    private double[] thresholdList;
    
    private Aligner aligner;
    
    public void setThresholdList(Reader r)
    		throws IOException
    	{
        BufferedReader br = new BufferedReader(r);
        List<Double> l = new ArrayList<Double>();
        for (String line = br.readLine(); line != null; line = br.readLine()) {
            l.add(new Double(line));
        }
        thresholdList = CollectTools.toDoubleArray(l);
    	}
    
    public void setAnnotations(File[] annoFiles) {
        this.annotations = annoFiles;
    }
    
    public void setGffAnnotations(File[] annoFiles) {
        this.gffAnnotations = annoFiles;
    }
    
    public void setMosaic(InputStream is) 
    		throws Exception
    {
        ObjectInputStream ois = new ObjectInputStream(is);
        mosaic = (MosaicSequenceBackground) ois.readObject();
        ois.close();
    }
    
    @Option(help="Name of Java class to use as alignment helper program", optional=true)
    public void setAlignerClass(String s) {
        this.alignerClass = s;
    }
    
    @Option(help="Cut-off threshold for showing motifs", optional=true)
    public void setScoreThreshold(double d) {
        this.scoreThreshold = d;
    }
    
    @Option(help="Species names", optional=true)
    public void setSpecies(String [] species) {
        this.species = species;
    }
    
    @Option(help="A FASTA file of sequences to view", optional=true)
    public void setSeqs(Reader[] r) throws Exception {
        seqs = new SequenceDB[r.length];
		    System.out.println("seq part "+alpha.getName());
        for (int s = 0; s < r.length; ++s) {
	        seqs[s] = new HashSequenceDB();
	        SequenceIterator si = SeqIOTools.readFasta(new BufferedReader(r[s]),alpha.getTokenization("token"));
		Sequence seq;
	        while (si.hasNext()) {
		    seq = si.nextSequence();
	            seqs[s].addSequence(seq);
	        }
        }
    }
    
    @Option(help="A file of motifs to show", optional=false)
    public void setMotifs(File f)
		throws Exception
	{
		if (f.getName().endsWith(".jos")) {
		    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(f));
		    List l = new ArrayList();
		    try {
		        while (true) {
		            l.add(ois.readObject());
		        }
		    } catch (Exception ex) {
		    }
		    motifs = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
		} else {
		    Motif[] d = MotifIOTools.loadMotifSetXML(new FileInputStream(f));
		    motifs = new WeightMatrix[d.length];
		    for (int m = 0; m < d.length; ++m) {
		        motifs[m] = d[m].getWeightMatrix();
		    }
		    System.out.println(alpha.getName());
		}
		alpha = (FiniteAlphabet) motifs[0].getAlphabet();
	}

    /**
     * 
     */
    public void main(String[] args) 
    		throws Exception
    {
        if (seqs == null) {
            throw new Exception("No sequences");
        }
        
        if (motifs == null) {
            throw new Exception("No motifs");
        }
        
        {
            File tmpCache = new File("me-cache");
            if (tmpCache.exists()) {
                cache = tmpCache;
            }
        }
        
        {
            for (int f = 0; f < annotations.length; ++f) {
                BufferedReader br = new BufferedReader(new FileReader(annotations[f]));
                for (String line = br.readLine(); line != null; line = br.readLine()) {
                    StringTokenizer toke = new StringTokenizer(line);
                    String seqName = toke.nextToken();
                    int min = Integer.parseInt(toke.nextToken());
                    int max = Integer.parseInt(toke.nextToken());
                    String type = toke.nextToken();
                    
                    Sequence seq = seqs[f].getSequence(seqName);
                    Feature.Template ft = new Feature.Template();
                    ft.type = type;
                    ft.source = "annotation";
                    ft.location = new RangeLocation(min, max);
                    ft.annotation = new SmallAnnotation();
                    seq.createFeature(ft);
                }
            }
        }
        
        {
            for (int f = 0; f < gffAnnotations.length; ++f) {
                GFFEntrySet gff = GFFTools.readGFF(gffAnnotations[f]);
                for (Iterator i = gff.lineIterator(); i.hasNext(); ) {
                    Object o = i.next();
                    if (o instanceof GFFRecord) {
                        GFFRecord record = (GFFRecord) o;
	                    String seqName = record.getSeqName();
	                    int min = record.getStart();
	                    int max = record.getEnd();
	                    String type = record.getFeature();
	                    
	                    Map gaga = record.getGroupAttributes();
	                    if (gaga.containsKey("Factor")) {
	                        // System.err.println("Got factor information");
	                        type = ((List) gaga.get("Factor")).get(0).toString();
	                    }
	                    
	                    Sequence seq = seqs[f].getSequence(seqName);
	                    Feature.Template ft = new Feature.Template();
	                    ft.type = type;
	                    ft.source = "annotation";
	                    ft.location = new RangeLocation(min, max);
	                    ft.annotation = new SmallAnnotation();
	                    seq.createFeature(ft);
                    }
                }
            }
        }
        
        aligner = getClass()
                .getClassLoader()
                .loadClass(alignerClass)
                .asSubclass(Aligner.class)
                .newInstance();
        
        JFrame explorerFrame = new JFrame("MotifExplorer");
        
        Box mainArea = Box.createVerticalBox();
        final JToolBar infoBar = new JToolBar();
        infoBar.setFloatable(false);
        infoBar.setBorderPainted(true);
        final JLabel statusBar = new JLabel("MotifExplorer");
        infoBar.add(statusBar);
        // final Action cancelAction = new AbstractAction("Cancel") {
        //     public void actionPerformed(ActionEvent e) {
                // TODO Auto-generated method stub
                
        //    }
        // };
        // infoBar.add(cancelAction);
        mainArea.add(infoBar);
        
        final MotifBeadPanel seqPanel = new MotifBeadPanel();
        mainArea.add(seqPanel);
        
        final JTree seqList;
        {
            Set allSeqIds = new TreeSet();
            for (int s = 0; s < seqs.length; ++s) {
                allSeqIds.addAll(seqs[s].ids());
            }
            
            DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
            for (Iterator idi = allSeqIds.iterator(); idi.hasNext(); ) {
                String id = (String) idi.next();
                
                // DefaultMutableTreeNode idNode = new DefaultMutableTreeNode(id);
                List<Integer> alignIdList = new ArrayList<Integer>();
                for (int s = 0; s < seqs.length; ++s) {
                    if (seqs[s].ids().contains(id)) {
                        // MutableTreeNode seqNode = new DefaultMutableTreeNode(species[s]);
                        // idNode.add(seqNode);
                        alignIdList.add(new Integer(s));
                    }
                }
                int[] alignIds = CollectTools.toIntArray(alignIdList);
                if (alignIds.length == 1) {
                    DefaultMutableTreeNode idNode = new DefaultMutableTreeNode(id);
                    rootNode.add(idNode);
                } else {
                    DefaultMutableTreeNode idNode = new DefaultMutableTreeNode(id);
                    for (int i = 0; i < alignIds.length; ++i) {
                        MutableTreeNode seqNode = new DefaultMutableTreeNode(species[alignIds[i]]);
                        idNode.add(seqNode);
                    }
	                for (int i = 0; i < alignIds.length; ++i) {
	                    for (int j = i + 1; j < alignIds.length; ++j) {
	                        MutableTreeNode alignNode = new DefaultMutableTreeNode(new AlignMemento(i, j));
	                        idNode.add(alignNode);
	                    }
	                }
	                rootNode.add(idNode);
                }
            }
            
            seqList = new JTree(rootNode);
            seqList.setRootVisible(false);
            
            seqList.addMouseListener(new MouseAdapter() {
                public void mouseClicked(MouseEvent mev) {
                    if (mev.getClickCount() == 2) {
                        TreePath tpath = seqList.getPathForLocation(mev.getPoint().x, mev.getPoint().y);
                        Object[] path = tpath.getPath();
                        if (path.length == 2) {
                            String id = (String) ((DefaultMutableTreeNode) path[1]).getUserObject();
                            for (int s = 0; s < seqs.length; ++s) {
                                if (seqs[s].ids().contains(id)) {
                                    try {
	    	                                seqPanel.setSequence(makeSequenceExplorer(seqs[s].getSequence(id)));
	    	                            } catch (Exception e) {
	    	                                e.printStackTrace();
	    	                            }
	    	                            break;
                                }
                            }
                        } else if (path.length == 3) {
                            String id = (String) ((DefaultMutableTreeNode) path[1]).getUserObject();
                            Object selector = ((DefaultMutableTreeNode) path[2]).getUserObject();
                            if (selector instanceof String) {
                                String speciesName = (String) selector;
	                            int speciesIdx = ArrayTools.indexOf(species, speciesName);
	                            try {
	                                seqPanel.setSequence(makeSequenceExplorer(seqs[speciesIdx].getSequence(id)));
	                            } catch (Exception e) {
	                                e.printStackTrace();
	                            }
                            } else if (selector instanceof AlignMemento) {
                                AlignMemento am = (AlignMemento) selector;
                                try {
                                    seqPanel.setSequence(makeAlignExplorer(id, am.species0, am.species1));
                                } catch (Exception e) {
	                                e.printStackTrace();
	                            }
                            }
                            
                            statusBar.setText(id + ":" + selector.toString());
                        }
                        
                    }
                }
            });
        }
        
        final JList motifList;
        {
            Vector<String> motifVector = new Vector<String>();
            for (int w = 0; w < motifs.length; ++w) {
                motifVector.add(WmTools.fluffyConsensus(motifs[w]) + " / " + WmTools.fluffyConsensus(WmTools.reverseComplement(motifs[w])));
            }
            motifList = new JList(motifVector);
            motifList.setCellRenderer(new MyCellRenderer());
            
            motifList.addMouseListener(new MouseAdapter() {
                public void mouseClicked(MouseEvent mev) {
                    if (mev.getClickCount() == 2) {
                        int index = motifList.locationToIndex(mev.getPoint());
                        try {
                            WMPanel.wmViewer(motifs[index], "Motif " + index);
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
            });
        }
        
        JSplitPane sideBar = new JSplitPane(JSplitPane.VERTICAL_SPLIT, motifList, new JScrollPane(seqList));
        JSplitPane mainLayout = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, sideBar, mainArea);
        explorerFrame.getContentPane().add(mainLayout);
        explorerFrame.pack();
        
        explorerFrame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent wev) {
                System.exit(0);
            }
        } );
        
        explorerFrame.setVisible(true);
    }

    /**
     * @param sequence
     * @throws Exception
     */
    private Sequence makeSequenceExplorer(Sequence sequence) throws Exception {
        Sequence mView = SequenceTools.view(sequence, sequence.getName() + "-me");
        annotateMotifs(mView, mView, StrandedFeature.UNKNOWN, "top");
        if (mosaic != null) {
            SymbolList msl = mView;
            int order = mosaic.getMosaicOrder();
            int classes = mosaic.getMosaicClasses();
            if (order > 1) {
                msl = SymbolListViews.orderNSymbolList(msl, order);
            }
            DP dp = mosaic.getBackgroundDP();
            
            SingleDPMatrix forwardMatrix = (SingleDPMatrix) dp.forwardMatrix(
                    new SymbolList[] {msl},
                    ScoreType.PROBABILITY
            );
            double score = forwardMatrix.getScore();
            SingleDPMatrix backwardMatrix = (SingleDPMatrix) dp.backwardMatrix(
                    new SymbolList[] {msl},
                    ScoreType.PROBABILITY
            );
            
            Matrix2D display = new SimpleMatrix2D(msl.length(), classes);
            {
                Pattern namePattern = Pattern.compile("patch([0-9]+)");
                
                for (int pos = 1; pos < msl.length(); ++pos) {
                    State[] states = forwardMatrix.states();
                    double[] fcol = forwardMatrix.scores[pos];
                    double[] bcol = backwardMatrix.scores[pos];
                    for (int s = 0; s < states.length; ++s) {
                        Matcher nameMatcher = namePattern.matcher(states[s].getName());
                        if (nameMatcher.matches()) {
                            int p = Integer.parseInt(nameMatcher.group(1));
                            display.set(pos - 1, p, Math.exp(fcol[s] + bcol[s] - score));
                        }
                    }  
                }
            }
            mView.getAnnotation().setProperty("mex.mosaic", display);
        }
        return mView;
    }
    
    /*
    
    private WeightMatrix wmToBm(WeightMatrix wm)
    		throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
        
        WeightMatrix bm = new SimpleWeightMatrix(wm.getAlphabet(), wm.columns(), DistributionFactory.DEFAULT);
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution wmCol = wm.getColumn(c);
            Distribution bmCol = bm.getColumn(c);
            
            double info = info(wmCol);
            for (Iterator si = alpha.iterator(); si.hasNext(); ) {
                Symbol s = (Symbol) si.next();
                bmCol.setWeight(s, info * wmCol.getWeight(s));
            }
        }
        return bm;
    }
    
    */
    
    private WeightMatrix wmToBm(WeightMatrix wm)
    		throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
        
        WeightMatrix bm = new SimpleWeightMatrix(wm.getAlphabet(), wm.columns(), NMSimpleDistribution.FACTORY);
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution wmCol = wm.getColumn(c);
            Distribution bmCol = bm.getColumn(c);
            
            for (Iterator si = alpha.iterator(); si.hasNext(); ) {
                Symbol s = (Symbol) si.next();
                double baseWeight = wmCol.getWeight(s);
                double moderatedWeight = (softness / 4) + baseWeight * (1.0 - softness);
                bmCol.setWeight(s, Math.log(moderatedWeight) / LOG_2);
            }
        }
        return bm;
    }
    
    private double info(Distribution dist)
    		throws Exception
    {
        double inf = 0;
        for (Iterator i = ((FiniteAlphabet) dist.getAlphabet()).iterator(); i.hasNext(); ) {
            Symbol s = (Symbol) i.next();
            double w = dist.getWeight(s);
            if (w > 0) {
                inf += w * Math.log(1.0 / w) / LOG_2;
            }
        }
        return 2.0 - inf;
    }
    
    private void annotateMotifs(SymbolList sl, Sequence seq, StrandedFeature.Strand strand, String posKey)
    		throws Exception
    {
        for (int w = 0; w < motifs.length; ++w) {
            double threshold = (thresholdList == null ? scoreThreshold : thresholdList[w]);
            WeightMatrix fwdMatrix = wmToBm(motifs[w]);
            WeightMatrix backMatrix = wmToBm(WmTools.reverseComplement(motifs[w]));
            int maxPos = sl.length() - motifs[w].columns() + 1;
            double maxScore = bmMaxScore(fwdMatrix);
            for (int pos = 1; pos <= maxPos; ++pos) {
                {
                    // double score = scoreWM(sl, motifs[w], pos);
                    // double bitsSubOptimal = -Math.log(score / wmMaxScore(motifs[w])) / LOG_2;
                    double bitsSubOptimal = maxScore - scoreBM(sl, fwdMatrix, pos);
                    
                    if (bitsSubOptimal <= threshold) {
                        StrandedFeature.Template template = new StrandedFeature.Template();
                        template.type = "motif" + w;
	                    template.source = "mica";
	                    template.location = new RangeLocation(pos, pos + motifs[w].columns() - 1);
	                    template.strand = strand == StrandedFeature.UNKNOWN ? StrandedFeature.POSITIVE : strand;
	                    template.annotation = new SmallAnnotation();
	                    template.annotation.setProperty("motif.score", new Double(threshold - bitsSubOptimal));
	                    template.annotation.setProperty("motif.draw_hint", posKey);
	                    seq.createFeature(template);
	                }
                }
                {
                    // double score = scoreWM(sl, backMatrix, pos);
                    // double bitsSubOptimal = -Math.log(score / wmMaxScore(motifs[w])) / LOG_2;
                    double bitsSubOptimal = maxScore - scoreBM(sl, backMatrix, pos);
                
                    if (bitsSubOptimal <= threshold) {
                        StrandedFeature.Template template = new StrandedFeature.Template();
                        template.type = "motif" + w;
	                    template.source = "mica";
	                    template.location = new RangeLocation(pos, pos + motifs[w].columns() - 1);
	                    template.strand = strand == StrandedFeature.UNKNOWN ? StrandedFeature.NEGATIVE : strand;
	                    template.annotation = new SmallAnnotation();
	                    template.annotation.setProperty("motif.score", new Double(threshold - bitsSubOptimal));
	                    template.annotation.setProperty("motif.draw_hint", posKey);
	                    seq.createFeature(template);
	                }
                }
            }
        }
    }
    
    private boolean isRealDNA(Symbol s) {
        if (s == DNATools.a()) {
            return true;
        } else if (s == DNATools.c()) {
            return true;
        } else if (s == DNATools.g()) {
            return true;
        } else if (s == DNATools.t()) {
            return true;
        }
        return false;
    }
    
    private void scanGaps(SymbolList align, Sequence seq, int part, StrandedFeature.Strand strand)
    		throws Exception
    {
        int gapStart = -1;
        int ungapStart = -1;
        for (int pos = 1; pos <= align.length(); ++pos) {
            BasisSymbol bs = (BasisSymbol) align.symbolAt(pos);
            Symbol s = (Symbol) bs.getSymbols().get(part);
            boolean dna = isRealDNA(s);
            if (!dna) {
                if (gapStart < 0) {
                    gapStart = pos;
                }
                if (ungapStart > 0) {
                    markUngap(seq, ungapStart, pos - 1, strand);
                    ungapStart = -1;
                }
            } else {
                if (ungapStart < 0) {
                    ungapStart = pos;
                }
                if (gapStart > 0) {
                    markGap(seq, gapStart, pos - 1, strand);
                    gapStart = -1;
                }
            }
        }
        
        if (gapStart > 0) {
            markGap(seq, gapStart, seq.length(), strand);
        }
        if (ungapStart > 0) {
            markUngap(seq, ungapStart, seq.length(), strand);
        }
    }
    
    private boolean isMatch(BasisSymbol sym)
    		throws Exception
    {
        Symbol gap = alpha.getGapSymbol();
        for (Iterator si = sym.getSymbols().iterator(); si.hasNext(); ) {
            if (!isRealDNA((Symbol) si.next())) {
                return false;
            }
        }
        return true;
    }
    
    private void markGap(Sequence seq, int start, int end, StrandedFeature.Strand strand) 
    		throws Exception
    {
        StrandedFeature.Template temp = new StrandedFeature.Template();
        temp.type = "gap";
        temp.source = "bjneedle";
        temp.location = new RangeLocation(start, end);
        temp.annotation = Annotation.EMPTY_ANNOTATION;
        temp.strand= strand;
        // System.err.println("Creating gap from " + temp.location.getMin() + " to " + temp.location.getMax());
        seq.createFeature(temp);
    }
    
    private void markUngap(Sequence seq, int start, int end, StrandedFeature.Strand strand)
    		throws Exception
    {
        StrandedFeature.Template temp = new StrandedFeature.Template();
        temp.type = "ungap";
        temp.source = "bjneedle";
        temp.location = new RangeLocation(start, end);
        temp.annotation = Annotation.EMPTY_ANNOTATION;
        temp.strand= strand;
        // System.err.println("Creating ungap from " + temp.location.getMin() + " to " + temp.location.getMax());
        seq.createFeature(temp);
    }
    
    private SymbolList extractAlignPart(SymbolList align, int part)
		throws Exception
	{
        List<Symbol> l = new ArrayList<Symbol>();
        
		for (int pos = 1; pos <= align.length(); ++pos) {
		    BasisSymbol bs = (BasisSymbol) align.symbolAt(pos);
		    Symbol s = (Symbol) bs.getSymbols().get(part);
		    l.add(s);
		}
		
		return new SimpleSymbolList(alpha, l);
	}
    
    private double seqId(SymbolList align)
		throws Exception
    {
        int matches = 0;
        for (int pos = 1; pos <= align.length(); ++pos) {
            BasisSymbol bs = (BasisSymbol) align.symbolAt(pos);
		    Symbol s1 = (Symbol) bs.getSymbols().get(0);
		    Symbol s2 = (Symbol) bs.getSymbols().get(1);
		    if (s1 == s2) {
		        ++matches;
		    }
        }
        
        return (1.0 * matches) / align.length();
    }
    
    private void markId(SymbolList niceAlign, Sequence mView, int start, int end) 
    		throws Exception
    {
        int alignLength = end - start + 1;
        if (alignLength >= 10) {
            int fragSize;
            {
                int frags = ((int) Math.ceil((1.0 * alignLength) / 60));
                fragSize = ((int) Math.floor((1.0 * alignLength) / frags));
            }
            
            int fragStart = start;
            while (fragStart < end) {
                int fragEnd = Math.min(end, fragStart + fragSize - 1);
                // allow for rounding errors.
                if (end - fragEnd < 10) {
                    fragEnd = end;
                }
                
                StrandedFeature.Template temp = new StrandedFeature.Template();
	    	        temp.type = "id";
	    	        temp.source = "bjneedle";
	    	        temp.location = new RangeLocation(fragStart, fragEnd);
	    	        temp.strand = StrandedFeature.UNKNOWN;
	    	        temp.annotation = new SmallAnnotation();
	    	        temp.annotation.setProperty("id", new Double(seqId(niceAlign.subList(fragStart, fragEnd))));
	    	        mView.createFeature(temp);
	    	        
	    	        fragStart = fragEnd + 1;
            }
        }
    }
    
    private SymbolList processAlignment(SymbolList rawAlign)
    		throws Exception
    {
        Alphabet alignAlpha = rawAlign.getAlphabet();
        Symbol alignGap = alignAlpha.getGapSymbol();
        List<Symbol> sl = new ArrayList<Symbol>();
        for (Iterator<Symbol> i = rawAlign.iterator(); i.hasNext(); ) {
            Symbol s = i.next();
            if (! (s instanceof BasisSymbol)) {
                continue;
            }
            if (s == alignGap) {
                continue;
            }
            sl.add(s);
        }
        SymbolList niceAlign = new SimpleSymbolList(alignAlpha, sl);
        return niceAlign;
    }
    
    private SymbolList getAlignment(String seqName, int sp0, int sp1)
    		throws Exception
    	{
        File cacheFile = null;
        if (cache != null) {
            cacheFile = new File(cache, "alicache-" + seqName + "-" + sp0 + "-" + sp1 + ".jos");
        }
        
        if (cacheFile != null && cacheFile.exists()) {
            System.err.println("Phew, found requested alignment in the cache :-)");
            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(cacheFile));
            SymbolList ali = (SymbolList) ois.readObject();
            ois.close();
            return processAlignment(ali);
        }
        
        Sequence seq0 = seqs[sp0].getSequence(seqName);
        Sequence seq1 = seqs[sp1].getSequence(seqName);
        
        System.err.println("Running alignment");
        SymbolList rawAlign = aligner.align(seq0, seq1);
        System.err.println("Done alignment");
                
        if (cacheFile != null) {
            System.err.println("Writing alignment back to cache.  It'll be faster next time!");
            ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(cacheFile));
            oos.writeObject(rawAlign);
            oos.close();
        }
        
        return processAlignment(rawAlign);
    }
    
    private Sequence makeAlignExplorer(String seqName, int sp0, int sp1) 
    	    throws Exception 
    	{
        SymbolList niceAlign = getAlignment(seqName, sp0, sp1);
        
        Sequence mView = new SimpleSequence(
                new DummySymbolList(DNATools.getDNA(), niceAlign.length()),
                seqName,
                seqName,
                Annotation.EMPTY_ANNOTATION
        );
        
        scanGaps(niceAlign, mView, 0, StrandedFeature.POSITIVE);
        scanGaps(niceAlign, mView, 1, StrandedFeature.NEGATIVE);
        
        {
            int matchStart = -1;
	        for (int pos = 1; pos < niceAlign.length(); ++pos) {
	            boolean match = isMatch((BasisSymbol) niceAlign.symbolAt(pos));
	            if (match && matchStart < 0) {
	                matchStart = pos;
	            } else if (!match && matchStart > 0) {
	                markId(niceAlign, mView, matchStart, pos - 1);
	                matchStart = -1;
	            }
	        }
	        if (matchStart > 0) {
	            markId(niceAlign, mView, matchStart, niceAlign.length());
	        }
	    }
        
        annotateMotifs(extractAlignPart(niceAlign, 0), mView, StrandedFeature.UNKNOWN, "top");
        annotateMotifs(extractAlignPart(niceAlign, 1), mView, StrandedFeature.UNKNOWN, "bottom");
        
        return mView;
    }
    
    private static double scoreWM(SymbolList sl, WeightMatrix wm, int pos)
    	   throws Exception
    {
        /*
        
	    double score = 1.0;
	    for (int col = 0; col < wm.columns(); ++col) {
	        Symbol s = sl.symbolAt(pos + col);
	        if (s instanceof AtomicSymbol) {
	            score *= wm.getColumn(col).getWeight(sl.symbolAt(pos + col));
	        } else {
	            return 0.0;
	        }
	    }
	    return score;
	    
	    */
        
	    double score = 1.0;
	    Symbol gap = sl.getAlphabet().getGapSymbol();
	    try {
		    for (int col = 0, scol = 0; col < wm.columns(); ++col, ++scol) {
		        Symbol s = sl.symbolAt(pos + scol);
		        while (s == gap) {
		            s = sl.symbolAt(pos + ++scol);
		        }
		        if (s instanceof AtomicSymbol) {
		            score *= wm.getColumn(col).getWeight(s);
		        } else {
		            return 0.0;
		        }
		    }
	    } catch (IndexOutOfBoundsException ex) {
	        return 0;
	    }
	    return score;
	}

    
    private static double wmMaxScore(WeightMatrix wm)
    		throws Exception
    {
	    double wmScore = 1.0;
	    for (int c = 0; c < wm.columns(); ++c) {
	        Distribution col = wm.getColumn(c);
	        double colScore = 0;
	        for (Iterator i = ((FiniteAlphabet) col.getAlphabet()).iterator(); i.hasNext(); ) {
	            colScore = Math.max(colScore, col.getWeight((Symbol) i.next()));
	        }
	        wmScore *= colScore;
	    }
	    return wmScore;
	}
    
    private static double scoreBM(SymbolList sl, WeightMatrix wm, int pos)
    		throws Exception
    {
        double score = 0.0;
	    Symbol gap = sl.getAlphabet().getGapSymbol();
	    try {
		    for (int col = 0, scol = 0; col < wm.columns(); ++col, ++scol) {
		        Symbol s = sl.symbolAt(pos + scol);
		        while (s == gap) {
		            s = sl.symbolAt(pos + ++scol);
		        }
		        if (s instanceof AtomicSymbol) {
		            score += wm.getColumn(col).getWeight(s);
		        } else {
		            return Double.NEGATIVE_INFINITY;
		        }
		    }
	    } catch (IndexOutOfBoundsException ex) {
	        return Double.NEGATIVE_INFINITY;
	    }
	    return score;
    }
    
    
    private static double bmMaxScore(WeightMatrix wm)
    		throws Exception
    {
        double wmScore = 0.0;
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution col = wm.getColumn(c);
            double colScore = Double.NEGATIVE_INFINITY;
            for (Iterator i = ((FiniteAlphabet) col.getAlphabet()).iterator(); i.hasNext(); ) {
                colScore = Math.max(colScore, col.getWeight((Symbol) i.next()));
            }
            wmScore += colScore;
        }
        return wmScore;
    }
}
