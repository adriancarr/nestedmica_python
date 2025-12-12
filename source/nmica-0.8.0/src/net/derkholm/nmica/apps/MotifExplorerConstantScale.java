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
import net.derkholm.nmica.seq.NMSimpleDistribution;
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.seq.align.Aligner;
import net.derkholm.nmica.utils.ArrayTools;
import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.utils.CollectTools;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.*;

import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.MutableTreeNode;
import javax.swing.tree.TreePath;
import java.awt.*;
import java.awt.event.*;
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
public class MotifExplorerConstantScale {
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
    private static final FiniteAlphabet DNA_ALPHA = DNATools.getDNA();
    
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
            Color.yellow
    } ;
    
    /**
     * @author thomas
     */
    private class MotifBeadPanel extends JPanel {
        Sequence seq = null;
	int isize = 500;
	int ilen = 0;
	double rmag = 4.0;
        
        public MotifBeadPanel() {
            super();
            setPreferredSize(new Dimension(isize, 220));
        }
        
        public void setSequence(Sequence seq) {
            this.seq = seq;
	    redraw();
        }

	public void setMag(double mag) {
	    this.rmag = mag;
	    redraw();
	}
        
        public void redraw() {
	    this.ilen=(int) (((float) seq.length())/rmag);
	    int isize=ilen + 40;
	    setPreferredSize(new Dimension(isize, 220));
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
            
            if (seq != null) {
		float rlen = (float) ilen;
                g.fillRect(20, 40, ilen, 10);
	            for (Iterator fi = seq.features(); fi.hasNext(); ) {
	                StrandedFeature sf = (StrandedFeature) fi.next();
	                
	                Matcher m = MOTIF_PATTERN.matcher(sf.getType());
	                if (m.matches()) {
	                    int w = Integer.parseInt(m.group(1));
	                    
		                int mid = (sf.getLocation().getMax() + sf.getLocation().getMin()) / 2;
		                int pos = 20 + (int) ((rlen * mid) / seq.length());
		                
		                Color c = Color.gray;
		                if (w < palette.length) {
		                    c = palette[w];
		                }
		                g.setColor(c);
		                
		                double score = ((Double) sf.getAnnotation().getProperty("motif.score")).doubleValue();
		                int height = (int) Math.floor(score * 10);
		                if (sf.getStrand() == StrandedFeature.POSITIVE) {
		                    g.fillRect(pos - 3, 40 - height, 6, height);
		                } else {
		                    g.fillRect(pos - 3, 50, 6, height);
		                }
	                } else if ("id".equals(sf.getType())) {
	                    int start = 20 + (int) ((rlen * sf.getLocation().getMin()) / seq.length());
	                    int end = 20 + (int) Math.ceil((rlen * sf.getLocation().getMax()) / seq.length());
	                    g.setColor(idColor(((Double) sf.getAnnotation().getProperty("id")).doubleValue()));
	                    g.fillRect(start, 40, Math.max(1, end-start), 10);
	                } else if ("gap".equals(sf.getType())) {
	                    Color c = Color.black;
	                    int start = 20 + (int) ((rlen * sf.getLocation().getMin()) / seq.length());
	                    int end = 20 + (int) ((rlen * sf.getLocation().getMax()) / seq.length());
	                    int height = 20;
	                    g.setColor(c);
		                if (sf.getStrand() == StrandedFeature.POSITIVE) {
		                    g.fillRect(start, 40 - height, Math.max(1, end - start), height);
		                } else {
		                    g.fillRect(start, 50, Math.max(1, end - start), height);
		                }
	                }
	            }

	            // this bit draws length arrow
	            g.setColor(Color.black);
		    int icen1 = (int) (((float) (ilen - 100))/2.0 + 20.0);
		    int icen2 = icen1 + 100;
			
	            g.drawLine(20, 130, icen1, 130);
	            g.drawLine(20, 130, 40, 120);
	            g.drawLine(20, 130, 40, 140);
	            
		    int iedge = ilen + 20;
	            g.drawLine(icen2, 130, iedge, 130);
	            g.drawLine(iedge, 130, ilen, 120);
	            g.drawLine(iedge, 130, ilen, 140);
	            
		    int icen3 = icen1 + 20;
	            g.drawString("" + seq.length() + "bp", icen3, 135);

	            if (seq.filter(new FeatureFilter.ByType("id")).countFeatures() > 0) {
		            g.drawString("%id", 5, 185);
		            int[] keyIds = new int[] {25, 40, 55, 70, 85, 100};
		            for (int k = 0; k < keyIds.length; ++k) {
		                double id = (1.0 * keyIds[k]) / 100;
		                Color c = idColor(id);
		                g.setColor(c);
		                g.fillRect(40 + k * 40, 160, 40, 40);
		                g.setColor(Color.black);
		                g.drawRect(40 + k * 40, 160, 40, 40);
		                g.drawString("" + keyIds[k] + "%", 45 + k * 40, 185);
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
            Color c = Color.gray;
            if (index < palette.length) {
                c = palette[index];
            }
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
    
    private Aligner aligner;
    
    public void setAlignerClass(String s) {
        this.alignerClass = s;
    }
    
    public void setScoreThreshold(double d) {
        this.scoreThreshold = d;
    }
    
    public void setSpecies(String [] species) {
        this.species = species;
    }
    
    public void setSeqs(Reader[] r) throws Exception {
        seqs = new SequenceDB[r.length];
        for (int s = 0; s < r.length; ++s) {
	        seqs[s] = new HashSequenceDB();
	        SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(r[s]));
	        while (si.hasNext()) {
	            seqs[s].addSequence(si.nextSequence());
	        }
        }
    }
    
    public void setMotifs(InputStream is)
    		throws Exception
    {
        ObjectInputStream ois = new ObjectInputStream(is);
        List l = new ArrayList();
        try {
            while (true) {
                l.add(ois.readObject());
            }
        } catch (Exception ex) {
        }
        motifs = (WeightMatrix[]) l.toArray(new WeightMatrix[0]);
    }
    
    public static void main(String[] args) throws Exception {
        MotifExplorerConstantScale app = new MotifExplorerConstantScale();
        CliTools.configureBean(app, args);
        app.run();
    }

    /**
     * 
     */
    private void run() 
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
        
        aligner = getClass()
                .getClassLoader()
                .loadClass(alignerClass)
                .asSubclass(Aligner.class)
                .newInstance();

	final JFrame explorerFrame = new JFrame("MotifExplorer");
        
        final Box mainArea = Box.createVerticalBox();
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

        mainArea.add(new JScrollPane(seqPanel));
        
	// add a magnification control menu
	JMenu magm = new JMenu("Magnification");
	final ButtonGroup group = new ButtonGroup();

	class TheEvent implements ActionListener
	{

	    TheEvent(){
		
	    }

	    public void actionPerformed(ActionEvent e){
		String tmag = group.getSelection().getActionCommand();
		double rmag = Double.parseDouble(tmag);
		seqPanel.setMag(rmag);
		explorerFrame.pack();
		//System.out.println(tmag);
		//JOptionPane.showMessageDialog(explorerFrame, tmag);
	    }
	}
    
	TheEvent mevent = new TheEvent();

	JRadioButtonMenuItem rbmi;
	
	rbmi = new JRadioButtonMenuItem("Variable",true);
	rbmi.setActionCommand("0");
	group.add(rbmi);
	rbmi.addActionListener(mevent);
	magm.add(rbmi);
	
	rbmi = new JRadioButtonMenuItem("2bp/pixel");
	rbmi.setActionCommand("2");
	group.add(rbmi);
	rbmi.addActionListener(mevent);
	magm.add(rbmi);
	
	rbmi = new JRadioButtonMenuItem("5bp/pixel");
	rbmi.setActionCommand("5");
	group.add(rbmi);
	rbmi.addActionListener(mevent);
	magm.add(rbmi);
	
	rbmi = new JRadioButtonMenuItem("10bp/pixel");
	rbmi.setActionCommand("10");
	group.add(rbmi);
	rbmi.addActionListener(mevent);
	magm.add(rbmi);
	
	JMenuBar menubar = new JMenuBar();
	menubar.add(magm);
	explorerFrame.setJMenuBar(menubar);

        final JTree seqList;
        {
            Set allSeqIds = new TreeSet();
            for (int s = 0; s < seqs.length; ++s) {
                allSeqIds.addAll(seqs[s].ids());
            }
            
            DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
            for (Iterator idi = allSeqIds.iterator(); idi.hasNext(); ) {
                String id = (String) idi.next();
                
                DefaultMutableTreeNode idNode = new DefaultMutableTreeNode(id);
                List<Integer> alignIdList = new ArrayList<Integer>();
                for (int s = 0; s < seqs.length; ++s) {
                    if (seqs[s].ids().contains(id)) {
                        MutableTreeNode seqNode = new DefaultMutableTreeNode(species[s]);
                        idNode.add(seqNode);
                        alignIdList.add(new Integer(s));
                    }
                }
                int[] alignIds = CollectTools.toIntArray(alignIdList);
                for (int i = 0; i < alignIds.length; ++i) {
                    for (int j = i + 1; j < alignIds.length; ++j) {
                        MutableTreeNode alignNode = new DefaultMutableTreeNode(new AlignMemento(i, j));
                        idNode.add(alignNode);
                    }
                }
                rootNode.add(idNode);
            }
            
            seqList = new JTree(rootNode);
            seqList.setRootVisible(false);
            
            seqList.addMouseListener(new MouseAdapter() {
                public void mouseClicked(MouseEvent mev) {
                    if (mev.getClickCount() == 2) {
                        TreePath tpath = seqList.getPathForLocation(mev.getPoint().x, mev.getPoint().y);
                        Object[] path = tpath.getPath();
                        if (path.length == 3) {
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
        annotateMotifs(mView, mView, StrandedFeature.UNKNOWN);
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
    
    private void annotateMotifs(SymbolList sl, Sequence seq, StrandedFeature.Strand strand)
    		throws Exception
    {
        for (int w = 0; w < motifs.length; ++w) {
            WeightMatrix fwdMatrix = wmToBm(motifs[w]);
            WeightMatrix backMatrix = wmToBm(WmTools.reverseComplement(motifs[w]));
            int maxPos = sl.length() - motifs[w].columns() + 1;
            double maxScore = bmMaxScore(fwdMatrix);
            for (int pos = 1; pos <= maxPos; ++pos) {
                {
                    // double score = scoreWM(sl, motifs[w], pos);
                    // double bitsSubOptimal = -Math.log(score / wmMaxScore(motifs[w])) / LOG_2;
                    double bitsSubOptimal = maxScore - scoreBM(sl, fwdMatrix, pos);
                    
                    if (bitsSubOptimal <= scoreThreshold) {
                        StrandedFeature.Template template = new StrandedFeature.Template();
                        template.type = "motif" + w;
	                    template.source = "mica";
	                    template.location = new RangeLocation(pos, pos + motifs[w].columns() - 1);
	                    template.strand = strand == StrandedFeature.UNKNOWN ? StrandedFeature.POSITIVE : strand;
	                    template.annotation = new SmallAnnotation();
	                    template.annotation.setProperty("motif.score", new Double(scoreThreshold - bitsSubOptimal));
	                    seq.createFeature(template);
	                }
                }
                {
                    // double score = scoreWM(sl, backMatrix, pos);
                    // double bitsSubOptimal = -Math.log(score / wmMaxScore(motifs[w])) / LOG_2;
                    double bitsSubOptimal = maxScore - scoreBM(sl, backMatrix, pos);
                
                    if (bitsSubOptimal <= scoreThreshold) {
                        StrandedFeature.Template template = new StrandedFeature.Template();
                        template.type = "motif" + w;
	                    template.source = "mica";
	                    template.location = new RangeLocation(pos, pos + motifs[w].columns() - 1);
	                    template.strand = strand == StrandedFeature.UNKNOWN ? StrandedFeature.NEGATIVE : strand;
	                    template.annotation = new SmallAnnotation();
	                    template.annotation.setProperty("motif.score", new Double(scoreThreshold - bitsSubOptimal));
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
        for (int pos = 1; pos <= align.length(); ++pos) {
            BasisSymbol bs = (BasisSymbol) align.symbolAt(pos);
            Symbol s = (Symbol) bs.getSymbols().get(part);
            boolean dna = isRealDNA(s);
            if (!dna && gapStart < 0) {
                gapStart = pos;
            } else if (dna && gapStart > 0) {
                markGap(seq, gapStart, pos - 1, strand);
                gapStart = -1;
            }
        }
        
        if (gapStart > 0) {
            markGap(seq, gapStart, seq.length(), strand);
        }
    }
    
    private boolean isMatch(BasisSymbol sym)
    		throws Exception
    {
        Symbol gap = DNA_ALPHA.getGapSymbol();
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
    
    private SymbolList extractAlignPart(SymbolList align, int part)
		throws Exception
	{
        List<Symbol> l = new ArrayList<Symbol>();
        
		for (int pos = 1; pos <= align.length(); ++pos) {
		    BasisSymbol bs = (BasisSymbol) align.symbolAt(pos);
		    Symbol s = (Symbol) bs.getSymbols().get(part);
		    l.add(s);
		}
		
		return new SimpleSymbolList(DNA_ALPHA, l);
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
        for (Iterator i = rawAlign.iterator(); i.hasNext(); ) {
            Symbol s = (Symbol) i.next();
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
        
        annotateMotifs(extractAlignPart(niceAlign, 0), mView, StrandedFeature.POSITIVE);
        annotateMotifs(extractAlignPart(niceAlign, 1), mView, StrandedFeature.NEGATIVE);
        
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
