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
import net.derkholm.nmica.seq.WmTools;
import net.derkholm.nmica.seq.motifxml.MotifHandler;
import net.derkholm.nmica.seq.motifxml.MotifWriter;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.ProteinTools.*;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.utils.stax.*;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;

import javax.swing.*;
import javax.xml.parsers.SAXParserFactory;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.SystemFlavorMap;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetAdapter;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
/**
 * @author Thomas Down
 */

@App(overview="View a set of motifs", generateStub=true)
public class MotifViewer {
    private boolean revComp;
    
    @Option(help="Should the program also show the reverse-complement of each motif", optional=true)
    public void setRevComp(boolean b) {
        this.revComp = b;
    }

    public void main(String[] args)
    		throws Exception
    {
        System.setProperty("apple.laf.useScreenMenuBar", "true");
        
        if (args.length > 0) {
            for (int i = 0; i < args.length; ++i) {
                MotifSet set = new MotifSet();
                set.loadFromFile(new File(args[i]));
                set.show();
            }
        } else {
            MotifSet blankSet = new MotifSet();
            blankSet.show();
        }
    }
    
    public static class Motif {
        private WeightMatrix weightMatrix;
        private String name;
        private double threshold;
        
        public Motif() {
            
        }
        
        public String getName() {
            return name;
        }
        public void setName(String name) {
            this.name = name;
        }
        public double getThreshold() {
            return threshold;
        }
        public void setThreshold(double threshold) {
            this.threshold = threshold;
        }
        public WeightMatrix getWeightMatrix() {
            return weightMatrix;
        }
        public void setWeightMatrix(WeightMatrix weightMatrix) {
            this.weightMatrix = weightMatrix;
        }
    }
    
    private static class XMotifHandler extends StAXContentHandlerBase {
        private Motif motif = new Motif();
        
        public Motif getMotif() {
            return motif;
        }
        
        public void startElement(String nsURI,
                String localName,
                String qName,
                Attributes attrs,
                DelegationManager dm
        )
                throws SAXException
        {
                    if ("weightmatrix".equals(localName)) {
                        dm.delegate(new MotifHandler());
                    } else if ("name".equals(localName)) {
                        dm.delegate(new StringElementHandlerBase() {
                            protected void setStringValue(String val) throws SAXException {
                                motif.setName(val);
                            }
                        });
                    }
        }
        
        public void endElement(String nsURI,
                String localName,
                String qName,
                StAXContentHandler delegate
        )
                throws SAXException
        {
                    if (delegate instanceof MotifHandler) {
                        motif.setWeightMatrix(((MotifHandler) delegate).getWeightMatrix());
                    }
       }
    }
    
    private static final DataFlavor WM_FLAVOR = new DataFlavor(DataFlavor.javaJVMLocalObjectMimeType, "WeightMatrix");
    
    private class MotifSet {
        private JFrame frame;
        
        private List<Motif> motifs = new ArrayList<Motif>();
        
        public MotifSet() {
            frame = new JFrame("Motif viewer");
            
            JMenuBar menubar = new JMenuBar();
            {
                JMenu menu = new JMenu("File");
                menu.add(new AbstractAction("New") {
                    public void actionPerformed(ActionEvent aev) {
                        if (motifs.size() != 0) {
                            MotifSet newMS = new MotifSet();
                            newMS.show();
                        }
                    }
                } );
                menu.add(new AbstractAction("Load") {
                    public void actionPerformed(ActionEvent aev) {
                        FileDialog filer = new FileDialog(frame, "Load motifs", FileDialog.LOAD);
	                      // fixme: should we be using filter.setVisible(true) instead?
                        filer.show();
                        File loadFile = new File(new File(filer.getDirectory()), filer.getFile());
                        
                        try {
	                        if (motifs.size() == 0) {
	                            loadFromFile(loadFile);
	                        } else {
	                            MotifSet newMS = new MotifSet();
	                            newMS.loadFromFile(loadFile);
	                            newMS.show();
	                        }
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                });
                menu.add(new AbstractAction("Save") {
                    public void actionPerformed(ActionEvent aev) {
                        FileDialog filer = new FileDialog(frame, "Save motifs", FileDialog.SAVE);
                        filer.show();
                        File saveFile = new File(new File(filer.getDirectory()), filer.getFile());
                        
                        try {
                            saveToFile(saveFile);
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                });
                menu.add(new AbstractAction("Export STUBB") {
                    public void actionPerformed(ActionEvent e) {
                        FileDialog filer = new FileDialog(frame, "Export motifs", FileDialog.SAVE);
                        filer.show();
                        File saveFile = new File(new File(filer.getDirectory()), filer.getFile());
                        
                        try {
                            stubbToFile(saveFile);
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                });
                menu.add(new AbstractAction("Export JOS") {
                    public void actionPerformed(ActionEvent e) {
                        FileDialog filer = new FileDialog(frame, "Export motifs", FileDialog.SAVE);
                        filer.show();
                        File saveFile = new File(new File(filer.getDirectory()), filer.getFile());
                        
                        try {
                            josToFile(saveFile);
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                });
                menu.add(new AbstractAction("Quit") {
                    public void actionPerformed(ActionEvent aev) {
                        System.exit(0);
                    }
                });
                
                menubar.add(menu);
            }
            
            frame.setJMenuBar(menubar);
            
            ((SystemFlavorMap) SystemFlavorMap.getDefaultFlavorMap()).addFlavorForUnencodedNative("WM", WM_FLAVOR);
            ((SystemFlavorMap) SystemFlavorMap.getDefaultFlavorMap()).addUnencodedNativeForFlavor(WM_FLAVOR, "WM");
            
            DropTarget target = new DropTarget(
                    frame,
                    new DropTargetAdapter() {
                        public void drop(DropTargetDropEvent dtde) {
                            System.err.println("User actions: " + dtde.getDropAction());
                            dtde.acceptDrop(dtde.getDropAction());
                            Transferable t = dtde.getTransferable();
                            try {
                                WeightMatrix wm = (WeightMatrix) t.getTransferData(WM_FLAVOR);
                                Motif m = new Motif();
                                m.setWeightMatrix(wm);
                                motifs.add(m);
                                dtde.dropComplete(true);
                                buildViewer();
                                frame.pack();
                                
                            } catch (UnsupportedFlavorException e) {
                                // TODO Auto-generated catch block
                                e.printStackTrace();
                            } catch (IOException e) {
                                // TODO Auto-generated catch block
                                e.printStackTrace();
                            }
                        }
                    }
            );
        }
        
        
        
        public void loadFromFile(File f)
        		throws Exception
        {
            List<Object> l = new ArrayList<Object>();
            try {
                ObjectInputStream ois = new ObjectInputStream(new FileInputStream(f));
                while (true) {
                    l.add(ois.readObject());
                }
            } catch (Exception ex) {
            }
            if (l.size() > 0) {
                motifs = new ArrayList<Motif>();
                for (int i = 0; i < l.size(); ++i) {
                    Motif m = new Motif();
                    m.setWeightMatrix((WeightMatrix) l.get(i));
                    m.setName("motif" + i);
                    motifs.add(m);
                }
            } else {
                Reader r = new FileReader(f);
                SAXParserFactory spf = SAXParserFactory.newInstance();
        	       spf.setNamespaceAware(true);
        	       XMLReader parser = spf.newSAXParser().getXMLReader();
        	       
        	       final List<Motif> motifList = new ArrayList<Motif>();
        	       StAXContentHandler handler = new StAXContentHandlerBase() {
        	           public void startElement(String nsURI,
        	                   String localName,
        	                   String qName,
        	                   Attributes attrs,
        	                   DelegationManager dm
        	          )
        	                   throws SAXException
        	          {
        	               if ("motif".equals(localName)) {
        	                           dm.delegate(new XMotifHandler());
        	               }
        	          }
        	           
        	          public void endElement(String nsURI,
        	                   String localName,
        	                   String qName,
        	                   StAXContentHandler delegate
        	          )
        	                   throws SAXException
        	          {
        	                       if (delegate instanceof XMotifHandler) {
        	                           motifList.add(((XMotifHandler) delegate).getMotif());
        	                       }
        	          }
        	       };
        	       parser.setContentHandler(new SAX2StAXAdaptor(handler));
        	       parser.parse(new InputSource(r));
        	       
        	       motifs = motifList;
            }
            
            frame.setTitle("Motif viewer: " + f.getName());
            buildViewer();
        }
        
        public void saveToFile(File f)
        		throws Exception
        {
            MotifWriter mw = new MotifWriter();
            
            PrintWriter pw = new PrintWriter(new FileWriter(f));
            XMLWriter xw = new PrettyXMLWriter(pw);
            xw.openTag("motifset");
            for (int i = 0; i < motifs.size(); ++i) {
                Motif m = motifs.get(i);
                
                xw.openTag("motif");
                  xw.openTag("name");
                  xw.print(m.getName());
                  xw.closeTag("name");
                  mw.writeMatrix(m.getWeightMatrix(), xw);
                xw.closeTag("motif");
            }
            xw.closeTag("motifset");
            pw.close();
            
            frame.setTitle("Motif viewer: " + f.getName());
        }
        
        public void stubbToFile(File f)
			throws Exception
		{
		    PrintWriter pw = new PrintWriter(new FileWriter(f));
		    for (int i = 0; i < motifs.size(); ++i) {
		        Motif m = motifs.get(i);
		        WeightMatrix wm = m.getWeightMatrix();
		        pw.println(">motif" + i + " " + wm.columns());
		        for (int c = 0; c < wm.columns(); ++c) {
			   //for (Iterator p=ProteinTools.getAlphabet().iterator();p.hasNext();) {
			   //	pw.print ((int) Math.round ( wm.getColumn(c).getWeight(p.next()) *100));
		           //     pw.print('\t');
			   //}
					
		        }
		        pw.println("<");
		    }
		    pw.close();
		}
        
        public void josToFile(File f)
        	   throws Exception
        {
            ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(f));
            for (int i = 0; i < motifs.size(); ++i) {
                Motif m = motifs.get(i);
                oos.writeObject(m.getWeightMatrix());
            }
            oos.close();
        }
        
        private void buildViewer() {
            Container content = frame.getContentPane();
            content.removeAll();
            
            Box b = Box.createVerticalBox();
            for (int i = 0; i < motifs.size(); ++i) {
                final Motif m = motifs.get(i);
                
                Box h = Box.createHorizontalBox();
                final JTextField label = new JTextField(m.getName());
                label.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aev) {
                        m.setName(label.getText());
                    }
                } );
                h.add(label);
                h.add(Box.createHorizontalStrut(50));
                WMPanel panel = new WMPanel(m.getWeightMatrix());
                h.add(panel);
                
                /*
                final DragSource dragSource = new DragSource();
                DragGestureListener dgl = new DragGestureListener() {

                    public void dragGestureRecognized(DragGestureEvent dge) {
                        dragSource.startDrag(
                                dge,
                                DragSource.DefaultCopyDrop,
                                new Transferable() {
                                    public DataFlavor[] getTransferDataFlavors() {
                                        return new DataFlavor[] {WM_FLAVOR};
                                    }

                                    public boolean isDataFlavorSupported(DataFlavor flavor) {
                                        return flavor == WM_FLAVOR;
                                    }

                                    public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException, IOException {
                                        if (flavor == WM_FLAVOR) {
                                            return m.getWeightMatrix();
                                        } else {
                                            throw new UnsupportedFlavorException(flavor);
                                        }
                                    }
                                    
                                },
                                new DragSourceAdapter() {
                                    public void dragDropEnd(DragSourceDropEvent ev) {
                                        System.err.println("Drag complete, success=" + ev.getDropSuccess());
                                    }
                                }
                        );
                    }
                    
                };
                dragSource.createDefaultDragGestureRecognizer(panel, DnDConstants.ACTION_COPY_OR_MOVE, dgl);
                
                */
                
                if (revComp) {
                     h.add(Box.createHorizontalStrut(50));
                     try {
                         h.add(new WMPanel(WmTools.reverseComplement(m.getWeightMatrix())));
                     } catch (IllegalAlphabetException iae) {
                     }
                }
                b.add(h);
            }
            content.add(b);
        }
        
        public void show() {
            frame.pack();
            frame.setVisible(true);
        }
    }
}
