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

package net.derkholm.nmica.motif;

import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.SAXParserFactory;

import net.derkholm.nmica.seq.motifxml.MotifWriter;

import org.biojava.bio.Annotation;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.SAX2StAXAdaptor;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;

/**
 * Utilities for reading and writing motifs.
 * 
 * @author thomas
 */
public class MotifIOTools {
	
	public static final String XMS_NS = "http://biotiffin.org/XMS/";
	
    private MotifIOTools() {
    }
    
    public static Motif[] loadWmJOS(InputStream is)
    		throws Exception
    {        
        List<WeightMatrix> l = new ArrayList<WeightMatrix>();
        try {
            ObjectInputStream ois = new ObjectInputStream(is);
            while (true) {
                l.add((WeightMatrix) ois.readObject());
            }
        } catch (Exception ex) {
        }
            List<Motif> motifs = new ArrayList<Motif>();
            for (int i = 0; i < l.size(); ++i) {
                Motif m = new Motif();
                m.setWeightMatrix((WeightMatrix) l.get(i));
                m.setName("motif" + i);
                motifs.add(m);
            }
        return motifs.toArray(new Motif[0]);
    }
    
    public static Motif[] loadMotifSetXML(InputStream is)
        throws Exception
    {
        return loadMotifSetXML(new InputSource(is));
    }
    
    public static Motif[] loadMotifSetXML(Reader is)
        throws Exception
    {
        return loadMotifSetXML(new InputSource(is));
    }
    
    private static Motif[] loadMotifSetXML(InputSource is)
    		throws Exception
    {
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
        parser.parse(is);
        
        return motifList.toArray(new Motif[0]);
    }
    
    public static void writeMotifSetXML(OutputStream os, Motif[] motifs)
		throws Exception
	{
    	writeMotifSetXML(os, motifs, Annotation.EMPTY_ANNOTATION);
	}
    
    public static void writeMotifSetXML(OutputStream os, Motif[] motifs, Annotation props)
    		throws Exception
    {
        MotifWriter mw = new MotifWriter();
        
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(os));
        XMLWriter xw = new PrettyXMLWriter(pw);
        xw.openTag("motifset");
        xw.attribute("xmlns", XMS_NS);
        writeAnnotation(xw, props);
        for (int i = 0; i < motifs.length; ++i) {
            Motif m = motifs[i];
            
            xw.openTag("motif");
              xw.openTag("name");
              xw.print(m.getName());
              xw.closeTag("name");
              if (m.getVersion() >= 0) {
            	  xw.openTag("version");
            	  xw.print("" + m.getVersion());
            	  xw.closeTag("version");
              }
              mw.writeMatrix(m.getWeightMatrix(), xw);
              double score = m.getThreshold();
              if (!Double.isNaN(score)) {
                  xw.openTag("threshold");
                  xw.print("" + score);
                  xw.closeTag("threshold");
              }
              writeAnnotation(xw, m.getAnnotation());
            xw.closeTag("motif");
        }
        xw.closeTag("motifset");
        pw.flush();
    }
    
    private static void writeAnnotation(XMLWriter xw, Annotation a)
    	throws Exception
    {
    	for (Iterator<Map.Entry<?,?>> pi = a.asMap().entrySet().iterator(); pi.hasNext(); ) {
            Map.Entry<?,?> prop = pi.next();
            xw.openTag("prop");
            xw.openTag("key");
            xw.print(prop.getKey().toString());
            xw.closeTag("key");
            xw.openTag("value");
            xw.print(prop.getValue().toString());
            xw.closeTag("value");
            xw.closeTag("prop");
        }
    }
}
