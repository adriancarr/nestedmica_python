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

import net.derkholm.nmica.seq.motifxml.MotifHandler;

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.DoubleElementHandlerBase;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * @author thomas
 */
 class XMotifHandler extends StAXContentHandlerBase {
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
                    } else if ("threshold".equals(localName)) {
                    	dm.delegate(new DoubleElementHandlerBase() {
							protected void setDoubleValue(double val) throws SAXException {
								motif.setThreshold(val);
							}
                    	});
                    } else if ("prop".equals(localName)) {
                        dm.delegate(new StAXContentHandlerBase() {
                            private String key;
                            private String value;
                            
                            public void startElement(String nsURI,
                                    String localName,
                                    String qName,
                                    Attributes attrs,
                                    DelegationManager dm
                            )
                                    throws SAXException
                            {
                                if ("key".equals(localName)) {
                                    dm.delegate(new StringElementHandlerBase() {
                                        protected void setStringValue(String val) {
                                            key = val;
                                        }
                                    });
                                } else if ("value".equals(localName)) {
                                    dm.delegate(new StringElementHandlerBase() {
                                        protected void setStringValue(String val) {
                                            value = val;
                                        }
                                    });
                                }
                            }
                            
                            public void endTree() 
                            		throws SAXException
                            {
                                try {
                                    motif.getAnnotation().setProperty(key, value);
                                } catch (ChangeVetoException ex) {
                                    throw new SAXException(ex);
                                }
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
