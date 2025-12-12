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

package net.derkholm.nmica.model.motif;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.XMLStreamWriter;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.dist.OrderNDistributionFactory;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * Utilities for reading and writing background models in XML format
 * 
 * @author thomas
 */

public class MosaicIO {
	public static final String MOSAIC_NS = "http://biotiffin.org/NestedMICA/mosaic-background/";
	
	private MosaicIO() {
	}
	
	public static void writeMosaic(XMLStreamWriter xw, Mosaic mosaic)
		throws XMLStreamException
	{
		FiniteAlphabet alpha = (FiniteAlphabet) mosaic.getAlphabet();
		int order = 1;
		{
			List al = alpha.getAlphabets();
			if (al.size() > 1) {
				order = al.size() + ((Alphabet) al.get(0)).getAlphabets().size() - 1;
				alpha = (FiniteAlphabet) al.get(al.size() - 1);
			}
		}
		
		xw.writeStartDocument("UTF-8", "1.0");
		xw.setDefaultNamespace(MOSAIC_NS);
		xw.writeStartElement(MOSAIC_NS, "mosaicBackground");
		xw.writeDefaultNamespace(MOSAIC_NS);
		  xw.writeStartElement(MOSAIC_NS, "alphabet");
		  xw.writeCharacters(alpha.getName());
		  xw.writeEndElement();
		  xw.writeStartElement(MOSAIC_NS, "order");
		  xw.writeCharacters("" + (order - 1));
		  xw.writeEndElement();
		  xw.writeStartElement(MOSAIC_NS, "transition");
		  xw.writeCharacters("" + mosaic.getTransition());
		  xw.writeEndElement();
		  for (Object pk : mosaic.getAnnotation().keys()) {
			  Object pv = mosaic.getAnnotation().getProperty(pk);
			  xw.writeStartElement(MOSAIC_NS, "prop");
				xw.writeStartElement(MOSAIC_NS, "key");
				xw.writeCharacters(pk.toString());
				xw.writeEndElement();
				xw.writeStartElement(MOSAIC_NS, "value");
				xw.writeCharacters(pv.toString());
				xw.writeEndElement();
			  xw.writeEndElement();
		  }
		  for (Distribution mc : mosaic.getDistributions()) {
			  xw.writeStartElement(MOSAIC_NS, "class");
			  if (order > 1) {
				  OrderNDistribution omc = (OrderNDistribution) mc;
				  for (Iterator<?> si = ((FiniteAlphabet) omc.getConditioningAlphabet()).iterator(); si.hasNext(); ) {
					  Symbol s = (Symbol) si.next();
					  xw.writeStartElement(MOSAIC_NS, "conditionedDistribution");
					  xw.writeAttribute("prefix", s.getName());
					  try {
						writeDistribution(xw, omc.getDistribution(s));
					  } catch (IllegalSymbolException e) {
						  throw new XMLStreamException(e);
					  }
					  xw.writeEndElement();
				  }
			  } else {
				  writeDistribution(xw, mc);
			  }
			  xw.writeEndElement();
		  }
		xw.writeEndElement();
		xw.writeEndDocument();
	}
	
	private static void writeDistribution(XMLStreamWriter xw, Distribution mc)
		throws XMLStreamException
	{
		for (Iterator<?> si = ((FiniteAlphabet) mc.getAlphabet()).iterator(); si.hasNext(); ) {
			  Symbol s = (Symbol) si.next();
			  xw.writeStartElement(MOSAIC_NS, "weight");
			  xw.writeAttribute("symbol", s.getName());
			  try {
				  xw.writeCharacters("" + mc.getWeight(s));
			  } catch (Exception ex) {
				  throw new XMLStreamException(ex);
			  }
			  xw.writeEndElement();
		  }
	}
	
	private static String readTextElement(XMLStreamReader xr) 
		throws XMLStreamException
	{
		StringBuilder sb = new StringBuilder();
		while (true) {
			final int evt = xr.next();
			if (evt == XMLStreamConstants.CHARACTERS || evt == XMLStreamConstants.CDATA || evt == XMLStreamConstants.SPACE) {
				sb.append(xr.getText());
			} else if (evt == XMLStreamConstants.START_ELEMENT) {
				throw new XMLStreamException("Child element found in unexpected location");
			} else if (evt == XMLStreamConstants.END_ELEMENT) {
				return sb.toString();
			}
		} 
	}
	
	private static void readIntoDistribution(XMLStreamReader xr, Alphabet alpha, Distribution dist)
		throws XMLStreamException, BioException, ChangeVetoException
	{
		while (true) {
			final int cevt = xr.nextTag();
			if (cevt == XMLStreamConstants.START_ELEMENT) {
				if ("weight".equals(xr.getLocalName())) {
					Symbol sym = alpha.getTokenization("name").parseToken(xr.getAttributeValue(0));
					dist.setWeight(sym, Double.parseDouble(readTextElement(xr)));
				} else if ("conditionedDistribution".equals(xr.getLocalName())) {
					Symbol prefix = ((OrderNDistribution) dist).getConditioningAlphabet().getTokenization("name").parseToken(xr.getAttributeValue(0));
					readIntoDistribution(xr, alpha, ((OrderNDistribution) dist).getDistribution(prefix));
				}
			} else if (cevt == XMLStreamConstants.END_ELEMENT) {
				return;
			}
		}
	}
	
	public static Mosaic readMosaic(XMLStreamReader xr)
		throws XMLStreamException
	{
		xr.nextTag();
		if (!MOSAIC_NS.equals(xr.getNamespaceURI()) || !"mosaicBackground".equals(xr.getLocalName())) {
			throw new XMLStreamException("Not a mosaicBackground");
		}
		
		double transition = Double.NaN;
		FiniteAlphabet alpha  = null;
		List<Distribution> dists = new ArrayList<Distribution>();
		Annotation props = new SmallAnnotation();
		int order = 0;
		
		while (true) {
			int evt = xr.nextTag();
			if (evt == XMLStreamConstants.START_ELEMENT) {
				if ("alphabet".equals(xr.getLocalName())) {
					alpha = (FiniteAlphabet) AlphabetManager.alphabetForName(readTextElement(xr));
				} else if ("transition".equals(xr.getLocalName())) {
					transition = Double.parseDouble(readTextElement(xr));
				} else if ("order".equals(xr.getLocalName())) {
					order = Integer.parseInt(readTextElement(xr));
				} else if ("prop".equals(xr.getLocalName())) {
					xr.nextTag();
					if (! "key".equals(xr.getLocalName())) {
						throw new XMLStreamException("<prop> missing required <key>");
					}
					String key = readTextElement(xr);
					xr.nextTag();
					if (! "value".equals(xr.getLocalName())) {
						throw new XMLStreamException("<prop> missing required <value>");
					}
					String value = readTextElement(xr);
					if (xr.nextTag() != XMLStreamConstants.END_ELEMENT) {
						throw new XMLStreamException("Unexpected tag in <prop> element");
					}
					try {
						props.setProperty(key, value);
					} catch (IllegalArgumentException e) {
					} catch (ChangeVetoException e) {
					}
				} else if ("class".equals(xr.getLocalName())) {
					try {
						Alphabet a = alpha;
						Distribution dist;
						if (order > 0) {
							a = AlphabetManager.getCrossProductAlphabet(Collections.nCopies(order + 1, alpha));
							dist = OrderNDistributionFactory.DEFAULT.createDistribution(a);
						} else {
							dist = DistributionFactory.DEFAULT.createDistribution(a);
						}
						// System.err.println(dist.getAlphabet().getName());
						readIntoDistribution(xr, alpha, dist);
						dists.add(dist);
					} catch (BioException e) {
						// e.printStackTrace();
						throw new XMLStreamException(e);
					} catch (ChangeVetoException e) {
					}
				}
			} else if (evt == XMLStreamConstants.END_ELEMENT) {
				break;
			}
		}
		
		return new Mosaic(dists.toArray(new Distribution[0]), transition, props);
	}
}
