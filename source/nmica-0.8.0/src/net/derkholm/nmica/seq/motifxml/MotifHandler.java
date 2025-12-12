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

package net.derkholm.nmica.seq.motifxml;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.DoubleElementHandlerBase;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Read in weight matrices which have been stored in XML format.
 * 
 * @author Thomas Down
 */
public class MotifHandler extends StAXContentHandlerBase {
    private WeightMatrix wm;
    private int depth = 0;
            
    public void startElement(String nsURI,
        String localName,
        String qName,
        Attributes attrs,
        DelegationManager dm
    )
        throws SAXException
    {
        if (depth == 0) {
            int columns = Integer.parseInt(attrs.getValue("columns"));
            String alphaName = attrs.getValue("alphabet");
            Alphabet alpha;
            try {
                alpha = AlphabetManager.alphabetForName(alphaName);
                wm = new SimpleWeightMatrix(alpha, columns, DistributionFactory.DEFAULT);
            } catch (Exception ex) {
                throw new SAXException("Couldn't find an alphabet named " + alphaName);
            }
        } else if (depth == 1) {
            if ("column".equals(localName)) {
                int pos = Integer.parseInt(attrs.getValue("pos"));
                dm.delegate(new DistributionHandler(wm.getColumn(pos)));
            }
        }
        ++depth;
    }

    public void endElement(String nsURI,
                           String localName,
			               String qName,
			               StAXContentHandler delegate)
                 throws SAXException
                 {
                     --depth;
                 }
                 
     public WeightMatrix getWeightMatrix() {
         return wm;
     }
     
     private static class DistributionHandler extends StAXContentHandlerBase {
         private Distribution dist;
         
         public DistributionHandler(Distribution dist) {
             super();
             this.dist = dist;
         }
         
         public void startElement(String nsURI,
            String localName,
            String qName,
            Attributes attrs,
            DelegationManager dm
         )
            throws SAXException
         {
             if ("weight".equals(localName)) {
                 String symName = attrs.getValue("symbol");
                 final Symbol sym;
                 try {
                     SymbolTokenization toke = dist.getAlphabet().getTokenization("name");
                     sym = toke.parseToken(symName);
                 } catch (Exception ex) {
                     throw new SAXException("Couldn't parse token " + symName + " in alphabet " + dist.getAlphabet().getName());
                 }
                 dm.delegate(new DoubleElementHandlerBase() {
                     protected void setDoubleValue(double val)
                         throws SAXException
                     {
                         try {
                             dist.setWeight(sym, val);
                         } catch (Exception ex) {
                             throw new SAXException("Couldn't set weight");
                         }
                     }
                 } );
             }
         }
     }
}
