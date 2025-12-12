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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.text.SimpleDateFormat;
import java.util.Date;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;
import net.derkholm.nmica.model.motif.SequenceBackground;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Convert old background model files to NestedMICA 0.8 format", generateStub=true)
@NMApp(launchName="nmconvertbg", vm=VirtualMachine.CLIENT)
public class ConvertBG {
	private File in;
	private File out;
	
	
	@Option(help="Old-format background model to process")
	public void setIn(File in) {
		this.in = in;
	}


	@Option(help="Output filename for new-format background model")
	public void setOut(File out) {
		this.out = out;
	}



	public void main(String[] args)
		throws Exception
	{
		Object oldBg;
		try {
			oldBg = new ObjectInputStream(new FileInputStream(in)).readObject();
		} catch (Exception ex) {
			throw new Exception("Coudln't read input file", ex);
		}
		
		if (oldBg instanceof SequenceBackground) {
			if (oldBg instanceof MosaicSequenceBackground) {
				MosaicSequenceBackground msb = (MosaicSequenceBackground) oldBg;
				Mosaic mosaic = new Mosaic(msb.getBackgroundDistributions(), msb.getBackgroundTransition());
		        mosaic.getAnnotation().setProperty("creator.name", "nmconvertbg");
		        mosaic.getAnnotation().setProperty("date", new SimpleDateFormat().format(new Date()));
		        mosaic.getAnnotation().setProperty("input", in.getName());
		        XMLOutputFactory factory = XMLOutputFactory.newInstance();
		        XMLStreamWriter xw = factory.createXMLStreamWriter(new FileOutputStream(out));
		        MosaicIO.writeMosaic(xw, mosaic);
		        xw.close();
			} else {
				System.err.println("Background models of type " + oldBg.getClass().getName() + " are not supported");
			}
		} else {
			System.err.println("Not a NestedMICA background model");
		}
	}
}
