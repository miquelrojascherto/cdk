/* Copyright (C) 2011  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.renderer.visitor;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import junit.framework.Assert;

import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.w3c.dom.Document;
import org.xml.sax.SAXException;

/**
 * @cdk.module  test-rendersvg
 * @cdk.githash
 *
 */
public class SVGGeneratorTest {

	private IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
	
	private StructureDiagramGenerator sdg = new StructureDiagramGenerator();
	
	public IMolecule layout(IMolecule molecule) {
		sdg.setMolecule(molecule);
		try {
			sdg.generateCoordinates();
		} catch (Exception e) {
			System.err.println(e);
		}
		return sdg.getMolecule();
	}

	@Test
	public void testEmptyModel() {
		IMolecule dummy = builder.newInstance(IMolecule.class);
        SVGGenerator svgGenerator = new SVGGenerator();
		List<IGenerator<IAtomContainer>> generators =
				new ArrayList<IGenerator<IAtomContainer>>();
		generators.add(new BasicSceneGenerator());
		generators.add(new BasicBondGenerator());
		BasicAtomGenerator atomGenerator = new BasicAtomGenerator();
		generators.add(atomGenerator);
		
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
		RendererModel model = renderer.getRenderer2DModel();
        renderer.paint(dummy, svgGenerator);
        String svg = svgGenerator.getResult();
		Assert.assertTrue(svg.length()>0);
		
		 // parse an XML document into a DOM tree
	    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	    dbf.setValidating(true);
	    DocumentBuilder parser = null;
		try {
			parser = dbf.newDocumentBuilder();
		} catch (ParserConfigurationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	    try {
	    	Document document = parser.parse (new ByteArrayInputStream (svg.getBytes()));
	    } catch (SAXException e) {
	        Assert.fail();
	    } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


}
