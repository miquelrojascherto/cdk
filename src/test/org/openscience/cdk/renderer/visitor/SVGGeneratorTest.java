/* Copyright (C) 2011,2013  Egon Willighagen <egonw@users.sf.net>
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

import java.awt.Color;
import java.awt.geom.AffineTransform;
import java.io.ByteArrayInputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.elements.TextElement;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.templates.MoleculeFactory;
import org.w3c.dom.Document;

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
	public void testEmptyModel() throws Exception {
		IMolecule dummy = builder.newInstance(IMolecule.class);
        SVGGenerator svgGenerator = new SVGGenerator();
		List<IGenerator<IAtomContainer>> generators =
				new ArrayList<IGenerator<IAtomContainer>>();
		generators.add(new BasicSceneGenerator());
		generators.add(new BasicBondGenerator());
		BasicAtomGenerator atomGenerator = new BasicAtomGenerator();
		generators.add(atomGenerator);
		
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
        renderer.paint(dummy, svgGenerator);
        String svg = svgGenerator.getResult();
		Assert.assertNotSame(0, svg.length());

		validateSVG(svg);
	}

	@Test
	public void testConstructor() {
		IDrawVisitor visitor = new SVGGenerator();
		Assert.assertNotNull(visitor);
	}

	@Test
	public void testSetFontManager() {
		IDrawVisitor visitor = new SVGGenerator();
		visitor.setFontManager(new AWTFontManager());
		// at least we now know it did not crash...
		Assert.assertNotNull(visitor);
	}

	@Test
	public void testSetRendererModel() {
		IDrawVisitor visitor = new SVGGenerator();
		visitor.setRendererModel(new RendererModel());
		// at least we now know it did not crash...
		Assert.assertNotNull(visitor);
	}

    @Test
	public void testVisit() throws Exception {
    	SVGGenerator svgGenerator = new SVGGenerator();
		svgGenerator.setFontManager(new AWTFontManager());
		svgGenerator.setTransform(new AffineTransform());
		svgGenerator.visit(new TextElement(2, 3, "Foo", Color.BLACK));
		// at least we now know it did not crash...
		Assert.assertNotNull(svgGenerator);
		String svg = svgGenerator.getResult();
		Assert.assertNotSame(0, svg.length());
		Assert.assertTrue(svg.contains("<text"));
		Assert.assertTrue(svg.contains("Foo"));

		validateSVG(svg);
	}

    @Test
	public void testRealWorldExample() throws Exception {
    	IMolecule molecule = MoleculeFactory.make124Triazole();
    	molecule = layout(molecule);
    	List<IGenerator<IAtomContainer>> generators =
				new ArrayList<IGenerator<IAtomContainer>>();
		generators.add(new BasicSceneGenerator());
		generators.add(new BasicBondGenerator());
		BasicAtomGenerator atomGenerator = new BasicAtomGenerator();
		generators.add(atomGenerator);
		
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
        
        SVGGenerator svgGenerator = new SVGGenerator();
		svgGenerator.setFontManager(new AWTFontManager());
		svgGenerator.setTransform(new AffineTransform());
        renderer.paint(molecule, svgGenerator);

		// at least we now know it did not crash...
		Assert.assertNotNull(svgGenerator);
		String svg = svgGenerator.getResult();
		Assert.assertNotSame(0, svg.length());

		validateSVG(svg);
	}

    private void validateSVG(String svg) throws Exception {
		 // parse an XML document into a DOM tree
	    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	    dbf.setFeature("http://xml.org/sax/features/namespaces", false);
	    dbf.setFeature("http://xml.org/sax/features/validation", false);
	    dbf.setFeature("http://apache.org/xml/features/nonvalidating/load-dtd-grammar", false);
	    dbf.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
	    DocumentBuilder parser = dbf.newDocumentBuilder();
	    Document document = parser.parse (new ByteArrayInputStream (svg.getBytes()));
	    Assert.assertNotNull(document);
    }

    @Test
    public void testPointTransformations() {
    	double originalX = 5.67;
    	double originalY = 1.23;
    	SVGGenerator svgGenerator = new SVGGenerator();
		svgGenerator.setTransform(new AffineTransform());
    	double[] transformed = svgGenerator.transformPoint(originalX, originalY);
    	double[] newOriginal = svgGenerator.invTransformPoint(transformed[0], transformed[1]);
    	Assert.assertEquals(originalX, newOriginal[0], 0.01);
    	Assert.assertEquals(originalY, newOriginal[1], 0.01);
    }
}
