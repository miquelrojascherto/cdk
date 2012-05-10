/*  Copyright (C) 2012  Klas Jönsson
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.tools;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomType;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

/**
 * A test class for the AtomTypeAwereSaturationChecker-class
 *
 * @author Klas J&ouml;nsson
 * @cdk.created 2012-04-18
 * @cdk.module  test-valency
 */
public class ATASaturationCheckerTest extends org.openscience.cdk.CDKTestCase {

	SaturationChecker satcheck = null;
	boolean standAlone = false;
	private static SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
	private AtomTypeAwareSaturationChecker atasc = new AtomTypeAwareSaturationChecker();

	
	@Before
    @Test
    public void setUp() throws Exception {
    	satcheck = new SaturationChecker();
    }
	
	/**
	 * A test of decideBondOrder(IAtomContainer) with a molecule we created
	 * from scratch.
	 * @throws Exception 
	 */
	@Test
	public void testASimpleCarbonRing() throws Exception {
		// First we create a simple carbon ring to play with... 
		IMolecule mol = new Molecule();
		IAtomType carbon = new AtomType(Elements.CARBON);
		
		IAtom a0 = new Atom("C");
		a0.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a0, carbon);
		IAtom a1 = new Atom("C");
		a1.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a1, carbon);
		IAtom a2 = new Atom("C");
		a2.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a2, carbon);
		IAtom a3 = new Atom("C");
		a3.setHybridization(IAtomType.Hybridization.SP2);

		AtomTypeManipulator.configureUnsetProperties(a3, carbon);
		IAtom a4 = new Atom("C");
		a4.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a4, carbon);
		IAtom a5 = new Atom("C");
		a5.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a5, carbon);
		
		mol.addAtom(a0);
		mol.addAtom(a1);
		mol.addAtom(a2);
		mol.addAtom(a3);
		mol.addAtom(a4);
		mol.addAtom(a5);
		
		IBond b0 = new Bond(a0,a1);
		b0.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b0);
		IBond b1 = new Bond(a1,a2);
		b1.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b1);
		IBond b2 = new Bond(a2,a3);
		b2.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b2);
		IBond b3 = new Bond(a3,a4);
		b3.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b3);
		IBond b4 = new Bond(a4,a5);
		b4.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b4);
		IBond b5 = new Bond(a5,a0);
		b5.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b5);

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		AtomTypeTools att = new AtomTypeTools();
		att.assignAtomTypePropertiesToAtom(mol, false);

		// ...then we send it to the method we want to test...
		atasc.decideBondOrder(mol);
		
		Assert.assertEquals(IBond.Order.DOUBLE, b0.getOrder());
		Assert.assertEquals(IBond.Order.SINGLE, b1.getOrder());
		Assert.assertEquals(IBond.Order.DOUBLE, b2.getOrder());
		Assert.assertEquals(IBond.Order.SINGLE, b3.getOrder());
		Assert.assertEquals(IBond.Order.DOUBLE, b4.getOrder()); 
		Assert.assertEquals(IBond.Order.SINGLE, b5.getOrder());
	}
	
	/**
	 * A test of decideBondOrder(IAtomContainer) with a molecule we created
	 * from a SMILES.
	 * @throws Exception 
	 */
	@Test
	public void testQuinone() throws Exception {
		IAtomContainer mol = sp.parseSmiles("O=c1ccc(=O)cc1");
		
		atasc.decideBondOrder(mol);
		
		Assert.assertTrue(mol.getAtom(1).getHybridization() == IAtomType.Hybridization.SP2);
		
		Assert.assertTrue(mol.getBond(0).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(0).getAtom(1).getSymbol().equals("O"));
		Assert.assertEquals(mol.getBond(0).getOrder(), IBond.Order.DOUBLE);
		
		Assert.assertTrue(mol.getBond(1).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(1).getAtom(1).getSymbol().equals("C"));
		Assert.assertEquals(mol.getBond(1).getOrder(), IBond.Order.SINGLE);

		Assert.assertTrue(mol.getBond(2).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(2).getAtom(1).getSymbol().equals("C"));
		Assert.assertEquals(mol.getBond(2).getOrder(), IBond.Order.DOUBLE);
		
		Assert.assertTrue(mol.getBond(3).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(3).getAtom(1).getSymbol().equals("C"));
		Assert.assertEquals(mol.getBond(3).getOrder(), IBond.Order.SINGLE);
		
		Assert.assertTrue(mol.getBond(4).getAtom(0).getSymbol().equals("O"));
		Assert.assertTrue(mol.getBond(4).getAtom(1).getSymbol().equals("C"));
		Assert.assertEquals(mol.getBond(4).getOrder(), IBond.Order.DOUBLE);
		
		Assert.assertTrue(mol.getBond(5).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(5).getAtom(1).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(5).getOrder() == IBond.Order.SINGLE);
		
		Assert.assertTrue(mol.getBond(6).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(6).getAtom(1).getSymbol().equals("C"));
		Assert.assertEquals(mol.getBond(6).getOrder(), IBond.Order.DOUBLE);
		
		Assert.assertTrue(mol.getBond(7).getAtom(0).getSymbol().equals("C"));
		Assert.assertTrue(mol.getBond(7).getAtom(1).getSymbol().equals("C"));
		Assert.assertEquals(mol.getBond(7).getOrder(), IBond.Order.SINGLE);
	
		Assert.assertEquals(mol.getBond(0).getAtom(0), mol.getBond(7).getAtom(1));
	}
	
	/**
	 * A test of decideBondOrder(IAtomContainer) with a simple carbon ring we 
	 * created from a SMILES.
	 * @throws CDKException 
	 */
	@Test
	public void testASimpleCarbonRing2() throws CDKException {
		IAtomContainer mol = sp.parseSmiles("c1ccccc1");
		
		atasc.decideBondOrder(mol);
		
		Assert.assertEquals(mol.getAtom(1).getHybridization(), IAtomType.Hybridization.SP2);
		Assert.assertEquals(mol.getBond(0).getOrder(), IBond.Order.DOUBLE);
		Assert.assertEquals(mol.getBond(1).getOrder(), IBond.Order.SINGLE);
		Assert.assertEquals(mol.getBond(2).getOrder(), IBond.Order.DOUBLE);
		Assert.assertEquals(mol.getBond(3).getOrder(), IBond.Order.SINGLE);
		Assert.assertEquals(mol.getBond(4).getOrder(), IBond.Order.DOUBLE);
		Assert.assertEquals(mol.getBond(5).getOrder(), IBond.Order.SINGLE);
	}
	
	/**
	 * This method tests the AtomTypeAwareSaturnationChecker with a large ring 
	 * system.
	 * @throws Exception
	 */
	@Test
    public void testALargeRingSystem() throws Exception {
        String smiles = "O=C1Oc6ccccc6(C(O)C1C5c2ccccc2CC(c3ccc(cc3)c4ccccc4)C5)";
        IAtomContainer mol = sp.parseSmiles(smiles);
        
        atasc.decideBondOrder(mol);
        
        // We should have 13 double bonds.
        int doubleBondCount = 0;
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) doubleBondCount++;
        }
        Assert.assertEquals(13, doubleBondCount);
	}
	
	/**
	 * This do the same as the method above, but with five other large ring 
	 * systems.
	 * @throws Exception
	 */
	@Test
    public void testLargeRingSystem1() throws Exception {
		// Should have 6 double bonds
		String smiles = "c1ccc2c(c1)CC4NCCc3cccc2c34";
        IAtomContainer mol = sp.parseSmiles(smiles);    
        atasc.decideBondOrder(mol);
        
        int doubleBondCount = 0; 
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) doubleBondCount++;
        }
        Assert.assertEquals(6, doubleBondCount);
	}
	
	@Test
	public void testLargeRingSystem2() throws Exception {
		// Should have 8 double bonds
		String smiles = "Oc1ccc(cc1)c1coc2c(c1=O)c(O)cc(c2)O";
        IAtomContainer mol = sp.parseSmiles(smiles);
        atasc.decideBondOrder(mol);
        
        int doubleBondCount = 0;
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) doubleBondCount++;
        }
        Assert.assertEquals(8, doubleBondCount);
	}
	
	@Test
	public void testADoubleRingWithANitrogenAtom() throws Exception {
		// Should have 8 double bonds
		String smiles = "c1ccn2cccc2c1";
        IAtomContainer mol = sp.parseSmiles(smiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(mol);
        
        IAtom nitrogen = mol.getAtom(3);

        atasc.decideBondOrder(mol);
        
        
        int doubleBondCount = 0, nSingleBonds = 0;
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) 
            	doubleBondCount++;
            if (bond.contains(nitrogen))
            	if (bond.getOrder().equals(Order.SINGLE))
            		nSingleBonds++;
        }
        Assert.assertEquals(4, doubleBondCount);
        Assert.assertEquals(3, nSingleBonds);
	}
	
	@Test
	public void testLargeRingSystem3() throws Exception {
		// Should have 17 double bonds
		String smiles = "O=C5C=C(O)C(N=Nc1ccc(cc1)Nc2ccccc2)=CC5(=NNc3ccc(cc3)Nc4ccccc4)";
        IAtomContainer mol = sp.parseSmiles(smiles);
        atasc.decideBondOrder(mol);
        
        int doubleBondCount = 0;
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) doubleBondCount++;
        }
        Assert.assertEquals(17, doubleBondCount);
	}
	
	@Test
	public void testLargeRingSystem4() throws Exception {
		// Should have 18 double bonds
		String smiles = "c1ccc(cc1)[Sn](c2ccccc2)(c3ccccc3)S[Sn](c4ccccc4)(c5ccccc5)c6ccccc6";
        IAtomContainer mol = sp.parseSmiles(smiles);
        atasc.decideBondOrder(mol);
        
        int doubleBondCount = 0;
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) doubleBondCount++;
        }
        Assert.assertEquals(18, doubleBondCount);
	}
	
	@Test
	public void testLargeRingSystem5() throws Exception {
		// Should have 24 double bonds
		String smiles = "O=C1c2ccccc2C(=O)c3c1ccc4c3[nH]c5c6C(=O)c7ccccc7C(=O)c6c8[nH]c9c%10C(=O)c%11ccccc%11C(=O)c%10ccc9c8c45";
		
        IAtomContainer mol = sp.parseSmiles(smiles);
        atasc.decideBondOrder(mol);
        
        int doubleBondCount = 0;
        for (IBond bond : mol.bonds()) {
            if (bond.getOrder().equals(Order.DOUBLE)) doubleBondCount++;
        }
        Assert.assertEquals(24, doubleBondCount);
	}
	
	/**
	 * From DeduceBondSystemToolTest
	 * 
	 * @throws Exception
	 */
	@Test
    public void testLargeBioclipseUseCase() throws Exception {
        String smiles = "COc1ccc2[C@@H]3[C@H](COc2c1)C(C)(C)OC4=C3C(=O)C(=O)C5=C4OC(C)(C)[C@@H]6COc7cc(OC)ccc7[C@H]56";
        IMolecule molecule = sp.parseSmiles(smiles);

        atasc.decideBondOrder(molecule);

        // we should have 14 double bonds
        int doubleBondCount = 0;
        for (int i = 0; i < molecule.getBondCount(); i++) {
            IBond bond = molecule.getBond(i);
            if (bond.getOrder() == Order.DOUBLE) doubleBondCount++;
        }
        Assert.assertEquals(10, doubleBondCount);
    }
	
	@Test
	public void testCyclobutadiene() throws CDKException {
		IAtomContainer mol = sp.parseSmiles("c1ccc1");

		atasc.decideBondOrder(mol);
		Assert.assertEquals(mol.getAtom(1).getHybridization(),
				IAtomType.Hybridization.SP2);
		Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(0).getOrder());
		Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(1).getOrder());
		Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(2).getOrder());
		Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(3).getOrder());
	}

    @Test
    public void testPyrrole() throws CDKException {
    	IAtomContainer mol = sp.parseSmiles("c1c[nH]cc1");

    	atasc.decideBondOrder(mol);
    	Assert.assertEquals(mol.getAtom(1).getHybridization(),
    			IAtomType.Hybridization.SP2);
    	Assert.assertTrue(mol.getBond(0).getFlag(CDKConstants.SINGLE_OR_DOUBLE));
    	Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(0).getOrder());
    	Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(1).getOrder());
    	Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(2).getOrder());
    	Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(3).getOrder());
    	Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(4).getOrder());
    }

    @Test
    public void testFurane() throws CDKException {
    	IAtomContainer mol = sp.parseSmiles("c1cocc1");
    	AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(mol);
    	
    	atasc.decideBondOrder(mol);
    	Assert.assertEquals(mol.getAtom(1).getHybridization(),
    			IAtomType.Hybridization.SP2);
    	Assert.assertTrue(mol.getBond(0).getFlag(CDKConstants.SINGLE_OR_DOUBLE));
    	Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(0).getOrder());
    	Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(1).getOrder());
    	//bond to oxygen
    	Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(2).getOrder());
    	//bond to oxygen
    	Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(3).getOrder());
    	Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(4).getOrder());
    }

    @Test
    public void testAnOtherDoubleRing() throws CDKException {
    	IAtomContainer mol = sp.parseSmiles("c1cccc2cccc2c1");

    	atasc.decideBondOrder(mol);
    	Assert.assertEquals(mol.getAtom(1).getHybridization(),
    			IAtomType.Hybridization.SP2);
    	int doubleBondCount = 0;
    	for (IBond bond : mol.bonds()) {
    		Assert.assertTrue(bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE));
    		if (bond.getOrder() == IBond.Order.DOUBLE) doubleBondCount++;
    	}
    	Assert.assertEquals(5, doubleBondCount);
    }
    
    @Test
    public void testAnOtherRingSystem() throws CDKException {
    	IAtomContainer mol = sp.parseSmiles("O=c2c1ccccc1c3ccccc23");

    	atasc.decideBondOrder(mol);
    	Assert.assertEquals(mol.getAtom(1).getHybridization(),
    			IAtomType.Hybridization.SP2);
    	int doubleBondCount = 0;
    	for (IBond bond : mol.bonds()) {
    		if (bond.getOrder() == IBond.Order.DOUBLE) doubleBondCount++;
    	}
    	Assert.assertEquals(7, doubleBondCount);
    }
    
	@Test
	public void testAzulene() throws CDKException {
		IAtomContainer mol = sp.parseSmiles("c12c(ccccc2)ccc1");
		atasc.decideBondOrder(mol);

    	int doubleBondCount = 0;
    	for (IBond bond : mol.bonds()) {
    		if (bond.getOrder() == IBond.Order.DOUBLE) doubleBondCount++;
    	}
    	Assert.assertEquals(5, doubleBondCount);
    	
		Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(0).getOrder());
	}
	
	@Test
	public void testButadiene() throws Exception {
        IMolecule mol = new Molecule();
		IAtomType carbon = new AtomType(Elements.CARBON);
		
		IAtom a0 = new Atom("C");
		a0.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a0, carbon);
		IAtom a1 = new Atom("C");
		a1.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a1, carbon);
		IAtom a2 = new Atom("C");
		a2.setHybridization(IAtomType.Hybridization.SP2);
		AtomTypeManipulator.configureUnsetProperties(a2, carbon);
		IAtom a3 = new Atom("C");
		a3.setHybridization(IAtomType.Hybridization.SP2);
		
		mol.addAtom(a0);
		mol.addAtom(a1);
		mol.addAtom(a2);
		mol.addAtom(a3);
		
		IBond b0 = new Bond(a0,a1);
		b0.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b0);
		IBond b1 = new Bond(a1,a2);
		b1.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b1);
		IBond b2 = new Bond(a2,a3);
		b2.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
		mol.addBond(b2);
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		AtomTypeTools att = new AtomTypeTools();
		att.assignAtomTypePropertiesToAtom(mol, false);
        
        atasc.decideBondOrder(mol);

		Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(0).getOrder());
		Assert.assertEquals(IBond.Order.SINGLE, mol.getBond(1).getOrder());
		Assert.assertEquals(IBond.Order.DOUBLE, mol.getBond(2).getOrder());
	}
}
