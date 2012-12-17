/* Copyright (C) 2008-2012  Egon Willighagen <egonw@users.sf.net>
 *               2008       Miguel Rojas-Cherto
 *               2011       Rajarshi Guha <rajarshi.guha@gmail.com>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA. 
 */
package org.openscience.cdk.atomtype;

import java.util.HashMap;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 * This class tests the perception of MM2 atom types, which uses
 * CDK atom type perception and mapping of CDK atom types to MM2
 * atom types.
 *
 * @cdk.module test-atomtype
 */
public class CDKBasedMM2AtomTypeMatcherTest extends AbstractMM2AtomTypeTest {

	private static Map<String, Integer> testedAtomTypes = new HashMap<String, Integer>();

	@Test public void testGetInstance_IChemObjectBuilder() {
		IAtomTypeMatcher matcher = CDKBasedMM2AtomTypeMatcher.getInstance(
		    SilentChemObjectBuilder.getInstance());
		Assert.assertNotNull(matcher);
	}
	
	@Test public void testFindMatchingAtomType_IAtomContainer_IAtom() throws Exception {
        IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        final IAtomType.Hybridization thisHybridization = IAtomType.Hybridization.SP3;
        atom.setHybridization(thisHybridization);
        mol.addAtom(atom);

        IAtomTypeMatcher matcher = CDKBasedMM2AtomTypeMatcher.getInstance(
        	SilentChemObjectBuilder.getInstance());
        IAtomType type = matcher.findMatchingAtomType(mol, atom);
        Assert.assertNotNull(type);
        Assert.assertEquals("Csp3", type.getAtomTypeName());
	}

	@Test public void testFindMatchingAtomType_IAtomContainer() throws Exception {
        IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        final IAtomType.Hybridization thisHybridization = IAtomType.Hybridization.SP3;
        atom.setHybridization(thisHybridization);
        mol.addAtom(atom);

        IAtomTypeMatcher matcher = CDKBasedMM2AtomTypeMatcher.getInstance(
        	SilentChemObjectBuilder.getInstance());
        IAtomType[] types = matcher.findMatchingAtomType(mol);
        Assert.assertNotNull(types);
        Assert.assertEquals(1, types.length);
        Assert.assertEquals("Csp3", types[0].getAtomTypeName());
	}

	@Test public void testCsp3() throws Exception {
        IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        final IAtomType.Hybridization thisHybridization = IAtomType.Hybridization.SP3;
        atom.setHybridization(thisHybridization);
        mol.addAtom(atom);

        String[] expectedTypes = {"Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
	}

    @Test public void testEthene() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        IAtom atom2 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addBond(0,1,CDKConstants.BONDORDER_DOUBLE);

        String[] expectedTypes = {"Csp2", "Csp2"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testImine() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        IAtom atom2 = new Atom("N");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addBond(0,1,CDKConstants.BONDORDER_DOUBLE);

        String[] expectedTypes = {"Csp2", "Nimine"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testPropyne() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        IAtom atom2 = new Atom("C");
        IAtom atom3 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addAtom(atom3);
        mol.addBond(0,1,CDKConstants.BONDORDER_TRIPLE);
        mol.addBond(2,1,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"Csp", "Csp", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testFormaldehyde() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("O");
        IAtom atom2 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addBond(0,1,CDKConstants.BONDORDER_DOUBLE);

        String[] expectedTypes = {"Osp2", "Ccarbonyl"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testMethanol() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("O");
        IAtom atom1 = new Atom("H");
        IAtom atom2 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom1);
        mol.addAtom(atom2);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(0,2,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"Osp3", "Halc", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testLithiumMethanoxide() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("O");
        IAtom atom2 = new Atom("C");
        IAtom atom3 = new Atom("Li");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addAtom(atom3);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(0,2,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"Osp3", "Csp3", "LI"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testHCN() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("N");
        IAtom atom2 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addBond(0,1,CDKConstants.BONDORDER_TRIPLE);

        String[] expectedTypes = {"Nsp", "Csp"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testNitromethane() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("N");
        IAtom atom2 = new Atom("O");
        IAtom atom3 = new Atom("O");
        IAtom atom4 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addAtom(atom3);
        mol.addAtom(atom4);
        mol.addBond(0,1,CDKConstants.BONDORDER_DOUBLE);
        mol.addBond(0,2,CDKConstants.BONDORDER_DOUBLE);
        mol.addBond(0,3,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"N2OX", "Osp2", "Osp2", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testMethylAmine() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("N");
        IAtom atom2 = new Atom("C");
        IAtom h1 = new Atom("H");
        IAtom h2 = new Atom("H");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addAtom(h1);
        mol.addAtom(h2);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(0,2,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(0,3,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"Nsp3", "Csp3", "HN", "HN"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testH2S() throws Exception {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        IAtom s = DefaultChemObjectBuilder.getInstance().newInstance(IAtom.class,"S");
        IAtom h1 = DefaultChemObjectBuilder.getInstance().newInstance(IAtom.class,"H");
        IAtom h2 = DefaultChemObjectBuilder.getInstance().newInstance(IAtom.class,"H");

        IBond b1 = DefaultChemObjectBuilder.getInstance().newInstance(IBond.class,s, h1, IBond.Order.SINGLE);
        IBond b2 = DefaultChemObjectBuilder.getInstance().newInstance(IBond.class,s, h2, IBond.Order.SINGLE);

        mol.addAtom(s);
        mol.addAtom(h1);
        mol.addAtom(h2);

        mol.addBond(b1);
        mol.addBond(b2);

        String[] expectedTypes = {"Ssp3", "HS", "HS"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testAmmonia() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("H");
        IAtom atom2 = new Atom("N");
        IAtom atom3 = new Atom("H");
        IAtom atom4 = new Atom("H");
        IAtom atom5 = new Atom("H");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        atom2.setFormalCharge(+1);
        mol.addAtom(atom3);
        mol.addAtom(atom4);
        mol.addAtom(atom5);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(1,2,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(1,3,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(1,4,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"HNplus", "Nsp3", "HNplus", "HNplus", "HNplus"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testTMS() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("C");
        IAtom atom2 = new Atom("Si");
        IAtom atom3 = new Atom("C");
        IAtom atom4 = new Atom("C");
        IAtom atom5 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addAtom(atom3);
        mol.addAtom(atom4);
        mol.addAtom(atom5);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(1,2,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(1,3,CDKConstants.BONDORDER_SINGLE);
        mol.addBond(1,4,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"Csp3", "SI", "Csp3", "Csp3", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testSaltsKations() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        
        IAtom atom = new Atom("Fe");
        atom.setFormalCharge(+2);
        mol.addAtom(atom);
        String[] expectedTypes = new String[]{"FE2"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        atom = new Atom("Fe");
        atom.setFormalCharge(+3);
        mol.addAtom(atom);
        expectedTypes = new String[]{"FE3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("Ni");
        atom.setFormalCharge(+2);
        mol.addAtom(atom);
        expectedTypes = new String[]{"NI2"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("Mg");
        atom.setFormalCharge(+2);
        mol.addAtom(atom);
        expectedTypes = new String[]{"MG"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("Co");
        atom.setFormalCharge(+2);
        mol.addAtom(atom);
        expectedTypes = new String[]{"CO2"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("Co");
        atom.setFormalCharge(+3);
        mol.addAtom(atom);
        expectedTypes = new String[]{"CO3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testFuran() throws Exception {
    	IAtomContainer furan = new AtomContainer();
    	furan.addAtom(new Atom("C"));
    	furan.addAtom(new Atom("C"));
    	furan.addAtom(new Atom("C"));
    	furan.addAtom(new Atom("C"));
    	furan.addAtom(new Atom("O"));
    	furan.addBond(0,1,CDKConstants.BONDORDER_DOUBLE);
    	furan.addBond(1,2,CDKConstants.BONDORDER_SINGLE);
    	furan.addBond(2,3,CDKConstants.BONDORDER_DOUBLE);
    	furan.addBond(3,4,CDKConstants.BONDORDER_SINGLE);
    	furan.addBond(4,0,CDKConstants.BONDORDER_SINGLE);
    	String[] expectedTypes = new String[]{
        	"Csp2","Csp2","Csp2","Csp2","Oar"
        };
        assertAtomTypes(testedAtomTypes, expectedTypes, furan);
    }

    @Test public void testMethane() throws Exception {
    	IAtomContainer mol = new AtomContainer();
        IAtom atom = new Atom("H");
        IAtom atom2 = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);

        String[] expectedTypes = {"H", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testHalogens() throws Exception {
        IAtomContainer mol = new AtomContainer();

        IAtom atom = new Atom("Cl");
        IAtom hydrogen = new Atom("C");
        mol.addAtom(atom);
        mol.addAtom(hydrogen);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        String[] expectedTypes = new String[]{"Cl", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("I");
        mol.addAtom(atom);
        mol.addAtom(hydrogen);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        expectedTypes = new String[]{"I", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("Br");
        mol.addAtom(atom);
        mol.addAtom(hydrogen);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        expectedTypes = new String[]{"Br", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);

        mol = new AtomContainer();
        atom = new Atom("F");
        mol.addAtom(atom);
        mol.addAtom(hydrogen);
        mol.addBond(0,1,CDKConstants.BONDORDER_SINGLE);
        expectedTypes = new String[]{"F", "Csp3"};
        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testNobleGases() throws Exception {
    	IAtomContainer mol = new AtomContainer();

    	mol.addAtom(new Atom("He"));
    	mol.addAtom(new Atom("Ne"));
    	mol.addAtom(new Atom("Ar"));
    	mol.addAtom(new Atom("Kr"));
    	mol.addAtom(new Atom("Xe"));

    	String[] expectedTypes = {"HE", "NE", "AR", "KR", "XE"}; 
    	assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testBoronTetraFluoride() throws Exception {
    	IAtomContainer mol = new AtomContainer();
    	mol.addAtom(new Atom("B")); mol.getAtom(0).setFormalCharge(-1);
    	mol.addAtom(new Atom("F"));
    	mol.addAtom(new Atom("F"));
    	mol.addAtom(new Atom("F"));
    	mol.addAtom(new Atom("F"));
    	mol.addBond(0,1,IBond.Order.SINGLE);
    	mol.addBond(0,2,IBond.Order.SINGLE);
    	mol.addBond(0,3,IBond.Order.SINGLE);
    	mol.addBond(0,4,IBond.Order.SINGLE);
    	
    	String[] expectedTypes = {"Btetrah", "F", "F", "F", "F"}; 
    	assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testBoron() throws Exception {
    	IAtomContainer mol = new AtomContainer();
    	mol.addAtom(new Atom("B"));
    	mol.addAtom(new Atom("F"));
    	mol.addAtom(new Atom("F"));
    	mol.addAtom(new Atom("F"));
    	mol.addBond(0,1,IBond.Order.SINGLE);
    	mol.addBond(0,2,IBond.Order.SINGLE);
    	mol.addBond(0,3,IBond.Order.SINGLE);

    	String[] expectedTypes = {"Btrig", "F", "F", "F"}; 
    	assertAtomTypes(testedAtomTypes, expectedTypes, mol);
    }

    @Test public void testForDuplicateDefinitions() {
    	super.testForDuplicateDefinitions();
    }

    /*
     * This method *must* be the last method in the class.
     */
    @Test public void countTestedAtomTypes() {
        super.countTestedAtomTypes(testedAtomTypes);
    }

}
