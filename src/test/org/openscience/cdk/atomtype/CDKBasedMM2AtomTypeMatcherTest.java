/* Copyright (C) 2008-2012  Egon Willighagen <egonw@users.sf.net>
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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
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
