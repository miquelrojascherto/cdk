/* Copyright (C) 2012  Egon Willighagen <egonw@users.sf.net>
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

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.AtomContainer;

/**
 * @cdk.module test-atomtype
 */
public class GeneratorAtomTypeMatcherTest {

	@Test
	public void testConstructor() {
        AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/dict/data/cdk-atom-types.owl", 
            new ChemObject().getBuilder()
        );
        IAtomType[] types = factory.getAtomTypes("C");
        List<IAtomType> typeList = new ArrayList<IAtomType>();
        for (IAtomType type : types) {
        	typeList.add(type);
        }
		GeneratorAtomTypeMatcher matcher = new GeneratorAtomTypeMatcher(typeList);
		IAtomContainer container = new AtomContainer();
		IAtom atom = new Atom("C");
		container.addAtom(atom);
		Assert.assertTrue(matcher.canBecomeValid(container, atom));
	}

	@Test
	public void testConstructorNull() {
		// nothing will match when a null atom type list is passed
		GeneratorAtomTypeMatcher matcher = new GeneratorAtomTypeMatcher(null);
		IAtomContainer container = new AtomContainer();
		IAtom atom = new Atom("C");
		container.addAtom(atom);
		Assert.assertFalse(matcher.canBecomeValid(container, atom));
	}

	@Test
	public void testPhosporus() {
        AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/dict/data/cdk-atom-types.owl", 
            new ChemObject().getBuilder()
        );
        IAtomType[] types = factory.getAllAtomTypes();
        List<IAtomType> typeList = new ArrayList<IAtomType>();
        for (IAtomType type : types) {
        	typeList.add(type);
        }
        GeneratorAtomTypeMatcher matcher = new GeneratorAtomTypeMatcher(typeList);
        IAtomContainer container = new AtomContainer();
        IAtom atom = new Atom("P");
        container.addAtom(atom);
        container.addAtom(new Atom("O"));
        container.addAtom(new Atom("O"));
        container.addAtom(new Atom("O"));
        container.addAtom(new Atom("O"));
        container.addBond(0, 1, Order.SINGLE);
        container.addBond(0, 2, Order.SINGLE);
        container.addBond(0, 3, Order.SINGLE);
        container.addBond(0, 4, Order.SINGLE);
        Assert.assertTrue(matcher.canBecomeValid(container, atom));
	}
}
