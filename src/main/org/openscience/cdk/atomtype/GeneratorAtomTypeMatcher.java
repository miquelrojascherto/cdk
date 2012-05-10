/* Copyright (C) 2012  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, version 2.1.
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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 * Support tool for structure generator that can test if a given atom can become
 * a valid atom type. For example, a carbon with three single bonds can become
 * the atom type C.sp2 and C.sp3 (CDK names, and only considering non-charged types).  
 *
 * @author         egonw
 * @cdk.module     atomtype
 * @cdk.githash
 */
@TestClass("org.openscience.cdk.atomtype.GeneratorAtomTypeMatcherTest")
public class GeneratorAtomTypeMatcher {

	private Map<String, List<IAtomType>> types;

	/**
	 * Constructs a new GeneratorAtomTypeMatcher using the atom types defined by the
	 * given AtomTypeFactory.
	 */
	@TestMethod("testConstructorNull,testConstructor")
	public GeneratorAtomTypeMatcher(List<IAtomType> types) {
		if (types == null) {
			this.types = Collections.emptyMap();
			return;
		}

		// make an index with types per element
		this.types = new HashMap<String, List<IAtomType>>();
		for (IAtomType type : types) {
			List<IAtomType> typesForElement = this.types.get(type.getSymbol());
			// create the list if it does not exist yet
			if (typesForElement == null) {
				typesForElement = new ArrayList<IAtomType>();
				this.types.put(type.getSymbol(), typesForElement);
			}
			// add the atom type to the list
			typesForElement.add(type);
		}
	}

	/**
	 * Decides if the given atom can become a valid atom type, by adding a neighbor or
	 * by increasing a bond order.
	 *
	 * @param container {@link IAtomContainer} container in which the atom participates
	 * @param atom      {@link IAtom} to be tested
	 * @return          true, if the atom can become a valid atom type
	 */
	@TestMethod("testPhosphorus")
	public boolean canBecomeValid(IAtomContainer container, IAtom atom) {
		// check all atom types for this atom
		List<IAtomType> typesForElement = this.types.get(atom.getSymbol());
		if (typesForElement == null || typesForElement.size() == 0) {
			// apparently there are not atom types defined for this element
			return false;
		}

		int currentNeighbors = container.getConnectedBondsCount(atom);
		IBond.Order maxBondOrder = container.getMaximumBondOrder(atom);

        for (IAtomType type : typesForElement) {
            if (currentNeighbors <= type.getFormalNeighbourCount() &&
            	!BondManipulator.isHigherOrder(maxBondOrder, type.getMaxBondOrder())) {
                return true; // the first hit is fine.
            }
        }

		return true;
	}

}

