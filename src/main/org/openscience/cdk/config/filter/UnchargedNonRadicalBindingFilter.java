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
package org.openscience.cdk.config.filter;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.interfaces.IAtomType;

/**
 * {@link IAtomTypeFilter} that only matches uncharged, non-radical, but binding atoms.
 *
 * @cdk.module core
 */
@TestClass("org.openscience.cdk.config.filter.UnchargedNonRadicalFilterTest")
public class UnchargedNonRadicalBindingFilter implements IAtomTypeFilter {

	/**
	 * Filter matches that indicates if the passed {@link IAtomType} passes the filter.
	 * 
	 * @param  type {@link IAtomType} to test.
	 * @return      true if the atom type passes the filter, false otherwise
	 */
	@TestMethod("testMatches")
	public boolean passes(IAtomType type) {
		// check the formal charge
		if (type.getFormalCharge() != null &&
			type.getFormalCharge() != 0) return false;

		// check for being a radical
		Object singleElectronCount = type.getProperty(CDKConstants.SINGLE_ELECTRON_COUNT);
		if (singleElectronCount != null &&
			(Integer)singleElectronCount != 0) return false;

		// check the number of neighbors
		if (type.getFormalNeighbourCount() != null &&
			type.getFormalNeighbourCount() == 0) return false;

		return true;
	}

}
