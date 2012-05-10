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

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.silent.SilentChemObjectBuilder;


/**
 * @cdk.module test-core
 */
public class UnchargedNonRadicalBindingFilterTest {

	@Test
	public void testMatches() {
		AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/dict/data/cdk-atom-types.owl", 
	        SilentChemObjectBuilder.getInstance()
		);
		IAtomType[] types = factory.getAllAtomTypes(new UnchargedNonRadicalBindingFilter());
		Assert.assertNotSame(0, types.length);
		for (IAtomType type : types) {
			// no atom type must be uncharged
			Integer charge = type.getFormalCharge();
			Assert.assertTrue(charge == null || charge == 0);
			// no atom type can be a radical
			Integer singleElectronCount = (Integer)type.getProperty(CDKConstants.SINGLE_ELECTRON_COUNT);
			Assert.assertTrue(singleElectronCount == null || singleElectronCount == 0);
		};
	}

}
