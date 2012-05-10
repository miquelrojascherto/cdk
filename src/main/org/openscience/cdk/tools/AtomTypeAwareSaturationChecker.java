/*  Copyright (C) 2012  Klas Jönsson
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
package org.openscience.cdk.tools;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 * This class tries to figure out the bond order of the bonds that has the flag 
 * <code>SINGLE_OR_DOUBLE</code> raised (i.e. set to <code>true</code>).<br>
 * The code is written with the assumption that the properties of the atoms in
 * the molecule has configured with the help of {@link AtomContainerManipulator}.
 * This class uses the {@link SaturationChecker} internally.
 * 
 * @author Klas J&ouml;nsson
 * @author Egon Willighagen
 * @cdk.created 2012-04-13
 * 
 * @cdk.keyword bond order
 * @cdk.module  valency
 */
@TestClass("org.openscience.cdk.tools.ATASaturationCheckerTest")
public class AtomTypeAwareSaturationChecker implements IValencyChecker,
IDeduceBondOrderTool {
	
	SaturationChecker sc;
	private static ILoggingTool logger =
			LoggingToolFactory.createLoggingTool(SaturationChecker.class);
	private IBond.Order oldBondOrder;
	private int startBond;

	/**
	 * Constructs an {@link AtomTypeAwareSaturationChecker} checker.
	 */
	public AtomTypeAwareSaturationChecker() {
		sc = new SaturationChecker();
	}
	
	/**
	 * This method decides the bond order on bonds that are in a ring and has 
	 * the <code>SINGLE_OR_DOUBLE</code>-flag raised.
	 *  
	 * @param atomContainer The molecule to investigate
	 * @throws CDKException 
	 */
	public void decideBondOrder(IAtomContainer atomContainer) throws CDKException {
		startBond = 0;
		decideBondOrder(atomContainer, 0);
	}
	
	/**
	 * This method decides the bond order on bonds that are in a ring and has 
	 * the <code>SINGLE_OR_DOUBLE</code>-flag raised.
	 *  
	 * @param atomContainer The molecule to investigate
	 * @param start The bond to start with
	 * @throws CDKException 
	 */
	private void decideBondOrder(IAtomContainer atomContainer, int start) throws CDKException {	
		for (int i = start; i < atomContainer.getBondCount(); i++) {
			checkBond(atomContainer, i);
		}
	}
	
	/**
	 * This method tries to set the bond order on the bond previous to the current 
	 * bond.
	 * 
	 * @param i The index of the current
	 * @param bond The current bond 
	 * @param atomContainer The molecule
	 * @return The new index, just to ensure it's updated
	 * @throws CDKException when no solution can be found
	 */
	private int backTrack(int i, IBond bond, IAtomContainer atomContainer) throws CDKException {
		bond.setOrder(IBond.Order.SINGLE);
		i--;
		if (i <= 0) {
			startBond++;
			for (IBond resetBond: atomContainer.bonds())
				if (resetBond.getFlag(CDKConstants.SINGLE_OR_DOUBLE))
					resetBond.setOrder(IBond.Order.SINGLE);
			if (startBond == atomContainer.getBondCount()-1 ) {
				/* If we end up here, then we have stepped thru the hole molecule
				 * as many times as it has bond each time starting on a the next 
				 * bond. And still we haven't found any solution. */
				throw new CDKException("Can't find any solution.");
			} else {
				decideBondOrder(atomContainer, startBond);
				if (startBond > 1) {
					/* If there's need to start further in we should go backwards
					 and try to set those bond as well */
					for (int j = startBond; j == 0; j--) {
						checkBond(atomContainer, j);
					}						
				}
			}
		}
		try {
			setMaxBondOrder(bond, atomContainer);
		} catch (CantDecideBondOrderException be) {
			// Need to go back one step more...
			i = backTrack(i, bond, atomContainer);
		}
		
		return i;
	}
	
	/**
	 * This method tries to set the bond order on the current bond.
	 * 
	 * @param atomContainer The molecule
	 * @param index The index of the current bond
	 * @throws CDKException when no suitable solution can be found
	 */
	private void checkBond(IAtomContainer atomContainer, int index) throws CDKException {
		IBond bond = atomContainer.getBond(index);
		if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE) && bond != null ) {
			try {
				oldBondOrder = bond.getOrder();
				bond.setOrder(IBond.Order.SINGLE);
				setMaxBondOrder(bond, atomContainer);
			} catch (CantDecideBondOrderException be) {
				logger.debug(be);
				index = backTrack(index, bond, atomContainer);
			} catch (CDKException e) {
				bond.setOrder(oldBondOrder);			
				logger.debug(e);
			}
		}
	}
	
	/**
	 * This method decides the highest bond order that the bond can have and set 
	 * it to that.
	 * 
	 * @param bond The bond to be investigated
	 * @param atomContainer The {@link IAtomContainer} that contains the bond
	 * @throws CDKEXception when the bond cannot be further increased
	 */
	private void setMaxBondOrder(IBond bond, IAtomContainer atomContainer) throws CDKException {
		if (bondOrderCanBeIncreased(bond, atomContainer)) {
			if (bond.getOrder() != IBond.Order.QUADRUPLE)
				bond.setOrder(BondManipulator.increaseBondOrder(bond.getOrder()));
			else
				throw new CDKException("Can't increase a quadruple bond!");
			}	
	}
	
	/**
	 * Check if the bond order can be increased. This method assumes that the 
	 * bond is between only two atoms.
	 * 
	 * @param bond The bond to check
	 * @param atomContainer The {@link IAtomContainer} that the bond belongs to
	 * @return True if it is possibly to increase the bond order
	 * @throws CDKException 
	 */
	public boolean bondOrderCanBeIncreased(IBond bond, IAtomContainer atomContainer) throws CDKException {
		boolean atom0isUnsaturated = false, atom1isUnsaturated = false;
		if (bondsUsed(bond.getAtom(0), atomContainer) < bond.getAtom(0).getBondOrderSum())
			atom0isUnsaturated = true;
		if (bondsUsed(bond.getAtom(1), atomContainer) < bond.getAtom(1).getBondOrderSum())
			atom1isUnsaturated = true;
		if (atom0isUnsaturated == atom1isUnsaturated)
			return atom0isUnsaturated;
		else {
			/*
			 * If one of the atoms is saturated and the other isn't, what do we do then?
			 * Look at the bonds on each side and decide from that...
			 */
			int myIndex = atomContainer.getBondNumber(bond);
			// If the previous bond is the reason it's no problem, so just move on...
			if (atomContainer.getBond(myIndex-1).getOrder() == IBond.Order.DOUBLE)
				return false;
			/*
			 * The only reason for trouble should now be that the next bond make 
			 * one of the atoms saturated, so lets throw an exception and reveres 
			 * until we can place a double bond and set it as single and continue  
			 */			
			if (isConnected(atomContainer.getBond(myIndex), atomContainer.getBond(0)))
				throw new CantDecideBondOrderException("Can't decide bond order of this bond");
			else {
				return false;
			}
		}
	}

	/**
	 * Look if any atoms in <code>bond1</code> also are in <code>bond2</code>
	 * and if so it conceder the bonds connected.
	 * @param bond1 The first bond
	 * @param bond2 The other bond
	 * @return True if any of  the atoms in <code>bond1</code> also are in <code>bond2</code>
	 */
	private boolean isConnected(IBond bond1, IBond bond2) {
		for (IAtom atom : bond1.atoms())
			if (bond2.contains(atom))
				return true;
		return false;
	}
	
	/**
	 * This method calculates the number of bonds that an <code>IAtom</code>
	 * can have.
	 * 
	 * @param atom The <code>IAtom</code> to be investigated
	 * @return The max number of bonds the <code>IAtom</code> can have
	 * @throws CDKException when the atom's valency is not set
	 */
	public double getMaxNoOfBonds(IAtom atom) throws CDKException {
		double noValenceElectrons = atom.getValency() == CDKConstants.UNSET ? -1 : atom.getValency();
		if (noValenceElectrons == -1) {
			throw new CDKException("Atom property not set: Valency");
		}
		// This will probably only work for group 13-18, and not for helium...
		return 8 - noValenceElectrons;
	}
	
	/**
	 * A small help method that count how many bonds an atom has, regarding 
	 * bonds due to its charge and to implicit hydrogens.
	 * 
	 * @param atom The atom to check
	 * @param atomContainer The atomContainer containing the atom
	 * @return The number of bonds that the atom has
	 */
	private double bondsUsed(IAtom atom, IAtomContainer atomContainer) {
		int bondsToAtom = 0;
		for (IBond bond : atomContainer.bonds())
			if (bond.contains(atom))
				bondsToAtom += BondManipulator.destroyBondOrder(bond.getOrder());
		int implicitHydrogens;
		if (atom.getImplicitHydrogenCount() == CDKConstants.UNSET || atom.getImplicitHydrogenCount() == null) {
			// Will probably only work with group 13-18, and not for helium...
			implicitHydrogens = (8 - atom.getValency()) - atom.getFormalNeighbourCount();
			String warningMessage = "Number of implicite hydrogens not set for atom " + atom.getAtomTypeName() 
					+ ". Estimated it to: " + implicitHydrogens;
			logger.warn(warningMessage);
		} else 
			implicitHydrogens = atom.getImplicitHydrogenCount();
	
		double charge;
		if (atom.getCharge() == CDKConstants.UNSET) 
			if (atom.getFormalCharge() == CDKConstants.UNSET) {
				charge = 0;
				String warningMessage = "Neither charge nor formal charge is set for atom " + atom.getAtomTypeName() 
						+ ". Estimate it to: 0";
				logger.warn(warningMessage);
			} else
				charge = atom.getFormalCharge();
		else
			charge = atom.getCharge();
		return bondsToAtom - charge + implicitHydrogens;
	}
	
	/* FIXME All or some of the methods below should probably be implemented 
	 * in this class and not as now be handled of the SaturateChecker class or
	 * if not used removed*/
	/** {@inheritDoc} */
	@Override
	public void saturate(IAtomContainer ac) throws CDKException {
		sc.saturate(ac);
	}

	/** {@inheritDoc} */
	@Override
	public boolean isSaturated(IAtomContainer ac) throws CDKException {
		return sc.isSaturated(ac);
	}

	/** {@inheritDoc} */
	@Override
	public boolean isSaturated(IAtom atom, IAtomContainer container)
			throws CDKException {
		return sc.isSaturated(atom, container);
	}
	
	/**
	 * This is a private exception thrown when it detects an error and needs to 
	 * start to back-trace.
	 * 
	 * @author Klas J&ouml;nsson
	 * 
	 * TODO Something about the suppressed warning?
	 */
	@SuppressWarnings("serial")
	private class CantDecideBondOrderException extends CDKException {

		public CantDecideBondOrderException(String message) {
			super(message);
		}
		
	}
	
}
