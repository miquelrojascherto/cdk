/* Copyright (C) 2007-2013  Egon Willighagen <egonw@users.sf.net>
 *                    2011  Nimish Gopal <nimishg@ebi.ac.uk>
 *                    2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *                    2011  Gilleain Torrance <gilleain.torrance@gmail.com>
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

import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.NoSuchAtomException;
import org.openscience.cdk.graph.SpanningTree;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IAtomType.Hybridization;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.tools.manipulator.BondManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 * Atom Type matcher that perceives atom types as defined in the CDK atom type list
 * <code>org/openscience/cdk/dict/data/cdk-atom-types.owl</code>. 
 * If there is not an atom type defined for the tested atom, then NULL 
 * is returned.
 *
 * @author         egonw
 * @cdk.created    2007-07-20
 * @cdk.module     core
 * @cdk.githash
 * @cdk.bug        1802998
 */
@TestClass("org.openscience.cdk.atomtype.CDKAtomTypeMatcherTest")
public class AtomTypeMatcher {

	public final static int REQUIRE_NOTHING = 1;
	public final static int REQUIRE_EXPLICIT_HYDROGENS = 2;
	
	private AtomTypeFactory factory;
	private int mode;
	private IAtomContainer atomContainer = null;
	private IBond.Order[] maxBondOrders = null;
	private int[] neighborCounts = null;
	private int[] piBondCounts = null;
	private int[] singleBondCounts = null;
	private SpanningTree st;
	private boolean[] isRingAtom = null;

    private AtomTypeMatcher(IAtomContainer container, int mode) {
    	atomContainer = container;
    	factory = AtomTypeFactory.getInstance(
			"org/openscience/cdk/dict/data/cdk-atom-types.owl",
			atomContainer.getBuilder()
		);
    	this.mode = mode;
    	cacheProperties();
    }
    
    @TestMethod("testGetInstance_IChemObjectBuilder")
    public static AtomTypeMatcher getInstance(IAtomContainer container) {
        return getInstance(container, REQUIRE_NOTHING);
    }

    @TestMethod("testGetInstance_IChemObjectBuilder_int")
    public static AtomTypeMatcher getInstance(IAtomContainer container, int mode) {
    	return new AtomTypeMatcher(container, mode);
    }
    
    @TestMethod("testFindMatchingAtomType_IAtomContainer")
    public IAtomType[] findMatchingAtomTypes() throws CDKException {
        IAtomType[] types = new IAtomType[atomContainer.getAtomCount()];
        int typeCounter = 0;
        for (IAtom atom : atomContainer.atoms()) {
            types[typeCounter] = findMatchingAtomType(atom);
            typeCounter++;
        }
        return types;
    }

    @TestMethod("testFindMatchingAtomType_IAtomContainer_IAtom")
    public IAtomType findMatchingAtomType(IAtom atom)
        throws CDKException {
        IAtomType type = null;
        if (atom instanceof IPseudoAtom) {
        	return factory.getAtomType("X");
        }
        Integer atomicNumberObj = atom.getAtomicNumber();
        String symbol = atom.getSymbol(); 
        if (atomicNumberObj == null) atomicNumberObj =  PeriodicTable.getAtomicNumber(symbol);
        if (atomicNumberObj == null) return null; // e.g. unknown element symbols
        int atomicNumber = atomicNumberObj.intValue();
        
        int atomNumber = this.atomContainer.getAtomNumber(atom);

        switch (atomicNumber) {
        case 1: type = perceiveHydrogens(atomNumber); break;
        case 2: type = perceiveHelium(atomNumber); break;
        case 3: type = perceiveLithium(atomNumber); break;
        case 4: type = perceiveBeryllium(atomNumber); break;
        case 5: type = perceiveBorons(atomNumber); break;
        case 6: type = perceiveCarbons(atomNumber); break;
        case 7: type = perceiveNitrogens(atomNumber); break;
        case 8: type = perceiveOxygens(atomNumber); break;
        case 9: type = perceiveFluors(atomNumber); break;
        case 10: type = perceiveNeon(atomNumber); break;
        case 11: type = perceiveSodium(atomNumber); break;
        case 12: type = perceiveMagnesium(atomNumber); break;
        case 13: type = perceiveAluminium(atomNumber); break;
        case 14: type = perceiveSilicon(atomNumber); break;
        case 15: type = perceivePhosphors(atomNumber); break;
        case 16: type = perceiveSulphurs(atomNumber); break;
        case 17: type = perceiveChlorine(atomNumber); break;
        case 18: type = perceiveArgon(atomNumber); break;
        case 19: type = perceivePotassium(atomNumber); break;
        case 20: type = perceiveCalcium(atomNumber); break;
        case 34: type = perceiveSelenium(atomNumber); break;
        case 35: type = perceiveBromine(atomNumber); break;
        case 52: type = perceiveTellurium(atomNumber); break;
        case 53: type = perceiveIodine(atomNumber); break;
        case 84: type = perceivePolodium(atomNumber); break;
        }
        if (type != null) return type;
        switch (symbol) {
        case "Cr": type = perceiveChromium(atomNumber); break;
        case "Mo": type = perceiveMolybdenum(atomNumber); break;
        case "Rb": type = perceiveRubidium(atomNumber); break;
        case "Cu": type = perceiveCopper(atomNumber); break;
        case "Ba": type = perceiveBarium(atomNumber); break;
        case "Ga": type = perceiveGallium(atomNumber); break;
        case "Ru": type = perceiveRuthenium(atomNumber); break;
        case "Zn": type = perceiveZinc(atomNumber); break;
        case "Ni": type = perceiveNickel(atomNumber); break;
        case "Gd": type = perceiveGadolinum(atomNumber); break;
        case "Ge": type = perceiveGermanium(atomNumber); break;
        case "Co": type = perceiveCobalt(atomNumber); break;
        case "Sb": type = perceiveAntimony(atomNumber); break;
        case "V": type = perceiveVanadium(atomNumber); break;
        case "Ti": type = perceiveTitanium(atomNumber); break;
        case "Sr": type = perceiveStrontium(atomNumber); break;
        case "Pb": type = perceiveLead(atomNumber); break;
        case "Tl": type = perceiveThallium(atomNumber); break;
        case "Pt": type = perceivePlatinum(atomNumber); break;
        case "Hg": type = perceiveMercury(atomNumber); break;
        case "Fe": type = perceiveIron(atomNumber); break;
        case "Ra": type = perceiveRadium(atomNumber); break;
        case "Au": type = perceiveGold(atomNumber); break;
        case "Ag": type = perceiveSilver(atomNumber); break;
        case "In": type = perceiveIndium(atomNumber); break;
        case "Pu": type = perceivePlutonium(atomNumber); break;
        case "Th": type = perceiveThorium(atomNumber); break;
        case "Mn": type = perceiveManganese(atomNumber); break;
        case "As": type = perceiveArsenic(atomNumber); break;
        case "Cd": type = perceiveCadmium(atomNumber); break;
        case "W": type = perceiveTungstun(atomNumber); break;
        case "Kr": type = perceiveKrypton(atomNumber); break;
        case "Xe": type = perceiveXenon(atomNumber); break;
        case "Rn": type = perceiveRadon(atomNumber); break;
        case "Sn": type = perceiveStannum(atomNumber); break;
        case "Sc": type = perceiveScandium(atomNumber); break;
        }
        return type;
    }
    
    private IAtomType perceiveGallium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
        if (!isCharged(atom) && maxBondOrder == IBond.Order.SINGLE &&
            getConnectedAtomsCount(atomNumber) <= 3) {
            IAtomType type = getAtomType("Ga");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() == 3) {
            IAtomType type = getAtomType("Ga.3plus");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }

    private IAtomType perceiveGermanium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
        if (!isCharged(atom) && maxBondOrder == IBond.Order.SINGLE &&
            getConnectedAtomsCount(atomNumber) <= 4) {
            IAtomType type = getAtomType("Ge");
            if (isAcceptable(atomNumber, type)) return type;
        }
        if (atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 3) {
            IAtomType type = getAtomType("Ge.3");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }

    private IAtomType perceiveSelenium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	int doublebondcount = countAttachedDoubleBonds(atomNumber);
    	if (atom.getFormalCharge() != CDKConstants.UNSET
    			&& atom.getFormalCharge() == 0) {
    		if (getConnectedAtomsCount(atomNumber) == 0) {
    			if (atom.getImplicitHydrogenCount() != null && atom.getImplicitHydrogenCount() == 0 ) {
    				IAtomType type = getAtomType("Se.2");
    				if (isAcceptable(atomNumber, type)) return type;
    			} else {
    				IAtomType type = getAtomType("Se.3");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		} else if (getConnectedAtomsCount(atomNumber) == 1) {

    			if (doublebondcount == 1) {
    				IAtomType type = getAtomType("Se.1");
    				if (isAcceptable(atomNumber, type)) return type;
    			} else if (doublebondcount == 0) {
    				IAtomType type = getAtomType("Se.3");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		} else if (getConnectedAtomsCount(atomNumber) == 2) {
    			if (doublebondcount == 0) {
    				IAtomType type = getAtomType("Se.3");
    				if (isAcceptable(atomNumber, type)) return type;
    			} else if (doublebondcount == 2) {
    				IAtomType type = getAtomType("Se.sp2.2");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		} else if (getConnectedAtomsCount(atomNumber) == 3) {
    			IAtomType type = getAtomType("Se.sp3.3");
    			if (isAcceptable(atomNumber, type)) return type;
    		} else if (getConnectedAtomsCount(atomNumber) == 4) {
    			if (doublebondcount == 2) {
    				IAtomType type = getAtomType("Se.sp3.4");
    				if (isAcceptable(atomNumber, type)) return type;
    			} else if (doublebondcount == 0) {
    				IAtomType type = getAtomType("Se.sp3d1.4");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		} else if (getConnectedAtomsCount(atomNumber) == 5) {
    			IAtomType type = getAtomType("Se.5");
    			if (isAcceptable(atomNumber, type)) return type;
    		}
    	} else if ((atom.getFormalCharge() != CDKConstants.UNSET && atom.getFormalCharge() == 4)
    			&& getConnectedAtomsCount(atomNumber) == 0) {
    		IAtomType type = getAtomType("Se.4plus");
    		if (isAcceptable(atomNumber, type)) return type;
    	} else if ((atom.getFormalCharge() != CDKConstants.UNSET && atom.getFormalCharge() == 1)
    			&& getConnectedAtomsCount(atomNumber) == 3) {
    		IAtomType type = getAtomType("Se.plus.3");
    		if (isAcceptable(atomNumber, type)) return type;
    	} else if ((atom.getFormalCharge() != CDKConstants.UNSET && atom.getFormalCharge() == -2)
    			&& getConnectedAtomsCount(atomNumber) == 0) {
    		IAtomType type = getAtomType("Se.2minus");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
        return null;
    }

    private IAtomType perceiveTellurium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
        if (!isCharged(atom) && maxBondOrder == IBond.Order.SINGLE && getConnectedAtomsCount(atomNumber) <= 2) {
            IAtomType type = getAtomType("Te.3");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() == 4) {
            if (getConnectedAtomsCount(atomNumber) == 0) {
                IAtomType type = getAtomType("Te.4plus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }

	private IAtomType perceiveBorons(int atomNumber)
		throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
	    IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	    if (atom.getFormalCharge() == -1 && 
	        maxBondOrder == IBond.Order.SINGLE &&
	        getConnectedAtomsCount(atomNumber) <= 4) {
	        IAtomType type = getAtomType("B.minus");
	        if (isAcceptable(atomNumber, type)) return type;
	    } else if (atom.getFormalCharge() == +3
                && getConnectedAtomsCount(atomNumber) == 4) {
            IAtomType type = getAtomType("B.3plus");
            if (isAcceptable(atomNumber, type)) return type;
	    } else if (getConnectedAtomsCount(atomNumber) <= 3) {
	        IAtomType type = getAtomType("B");
	        if (isAcceptable(atomNumber, type)) return type;
	    }
    	return null;
    }

    private IAtomType perceiveBeryllium(int atomNumber)
    	throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
		if (atom.getFormalCharge() == -2 &&
		    getMaximumBondOrder(atomNumber) == IBond.Order.SINGLE &&
		    getConnectedAtomsCount(atomNumber) <= 4) {
		    IAtomType type = getAtomType("Be.2minus");
		    if (isAcceptable(atomNumber, type)) return type;
		} else if (atom.getFormalCharge() == 0 &&
                   getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Be.neutral");
            if (isAcceptable(atomNumber, type)) return type;
        }
		return null;
    }

    private IAtomType perceiveCarbonRadicals(int atomNumber) throws CDKException {
        if (getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("C.radical.planar");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (getConnectedAtomsCount(atomNumber) <= 3) {
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (maxBondOrder == IBond.Order.SINGLE) {
                IAtomType type = getAtomType("C.radical.planar");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (maxBondOrder == IBond.Order.DOUBLE) {
                IAtomType type = getAtomType("C.radical.sp2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (maxBondOrder == IBond.Order.TRIPLE) {
                IAtomType type = getAtomType("C.radical.sp1");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }
    
	private IAtomType perceiveCarbons(int atomNumber)
    	throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
		boolean isCharged = isCharged(atom);
	    // if hybridization is given, use that
	    if (hasOneSingleElectron(atomNumber)) {
	        return perceiveCarbonRadicals(atomNumber);
	    } else if (hasHybridization(atom) && !isCharged) {
	        if (atom.getHybridization() == Hybridization.SP2) {
	            IAtomType type = getAtomType("C.sp2");
	            if (isAcceptable(atomNumber, type)) return type;
	        } else if (atom.getHybridization() == Hybridization.SP3) {
	            IAtomType type = getAtomType("C.sp3");
	            if (isAcceptable(atomNumber, type)) return type;
	        } else if (atom.getHybridization() == Hybridization.SP1) {
	        	IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	        	if (maxBondOrder == Order.TRIPLE) {
		            IAtomType type = getAtomType("C.sp");
		            if (isAcceptable(atomNumber, type)) return type;
	        	} else {
	        		IAtomType type = getAtomType("C.allene");
	        		if (isAcceptable(atomNumber, type)) return type;
	        	}
	        }
	    } else if (atom.getFlag(CDKConstants.ISAROMATIC)) {
	        IAtomType type = getAtomType("C.sp2");
	        if (isAcceptable(atomNumber, type)) return type;
	    } else if (hasOneOrMoreSingleOrDoubleBonds(atomNumber)) {
	        IAtomType type = getAtomType("C.sp2");
	        if (isAcceptable(atomNumber, type)) return type;
	    } else if (isCharged) {
	    	int charge = atom.getFormalCharge().intValue();
	        if (charge == 1) {
	            if (getConnectedAtomsCount(atomNumber) == 0) {
	                IAtomType type = getAtomType("C.plus.sp2");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else {
	                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	                if (maxBondOrder == CDKConstants.BONDORDER_TRIPLE) {
	                    IAtomType type = getAtomType("C.plus.sp1");
	                    if (isAcceptable(atomNumber, type)) return type;
	                } else if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
	                    IAtomType type = getAtomType("C.plus.sp2");
	                    if (isAcceptable(atomNumber, type)) return type;
	                } else if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
	                    IAtomType type = getAtomType("C.plus.planar");
	                    if (isAcceptable(atomNumber, type)) return type;
	                } 
	            }
	        } else if (charge == -1) {
	            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	            if (maxBondOrder == CDKConstants.BONDORDER_SINGLE &&
	                    getConnectedAtomsCount(atomNumber) <= 3) {
	                if (isRingAtom(atomNumber) && bothNeighborsAreSp2(atomNumber)) {
	                    IAtomType type = getAtomType("C.minus.planar");
	                    if (isAcceptable(atomNumber, type)) return type;
	                }
	                IAtomType type = getAtomType("C.minus.sp3");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE &&
	                    getConnectedAtomsCount(atomNumber) <= 3) {
	                IAtomType type = getAtomType("C.minus.sp2");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else if (maxBondOrder == CDKConstants.BONDORDER_TRIPLE &&
	                    getConnectedAtomsCount(atomNumber) <= 1) {
	                IAtomType type = getAtomType("C.minus.sp1");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	        }
	        return null;
	    } else if (getConnectedAtomsCount(atomNumber) > 4) {
	        // FIXME: I don't perceive carbons with more than 4 connections yet
	        return null;
	    } else { // OK, use bond order info
	        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	        if (maxBondOrder == IBond.Order.QUADRUPLE) {
	            // WTF??
	            return null;
	        } else if (maxBondOrder == CDKConstants.BONDORDER_TRIPLE) {
	            IAtomType type = getAtomType("C.sp");
	            if (isAcceptable(atomNumber, type)) return type;
	        } else if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
	            // OK, one or two double bonds?
	            int doubleBondCount = countAttachedDoubleBonds(atomNumber);
	            if (doubleBondCount == 2) {
	                IAtomType type = getAtomType("C.allene");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else if (doubleBondCount == 1) {
	                IAtomType type = getAtomType("C.sp2");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	        } else {
	            if (hasAromaticBond(atomNumber)) {
	                IAtomType type = getAtomType("C.sp2");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	            IAtomType type = getAtomType("C.sp3");
	            if (isAcceptable(atomNumber, type)) return type;
	        }
	    }
    	return null;
    }

	private void cacheProperties() {
		int atomCount = this.atomContainer.getAtomCount();
		maxBondOrders = new IBond.Order[atomCount]; // all null, which is always lowest
		neighborCounts = new int[atomCount]; // all 0, which sounds like a good default
		piBondCounts = new int[atomCount];
		singleBondCounts = new int[atomCount];

		for (IBond bond : atomContainer.bonds()) {
			IBond.Order order = bond.getOrder();
			for (IAtom bondAtom : bond.atoms()) {
				int bondAtomCount = atomContainer.getAtomNumber(bondAtom);
				neighborCounts[bondAtomCount]++; // assume no bond double, etc
				if (order == IBond.Order.DOUBLE) {
					piBondCounts[bondAtomCount]++;
				} else if (order == IBond.Order.SINGLE) {
					singleBondCounts[bondAtomCount]++;
				}
				IBond.Order prevOrder = maxBondOrders[bondAtomCount];
				if (prevOrder == null || BondManipulator.isHigherOrder(order, prevOrder))
					maxBondOrders[bondAtomCount] = order;
			}
		}
	}
	
    private int getConnectedAtomsCount(int atomNumber) {
		return neighborCounts[atomNumber];
	}

    private Order getMaximumBondOrder(int atomNumber) {
		return maxBondOrders[atomNumber];
	}

	private boolean hasOneOrMoreSingleOrDoubleBonds(int atomNumber) {
    	for (IBond bond : atomContainer.getConnectedBondsList(atomContainer.getAtom(atomNumber))) {
    		if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE)) return true;
    	}
		return false;
	}

	private boolean hasOneSingleElectron(int atomNumber) {
	    Iterator<ISingleElectron> singleElectrons = atomContainer.singleElectrons().iterator();
	    while (singleElectrons.hasNext()) {
	    	if (singleElectrons.next().contains(atomContainer.getAtom(atomNumber))) return true;
	    }
	    return false;
    }

    private int countSingleElectrons(int atomNumber) {
	    Iterator<ISingleElectron> singleElectrons = atomContainer.singleElectrons().iterator();
	    int count = 0;
	    while (singleElectrons.hasNext()) {
	    	if (singleElectrons.next().contains(atomContainer.getAtom(atomNumber))) count++;
	    }
	    return count;
    }

    private IAtomType perceiveOxygenRadicals(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() == 0) {
            if (getConnectedAtomsCount(atomNumber) <= 1) {
                IAtomType type = getAtomType("O.sp3.radical");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (atom.getFormalCharge() == +1) {
            if (getConnectedAtomsCount(atomNumber) == 0) {
                IAtomType type = getAtomType("O.plus.radical");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (getConnectedAtomsCount(atomNumber) <= 2) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == IBond.Order.SINGLE) {
                    IAtomType type = getAtomType("O.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (maxBondOrder == IBond.Order.DOUBLE) {
                    IAtomType type = getAtomType("O.plus.sp2.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            }
        }
        return null;
    }
    
    private boolean isCharged(IAtom atom) {
        return (atom.getFormalCharge() != CDKConstants.UNSET && atom.getFormalCharge() != 0);
    }
    
    private boolean hasHybridization(IAtom atom) {
        return atom.getHybridization() != CDKConstants.UNSET;
    }
    
	private IAtomType perceiveOxygens(int atomNumber) throws CDKException {
	    if (hasOneSingleElectron(atomNumber)) {
	        return perceiveOxygenRadicals(atomNumber);
	    }
    	IAtom atom = atomContainer.getAtom(atomNumber);
	    
	    // if hybridization is given, use that
	    if (hasHybridization(atom) && !isCharged(atom)) {
	        if (atom.getHybridization() == Hybridization.SP2) {
	            int connectedAtomsCount = getConnectedAtomsCount(atomNumber);
	            if (connectedAtomsCount == 1) {
	                if (isCarboxylate(atomNumber)) {
	                    IAtomType type = getAtomType("O.sp2.co2");
	                    if (isAcceptable(atomNumber, type)) return type;    				        
	                } else {
	                    IAtomType type = getAtomType("O.sp2");
	                    if (isAcceptable(atomNumber, type)) return type;
	                }
	            } else if (connectedAtomsCount == 2) {
	                IAtomType type = getAtomType("O.planar3");
	                if (isAcceptable(atomNumber, type)) return type;
	            }    				
	        } else if (atom.getHybridization() == Hybridization.SP3) {
	            IAtomType type = getAtomType("O.sp3");
	            if (isAcceptable(atomNumber, type)) return type;
	        } else if (atom.getHybridization() == Hybridization.PLANAR3) {
	            IAtomType type = getAtomType("O.planar3");
	            if (isAcceptable(atomNumber, type)) return type;
	        }
	    } else if (isCharged(atom)) {
	        if (atom.getFormalCharge() == -1 &&
	                getConnectedAtomsCount(atomNumber) <= 1) {
	            if (isCarboxylate(atomNumber)) {
	                IAtomType type = getAtomType("O.minus.co2");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else {
	                IAtomType type = getAtomType("O.minus");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	        } else if (atom.getFormalCharge() == -2 &&
	                getConnectedAtomsCount(atomNumber) == 0) {
	            IAtomType type = getAtomType("O.minus2");
	            if (isAcceptable(atomNumber, type)) return type;
	        } else if (atom.getFormalCharge() == +1) {
	            if (getConnectedAtomsCount(atomNumber) == 0) {
	                IAtomType type = getAtomType("O.plus");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	            if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
	                IAtomType type = getAtomType("O.plus.sp2");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else if (maxBondOrder == CDKConstants.BONDORDER_TRIPLE) {
	                IAtomType type = getAtomType("O.plus.sp1");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else {
	                IAtomType type = getAtomType("O.plus");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	        }
	        return null;
	    } else if (getConnectedAtomsCount(atomNumber) > 2) {
	        // FIXME: I don't perceive carbons with more than 4 connections yet
	        return null;
	    } else if (getConnectedAtomsCount(atomNumber) == 0) {
	        IAtomType type = getAtomType("O.sp3");
	        if (isAcceptable(atomNumber, type)) return type;
	    } else { // OK, use bond order info
	        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
	        if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
	            if (isCarboxylate(atomNumber)) {
	                IAtomType type = getAtomType("O.sp2.co2");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else {
	                IAtomType type = getAtomType("O.sp2");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	        } else if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
	            int explicitHydrogens = countExplicitHydrogens(atomNumber);
	            int connectedHeavyAtoms = getConnectedAtomsCount(atomNumber) - explicitHydrogens; 
	            if (connectedHeavyAtoms == 2) {
	                // a O.sp3 which is expected to take part in an aromatic system
	                if (isRingAtom(atomNumber) && bothNeighborsAreSp2(atomNumber)) {
	                    IAtomType type = getAtomType("O.planar3");
	                    if (isAcceptable(atomNumber, type)) return type;
	                }
	                IAtomType type = getAtomType("O.sp3");
	                if (isAcceptable(atomNumber, type)) return type;
	            } else {
	                IAtomType type = getAtomType("O.sp3");
	                if (isAcceptable(atomNumber, type)) return type;
	            }
	        }
	    }
    	return null;
    }

    private boolean isCarboxylate(int atomNumber) {
        // assumes that the oxygen only has one neighbor (C=O, or C-[O-])
        List<IAtom> neighbors = atomContainer.getConnectedAtomsList(atomContainer.getAtom(atomNumber));
        if (neighbors.size() != 1) return false;
        IAtom carbon = neighbors.get(0);
        if (!"C".equals(carbon.getSymbol())) return false;
        
        int oxygenCount = 0;
        int singleBondedNegativeOxygenCount = 0;
        int doubleBondedOxygenCount = 0;
        for (IBond cBond : atomContainer.getConnectedBondsList(carbon)) {
            IAtom neighbor = cBond.getConnectedAtom(carbon);
            if ("O".equals(neighbor.getSymbol())) {
                oxygenCount++;
                IBond.Order order = cBond.getOrder();
                Integer charge = neighbor.getFormalCharge();
                if (order == IBond.Order.SINGLE && charge != null && charge == -1) {
                    singleBondedNegativeOxygenCount++;
                } else if (order == IBond.Order.DOUBLE) {
                    doubleBondedOxygenCount++;
                }
            }
        }
        return (oxygenCount == 2) && (singleBondedNegativeOxygenCount == 1) && (doubleBondedOxygenCount == 1);
    }

    private boolean atLeastTwoNeighborsAreSp2(int atomNumber) {
    	int count = 0;
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	Iterator<IAtom> atoms = atomContainer.getConnectedAtomsList(atom).iterator();
    	while (atoms.hasNext() && (count < 2)) {
    		IAtom nextAtom = atoms.next();
    		if (!nextAtom.getSymbol().equals("H")) {
    			if (nextAtom.getHybridization() != CDKConstants.UNSET &&
    				nextAtom.getHybridization() == Hybridization.SP2) {
    				// OK, it's SP2
    				count++;
    			} else if (piBondCounts[atomContainer.getAtomNumber(nextAtom)] > 0) {
    				// OK, it's SP2
    				count++;
    			} else if (atomContainer.getBond(atom, nextAtom).getFlag(CDKConstants.ISAROMATIC)) {
                    // two aromatic bonds indicate sp2
                    count++;
                } // OK, not SP2
    		}
    	}
    	return count >= 2;
    }

    private boolean bothNeighborsAreSp2(int atomNumber) {       
    	return atLeastTwoNeighborsAreSp2(atomNumber);
    }

    private IAtomType perceiveNitrogenRadicals(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (getConnectedAtomsCount(atomNumber) >= 1 &&
                getConnectedAtomsCount(atomNumber) <= 2) {
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (atom.getFormalCharge() != CDKConstants.UNSET &&
                    atom.getFormalCharge() == +1) {
                if (maxBondOrder == IBond.Order.DOUBLE) {
                    IAtomType type = getAtomType("N.plus.sp2.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (maxBondOrder == IBond.Order.SINGLE) {
                    IAtomType type = getAtomType("N.plus.sp3.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (atom.getFormalCharge() == CDKConstants.UNSET ||
                    atom.getFormalCharge() == 0) {
                if (maxBondOrder == IBond.Order.SINGLE) {
                    IAtomType type = getAtomType("N.sp3.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (maxBondOrder == IBond.Order.DOUBLE) {
                    IAtomType type = getAtomType("N.sp2.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            }
        } else {
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (atom.getFormalCharge() != CDKConstants.UNSET &&
                    atom.getFormalCharge() == +1 && maxBondOrder == IBond.Order.SINGLE) {
                IAtomType type = getAtomType("N.plus.sp3.radical");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }
    private IAtomType perceiveMolybdenum(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 4) {
                IAtomType type = getAtomType("Mo.4");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
            IAtomType type1 = getAtomType("Mo.metallic");
            if (isAcceptable(atomNumber, type1)) {
                return type1;
            }
        }
        return null;
    }
    private IAtomType perceiveNitrogens(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
       // if hybridization is given, use that
        if (hasOneSingleElectron(atomNumber)) {
            return perceiveNitrogenRadicals(atomNumber);
        } else if (hasHybridization(atom) && !isCharged(atom)) {
            if (atom.getHybridization() == Hybridization.SP1) {
                int neighborCount = getConnectedAtomsCount(atomNumber);
                if (neighborCount > 1) {
                    IAtomType type = getAtomType("N.sp1.2");
                    if (isAcceptable(atomNumber, type)) return type;
                } else {
                    IAtomType type = getAtomType("N.sp1");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (atom.getHybridization() == Hybridization.SP2) {
            	if (isAmide(atomNumber)) {
                    IAtomType type = getAtomType("N.amide");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (isThioAmide(atomNumber)) {
                    IAtomType type = getAtomType("N.thioamide");
                    if (isAcceptable(atomNumber, type)) return type;
                }
                // but an sp2 hyb N might N.sp2 or N.planar3 (pyrrole), so check for the latter
            	int neighborCount = getConnectedAtomsCount(atomNumber);
            	if (neighborCount == 4 &&
            	    IBond.Order.DOUBLE == getMaximumBondOrder(atomNumber)) {
            	    IAtomType type = getAtomType("N.oxide");
                    if (isAcceptable(atomNumber, type)) return type;
            	} else
            	if (neighborCount > 1 && bothNeighborsAreSp2(atomNumber)) {
            		IRing ring = getRing(atomNumber);
            		int ringSize = ring == null ? 0 : ring.getAtomCount();
            		if (ring != null && ringSize > 0) {
            			if (neighborCount == 3) {
                            IBond.Order maxOrder = getMaximumBondOrder(atomNumber);
                            if (maxOrder == IBond.Order.DOUBLE) {
                                IAtomType type = getAtomType("N.sp2.3");
                                if (isAcceptable(atomNumber, type)) return type;
                            } else if (maxOrder == IBond.Order.SINGLE) {
                                IAtomType type = getAtomType("N.planar3");
                                if (isAcceptable(atomNumber, type)) return type;
                            }
            			} else if (neighborCount == 2) {
            				IBond.Order maxOrder = getMaximumBondOrder(atomNumber);
            				if (maxOrder == IBond.Order.SINGLE) {
            				    if (atom.getImplicitHydrogenCount() != CDKConstants.UNSET && atom.getImplicitHydrogenCount() == 1) {
            						IAtomType type = getAtomType("N.planar3");
            						if (isAcceptable(atomNumber, type)) return type;
            					} else {
            						IAtomType type = getAtomType("N.sp2");
            						if (isAcceptable(atomNumber, type)) return type;
            					}
            				} else if (maxOrder == IBond.Order.DOUBLE) {
            					IAtomType type = getAtomType("N.sp2");
            					if (isAcceptable(atomNumber, type)) return type;
            				}
            			}
            		}
            	}
                IAtomType type = getAtomType("N.sp2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getHybridization() == Hybridization.SP3) {
                IAtomType type = getAtomType("N.sp3");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getHybridization() == Hybridization.PLANAR3) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (getConnectedAtomsCount(atomNumber) == 3 &&
                        maxBondOrder == CDKConstants.BONDORDER_DOUBLE &&
                        countAttachedDoubleBonds(atomNumber, "O") == 2) {
                    IAtomType type = getAtomType("N.nitro");
                    if (isAcceptable(atomNumber, type)) return type;
                }
                IAtomType type = getAtomType("N.planar3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (isCharged(atom)) {
            if (atom.getFormalCharge() == 1) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == CDKConstants.BONDORDER_SINGLE ||
                        getConnectedAtomsCount(atomNumber) == 0) {
                    if (atom.getHybridization() == IAtomType.Hybridization.SP2) {
                        IAtomType type = getAtomType("N.plus.sp2");
                        if (isAcceptable(atomNumber, type)) return type;
                    }
                    IAtomType type = getAtomType("N.plus");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
                    int doubleBonds= countAttachedDoubleBonds(atomNumber);
                    if (doubleBonds == 1) {
                        IAtomType type = getAtomType("N.plus.sp2");
                        if (isAcceptable(atomNumber, type)) return type;
                    } else if (doubleBonds == 2) {
                        IAtomType type = getAtomType("N.plus.sp1");
                        if (isAcceptable(atomNumber, type)) return type;
                    }
                } else if (maxBondOrder == CDKConstants.BONDORDER_TRIPLE) {
                    if (getConnectedAtomsCount(atomNumber) == 2) {
                        IAtomType type = getAtomType("N.plus.sp1");
                        if (isAcceptable(atomNumber, type)) return type;
                    }
                }
            } else if (atom.getFormalCharge() == -1) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
                    if (getConnectedAtomsCount(atomNumber) >= 2 &&
                    		bothNeighborsAreSp2(atomNumber) &&
                    		isRingAtom(atomNumber)) {
                        IAtomType type = getAtomType("N.minus.planar3");
                        if (isAcceptable(atomNumber, type)) return type;
                    } else if (getConnectedAtomsCount(atomNumber) <= 2) {
                        IAtomType type = getAtomType("N.minus.sp3");
                        if (isAcceptable(atomNumber, type)) return type;
                    }
                } else if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
                    if (getConnectedAtomsCount(atomNumber) <= 1) {
                        IAtomType type = getAtomType("N.minus.sp2");
                        if (isAcceptable(atomNumber, type)) return type;
                    }
                }
            }
        } else if (getConnectedAtomsCount(atomNumber) > 3) {
            if (getConnectedAtomsCount(atomNumber) == 4 &&
                countAttachedDoubleBonds(atomNumber) == 1) {
                IAtomType type = getAtomType("N.oxide");
                if (isAcceptable(atomNumber, type)) return type;
            }
            return null;
        } else if (getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("N.sp3");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (hasOneOrMoreSingleOrDoubleBonds(atomNumber)) {
        	int connectedAtoms = getConnectedAtomsCount(atomNumber) +
        		(atom.getImplicitHydrogenCount() == CDKConstants.UNSET
        		    ? 0
        			: atom.getImplicitHydrogenCount());
        	if (connectedAtoms == 3) {
            	IAtomType type = getAtomType("N.planar3");
            	if (isAcceptable(atomNumber, type)) return type;
        	}
        	IAtomType type = getAtomType("N.sp2");
        	if (isAcceptable(atomNumber, type)) return type;
        } else { // OK, use bond order info
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
                if (isAmide(atomNumber)) {
                    IAtomType type = getAtomType("N.amide");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (isThioAmide(atomNumber)) {
                    IAtomType type = getAtomType("N.thioamide");
                    if (isAcceptable(atomNumber, type)) return type;
                }
                int explicitHydrogens = countExplicitHydrogens(atomNumber);
                int connectedHeavyAtoms = getConnectedAtomsCount(atomNumber) - explicitHydrogens;
                if (connectedHeavyAtoms == 2) {
                	List<IBond> bonds = atomContainer.getConnectedBondsList(atom);
                    if (bonds.get(0).getFlag(CDKConstants.ISAROMATIC) &&
                            bonds.get(1).getFlag(CDKConstants.ISAROMATIC)) {
                        Integer hCount = atom.getImplicitHydrogenCount();
                        if (hCount == CDKConstants.UNSET || hCount == 0) {
                            if (getMaximumBondOrder(atomNumber) == CDKConstants.BONDORDER_SINGLE &&
                                    isSingleHeteroAtom(atomNumber)) {
                                IAtomType type = getAtomType("N.planar3");
                                if (isAcceptable(atomNumber, type)) return type;
                            } else {
                                IAtomType type = getAtomType("N.sp2");
                                if (isAcceptable(atomNumber, type)) return type;
                            }
                        } else if (hCount == 1) {
                            IAtomType type = getAtomType("N.planar3");
                            if (isAcceptable(atomNumber, type))
                                return type;
                        }
                	} else if (bothNeighborsAreSp2(atomNumber) && isRingAtom(atomNumber)) {
                		// a N.sp3 which is expected to take part in an aromatic system
                		IAtomType type = getAtomType("N.planar3");
                		if (isAcceptable(atomNumber, type)) return type;
                	} else {
                		IAtomType type = getAtomType("N.sp3");
                		if (isAcceptable(atomNumber, type)) return type;
                	}
                } else if (connectedHeavyAtoms == 3) {
                	if (bothNeighborsAreSp2(atomNumber) && isRingAtom(atomNumber)) {
                		IAtomType type = getAtomType("N.planar3");
                		if (isAcceptable(atomNumber, type)) return type;
                	}
                	IAtomType type = getAtomType("N.sp3");
                	if (isAcceptable(atomNumber, type)) return type;
                } else if (connectedHeavyAtoms == 1) {
                    IAtomType type = getAtomType("N.sp3");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (connectedHeavyAtoms == 0) {
                    IAtomType type = getAtomType("N.sp3");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
                if (getConnectedAtomsCount(atomNumber) == 3 &&
                        countAttachedDoubleBonds(atomNumber, "O") == 2) {
                    IAtomType type = getAtomType("N.nitro");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (getConnectedAtomsCount(atomNumber) == 3 &&
                        countAttachedDoubleBonds(atomNumber) > 0) {
                    IAtomType type = getAtomType("N.sp2.3");
                    if (isAcceptable(atomNumber, type)) return type;
                }
                IAtomType type = getAtomType("N.sp2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (maxBondOrder == CDKConstants.BONDORDER_TRIPLE) {
                int neighborCount = getConnectedAtomsCount(atomNumber);
                if (neighborCount > 1) {
                    IAtomType type = getAtomType("N.sp1.2");
                    if (isAcceptable(atomNumber, type)) return type;
                } else {
                    IAtomType type = getAtomType("N.sp1");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            }
        }
    	return null;
    }

    /**
     * Determines whether the bonds (up to two spheres away) are only to non
     * hetroatoms. Currently used in N.planar3 perception of (e.g. pyrrole).
     *
     * @param atom an atom to test
     * @param container container of the atom
     *
     * @return whether the atom's only bonds are to hetroatoms
     * @see #perceiveNitrogens(org.openscience.cdk.interfaces.IAtomContainer, org.openscience.cdk.interfaces.IAtom)
     */
    private boolean isSingleHeteroAtom(int atomNumber) {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        List<IAtom> connected = atomContainer.getConnectedAtomsList(atom);

        for (IAtom atom1 : connected) {

            boolean aromatic = atomContainer.getBond(atom, atom1).getFlag(CDKConstants.ISAROMATIC);

            // ignoring non-aromatic bonds
            if(!aromatic)
                continue;

            // found a hetroatom - we're not a single hetroatom
            if(!"C".equals(atom1.getSymbol()))
                return false;

            // check the second sphere
            for (IAtom atom2 : atomContainer.getConnectedAtomsList(atom1)) {

                if (atom2 != atom
                   && atomContainer.getBond(atom1, atom2).getFlag(CDKConstants.ISAROMATIC)
                   && !"C".equals(atom2.getSymbol())) {
                        return false;
                }

            }

        }

        return true;

    }

    private void cacheRingProperties() {
    	int atomCount = atomContainer.getAtomCount();
    	st = new SpanningTree(atomContainer);
    	isRingAtom = new boolean[atomCount];
    	if (st.getBondsCyclicCount() != 0)
    		for (int i=0; i<atomCount; i++) {
    			isRingAtom[i] = st.getCyclicFragmentsContainer().contains(atomContainer.getAtom(i));
    		}
    }
    
    private boolean isRingAtom(int atomNumber) {
    	if (st == null) cacheRingProperties();
        return isRingAtom[atomNumber];
    }

    private IRing getRing(int atomNumber) {
    	if (st == null) cacheRingProperties();
    	if (!isRingAtom[atomNumber]) return null;
    	try {
    		IRingSet set = st.getAllRings();
    		IAtom atom = atomContainer.getAtom(atomNumber);
    		for (int i=0; i<set.getAtomContainerCount(); i++) {
    			IRing ring = (IRing)set.getAtomContainer(i);
    			if (ring.contains(atom)) {
    				return ring;
    			}
    		}
    	} catch (NoSuchAtomException exception) {
    		return null;
    	}
    	return null;
    }

    private boolean isAmide(int atomNumber) {
    	List<IAtom> neighbors = atomContainer.getConnectedAtomsList(atomContainer.getAtom(atomNumber));
    	for (IAtom neighbor : neighbors) {
    		if (neighbor.getSymbol().equals("C")) {
    			if (countAttachedDoubleBonds(atomNumber, "O") == 1) return true;
    		}
    	}
    	return false;
    }

    private boolean isThioAmide(int atomNumber) {
        List<IAtom> neighbors = atomContainer.getConnectedAtomsList(atomContainer.getAtom(atomNumber));
        for (IAtom neighbor : neighbors) {
            if (neighbor.getSymbol().equals("C")) {
                if (countAttachedDoubleBonds(atomNumber, "S") == 1) return true;
            }
        }
        return false;
    }

    private int countExplicitHydrogens(int atomNumber) {
    	int count = 0;
        for (IAtom aAtom : atomContainer.getConnectedAtomsList(atomContainer.getAtom(atomNumber))) {
            if (aAtom.getSymbol().equals("H")) {
                count++;
            }
        }
    	return count;
    }
    
    private IAtomType perceiveIron(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if ("Fe".equals(atom.getSymbol())) {
            if (hasOneSingleElectron(atomNumber)) {
                // no idea how to deal with this yet
                return null;
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == 0)) {
                IAtomType type = getAtomType("Fe.metallic");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
                int neighbors = getConnectedAtomsCount(atomNumber);
                if (neighbors == 2) {
                    IAtomType type5 = getAtomType("Fe.2");
                    if (isAcceptable(atomNumber, type5)) {
                        return type5;
                    }
                } else if (neighbors == 3) {
                    IAtomType type6 = getAtomType("Fe.3");
                    if (isAcceptable(atomNumber, type6)) {
                        return type6;
                    }
                } else if (neighbors == 4) {
                    IAtomType type7 = getAtomType("Fe.4");
                    if (isAcceptable(atomNumber, type7)) {
                        return type7;
                    }
                } else if (neighbors == 5) {
                    IAtomType type8 = getAtomType("Fe.5");
                    if (isAcceptable(atomNumber, type8)) {
                        return type8;
                    }
                } else if (neighbors == 6) {
                    IAtomType type9 = getAtomType("Fe.6");
                    if (isAcceptable(atomNumber, type9)) {
                        return type9;
                    }
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == 2)) {
                int neighbors = getConnectedAtomsCount(atomNumber);
                if (neighbors <= 1) {
                    IAtomType type = getAtomType("Fe.2plus");
                    if (isAcceptable(atomNumber, type)) {
                        return type;
                    }
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == 1)) {
                int neighbors = getConnectedAtomsCount(atomNumber);

                if (neighbors == 2) {
                    IAtomType type0 = getAtomType("Fe.plus");
                    if (isAcceptable(atomNumber, type0)) {
                        return type0;
                    }
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == 3)) {
                IAtomType type1 = getAtomType("Fe.3plus");
                if (isAcceptable(atomNumber, type1)) {
                    return type1;
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == -2)) {
                IAtomType type2 = getAtomType("Fe.2minus");
                if (isAcceptable(atomNumber, type2)) {
                    return type2;
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == -3)) {
                IAtomType type3 = getAtomType("Fe.3minus");
                if (isAcceptable(atomNumber, type3)) {
                    return type3;
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == -4)) {
                IAtomType type4 = getAtomType("Fe.4minus");
                if (isAcceptable(atomNumber, type4)) {
                    return type4;
                }
            }
        }
        return null;
    }


    private IAtomType perceiveMercury(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if ("Hg".equals(atom.getSymbol())) {
            if (hasOneSingleElectron(atomNumber)) {
                // no idea how to deal with this yet
                return null;
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == -1)) {
                IAtomType type = getAtomType("Hg.minus");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == 2)) {
                IAtomType type = getAtomType("Hg.2plus");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == +1)) {
                int neighbors = getConnectedAtomsCount(atomNumber);
                if (neighbors <= 1) {  
                    IAtomType type = getAtomType("Hg.plus");
                    if (isAcceptable(atomNumber, type)) {
                        return type;
                    }
                }
            } else if ((atom.getFormalCharge() != null
                    && atom.getFormalCharge() == 0)) {
                int neighbors = getConnectedAtomsCount(atomNumber);
                if (neighbors == 2) {
                    IAtomType type = getAtomType("Hg.2");
                    if (isAcceptable(atomNumber, type)) {
                        return type;
                    }
                } else if (neighbors == 1) {
                    IAtomType type = getAtomType("Hg.1");
                    if (isAcceptable(atomNumber, type)) {
                        return type;
                    }
                } else if (neighbors == 0) {
                    IAtomType type = getAtomType("Hg.metallic");
                    if (isAcceptable(atomNumber, type)) {
                        return type;
                    }
                }
            }
        }
        return null;
    }

    private IAtomType perceiveSulphurs(int atomNumber)
    throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        List<IBond> neighbors = atomContainer.getConnectedBondsList(atom);
        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
        int neighborcount = neighbors.size();
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if (atom.getHybridization() != CDKConstants.UNSET &&
                   atom.getHybridization() == Hybridization.SP2 &&
                   atom.getFormalCharge() != CDKConstants.UNSET &&
                   atom.getFormalCharge() == +1) {
            if (neighborcount == 3) {
                IAtomType type = getAtomType("S.inyl.charged");
                if (isAcceptable(atomNumber, type)) return type;
            } else {
                IAtomType type = getAtomType("S.plus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET &&
                   atom.getFormalCharge() != 0) {

            if (atom.getFormalCharge() == -1 &&
                    neighborcount == 1) {
                IAtomType type = getAtomType("S.minus");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() == +1 &&
                    neighborcount == 2) {
                IAtomType type = getAtomType("S.plus");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() == +1 &&
                    neighborcount == 3) {
                IAtomType type = getAtomType("S.inyl.charged");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() == +2 &&
                    neighborcount == 4) {
                IAtomType type = getAtomType("S.onyl.charged");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() == -2 &&
                    neighborcount == 0) {
                IAtomType type = getAtomType("S.2minus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 0) {
            if (atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == 0) {
                IAtomType type = getAtomType("S.3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 1) {
            if (atomContainer.getConnectedBondsList(atom).get(0).getOrder() == CDKConstants.BONDORDER_DOUBLE) {
                IAtomType type = getAtomType("S.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atomContainer.getConnectedBondsList(atom).get(0).getOrder() == CDKConstants.BONDORDER_SINGLE) {
                IAtomType type = getAtomType("S.3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 2) {
            if (isRingAtom(atomNumber) && bothNeighborsAreSp2(atomNumber)) {
                if (countAttachedDoubleBonds(atomNumber) == 2) {
                    IAtomType type = getAtomType("S.inyl.2");
                    if (isAcceptable(atomNumber, type)) return type;
                } else {
                    IAtomType type = getAtomType("S.planar3");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (countAttachedDoubleBonds(atomNumber, "O") == 2) {
                IAtomType type = getAtomType("S.oxide");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (countAttachedDoubleBonds(atomNumber) == 2) {
                IAtomType type = getAtomType("S.inyl.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (countAttachedDoubleBonds(atomNumber) <= 1) {
                IAtomType type = getAtomType("S.3");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (countAttachedDoubleBonds(atomNumber) == 0
                    && countAttachedSingleBonds(atomNumber) == 2) {
                IAtomType type = getAtomType("S.octahedral");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 3) {
            int doubleBondedAtoms = countAttachedDoubleBonds(atomNumber);
            if (doubleBondedAtoms == 1) {
                IAtomType type = getAtomType("S.inyl");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (doubleBondedAtoms == 3) {
                IAtomType type = getAtomType("S.trioxide");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (doubleBondedAtoms == 0) {
                IAtomType type = getAtomType("S.anyl");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 4) {
            // count the number of double bonded oxygens
            int doubleBondedOxygens = countAttachedDoubleBonds(atomNumber, "O");
            int doubleBondedNitrogens = countAttachedDoubleBonds(atomNumber, "N");
            int doubleBondedSulphurs = countAttachedDoubleBonds(atomNumber, "S");
            int countAttachedDoubleBonds = countAttachedDoubleBonds(atomNumber);

            if (doubleBondedOxygens + doubleBondedNitrogens == 2) {
                IAtomType type = getAtomType("S.onyl");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (doubleBondedSulphurs == 1
                    && doubleBondedOxygens == 1) {
                IAtomType type = getAtomType("S.thionyl");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
                IAtomType type = getAtomType("S.anyl");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (doubleBondedOxygens == 1) {
                IAtomType type = getAtomType("S.sp3d1");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (countAttachedDoubleBonds == 2
                    && maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
                IAtomType type = getAtomType("S.sp3.4");
                if (isAcceptable(atomNumber, type)) return type;
            }

        } else if (neighborcount == 5) {

            if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {

                IAtomType type = getAtomType("S.sp3d1");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
                IAtomType type = getAtomType("S.octahedral");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 6) {
            if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
                IAtomType type = getAtomType("S.octahedral");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }

    private IAtomType perceivePhosphors(int atomNumber)
    throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        List<IBond> neighbors = atomContainer.getConnectedBondsList(atom);
        int neighborcount = neighbors.size();
        IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
        if (countSingleElectrons(atomNumber) == 3) {
        	IAtomType type = getAtomType("P.se.3");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if (neighborcount == 0) {
            if (atom.getFormalCharge() == null ||
                atom.getFormalCharge().intValue() == 0) {
                IAtomType type = getAtomType("P.ine");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 1) {
            if (atom.getFormalCharge() == null ||
                atom.getFormalCharge().intValue() == 0) {
                IAtomType type = getAtomType("P.ide");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 3) {
        	int doubleBonds = countAttachedDoubleBonds(atomNumber);
            if (atom.getFormalCharge() != null &&
                atom.getFormalCharge().intValue() == 1) {
                IAtomType type = getAtomType("P.anium");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (doubleBonds == 1) {
            	IAtomType type = getAtomType("P.ate");
                if (isAcceptable(atomNumber, type)) return type;
            } else {
                IAtomType type = getAtomType("P.ine");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 2) {
            if (maxBondOrder == CDKConstants.BONDORDER_DOUBLE) {
                if (atom.getFormalCharge() != null &&
                    atom.getFormalCharge().intValue() == 1) {
                    IAtomType type = getAtomType("P.sp1.plus");
                    if (isAcceptable(atomNumber, type)) return type;
                } else {
                    IAtomType type = getAtomType("P.irane");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (maxBondOrder == CDKConstants.BONDORDER_SINGLE) {
                IAtomType type = getAtomType("P.ine");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 4) {
            // count the number of double bonded oxygens
            int doubleBonds = countAttachedDoubleBonds(atomNumber);
            if (atom.getFormalCharge() == 1 && doubleBonds == 0) {
                IAtomType type = getAtomType("P.ate.charged");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (doubleBonds == 1){
                IAtomType type = getAtomType("P.ate");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 5) {
            if (atom.getFormalCharge() == null ||
                atom.getFormalCharge().intValue() == 0) {
                IAtomType type = getAtomType("P.ane");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
    	return null;
    }
    
    private IAtomType perceiveHydrogens(int atomNumber)
    throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        int neighborcount = getConnectedAtomsCount(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            if ((atom.getFormalCharge() == CDKConstants.UNSET || atom.getFormalCharge() == 0) &&
                    neighborcount == 0) {
                IAtomType type = getAtomType("H.radical");
                if (isAcceptable(atomNumber, type)) return type;
            }
            return null;
        } else if (neighborcount == 2) {
            // FIXME: bridging hydrogen as in B2H6
            return null;
        } else if (neighborcount == 1) {
            if (atom.getFormalCharge() == CDKConstants.UNSET || atom.getFormalCharge() == 0) {
                IAtomType type = getAtomType("H");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 0) {
            if (atom.getFormalCharge() == CDKConstants.UNSET || atom.getFormalCharge() == 0) {
                IAtomType type = getAtomType("H");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() == 1){
                IAtomType type = getAtomType("H.plus");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() == -1){
                IAtomType type = getAtomType("H.minus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
    	return null;
    }

    private IAtomType perceiveLithium(int atomNumber)
    	throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        int neighborcount = getConnectedAtomsCount(atomNumber);
        if (neighborcount == 1) {
            if (atom.getFormalCharge() == CDKConstants.UNSET ||
                    atom.getFormalCharge() == 0) {
                IAtomType type = getAtomType("Li");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (neighborcount == 0) {
            if (atom.getFormalCharge() == CDKConstants.UNSET
                    || atom.getFormalCharge() == 0) {
                IAtomType type = getAtomType("Li.neutral");
                if (isAcceptable(atomNumber, type)) return type;
            }
            if (atom.getFormalCharge() == CDKConstants.UNSET
                    || atom.getFormalCharge() == +1) {
                IAtomType type = getAtomType("Li.plus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
    	return null;
    }

    private IAtomType perceiveFluors(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		if (getConnectedAtomsCount(atomNumber) == 0) {
    			if (atom.getFormalCharge() != CDKConstants.UNSET &&
    					atom.getFormalCharge() == +1) {
    				IAtomType type = getAtomType("F.plus.radical");
    				if (isAcceptable(atomNumber, type)) return type;
    			} else if (atom.getFormalCharge() == CDKConstants.UNSET ||
    					atom.getFormalCharge() == 0) {
    				IAtomType type = getAtomType("F.radical");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		} else if (getConnectedAtomsCount(atomNumber) <= 1) {
    			IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
    			if (maxBondOrder == IBond.Order.SINGLE) {
    				IAtomType type = getAtomType("F.plus.radical");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		}
    		return null;
    	} else if (atom.getFormalCharge() != CDKConstants.UNSET &&
    			atom.getFormalCharge() != 0) {
    		if (atom.getFormalCharge() == -1) {
    			IAtomType type = getAtomType("F.minus");
    			if (isAcceptable(atomNumber, type)) return type;
    		} else if (atom.getFormalCharge() == 1) {
    			IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
    			if (maxBondOrder == IBond.Order.DOUBLE) {
    				IAtomType type = getAtomType("F.plus.sp2");
    				if (isAcceptable(atomNumber, type)) return type;
    			}else if (maxBondOrder == IBond.Order.SINGLE){
    				IAtomType type = getAtomType("F.plus.sp3");
    				if (isAcceptable(atomNumber, type)) return type;
    			}
    		}
    	} else if (getConnectedAtomsCount(atomNumber) == 1 ||
    			getConnectedAtomsCount(atomNumber) == 0) {
    		IAtomType type = getAtomType("F");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }

    private IAtomType perceiveArsenic(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +1
                && getConnectedAtomsCount(atomNumber) <= 4)) {
            IAtomType type = getAtomType("As.plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 4) {
                IAtomType type = getAtomType("As.5");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
            if (neighbors == 2) {
                IAtomType type = getAtomType("As.2");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
            IAtomType type = getAtomType("As");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +3)) {
            IAtomType type = getAtomType("As.3plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -1)) {
            IAtomType type = getAtomType("As.minus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        }
        return null;
    }   
     
    private IAtomType perceiveThorium(int atomNumber)
            throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if ("Th".equals(atom.getSymbol())) {
            if (atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 0) {
                IAtomType type = getAtomType("Th");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
        }
        return null;
    }

    private IAtomType perceiveRubidium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +1) {
            IAtomType type = getAtomType("Rb.plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0) {
            IAtomType type = getAtomType("Rb.neutral");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        }
        return null;
    }
    private IAtomType perceiveTungstun(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		IAtomType type = getAtomType("W.metallic");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveCopper(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +2)) {
            IAtomType type = getAtomType("Cu.2plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 1) {
                IAtomType type = getAtomType("Cu.1");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            } else {
                IAtomType type01 = getAtomType("Cu.metallic");
                if (isAcceptable(atomNumber, type01)) {
                    return type01;
                }
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +1) {
            IAtomType type02 = getAtomType("Cu.plus");
            if (isAcceptable(atomNumber, type02)) {
                return type02;
            }
        }
        return null;
    }
    private IAtomType perceiveBarium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 2)) {
            IAtomType type = getAtomType("Ba.2plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        }
        return null;
    }
    private IAtomType perceiveAluminium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 3) {
            int connectedBondsCount = getConnectedAtomsCount(atomNumber);
            if (connectedBondsCount == 0) {
                IAtomType type = getAtomType("Al.3plus");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0
                && getConnectedAtomsCount(atomNumber) == 3) {
            IAtomType type = getAtomType("Al");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -3
                && getConnectedAtomsCount(atomNumber) == 6) {
            IAtomType type = getAtomType("Al.3minus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        }
        return null;
    }
    private IAtomType perceiveZinc(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if (getConnectedAtomsCount(atomNumber) == 0
                && (atom.getFormalCharge() != null
                && atom.getFormalCharge() == 0)) {
            IAtomType type = getAtomType("Zn.metallic");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (getConnectedAtomsCount(atomNumber) == 0
                && (atom.getFormalCharge() != null
                && atom.getFormalCharge() == 2)) {
            IAtomType type = getAtomType("Zn.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (getConnectedAtomsCount(atomNumber) == 1
                && (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)) {
            IAtomType type = getAtomType("Zn.1");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (getConnectedAtomsCount(atomNumber) == 2
                && (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)) {
            IAtomType type = getAtomType("Zn");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    private IAtomType perceiveChromium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0
                && getConnectedAtomsCount(atomNumber) == 6) {
            IAtomType type = getAtomType("Cr");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0
                && getConnectedAtomsCount(atomNumber) == 4) {
            IAtomType type = getAtomType("Cr.4");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 6
                && getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Cr.6plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0
                && getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Cr.neutral");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ("Cr".equals(atom.getSymbol())) {
            if (atom.getFormalCharge() != CDKConstants.UNSET
                    && atom.getFormalCharge() == 3
                    && getConnectedAtomsCount(atomNumber) == 0) {
                IAtomType type = getAtomType("Cr.3plus");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
        }
        return null;
    }
    private IAtomType perceivePolodium(int atomNumber) throws CDKException {
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if (getConnectedAtomsCount(atomNumber) == 2) {
    		IAtomType type = getAtomType("Po");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveStannum(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
    			atom.getFormalCharge() == 0 &&
    			getConnectedAtomsCount(atomNumber) <= 4)) {
    		IAtomType type = getAtomType("Sn.sp3");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveScandium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (atom.getFormalCharge() != CDKConstants.UNSET &&
    			atom.getFormalCharge() == -3 &&
    			getConnectedAtomsCount(atomNumber) == 6) {
    		IAtomType type = getAtomType("Sc.3minus");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveNickel(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +2)) {
            IAtomType type = getAtomType("Ni.2plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)
                && getConnectedAtomsCount(atomNumber) == 2) {
            IAtomType type = getAtomType("Ni");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)
                && getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Ni.metallic");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 1)
                && getConnectedAtomsCount(atomNumber) == 1) {
            IAtomType type = getAtomType("Ni.plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        }
        return null;
    }
    private IAtomType perceiveHelium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		IAtomType type = getAtomType("He");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveNeon(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		IAtomType type = getAtomType("Ne");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveArgon(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		IAtomType type = getAtomType("Ar");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveKrypton(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		IAtomType type = getAtomType("Kr");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }
    private IAtomType perceiveXenon(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		if (getConnectedAtomsCount(atomNumber) == 0) {
    			IAtomType type = getAtomType("Xe");
    			if (isAcceptable(atomNumber, type)) return type;
    		} else {
    			IAtomType type = getAtomType("Xe.3");
    			if (isAcceptable(atomNumber, type)) return type;
    		}
    	}
    	return null;
    }
    private IAtomType perceiveRadon(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (hasOneSingleElectron(atomNumber)) {
    		// no idea how to deal with this yet
    		return null;
    	} else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
    			atom.getFormalCharge() == 0)) {
    		IAtomType type = getAtomType("Rn");
    		if (isAcceptable(atomNumber, type)) return type;
    	}
    	return null;
    }

    private IAtomType perceiveSilicon(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0) {
            if (getConnectedAtomsCount(atomNumber) == 2) {
                IAtomType type = getAtomType("Si.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (getConnectedAtomsCount(atomNumber) == 3) {
                IAtomType type = getAtomType("Si.3");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (getConnectedAtomsCount(atomNumber) == 4) {
                IAtomType type = getAtomType("Si.sp3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -2) {
            IAtomType type = getAtomType("Si.2minus.6");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveManganese(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != null
                && atom.getFormalCharge() == 0)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 2) {
                IAtomType type02 = getAtomType("Mn.2");
                if (isAcceptable(atomNumber, type02)) return type02;
            } else if (neighbors == 0) {
                IAtomType type03 = getAtomType("Mn.metallic");
                if (isAcceptable(atomNumber, type03)) return type03;
            }
        } else if ((atom.getFormalCharge() != null
                && atom.getFormalCharge() == +2)) {
            IAtomType type = getAtomType("Mn.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != null
                && atom.getFormalCharge() == +3)) {
            IAtomType type = getAtomType("Mn.3plus");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveSodium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 1)) {
            IAtomType type = getAtomType("Na.plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() == CDKConstants.UNSET
                || atom.getFormalCharge() == 0)
                && getConnectedAtomsCount(atomNumber) == 1) {
            IAtomType type = getAtomType("Na");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)
                && getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Na.neutral");
            if (isAcceptable(atomNumber, type)) return type;
        } 
        return null;
    }
    
    private IAtomType perceiveIodine(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            if (getConnectedAtomsCount(atomNumber) == 0) {
                if (atom.getFormalCharge() != CDKConstants.UNSET &&
                    atom.getFormalCharge() == +1) {
                    IAtomType type = getAtomType("I.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (atom.getFormalCharge() == CDKConstants.UNSET ||
                           atom.getFormalCharge() == 0) {
                    IAtomType type = getAtomType("I.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (getConnectedAtomsCount(atomNumber) <= 1) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == IBond.Order.SINGLE) {
                    IAtomType type = getAtomType("I.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            }
            return null;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET && 
               atom.getFormalCharge() != 0) {
      if (atom.getFormalCharge() == -1) {
          if (getConnectedAtomsCount(atomNumber) == 0) {
              IAtomType type = getAtomType("I.minus");
              if (isAcceptable(atomNumber, type)) return type;
          } else {
              IAtomType type = getAtomType("I.minus.5");
              if (isAcceptable(atomNumber, type)) return type;
          }
            } else if (atom.getFormalCharge() == 1) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == IBond.Order.DOUBLE) {
                    IAtomType type = getAtomType("I.plus.sp2");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (maxBondOrder == IBond.Order.SINGLE){
                    IAtomType type = getAtomType("I.plus.sp3");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            }
        } else if (getConnectedAtomsCount(atomNumber) == 3) {
            int doubleBondCount = countAttachedDoubleBonds(atomNumber);
            if (doubleBondCount == 2) {
                IAtomType type = getAtomType("I.5");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (atom.getFormalCharge() != CDKConstants.UNSET
                    && atom.getFormalCharge() == 0) {
                IAtomType type = getAtomType("I.sp3d2.3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (getConnectedAtomsCount(atomNumber) == 2) {
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (maxBondOrder == IBond.Order.DOUBLE) {
                IAtomType type = getAtomType("I.3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (getConnectedAtomsCount(atomNumber) == 1 ||
                getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("I");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveRuthenium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0) {
            IAtomType type = getAtomType("Ru.6");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -2) {
            IAtomType type = getAtomType("Ru.2minus.6");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -3) {
            IAtomType type = getAtomType("Ru.3minus.6");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceivePotassium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == +1)) {
            IAtomType type = getAtomType("K.plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() == CDKConstants.UNSET
                || atom.getFormalCharge() == 0) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 1) {
                IAtomType type = getAtomType("K.neutral");
                if (isAcceptable(atomNumber, type)) return type;
            }
            IAtomType type = getAtomType("K.metallic");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceivePlutonium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
       if (atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Pu");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveCadmium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +2)) {
            IAtomType type = getAtomType("Cd.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)) {
            if (getConnectedAtomsCount(atomNumber) == 0) {
                IAtomType type = getAtomType("Cd.metallic");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (getConnectedAtomsCount(atomNumber) == 2) {
                IAtomType type = getAtomType("Cd.2");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }
    
    private IAtomType perceiveIndium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 3) {
            IAtomType type = getAtomType("In.3");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() == 3 && getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("In.3plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 1) {
            IAtomType type = getAtomType("In.1");
            if (isAcceptable(atomNumber, type)) return type;
        } else {
            IAtomType type = getAtomType("In");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveChlorine(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            if (getConnectedAtomsCount(atomNumber) > 1) {
                if (atom.getFormalCharge() != CDKConstants.UNSET
                        && atom.getFormalCharge() == +1) {
                    IAtomType type = getAtomType("Cl.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (getConnectedAtomsCount(atomNumber) == 1) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == IBond.Order.SINGLE) {
                    IAtomType type = getAtomType("Cl.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (getConnectedAtomsCount(atomNumber) == 0
                    && (atom.getFormalCharge() == CDKConstants.UNSET
                    || atom.getFormalCharge() == 0)) {
                IAtomType type = getAtomType("Cl.radical");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (atom.getFormalCharge() == CDKConstants.UNSET
                || atom.getFormalCharge() == 0) {
            int neighborcount = getConnectedAtomsCount(atomNumber);
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);

            if (maxBondOrder == IBond.Order.DOUBLE) {
                int neighbor = getConnectedAtomsCount(atomNumber);
                if (neighbor == 2) {
                    IAtomType type = getAtomType("Cl.2");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (neighbor == 3) {
                    IAtomType type = getAtomType("Cl.chlorate");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (neighbor == 4) {
                    IAtomType type = getAtomType("Cl.perchlorate");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (neighborcount <= 1) {
                IAtomType type = getAtomType("Cl");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -1)) {
            IAtomType type = getAtomType("Cl.minus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET && atom.getFormalCharge() == 1) {
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (maxBondOrder == IBond.Order.DOUBLE) {
                IAtomType type = getAtomType("Cl.plus.sp2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (maxBondOrder == IBond.Order.SINGLE) {
                IAtomType type = getAtomType("Cl.plus.sp3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == +3) && getConnectedAtomsCount(atomNumber) == 4) {
            IAtomType type = getAtomType("Cl.perchlorate.charged");
            if (isAcceptable(atomNumber, type)) return type;
        } else {
            int doubleBonds = countAttachedDoubleBonds(atomNumber);
            if (getConnectedAtomsCount(atomNumber) == 3
                    && doubleBonds == 2) {
                IAtomType type = getAtomType("Cl.chlorate");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (getConnectedAtomsCount(atomNumber) == 4
                    && doubleBonds == 3) {
                IAtomType type = getAtomType("Cl.perchlorate");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }
    
    private IAtomType perceiveSilver(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 1) {
                IAtomType type = getAtomType("Ag.1");
                if (isAcceptable(atomNumber, type)) return type;
            }
            IAtomType type = getAtomType("Ag.neutral");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 1)) {
            IAtomType type = getAtomType("Ag.plus");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveGold(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        }
        int neighbors = getConnectedAtomsCount(atomNumber);
        if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0) && neighbors == 1) {
            IAtomType type = getAtomType("Au.1");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveRadium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)) {
            IAtomType type = getAtomType("Ra.neutral");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveCalcium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if ("Ca".equals(atom.getSymbol())) {
            if (hasOneSingleElectron(atomNumber)) {
                // no idea how to deal with this yet
                return null;
            } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                    && atom.getFormalCharge() == 2 && getConnectedAtomsCount(atomNumber) == 0)) {
                IAtomType type = getAtomType("Ca.2plus");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                    && atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 2)) {
                IAtomType type = getAtomType("Ca.2");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                    && atom.getFormalCharge() == 0 && getConnectedAtomsCount(atomNumber) == 1)) {
                IAtomType type = getAtomType("Ca.1");
                if (isAcceptable(atomNumber, type)) {
                    return type;
                }
            }
        }
        return null;
    }
    
    private IAtomType perceivePlatinum(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == +2)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 4) {
                IAtomType type = getAtomType("Pt.2plus.4");
                if (isAcceptable(atomNumber, type)) return type;
            } else {
                IAtomType type = getAtomType("Pt.2plus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
                atom.getFormalCharge() == 0)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 2) {
                IAtomType type = getAtomType("Pt.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 4) {
                IAtomType type = getAtomType("Pt.4");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 6) {
                IAtomType type = getAtomType("Pt.6");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }
    
    private IAtomType perceiveAntimony(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                    atom.getFormalCharge() == 0 &&
                    getConnectedAtomsCount(atomNumber) == 3)) {
            IAtomType type = getAtomType("Sb.3");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET && 
                    atom.getFormalCharge() == 0 &&
                    getConnectedAtomsCount(atomNumber) == 4)) {
            IAtomType type = getAtomType("Sb.4");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveGadolinum(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET &&
            atom.getFormalCharge() == +3 &&
            getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Gd.3plus");
            if (isAcceptable(atomNumber, type)) {
                return type;
            }
        }
        return null;
    }

    private IAtomType perceiveMagnesium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                    atom.getFormalCharge() == 0)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 4) {
                IAtomType type = getAtomType("Mg.neutral");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 2) {
                IAtomType type = getAtomType("Mg.neutral.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 1) {
                IAtomType type = getAtomType("Mg.neutral.1");
                if (isAcceptable(atomNumber, type)) return type;
            } else {
                IAtomType type = getAtomType("Mg.neutral");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == +2)) {
            IAtomType type = getAtomType("Mg.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveThallium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET &&
            atom.getFormalCharge() == +1 &&
            getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Tl.plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET &&
                   atom.getFormalCharge() == 0 &&
                   getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Tl");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET &&
                   atom.getFormalCharge() == 0 &&
                   getConnectedAtomsCount(atomNumber) == 1) {
            IAtomType type = getAtomType("Tl.1");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveLead(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET &&
            atom.getFormalCharge() == 0 &&
            getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Pb.neutral");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET &&
                   atom.getFormalCharge() == 2 &&
                   getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Pb.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET &&
                   atom.getFormalCharge() == 0 &&
                   getConnectedAtomsCount(atomNumber) == 1) {
            IAtomType type = getAtomType("Pb.1");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveStrontium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 2)) {
            IAtomType type = getAtomType("Sr.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveTitanium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET &&
            atom.getFormalCharge() == -3 &&
            getConnectedAtomsCount(atomNumber) == 6) {
            IAtomType type = getAtomType("Ti.3minus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
                    atom.getFormalCharge() == 0) &&
                   getConnectedAtomsCount(atomNumber) == 4) {
            IAtomType type = getAtomType("Ti.sp3");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == 0)
                && getConnectedAtomsCount(atomNumber) == 2) {
            IAtomType type = getAtomType("Ti.2");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveVanadium(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == -3 &&
                getConnectedAtomsCount(atomNumber) == 6) {
            IAtomType type = getAtomType("V.3minus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() != CDKConstants.UNSET
                && atom.getFormalCharge() == -3
                && getConnectedAtomsCount(atomNumber) == 4) {
            IAtomType type = getAtomType("V.3minus.4");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private IAtomType perceiveBromine(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            if (getConnectedAtomsCount(atomNumber) == 0) {
                if (atom.getFormalCharge() != CDKConstants.UNSET &&
                        atom.getFormalCharge() == +1) {
                    IAtomType type = getAtomType("Br.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                } else if (atom.getFormalCharge() == CDKConstants.UNSET ||
                        atom.getFormalCharge() == 0) {
                    IAtomType type = getAtomType("Br.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            } else if (getConnectedAtomsCount(atomNumber) <= 1) {
                IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
                if (maxBondOrder == IBond.Order.SINGLE) {
                    IAtomType type = getAtomType("Br.plus.radical");
                    if (isAcceptable(atomNumber, type)) return type;
                }
            }
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == -1)) {
            IAtomType type = getAtomType("Br.minus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (atom.getFormalCharge() == 1) {
            IBond.Order maxBondOrder = getMaximumBondOrder(atomNumber);
            if (maxBondOrder == IBond.Order.DOUBLE) {
                IAtomType type = getAtomType("Br.plus.sp2");
                if (isAcceptable(atomNumber, type)) return type;
            }else if (maxBondOrder == IBond.Order.SINGLE){
                IAtomType type = getAtomType("Br.plus.sp3");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if (getConnectedAtomsCount(atomNumber) == 1 ||
                getConnectedAtomsCount(atomNumber) == 0) {
            IAtomType type = getAtomType("Br");
            if (isAcceptable(atomNumber, type)) return type;
        } else if (getConnectedAtomsCount(atomNumber) == 3) {
            IAtomType type = getAtomType("Br.3");
            if (isAcceptable(atomNumber, type)) return type;
        }
        return null;
    }
    
    private int countAttachedDoubleBonds(int atomNumber, String symbol) {
        return countAttachedBonds(atomNumber, IBond.Order.DOUBLE, symbol);
    }
    
    private IAtomType perceiveCobalt(int atomNumber) throws CDKException {
    	IAtom atom = atomContainer.getAtom(atomNumber);
        if (hasOneSingleElectron(atomNumber)) {
            // no idea how to deal with this yet
            return null;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == +2)) {
            IAtomType type = getAtomType("Co.2plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() != CDKConstants.UNSET &&
                atom.getFormalCharge() == +3)) {
            IAtomType type = getAtomType("Co.3plus");
            if (isAcceptable(atomNumber, type)) return type;
        } else if ((atom.getFormalCharge() == CDKConstants.UNSET ||
                atom.getFormalCharge() == 0)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 2) {
                IAtomType type = getAtomType("Co.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 4) {
                IAtomType type = getAtomType("Co.4");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 6) {
                IAtomType type = getAtomType("Co.6");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 1) {
                IAtomType type = getAtomType("Co.1");
                if (isAcceptable(atomNumber, type)) return type;
            } else {
                IAtomType type = getAtomType("Co.metallic");
                if (isAcceptable(atomNumber, type)) return type;
            }
        } else if ((atom.getFormalCharge() != null
                && atom.getFormalCharge() == +1)) {
            int neighbors = getConnectedAtomsCount(atomNumber);
            if (neighbors == 2) {
                IAtomType type = getAtomType("Co.plus.2");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 4) {
                IAtomType type = getAtomType("Co.plus.4");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 1) {
                IAtomType type = getAtomType("Co.plus.1");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 6) {
                IAtomType type = getAtomType("Co.plus.6");
                if (isAcceptable(atomNumber, type)) return type;
            } else if (neighbors == 5) {
                IAtomType type = getAtomType("Co.plus.5");
                if (isAcceptable(atomNumber, type)) return type;
            } else {
                IAtomType type = getAtomType("Co.plus");
                if (isAcceptable(atomNumber, type)) return type;
            }
        }
        return null;
    }

	private int countAttachedDoubleBonds(int atomNumber) {
    	return piBondCounts[atomNumber];
    }
    
    private int countAttachedSingleBonds(int atomNumber) {
        return singleBondCounts[atomNumber];
    }

    private boolean hasAromaticBond(int atomNumber) {
        List<IBond> neighbors = atomContainer.getConnectedBondsList(atomContainer.getAtom(atomNumber));
        for (IBond bond : neighbors) {
            if (bond.getFlag(CDKConstants.ISAROMATIC)) return true;
        }
        return false;
    }

    /**
     * Count the number of doubly bonded atoms.
     *
     * @param container the molecule in which to look
     * @param atom the atom being looked at
     * @param order the desired bond order of the attached bonds 
     * @param symbol If not null, then it only counts the double bonded atoms which
     *               match the given symbol.
     * @return the number of doubly bonded atoms
     */
    private int countAttachedBonds(int atomNumber, IBond.Order order, String symbol) {
    	// count the number of double bonded oxygens
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	List<IBond> neighbors = atomContainer.getConnectedBondsList(atom);
    	int neighborcount = neighbors.size();
    	int doubleBondedAtoms = 0;
    	for (int i=neighborcount-1;i>=0;i--) {
            IBond bond =  neighbors.get(i);
    		if (bond.getOrder() == order) {
    			if (bond.getAtomCount() == 2 && bond.contains(atom)) {
    				if (symbol != null) {
    				    // if other atom is a sulphur
    					if ((bond.getAtom(0) != atom &&
    					     bond.getAtom(0).getSymbol().equals(symbol)) ||
    						(bond.getAtom(1) != atom &&
    						 bond.getAtom(1).getSymbol().equals(symbol))) {
    						doubleBondedAtoms++;
    					}
    				} else {
    					doubleBondedAtoms++;
    				}
    			}
    		}
    	}
    	return doubleBondedAtoms;
    }

    private IAtomType getAtomType(String identifier) throws CDKException {
    	IAtomType type = factory.getAtomType(identifier);
    	type.setValency((Integer)type.getProperty(CDKConstants.PI_BOND_COUNT) +
    			        type.getFormalNeighbourCount());
    	return type;
    }
    
    private boolean isAcceptable(int atomNumber, IAtomType type) {
    	IAtom atom = atomContainer.getAtom(atomNumber);
    	if (mode == REQUIRE_EXPLICIT_HYDROGENS) {
    		// make sure no implicit hydrogens were assumed
    		int actualContainerCount = atomContainer.getConnectedAtomsCount(atom);
    		int requiredContainerCount = type.getFormalNeighbourCount();
    		if (actualContainerCount != requiredContainerCount)
    			return false;
    	} else if (atom.getImplicitHydrogenCount() != CDKConstants.UNSET) {
    		// confirm correct neighbour count
    		int connectedAtoms = atomContainer.getConnectedAtomsCount(atom);
    		int hCount = atom.getImplicitHydrogenCount();
    		int actualNeighbourCount =  connectedAtoms + hCount;
    		int requiredNeighbourCount = type.getFormalNeighbourCount();
    		if (actualNeighbourCount > requiredNeighbourCount)
    			return false;
    	}

    	// confirm correct bond orders
        IBond.Order typeOrder = type.getMaxBondOrder(); 
    	if (typeOrder != null) {
    		for (IBond bond : atomContainer.getConnectedBondsList(atom)) {
    			IBond.Order order = bond.getOrder();
    			if (order != CDKConstants.UNSET && order != IBond.Order.UNSET) {
    				if (BondManipulator.isHigherOrder(order, typeOrder)) return false;
    			} else if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE)) {
    				if (typeOrder != IBond.Order.SINGLE &&
        				typeOrder != IBond.Order.DOUBLE) return false;
    			} else {
    				return false;
    			}
    		}
    	}
    		
    	// confirm correct valency
    	if (type.getValency() != CDKConstants.UNSET && atomContainer.getBondOrderSum(atom) > type.getValency())
    		return false;

    	// confirm correct formal charge
        if (atom.getFormalCharge() != CDKConstants.UNSET &&
            !atom.getFormalCharge().equals(type.getFormalCharge()))
    		return false;

    	return true;
    }
    
    private boolean isHueckelNumber(int electronCount) {
        return (electronCount % 4 == 2) && (electronCount >= 2);
    }
    
}

