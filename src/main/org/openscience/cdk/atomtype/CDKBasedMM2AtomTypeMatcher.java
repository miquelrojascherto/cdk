/* Copyright (C) 2008-2012  Egon Willighagen <egonw@users.sf.net>
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

import java.io.InputStream;
import java.util.Hashtable;
import java.util.Map;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.atomtype.mapper.AtomTypeMapper;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

/**
 * Atom Type matcher for MM2 atom types. It uses the {@link CDKAtomTypeMatcher}
 * for perception and then maps CDK to MM2 atom types.
 *
 * @cdk.module     atomtype
 * @cdk.githash
 * @cdk.keyword    atom type, MM2
 * @cdk.keyword    MM2
 */
@TestClass("org.openscience.cdk.atomtype.SybylAtomTypeMatcherTest")
public class CDKBasedMM2AtomTypeMatcher implements IAtomTypeMatcher {

    private final static String SYBYL_ATOM_TYPE_LIST = "org/openscience/cdk/dict/data/mm2-atom-types.owl";
    private final static String CDK_TO_MM2_MAP = "org/openscience/cdk/dict/data/cdk-mm2-mappings.owl";

	private AtomTypeFactory factory;
	private CDKAtomTypeMatcher cdkMatcher;
	private AtomTypeMapper mapper;

    private static Map<IChemObjectBuilder,CDKBasedMM2AtomTypeMatcher>
        factories = new Hashtable<IChemObjectBuilder,CDKBasedMM2AtomTypeMatcher>(1);

    private CDKBasedMM2AtomTypeMatcher(IChemObjectBuilder builder) {
        InputStream stream = this.getClass().getClassLoader().getResourceAsStream(SYBYL_ATOM_TYPE_LIST);
        factory = AtomTypeFactory.getInstance(stream, "owl", builder);
        cdkMatcher = CDKAtomTypeMatcher.getInstance(builder);
        InputStream mapStream = this.getClass().getClassLoader().getResourceAsStream(CDK_TO_MM2_MAP);
        mapper = AtomTypeMapper.getInstance(CDK_TO_MM2_MAP, mapStream);
    }

    @TestMethod("testGetInstance_IChemObjectBuilder")
    public static CDKBasedMM2AtomTypeMatcher getInstance(IChemObjectBuilder builder) {
    	if (!factories.containsKey(builder))
    		factories.put(builder, new CDKBasedMM2AtomTypeMatcher(builder));
    	return factories.get(builder);
    }

    @TestMethod("testFindMatchingAtomType_IAtomContainer")
    public IAtomType[] findMatchingAtomType(IAtomContainer atomContainer) throws CDKException {
        for (IAtom atom : atomContainer.atoms()) {
            IAtomType type = cdkMatcher.findMatchingAtomType(atomContainer, atom);
            atom.setAtomTypeName(type == null ? null : type.getAtomTypeName());
            atom.setHybridization(type == null ? null : type.getHybridization());
        }
        CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
        IAtomType[] types = new IAtomType[atomContainer.getAtomCount()];
        int typeCounter = 0;
        for (IAtom atom : atomContainer.atoms()) {
            String mappedType = mapCDKToMM2Type(atom);
            if (mappedType == null) {
                types[typeCounter] = null;
            } else {
                types[typeCounter] = factory.getAtomType(mappedType);
            }
            typeCounter++;
        }
        return types;
    }

    /**
     * MM2 atom type perception for a single atom.
     */
    @TestMethod("testFindMatchingAtomType_IAtomContainer_IAtom")
    public IAtomType findMatchingAtomType(IAtomContainer atomContainer, IAtom atom)
        throws CDKException {
        IAtomType type = cdkMatcher.findMatchingAtomType(atomContainer, atom);
        if (type == null) return null;
        else atom.setAtomTypeName(type.getAtomTypeName());
        String mappedType = mapCDKToMM2Type(atom);
        if (mappedType == null) return null;
        return factory.getAtomType(mappedType);
    }

    private String mapCDKToMM2Type(IAtom atom) {
        String typeName = atom.getAtomTypeName();
        if (typeName == null) return null;
        String mappedType = mapper.mapAtomType(typeName);
        // the next lines should perceive more details atom types
        return mappedType;
    }

}

