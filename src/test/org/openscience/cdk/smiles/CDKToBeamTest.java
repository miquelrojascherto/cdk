/*
 * Copyright (c) 2013 European Bioinformatics Institute (EMBL-EBI)
 *                    John May <jwmay@users.sf.net>
 *  
 * Contact: cdk-devel@lists.sourceforge.net
 *  
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version. All we ask is that proper credit is given
 * for our work, which includes - but is not limited to - adding the above 
 * copyright notice to the beginning of your source code files, and to any
 * copyright notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 U
 */

package org.openscience.cdk.smiles;

import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.PseudoAtom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.templates.TestMoleculeFactory;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.beam.Graph;
import uk.ac.ebi.beam.Element;
import uk.ac.ebi.beam.Functions;

import java.util.Map;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation.OPPOSITE;
import static org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation.TOGETHER;
import static org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo.ANTI_CLOCKWISE;
import static org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo.CLOCKWISE;

/**
 * Unit tests for converting CDK IAtomContainer's to the grins object module.
 * For clarity often the SMILES output is verified if a test fails it could be
 * the Grins output changed and there was not a problem with the conversion.
 *
 * @author John May
 * @cdk.module test-smiles
 */
public class CDKToBeamTest {

    @Test(expected = NullPointerException.class)
    public void noImplicitHCount() throws Exception {
        new CDKToBeam().toBeamAtom(new Atom("C"));
    }

    @Test(expected = NullPointerException.class)
    public void noSymbol() throws Exception {
        new CDKToBeam().toBeamAtom(new Atom());
    }

    @Test
    public void unknownSymbol() throws Exception {
        IAtom a = new Atom("ALA");
        a.setImplicitHydrogenCount(0);
        assertThat(new CDKToBeam().toBeamAtom(a).element(),
                   is(Element.Unknown));
    }

    @Test
    public void unknownSymbol_Pseudo() throws Exception {
        IAtom a = new PseudoAtom("R1");
        a.setImplicitHydrogenCount(0);
        assertThat(new CDKToBeam().toBeamAtom(a).element(),
                   is(Element.Unknown));
    }

    @Test
    public void methane_Atom() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(4);
        assertThat(new CDKToBeam().toBeamAtom(a).element(),
                   is(Element.Carbon));
        assertThat(new CDKToBeam().toBeamAtom(a).hydrogens(),
                   is(4));
    }

    @Test
    public void water_Atom() throws Exception {
        IAtom a = new Atom("O");
        a.setImplicitHydrogenCount(2);
        assertThat(new CDKToBeam().toBeamAtom(a).element(),
                   is(Element.Oxygen));
        assertThat(new CDKToBeam().toBeamAtom(a).hydrogens(),
                   is(2));
    }

    @Test public void chargedAtom() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(0);
        for (int chg = -10; chg < 10; chg++) {
            a.setFormalCharge(chg);
            assertThat(new CDKToBeam().toBeamAtom(a).charge(),
                       is(chg));
        }
    }

    @Test public void aliphaticAtom() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(0);
        assertFalse(new CDKToBeam().toBeamAtom(a).aromatic());
    }

    @Test public void aromaticAtom() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(0);
        a.setFlag(CDKConstants.ISAROMATIC, true);
        assertTrue(new CDKToBeam().toBeamAtom(a).aromatic());
    }

    @Test public void unspecifiedIsotope() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(0);
        assertThat(new CDKToBeam().toBeamAtom(a).isotope(),
                   is(-1));
    }

    @Test public void specifiedIsotope() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(0);
        a.setMassNumber(13);
        assertThat(new CDKToBeam().toBeamAtom(a).isotope(),
                   is(13));
    }
    
    @Test public void defaultIsotope() throws Exception {
        IAtom a = new Atom("C");
        a.setImplicitHydrogenCount(0);
        a.setMassNumber(12);
        assertThat(new CDKToBeam().toBeamAtom(a).isotope(),
                   is(-1));
    }

    @SuppressWarnings("unchecked")
    @Test(expected = IllegalArgumentException.class)
    public void unsetBondOrder() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v, IBond.Order.UNSET);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        new CDKToBeam().toBeamEdge(b, mock);
    }

    @SuppressWarnings("unchecked")
    @Test(expected = NullPointerException.class)
    public void undefBondOrder() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v, null);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        new CDKToBeam().toBeamEdge(b, mock);
    }

    @SuppressWarnings("unchecked")
    @Test(expected = IllegalArgumentException.class)
    public void tooFewAtoms() throws Exception {
        IBond b = new Bond(new IAtom[]{mock(IAtom.class)});
        new CDKToBeam().toBeamEdge(b, mock(Map.class));
    }

    @SuppressWarnings("unchecked")
    @Test(expected = IllegalArgumentException.class)
    public void tooManyAtoms() throws Exception {
        IBond b = new Bond(new IAtom[]{mock(IAtom.class),
                                       mock(IAtom.class),
                                       mock(IAtom.class)});
        new CDKToBeam().toBeamEdge(b, mock(Map.class));
    }

    @SuppressWarnings("unchecked")
    @Test public void singleBond() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        CDKToBeam c2g = new CDKToBeam();
        assertThat(c2g.toBeamEdge(b, mock),
                   is(uk.ac.ebi.beam.Bond.SINGLE.edge(0, 1)));
    }

    @SuppressWarnings("unchecked")
    @Test public void aromaticBond() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v);
        b.setFlag(CDKConstants.ISAROMATIC, true);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        CDKToBeam c2g = new CDKToBeam();
        assertThat(c2g.toBeamEdge(b, mock),
                   is(uk.ac.ebi.beam.Bond.AROMATIC.edge(0, 1)));
    }

    @SuppressWarnings("unchecked")
    @Test public void doubleBond() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v, IBond.Order.DOUBLE);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        CDKToBeam c2g = new CDKToBeam();
        assertThat(c2g.toBeamEdge(b, mock),
                   is(uk.ac.ebi.beam.Bond.DOUBLE.edge(0, 1)));
    }

    @SuppressWarnings("unchecked")
    @Test public void tripleBond() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v, IBond.Order.TRIPLE);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        CDKToBeam c2g = new CDKToBeam();
        assertThat(c2g.toBeamEdge(b, mock),
                   is(uk.ac.ebi.beam.Bond.TRIPLE.edge(0, 1)));
    }

    @SuppressWarnings("unchecked")
    @Test public void quadrupleBond() throws Exception {
        IAtom u = mock(IAtom.class);
        IAtom v = mock(IAtom.class);
        IBond b = new Bond(u, v, IBond.Order.QUADRUPLE);
        Map mock = mock(Map.class);
        when(mock.get(u)).thenReturn(0);
        when(mock.get(v)).thenReturn(1);
        CDKToBeam c2g = new CDKToBeam();
        assertThat(c2g.toBeamEdge(b, mock),
                   is(uk.ac.ebi.beam.Bond.QUADRUPLE.edge(0, 1)));
    }

    @Test public void adeneine() throws Exception {
        Graph g = convert(TestMoleculeFactory.makeAdenine());
        assertThat(g.toSmiles(),
                   is("[C]-1-2=[C](-[N]=[CH]-[N]=[C]1-[NH2])-[NH]-[CH]=[N]2"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("C12=C(N=CN=C1N)NC=N2"));
    }

    @Test public void benzene_kekule() throws Exception {
        Graph g = convert(TestMoleculeFactory.makeBenzene());
        assertThat(g.toSmiles(),
                   is("[CH]=1-[CH]=[CH]-[CH]=[CH]-[CH]1"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("C=1C=CC=CC1"));
    }

    @Test public void benzene() throws Exception {
        IAtomContainer ac = TestMoleculeFactory.makeBenzene();
        Graph g = convert(ac, true, true);
        assertThat(g.toSmiles(),
                   is("[cH]:1:[cH]:[cH]:[cH]:[cH]:[cH]1"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("c1ccccc1"));
    }

    @Test public void imidazole_kekule() throws Exception {
        Graph g = convert(TestMoleculeFactory
                                          .makeImidazole(), false, true);
        assertThat(g.toSmiles(),
                   is("[CH]=1-[NH]-[CH]=[N]-[CH]1"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("C=1NC=NC1"));
    }

    @Test public void imidazole() throws Exception {
        Graph g = convert(TestMoleculeFactory
                                          .makeImidazole(), true, true);
        assertThat(g.toSmiles(),
                   is("[cH]:1:[nH]:[cH]:[n]:[cH]1"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("c1[nH]cnc1"));
    }

    @Test public void C13_isomeric() throws Exception {
        IAtomContainer ac = new AtomContainer();
        IAtom a = new Atom("C");
        a.setMassNumber(13);
        ac.addAtom(a);
        Graph g = convert(ac);
        assertThat(g.atom(0).isotope(), is(13));
        assertThat(g.toSmiles(), is("[13CH4]"));
    }

    @Test public void C13_nonIsomeric() throws Exception {
        IAtomContainer ac = new AtomContainer();
        IAtom a = new Atom("C");
        a.setMassNumber(13);
        ac.addAtom(a);
        Graph g = convert(ac, false, false); // non-isomeric
        assertThat(g.atom(0).isotope(), is(-1));
        assertThat(g.toSmiles(), is("[CH4]"));
    }

    @Test public void azanium() throws Exception {
        IAtomContainer ac = new AtomContainer();
        IAtom a = new Atom("N");
        a.setFormalCharge(+1);
        ac.addAtom(a);
        Graph g = convert(ac);
        assertThat(g.atom(0).charge(), is(+1));
        assertThat(g.toSmiles(), is("[NH4+]"));
    }

    @Test public void oxidanide() throws Exception {
        IAtomContainer ac = new AtomContainer();
        IAtom a = new Atom("O");
        a.setFormalCharge(-1);
        ac.addAtom(a);
        Graph g = convert(ac);
        assertThat(g.atom(0).charge(), is(-1));
        assertThat(g.toSmiles(), is("[OH-]"));
    }

    @Test public void oxidandiide() throws Exception {
        IAtomContainer ac = new AtomContainer();
        IAtom a = new Atom("O");
        a.setFormalCharge(-2);
        ac.addAtom(a);
        Graph g = convert(ac);
        assertThat(g.atom(0).charge(), is(-2));
        assertThat(g.toSmiles(), is("[O-2]"));
    }

    /**
     * (E)-1,2-difluoroethene
     *
     * @cdk.inchi InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1+
     */
    @Test public void e_1_2_difluoroethene() throws Exception {

        IAtomContainer ac = new AtomContainer();
        ac.addAtom(new Atom("F"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("F"));
        ac.addBond(0, 1, SINGLE);
        ac.addBond(1, 2, DOUBLE);
        ac.addBond(2, 3, SINGLE);

        ac.addStereoElement(new DoubleBondStereochemistry(ac.getBond(1),
                                                          new IBond[]{
                                                                  ac.getBond(0),
                                                                  ac.getBond(2)
                                                          },
                                                          OPPOSITE));
        Graph g = convert(ac);
        assertThat(g.toSmiles(),
                   is("[F]/[CH]=[CH]/[F]"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("F/C=C/F"));
    }

    /**
     * (Z)-1,2-difluoroethene
     *
     * @cdk.inchi InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1-
     */
    @Test public void z_1_2_difluoroethene() throws Exception {

        IAtomContainer ac = new AtomContainer();
        ac.addAtom(new Atom("F"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("F"));
        ac.addBond(0, 1, SINGLE);
        ac.addBond(1, 2, DOUBLE);
        ac.addBond(2, 3, SINGLE);

        ac.addStereoElement(new DoubleBondStereochemistry(ac.getBond(1),
                                                          new IBond[]{
                                                                  ac.getBond(0),
                                                                  ac.getBond(2)
                                                          },
                                                          TOGETHER));
        Graph g = convert(ac);
        assertThat(g.toSmiles(),
                   is("[F]/[CH]=[CH]\\[F]"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("F/C=C\\F"));
    }

    /**
     * (2R)-butan-2-ol
     *
     * @cdk.inchi InChI=1/C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3/t4-/s2
     */
    @Test public void _2R_butan_2_ol() throws Exception {

        IAtomContainer ac = new AtomContainer();
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("O"));
        ac.addAtom(new Atom("H"));
        ac.addBond(0, 1, SINGLE);
        ac.addBond(1, 2, SINGLE);
        ac.addBond(2, 3, SINGLE);
        ac.addBond(2, 4, SINGLE);
        ac.addBond(2, 5, SINGLE);

        ac.addStereoElement(new TetrahedralChirality(ac.getAtom(2),
                                                     new IAtom[]{
                                                             ac.getAtom(1), // C-C
                                                             ac.getAtom(3), // C
                                                             ac.getAtom(4), // O
                                                             ac.getAtom(5), // H
                                                     },
                                                     CLOCKWISE));

        Graph g = convert(ac);
        assertThat(g.toSmiles(),
                   is("[CH3]-[CH2]-[C@@](-[CH3])(-[OH])-[H]"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("CC[C@@](C)(O)[H]"));
    }

    /**
     * (2S)-butan-2-ol
     *
     * @cdk.inchi InChI=1/C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3/t4-/s2
     */
    @Test public void _2S_butan_2_ol() throws Exception {

        IAtomContainer ac = new AtomContainer();
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("C"));
        ac.addAtom(new Atom("O"));
        ac.addAtom(new Atom("H"));
        ac.addBond(0, 1, SINGLE);
        ac.addBond(1, 2, SINGLE);
        ac.addBond(2, 3, SINGLE);
        ac.addBond(2, 4, SINGLE);
        ac.addBond(2, 5, SINGLE);

        ac.addStereoElement(new TetrahedralChirality(ac.getAtom(2),
                                                     new IAtom[]{
                                                             ac.getAtom(1), // C-C
                                                             ac.getAtom(3), // C
                                                             ac.getAtom(4), // O
                                                             ac.getAtom(5), // H
                                                     },
                                                     ANTI_CLOCKWISE));

        Graph g = convert(ac);
        assertThat(g.toSmiles(),
                   is("[CH3]-[CH2]-[C@](-[CH3])(-[OH])-[H]"));
        assertThat(Functions.collapse(g).toSmiles(),
                   is("CC[C@](C)(O)[H]"));
    }
    
    static Graph convert(IAtomContainer ac) throws Exception {
        return convert(ac, false, true);
    }

    static Graph convert(IAtomContainer ac,
                                 boolean aromatic,
                                 boolean isomeric) throws
                                                   Exception {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
        CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance())
                        .addImplicitHydrogens(ac);
        if (aromatic)
            CDKHueckelAromaticityDetector.detectAromaticity(ac);
        return new CDKToBeam(isomeric).toBeamGraph(ac);
    }
}
