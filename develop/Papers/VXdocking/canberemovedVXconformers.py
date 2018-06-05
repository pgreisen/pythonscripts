from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms as rdmt
import numpy as np
# from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Geometry import rdGeometry as geom


class VXconformers:


    def __init__(self):
        # what is this value???
        # distance between P-N
        # remove self interaction
        self.distance_threshold = 3.0

        self.p_id = None
        self.s_id  = None
        self.pc_id  = None
        self.nh_id  = None
        self.po_id  = None
        self.poe_id  = None



    def get_conformer_rmsd(self,mol):
        """
        Calculate conformer-conformer RMSD.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        rmsd = np.zeros((mol.GetNumConformers(), mol.GetNumConformers()),
                        dtype=float)
        for i, ref_conf in enumerate(mol.GetConformers()):
            for j, fit_conf in enumerate(mol.GetConformers()):
                if i >= j:
                    continue
                rmsd[i, j] = AllChem.GetBestRMS(mol, mol, ref_conf.GetId(),
                                                fit_conf.GetId())
                rmsd[j, i] = rmsd[i, j]
        return rmsd

    def get_conformer_energies(self, mol):
        """
        Calculate conformer energies.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        Returns
        -------
        energies : array_like
            Minimized conformer energies.
        """
        energies = []
        for conf in mol.GetConformers():
            ff = self.get_molecule_force_field(mol, conf_id=conf.GetId())
            energy = ff.CalcEnergy()
            energies.append(energy)
        energies = np.asarray(energies, dtype=float)
        return energies


    def set_new_bond_distance(self, mol, dist, atm1, atm2, confid=0):
        """
        Set a new bond distance between two atoms
        ----------
        mol : RDKit molecule
        dist : float - new bond distance
        atm1 : int - idx of atom 1
        atm2 : int - idx of atom 2
        confid : int - conformer from ensemble
        Returns
        ----------
        Mol with new bond length between atom 1 and atom 2

        """
        conf = mol.GetConformer(confid)
        rdmt.SetBondLength(conf, atm1, atm2, dist)
        return mol


    def set_new_angle(self,mol, angle, atm1, atm2, atm3, confid=0):
        """
        Set a new angle between three atoms
        ----------
        mol : RDKit molecule
        angle : float - new angle in degrees
        atm1 : int - idx of atom 1
        atm2 : int - idx of atom 2
        atm3 : int - idx of atom 3
        confid : int - conformer from ensemble
        Returns
        ----------
        Mol with new angle between atom 1, atom 2, atom 3

        """
        conf = mol.GetConformer(confid)
        rdmt.SetAngleDeg(conf, atm1, atm2, atm3, angle)
        return mol


    def generate_molecule(self, name, smiles):
        """
        Generate the 3D molecular structure based on input SMILES
        ----------
        name : name of molecule
        smiles: SMILES of molecule
        Returns
        ----------
        Mol

        """
        LIGAND_NAME = name
        m = Chem.MolFromSmiles(smiles)
        m_h = Chem.AddHs(m)
        # Embeed the geometry
        AllChem.EmbedMolecule(m_h, AllChem.ETKDG())
        # Setting name of molecule
        m_h.SetProp("_Name", LIGAND_NAME)
        return m_h



    def get_conformers(self,mol,nr=500,rmsthreshold=0.1):
        """
        Generate 3D conformers of molecule using CSD-method
        ----------
        mol : RKdit molecule
        nr : integer, number of conformers to be generate
        rmsthreshold : float, prune conformers that are less rms away from another conf
        Returns
        ----------
        List of new conformation IDs
        """
        # Generate conformers on the CSD-method
        return AllChem.EmbedMultipleConfs(mol, numConfs=nr,useBasicKnowledge=True,\
                                          pruneRmsThresh=rmsthreshold,useExpTorsionAnglePrefs=True)


    def get_index_atms(self,mol):
        """
        get the index of atoms to set geometry
        methyl bonded to phosphorus
        ether bonded to phosphorus
        phosphoryl bonded to phosphorus
        sulfur atom bonded to phosphorus
        hydrogen bonded to nitrogen
        ----------
        mol: RKdit mol
        Returns
        ----------
        index of atoms

        """
        chiral_atm = {}
        o_groups = []
        for atom in mol.GetAtoms():
            tmp_atm = atom.GetSymbol()
            bonds = atom.GetBonds()
            if (tmp_atm == 'P'):
                self.p_id = atom.GetIdx()
                for bond in bonds:
                    target = bond.GetOtherAtom(atom)
                    # methyl group
                    if (target.GetSymbol() == 'C'):
                        self.pc_id = target.GetIdx()
                    elif (target.GetSymbol() == 'O'):
                        if (str(bond.GetBondType()) == 'DOUBLE'):
                            self.po_id = target.GetIdx()
                        else:
                            self.poe_id = target.GetIdx()
            elif (tmp_atm == 'S'):
                self.s_id = atom.GetIdx()
            elif (tmp_atm == 'N'):
                for bond in bonds:
                    target = bond.GetOtherAtom(atom)
                    if (target.GetSymbol() == 'H'):
                        self.nh_id = atom.GetIdx()
            else:
                for bond in bonds:
                    target = bond.GetOtherAtom(atom)
        return self.p_id, self.s_id, self.pc_id, self.nh_id, self.po_id, self.poe_id

    def get_ts_geom(self, m, p_s_dist, p_id, s_id, pc_id, nh_id, po_id, poe_id):
        '''


        '''
        plane_angle = 120
        off_plane = 90
        # P-S bond distance is set here
        m = self.set_new_bond_distance(m, p_s_dist, p_id, s_id)
        # setting o-ethyl group
        m = self.set_new_angle(m, 90, s_id, p_id, poe_id)
        # setting methyl-group
        m = self.set_new_angle(m, 90, s_id, p_id, pc_id)
        # setting phosphoryl
        m = self.set_new_angle(m, 90, s_id, p_id, po_id)
        # Setting plane angle to 120 degrees
        m = self.set_new_angle(m, 120, poe_id, p_id, pc_id)
        m = self.set_new_angle(m, 120, po_id, p_id, pc_id)
        # conformers
        conf = m.GetConformer(0)
        # Proton sulfur distance
        dst1 = rdmt.GetBondLength(conf, po_id, nh_id)
        if (dst1 > self.distance_threshold):
            return m
        else:
            return

    def get_transition_state_geometry(self, p_s_dist, p_id, s_id, pc_id, nh_id, po_id, poe_id ):
        '''
        ----------
        dist: float, bond breaking distance

        Returns
        ----------
        index of atoms
        '''
        import os
        fls = os.listdir('.')
        s_isomers = []
        r_isomers = []
        k = True

        for fl in fls:
            if (fl.startswith('prune')):
                m = Chem.MolFromMolFile(fl, True, False)
                if( m is None):
                    continue

                Chem.AssignAtomChiralTagsFromStructure(m)
                tmp = Chem.FindMolChiralCenters(m, includeUnassigned=True)

                if (tmp[0][1] == 'S'):
                    m = self.get_ts_geom(m, p_s_dist, p_id, s_id, pc_id, nh_id, po_id, poe_id)
                    if (m != None):
                        s_isomers.append(m)
                elif (tmp[0][1] == 'R'):
                    m = self.get_ts_geom(m, p_s_dist, p_id, s_id, pc_id, nh_id, po_id, poe_id)

                    if (m != None):
                        r_isomers.append(m)
                else:
                    print "Error"
        return s_isomers, r_isomers

    def write_aligned_to_file(self, list_of_confs, atoms_to_match=(8, 9, 10, 11), filename='Aligned.sdf'):
        aligned = None
        aligned = list_of_confs[0]
        tmp_length_ = len(list_of_confs)
        for i in range(tmp_length_):
            aligned.AddConformer(list_of_confs[i].GetConformer(0), assignId=True)
        conf_ids = [conf.GetId() for conf in aligned.GetConformers()]
        rmslst = []
        AllChem.AlignMolConformers(aligned, atomIds=atoms_to_match, RMSlist=rmslst)

        for i in conf_ids:
            writer3 = Chem.SDWriter(filename + "_" + str(i) + ".sdf")
            writer3.write(aligned, confId=i)
        return aligned


    def get_index_of_nucleophile(self, combo):
        '''

        '''
        o_hydroxide = 0
        for atom in combo.GetAtoms():
            bonds = atom.GetBonds()
            for bond in bonds:
                target = bond.GetOtherAtom(atom)
                if (target.GetSymbol() == 'H' and atom.GetSymbol() == 'O'):
                    o_hydroxide = atom.GetIdx()
                    h_hydroxide = target.GetIdx()
        return o_hydroxide, h_hydroxide


    # setup and generate transition states
    def add_hydroxide_to_transition_state_model(self, nucleophile, p_id, s_id, conformer, bond_increase=1.9):
        '''


        '''
        # get directional vector from phosphorus and sulfur atom
        p_vec = np.array(conformer.GetConformer(0).GetAtomPosition(p_id))
        s_vec = np.array(conformer.GetConformer(0).GetAtomPosition(s_id))
        sp_vec = (p_vec - s_vec) / np.linalg.norm(p_vec - s_vec) * bond_increase + p_vec
        tmp_ = geom.Point3D(sp_vec[0], sp_vec[1], sp_vec[2])
        nucleophile.GetConformer().SetAtomPosition(0, tmp_)
        # Hydrogen is removed and add later
        nucleophile = Chem.RemoveHs(nucleophile)
        nucleophile = Chem.AddHs(nucleophile, addCoords=True)
        # combine molecule
        combo = Chem.CombineMols(conformer, nucleophile)
        # get new index of oxygen of hydroxide
        o_hydroxide, h_hydroxide = self.get_index_of_nucleophile(combo)
        edcombo = Chem.EditableMol(combo)
        edcombo.AddBond(p_id, o_hydroxide, order=Chem.rdchem.BondType.SINGLE)
        back = edcombo.GetMol()
        # rdmt.SetAngleDeg(back.GetConformer(0),p_id,o_hydroxide, h_hydroxide,180)
        return back





