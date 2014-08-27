import sys,os

class get_atoms_binuclear:

    # Requires file
    # returns file object
    def read_file(self,filename):
        fl = open(file,'r')
        return fl
    
    # Only heteroatoms
    # Requires filename and list of atoms
    # Returns atoms from pdb file
    def get_hetatoms(self,filename,atoms):
        fl = open(filename,'r')
        # List of lines
        atm = {}
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:6] == 'HETATM':
                i = line.split()[2]
                # modify for carbamated lysine residue
                carbamated_resname = str(line[17:20]).strip()
                #                if i in atoms and carbamated_resname != 'KCX':
                if i in atoms and carbamated_resname != 'KCX' or carbamated_resname != 'FMT':
                    atm[i] = line
        return atm


    def get_hetatoms_pdb(self,filename,atoms):
        tmp_het = []
        fl = open(filename,'r')
        atmnr = []
        atmname = []
        print atoms
        for elem in atoms:
            # Need to split into atomnr and name
            tmp_atmnr,tmp_atmname = elem.split()
            atmnr.append(str(tmp_atmnr).strip())
            atmname.append(tmp_atmname[0:2])
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:6] == 'HETATM':

                i = line[7:11]
                j = line[13:15]
                j = j.strip()
                i = i.strip()
                carbamated_resname = str(line[17:20]).strip()
                if i in atmnr:
                    # if j in atmname and carbamated_resname != 'KCX':
                    if i in atoms and carbamated_resname != 'KCX' or carbamated_resname != 'FMT':
                        tmp_het.append(line)
        return tmp_het


    # Requires pdb-file and metal string
    # Return list with metal ions
    def get_metalion_pdb(self,filename,atoms):
        metal_coor = []
        fl = open(filename,'r')
        # List of lines
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:6] == 'HETATM':
                i = line.split()[2]
                metal_resname = str(line[17:20]).strip()
                if i in atoms and metal_resname in atoms:
                    metal_coor.append(line)
        return metal_coor

    # Requires filename, resname, resid, atom
    # Returns line in pdb file
    def get_lig_atoms(self,filename,resn,resid,atom):
        fl = open(filename,'r')
        # List of lines
        atm = {}
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:4] == 'ATOM':
                lst = line.split()
                # changed 06-01-2010
                #if lst[5] == resid:
                if str(line[23:26]).strip() == resid:
                    #if lst[3] == resn:
                    if line[17:20] == resn:
                        # if lst[2] == atom:
                        if str(line[12:16]).strip() == atom:
                            return line                            

    # Requires list of atoms
    # Write cry file with atoms
    def write_atoms_file(self,atoms):
        wr = open('cry.pdb','w')
        #print 'all atoms given',atoms
        for atm in atoms:
            print 'Debug atoms in crystal structure\n',atm
            wr.write(atm)


    def main():
        fl1 = sys.argv[1]
        atoms = []
        position = 0 
        for arg in sys.argv[2:]:
            atoms.append(str(arg))
        ln = get_atoms_file(fl1,atoms)
        write_atoms_file(atoms,ln)
