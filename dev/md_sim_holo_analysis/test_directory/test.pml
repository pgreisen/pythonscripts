hide everything
show cartoon
show sticks, het
create protein_rmsf, all and not het
create torsional_angles, resi 
create rmsd_residues, resi 102+104+41+37+73+88+135+20+45+77+119+137+117+57+24+133+90+
show sticks, torsional_angles
show sticks, rmsd_residues
color red, torsional_angles and name c*
color yellow, rmsd_residues and name c*
inFile = open("bfactor_data", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
alter protein_rmsf, b=0.0
alter protein_rmsf and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, protein_rmsf and n. CA
create p_rmsf, protein_rmsf
spectrum b, rainbow_rev, p_rmsf and n. CA
cartoon putty, p_rmsf
inFile = open("cc_vector_102.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_102, protein_rmsf
alter cc_residue_102, b=0.0
alter cc_residue_102 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_102 and n. CA
