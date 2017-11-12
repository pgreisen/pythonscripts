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
inFile = open("cc_vector_104.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_104, protein_rmsf
alter cc_residue_104, b=0.0
alter cc_residue_104 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_104 and n. CA
inFile = open("cc_vector_41.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_41, protein_rmsf
alter cc_residue_41, b=0.0
alter cc_residue_41 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_41 and n. CA
inFile = open("cc_vector_37.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_37, protein_rmsf
alter cc_residue_37, b=0.0
alter cc_residue_37 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_37 and n. CA
inFile = open("cc_vector_73.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_73, protein_rmsf
alter cc_residue_73, b=0.0
alter cc_residue_73 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_73 and n. CA
inFile = open("cc_vector_88.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_88, protein_rmsf
alter cc_residue_88, b=0.0
alter cc_residue_88 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_88 and n. CA
inFile = open("cc_vector_135.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_135, protein_rmsf
alter cc_residue_135, b=0.0
alter cc_residue_135 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_135 and n. CA
inFile = open("cc_vector_20.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_20, protein_rmsf
alter cc_residue_20, b=0.0
alter cc_residue_20 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_20 and n. CA
inFile = open("cc_vector_45.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_45, protein_rmsf
alter cc_residue_45, b=0.0
alter cc_residue_45 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_45 and n. CA
inFile = open("cc_vector_77.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_77, protein_rmsf
alter cc_residue_77, b=0.0
alter cc_residue_77 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_77 and n. CA
inFile = open("cc_vector_119.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_119, protein_rmsf
alter cc_residue_119, b=0.0
alter cc_residue_119 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_119 and n. CA
inFile = open("cc_vector_137.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_137, protein_rmsf
alter cc_residue_137, b=0.0
alter cc_residue_137 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_137 and n. CA
inFile = open("cc_vector_117.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_117, protein_rmsf
alter cc_residue_117, b=0.0
alter cc_residue_117 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_117 and n. CA
inFile = open("cc_vector_57.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_57, protein_rmsf
alter cc_residue_57, b=0.0
alter cc_residue_57 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_57 and n. CA
inFile = open("cc_vector_24.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_24, protein_rmsf
alter cc_residue_24, b=0.0
alter cc_residue_24 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_24 and n. CA
inFile = open("cc_vector_133.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
create cc_residue_133, protein_rmsf
alter cc_residue_133, b=0.0
alter cc_residue_133 and n. CA, b=stored.newB.pop(0)
spectrum b, rainbow_rev, cc_residue_133 and n. CA
inFile = open("cc_vector_90.dat", 'r')
stored.newB = []
for line in inFile.readlines(): stored.newB.append( float(line) )
inFile.close()
hide everything, elem h
