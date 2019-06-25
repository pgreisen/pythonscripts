import sys, shutil

# map residues

res_map = {
    'C': 'CYS',
    'H': 'HIS',
    'D': 'ASP',
    'E': 'GLU',
    'K': 'LYS',
    'S': 'SER',
    'R': 'ARG',
    'T': 'THR',
    'N': 'ASN',
    'Q': 'GLN',
    'Y': 'TYR',
    'M': 'MET',
    'G': 'GLY',
    'W': 'TRP',
    'P': 'PRO',
    'A': 'ALA',
    'V': 'VAL',
    'I': 'ILE',
    'L': 'LEU',
    'F': 'PHE'}

# xml template
def generate_xml_file(muttask, mutmover, mutprotocol):
    xml = '''
<ROSETTASCRIPTS>

  <TASKOPERATIONS>

    <InitializeFromCommandline name="init"/>
    <LimitAromaChi2 name="limchi2"/>
    <RestrictToRepacking name="repack_only"/>
    <DesignAround name="da" design_shell="8.0" resnums="'''+muttask+'''" repack_shell="6.0" allow_design="0"/>

  </TASKOPERATIONS>

 <SCOREFXNS>
    <ScoreFunction name="hard_rep" weights="ref2015">

    </ScoreFunction>
 </SCOREFXNS>


  <MOVERS>
  
  '''+mutmover+'''

  # Minimization of complex - no design allowed
  <TaskAwareMinMover name="min" bb="0" chi="1" jump="1" scorefxn="hard_rep" task_operations="init"/>
  # Packing of rotamers making sure no aromatic with chi2 of 90 degrees
  <PackRotamersMover name="repack" task_operations="init,limchi2,repack_only,da"/>
  
  <ParsedProtocol name="min_repack_min">
    <Add mover="min"/>
    <Add mover="repack"/>
    <Add mover="min"/>
  </ParsedProtocol>

  </MOVERS>


  <PROTOCOLS>

    # Insert mutations
    '''+mutprotocol+'''
    <Add mover_name="min_repack_min"/>

  </PROTOCOLS>

</ROSETTASCRIPTS>

'''
    return xml


def generate_mover(id, resid_chains,resname):
    return "<MutateResidue name=\"mr"+id+"\" target=\""+resid_chains+"\" new_res=\""+resname+"\"/>"

def generate_protocol(id):
    return "<Add mover_name=\"mr"+id+"\"/>"

# default chain
chain = "A"

# Grep variants to generate
variants = {}
with open("tot.fasta", 'r') as f:
    for line in f:
        if(line[0] == '>'):
            name_ = ""
            task_ = ""
            mut_ = ""
            protocol_ = ""
            tmpvar = line[1:].strip().split('_')
            for i in tmpvar[0:-1]:
                name_ += i+','
                id_ = i[1:-1]
                resid_chain_ = id_+chain
                resname_ = res_map[i[-1]]
                task_ += resid_chain_+","
                mut_ += generate_mover(id_, resid_chain_, resname_) + "\n"
                protocol_ += generate_protocol(id_) + "\n"
            variants[name_[:-1]] = (task_[:-1], mut_, protocol_)

for i in variants.keys():
    mt = variants[i][0]
    mm  = variants[i][1]
    mp =  variants[i][2]
    xml = generate_xml_file(mt, mm, mp)
    with open(i+'.xml', 'w') as f:
        f.write(xml)
    assert 1 == 0



