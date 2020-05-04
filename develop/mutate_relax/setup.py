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

 <FILTERS>

    <Ddg name="ddg" scorefxn="hard_rep" jump="1" repack="1" repeats="5" threshold="-7.5" confidence="1"/>
    <Sasa name="sasa" threshold="500" confidence="1"/>
    <BuriedUnsatHbonds name="buriedUnsatBonds" scorefxn="hard_rep" jump_number="1" cutoff="9" confidence="0"/>
    <ShapeComplementarity name="shapeComplementarity" jump="1" verbose="1" min_sc="0.5" write_int_area="0" confidence="1"/>
    <PackStat name="packstat" repeats="5" threshold="0.59" confidence="0"/>
    <ScoreType name="hpatch" scorefxn="hard_rep" score_type="hpatch" threshold="35" confidence="0"/>
    <ScoreType name="lr_elec" scorefxn="hard_rep" score_type="fa_elec" threshold="1200" confidence="0"/>
    <ScoreType name="total_score" scorefxn="hard_rep" score_type="total_score" threshold="0" confidence="0"/>
           
  </FILTERS>



  <PROTOCOLS>

    # Insert mutations
    '''+mutprotocol+'''
    <Add mover_name="min_repack_min"/>
    <Add filter="shapeComplementarity"/>
    <Add filter="ddg"/>
    <Add filter="total_score"/>
    <Add filter="sasa"/>
    <Add filter="buriedUnsatBonds"/>
    <Add filter="packstat"/>
    <Add filter="hpatch"/>
    <Add filter="lr_elec"/>

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
            # _denovo.pdb
            tmpvar1 = line[1:].strip().split(',')
            tmpvar = tmpvar1[0].split("#")
            # insert pdb name into xml-name
            name_ += tmpvar1[-1].replace("_denovo.pdb","")+"_"
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
    filename = i.replace(',','_')
    with open(filename+'.xml', 'w') as f:
        f.write(xml)



