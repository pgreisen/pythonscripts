#!/usr/bin/python
from optparse import OptionParser
from math import sqrt
from numpy import *
import os,shutil



class InsertMutation:

    def __init__(self):
        self.filename = "",
        self.off_seq = 0,
        self.threshold = -1



    def get_template(self,residuenumber,newaminoacid):
        print str(int(residuenumber)+self.off_seq)
        template = '''<ROSETTASCRIPTS>
  <TASKOPERATIONS>
    <InitializeFromCommandline name=init/>
    <LimitAromaChi2 name=limchi2/>
    <RestrictToRepacking name=repack_only/>
    <DesignAround name=da resnums="%%resnum%%" repack_shell=6.0 design_shell=6.0 allow_design=0/>
  </TASKOPERATIONS>

  <SCOREFXNS>
    <enzdes weights=talaris2013.wts/>
  </SCOREFXNS>

  <MOVERS>
    <MutateResidue name=mr1 target='''+str(int(residuenumber)+self.off_seq)+'''A new_res='''+self.get_single_residue(newaminoacid)+''' />
    # Minimization of complex - no design allowed
    <TaskAwareMinMover name=min bb=0 chi=1 jump=1 scorefxn=enzdes task_operations=init/>
    # Packing of rotamers making sure no aromatic with chi2 of 90 degrees
    <PackRotamersMover name=repack task_operations=init,limchi2,repack_only,da/>

    <ParsedProtocol name=min_repack_min>
      <Add mover=min/>
      <Add mover=repack/>
      <Add mover=min/>
    </ParsedProtocol>

  </MOVERS>

  <PROTOCOLS>
    <Add mover_name=mr1/>
    <Add mover_name=min_repack_min/>
  </PROTOCOLS>

</ROSETTASCRIPTS>

        '''
        return template

    def write_to_file(self,template,name):
        with open(name,'w') as f:
            for line in template:
                f.write(line)


    # Map of amino acids necessary for
    # Rosetta internal naming
    # Requires residue name
    # Return rosetta 1-letter code
    def get_single_residue(self,resn):

        res_map = {
            'C' : 'CYS',
            'H' : 'HIS',
            'D' : 'ASP',
            'E' : 'GLU',
            'K' : 'LYS',
            'S' : 'SER',
            'R' : 'ARG',
            'T' : 'THR',
            'N' : 'ASN',
            'Q' : 'GLN',
            'Y' : 'TYR',
            'M' : 'MET',
            'G' : 'GLY',
            'W' : 'TRP',
            'P' : 'PRO',
            'A' : 'ALA',
            'V' : 'VAL',
            'I' : 'ILE',
            'L' : 'LEU',
            'F' : 'PHE'
            }
        return res_map[resn]



    def main(self):

        parser = OptionParser()
        parser.add_option('-f', dest='filename', help='')
        parser.add_option('-o',dest='off_seq',help='',type=int)

        (options, args) = parser.parse_args()
        self.filename = options.filename
        self.off_seq = options.off_seq

        with open(self.filename,'r') as f:
            for line in f:
                tmpline = line.split(',')
                # print tmpline
                if( float( tmpline[3]) >= self.threshold ):
                    if ( tmpline[0][-1] == '*' ):
                        continue
                    else:
                        tmp_template = self.get_template(tmpline[0][1:-1], tmpline[0][-1])

                        tmp_name = tmpline[0][0]+str( int(tmpline[0][1:-1]) +self.off_seq ) +tmpline[0][-1]

                        self.write_to_file(tmp_template, "mutate_residue_"+tmp_name+".xml")


if __name__ == "__main__":
    run = InsertMutation()
    run.main()
