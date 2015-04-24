import os,sys, argparse
from numpy import *
from numpy import linalg as LA
from collections import defaultdict


class FastaManipulation:



    def __init__(self):
        self.fastafiles = {}
        self.native = ""
        self.design = ""
        self.first = True
        self.chain = 'A'

        self.mutations = {}


    def set_sequences(self, file):
        dummy = 1
        with open(file,'r') as f:
            for line in f:
                if( dummy == 1 ):
                    self.design = line
                elif( dummy == 2 ):
                    self.native = line
                dummy += 1
        print self.native
        print self.design

    def amino_acids(self, residue):
        atms = {
            'H' :'HIS',
            'E' :'GLU',
            'D' :'ASP',
            'C' :'CYS',
            'K' :'LYS',
            'T' :'THR',
            'S' :'SER',
            'R' :'ARG',
            'Q' :'GLN',
            'N' :'ASN',
            'Y' :'TYR',
            'M' :'MET',
            'P' :'PRO',
            'A' :'ALA',
            'V' :'VAL',
            'L' :'LEU',
            'I' :'ILE',
            'F' :'TYR',
            'W' :'TRP',
            'G' :'GLY',
            '-' : '-'
        }
        return atms[residue]


    def get_offset(self):
        """
        Get off set between sequences

        """
        return 0



    def get_diff_seq(self):
        substitution_string = []
        a_string = ""
        design_around_task = ""
        protocol = ""
        dummy = 0
        assert len( self.native ) == len( self.design )
        for i in range( len(self.native) ):
            if (self.native[i] != self.design[i]) :
                if(self.amino_acids(self.design[i]) == '-'):
                    continue
                a = "<MutateResidue name=mr"+str(dummy)+" target="+str(i)+"A new_res="+self.amino_acids(self.design[i] )+"/>"
                substitution_string.append(a)
                design_around_task = design_around_task+str(i)+self.chain+","
                protocol += "<Add mover_name=mr"+str(dummy)+"/>\n"
                dummy +=1


        design_around_task = design_around_task[0:-1]
        b = "<DesignAround name=da design_shell=8.0 resnums="+design_around_task+" repack_shell=6.0 allow_design=0/>"

        for i in substitution_string:
            a_string = a_string+i+"\n"

        return a_string, b, protocol


    def generate_xml(self,a,b,c):
        template = '''
        <ROSETTASCRIPTS>

  <TASKOPERATIONS>

    <InitializeFromCommandline name=init/>
    <LimitAromaChi2 name=limchi2/>
    <RestrictToRepacking name=repack_only/>
    '''+b+'''

  </TASKOPERATIONS>

  <SCOREFXNS>

    <enzdes weights=talaris2013.wts/>

  </SCOREFXNS>

  <MOVERS>

    '''+a+'''

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

    # Insert mutations
    '''+c+'''

    <Add mover_name=min_repack_min/>

  </PROTOCOLS>

</ROSETTASCRIPTS>

        '''
        return template



    def write_to_file(self, template):
        with open("mutate_residues.xml",'w') as f:
            for line in template:
                f.write(line)




    def main(self):

        parser = argparse.ArgumentParser(description="Manipulate the fasta-sequences")
        # get the initial rosetta design as input
        parser.add_argument("-f", dest="fastafile", help="An aligned fasta file ")

        input_variables = parser.parse_args()

        self.set_sequences( input_variables.fastafile )

        a,b,c = self.get_diff_seq()

        template = self.generate_xml( a,b,c)
        self.write_to_file( template )


if __name__ == "__main__":
    run = FastaManipulation()
    run.main()
