import os, shutil


class SetupFFTdockProtocol():

    def __init__(self):
        self.side1 = "A"
        self.side2 = "B"
        self.rot_step = 60
        self.nmodels = 10
        self.trans_step = 1.0
        self.electrostatics = "false"


    def get_xml_file(self,side1,side2,rot_step,nmodels,trans_step,electrostatics):
        template="""
<dock_design>
<TASKOPERATIONS>
</TASKOPERATIONS>
<SCOREFXNS>
</SCOREFXNS>
<FILTERS>
</FILTERS>
<MOVERS>
    <FFTdockallatoms name=dock side1="""+str(side1)+""" side2="""+str(side2)+""" rot_step="""+str(rot_step)+""" nmodels="""+str(nmodels)+""" electrostatics="""+str(electrostatics)+""" density=true scale_radii=1 trans_step="""+str(trans_step)+"""/>
</MOVERS>
<APPLY_TO_POSE>
</APPLY_TO_POSE>
<PROTOCOLS>
    <Add mover=dock/>
</PROTOCOLS>
</dock_design>

        """
        return template


    def get_chain_id_and_pdb_length(self,pdbfile):
        chain = ""
        length_of_protein = 0
        with open(pdbfile,'r') as f:
            for line in f:
                if( line[0:4] == "ATOM" and length_of_protein == 0 ):
                    chain = line[21:22]
                if(line[0:4] == "ATOM" and line[13:15] == "CA" ):
                    length_of_protein += 1
        return chain, length_of_protein


    def set_chain_id(self, pdbfile,newchain):
        tmppdbfile = open("tmp.pdb",'w')
        with open(pdbfile, 'r') as f:
            for line in f:
                if line[0:4] == "ATOM":
                    newline = str(line[0:21])+str(newchain)+str(line[22:])
                    tmppdbfile.write(newline)
        tmppdbfile.close()
        shutil.move("tmp.pdb", pdbfile)

    def write_xml( self, xml_template ):
        with open("dock.xml",'w') as f:
            for line in xml_template:
                f.write(line)


    def main(self):
        path = './'
        dirs = os.listdir( path )
        # loop over all the directories
        for dir in dirs:
            if ( os.path.isdir ( dir ) ):

                # chdir
                os.chdir( dir )
                pdbfiles = os.listdir( path )
                contains_pdbs = False
                pdbs = []
                chain = {}
                length_of_protein = {}
                # loop over all the pdbfiles
                for pdbfile in pdbfiles:
                    if( os.path.isfile(pdbfile) and pdbfile.endswith("pdb")):
                        pdbs.append( pdbfile )
                        tmp = self.get_chain_id_and_pdb_length(pdbfile)
                        print "THE CHAIN OF THE PROTEIN IS DETERMINED TO BE ", tmp
                        length_of_protein[ pdbfile ] = tmp[1]
                        chain[pdbfile] = tmp[0]
                        contains_pdbs = True

                if( contains_pdbs ):
                    if( len(pdbs) > 2 ):
                        print "More than two pdb files", len(pdbfiles), dir
                    if( len(pdbs) == 2):

                        if(chain[pdbs[0]] == chain[pdbs[1]]):
                            print "These pdbfiles contain the chain ID", chain[pdbs[0]], chain[pdbs[1]]

                            if(chain[pdbs[0]] == "A"):
                                self.set_chain_id( pdbs[0],'B')
                                chain[ pdbs[0]] = 'B'
                            else:
                                self.set_chain_id( pdbs[0],'A')
                                chain[ pdbs[0]] = 'A'

                    if( length_of_protein[pdbs[0] ] > length_of_protein[pdbs[1] ] ):
                        xml_template = self.get_xml_file(chain[ pdbs[0]], chain[ pdbs[1]], self.rot_step, self.nmodels, self.trans_step, self.electrostatics)

                    else:
                        xml_template = self.get_xml_file(chain[ pdbs[1]], chain[ pdbs[0]], self.rot_step, self.nmodels, self.trans_step, self.electrostatics)
                    self.write_xml( xml_template )
                os.chdir("../")


if __name__ == '__main__':
    run = SetupFFTdockProtocol()
    run.main()
