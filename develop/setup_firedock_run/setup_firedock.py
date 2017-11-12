import os,sys


class setup_firedock:


    def __init__(self):
        pass


    def get_position_file(self,inputfile):
        pass


    def write_flexible_ligands(self,list_of_position):

        with open("rec_flex_residues.txt","w") as position_file:

        for i in list:
            position_file.write(str(i) + "A\n")

        position_file.close()


    @property
    def get_firedock_template_file(self,proteinfilename,ligandfilename):

        template = '''
    ### I/O
    receptorPDBFileName '''+proteinfilename+'''
    ligandPDBFileName '''+ligandfilename+'''
    ## reference for rmsd calculations
    # templateLigandPDBFileName 1ACB_l.pdb.CHB.pdb
    # transformations for refinement
    #transFileName 1ACB.zdock.trans
    ## libraries files
    rotamerLibFile  /work/greisen/programs/FireDock/lib/bbdep02.May.sortlib
    protLib /work/greisen/programs/FireDock/lib/chem.lib
    # pdbConventionFile /specific/a/netapp2/vol/private/mol/bluvsh/FireDock123/lib/Names.CHARMM.db
    ## output file
    energiesOutFileName firedock.result.out

    ### Input Options
    ## these residues will be fixed unless appear also in receptorFlexibleResiduesFile
    # receptorFixedResiduesFile rec_fixed_residues.txt
    # ligandFixedResiduesFile lig_fixed_residues.txt
    ## defines for which chains to build surface residues rotamers. If empty than will build for all chains
    ## For bound receptor, you should uncomment it
    # flexibleReceptorChains @
    ## For bound ligand, you should uncomment it
    # flexibleLigandChains @
    ## Residues specified here will be flexible
    ## The residue can be flexible if its chain was specified in flexibleReceptorChains or if flexibleReceptorChains is empty
    ## PG notes file generated above
    receptorFlexibleResiduesFile rec_flex_residues.txt
    # ligandFlexibleResiduesFile lig_flex_residues.txt

    ### Output Options
    # to output refined complexes
    printRefinedComplexes 0
    ## 1 - only energy caclulaltion is performed without refinement (run only with flag toAMPL)
    onlyEnergyCalculation 0

    ### side-chain optimization
    # 1 - only clashing residues are flexible, 0 - all residues are flexible
    receptorOnlyClashesMovable 1
    ligandOnlyClashesMovable 1
    # 0 - small rotamer set, 1 - extended rotamer set
    extraRotamers 0

    ### rigid-body optimization
    # num of MC cycles for rigid-body minimization (if 0 - no RBM)
    rigidBodyMinimizationCycles 50

    ### weights for energy score (Others)
    attrVdWWeight   1.5
    repVdWWeight    0.8
    ACEWeight       1.6
    attrElWeight    0.21
    repElWeight     0.21
    l_attrElWeight  0.0
    l_repElWeight   0.69
    HBWeight        1.2
    pipiWeight      1.0
    catpiWeight     0.7
    aliphWeight     2.5
    insidenessWeight 1.55
    confProbWeight  0.0
    radiiScaling 0.8

    '''


        return template






    def main(self):

        inputfile = sys.argv[1]
        # converts rosetta position into FireDockPositions
        initial_position_file = get_position_file(inputfile)





if __name__ == "__main__":
   setup = setup_firedock()
   setup_firedock.main()