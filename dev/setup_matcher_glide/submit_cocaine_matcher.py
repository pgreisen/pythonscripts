import sys, shutil, os, subprocess

'''

Requires a directory with pdb files. Loops over all the pdbs and makes a directory
and submit the job.

Requires: User needs to set the path to the different parameters such as parameters, flags,
xml, database, and executables. 

and is simply executed like this:

python submit_glide.py


'''
# defining the xml-file
xml="enzdes_w_cst_cocaine.xml"
params = "CO2.fa.params"
cstfile = "CO2designDS.cst"



def write_xml_file_to_disk(cstfile_name):

    template='''Design using cst-files to keep certain interactions in the
beginning.
With talaris I am not using the soft-scorefunction at the moment.
<ROSETTASCRIPTS>
      <SCOREFXNS>
        <enzdes weights=talaris2013>
          <Reweight scoretype=res_type_constraint weight=1.0/>
        </enzdes>
        <soft weights=ligand_soft_rep>
          <Reweight scoretype=res_type_constraint weight=1.0/>
        </soft>
      </SCOREFXNS>
      <TASKOPERATIONS>
        <InitializeFromCommandline name=init/>
        # Design in the interface between ligand and protein - distance based on ligand
     <DetectProteinLigandInterface name=shelldesign design=1 cut1=0.0 cut2=8.0 cut3=10.0 cut4=12.0 catres_interface=0/>
        # No design allowed in the interface between protein and ligand. Distances based on ligand.
        <DetectProteinLigandInterface name=shellnodesign design=0 cut1=6.0 cut2=8.0 cut3=10.0 cut4=12.0 catres_interface=0/>

        <DetectProteinLigandInterface name=dpli_nodesign design=0 cut1=6.0 cut2=8.0 cut3=10.0 cut4=12.0 catres_interface=0/>

        # Linear access to memory
        <SetIGType name=linmem_ig lin_mem_ig=1/>
        # Only use Trp if native or present in the seq
        # Disallow if non-native residue
        <DisallowIfNonnative name=dis_allow resnum=0 disallow_aas=W/>
        # Limit aromatic residues with angle of 90 degrees
        <LimitAromaChi2 name=limchi2/>
        # Upweight the interface energy between ligand and protein
        <ProteinLigandInterfaceUpweighter name=up interface_weight=1.2/>
        # Make sure that the catalytic residues are not designed but allowed to repack
        <SetCatalyticResPackBehavior name=catres fix_catalytic_aa=0/>
        # Catalytic residues are not designed and not allowed to get repack
        <SetCatalyticResPackBehavior name=fixcat fix_catalytic_aa=1/>
          # Interface between protein and ligand with design
        # shells of design
        # cut1 is set to 0 such that the Ca-Cb vector is only designed
        <DetectProteinLigandInterface name=dpli cut1=0.0 cut2=8.0 cut3=10.0 cut4=12.0 design=1/>
        # No design only repacking of residues
        <RestrictToRepacking name=repack_only/>

      </TASKOPERATIONS>


      <MOVERS>

        # Minimization of complex - no design allowed
        <TaskAwareMinMover name=min bb=0 chi=1 jump=1 scorefxn=enzdes task_operations=shellnodesign,init/>
        # Packing of rotamers making sure no aromatic with chi2 of 90 degrees
        <PackRotamersMover name=repack task_operations=init,shellnodesign,limchi2/>

        <ParsedProtocol name=min_repack_min>
          <Add mover=min/>
          <Add mover=repack/>
          <Add mover=min/>
        </ParsedProtocol>
         # Add cst to the pose - make sure cst-file in same directory
        # or absolute path
        <AddOrRemoveMatchCsts name=cstadd cstfile="'''+cstfile_name+'''" cst_instruction=add_new/>
        # Remove cst
        <AddOrRemoveMatchCsts name=cstrem cstfile="'''+cstfile_name+'''" cst_instruction=remove/>
        # Add cst back on - remember to specify absolute path to cst
        # file or have a copy of file in directory
        <AddOrRemoveMatchCsts name=fincstadd cstfile="'''+cstfile_name+'''" cst_instruction=add_pregenerated/>
      </MOVERS>

     <FILTERS>
        # Interface energy between protein and ligand
        <LigInterfaceEnergy name=b_lie scorefxn=enzdes jump_number=1 confidence=1/>
        <Delta name=delta_ife upper=1 range=0 filter=b_lie/>
        # Ligand interface where the cst_score is NOT substrated hence
        # bad cst high score
        <LigInterfaceEnergy name=lig_cst  scorefxn=enzdes include_cstE=1 jump_number=1 energy_cutoff=-8.08/>
        # Ligand interface energy - make sure decent score
        <LigInterfaceEnergy name=lig_ife scorefxn=enzdes jump_number=1 energy_cutoff=-8.08/>
        # Delta filter for IFE with cst
        <Delta name=delta_ife_cst upper=1 range=1 filter=lig_cst/>
        # Computes the cst-score for the pose - the energy cutoff
        # should scale with number of csts i use 0.5 per cst
        # That might be a little to high
        <EnzScore name=cst_score  scorefxn=enzdes whole_pose=1 score_type=cstE energy_cutoff=1.0/>
        # chain=0 means compute for the complex - the threshold value
        # 0.6 is based on empirical observation from Nobu
        <PackStat name=pstat chain=0 threshold=0.6/>
        # Make sure AtomA closer to center of mass compared with AtomB
        <DiffAtomCenterOfMass name="com_ligand" AtomA="C6" AtomB="C12"/>
        # 1 is a little lenient
        RepackWithoutLigand name=rwl scorefxn=enzdes rms_threshold=1
        The calculations is run 10 times to get a semi-reasonable value, the threshold is
        set quite lenient.
        <Ddg name=ddg_complex scorefxn=enzdes threshold=-5 jump=1 repeats=10 relax_mover=min_repack_min/>
        <DSasa name=dsasa lower_threshold=0.90/>
        <ShapeComplementarity name=sc min_sc=0.65 jump=1/>
        # Make sure that the main interaction to the ligand
        # is place on seconadry structure elements
        <SecondaryStructureDensity name=ssd threshold=0.8 shell=8/>
        </FILTERS>

   <MOVERS>

        # Adding a bonus to the native sequence
        # Not used in this run
        <FavorSequenceProfile name=fsp scaling=prob use_starting=1 matrix=IDENTITY weight=1.5 scorefxns=enzdes/>

        # Design with hard minimization / soft repacking with the
        # following task operations
        # Fix catalytic residues (fixcat )
        # Interaction between prot/lig upweighted 1.5 ( up )
        # Trp not allowed unless it is native ( dis_allow )
        #
        # With a generic monte carlo using the cst energy
        <EnzRepackMinimize name=desmincst design=1 cst_opt=1 repack_only=0 scorefxn_minimize=enzdes scorefxn_repack=soft_rep minimize_rb=1 minimize_sc=1 minimize_bb=0 cycles=1 minimize_lig=1 min_in_stages=0 backrub=0 task_operations=init,dpli,limchi2,linmem_ig,up,dis_allow,fixcat/>
        # Only accept lower energies temperature=0
        <GenericMonteCarlo name=gmc_desmincst mover_name=desmincst filter_name=lig_cst trials=5 sample_type=low temperature=0 drift=0/>

        # Same as above but with hard repack instead of soft.
        <EnzRepackMinimize name=desmincst_hard design=1 cst_opt=1 repack_only=0 scorefxn_minimize=enzdes scorefxn_repack=enzdes minimize_rb=1 minimize_sc=1 minimize_bb=0 cycles=1 minimize_lig=1 min_in_stages=0 backrub=0 task_operations=init,dpli,limchi2,linmem_ig,up,dis_allow,fixcat/>
        # Only accept lower energies temperature=0
        <GenericMonteCarlo name=gmc_desmincst_hard mover_name=desmincst_hard filter_name=lig_cst trials=5 sample_type=low temperature=0 drift=0/>
                # Design with hard minimization / soft repacking with the
        # and no optimization of cst-values
        # following task operations
        # Fix catalytic residues but able to repacked (fixcat )
        # Interaction between prot/lig upweighted 1.5 ( up )
        # Trp not allowed unless it is native ( dis_allow )
        #
        # For the ligand interface energy with cst we now use a
        # a delta filter such that this score is not getting worse
        # than 1 REU
        #
        # With a generic monte carlo using the cst energy
        <EnzRepackMinimize name=desmin_repack_cst design=1 repack_only=0 scorefxn_minimize=enzdes scorefxn_repack=soft_rep minimize_rb=1 minimize_sc=1 minimize_bb=0 cycles=1 minimize_lig=1 min_in_stages=0 backrub=0 task_operations=init,dpli,limchi2,linmem_ig,up,dis_allow,catres/>
        # Only accept lower energies temperature=0
        <GenericMonteCarlo name=gmc_desmin_repack_cst mover_name=desmin_repack_cst filter_name=delta_ife_cst trials=5 sample_type=low temperature=0 drift=0/>
         # Same as above with hard repacking
        <EnzRepackMinimize name=desmin_repack_cst_hard design=1 repack_only=0 scorefxn_minimize=enzdes scorefxn_repack=enzdes minimize_rb=1 minimize_sc=1 minimize_bb=0 cycles=1 minimize_lig=1 min_in_stages=0 backrub=0 task_operations=init,dpli,limchi2,linmem_ig,up,dis_allow,catres/>
        # Only accept lower energies temperature=0
        <GenericMonteCarlo name=gmc_desmin_repack_cst_hard mover_name=desmin_repack_cst filter_name=delta_ife_cst trials=5 sample_type=low temperature=0 drift=0/>


        # Hard minimization / hard repack without design and no backbone
         # minimization
        <EnzRepackMinimize name=rpkmin repack_only=1 design=0 scorefxn_minimize=enzdes scorefxn_repack=enzdes minimize_rb=1 minimize_sc=1 minimize_bb=0 cycles=3 task_operations=init,dpli_nodesign,limchi2,linmem_ig,fixcat/>

        # Hard minimization / hard repack without design with backbone
        # minimization
        <EnzRepackMinimize name=min_w_bb repack_only=0 design=0 scorefxn_minimize=enzdes scorefxn_repack=enzdes minimize_rb=1 minimize_sc=1 minimize_bb=1 cycles=3 task_operations=init,dpli_nodesign,limchi2,linmem_ig,fixcat/>
     # Hard minimization / hard repack without design and no backbone
        # minimization
        <EnzRepackMinimize name=rpkmin_nocst repack_only=1 design=0 scorefxn_minimize=enzdes scorefxn_repack=enzdes minimize_rb=1 minimize_sc=1 minimize_bb=0 cycles=3 task_operations=init,dpli_nodesign,limchi2,linmem_ig/>

      </MOVERS>


      <PROTOCOLS>
              <Add filter=com_ligand/>
              <Add mover_name=cstadd/>
              <Add mover_name=gmc_desmincst/>
              <Add mover_name=gmc_desmincst_hard/>
              <Add mover_name=gmc_desmin_repack_cst/>
              <Add mover_name=gmc_desmin_repack_cst_hard/>
              <Add mover_name=rpkmin/>
              <Add mover_name=min_w_bb/>
              <Add filter_name=cst_score/>
            # Remove cst
              <Add mover_name=cstrem/>
              <Add mover_name=rpkmin_nocst/>


              # Only output poses which passes the filters
              <Add filter=lig_ife/>
              <Add filter=pstat/>
              <Add filter=dsasa/>
              <Add filter=sc/>
              <Add filter=ddg_complex/>
              Add filter=rwl

      </PROTOCOLS>

</ROSETTASCRIPTS>
    '''
    tmpfile = open("enzdes_w_cst_cocaine.xml", 'w')
    tmpfile.write(template)
    tmpfile.close()

    # cstfile = "CO2designDS.cst"



def write_wrapper(pdbname):
    wrapper = '''
    #!/bin/bash

    tar -xvf database.tgz

    ./rosetta_scripts.static.linuxgccrelease -overwrite @flags -parser:protocol '''+str(xml)+''' -out:prefix ${1} -extra_res_fa '''+params+''' -database database -s '''+str(pdbname)+''' > run_${1}.txt

    RETVAL=$?

    echo rosetta retval is $RETVAL
    if [ "$RETVAL" != "0" ]; then
    echo exiting with $RETVAL
    exit $RETVAL
    fi

    echo checking log file for "jobs considered"
    grep -q \'jobs considered\' run_${1}.txt

    RETVAL=$?
    echo grep retval is $RETVAL
    
    if [ "$RETVAL" != "0" ]; then
    echo NOT found, exiting with 1
    exit 1
    fi
    
    echo string found, exiting with 0
    exit 0
    '''
    return wrapper



def update_condor_script(pdbfile):
    MAINDIR = '/home/greisen/files/'
    PTH = os.path.abspath('./')
    print 'The path is ', PTH
    template = '''
    notification=Never
    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT
    on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
    
    transfer_input_files = /home/greisen/database.tgz, '''+str(PTH)+'/'+xml+''', '''+MAINDIR+'''flags, '''+MAINDIR+params+''', '''+str(PTH)+'''/'''+str(pdbfile)+''', '''+str(PTH)+'''/run_wrapper.sh, /home/greisen/files/'''+cstfile+''', /home/greisen/rosetta_source_bin_01_11_2013_upweight_hb/rosetta_scripts.static.linuxgccrelease

    Executable = run_wrapper.sh
    universe = vanilla
    copy_to_spool = false

    Error = err.$(Process)
    Output = out.$(Process)
    Log = condor_log.txt

    Arguments = $(Process)
    queue 1
    
    '''

    return template

def write_template_to_file(template,name='condor.submit'):
    cndr = open(name,'w')
    for line in template:
        cndr.write(line)

PTH = os.path.abspath('./')
# Get the files in the directory
pdbfiles = os.listdir(PTH)
i = 1
for pdbfile in pdbfiles:

    if(os.path.isfile(pdbfile) and str(pdbfile[-3:]) == 'pdb'):

        shutil.os.mkdir(str(i))

        shutil.move(pdbfile,str(i))

        os.chdir(str(i))
        write_xml_file_to_disk(cstfile)
        # Condor
        condor_template = update_condor_script(pdbfile)

        write_template_to_file(condor_template,name='condor.submit')
        
        # Wrapper
        wrapper_template = write_wrapper(pdbfile)

        write_template_to_file(wrapper_template,name='run_wrapper.sh')

        #subprocess.Popen('condor_submit condor.submit',shell=True)

        os.chdir('../')

        i = i +1
