import os, shutil, subprocess


class get_native_protein_and_setup:

    global parameterfile
    parameterfile = "CDL.params"


    global native_pdb_id
    native_pdb_id = 1


    def get_pdb_id(self, pdbfile ):

        global native_pdb_id

        return  pdbfile.split('_')[native_pdb_id]



    def copy_native_pdb(self,pdbname):
        assert len( pdbname ) == 4

        pdbname = pdbname.lower()


        native_pdb = pdbname+'_nohet_1_relax.pdb'

        pdb_path = '/lab/shared/scaffolds/'+pdbname[1:3]+'/'+pdbname+'/'+native_pdb

        move_files = 'cp '+pdb_path+' .'

        subprocess.Popen(move_files,shell=True).wait()

        return native_pdb

    def run(self,native,design,parameterpath):

        exe = "sh /work/greisen/files/post_analysis_holo_protein/run.sh "+str(design)+" "+str(native)+" "+parameterpath+"/"+parameterfile+" >log&"

        subprocess.Popen(exe,shell=True).wait()


    def main(self):

        dummy = 1

        path = './'
        files = os.listdir(path)

        parameterpath = os.path.abspath('./')

        print parameterpath

        for pdbfile in files:

            if(pdbfile.endswith(".pdb")):

                native_pdbid = self.get_pdb_id(pdbfile)

                dest = str( dummy )

                os.mkdir( dest )

                shutil.move(pdbfile, dest )

                os.chdir( dest )

                # copy native pdb to directory
                native_pdb = self.copy_native_pdb( native_pdbid )

                # run 
                self.run(native_pdb, pdbfile, parameterpath )

                os.chdir('../')

                dummy += 1



if __name__ == "__main__":
   pos = get_native_protein_and_setup()
   pos.main()