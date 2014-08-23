import matplotlib,os,csv,subprocess
# matplotlib.use('macosx')
from AnalyseMutations import *
from pylab import *

'''

The class can be tested by executing the following:

python write_cvs_file_docking_results.py logfile score.sc

The first input is a file that contains the scores from a docking run.

NOT FINISHED

'''



class AnalysisOfRosettaDocking:


    def get_file(self,filename):
        tmpfile = open(filename,'r')
        fl = tmpfile.readlines()
        return fl

    # seven parameters
    #         self.write_csv_file( designnames ,number_of_mutations, sc, ife, pstat, rwl, dsasa, ddg_complex  )
    def write_csv_file(self,designname, number_of_mutations, sc, ife, pstat, rwl, dsasa, ddg_complex ):
        with open('design_protocol.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["Name","Number of mutations", "SC", "IFE", "Pstat", "rwl", "Dsasa","ddG"])
            
            number_of_tags = len(designname)
            print len(designname), len(number_of_mutations), len(sc), len(ife), len(pstat), len(dsasa), len(rwl)
            i = 0
            for i in range(number_of_tags):
                spamwriter.writerow([designname[i],number_of_mutations[i], sc[i], ife[i], pstat[i],rwl[i],dsasa[i],ddg_complex[i] ] )
                i += 1

    def plot(self,rmsd,ife):
        fig = figure()
        axis = fig.add_subplot(111)
        axis.set_xlim(0,9)
        axis.set_ylim(min(ife),0)
        plot(rmsd,ife,'.')
        # xlabel('rmsd')
        # ylabel('Interface energy')
        # savefig(title+'.png')
        show()



    def write_score_term_to_csv_file(self,filename,scoreterm):

        scores = []
        allow = False

        for line in filename:

            if( line[0:4] == "SCOR" and allow == True):
                scores.append( float (line.split()[place_of_score_term]) )

            if( line[0:4] == "SCOR"):
                dummy = -1
                tmp= line.split()
                for i in tmp:
                    dummy +=1
                    if i.strip() == scoreterm:
                        print i.strip()
                        place_of_score_term = dummy
                        allow = True
            else:
                continue
        return scores


    def write_single_column_to_csv_file(self,list_of_values):
        with open('docking_data_ddg.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["ddG"])
            i = 0
            for i in range( len( list_of_values ) ):
                spamwriter.writerow([list_of_values[i]])
                i += 1


    def plot_rmsd_versus_energy(self,rmsd,energy):
        print len(rmsd),len(energy)
        assert len(rmsd) == len(energy)
        fig = figure()
        axis = fig.add_subplot(111)
        axis.set_xlim(0,9)
        axis.set_ylim(min(energy),0)
        plot(rmsd,energy,'.')
        show()

    def get_number_of_mutations(self, nativepdb, designpdb):
        
        gnm = AnalyseMutations()

        a = gnm.get_seq_pdbfile(nativepdb)
        b = gnm.get_seq_pdbfile(designpdb)

        number_of_mutations = gnm.get_number_of_mutations(a,b)

        return number_of_mutations



    def get_scores_from_design(self,designfile):
        pdbfile = open(designfile,'r')
        sc = float(NaN)
        pstat = float(NaN)
        rwl = float(NaN)
        ife = float(NaN)
        dsasa = float(NaN)
        ddg_complex = float(NaN)
        for line in pdbfile:
            if(len(line) > 4 ):
                
                if( line[0:4] == "ATOM"):
                    continue
                elif( line[0:4] == "HETA"):
                    continue
                else:
                    tmp = line.split()

                    if(tmp[0] == "sc" or tmp[0] == "SC"):
                        sc = float( tmp[1] )

                    elif( tmp[0] == "pstat" ):
                        pstat = float( tmp[1] )
                    
                    elif( tmp[0] == "rwl" ):
                        rwl = float( tmp[1] )

                    elif( tmp[0] == "ddg_complex" ):
                        ddg_complex = float( tmp[1] )

                    elif( tmp[0] == "ife" or tmp[0] == "interfE"):
                        ife = float( tmp[1] )
        pdbfile.close()
        return sc, pstat, rwl, ife, dsasa, ddg_complex 

    def get_statistics(self,list):
        return round( mean(list) ,3) , round( sqrt(var(list)) , 3 )

    def main(self):
        path = './'
        # Loop over the directories
        dirs = os.listdir( path )

        sc = []
        pstat = []
        rwl = []
        ife = []
        dsasa = []
        ddg_complex = []
        designnames = []
        number_of_mutations = []
        
        for dir in dirs:
            if ( os.path.isdir ( dir )):
                os.chdir( dir )
                
                files = os.listdir( path )

                for file in files:
                    if( os.path.isfile( file ) and file[0] == "0" and file.endswith(".pdb") ):
                        design = file
                    elif ( os.path.isfile( file ) and file.endswith(".pdb") ):
                        native = file
                
                tmp_number_of_mutations = self.get_number_of_mutations(native, design)
                

                tmp_sc, tmp_pstat, tmp_rwl, tmp_ife, tmp_dsasa, tmp_ddg_complex = self.get_scores_from_design( design )
        
                sc.append( tmp_sc )
                pstat.append( tmp_pstat )
                rwl.append( tmp_rwl )
                dsasa.append( tmp_dsasa )
                ddg_complex.append( tmp_ddg_complex )
                number_of_mutations.append( tmp_number_of_mutations )                
                designnames.append( design )
                ife.append( tmp_ife )

                
                os.chdir( "../" )
        
        self.write_csv_file( designnames ,number_of_mutations, sc, ife, pstat, rwl, dsasa, ddg_complex  )

        self.write_statistics( number_of_mutations, sc, ife, pstat, rwl, dsasa, ddg_complex)


if __name__ == "__main__":
    run = AnalysisOfRosettaDocking()
    run.main()
