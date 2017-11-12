import matplotlib,os,csv,subprocess
# matplotlib.use('macosx')
from pylab import *
# from operator import itemgetter
from operator import *
from collections import OrderedDict

'''

The class can be tested by executing the following:

python write_cvs_file_docking_results.py logfile score.sc

The first input is a file that contains the scores from a docking run.

The second file is a rescore of the first with different weights.

plots the rmsd versus the energy in the second file

'''

class AnalysisOfRosettaDocking:


    def get_file(self,filename):
        tmpfile = open(filename,'r')
        fl = tmpfile.readlines()
        return fl

    def get_rmsd_ife_tag(self,file):
        energy = []
        rmsd = []
        tag = []
        i = 1
        for line in file:
            if i > 1 and len(line.split()) == 64:
                tmp1 = float(line.split()[51])
                tmp2 = float(line.split()[52])
                tmp3 = str(line.split()[0])
                energy.append(tmp1)
                rmsd.append(tmp2)
                tag.append(tmp3)
            i = 1+i

        return tag,rmsd,energy

    def get_top_solutions_from_docking(self,tag,energy,parameter):
        
        write_top_solutions_to_file = open("get_top_solutions.sh",'w')
        template_line = "~/git_rosetta/rosetta_master/Rosetta/main/source/bin/extract_atomtree_diffs.linuxgccrelease -database ~/git_rosetta/rosetta_master/Rosetta/main/database -extra_res_fa "+parameter+" -s silent.out -tags "

        top_solutions =  {}
        assert len(tag) == len(energy)
        
        for i in range( len(tag) ):
            top_solutions[tag[i]] = energy[i]

        #        sorted_x = sorted(top_solutions.iteritems(), key=itemgetter(1))
        d_sorted_by_value = OrderedDict(sorted(top_solutions.items(), key=lambda x: x[1]))
        dummy = 0
        for k,v in d_sorted_by_value.iteritems():
            if dummy <= 10:
                write_top_solutions_to_file.write(template_line+k+'\n')
            else:
                break
            dummy = dummy + 1
        write_top_solutions_to_file.close()
            


    def write_csv_file(self,tag,rmsd,energy):
        with open('docking_data.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["Tag","RMSD","IFE"])
            number_of_tags = len(tag)
            assert len(tag) == len(rmsd) == len(energy)
            i = 0
            for i in range(number_of_tags):
                spamwriter.writerow([tag[i],rmsd[i],energy[i]])
                i += 1

    def plot(self,rmsd,ife):
        fig = figure()
        axis = fig.add_subplot(111)
        axis.set_xlim(0,20)
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

    def main(self):

        inputfile = sys.argv[1]

        parameterfile = sys.argv[2]

        fl = self.get_file( inputfile )
        tag,rmsd,ife = self.get_rmsd_ife_tag( fl )
        self.write_csv_file(tag,rmsd,ife)
        self.plot(rmsd,ife)

        self.get_top_solutions_from_docking(tag,ife,parameterfile)


        # get ddg
        #inputfile2 = sys.argv[2]
        #fl2 = self.get_file( inputfile2 )

        #scoreterm = "ddg"

        #list_with_values = self.write_score_term_to_csv_file(fl2,scoreterm)

        #self.write_single_column_to_csv_file(list_with_values)

        #self.plot_rmsd_versus_energy(rmsd,list_with_values)

if __name__ == "__main__":
    run = AnalysisOfRosettaDocking()
    run.main()
