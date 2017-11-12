import subprocess,sys
from xml.dom import minidom


class GetPfam:

    def get_pfam_id(self,pdbid):
        # curl is a bash command equivalent to wget
        exe = "curl http://pdb.org/pdb/rest/hmmer?structureId="+str(pdbid)+" >pfam_tmp.xml "
        subprocess.Popen(exe,shell=True).wait()


    def __init__(self):
        # test if programs are present
        pass


    def read_pfam_xml(self):
        xmltest = minidom.parse("pfam_tmp.xml")
        # get the outer tag
        itemlist = xmltest.getElementsByTagName('pfamHit')
        for s in itemlist :
            pfamname =  s.attributes['pfamName'].value
        tmp_pfam = open("pfam_id.dat","w")
        tmp_pfam.write("pfam name: "+pfamname)
        tmp_pfam.close()

    def main(self):
        pdbid = sys.argv[1]
        print "PDB ID is ",pdbid
        self.get_pfam_id(pdbid)
        self.read_pfam_xml()

if __name__ == "__main__":
   gp = GetPfam()
   gp.main()
