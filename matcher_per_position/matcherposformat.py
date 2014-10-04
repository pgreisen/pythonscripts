
class matcherposformat:
    '''


    '''
    
    # @requires integer how many geometrical constraints are necessary, name for position file
    # @return string with number of constraints, geometrical constraint, positions on scaffold
    
    def set_ncst_string(self,positions,ncst=0,posfile='new_posfile.pos'):
        numberofgeometrical = len(positions)
        if(ncst == 0):
            str1 = 'N_CST '+str(numberofgeometrical)+'\n'
        else:
            str1  = 'N_CST '+str(ncst)+'\n'

        posfile = open(posfile,'w')
        posfile.write(str1)
        dummy = 0
        for i in range(numberofgeometrical):
            tmp_line = str(i+1)+':'
            for j in positions[i]:
                tmp_line = tmp_line+' '+str(j)
            posfile.write(tmp_line+'\n')
            dummy += 1

        for k in range(ncst-1):
            dummy += 1
            posfile.write(str(dummy)+': ALL\n')
