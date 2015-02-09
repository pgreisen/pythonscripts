
class matcherposformat:
    '''


    '''
    
    # @requires integer how many geometrical constraints are necessary, name for position file
    # @return string with number of constraints, geometrical constraint, positions on scaffold
    
    def set_ncst_string(self,positions,posfile='new_posfile.pos'):
        numberofgeometrical = len(positions)
        str1 = 'N_CST '+str(numberofgeometrical)+'\n'
        posfile = open(posfile,'w')
        posfile.write(str1)
        for i in range(numberofgeometrical):
            tmp_line = str(i+1)+':'
            for j in positions[i]:
                tmp_line = tmp_line+' '+str(j)
            posfile.write(tmp_line+'\n')
    
    
