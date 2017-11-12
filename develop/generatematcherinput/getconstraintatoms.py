#!/usr/bin/env python

'''

@requires cstfile in list

@return number of constraints, atoms to be constrainted


'''
# Used to identify numbers of constraints
NCST = 'CST::BEGIN'
# Ligand atom names
ATM = 'atom_name:'
LEN_ATM = len(ATM)

class getconstraintatoms:


    def getligandatomnames(self,cstatoms,line):
        # Assuming first postion is name for match
        cstatoms.append(line.split()[0])
        return cstatoms



    def getnumberconstraints(self,cstfile):
        ncst = 0
        cstatoms = []
        tmp = open(cstfile,'r')
        for line in tmp:
            tmp_pos = line.find(NCST)
            tmp_pos2 = line.find(ATM)
            if(tmp_pos == 0):
                ncst = ncst +1
            elif(tmp_pos2 != -1):
                # We remove the matched string
                # as the atom names should be there
                pos_string = tmp_pos2+LEN_ATM
                self.getligandatomnames(cstatoms,line[pos_string:])
        tmp.close()
        
        return ncst,cstatoms
