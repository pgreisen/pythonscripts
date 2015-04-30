
class GeneratePSSM:


    def __init__(self):
        self.fastasequence = []
        self.aa = ['A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V']
        self.insert_text = "\n PSSM for design\n"
        self.padding_4 = "    "
        self.padding_5 = "   "
        self.padding_6 = "  "
        self.padding_7 = " "
        self.padding = "  "
        #self.padding_new = "    "
        self.padding_new = "  "
        self.number_of_amino_acids = 20
        self.beneficial_score  = 3
        self.g1_g10_score = 2
        self.neutral_score = 1
        self.deleterious_score = -1

        # off set between sequences +1 due to counting in python
        self.offset = 35

        # beneficial substitutions
        self.beneficial = {
            'A_203' : ['L'],
            'G_272' : ['E'],
            'Y_309' : ['W']
            }

        self.g1_g10 = {
            'C_59' : ['V','F','M'],
            'S_61' : ['G'],
            'K_77' : ['A','T'],
            'A_80' : ['V'],
            'V_99' : ['I'],
            'I_106' : ['A','V','C'],
            'S_111' : ['R'],
            'L_113' : ['M'],
    'L_130' : ['M'],
    'W_131' : ['H','F'],
    'E_132' : ['A','D','E'],
    'L_136' : ['E','P'],
    'M_138' : ['F'],
    'T_147' : ['Y'],
    'A_171' : ['S'],
    'T_172' : ['I'],
    'T_173' : ['N','Q','K','W','Y'],
    'T_177': ['N'],
    'F_179' : ['M'],
    'R_185' : ['K'],
    'A_203' : ['L','F','V','W','I','M','S','F','D','R','N','H'],
    'A_204' : ['G'],
    'G_208' : ['D'],
    'I_228' : ['V'],
    'D_233' : ['N','M','L','R','V','G','A','S','T'],
    'H_254' : ['N','Q','R','A','G','S','T','D'],
    'I_255' : ['M','L'],
    'P_256' : ['G'],
    'H_257' : ['Y'],
    'S_258' : ['N','E'],
    'S_267' : ['M','R'],
    'A_270' : ['E','T','D','I','M','S'],
    'L_271' : ['Y','P','F','V'],
    'I_274' : ['S','N'],
    'W_302' : ['V'],
    'L_303' : ['T','V'],
    'F_306' : ['H','I'],
    'S_307' : ['D'],
    'S_308' : ['M','G'],
    'Y_309' : ['S','W'],
    'V_310' : ['T','A','E','L'],
    'S_319' : ['R'],
    'P_342' : ['S'],
    'T_350' : ['M','L','F']
}


        self.deleterious = {
    'C_59' : ['T'],
    'G_60' : ['L','E'],
    'F_72' : ['E'],
    'F_73' : ['L','W'],
    'D_100' : ['E'],
    'V_101' : ['C'],
    'S_102' : ['L'],
    'I_106' :['M'],
    'W_131' : ['L','V', 'Y'],
    'E_132' : ['A'],
    'D_133' : ['K'],
    'M_138' : ['V'],
    'F_179' : ['D'],
    'Q_212' : ['E'],
    'D_233' : ['I'],
    'H_254' : ['D'],
    'I_255' : ['S'],
    'H_257' : ['Y','W'],
    'N_300' : ['H','Q'],
    'L_303' : ['F'],
    'S_307' : ['D'],
    'S_308' : ['E','N','D'],
    'Y_309' : ['R'],
    'M_317' : ['K']
}
        self.neutral = {
    'I_106' : ['L','V','A','C'],
    'L_271' : ['F','D','C', 'V', 'L', 'M', 'N', 'S'],
    'G_273' : ['L', 'D', 'T', 'V', 'S'],
    'F_306' : ['W'],
    'S_308' : ['T'],
    'M_317' : ['L']
}



    def set_fasta_sequence(self,fastafile):
        """
        Takes a fasta-file as input and sets the fasta sequence in the constructor.
        @requires a fasta file
        @sets the class variable self.fastasequence

        """
        with open(fastafile,'r') as f:
            for line in f:
                if(line[0] != '>'):
                    tmp = list( line.rstrip() )
                    for aa in tmp:
                        self.fastasequence.append(aa)


    def get_scores(self,position):
        """
        here the different scores are inserted into a string
        @requires position
        @returns string with different scores
        """
        dummy_array = [0]*20
        # check beneficial
        for key in self.beneficial:
            aa,pos = key.split('_')
            if( pos == str( position + self.offset) ):
                print pos, position - self.offset
                # # list.index('a')
                # print self.aa.index(aa), pos,aa
                # loop over list and insert values
                for value in self.beneficial[key]:
                    dummy_array[ self.aa.index(  value )  ] = self.beneficial_score
        for key in self.g1_g10:
            aa,pos = key.split('_')
            if( pos == str( position + self.offset) ):
                for value in self.g1_g10[key]:
                    dummy_array[ self.aa.index(  value )  ] = self.g1_g10_score
        # self.deleterious
        for key in self.deleterious:
            aa,pos = key.split('_')
            if( pos == str( position + self.offset) ):
                for value in self.deleterious[key]:
                    dummy_array[ self.aa.index(  value )  ] = self.deleterious_score
        for key in self.neutral:
            aa,pos = key.split('_')
            if( pos == str( position + self.offset) ):
                for value in self.neutral[key]:
                    dummy_array[ self.aa.index(  value )  ] = self.neutral_score

        return dummy_array


    def write_to_file(self):

        with open("pssm.dat",'w') as f:
            f.write( self.insert_text )
            f.write("           ")
            for aa in self.aa:
                f.write(aa+self.padding)
            f.write("\n")
            # loop over fasta sequence
            for position in range( len(self.fastasequence) ):

                pos_aa = str(position+1)+' '+self.fastasequence[position]

                if( len(pos_aa) == 3):
                    f.write(self.padding_4+pos_aa+self.padding_new)

                elif( len(pos_aa) == 4 ):
                    f.write(self.padding_5+pos_aa+self.padding_new)

                elif( len(pos_aa) == 5 ):
                    f.write(self.padding_6+pos_aa+self.padding_new)

                elif( len(pos_aa) == 6 ):
                    f.write(self.padding_7+pos_aa+self.padding_new)
                # print position
                values = self.get_scores(position)

                assert len( values ) == self.number_of_amino_acids

                # insert values into martix
                for i in values:
                    if(len(str(i)) == 1):
                        f.write('  '+str(i))

                    elif(len(str(i)) == 2):

                        f.write(' '+str(i))


                f.write('\n')


gs = GeneratePSSM()
gs.set_fasta_sequence("A53.fasta")
gs.write_to_file()

##print gs.fastasequence, len(gs.fastasequence)

# list.index('a')


# beneficial substitutions
beneficial = {
    'A_203' : ['L'],
    'G_272' : ['E'],
    'Y_309' : ['W']
}

g1_g10 = {
    'C_59' : ['V','F','M'],
    'S_61' : ['G'],
    'K_77' : ['A','T'],
    'A_80' : ['V'],
    'V_99' : ['I'],
    'I_106' : ['A','V','C'],
    'S_111' : ['R'],
    'L_113' : ['M'],
    'L_130' : ['M'],
    'W_131' : ['H','F'],
    'E_132' : ['A','D','E'],
    'L_136' : ['E','P'],
    'M_138' : ['F'],
    'T_147' : ['Y'],
    'A_171' : ['S'],
    'T_172' : ['I'],
    'T_173' : ['N','Q','K','W','Y'],
    'T_177': ['N'],
    'F_179' : ['M'],
    'R_185' : ['K'],
    'A_203' : ['L','F','V','W','I','M','S','F','D','R','N','H'],
    'A_204' : ['G'],
    'G_208' : ['D'],
    'I_228' : ['V'],
    'D_233' : ['N','M','L','R','V','G','A','S','T'],
    'H_254' : ['N','Q','R','A','G','S','T','D'],
    'I_255' : ['M','L'],
    'P_256' : ['G'],
    'H_257' : ['Y'],
    'S_258' : ['N','E'],
    'S_267' : ['M','R'],
    'A_270' : ['E','T','D','I','M','S'],
    'L_271' : ['Y','P','F','V'],
    'I_274' : ['S','N'],
    'W_302' : ['V'],
    'L_303' : ['T','V'],
    'F_306' : ['H','I'],
    'S_307' : ['D'],
    'S_308' : ['M','G'],
    'Y_309' : ['S','W'],
    'V_310' : ['T','A','E','L'],
    'S_319' : ['R'],
    'P_342' : ['S'],
    'T_350' : ['M','L','F']
}


neutral = {
    'I_106' : ['L','V','A','C'],
    'L_271' : ['F','D','C', 'V', 'L', 'M', 'N', 'S'],
    'G_273' : ['L', 'D', 'T', 'V', 'S'],
    'F_306' : ['W'],
    'S_308' : ['T'],
    'M_317' : ['L']
}



deleterious = {
    'C_59' : ['T'],
    'G_60' : ['L','E'],
    'F_72' : ['E'],
    'F_73' : ['L','W'],
    'D_100' : ['E'],
    'V_101' : ['C'],
    'S_102' : ['L'],
    'I_106' :['M'],
    'W_131' : ['L','V', 'Y'],
    'E_132' : ['A'],
    'D_133' : ['K'],
    'M_138' : ['V'],
    'F_179' : ['D'],
    'Q_212' : ['E'],
    'D_233' : ['I'],
    'H_254' : ['D'],
    'I_255' : ['S'],
    'H_257' : ['Y','W'],
    'N_300' : ['H','Q'],
    'L_303' : ['F'],
    'S_307' : ['D'],
    'S_308' : ['E','N','D'],
    'Y_309' : ['R'],
    'M_317' : ['K']
}

beneficial_score  = 3
g1_g10_score = 2
neutral_score = 1
deleterious_score = -1

dd_dec_seq = ['10','11','12','13','14','15','16','17','18','19']

dec_seq = ['1','2','3','4','5','6','7','8','9']
clean_matrix = []

with open("pte.pssm",'r') as f:
    for line in f:
        tmp = line.split()
        # print tmp

        if ( len(tmp) > 10):
            # replace characters
            # str.replace
            tmp_insert = line[10:70]

            for dd_dec in dd_dec_seq:
                tmp_insert = tmp_insert.replace(dd_dec,' 0',20)

            for dec in dec_seq:
                tmp_insert = tmp_insert.replace(dec,'0',20)

            new_line = line[0:10]+tmp_insert+line[71:]
            clean_matrix.append( new_line )
        else:
            clean_matrix.append( line )

with open("new_pte.pssm",'w') as f:
    for line in clean_matrix:
        f.write(line)

