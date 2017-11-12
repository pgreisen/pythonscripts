import sys
def write_parameter_file(residuenumber,temperature):
    template = '''
100ps MD with res on protein
&cntrl
imin   = 0,
irest  = 0,
ntx    = 1,
ntb    = 1,
cut    = 10.0,
ntr    = 1,
ntc    = 2,
ntf    = 2,
tempi  = 0.0,
temp0  = '''+str(temperature)+''',
ntt    = 3,
gamma_ln = 1.0,
nstlim = 50000, dt = 0.002
ntpr = 5000, ntwx = 5000, ntwr = 5000
/
Keep protein fixed with weak restraints
50.0
RES 1 '''+str(residuenumber)+'''
END
END
'''
    return template

def get_number_of_residues(pdbfile):
    tmpfile = open(pdbfile,'r')
    max_number = 0
    for line in tmpfile:
        try:
            query_number = int( line[23:26].strip() )

            if( query_number > max_number):
                max_number = int(line[23:26])
        except:
            pass
    return max_number

def write_out_parameterfile(list):
    tmpfile = open("heat.in",'w')
    for line in list:
        tmpfile.write(line)


def main():
    pdbfile = sys.argv[1]
    temperature = sys.argv[2]
    max_number_residues = get_number_of_residues(pdbfile)
    print "The maximum number of residues in pdb is ", max_number_residues

        
    list_template = write_parameter_file(max_number_residues,temperature)

    write_out_parameterfile(list_template)


if __name__ == "__main__":
        main()
