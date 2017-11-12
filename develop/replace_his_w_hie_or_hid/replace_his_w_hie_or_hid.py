import sys

def read_file(input_pdb):
    tmp_file = open(input_pdb,'r')
    file = tmp_file.readlines()
    tmp_file.close()
    return file


input_pdb = sys.argv[1]


pdb_file = read_file(input_pdb)

dict_his = {
    "28" : "HIE",
    "87" : "HIE",
    "107": "HIE",
    "110": "HID",
    "117": "HID",
    "118": "HID",
    "119": "HIE",
    "141": "HID"
}

new_file = open('new_file.pdb','w')

for line in pdb_file:

    if(line[17:19] == "HI"):

        for k,v in dict_his.iteritems():

            if k == line[23:26].strip():
                
                line = str(line[0:17])+v+str(line[20:])
    new_file.write(line)

