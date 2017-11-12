import sys

def read_file(input_pdb):
    xvalue = []
    yvalue = []
    tmp_file = open(input_pdb,'r')
    for line in tmp_file:
        tmpline = line.split()
        xvalue.append( tmpline[0] )
        yvalue.append( float( tmpline[1] ) )

    tmp_file.close()
    return xvalue,yvalue

def write_to_files(x,y,inputname):
    with open("sorted_"+inputname, 'w') as f:
        for f1, f2 in zip(x, y):
            f.write(f1+"  "+str(f2)+"\n")

        
    
def main():
    # test()
    inputfile = sys.argv[1]
    x,y = read_file( inputfile )
    ##print y
    #    x, y = (list(tmp) for tmp in zip(*sorted(zip(x, y), key=lambda pair: pair[1])))
    x, y = (list(tmp) for tmp in zip(*sorted(zip(x, y), key=lambda pair: pair[1] , reverse=True) ))
    ##print y
    write_to_files(x,y,inputfile)
    
def test():
    list1 = [1, 2, 5, 4, 4, 3, 6]
    list2 = [3, 2, 1, 2, 1, 7, 8]
    print "Initial input"
    print "list1 ", list1
    print "list2 ", list2
    list1, list2 = (list(x) for x in zip(*sorted(zip(list1, list2), key=lambda pair: pair[0])))
    print list1
    # [1, 2, 3, 4, 4, 5, 6]
    print list2
    # [3, 2, 7, 2, 1, 1, 8]

if __name__ == "__main__":
   main()
