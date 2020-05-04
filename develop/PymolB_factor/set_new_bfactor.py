

#cat 1XQ8.pdb |grep CA

# load the protein
cmd.load("protA.pdb")

# open the file of new values (just 1 column of numbers, one for each alpha carbon)
inFile = open("newBFactors", 'r')

# create the global, stored array
stored.newB = []

# read the new B factors from file
for line in inFile.readlines(): stored.newB.append( float(line) )

# close the input file
inFile.close()

# clear out the old B Factors
alter protA, b=0.0

# update the B Factors with new properties
alter protA and n. CA, b=stored.newB.pop(0)

# color the protein based on the new B Factors of the alpha carbons
cmd.spectrum("b", "protA and n. CA")
