tmpObj="new"
cmd.create( tmpObj, "all and name ca and ss s");
lrr = []
# works
#cmd.iterate("all and name ca and ss s","print resi"  )
cmd.iterate("all and name ca and ss s","lrr.append( resi )"  )

#cmd.iterate("all and name ca and ss s",print "distance resi and name ca, resn TSV and name O8"  )

print lrr

for i in lrr:
    print i
    print cmd.distance("dst","resi %s and name ca" %i, "resn tsv and name o8")


print lrr

