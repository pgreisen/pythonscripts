def combinations(string):
     yield ''
     # (0, thing[0]), (1, thing[1]), (2, thing[2]), and so forth
     for i, d in enumerate(string):
          print "i,d",i,d
          for comb in combinations(string[i+1:]):
               yield d + comb


a= combinations("1234")
for i in a:
     print i
