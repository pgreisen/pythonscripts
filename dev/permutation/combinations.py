def combinations(string):
     yield ''
     for i, d in enumerate(string):
         for comb in combinations(string[i+1:]):
             yield d + comb


a= combinations("1234")
for i in a:
     print i
