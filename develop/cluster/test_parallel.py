import itertools
import multiprocessing

#Generate values for each parameter
a = range(3)
b = range(3)
c = range(3)
d = range(3)

#Generate a list of tuples where each tuple is a combination of parameters.
#The list will contain all possible combinations of parameters.
paramlist = list(itertools.product(a,b,c,d))

#A function which will process a tuple of parameters
def func(params):
  a = params[0]
  b = params[1]
  c = params[2]
  d = params[3]
  return (a*b*c*d,"Helloe",a,b,c,d)

#Generate processes equal to the number of cores
pool = multiprocessing.Pool()

#Distribute the parameter sets evenly across the cores
res  = pool.map(func,paramlist)

print(res)
