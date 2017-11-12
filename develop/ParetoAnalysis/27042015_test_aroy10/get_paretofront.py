from pylab import *
import sys,os,  argparse
from collections import defaultdict

'''
Method to take two equally-sized lists and return just the elements which lie 
on the Pareto frontier, sorted into order.
Default behaviour is to find the maximum for both X and Y, but the option is
available to specify maxX = False or maxY = False to find the minimum for either
or both of the parameters.


def pareto_frontier(Xs, Ys, maxX = True, maxY = True):
# Sort the list in either ascending or descending order of X
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
# Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]    
# Loop through the sorted list
    for pair in myList[1:]:
        if maxY: 
            if pair[1] >= p_front[-1][1]: # Look for higher values of Y
                p_front.append(pair) # and add them to the Pareto frontier
        else:
            if pair[1] <= p_front[-1][1]: # Look for lower values of Y
                p_front.append(pair) # and add them to the Pareto frontier
# Turn resulting pairs back into a list of Xs and Ys
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY


import matplotlib.pyplot as plt   
Xs, Ys = # get your data from somewhere to go here
# Find lowest values for cost and highest for savings
p_front = pareto_frontier(Xs, Ys, maxX = False, maxY = True) 
# Plot a scatter graph of all results
plt.scatter(Xs, Ys)
# Then plot the Pareto frontier on top
plt.plot(p_front[0], p_front[1])
plt.show()



'''


class get_union:



    def __init__(self):
        self.mutations_positions = defaultdict(list)
        self.maxX = True
        self.maxY = True
        self.pdbs_on_front = []
        self.iterations = 3
        self.px = []
        self.py = []

    '''
    def pareto_frontier_multi(myArray):
        # Sort on first dimension
        myArray = myArray[myArray[:,0].argsort()]
        # Add first row to pareto_frontier
        pareto_frontier = myArray[0:1,:]
        # Test next row against the last row in pareto_frontier
        for row in myArray[1:,:]:
            if sum([row[x] >= pareto_frontier[-1][x]
                for x in range(len(row))]) == len(row):
                # If it is better on all features add the row to pareto_frontier
                pareto_frontier = np.concatenate((pareto_frontier, [row]))
        return pareto_frontier
    '''



    def pareto_frontier(self,Xs, Ys, maxX = True, maxY = False):
        # Sort the list in either ascending or descending order of X
        myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)

        # Start the Pareto frontier with the first value in the sorted list
        p_front = [myList[0]]
        # Loop through the sorted list
        for pair in myList[1:]:
            if maxY:
                if pair[1] >= p_front[-1][1]: # Look for higher values of Y
                    p_front.append(pair) #  and add them to the Pareto frontier
            else:
                if pair[1] <= p_front[-1][1]: # Look for lower values of Y
                    p_front.append(pair) #  and add them to the Pareto frontier
        # Turn resulting pairs back into a list of Xs and Ys
        p_frontX = [pair[0] for pair in p_front]
        p_frontY = [pair[1] for pair in p_front]
        return p_frontX, p_frontY


    def get_data(self,datafile,scale=1):
        data = {}
        dt = []
        with open(datafile,'r') as f:
            for line in f:
                tmp = line.split()

                if( len(tmp) == 2 ):
                    key = format( float(tmp[1]),'.3f' )
                    data[ key ] = tmp[0].split(':')[0]
                    dt.append( key )        

        return data, dt

    def get_pdbs(self):
        with open("get_pareto_pdbs.sh","w") as f:
            f.write("mkdir pareto_pdbs;\n")
            f.write("for i in ")
            for i in self.pdbs_on_front:
                f.write(i+" ")
            f.write(";\n")
            f.write("do\n cp $i pareto_pdbs/;\ndone;" )





    def main(self):

        parser = argparse.ArgumentParser(description="Compute the Pareto front between the list of data - default is maximum for both X and Y ( multiple dimensions are under construction")
        # get the initial rosetta design as input
        parser.add_argument("--data1", dest="data1", help="Data set 1" )

        parser.add_argument("--data2", dest="data2", help="Data set 2" )

        parser.add_argument("--xlabel", dest="xlabel", help="xlabel", default="SC" )

        parser.add_argument("--ylabel", dest="ylabel", help="ylabel", default="IFE" )

        # number of rounds to take pareto front
        parser.add_argument("--iterations", dest="iterations", help="Number of iterations that we will remove the pareto front and recompute a new one", type=int, default=3 )

        input_variables = parser.parse_args()

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        datafile_1 = input_variables.data1
        datafile_2 = input_variables.data2
        # two hash table
        x,x_array = self.get_data( datafile_1 )
        y,y_array = self.get_data( datafile_2 )

        for iteration in range(self.iterations):
            # tmp variables
            px,py = self.pareto_frontier(x_array, y_array, True, False )

            for i,j in zip(px, py):
                # print i,j
                self.pdbs_on_front.append( x[str(i)] )
                x.pop( str(i) )
                x_array.remove( i )
                y_array.remove( j )

            self.px = self.px + px
            self.py = self.py + py


        scatter(x_array,y_array)
        plot(self.px,self.py,'^',color='red')
        xlabel( input_variables.xlabel)
        ylabel( input_variables.ylabel)
        show()

        self.get_pdbs()



if __name__ == "__main__":
    run = get_union()
    run.main()

        

#keys_a = set(dict_a.keys())
#keys_b = set(dict_b.keys())
#intersection = keys_a & keys_b # '&' operator is used for set intersection
