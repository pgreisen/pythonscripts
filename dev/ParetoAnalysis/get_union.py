from pylab import *
import sys,os
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
        self.residue_numbers_close_to_metal_site = [199, 254, 101,60,253]
        self.mutations_positions = defaultdict(list)
        self.top_mutations_to_remove = ["T_54_L", "T_54_M", "T_350_L","T_350_M", "G_254_R", "A_49_M", "T_350_F", "S_231_L", "N_321_R", "T_45_P","T_97_M", "D_253_W", "T_199_W", "T_199_F","L_249_W", "G_60_W", "G_251_H", "S_258_W", "G_254_R", "D_253_F", "D_253_T", "D_253_Y", "L_249_F", "T_350_Y" ]

    def pareto_frontier(self,Xs, Ys, maxX = True, maxY = True):
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


    def get_data(self,datafile):
        data = {}
        with open(datafile,'r') as f:
            for line in f:
                tmp = line.split(',')

                if( tmp[0] in self.top_mutations_to_remove):
                    continue

                if( len(tmp) == 2 ):
                    native, residue_number, mutation  = tmp[0].split('_')
                    # for some reason FoldX adds an o for histidine 
                    # check this !
                    if( mutation == 'o' ):
                        mutation = 'H'
                    if( mutation == 'e' ):
                        mutation = 'H'
                    # for self mutations it is not a capital letter
                    if ( mutation.isupper() == False ):
                        print "Self mutation", native, residue_number, mutation
                        mutation = mutation.upper()
                    ##print native, residue_number, mutation
                    if( int( residue_number) not in self.residue_numbers_close_to_metal_site):

                        self.mutations_positions[residue_number].append(mutation)

                        data[tmp[0] ] = tmp[1]
        return data

    def get_intersection_values(self,intersection,dt1,dt2):
        values = {}
        for i in intersection:
            values[ dt1[i] ] = i+"Rosetta"
            values[ dt2[i] ] = i+"FoldX"
        return values

    def main(self):
        datafile_1 = sys.argv[1]
        datafile_2 = sys.argv[2]

        dt1 = self.get_data( datafile_1 )
        dt2 = self.get_data( datafile_2 )

        keys_a = set(dt1.keys())
        keys_b = set(dt2.keys())
        intersection = keys_a & keys_b # '&' operator is used for set intersection
        ##print intersection

        x = []
        y = []
        labels = []
        for i in intersection:
            x.append( dt1[i]  )
            y.append( dt2[i]  )
            labels.append(i)

        scatter(x,y,)
        xlabel("ddG_FoldX")
        ylabel("ddG_Rosetta")

        for label, i, j in zip(labels, x, y):     
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        show()

        px,py = self.pareto_frontier(x, y )
        #print px
        #print py

        values = self.get_intersection_values(intersection,dt1,dt2)
        for i in px:
            print values[i]
        for i in py:
            print values[i]

        protein_position = []
        number_of_mutations = []

        for pos in self.mutations_positions:

            protein_position.append( pos )
            tmp_set = len( set ( self.mutations_positions[pos] ) )
            number_of_mutations.append( tmp_set )

            print tmp_set, pos, self.mutations_positions[pos]

        plot(protein_position, number_of_mutations,'o')
        show()


if __name__ == "__main__":
    run = get_union()
    run.main()

        

#keys_a = set(dict_a.keys())
#keys_b = set(dict_b.keys())
#intersection = keys_a & keys_b # '&' operator is used for set intersection
