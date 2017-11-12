import numpy as np


def pareto_frontier(Xs, Ys, maxX = False, maxY = False):
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    p_front = [myList[0]]    
    for pair in myList[1:]:
        if maxY: 
            if pair[0] == p_front[-1][0] and pair[1] >= p_front[-1][1]:
                p_front.pop(-1)
                p_front.append(pair)
            elif pair[1] >= p_front[-1][1]:
                p_front.append(pair)
        else:
            if pair[0] == p_front[-1][0] and pair[1] <= p_front[-1][1]:
                p_front.pop(-1)
                p_front.append(pair)
            elif pair[1] <= p_front[-1][1]:
                p_front.append(pair)
            
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY

def pareto_frontier(Xs, Ys, maxX = False, maxY = False):
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    p_front = [myList[0]]    
    for pair in myList[1:]:
        if maxY: 
            if pair[0] == p_front[-1][0] and pair[1] >= p_front[-1][1]:
                p_front.pop(-1)
                p_front.append(pair)
            elif pair[1] >= p_front[-1][1]:
                p_front.append(pair)
        else:
            if pair[0] == p_front[-1][0] and pair[1] <= p_front[-1][1]:
                p_front.pop(-1)
                p_front.append(pair)
            elif pair[1] <= p_front[-1][1]:
                p_front.append(pair)
            
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY

 
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


def pareto_frontier_n(Xs, Ys, maxX = False, maxY = False):
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    p_front = [myList[0]]    
    for pair in myList[1:]:
        if maxY: 
            if pair[0] == p_front[-1][0] and pair[1] >= p_front[-1][1]:
                p_front.pop(-1)
                p_front.append(pair)
            elif pair[1] >= p_front[-1][1]:
                p_front.append(pair)
        else:
            if pair[0] == p_front[-1][0] and pair[1] <= p_front[-1][1]:
                p_front.pop(-1)
                p_front.append(pair)
            elif pair[1] <= p_front[-1][1]:
                p_front.append(pair)
            
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY
 
myArray = np.array([[1,1,1],[2,2,2],[4,4,4],[6,5,3],[1,1,3],[9,3,2],[5,3,5]])
print pareto_frontier_multi(myArray)
