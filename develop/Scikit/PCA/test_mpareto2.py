def simple_cull(inputPoints, dominates):
    paretoPoints = set()
    candidateRowNr = 0
    dominatedPoints = set()
    while len(inputPoints):
        candidateRow = inputPoints[candidateRowNr]
        inputPoints.remove(candidateRow)
        rowNr = 0
        nonDominated = True
        while len(inputPoints) != 0 and rowNr < len(inputPoints):
            row = inputPoints[rowNr]
            if dominates(candidateRow, row):
                inputPoints.remove(row)
                dominatedPoints.add(tuple(row))
            elif dominates(row, candidateRow):
                nonDominated = False
                dominatedPoints.add(tuple(candidateRow))
                rowNr += 1
            else:
                rowNr += 1
                
        if nonDominated:
            # add the non-dominated point to the Pareto frontier
            paretoPoints.add(tuple(candidateRow))

    return paretoPoints, dominatedPoints

def dominates(row, anotherRow):
    return sum([row[x] >= anotherRow[x] for x in range(len(row))]) == len(row) # maximization domination

def testPoints(points):
    paretoPoints, dominatedPoints = simple_cull(list(points), dominates)
    print("\nInput points")
    print(points)
    print("*"*8 + " non-dominated answers " + ("*"*8))
    print(paretoPoints)
    print("*"*8 + " dominated answers " + ("*"*8))
    print (dominatedPoints)
    
testPoints([[1,1,1],[2,2,2],[3,3,3],[4,4,4],[3,5,6],[1,3,7],[4,6,2]])

testPoints([[1,2],[2,1]])

testPoints([[1],[2],[2],[1],[2],[3],[7]])
