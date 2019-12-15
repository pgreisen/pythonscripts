class ParetoAnalysis:

    def __init__(self):
        print("initialized again")

    def pareto_frontier(self, X, Y, maxX=True, maxY=True, indices=False):
        """Determine Pareto frontier, returns list of sorted points.
        Args:
          X, Y: data.
          maxX, maxY: (bool) whether to maximize or minimize along respective
            coordinate.
        """
        assert len(X) == len(Y)
        if len(X) == 0:
            return []

        xx = -1 if maxX else +1
        yy = -1 if maxY else +1
        key = lambda xy: (xy[0] * xx, xy[1] * yy)  # need to break ties
        a = sorted(zip(X, Y, range(len(X))), key=key)
        frontier = []
        lastx = float('inf') * xx
        lasty = float('inf') * yy
        for xy in a:
            x, y, i = xy
            if yy * y <= yy * lasty:
                if lastx != x:
                    frontier.append(xy)
                lasty = y
                lastx = x

        x, y, i = zip(*frontier)
        return x, y, i

    def getIndexes(self, dfObj, value):
        ''' Get index positions of value in dataframe i.e. dfObj.'''

        listOfPos = list()
        # Get bool dataframe with True at positions where the given value exists
        result = dfObj.isin([value])
        # Get list of columns that contains the value
        seriesObj = result.any()
        columnNames = list(seriesObj[seriesObj == True].index)
        # Iterate over list of columns and fetch the rows indexes where value exists
        for col in columnNames:
            rows = list(result[col][result[col] == True].index)
            for row in rows:
                listOfPos.append(row)
        # Return a list of tuples indicating the positions of value in the dataframe
        return listOfPos


    def get_multiple_iterations_on_pareto_front(self, dfinit, xlabel, ylabel, maxX=True,maxY=True,  number_of_rounds=2):
        df = dfinit.copy()
        index_of_hits = []
        for iterations in range(0,number_of_rounds):
            x_,y_,i_ = self.pareto_frontier(df[xlabel],df[ylabel],maxX,maxY)
            # assumes that one value is enough
            for xtmp in x_:
                lst = self.getIndexes(dfinit, xtmp)
                for i in lst:
                    index_of_hits.append(i)
                    df.drop(index=i,inplace=True)
        return index_of_hits

