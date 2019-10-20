import numpy as np
import pandas as pd


#class Pareto(object):
class Pareto:

    def __init__(self, df, xcol, ycol, ax=None):
        self.df = df
        self.xcol = xcol
        self.ycol = ycol
        self.ax = ax
        self.frontier = self.pareto_ix(df[xcol], df[ycol])
        assert np.isfinite(df[xcol]).all() and np.isfinite(df[ycol]).all(), \
            'Pareto: `DataFrame` contains non-finite values.'

    def scatter(self, **plot_kwargs):
        if self.ax is None:
            self.ax = pl.gca()
        kwargs = dict(lw=0)
        kwargs.update(plot_kwargs)
        self.ax.scatter(self.df[self.xcol], self.df[self.ycol], **kwargs)

    def show_frontier(self, **plot_kwargs):
        if self.ax is None:
            self.ax = pl.gca()
        show_frontier(self.df[self.xcol], self.df[self.ycol], ax=self.ax, **plot_kwargs)

    def lookup_x(self, x):
        """
        Find Pareto point constrained by x.
            argmax    p.y
         p: p.x <= x
        """
        s = self.df
        s = s[s[self.xcol] <= x]    # filter
        if s.empty:
            return np.nan
        yy = s[self.ycol].argmax()
        return s.ix[yy][self.ycol]

    def lookup_y(self, y):
        """
        Find Pareto point constrained by y.
            argmax    p.x
         p: p.x >= y
        """
        s = self.df
        s = s[s[self.ycol] >= y]    # filter
        if s.empty:
            return np.nan
        xx = s[self.xcol].argmin()
        return s.ix[xx][self.xcol]


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
        key = lambda xy: (xy[0] * xx, xy[1] * yy)   # need to break ties

        a = sorted(zip(X, Y, range(len(X))), key=key)
        frontier = []
        lastx = float('inf') * xx
        lasty = float('inf') * yy
        for xy in a:
            x,y,i = xy
            if yy*y <= yy*lasty:
                if lastx != x:
                    frontier.append(xy)
                lasty = y
                lastx = x

        x,y,i = zip(*frontier)
        return list(i) if indices else zip(x,y)

    def pareto_ix(self, X, Y, *a, **kw):
        "Determine Pareto frontier, returns list of indices"
        return self.pareto_frontier(X, Y, *a, **kw)
