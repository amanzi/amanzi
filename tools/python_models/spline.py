import numpy as np

class Spline(object):
    """Forms a cublic spline on an interval given values and derivatives at the endpoints of that interval."""
    def __init__(self, x1, y1, dy1, x2, y2, dy2):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.dy1 = dy1
        self.dy2 = dy2

    def T(self, x):
        return (x - self.x1) / (self.x2 - self.x1)

    def Value(self, x):
        t = self.T(x);
        return (1-t)**2 * ((1+2*t) * self.y1 + t * (self.x2 - self.x1) * self.dy1) \
          + t**2 * ((3-2*t) * self.y2 + (t-1) * (self.x2 - self.x1) * self.dy2)
        
    def Derivative(self, x):
        t = self.T(x)
        dtdx = 1./(self.x2 - self.x1)
        dydt = (6*t**2 - 6*t)* self.y1 \
          + (3*t**2 - 4*t + 1) * (self.x2 - self.x1) * self.dy1 \
          + (-6*t**2 + 6*t) * self.y2 \
          + (3*t**2 - 2*t) * (self.x2 - self.x1) * self.dy2
        return dydt * dtdx;

class SplineList(object):
    def __init__(self, x, y, dy, left_boundary=None, right_boundary=None, monotonicity_preserving=False):
        """Spline using arrays of x, y, and y'.

        Values outside of [x[0], x[-1]] depend upon left and right boundary options:
          "linear" (default)  | linear extrapolation using y' at the boundary
          "constant"          | gives the boundary value
          "extrapolate"       | simply continues the endmost spline
          "nan"               | returns NaN
        """

        assert len(x) > 1
        assert len(x) == len(y) == len(dy)
        assert (x[1:] - x[:-1]).min() > 0.

        if left_boundary is None:
            left_boundary = "linear"
        if right_boundary is None:
            right_boundary = "linear"
        self._left = left_boundary
        self._right = right_boundary
        
        if monotonicity_preserving:
            # ensure monotonicity of the spline, using the deBoor formula e.g. Hyman '83 monotonicity preserving spline
            S_half = (y[1:] - y[:-1])/(x[1:] - x[:-1])
            
            if S_half.mean() > 0:
                if (S_half.min() < 0.0):
                    raise RuntimeError("Spline: monotonicity_preserving requested but the data is not monotonic")
                dy_maxs = np.zeros(dy.shape,'d')
                dy_maxs[0] = np.inf
                dy_maxs[-1] = np.inf
                dy_maxs[1:-1] = 3*np.minimum(S_half[1:], S_half[:-1])
                dy = np.minimum(dy_maxs, dy)
                
            elif S_half.mean() < 0:
                if (S_half.max() > 0.0):
                    raise RuntimeError("Spline: monotonicity_preserving requested but the data is not monotonic")
                dy_mins = np.zeros(dy.shape,'d')
                dy_mins[0] = -np.inf
                dy_mins[-1] = -np.inf
                dy_mins[1:-1] = 3*np.maximum(S_half[1:], S_half[:-1])                   
                dy2 = np.maximum(dy_mins, dy)
                print "Spline slopes:", x[0:5], y[0:5], dy[0:5], dy2[0:5]
                dy = dy2
                    
            else:
                raise RuntimeError("Spline: monotonicity_preserving requested but the data is not monotonic")

        self._x = x
        self._splines = [Spline(x[i], y[i], dy[i], x[i+1], y[i+1], dy[i+1]) for i in range(len(x)-1)]

    def _index(self,x):
        i_list = np.where(x < self._x)[0]
        if len(i_list) == 0:
            # off the right edge
            if self._right == "extrapolate":
                return len(self._x)-2
            else:
                return len(self._x)-1
        elif i_list[0] == 0:
            # off the left edge
            if self._left == "extrapolate":
                return 0
            else:
                return -1
        else:
            return i_list[0]-1
        
    def Value(self, x):
        i = self._index(x)
        try:
            return self._splines[i].Value(x)
        except IndexError:
            if i == 0:
                # left boundary
                if self._left == "constant":
                    return self._splines[0].Value(self._x[0])
                elif self._left == "linear":
                    return self._splines[0].Value(self._x[0]) - (x - self._x[0])*self._splines[0].Derivative(self._x[0])
                else:
                    return np.nan
            elif i == len(self._splines):
                # right boundary
                if self._right == "constant":
                    return self._splines[-1].Value(self._x[-1])
                elif self._right == "linear":
                    return self._splines[-1].Value(self._x[-1]) + (x - self._x[-1])*self._splines[-1].Derivative(self._x[-1])
                else:
                    return np.nan
            else:
                raise IndexError("invalid index?  internal error! Asked for %d in array of length %d with value %g"%(i, len(self._splines), x))

    def Derivative(self, x):
        i = self._index(x)
        try:
            return self._splines[i].Derivative(x)
        except IndexError:
            if i == 0:
                # left boundary
                if self._left == "constant":
                    return 0.
                elif self._left == "linear":
                    return self._splines[0].Derivative(self._x[0])
                else:
                    return np.nan
            elif i == len(self._splines):
                # right boundary
                if self._right == "constant":
                    return 0.
                elif self._right == "linear":
                    return self._splines[-1].Derivative(self._x[-1])
                else:
                    return np.nan
            else:
                raise IndexError("invalid index?  internal error!")
            

    def __call__(self, x):
        return self.Value(x)        

  

if __name__ == "__main__":
    import numpy as np
    from matplotlib import pyplot as plt
    m = (-1.3 - 2.2) / (1.3 - 0.1)
    s = Spline(0.1,2.2,m, 1.3, -1.3, 0.5*m, True)
    x = np.linspace(0.1, 1.3, 10000)
    y = np.array([s.Value(lcv) for lcv in x])
    plt.figure()
    plt.subplot(211)
    plt.plot(x,y,'b')

    plt.subplot(212)
    dy = np.array([s.Derivative(lcv) for lcv in x])
    dy_num = (y[1:] - y[:-1]) / (x[1:] - x[:-1])
    xmean = (x[1:] + x[:-1])/2.
    plt.plot(x, dy,'b')
    plt.plot(xmean, dy_num,'r')    
    
    plt.show()
