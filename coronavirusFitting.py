import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from odeSolve import odeSolve
from scipy.optimize import least_squares


def sir(t,x,b,k):
    s = x[0]
    i = x[1]

    dsdt = -b*s*i
    drdt = k*i
    didt = -drdt - dsdt

    return np.array([dsdt,didt,drdt])

def sirFit(b,t,x):
    fun = lambda t,x: sir(t,x,b[0],b[1])
    ode = odeSolve(fun=fun, solver='rk45')
    ode.y0 = np.array([1.0,x[0],0.0])
    ode.tspan = [t[0], t[-1]]
    ode.solve()

    xx = np.interp(t,ode.sol[:,0],ode.sol[:,2])

    return x - xx

if __name__ == "__main__":
    # tVals = np.array([17,19, 20, 21, 22, 24,  25,  26])
    # xVals = np.array([41,62,201,291,440,830,1287,1975])/7713468000.0
    cVals = np.array([278.0, 326.0, 547.0, 639.0, 916.0, 2000.0, 2700.0, 4400.0, 6000.0, 7700.0, 8100.0])
    oVals = np.array([4.0, 6.0, 8.0, 14.0, 25.0, 40.0, 57.0, 64.0, 87.0, 105.0, 112.0,])
    xVals = (oVals + cVals)/7713468000.0
    tVals = [k for k in range(len(xVals))]

    fun = lambda b: sirFit(b,tVals,xVals)
    res = least_squares(fun, np.array([.5, 1.0/3.0]),verbose=2)
    print(res)

    fun = lambda t,x: sir(t,x,res.x[0],res.x[1])
    ode = odeSolve(fun=fun, solver='rk45')
    ode.y0 = np.array([1.0,xVals[0],0.0])
    ode.tspan = [tVals[0], 365.0]
    ode.solve()
    plt.plot(tVals,xVals,'k*')
    plt.plot(ode.sol[:,0], ode.sol[:,2])
    plt.show()
    ode.plot()