import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from tqdm import tqdm
from odeSolve import odeSolve
from scipy.optimize import fmin
from scipy.optimize import minimize
from scipy.optimize import least_squares

TOTAL_POPULATION = 7530000000.0


def sir(t,x,cS,cR,cD,tS,tR,tD):
    s = x[0]
    i = x[1]

    dsdt = -cS*s*i/(t/tS+1.0)
    drdt = cR*i*(t/tR+1.0)
    dddt = cD*i/(t/tD+1.0)
    didt = -drdt - dsdt - dddt

    return np.array([dsdt,didt,drdt,dddt])

def sirFit(b,t,iV,dV,rV):
    fun = lambda t,x: sir(t,x,b[0],b[1],b[2],b[3],b[4],b[5])
    ode = odeSolve(fun=fun, solver='rk45')
    ode.maxDt = 0.25
    ode.y0 = np.array([1.0-iV[0],iV[0],0.0,0.0])
    ode.tspan = [t[0], t[-1]]
    ode.solve()

    iVV = np.interp(t,ode.sol[:,0],ode.sol[:,2])
    rVV = np.interp(t,ode.sol[:,0],ode.sol[:,3])
    dVV = np.interp(t,ode.sol[:,0],ode.sol[:,4])

    res = (iV-iVV)/iV[-1]
    res = np.concatenate((res,(dV-dVV)/dV[-1]), axis = 0)
    res = np.concatenate((res,(rV-rVV)/rV[-1]), axis = 0)

    return 100*res

if __name__ == "__main__":
    #                    20     21     22     23     24      25      26      27      28      29      30
    iVals = np.array([282.0, 332.0, 555.0, 653.0, 941.0, 2019.0, 2794.0, 4473.0, 6057.0, 7783.0, 9776.0])/TOTAL_POPULATION
    rVals = np.array([  0.0,   0.0,   0.0,  30.0,  36.0,   49.0,   54.0,   63.0,  111.0,  133.0,  187.0])/TOTAL_POPULATION
    dVals = np.array([  0.0,   0.0,   0.0,  18.0,  26.0,   56.0,   80.0,  107.0,  133.0,  170.0,  213.0])/TOTAL_POPULATION
    
    tVals = [k for k in range(len(iVals))]

    fun = lambda b: sirFit(b,tVals,iVals,dVals,rVals)
    minCost = np.inf
    ub = np.array([1.0,1.0,1.0,365.0,365.0,365.0])
    for k in tqdm(range(50)):
        b0 = np.random.rand()*ub
        curRes = least_squares(
            fun=fun, 
            x0=b0,
            bounds=(0.0,ub),
            verbose=0,
            loss='linear',
            xtol=5.0e-16
        )
        b = curRes.x
        cost = curRes.cost

        if cost < minCost:
            print("   Min Cost =", cost)
            minCost = cost
            resB = b

    fun = lambda t,x: sir(t,x,resB[0],resB[1],resB[2],resB[3],resB[4],resB[5])
    ode = odeSolve(fun=fun, solver='rk45')
    ode.dt = 0.25

    ode.y0 = np.array([1.0-iVals[0],iVals[0],0.0,0.0])
    ode.tspan = [tVals[0], 2*365]
    # ode.tspan = [tVals[0], tVals[-1]+0]
    ode.solve()
    print("=======================================")
    print("========== Ending Statistics ==========")
    print("             Dead: %10d" % int(TOTAL_POPULATION*ode.sol[-1,4]))
    print("        Recovered: %10d" % int(TOTAL_POPULATION*ode.sol[-1,3]))
    print("       Unaffected: %10d" % int(TOTAL_POPULATION*ode.sol[-1,1]))
    print("   Still Infected: %10d" % int(TOTAL_POPULATION*ode.sol[-1,2]))
    print("   Coeficients")
    print("   cS:", resB[0])
    print("   cR:", resB[1])
    print("   cD:", resB[2])
    print("   tS:", resB[3])
    print("   tR:", resB[4])
    print("   tD:", resB[5])
    plt.plot(tVals,100*iVals,'b*')
    plt.plot(ode.sol[:,0], 100*ode.sol[:,2],'b-')
    plt.plot(tVals,100*rVals,'r*')
    plt.plot(ode.sol[:,0], 100*ode.sol[:,3],'r-')
    plt.plot(tVals,100*dVals,'k*')
    plt.plot(ode.sol[:,0], 100*ode.sol[:,4],'k-')
    plt.ylabel("Percent Population (%)")
    plt.xlabel("Time (Days)")
    plt.show()
    ode.plot()