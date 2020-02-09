import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from tqdm import tqdm
from datetime import datetime
from odeSolve import odeSolve
from scipy.optimize import fmin
from scipy.optimize import minimize
from scipy.optimize import least_squares

TOTAL_POPULATION = 7530000000.0
NUMRUNS = 100

def sir(t,x,cS,cR,cD,tS,tR,tD):
    s = x[0]
    i = x[1]
    if s < 0.0:
        s = 0.0
    if s > 1.0:
        s = 1.0
    if i < 0.0:
        i = 0.0
    if i > 1.0:
        i = 0.0

    # dsdt = -(cS - tS*t)*s*i
    # drdt =  (cR + tR*t)*i
    # dddt =  (cD - tD*t)*i

    dsdt = -cS*np.exp(-tS*t)*s*i
    drdt =  cR*np.exp( tR*t)*i
    dddt =  cD*np.exp(-tD*t)*i

    if dddt < 0.0:
        dddt = 0.0

    didt = -drdt - dsdt - dddt
    return np.array([dsdt,didt,drdt,dddt])

def getRes(t,X,iV,dV,rV):
    iVV = np.interp(t,X[:,0],X[:,2])
    rVV = np.interp(t,X[:,0],X[:,3])
    dVV = np.interp(t,X[:,0],X[:,4])
    res = (iV-iVV)/iV[-1]
    res = np.concatenate((res,(dV-dVV)/dV[-1]), axis = 0)
    res = np.concatenate((res,(rV-rVV)/rV[-1]), axis = 0)
    return res

def getStats(t,X,iV,dV,rV):
    iVV = np.interp(t,X[:,0],X[:,2])
    rVV = np.interp(t,X[:,0],X[:,3])
    dVV = np.interp(t,X[:,0],X[:,4])
    Yp = np.concatenate((iVV,rVV,dVV), axis = 0)
    Y = np.concatenate((iV,rV,dV), axis = 0)
    Ybar = np.mean(Y)
    res = Yp-Y
    SSE = np.sum(res*res)
    SST = np.sum((Y-Ybar)*(Y-Ybar))
    n = len(Y)
    p = 6.0
    R2 = 1.0 - SSE/SST
    R2adj = 1.0 - SSE*(n-1.0)/(SST*(n-p-1.0))
    return SSE, SST, R2, R2adj

def sirFit(b,t,iV,dV,rV):
    fun = lambda t,x: sir(t,x,b[0],b[1],b[2],b[3],b[4],b[5])
    ode = odeSolve(fun=fun, solver='scipy45')
    ode.rtol = 1.0e-12
    ode.atol = 1.0e-14
    ode.y0 = np.array([1.0-iV[0]-rV[0]-dV[0],iV[0],rV[0],dV[0]])
    ode.tspan = [t[0], t[-1]]
    ode.solve()
    return 100*getRes(t,ode.sol,iV,dV,rV)

def save(ode, resB):
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
    SSE, SST, R2, R2adj = getStats(tVals, ode.sol, iVals, dVals, rVals)
    print("     SSE =",SSE)
    print("     SST =",SST)
    print("      r2 =", R2)
    print("   r2adj =", R2adj)

    now = datetime.now()
    filename = "results_" + now.strftime("%Y_%m_%d_%H_%M_%S") + ".txt"
    with open(filename,'w') as f:
        printWrite(f, "===================================")
        printWrite(f, "======== Ending Statistics ========")
        printWrite(f, "             Dead: %10d" % int(TOTAL_POPULATION*ode.sol[-1,4]))
        printWrite(f, "        Recovered: %10d" % int(TOTAL_POPULATION*ode.sol[-1,3]))
        printWrite(f, "       Unaffected: %10d" % int(TOTAL_POPULATION*ode.sol[-1,1]))
        printWrite(f, "   Still Infected: %10d" % int(TOTAL_POPULATION*ode.sol[-1,2]))
        printWrite(f, "")
        printWrite(f, "   Coeficients")
        printWrite(f, "   cS: %e" % resB[0])
        printWrite(f, "   cR: %e" % resB[1])
        printWrite(f, "   cD: %e" % resB[2])
        printWrite(f, "   tS: %e" % resB[3])
        printWrite(f, "   tR: %e" % resB[4])
        printWrite(f, "   tD: %e" % resB[5])
        printWrite(f,"")
        printWrite(f, "   Fit Statistics")
        printWrite(f, "     SSE = %e" % SSE)
        printWrite(f, "     SST = %e" % SST)
        printWrite(f, "      r2 = %e" % R2)
        printWrite(f, "   r2adj = %e" % R2adj)
        printWrite(f, "===================================")

def printWrite(f, curStr):
    print(curStr)
    f.write(curStr + "\n")

if __name__ == "__main__":
	# #                Jan 23     24      25      26      27      28      29      30       31    Feb 1        2        3        4        5        6
    # iVals = np.array([653.0, 941.0, 2019.0, 2794.0, 4473.0, 6057.0, 7783.0, 9776.0, 11377.0, 14549.0, 17295.0, 20588.0, 24503.0, 24630.0, 30877.0])/TOTAL_POPULATION
    # rVals = np.array([ 30.0,  36.0,   49.0,   54.0,   63.0,  111.0,  133.0,  187.0,   252.0,   340.0,   487.0,   644.0,   899.0,  1029.0,  1499.0])/TOTAL_POPULATION
    # dVals = np.array([ 18.0,  26.0,   56.0,   80.0,  107.0,  133.0,  170.0,  213.0,   259.0,   305.0,   362.0,   426.0,   492.0,   494.0,   636.0])/TOTAL_POPULATION
    
    vals = np.genfromtxt("inputData.csv", delimiter=",")
    iVals = vals[:,1]/TOTAL_POPULATION
    rVals = vals[:,2]/TOTAL_POPULATION
    dVals = vals[:,3]/TOTAL_POPULATION

    tVals = [k for k in range(len(iVals))]

    fun = lambda b: sirFit(b,tVals,iVals,dVals,rVals)
    minCost = np.inf
    lb = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
    ub = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    for k in tqdm(range(NUMRUNS)):
        b0 = np.random.rand()*(ub-lb) + lb
        curRes = least_squares(
            fun=fun, 
            x0=b0,
            bounds=(lb,ub),
            verbose=0,
            loss='linear',
            xtol=5.0e-16,
            ftol=5.0e-16
        )
        b = curRes.x
        cost = curRes.cost

        if cost < minCost:
            print("   Min Cost =", cost)
            minCost = cost
            resB = b
            # break

    fun = lambda t,x: sir(t,x,resB[0],resB[1],resB[2],resB[3],resB[4],resB[5])
    ode = odeSolve(fun=fun, solver='scipy45')
    ode.rtol = 1.0e-12
    ode.atol = 1.0e-14
    ode.y0 = np.array([1.0-iVals[0]-rVals[0]-dVals[0],iVals[0],rVals[0],dVals[0]])
    ode.tspan = [tVals[0], 0.25*365]
    # ode.tspan = [tVals[0], tVals[-1]+1]
    ode.solve()
    save(ode,resB)
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
