import time
import numpy as np
import scipy as sp

from tqdm import tqdm
from datetime import datetime
from odeSolve import odeSolve
from scipy.optimize import least_squares

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

def isZero(t,x):
    if x[1] < 0.0:
        return -1.0
    else:
        return  1.0

def save(ode, resB, totalPopulation, tVals, iVals, rVals, dVals, region, write2file=True):
    SSE, SST, R2, R2adj = getStats(tVals, ode.sol, iVals, dVals, rVals)
    
    if write2file:
        reg = region.replace(' ','_')
        now = datetime.now()
        filename = "results/" + reg + now.strftime("_%Y_%m_%d_%H_%M_%S") + ".txt"
        with open(filename,'w') as f:
            printWrite(f, "===================================", write2file)
            printWrite(f, "======== Ending Statistics ========", write2file)
            printWrite(f, "             Dead: %10d" % int(totalPopulation*ode.sol[-1,4]), write2file)
            printWrite(f, "        Recovered: %10d" % int(totalPopulation*ode.sol[-1,3]), write2file)
            printWrite(f, "       Unaffected: %10d" % int(totalPopulation*ode.sol[-1,1]), write2file)
            printWrite(f, "   Still Infected: %10d" % int(totalPopulation*ode.sol[-1,2]), write2file)
            printWrite(f, "", write2file)
            printWrite(f, "   Coeficients", write2file)
            printWrite(f, "   cS: %e" % resB[0], write2file)
            printWrite(f, "   cR: %e" % resB[1], write2file)
            printWrite(f, "   cD: %e" % resB[2], write2file)
            printWrite(f, "   tS: %e" % resB[3], write2file)
            printWrite(f, "   tR: %e" % resB[4], write2file)
            printWrite(f, "   tD: %e" % resB[5], write2file)
            printWrite(f,"", write2file)
            printWrite(f, "   Fit Statistics", write2file)
            printWrite(f, "     SSE = %e" % SSE, write2file)
            printWrite(f, "     SST = %e" % SST, write2file)
            printWrite(f, "      r2 = %e" % R2, write2file)
            printWrite(f, "   r2adj = %e" % R2adj, write2file)
            printWrite(f, "===================================", write2file)
    else:
        f = None
        printWrite(f, "===================================", write2file)
        printWrite(f, "======== Ending Statistics ========", write2file)
        printWrite(f, "             Dead: %10d" % int(totalPopulation*ode.sol[-1,4]), write2file)
        printWrite(f, "        Recovered: %10d" % int(totalPopulation*ode.sol[-1,3]), write2file)
        printWrite(f, "       Unaffected: %10d" % int(totalPopulation*ode.sol[-1,1]), write2file)
        printWrite(f, "   Still Infected: %10d" % int(totalPopulation*ode.sol[-1,2]), write2file)
        printWrite(f, "", write2file)
        printWrite(f, "   Coeficients", write2file)
        printWrite(f, "   cS: %e" % resB[0], write2file)
        printWrite(f, "   cR: %e" % resB[1], write2file)
        printWrite(f, "   cD: %e" % resB[2], write2file)
        printWrite(f, "   tS: %e" % resB[3], write2file)
        printWrite(f, "   tR: %e" % resB[4], write2file)
        printWrite(f, "   tD: %e" % resB[5], write2file)
        printWrite(f,"", write2file)
        printWrite(f, "   Fit Statistics", write2file)
        printWrite(f, "     SSE = %e" % SSE, write2file)
        printWrite(f, "     SST = %e" % SST, write2file)
        printWrite(f, "      r2 = %e" % R2, write2file)
        printWrite(f, "   r2adj = %e" % R2adj, write2file)
        printWrite(f, "===================================", write2file)

def printWrite(f, curStr, write2file):
    print(curStr)
    if write2file:
        f.write(curStr + "\n")

def runOne(fun, lb, ub, doeVals, numRuns):
    minCost = np.inf
    for k in tqdm(range(numRuns)):
        b0 = doeVals[k,:]*(ub-lb) + lb
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
    return resB
