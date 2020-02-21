from infectionFitting import *
    
totalPopulation = 7530000000.0
numRuns = 15
infectionName = "coronavirus"

vals = np.genfromtxt("inputDataCoronaVirus.csv", delimiter=",")
iVals = vals[:,1]/totalPopulation
rVals = vals[:,2]/totalPopulation
dVals = vals[:,3]/totalPopulation

tVals = [k for k in range(len(iVals))]

fun = lambda b: sirFit(b,tVals,iVals,dVals,rVals)
lb = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
ub = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
doeVals = pyDOE.lhs(len(lb),numRuns,'center')

# Run on multiple cores...
resB = runOne(fun, lb, ub, doeVals, numRuns)

fun = lambda t,x: sir(t,x,resB[0],resB[1],resB[2],resB[3],resB[4],resB[5])
ode = odeSolve(fun=fun, solver='scipy45')
ode.rtol = 1.0e-12
ode.atol = 1.0e-14
ode.y0 = np.array([1.0-iVals[0]-rVals[0]-dVals[0],iVals[0],rVals[0],dVals[0]])
ode.tspan = [tVals[0], 0.25*365]
ode.solve()
save(ode, resB, totalPopulation, tVals, iVals, rVals, dVals, infectionName)
plt.plot(tVals,100*iVals,'b*')
plt.plot(ode.sol[:,0], 100*ode.sol[:,2],'b-')
plt.plot(tVals,100*rVals,'r*')
plt.plot(ode.sol[:,0], 100*ode.sol[:,3],'r-')
plt.plot(tVals,100*dVals,'k*')
plt.plot(ode.sol[:,0], 100*ode.sol[:,4],'k-')
plt.ylabel("Percent Population (%)")
plt.xlabel("Time (Days)")
plt.show()
# ode.plot()
