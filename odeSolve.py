import numpy as np
import matplotlib.pyplot as plt

class odeSolve():
    def __init__(self, fun=None, solver='rk45'):
        self.solver = solver
        self.tspan = []
        self.sol = None
        self.fun = fun
        self.y0 = None
        self.numSteps = None
        self.rtol = 1.0e-3
        self.atol = 1.0e-6
        self.dt = 1.0e-6
        self.maxDt = np.inf
        self.failed = 0
        self.timeStep = []

    def solve(self):
        if self.solver == 'euler':
            self.sol = self.euler()
        elif self.solver == 'mod-euler':
            self.sol = self.mod_euler()
        elif self.solver == 'rk4':
            self.sol = self.rk4()
        elif self.solver == 'rk12':
            self.sol = self.rk12()
        elif self.solver == 'rk45':
            self.sol = self.rk45()
        elif self.solver == "scipy45":
            self.sol = self.sp45()
        else:
            print("ERROR: Solver not found. Current Solvers Include")
            print(
                [
                'euler',
                'mod-euler',
                'rk4',
                'rk12',
                'rk45',
                'scipy45'
                ]
            )
    
    def sp45(self):
        from scipy.integrate import solve_ivp as ode
        
        sol = ode(
            fun=self.fun,
            t_span=self.tspan,
            y0=self.y0,
            method='RK45',
            max_step=self.maxDt,
            rtol=self.rtol,
            atol=self.atol
        )

        count = len(sol.t)
        numVars = len(self.y0)
        out = np.zeros([count,numVars+1])
        out[:,0] = sol.t
        out[:,1:] = np.transpose(sol.y)
        dt = np.gradient(sol.t)
        self.timeStep = dt
        return out

    def rk45(self):
        t0 = self.tspan[0]
        tend = self.tspan[1]

        yn = np.array([self.y0])
        tn = np.array([t0])
        dt = self.dt
        self.timeStep.append(dt)

        oneOnFive = 1.0/5.0

        numVars = len(self.y0)

        done = False
        count = 0
        while not done:
            if tn[-1] + dt > tend:
                dt = tend - tn[-1]
                done = True
            
            y = yn[-1,:]
            t = tn[-1]

            k1 = dt*self.fun(t, y)
            k2 = dt*self.fun(t+0.25*dt, y + 0.25*k1)
            k3 = dt*self.fun(t+0.375*dt, y + (3.0*k1+9.0*k2)/32.0)
            k4 = dt*self.fun(t+12.0/13.0*dt, y + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
            k5 = dt*self.fun(t+dt, y + (439.0*k1 - 8.0*k2 + 3680.0*k3 - 845.0*k4)/4104.0)
            k6 = dt*self.fun(t+0.5*dt, y - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)

            y1 = y + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 - 9.0*k5/50.0 + 2.0*k6/55.0
            y2 = y + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0
            fE = y2-y1
            err = dt * np.linalg.norm(fE/(np.max([np.max([np.abs(y),np.abs(y2)]),5.0e-16])),ord=np.inf)

            if err < self.rtol: 
                count += 1
                z = np.zeros([count+1,numVars])
                z[:-1,:] = yn
                z[-1,:] = y2
                yn = z
                tn = np.append(tn,tn[-1]+dt)
                self.timeStep.append(dt)
            else:
                done = False
                self.failed += 1
            dt = 0.8 * dt * np.power((self.rtol/err), oneOnFive)
            if dt > self.maxDt:
                dt = self.maxDt
        
        out = np.zeros([count+1,numVars+1])
        out[:,0] = tn
        out[:,1:] = yn
        return out

    def rk12(self):
        t0 = self.tspan[0]
        tend = self.tspan[1]

        yn = np.array([self.y0])
        tn = np.array([t0])
        dt = self.dt
        self.timeStep.append(dt)

        numVars = len(self.y0)

        done = False
        count = 0
        while not done:
            if tn[-1] + dt > tend:
                dt = tend - tn[-1]
                done = True
            
            k1 = dt*self.fun(tn[-1], yn[-1,:])
            k2 = dt*self.fun(tn[-1]+dt, yn[-1,:] + k1)

            y1 = yn[-1,:] + k1
            y2 = yn[-1,:] + 0.5*(k1 + k2)

            err = np.sqrt(np.sum(np.power(y2-y1,2.0)))

            if err < self.atol:
                count += 1
                z = np.zeros([count+1,numVars])
                z[:-1,:] = yn
                z[-1,:] = y2
                yn = z
                tn = np.append(tn,tn[-1]+dt)
                self.timeStep.append(dt)
            else:
                done = False
            dt = 0.9*dt*np.power(self.atol/(1.0*err+1.0e-15),1/2)
            if dt > self.maxDt:
                dt = self.maxDt
        
        out = np.zeros([count+1,numVars+1])
        out[:,0] = tn
        out[:,1:] = yn
        return out


    
    def rk4(self):
        t0 = self.tspan[0]
        tend = self.tspan[1]
        if self.numSteps is None:
            print("Need number of steps for euler solver")
        out = np.zeros([self.numSteps+1,1+np.size(self.y0)])
        out[0,1:] = self.y0
        out[0,0] = t0
        dt = (tend-t0)/self.numSteps
        self.timeStep.append(dt)

        for k in range(self.numSteps):
            tn = out[k,0]
            yn = out[k,1:]
            k1 = dt * self.fun(tn,yn)
            k2 = dt * self.fun(tn+0.5*dt, yn + 0.5*k1)
            k3 = dt * self.fun(tn+0.5*dt, yn + 0.5*k2)
            k4 = dt * self.fun(tn+dt, yn + k3)
            out[k+1,1:] = out[k,1:] + (k1+2.0*k2+2.0*k3+k4)/6.0
            out[k+1,0] = out[k,0] + dt
            self.timeStep.append(dt)

        return out


    def mod_euler(self):
        t0 = self.tspan[0]
        tend = self.tspan[1]
        if self.numSteps is None:
            print("Need number of steps for euler solver")
        
        out = np.zeros([self.numSteps+1,1+np.size(self.y0)])
        out[0,1:] = self.y0
        out[0,0] = t0
        dt = (tend-t0)/self.numSteps
        self.timeStep.append(dt)

        for k in range(self.numSteps):
            k1 = out[k,1:] + dt*self.fun(out[k,0],out[k,1:])
            out[k+1,1:] = out[k,1:] + 0.5*dt*(self.fun(out[k,0],out[k,1:]) + self.fun(out[k,0]+dt,k1))
            out[k+1,0] = out[k,0] + dt
            self.timeStep.append(dt)

        return out

    def euler(self):
        t0 = self.tspan[0]
        tend = self.tspan[1]
        if self.numSteps is None:
            print("Need number of steps for euler solver")
        
        out = np.zeros([self.numSteps+1,1+np.size(self.y0)])
        out[0,1:] = self.y0
        out[0,0] = t0
        dt = (tend-t0)/self.numSteps
        self.timeStep.append(dt)
        for k in range(self.numSteps):
            out[k+1,1:] = out[k,1:] + dt*self.fun(out[k,0],out[k,1:])
            out[k+1,0] = out[k,0] + dt
            self.timeStep.append(dt)
        
        return out
    
    def plot(self):
        plt.subplot(211)
        sol = self.sol
        if sol is None:
            print("No Solution available...")
            return
        t = sol[:,0]
        plt.plot(t,self.timeStep)
        plt.ylabel("Time Step Size")

        plt.subplot(212)
        y = sol[:,1:]
        plt.plot(t,y,'-o')
        ylabelString = []
        for k in range(len(self.y0)):
            ylabelString.append("x[%i]" % k)
        plt.xlabel("Time")
        plt.legend(ylabelString)
        plt.show()

