
import sys
import numpy as np
import pandas as pd
import datetime as dt
from infectionFitting import *
    
totalPopulation = 7530000000.0
numRuns = 15
infectionName = "coronavirus"



startDate = dt.date(2020,1,22)
endDate = dt.datetime.now() + dt.timedelta(days=-1)
endDate = dt.date(endDate.year,endDate.month,endDate.day)
baseURL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/'
# gg = df.loc[df['Country/Region']=='Mainland China']
# region = df['Country/Region'].unique()

def getRegions():
    date = endDate + dt.timedelta(days=-3)
    mdy = getDateStr(date)
    url = baseURL + mdy
    print(url)
    df = pd.read_csv(url,parse_dates=[2])
    regions = df['Country/Region'].unique()
    return regions

def getData(date):
    mdy = getDateStr(date)
    print(mdy)
    url = baseURL + mdy
    try:
        df = pd.read_csv(url,parse_dates=[2])
    except:
        df = None
    return df

def getDateStr(date):
    m = date.month
    d = date.day
    if m < 10:
        mdy = "0%s-"%m
    else: 
        mdy = "%s-"%m
    if d < 10:
        mdy += "0%s-2020.csv"%d
    else:
        mdy += "%s-2020.csv"%d
    return mdy

def createTimeData(region):
    infected = []
    recoverd = []
    deaths = []
    days = []

    posRegions = getRegions()
    if region not in posRegions:
        print("Region: %s, is not a posible region..." % region)
        print("Possible Regions:\n")
        print(posRegions)
        exit()

    k = 0
    while True:
        curDate = startDate + dt.timedelta(days=k)
        k+=1
        if curDate > endDate:
            break
        df = getData(curDate)
        if df is None:
            continue
        if region == "all":
            break
        else:
            gg = df.loc[df['Country/Region']==region]
            infected.append(np.sum(gg['Confirmed']))
            recoverd.append(np.sum(gg['Recovered']))
            deaths.append(np.sum(gg['Deaths']))
            days.append(k)
    
    infected = np.array(infected)
    recoverd = np.array(recoverd)
    deaths = np.array(deaths)
    days = np.array(days)

    ids = []
    temp = np.where(infected>0)
    ids.append(temp[0])
    temp = np.where(recoverd>0)
    ids.append(temp[0])
    temp = np.where(deaths>0)
    ids.append(temp[0])
    minIDlen = 1.0e10
    for idc in ids:
        if len(idc) < minIDlen:
            minIDlen = len(idc)
            minID = idc

    return infected[minID],recoverd[minID],deaths[minID],days[minID]

def printGrowthFactor(I):
    deltaI = I[1:]-I[0:-1]
    print(deltaI[1:]/deltaI[0:-1])


if __name__ == "__main__":
    if len(sys.argv) > 1:
        region = sys.argv[1]
        if region == 'getRegions':
            regions = getRegions()
            print(regions)
            exit()
    else:
        region = 'US'
    
    I,R,D,T = createTimeData(region)
    printGrowthFactor(I)

    iVals = I/totalPopulation
    rVals = R/totalPopulation
    dVals = D/totalPopulation

    tVals = T

    fun = lambda b: sirFit(b,tVals-tVals[0],iVals,dVals,rVals)
    lb = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
    ub = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    doeVals = pyDOE.lhs(len(lb),numRuns,'center')

    # Run on multiple cores...
    resB = runOne(fun, lb, ub, doeVals, numRuns)

    fun = lambda t,x: sir(t,x,resB[0],resB[1],resB[2],resB[3],resB[4],resB[5])
    ode = odeSolve(fun = fun, solver = 'scipy45', events = isZero)
    ode.rtol = 1.0e-12
    ode.atol = 1.0e-14
    ode.y0 = np.array([1.0-iVals[0]-rVals[0]-dVals[0],iVals[0],rVals[0],dVals[0]])
    ode.tspan = [0, 10.0*365] # Simulate until no one is still infected...
    ode.solve()
    save(ode, resB, totalPopulation, tVals-tVals[0], iVals, rVals, dVals, infectionName, False)
    plt.plot(tVals-tVals[0],totalPopulation*iVals,'b*')
    plt.plot(ode.sol[:,0], (totalPopulation*ode.sol[:,2]),'b-',label="Infected")
    plt.plot(tVals-tVals[0],totalPopulation*rVals,'r*')
    plt.plot(ode.sol[:,0], (totalPopulation*ode.sol[:,3]),'r-',label="Recovered")
    plt.plot(tVals-tVals[0],totalPopulation*dVals,'k*')
    plt.plot(ode.sol[:,0], (totalPopulation*ode.sol[:,4]),'k-',label="Deaths")
    plt.ylabel("Number of People")
    plt.xlabel("Time (Days)")
    plt.legend()
    plt.tight_layout()
    plt.show()