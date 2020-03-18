
import sys
import numpy as np
import pandas as pd
import datetime as dt
from infectionFitting import *
    
global totalPopulation
totalPopulation = 7530000000.0
numRuns = 15
infectionName = "coronavirus"



startDate = dt.date(2020,1,22)
endDate = dt.datetime.now() + dt.timedelta(days=-1)
endDate = dt.date(endDate.year,endDate.month,endDate.day)
baseURL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/'

def getRegions():
    date = endDate + dt.timedelta(days=-3)
    mdy = getDateStr(date)
    url = baseURL + mdy
    df = pd.read_csv(url,parse_dates=[2])
    regions = df['Country/Region'].unique()
    return regions

def getData(date):
    mdy = getDateStr(date)
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

def getTotalPopulation(region):
    global totalPopulation
    if region == 'world':
        return
    elif region == 'US':
        totalPopulation = 327200000.0
    elif region == 'Italy':
        totalPopulation = 60480000.0
    elif region == 'China':
        totalPopulation = 1386000000.0
    elif region == 'Mainland China':
        totalPopulation = 1386000000.0
    elif region == 'Panama':
        totalPopulation = 4099000.0
    elif region == 'Australia':
        totalPopulation = 24600000.0
    elif region == 'Croatia':
        totalPopulation = 4076000.0
    elif region == 'Burkina Faso':
        totalPopulation = 19190000.0
    elif region == 'Albania':
        totalPopulation = 2877000.0
    elif region == 'Canada':
        totalPopulation = 37590000.0
    elif region == 'United Kingdom':
        totalPopulation = 66440000.0
    elif region == 'Philippines':
        totalPopulation = 104900000.0
    elif region == 'Nepal':
        totalPopulation = 29300000.0
    elif region == 'Sri Lanka':
        totalPopulation = 21440000.0
    elif region == 'United Arab Emirates':
        totalPopulation = 9400000.0
    elif region == 'Thailand':
        totalPopulation = 69040000.0
    elif region == 'Iran':
        totalPopulation = 81160000.0
    elif region == 'Iraq':
        totalPopulation = 38270000.0
    elif region == 'Mexico':
        totalPopulation = 129200000.0
    elif region == 'Germany':
        totalPopulation = 82790000.0
    elif region == 'Turkey':
        totalPopulation = 80810000.0
    elif region == 'France':
        totalPopulation = 66890000.0
    elif region == 'Korea, South':
        totalPopulation = 51470000.0
    elif region == 'India':
        totalPopulation = 1339000000.0
    else:
        print("WARNING: Using world population...")




def createTimeData(region):
    infected = []
    recoverd = []
    deaths = []
    days = []

    if not (region == 'world' or region == 'all'):
        posRegions = getRegions()
        if region not in posRegions:
            print("Region: %s, is not a posible region..." % region)
            print("Possible Regions:\n")
            print(posRegions)
            exit()
    
    getTotalPopulation(region)

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
        elif region == "world":
            infected.append(np.sum(df['Confirmed']))
            recoverd.append(np.sum(df['Recovered']))
            deaths.append(np.sum(df['Deaths']))
            days.append(k)
        else:
            if region == 'China':
                gg = df.loc[df['Country/Region']=='China']
                if len(gg) == 0:
                    gg = gg = df.loc[df['Country/Region']=='Mainland China']
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

    if minIDlen > 3:
        return infected[minID],recoverd[minID],deaths[minID],days[minID]
    else:
        print("Not enough data... Exiting...")
        exit()

def printGrowthFactor(I):
    deltaI = I[1:]-I[0:-1]
    print(deltaI[1:]/deltaI[0:-1])


def printRawData(I,R,D,T):
    data = np.transpose(np.array([T,I,R,D]))
    print("Days, Infected, Recovered, Deaths")
    print(data)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        region = sys.argv[1]
        if region == 'getRegions':
            regions = getRegions()
            print(regions)
            exit()
    else:
        region = 'China'
    
    print("\nGetting Data")
    I,R,D,T = createTimeData(region)
    
    print("\nRaw Data")
    printRawData(I,R,D,T)
    print("\nGrowth Factor")
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
    print("\nFitting Data")
    resB = runOne(fun, lb, ub, doeVals, numRuns)

    fun = lambda t,x: sir(t,x,resB[0],resB[1],resB[2],resB[3],resB[4],resB[5])
    ode = odeSolve(fun = fun, solver = 'scipy45', events = isZero)
    ode.rtol = 1.0e-12
    ode.atol = 1.0e-14
    ode.y0 = np.array([1.0-iVals[0]-rVals[0]-dVals[0],iVals[0],rVals[0],dVals[0]])
    ode.tspan = [0, 100.0*365] # Simulate until no one is still infected...
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