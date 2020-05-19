
import sys
import pyDOE
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import infectionFitting as infect
# from infectionFitting import *
    
global totalPopulation
totalPopulation = 327200000.0
numRuns = 25
infectionName = "coronavirus"



startDate = dt.date(2020,3,1)
endDate = dt.datetime.now() + dt.timedelta(days=-1)
endDate = dt.date(endDate.year,endDate.month,endDate.day)
baseURL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/'

def getRegions():
    date = endDate + dt.timedelta(days=-3)
    mdy = getDateStr(date)
    url = baseURL + mdy
    try:
        df = pd.read_csv(url,parse_dates=[2])
    except:
        print(url)
        df = pd.read_csv(url)
        print(df)
        exit()
    cr = 'Province_State'
    try:
        df[cr]
    except:
        cr = 'Province_State'
    regions = df[cr].unique()
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
    if region == 'US':
        return
    elif region == 'Iowa':
        totalPopulation = 3156000.0
    elif region == 'Wyoming':
        totalPopulation = 577737.0
    elif region == 'Idaha':
        totalPopulation = 1754000.0
    elif region == 'New York':
        totalPopulation = 8623000.0
    elif region == 'Florida':
        totalPopulation = 21300000.0
    elif region == 'South Dakota':
        totalPopulation = 882235.0
    elif region == 'California':
        totalPopulation = 39560000.0
    elif region == 'Louisiana':
        totalPopulation = 4660000.0
    elif region == 'Minnesota':
        totalPopulation = 5640000.0
    elif region == 'Montana':
        totalPopulation = 1069000.0
    elif region == 'Nebraska':
        totalPopulation = 1934000.0
    elif region == 'Idaho':
        totalPopulation = 1787000.0
    elif region == 'Texas':
        totalPopulation = 29000000.0
    else:
        print("WARNING: Using US population...")




def createTimeData(region):
    infected = []
    deaths = []
    days = []

    if not (region == 'US' or region == 'all'):
        posRegions = getRegions()
        if region not in posRegions:
            print("Region: %s, is not a posible region..." % region)
            print("Possible Regions:\n")
            print(posRegions)
            exit()
    
    getTotalPopulation(region)

    k = 0
    gotData = False
    while True:
        curDate = startDate + dt.timedelta(days=k)
        k+=1
        if curDate > endDate:
            break
        df = getData(curDate)
        cr = 'Province_State'
        try:
            df[cr]
        except:
            cr = 'Province_State'
        try:
            if df is None:
                continue
            if region == "all":
                break
            elif region == "US":
                infected.append(np.sum(df['Confirmed']))
                deaths.append(np.sum(df['Deaths']))
                days.append(k)
            else:
                gg = df.loc[df[cr]==region]
                numInfected = np.sum(gg['Confirmed'])
                numDeaths = np.sum(gg['Deaths'])
                if numInfected*numDeaths > 0:
                    infected.append(np.sum(gg['Confirmed']))
                    deaths.append(np.sum(gg['Deaths']))
                    days.append(k)
            gotData = True
        except:
            if gotData:
                print(df)
    
    infected = np.array(infected)
    deaths = np.array(deaths)
    days = np.array(days)

    ids = []
    temp = np.where(infected>0)
    ids.append(temp[0])
    temp = np.where(deaths>0)
    ids.append(temp[0])
    minIDlen = 1.0e10
    for idc in ids:
        if len(idc) < minIDlen:
            minIDlen = len(idc)
            minID = idc

    if minIDlen > 3:
        return infected[minID],deaths[minID],days[minID]
    else:
        print("Not enough data... Exiting...")
        print(days,infected,deaths)
        exit()

def printGrowthFactor(I):
    deltaI = I[1:]-I[0:-1]
    print(deltaI[1:]/deltaI[0:-1])


def printRawData(I,D,T):
    data = np.transpose(np.array([T,I,D]))
    print("Days, Infected, Deaths")
    print(data)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        region = sys.argv[1]
        if region == 'getRegions':
            regions = getRegions()
            print(regions)
            exit()
    else:
        region = 'Iowa'
    
    print("\nGetting Data")
    I,D,T = createTimeData(region)
    
    print("\nRaw Data")
    printRawData(I,D,T)
    print("\nGrowth Factor")
    printGrowthFactor(I)

    iVals = I/totalPopulation
    dVals = D/totalPopulation

    tVals = T

    fun = lambda b: infect.sirFit(b,tVals-tVals[0],iVals,dVals)
    lb = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
    ub = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    doeVals = pyDOE.lhs(len(lb),numRuns,'center')

    # Run on multiple cores...
    print("\nFitting Data")
    resB = infect.runOne(fun, lb, ub, doeVals, numRuns)

    fun = lambda t,x: infect.sir(t,x,resB[0],resB[1],resB[2],resB[3],resB[4],resB[5])
    ode = infect.odeSolve(fun = fun, solver = 'scipy45', events = infect.isZero)
    ode.rtol = 1.0e-12
    ode.atol = 1.0e-14
    ode.y0 = np.array([1.0-iVals[0]-dVals[0],iVals[0],dVals[0]])
    ode.tspan = [0, 100.0*365] # Simulate until no one is still infected...
    ode.solve()
    infect.save(ode, resB, totalPopulation, tVals-tVals[0], iVals, dVals, region, True)
    plt.plot(tVals-tVals[0],totalPopulation*iVals,'b*')
    plt.plot(ode.sol[:,0], (totalPopulation*ode.sol[:,2]),'b-',label="Infected")
    plt.plot(tVals-tVals[0],totalPopulation*dVals,'k*')
    plt.plot(ode.sol[:,0], (totalPopulation*ode.sol[:,3]),'k-',label="Deaths")
    plt.ylabel("Number of People")
    plt.xlabel("Time (Days)")
    plt.legend()
    plt.tight_layout()
    plt.show()