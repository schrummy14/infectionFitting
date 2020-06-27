
import sys
import pyDOE
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import infectionFitting as infect
# from infectionFitting import *
    
global totalPopulation
totalPopulation = 3156000.0
numRuns = 5
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
    df = df[df[cr]=='Iowa']
    cr = 'Admin2'
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
    y = date.year
    # print(date.day,date.month,date.year)
    if m < 10:
        mdy = "0%s-"%m
    else: 
        mdy = "%s-"%m
    # if d < 10:
    #     mdy += "0%s-2020.csv"%d
    # else:
    #     mdy += "%s-2020.csv"%d

    if d < 10:
        mdy += "0%s-"%d
    else:
        mdy += "%s-"%d

    mdy += "%s.csv"%y

    # print(mdy)

    return mdy

def getTotalPopulation(region):
    global totalPopulation
    if region == 'all':
        return
    elif region == 'Story':
        totalPopulation = 97117.0
    elif region == 'Sac':
        totalPopulation = 9719.0
    elif region == 'Polk':
        totalPopulation = 490161.0
    elif region == 'Linn':
        totalPopulation = 226706.0
    elif region == 'Wapello':
        totalPopulation = 34969.0
    else:
        print("WARNING: Using State of Iowa Population...")




def createTimeData(region):
    infected = []
    deaths = []
    days = []

    infectedAll = []
    deathsAll = []
    daysAll = []

    if not (region == 'all'):
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
            df = df[df[cr]=='Iowa']
            cr = 'Admin2'
        except:
            cr = 'Province_State'
        try:
            if df is None:
                continue
            if region == "all":
                infected.append(np.sum(df['Confirmed']))
                deaths.append(np.sum(df['Deaths']))
                days.append(k)
            else:
                gg = df.loc[df[cr]==region]
                numInfected = np.sum(gg['Confirmed'])
                numDeaths = np.sum(gg['Deaths'])
                infectedAll.append(numInfected)
                deathsAll.append(numDeaths)
                daysAll.append(k)
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
        daysAll = np.array(daysAll)
        infectedAll = np.array(infectedAll)
        deathsAll = np.array(deathsAll)
        printFlag = False
        print('Days, Infected, Deaths')
        for k in range(len(infectedAll)):
            if infectedAll[k] > 0:
                printFlag = True
            if printFlag:
                print("%d, %d, %d" % (daysAll[k],infectedAll[k],deathsAll[k]))
                # print(daysAll[k],', ',infectedAll[k],', ',deathsAll[k])
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