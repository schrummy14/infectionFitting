
import sys
import pyDOE
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import infectionFitting as infect

infectionName = "coronavirus"

startDate = dt.date(2020,3,1)
endDate = dt.datetime.now() + dt.timedelta(days=-1)
endDate = dt.date(endDate.year,endDate.month,endDate.day)
baseURL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/'

def getData(date):
    mdy = getDateStr(date)
    df = infect.getData(mdy)
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


def printGrowthFactor(I):
    deltaI = I[1:]-I[0:-1]
    print(deltaI[1:]/deltaI[0:-1])


def createTimeDataCountry():
    region = sys.argv[1]
    infected = []
    deaths = []
    days = []

    posRegions = getRegionsCountry()
    if region not in posRegions:
        print("Region: %s, is not a posible region..." % region)
        print("Possible Regions:\n")
        print(np.sort(posRegions))
        exit()

    k = 0
    while True:
        curDate = startDate + dt.timedelta(days=k)
        k+=1
        if curDate > endDate:
            break
        df = getData(curDate)
        cr = 'Country/Region'
        try:
            df[cr]
        except:
            cr = 'Country_Region'
        try:
            if df is None:
                continue
            if region == "all":
                break
            elif region == "world":
                infected.append(np.sum(df['Confirmed']))
                deaths.append(np.sum(df['Deaths']))
                days.append(k)
            else:
                if region == 'China':
                    gg = df.loc[df[cr]=='China']
                    if len(gg) == 0:
                        gg = gg = df.loc[df[cr]=='Mainland China']
                else:
                    gg = df.loc[df[cr]==region]
                numInfected = np.sum(gg['Confirmed'])
                numDeaths = np.sum(gg['Deaths'])
                if numInfected*numDeaths > 0:
                    infected.append(np.sum(gg['Confirmed']))
                    deaths.append(np.sum(gg['Deaths']))
                    days.append(k)
        except:
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
        infect.showDailyChange(daysAll,infectedAll,deathsAll)
        exit()

def getRegionsCountry():
    date = endDate + dt.timedelta(days=-3)
    mdy = getDateStr(date)
    try:
        df = pd.read_csv(mdy,parse_dates=[2])
    except:
        url = baseURL + mdy
        df = pd.read_csv(url,parse_dates=[2])
    cr = 'Country/Region'
    try:
        df[cr]
    except:
        cr = 'Country_Region'
    regions = df[cr].unique()
    return regions

def createTimeDataState():
    
    
    infected = []
    deaths = []
    days = []
    region = sys.argv[2]

    posRegions = getRegionsState()
    if region not in posRegions:
        print("Region: %s, is not a posible region..." % region)
        print("Possible Regions:\n")
        print(np.sort(posRegions))
        exit()
    

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
        infect.showDailyChange(daysAll,infectedAll,deathsAll)
        exit()

def getRegionsState():
    date = endDate + dt.timedelta(days=-3)
    mdy = getDateStr(date)
    url = baseURL + mdy
    try:
        df = pd.read_csv(mdy,parse_dates=[2])
    except:
        url = baseURL + mdy
        df = pd.read_csv(url,parse_dates=[2])

    cr = 'Country/Region'
    try:
        df[cr]
    except:
        cr = 'Country_Region'
    df = df[df[cr]==sys.argv[1]]
    regions = df['Province_State'].unique()
    return regions

def createTimeDataCounty():
    infected = []
    deaths = []
    days = []

    infectedAll = []
    deathsAll = []
    daysAll = []

    state = sys.argv[2]
    region = sys.argv[3]

    if not (region == 'all'):
        posRegions = getRegionsCounty()
        if region not in posRegions:
            print("Region: %s, is not a posible region..." % region)
            print("Possible Regions:\n")
            print(np.sort(posRegions))
            exit()

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
            df = df[df[cr]==state]
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
        infect.showDailyChange(daysAll,infectedAll,deathsAll)
        exit()


def getRegionsCounty():
    if not sys.argv[1] == 'US':
        print("Currently only works with the US")
        exit()
    date = endDate + dt.timedelta(days=-3)
    mdy = getDateStr(date)
    try:
        df = pd.read_csv(mdy,parse_dates=[2])
    except:
        url = baseURL + mdy
        df = pd.read_csv(url,parse_dates=[2])
    cr = 'Province_State'
    try:
        df[cr]
    except:
        cr = 'Province_State'
    df = df[df[cr]==sys.argv[2]]
    cr = 'Admin2'
    regions = df[cr].unique()
    return regions

def printRawData(I,D,T):
    data = np.transpose(np.array([T,I,D]))
    print("Days, Infected, Deaths")
    print(data)

if __name__ == "__main__":

    numArg = len(sys.argv)
    print(sys.argv)
    if numArg == 1:
        print("Error: Usage-> plotCorona.py [country] [state {optional}] [county {optional}]")
        print("               plotCorona.py [country {optional}] [state {optional}] getRegions")
        exit()
    
    if sys.argv[-1] == 'getRegions':
        if numArg == 2:
            print(np.sort(getRegionsCountry()))
        elif numArg == 3:
            print(np.sort(getRegionsState()))
        else:
            print(np.sort(getRegionsCounty()))
        exit()

    country = sys.argv[1]
    if numArg == 2:
        I,D,T = createTimeDataCountry()
    elif numArg == 3:
        I,D,T = createTimeDataState()
    else:
        I,D,T = createTimeDataCounty()
        
    
    print("\nRaw Data")
    printRawData(I,D,T)
    print("\nInfection Growth Factor")
    printGrowthFactor(I)
    print("\nDeath Growth Factor")
    printGrowthFactor(D)

    infect.showDailyChange(T,I,D)