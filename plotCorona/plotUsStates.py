
import sys
import numpy as np
import pandas as pd
import plotCorona as pc
import infectionFitting as infect
import matplotlib.pyplot as plt

# python37.exe .\plotUsStates.py Arizona Florida Georgia Illinois

global popData
popData = None

def getUSstates():
    return ['Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado',
 'Connecticut', 'Delaware', 'Diamond Princess', 'District of Columbia',
 'Florida', 'Georgia', 'Grand Princess', 'Guam', 'Hawaii', 'Idaho', 'Illinois',
 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland',
 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana',
 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York',
 'North Carolina', 'North Dakota', 'Northern Mariana Islands', 'Ohio',
 'Oklahoma', 'Oregon', 'Pennsylvania', 'Puerto Rico', 'Rhode Island', 
 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah',
 'Vermont', 'Virgin Islands', 'Virginia', 'Washington', 'West Virginia',
 'Wisconsin', 'Wyoming']

def getPopulation(state):
    global popData
    if popData is None:
        print("Reading in Population Data")
        popData = pd.read_csv('nst-est2019-alldata.csv')
    gg = popData.loc[popData['NAME']==state]
    return float(gg['POPESTIMATE2019'])

if __name__=="__main__":
    narg = len(sys.argv)
    if narg > 1:
        states = []
        for k in range(1,narg):
            states.append(sys.argv[k])
    else:
        states = getUSstates()
    if states[-1] == 'byPopulation':
        byPop = True
        states.remove(states[-1])
    else:
        byPop = False
    fig, axes = plt.subplots(1,2,figsize=(14,6))
    for state in states:
        print("Plotting State:", state)
        if state == 'United States':
            sys.argv = ['temp','US']
            I,D,T = pc.createTimeDataCountry()
        else:
            sys.argv = ['temp','US',state]
            I,D,T = pc.createTimeDataState()

        T-=T[0]
        dI = np.zeros(len(I))
        dI[1:] = I[1:] - I[:-1]
        mdI = infect.running_mean(dI,14)

        dD = np.zeros(len(D))
        dD[1:] = D[1:] - D[:-1]
        mdD = infect.running_mean(dD,14)

        if byPop:
            statePop = getPopulation(state)
            mdI = 100.0*mdI/statePop
            mdD = 100.0*mdD/statePop
        axes[0].plot(T,mdI,label=state)
        axes[1].plot(T,mdD,label=state)

    if byPop:
        axes[0].set(xlabel="Days Since First Death",ylabel="Percent Change in Infections")
        axes[1].set(xlabel="Days Since First Death",ylabel="Percent Change in Deaths")
    else:
        axes[0].set(xlabel="Days Since First Death",ylabel="Change in Infections")
        axes[1].set(xlabel="Days Since First Death",ylabel="Change in Deaths")
    axes[0].legend(loc=2, ncol=2)
    axes[1].legend(loc=2, ncol=2)
    plt.tight_layout()
    plt.show()
