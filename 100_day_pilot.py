### 100-day pilot

##packages
import math
import numpy as np 
import matplotlib.pyplot as pl 
from scipy.integrate import solve_ivp

##initial states, time span, and parameters
timespan = 100
x0 = [1, 0, 0, 0] + [0]*4*timespan
outT = np.linspace(0, timespan, (timespan+1)*10)
parameters = [0.25, 0.15]

##set up the model
def quarantine(t, x, param):
    #assign state variables
    I1 = x[0]
    I2 = x[1]
    I3 = x[2]
    R  = x[3] 
    #assign parameters
    progressRate = parameters[0]
    traceRate    = parameters[1] 
    #equations
    dxdt = np.zeros(len(x))
    dxdt[0] = - progressRate*I1 - traceRate*I1 
    dxdt[1] = + progressRate*I1 - progressRate*I2 - traceRate*I2 
    dxdt[2] = + progressRate*I2 - progressRate*I3 - traceRate*I3
    dxdt[3] = + progressRate*I3 - traceRate*R
    #for loop
    date = math.ceil(t) #0.1->1, true date
    for i in range(date):
        DAY = i + 1 #date corresponding to set of equations numbers 
        if DAY < date:
            if (date - DAY) <= 15: #escape in the middle of the quarantine
                escapeRate = (date - DAY)/10 #escape rate increases across time
                QI1 = x[DAY*4]
                QI2 = x[DAY*4+1]
                QI3 = x[DAY*4+2]
                QR  = x[DAY*4+3]
                dxdt[DAY*4]   =                    - progressRate*QI1 - escapeRate*QI1 
                dxdt[DAY*4+1] = + progressRate*QI1 - progressRate*QI2 - escapeRate*QI2 
                dxdt[DAY*4+2] = + progressRate*QI2 - progressRate*QI3 - escapeRate*QI3
                dxdt[DAY*4+3] = + progressRate*QI3                    - escapeRate*QR

                dxdt[0] = dxdt[0] + escapeRate*QI1 
                dxdt[1] = dxdt[1] + escapeRate*QI2 
                dxdt[2] = dxdt[2] + escapeRate*QI3 
                dxdt[3] = dxdt[3] + escapeRate*QR 
            elif (date - DAY) == 16: #all stop quarantine after 15 days
                x[0] = x[0] + x[DAY*4]
                x[1] = x[1] + x[DAY*4+1]
                x[2] = x[2] + x[DAY*4+2]
                x[3] = x[3] + x[DAY*4+3]
                x[DAY*4] = 0 
                x[DAY*4+1] = 0
                x[DAY*4+2] = 0
                x[DAY*4+3] = 0
                dxdt[DAY*4] = 0
                dxdt[DAY*4+1] = 0
                dxdt[DAY*4+2] = 0
                dxdt[DAY*4+3] = 0
            else: #no change after 15 days 
                dxdt[DAY*4] = 0
                dxdt[DAY*4+1] = 0
                dxdt[DAY*4+2] = 0
                dxdt[DAY*4+3] = 0         
        elif DAY == date: #go to quarantine if being traced 
            dxdt[DAY*4]   = + traceRate*I1
            dxdt[DAY*4+1] = + traceRate*I2 
            dxdt[DAY*4+2] = + traceRate*I3 
            dxdt[DAY*4+3] = + traceRate*R
    return dxdt

##simulate the model
timerange = (0, timespan)
out = solve_ivp(quarantine, timerange, x0, args=(parameters, ), dense_output=True)

##figures
x = out.sol(outT).T
I1 = x[:, 0]
I2 = x[:, 1]
I3 = x[:, 2]
R  = x[:, 3]

pl.figure(1)
pl.plot(outT, I1, label='I1')
pl.plot(outT, I2, label='I2')
pl.plot(outT, I3, label='I3')
pl.plot(outT, R,  label='R')
pl.legend(loc='upper right')
pl.show()

DAY1QI1 = x[:, 4] #those traced on day 1
DAY1QI2 = x[:, 5]
DAY1QI3 = x[:, 6]
DAY1QR  = x[:, 7]

pl.figure(2)
pl.plot(outT, DAY1QI1, label='DAY1QI1')
pl.plot(outT, DAY1QI2, label='DAY1QI2')
pl.plot(outT, DAY1QI3, label='DAY1QI3')
pl.plot(outT, DAY1QR,  label='DAY1QR')
pl.legend(loc='upper right')
pl.show()

