# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:55:50 2015

@author: Xiaoping
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from pylab import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas 
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

ConcentrationH2O = 3 # percent, absolute
ConcentrationCH4 = 10 #ppm
PlotRange = [3.3, 3.33]
LrtzFilter = [3.315, 0.001/2]
N = 1000

def IntrapolateData(data, PlotRange, N):
    xi = linspace(min(PlotRange), max(PlotRange), N);
    x = data[:,0]
    y = data[:,1]
    order = np.argsort(x)
    ius = InterpolatedUnivariateSpline(x[order], y[order]);
    yi = ius(xi);
    return [xi, yi];
    
def lorentzian_profile(kappa, gamma, kappa0, height):
    L = height * gamma**2 / ((kappa - kappa0)**2 + gamma**2)
    return(L);

def createFilter(x,lorentzianPara):
    kappa0 = lorentzianPara[0];
    gamma = lorentzianPara[1];
    y = np.zeros(len(x));
    for i in range(0, len(x)-1):
        y[i] = lorentzian_profile(x[i],gamma,kappa0, height = 1)
    return y;

def FilteredAbsorption(y, yFilter):
    yi = np.zeros(len(y));
    for i in range(0, len(y)-1):
        yi[i] = y[i] * yFilter[i]
    return yi;
    
DataCH4 = np.loadtxt('./data/CH4.dat');
DataH2O = np.loadtxt('./data/H2O.dat');
[xCH4,yCH4] = IntrapolateData(DataCH4, PlotRange, N);
[xH2O,yH2O] = IntrapolateData(DataH2O, PlotRange, N);
yH2O = yH2O * ConcentrationH2O/3;
yCH4 = yCH4 * ConcentrationCH4/100;

yFilter = createFilter(xH2O,LrtzFilter)
yCH4f = FilteredAbsorption(yCH4, yFilter);
yH2Of = FilteredAbsorption(yH2O, yFilter);
print(['CH4/H2O ratio',sum(yCH4f)/sum(yH2Of)])
print(['Signal Level',sum(yCH4f*yFilter)/sum(yFilter)])
print('noise level depends on electronic noise and shot noise, please calculate number of photons')

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(xCH4,yCH4,label = 'CH4')
ax1.plot(xH2O,yH2O,label = 'H2O')
ax1.plot(xH2O,yFilter*max(yH2O),label = 'Filter')
ax1.legend(bbox_to_anchor=(1, 1))

ax2 = fig.add_subplot(212)
ax2.plot(xCH4,yCH4f,label = 'CH4')
ax2.plot(xH2O,yH2Of,label = 'H2O')
ax2.legend(bbox_to_anchor=(1, 1))

plt.show()
