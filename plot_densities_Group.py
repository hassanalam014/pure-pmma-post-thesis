# Date: April 2017
#
# Description: The purpose of this file is to plot Polystyrene (PS) density information based on experiment and theory for comparison.
#

import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from matplotlib.ticker import AutoMinorLocator
from all_p_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from findVectors import findVectors
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity

#Setting font size
axis_size = 20
title_size = 20
size = 14
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.pdf'
img_dpi = None
output_folder = 'plot_density'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Defining linetype
T402_line = '-'
T452_line = '--'
T503_line = ':'

#General line properties.
linewidth = 3
markersize = 8

#Setting which set of parameters to use for calculation.
param_set = 'Self'

if param_set == 'GLi':
	title = 'Density of PS using GLi 2007 SL Parameters'
	Pstar = GLi_Pstar
	Tstar = GLi_Tstar
	Rstar = GLi_Rstar
elif param_set == 'Self':
	title = 'Density of PS using our own SL Parameters'
	Pstar = Self_Pstar
	Tstar = Self_Tstar
	Rstar = Self_Rstar

#Initializing the array of densities.
R0 = npy.linspace(0.01,0.99*Rstar,300)

gamma,vh,epsilon = calculateNewMolecularParameters(Pstar,Tstar,Rstar,M0[0])
vh = vh/NA
epsilon = epsilon/NA
print('The molecular parameters are: gamma = {}, vh = {}, and epsilon = {}.'.format(gamma,vh,epsilon))

Pmin = min(P0)
Pmax = max(P0)
Tmin = min(T0)
Tmax = max(T0)
print('The pressure range is {}-{}MPa and the temperature range is {}-{}K.'.format(Pmin,Pmax,Tmin,Tmax))

#==============================================================================================================
#Calculating Isotherms.
#==============================================================================================================

temp = ['402K','412K','422K','432K','442K','452K','463K','473K','482K','492K','503K','513K','524K']

for i in range(0,len(temp)):
	exec "result = calculatePressure(T0_%s[0],R0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (temp[i],temp[i])
	exec "T%s_P = result[0]" % (temp[i])
	exec "vector_%s = findVectors(T%s_P,R0,P0_%s,R0_%s)" % (temp[i],temp[i],temp[i],temp[i])

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#P versus R plots.
figPUREPS=plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()
plt.plot(T402K_P,R0,'k',lw=linewidth,ls=T402_line,label='PS at 402K theory')
plt.plot(T452K_P,R0,'k',lw=linewidth,ls=T452_line,label='452K theory')
plt.plot(T503K_P,R0,'k',lw=linewidth,ls=T503_line,label='453K theory')
#plt.plot(T524K_P,R0,'y')
plt.plot(P0_402K,R0_402K,'ok',ms=markersize,label='402K experiment')
plt.plot(P0_452K,R0_452K,'^k',ms=markersize,label='452K experiment')
plt.plot(P0_503K,R0_503K,'sk',ms=markersize,label='503K experiment')
#plt.plot(P0_524K,R0_524K,'*y')
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Density $\rho$ ($g/cm^3$)',fontsize=axis_size)
plt.axis([0,210,0.9,1.11])
plt.legend(loc=4,fontsize=size,numpoints=1)
#minorLocator = AutoMinorLocator()
#ax.xaxis.set_minor_locator(minorLocator)
#plt.tick_params(which='both', width=1)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4)
#minorLocator = AutoMinorLocator()
#ax.yaxis.set_minor_locator(minorLocator)
#plt.tick_params(which='both', width=1)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4)
figPUREPS.savefig('./'+output_folder+r'\pure_PS_density'+img_extension,dpi=img_dpi)

plt.show()
