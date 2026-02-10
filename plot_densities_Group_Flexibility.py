# Date: April 2017
#
# Description: The purpose of this file is to plot Polystyrene (PS) density information based on experiment and theory for comparison.
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
#from matplotlib.ticker import AutoMinorLocator
from all_p_params import *
from loadSpecificHeatExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
#from findVectors import findVectors
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
from wrapperFlexibilityFunctions import calculateSpecificHeat

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
M179_line = '-'
M240_line = '--'
M258_line = ':'

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

Pstarstar=Self_Pstarstar
Tstarstar=Self_Tstarstar
Rstarstar=Self_Rstarstar

#Initializing the array of densities.
T0 = npy.linspace(290,500,300)
#R0 = npy.zeros(len(T0))

# gamma,vh,epsilon = calculateNewMolecularParameters(Pstar,Tstar,Rstar,M0[0])
# vh = vh/NA
# epsilon = epsilon/NA
# print('The molecular parameters are: gamma = {}, vh = {}, and epsilon = {}.'.format(gamma,vh,epsilon))

Pmin = min(P0)
Pmax = max(P0)
Tmin = min(T0)
Tmax = max(T0)
print('The pressure range is {}-{}MPa and the temperature range is {}-{}K.'.format(Pmin,Pmax,Tmin,Tmax))

#==============================================================================================================
#Calculating Cp for various Molar Mass chains:
#==============================================================================================================

molarmass = ['179K','240K','258K']							#K=kilo; not Kelvin

for i in range(0,len(molarmass)):
	exec "R0=calculateDensity(P0_%s[0],T0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (molarmass[i],molarmass[i])
	exec "R0=R0[1]"	#Because caclulateDensity returns nested list whose 2nd row is list of R0 values
	exec "result = calculateSpecificHeat(P0_%s[0],T0,R0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)" % (molarmass[i],molarmass[i])
	exec "M%s_C = result[0]" % (molarmass[i])
	#exec "vector_%s = findVectors(T%s_P,R0,P0_%s,R0_%s)" % (temp[i],temp[i],temp[i],temp[i])

arrow_ls = 'dashdot'
show_arrows = True
#print M179K_C
#==================================================================================
#P versus R plots.
figPUREPS=plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()
plt.plot(T0,M179K_C,'k',lw=linewidth,ls=M179_line,label='PS at 179kilo theory')
plt.plot(T0,M240K_C,'k',lw=linewidth,ls=M240_line,label='240kilo theory')
plt.plot(T0,M258K_C,'k',lw=linewidth,ls=M258_line,label='258kilo theory')
#plt.plot(T524K_P,R0,'y')
plt.plot(T0_179K,C0_179K,'ok',ms=markersize,label='179kilo experiment')
plt.plot(T0_240K,C0_240K,'^k',ms=markersize,label='240kilo experiment')
plt.plot(T0_258K,C0_258K,'sk',ms=markersize,label='258kilo experiment')
#plt.plot(P0_524K,R0_524K,'*y')
plt.xlabel('Temperature T (K)',fontsize=axis_size)
plt.ylabel(r'Specific Heat $\rho$ ($J/g.K$)',fontsize=axis_size)
#plt.axis([290,500,0,3])
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
figPUREPS.savefig('./'+output_folder+r'\pure_PS_specificHeat'+img_extension,dpi=img_dpi)

plt.show()
