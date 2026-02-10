# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
# from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
# from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
# from wrapperFlexibilityFunctions import calculateSpecificHeat
# from isListOrNpyArray import *
from Parameters_of_Different_Polymers import *
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
from optimizeResidualFunctions import pureEOSResidual,pureChemicalPotentialResidual
from loadSpecificHeatExperimentalData import *
from sympy import Symbol, nsolve
import sympy
import mpmath
import cmath

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

### Programs My Theory

def glassTransitionDerivative_for_eps2(epsilon_2,P,T,R,M,x,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))

	# print dS_dT_Own_Criteria_1
	res=Own_Criteria_1

	return res

def eps_2(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
		
	for i in range(0,50000,500):
		
		epsilon_2=0.0
		try:
			epsilon_2 = bisect(glassTransitionDerivative_for_eps2,i,i+500,args=(P,T,R,M,x,Rratio,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get epsilon_2")
			pass
		if epsilon_2!=0.0:
			print 'Hurry! epsilon_2_dependent is:', epsilon_2, 'Rratio_independent is:', Rratio
			break
	
	if epsilon_2==0.0:
		print 'Program Failed to get value of epsilon_2 at Rratio=', Rratio
		# epsilon_2=50000

	return epsilon_2


def glassTransitionDerivative_for_Rratio(Rratio,P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
	# print dS_dT_Own_Criteria_1
	res=Own_Criteria_1

	return res

def Rrat(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	for i in range(0,100,1):
		
		Rratio=0.0
		try:
			Rratio = bisect(glassTransitionDerivative_for_Rratio,i,i+1,args=(P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get Rratio")
			pass
		if Rratio!=0.0:
			print 'Hurry! Rratio_dependent is:', Rratio, 'epsilon_2_independent is:', epsilon_2
			break
	
	if Rratio==0.0:
		print 'Program Failed to get value of Rratio at epsilon_2=',epsilon_2
		# Rratio=50000

	return Rratio


P = P_atm
T=Tg_atm
M=M_infinity
R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

Ptilde=P/Pstar
Ttilde=T/Tstar
Rtilde=R/Rstar
dPtilde_dT=dP_dT_atm/Pstar
dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar

Vratio=1.0
#####################################################################################################

######################################################################################################
epsilon_2_independent = npy.linspace(0,10000,20)
Rratio_dependent=npy.zeros(len(epsilon_2_independent))
Rratio_independent = npy.linspace(0.0,100,20)
epsilon_2_dependent=npy.zeros(len(Rratio_independent))

x=0.32

for i in range(0,len(Rratio_independent)):
	epsilon_2_dependent[i]=eps_2(P,T,R,M,x=x,Rratio=Rratio_independent[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	if epsilon_2_dependent[i]==0:
		Rratio_independent[i]=0

for i in range(0,len(epsilon_2_independent)):
	Rratio_dependent[i]=Rrat(P,T,R,M,x=x,epsilon_2=epsilon_2_independent[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	if Rratio_dependent[i]==0:
		epsilon_2_independent[i]=0

print Rratio_dependent
print epsilon_2_independent

print epsilon_2_dependent
print Rratio_independent
###############################################################################

###############################################################################

#Setting font size
axis_size = 20
title_size = 20
size = 14
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'plots'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#General line properties.
linewidth = 1
markersize = 3

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

# plt.plot(epsilon_2_independent,Rratio_dependent,'k',color='g',lw=linewidth,ls='-',label='eps_indep.')
# plt.plot(epsilon_2_dependent,Rratio_independent,'k',color='r',lw=linewidth,ls='-',label='Rratio_indep')
plt.plot(epsilon_2_independent,Rratio_dependent,'sk',color='g',ms=markersize,label='eps_indep.')
plt.plot(epsilon_2_dependent,Rratio_independent,'sk',color='r',ms=markersize,label='Rratio_indep')

# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')

plt.xlabel('epsilon_2',fontsize=axis_size)
plt.ylabel(r'Rratio',fontsize=axis_size)
# plt.axis([0,50,350,400])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()
