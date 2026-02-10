# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
# from loadSpecificHeatExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from findVectors import findVectors
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
from isListOrNpyArray import *
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
import winsound
import pyttsx		#Text to speech
from lazyme.string import palette, highlighter, formatter, color_print		#For command text color and underline etc.
# color_print('abc', color='red', underline=True, bold=True, highlight='white', faint=True, reverse=True)

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def glassTransition_for_Rratio(Rratio,P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
	res=Own_Criteria_1

	return res

def Rrat(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	for i in range(0,100,1):
		
		Rratio=0.0
		try:
			Rratio = bisect(glassTransition_for_Rratio,i,i+1,args=(P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar))
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

def glassTransitionCriteria(T,P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar

	# MY Theory:
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
	res=Own_Criteria_1

	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg

def ResidualArray(params,P,Tg):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	M = params['M'].value
	epsilon_2 = params['epsilon_2'].value
	Vratio = params['Vratio'].value
	x = params['x'].value
	Tg_atm = params['Tg_atm'].value

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	Rg_atm=density(P_atm,Tg_atm,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Rratio=Rrat(P_atm,Tg_atm,Rg_atm,M,x=x,epsilon_2=epsilon_2,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio,'x':x,'Tg_atm':Tg_atm}
	
	print 'epsilon_2=', epsilon_2
	print 'Rratio	=', Rratio
	print 'x	=', x
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],**kwargs)
		residual[j] = abs((Tg[j]-Tg_calculated))
	
	return residual

# engine = pyttsx.init()	## Text to Speech
# engine.say('The program has started')
# engine.runAndWait()

P = P_atm
T=Tg_atm
M=M_infinity
# R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Rg_atm=density(P_atm,Tg_atm,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Vratio=1.0
r = (Pstar*M)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

T_last=Tg_exp[-1]		#Picking last experimental value
P_last=Pg_exp[-1]		#Picking last experimental value

# Ptilde=P/Pstar
# Ttilde=T/Tstar
# Rtilde=R/Rstar
# dPtilde_dT=dP_dT_atm/Pstar
# dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar

######################################################################################################

######################################################################################################

r = (Pstar*M)/(kB*Tstar*Rstar)

R_atm=density(P_atm,Tg_atm,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Ptilde_atm=P_atm/Pstar
Ttilde_atm=Tg_atm/Tstar
Rtilde_atm=R_atm/Rstar

R_last=density(P_last,T_last,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Ptilde_last=P_last/Pstar
Ttilde_last=T_last/Tstar
Rtilde_last=R_last/Rstar

x_min_before=0.32
x_min_after=0.32
x= npy.linspace(x_min_before,x_min_after,1)
Rratio_array=npy.zeros(len(x))
epsilon_2_array=npy.zeros(len(x))


for i in range(0,len(x)):

	#Simultaneous Equation Solver
	mpmath.mp.dps = 15
	Vratio=1.0
	Rratio = Symbol('Rratio')
	epsilon_2 = Symbol('epsilon_2')

	Own_Criteria_1_atm=(Pstar/(Rstar*Tstar))*(-((1-Rtilde_atm)*(ln(1-Rtilde_atm))/Rtilde_atm)-((ln(Rtilde_atm))/r)+((1/Ttilde_atm)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*Tg_atm)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*Tg_atm))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*Tg_atm))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))
	Own_Criteria_1_last=(Pstar/(Rstar*Tstar))*(-((1-Rtilde_last)*(ln(1-Rtilde_last))/Rtilde_last)-((ln(Rtilde_last))/r)+((1/Ttilde_last)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T_last)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T_last))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T_last))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))

	answer=nsolve((Own_Criteria_1_atm, Own_Criteria_1_last), (Rratio, epsilon_2), (1.00, 7483.0),verify=False)
	print 'last successful value of x is:', x[i]
	Rratio_array[i]=answer[0]
	epsilon_2_array[i]=answer[1]

print Rratio_array
print epsilon_2_array

Tstarstar=npy.zeros(len(epsilon_2_array))
Tratio=npy.zeros(len(epsilon_2_array))

for i in range(0,len(epsilon_2_array)):

	Tstarstar[i]=epsilon_2_array[i]/kB
	Tratio[i]=Tstarstar[i]/Tstar

#Initializing the array of densities.
P = npy.linspace(0.101325,700,20)
R=npy.zeros(len(P))
Tg_calculated=npy.zeros(len(P))


# for i in range(0,len(P0)):
# 	R0[i]=density(P0[i],T0[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)	
	
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
output_folder = 'plot_Rratio_min_Trick'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#General line properties.
linewidth = 1
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

###################################################################################
###    Tg VERSUS P Plots  #### 
###################################################################################
plot_color=['r','k','g','m']#,'y','r','k']

for j in range(0,len(epsilon_2_array)):
	for i in range(0,len(P)):
		Tg_calculated[i]=glassTemp(P[i],M=M,x=x[j],Rratio=Rratio_array[j],Tratio=Tratio[j],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
		print P[i]

	exec "plt.plot(P,Tg_calculated,'k',color='%s',lw=linewidth,ls='-',label='x= %s')" %(plot_color[j],x[j])

# plt.plot(P_line,T_line,'k',color='r',lw=linewidth,ls='-',label='Pure PMMA Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Exp. Data')
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
# plt.axis([0.0,300,300,500])
#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)
#############################################################################################################

#############################################################################################################
# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
# plt.axis([0,165,370,430])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()
