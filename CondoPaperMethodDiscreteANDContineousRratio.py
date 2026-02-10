# Date: 2019
#
# Description: The purpose of this file is to plot Polystyrene (PS) Thermodynamics Properties
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
import winsound    # Play a beep sound
#from matplotlib.ticker import AutoMinorLocator
# from all_p_params import *
# from loadSpecificHeatExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from findVectors import findVectors
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
# from wrapperFlexibilityFunctions import calculateSpecificHeat
from isListOrNpyArray import *
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
from optimizeResidualFunctions import pureEOSResidual,pureChemicalPotentialResidual
# from loadSpecificHeatExperimentalData import *
from sympy import Symbol, nsolve
import sympy
import mpmath
from Parameters_of_Different_Polymers import *


def TgCriteriaForRratio(Rratio,P,T,M,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar

	#Condo Theory:
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))

	res=S_condo

	return res

def Rrat(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Rratio = bisect(TgCriteriaForRratio,0.4,100,args=(P,T,M,Tratio,Vratio,Pstar,Tstar,Rstar))
	return Rratio


def ResidualArrayContineous(params,P,T):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	M = params['M'].value
	epsilon_2 = params['epsilon_2'].value
	Vratio = params['Vratio'].value
	Tg_atm = params['Tg_atm'].value
	P_atm = params['P_atm'].value

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	Rratio=Rrat(P_atm,Tg_atm,M,Tratio=Tratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	print epsilon_2
	print Rratio
	# kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Vratio':Vratio}
	# Rratio=Rrat(P_atm,Tg_atm,M,**kwargs)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio}
	
	# print Rratio
	# print epsilon_2
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],M,**kwargs)
		residual[j] = abs((T[j]-Tg_calculated))/T[j]
	
	return residual


def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R


def TgCriteriaForEpsilon(epsilon_2,P,T,M,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Pratio=Tratio/Vratio
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar

	#Condo Theory:
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))

	res=S_condo

	return res

def eps(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	# epsilon_2 = bisect(TgCriteriaForEpsilon,1500,20000,args=(P,T,M,Rratio,Vratio,Pstar,Tstar,Rstar))
		
	for i in range(0,50000,500):
		
		epsilon_2=0.0
		try:
			epsilon_2 = bisect(TgCriteriaForEpsilon,i,i+500,args=(P,T,M,Rratio,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get epsilon_2")
			pass
		if epsilon_2!=0.0:
			print 'Hurry! epsilon_2_dependent is:', epsilon_2, 'Rratio_independent is:', Rratio
			break
	
	if epsilon_2==0.0:
		print 'Program Failed to get value of epsilon_2 against Rratio_independent=', Rratio
		# epsilon_2=5000

	return epsilon_2


def glassTransitionCriteria(T,P,M,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar

	#Condo Theory:
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))
	# dS_dTg_condo=(1+(ln(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(ln(1-Rtilde)+Rtilde))-(((((Rratio)*exp(-epsilon_2/(kB*T)))/(1+((Rratio)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(((1-(((Rratio)*exp(-epsilon_2/(kB*T)))/(1+((Rratio)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)
	# dS_dTg_condo_again=(((r-2)/r)*((epsilon_2/(kB*T))**2)*(((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))*(1-((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T)))))))/T)*((Rtilde)**2)*(((Ttilde/Rtilde)*((1/r)+(Rtilde/(1-Rtilde))))-2))+((((ln(1-Rtilde))/(Rtilde))+1-(1/r))*((dPtilde_dT)+((1/Tstar)*((ln(1-Rtilde))+((1-(1/r))*Rtilde)))))

	res=S_condo

	return res

def glassTemp(P,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg

def ResidualArray(P,T,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	# kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Vratio':Vratio}
	# Rratio=Rrat(P_atm,Tg_atm,M,**kwargs)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio}
	
	# print Rratio
	# print epsilon_2
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],M,**kwargs)
		residual[j] = abs((T[j]-Tg_calculated))/T[j]
	
	return residual


Program_Running_For=['BPP Hollander 15kilo']

Pick_List_Element = Program_Running_For[0]
Divide_List_Picked_Element = Pick_List_Element.split()

print(Divide_List_Picked_Element)

Polymer_Type=Divide_List_Picked_Element[0]
Reference=Divide_List_Picked_Element[1]
Polymer_Weight=Divide_List_Picked_Element[2]
# class Polymer_Type

kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

M=M_infinity
Vratio=1.0
# dP_dT=4.32  #PS MY:2.86	#PS CONDO:3.164	#PMMA MY:4.314	#PMMA CONDO Actual Value in Paper = 5.12=Wrong??
# dTg_dP_atm_condo = 0.236         		 #Ref[53] Condo Value=0.236, 
# dTg_dP_atm_extreme_ends = 0.2318         #Line passing through extreme ends of exp data
dTg_dP_atm_extreme_ends = (T_upper-Tg_exp[0])/(P_upper-Pg_exp[0])         #Line passing through extreme ends of only linear range of exp data
# dP_dT_condo=1/dTg_dP_atm_condo
dP_dT_extreme_ends=1/dTg_dP_atm_extreme_ends
# dPtilde_dT=dP_dT/Pstar

P_line = npy.linspace(0.101325,P_upper,10)
Tg_line = npy.zeros(len(P_line))

#Ideal Experimental Straight Line Data for Fitting:
for i in range(0,len(P_line)):
	Tg_line[i]=((P_line[i]-P_atm)/dP_dT_extreme_ends)+Tg_atm		#Extreme End Line Passing Line
	# Tg_line[i]=((P_line[i]-P_atm)/dP_dT_condo)+Tg_atm				#Condo Ref[53] Slope Line
	# Tg_line[i] = 0.22774*P_line[i] + 378.57256						#Free Regression Line
	# Tg_line[i] = 0.23239*P_line[i] + 377.97933						#Tg_atm Fixed Regression Line
'''
Rratio=[2.0,3.0,4.0,5.0]
epsilon_2_list = npy.zeros(len(Rratio))
residual_list=npy.zeros(len(Rratio))

for j in range(0,len(Rratio)):
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Rratio':Rratio[j],'Vratio':Vratio}
	epsilon_2_list[j]=eps(P_atm,Tg_atm,M,**kwargs)

for j in range(0,len(Rratio)):
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Rratio':Rratio[j],'Vratio':Vratio,'epsilon_2':epsilon_2_list[j],'Tg_atm':Tg_atm,'P_atm':P_atm}
	residual_array=ResidualArray(P_line,Tg_line,**kwargs)
	SumSquareError=0

	for i in residual_array:
		SumSquareError+=i*i

	residual_list[j]=SumSquareError

print Rratio
print epsilon_2_list
print residual_list

residual_min=min(residual_list)
index_min=npy.argmin(residual_list)
epsilon_2_min=epsilon_2_list[index_min]
Rratio_min=Rratio[index_min]

print Rratio_min
print epsilon_2_min

########################################################################
'''
'''
########################################################################

#Fitting Idealized Experimental Straight Line Data:
params = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#						(Name,			Value,		        Vary?,	Min,		Max,	Expr)
params.add_many((		'epsilon_2',	7000.0,		    	True,	3000.0,		20000,	None),				
				(		'Vratio',		Vratio,				False,	0.0,		1.0,	None),				
				(		'M',			M,					False,	0.0,		None,	None),				
				(		'P_atm',		P_atm,				False,	0.0,		None,	None),
				(		'Tg_atm',		Tg_atm,				False,	0.0,		None,	None),
				(		'Pstar',		Pstar,				False,	0.0,		None,	None),
				(		'Tstar',		Tstar,				False,	0.0,		None,	None),
				(		'Rstar',		Rstar,				False,	0.0,		None,	None))

#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
#These need to be added to the experimental datapints to find the fitted specific heat.
fit = minimize(ResidualArrayContineous,params,args=(P_line,Tg_line))

#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
report_fit(fit.params)

if 'epsilon_2' in fit.params:
	epsilon_2 = fit.params['epsilon_2'].value

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Vratio':Vratio}
Rratio=Rrat(P_atm,Tg_atm,M,**kwargs)

print 'Rratio is', Rratio
print 'epsilon_2 is', epsilon_2

#Play a Beep Sound
duration = 1000  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)

##########################################################################################
'''
##########################################################################################
#For Plot of Tg Versus Pressure:
#							#MyOwnCriteria_1 (Best fit)				#Condo Theory			#My Theory d2S/dT2|p=0 with v!=v_0	##My Theory dCp/dT|p=0 with v!=v_0
Rratio=3			#2.3984278705790363									#4.915257075699675 		#0.3411518312632549					#0.04394803
epsilon_2=4354.015059

#7020.401851731026			#7235.098320856251 		#10818.738735660516					#39795.7441
# Vratio=0.18740252			#0.37374119657062027								#----------------		#1.0143895679458907					#40.2622941

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar

#Initializing the array of densities.
P=npy.linspace(0.101325,500,20)
R=npy.zeros(len(P))
Tg_calculated=npy.zeros(len(P))

for i in range(0,len(P)):
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio}
	Tg_calculated[i] = glassTemp(P[i],M,**kwargs)
	print P[i]

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
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

plt.plot(P,Tg_calculated,'k',color='g',lw=linewidth,ls='-',label='Therotical Curve')
plt.plot(P_line,Tg_line,'k',color='r',lw=linewidth,ls='-',label='Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Exp. Data')

# plt.axhline(y=243.5,lw=0.5,color='k', linestyle='-.')
# plt.axvline(x=0.101325,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()
