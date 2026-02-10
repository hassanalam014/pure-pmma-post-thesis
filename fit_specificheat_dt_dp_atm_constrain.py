from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from loadSpecificHeatExperimentalData import *
# from loadExperimentalData import Pc0,Tc0,Rc0
# from calculatePureVariables import density
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
# from calculateSimpleFlexibilityResidual import calculatePureSpecificHeatResidual,calculatePureSpecificHeatResidual2
# from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
# from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
# from isListOrNpyArray import *
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
from optimizeResidualFunctions import pureEOSResidual,pureChemicalPotentialResidual
from Parameters_of_Different_Polymers import *


def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def OwnCriteria1(x,P,T,R,M,Pratio,Tratio,Rratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))

	res=Own_Criteria_1

	return res

def x_atm(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	epsilon_2 = bisect(OwnCriteria1,0.01,0.99,args=(P,T,R,M,Pratio,Tratio,Rratio,Pstar,Tstar,Rstar))
	
	return epsilon_2



def glassTransitionDerivative(epsilon_2,P,T,R,M,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	dS_dT_Own_Criteria_1=(Vratio*Ttilde*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((1/Rtilde)*(dPtilde_dTtilde)))/(Rtilde*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))+(((Vratio*epsilon_2/(kB*T))**2)*((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))*(1-((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))))
	# print 'dS_dT_Own_Criteria_1 sign at epsilon_2=',epsilon_2, ' is:', dS_dT_Own_Criteria_1

	res=dS_dT_Own_Criteria_1

	return res

def eps_2(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	
	for i in range(0,50000,500):
		
		epsilon_2=0.0
		try:
			epsilon_2 = bisect(glassTransitionDerivative,i,i+500,args=(P,T,R,M,Rratio,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get epsilon_2")
			pass
		if epsilon_2!=0.0:
			break
	
	if epsilon_2==0.0:
		print 'Program Failed to get value of epsilon_2'
		epsilon_2=50000

	return epsilon_2


def specificHeat(P,T,R,M,fit_type,**kwargs):    
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	if 'Fit_above_Tg' in fit_type:
		
		Tstarstar=epsilon_2/kB
		Tratio=Tstarstar/Tstar

		Pratio=Tratio/Vratio

		Tstarstar=Tratio*Tstar
		Pstarstar=Pratio*Pstar
		Rstarstar=Rratio*Rstar
		
		C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
		C2=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
		C3=A+B*T
		C=C1+C2+C3

	if 'Fit_below_Tg' in fit_type:
		C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
		C3=A+B*T
		C=C1+C3

	return	C

def specificHeatResidualArray(params,C,P,T,R,M,I,fit_type):
	
	A = params['A'].value
	B = params['B'].value
	Rratio = params['Rratio'].value
	Vratio = params['Vratio'].value
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	Tg_atm = params['Tg_atm'].value
	dP_dT_atm = params['dP_dT_atm'].value
	Mg=M[0]
	Rg_atm=density(P_atm,Tg_atm,Mg,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	epsilon_2=eps_2(P_atm,Tg_atm,Rg_atm,Mg,Rratio=Rratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'Rratio':Rratio,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}

	if 'Fit_below_Tg' in fit_type:
		print 'A is:',A
		print 'B is:',B
		print '----------------------'

	if 'Fit_above_Tg' in fit_type:
		print 'epsilon_2 is:',epsilon_2
		print 'Rratio is:',Rratio
		print '-----------------------'

	residual=npy.zeros(len(C))

	for j in range(0,len(C)):
		C_calculated = specificHeat(P[j],T[j],R[j],M[j],fit_type,**kwargs)
		residual[j] = (C[j]-C_calculated)

	return residual

# P = P_atm
# M=M_infinity
# R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
# r = (Pstar*M)/(kB*Tstar*Rstar)
T=Tg_atm
dP_dT_atm=1/dTg_dP_atm
Vratio=1.0
Mg=M0_complete_Tg[0]

# Ptilde=P/Pstar
# Ttilde=T/Tstar
# Rtilde=R/Rstar
# # vtilde=1/Rtilde	
dPtilde_dT=dP_dT_atm/Pstar
dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar
##########################################################################################################

##########################################################################################################

R0_complete_Tg=npy.zeros(len(C0_complete_Tg))
for j in range(0,len(C0_complete_Tg)):
	R0_complete_Tg[j]=density(P0_complete_Tg[j],T0_complete_Tg[j],M0_complete_Tg[j],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

R0_below_Tg=npy.zeros(len(T0_below_Tg))
for j in range(0,len(T0_below_Tg)):
	R0_below_Tg[j]=density(P0_below_Tg[j],T0_below_Tg[j],M0_below_Tg[j],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

R0_at_Tg=npy.zeros(len(T0_at_Tg))
for j in range(0,len(T0_at_Tg)):
	R0_at_Tg[j]=density(P0_at_Tg[j],T0_at_Tg[j],M0_at_Tg[j],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

R0_above_Tg=npy.zeros(len(T0_above_Tg))
for j in range(0,len(T0_above_Tg)):
	R0_above_Tg[j]=density(P0_above_Tg[j],T0_above_Tg[j],M0_above_Tg[j],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

R0_excluding_Tg=npy.zeros(len(T0_excluding_Tg))
for j in range(0,len(T0_excluding_Tg)):
	R0_excluding_Tg[j]=density(P0_excluding_Tg[j],T0_excluding_Tg[j],M0_excluding_Tg[j],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)


#################################################################
'''
#################################################################

print T0_below_Tg
print 'First Fitting Base Curve'

#Fitting Data to the base curve below glass transition:
params_below_Tg = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#						(Name,			Value,		        Vary?,	Min,	Max,	Expr)
params_below_Tg.add_many(('Rratio',		1.64,			    False,	0.0,	None,	None),
				(		'A',			0.1,				True,	0,		None,	None),
				(		'B',			0.00258531,			True,	0,		None,	None),
				(		'Vratio',		Vratio,				False,	0,		None,	None),
				(		'Pstar',		Pstar,				False,	0,		None,	None),
				(		'Tstar',		Tstar,				False,	0,		None,	None),
				(		'Rstar',		Rstar,				False,	0,		None,	None),
				(		'Tg_atm',		Tg_atm,		        False,	0,		None,	None),
				(		'dP_dT_atm',	dP_dT_atm,	        False,	0,		None,	None))

#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
#These need to be added to the experimental datapints to find the fitted specific heat.
fit = minimize(specificHeatResidualArray,params_below_Tg,args=(C0_below_Tg,P0_below_Tg,T0_below_Tg,R0_below_Tg,M0_below_Tg,I0_below_Tg,'Fit_below_Tg'))

#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
report_fit(fit.params)

if 'A' in fit.params and 'B' in fit.params:
	A = fit.params['A'].value
	B = fit.params['B'].value
	#kwargs = {'A':A,'B':B}

######################################################################
'''
######################################################################

A=3.6870e-11	#PVME=0.35758515		#PMMA 140kilo=3.6870e-11	#PMMA 44kilo = 3.6870e-11		#PS=3.6870e-11
B=0.00372975	#PVME=0.00280566		#PMMA 140kilo=0.00390689	#PMMA 44kilo = 0.00551656		#PS=0.00372975

#Fitting above to the base curve above glass transition:
params_above_Tg = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters thabove will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#						(Name,			Value,		        Vary?,	Min,	Max,		Expr)
params_above_Tg.add_many(('Rratio',		1.64,				True,	0.0,	10000,		None),
				(		'A',			A,					False,	0,		None,		None),
				(		'B',			B,					False,	0,		None,		None),
				(		'Vratio',		Vratio,				False,	0,		None,		None),
				(		'Pstar',		Pstar,				False,	0,		None,		None),
				(		'Tstar',		Tstar,				False,	0,		None,		None),
				(		'Rstar',		Rstar,				False,	0,		None,		None),
				(		'Tg_atm',		Tg_atm,		        False,	0,		None,		None),
				(		'dP_dT_atm',	dP_dT_atm,	        False,	0,		None,		None))

#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
#These need to be added to the experimental daboveapints to find the fitted specific heabove.
fit = minimize(specificHeatResidualArray,params_above_Tg,args=(C0_above_Tg,P0_above_Tg,T0_above_Tg,R0_above_Tg,M0_above_Tg,I0_above_Tg,'Fit_above_Tg'))

#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
report_fit(fit.params)
	
if 'Rratio' in fit.params:

	Rratio = fit.params['Rratio'].value
	#kwargs = {'A':A,'B':B}

Rg_atm=density(P_atm,Tg_atm,Mg,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
epsilon_2=eps_2(P_atm,Tg_atm,Rg_atm,Mg,Rratio=Rratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Pratio=Tratio/Vratio
x=x_atm(P_atm,Tg_atm,Rg_atm,Mg,Pratio=Pratio,Tratio=Tratio,Rratio=Rratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
print 'x is:', x
print 'Rratio is:', Rratio
print 'epsilon_2 is:', epsilon_2

# print C0_above_Tg
# print P0_above_Tg 
# print T0_above_Tg 
# print R0_above_Tg
# print M0_above_Tg
# print I0_above_Tg

