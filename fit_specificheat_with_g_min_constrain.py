from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from loadSpecificHeatExperimentalData import *
from loadExperimentalData import Pc0,Tc0,Rc0
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


def OwnCriteria1(Rratio,epsilon_2,x,Vratio,Tg_atm,dP_dT_atm,Mg,Pstar,Tstar,Rstar):  
		
	# for key,value in kwargs.items():
	# 	exec "%s=%s" % (key,value)

	M=Mg
	T=Tg_atm
	P=P_atm
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	Pratio=Tratio/Vratio

	# Rratio = Symbol('Rratio')
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
	# answer= solve(Own_Criteria_1, Rratio)
	# print 'Rratio is:', answer

	res=Own_Criteria_1

	return res

def Rratio_From_Own_Criteria_1(**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	# Tstarstar=epsilon_2/kB
	# Tratio=Tstarstar/Tstar
	
	# print 'The values are:'
	# print P
	# print T
	# print R
	# print M
	# print Pstar
	# print Tstar
	# print Rstar
	# print dP_dT
	# print Pratio
	# print Tratio
	# print Vratio
	# print Ptilde
	# print Ttilde
	# print Rtilde
	# print r
	# print F
	print 'Problematic Value is:',epsilon_2
	# print Rratio
	Rratio = bisect(OwnCriteria1,0.3,10,args=(epsilon_2,x,Vratio,Tg_atm,dP_dT_atm,Mg,Pstar,Tstar,Rstar))
	# print 'Rratio is', Rratio

	return Rratio


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
		# print 'epsilon_2 is:',epsilon_2
		# print 'Rratio is:',Rratio

	if 'Fit_below_Tg' in fit_type:
		C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
		C3=A+B*T
		C=C1+C3
		# print 'A is:',A
		# print 'B is:',B

	return	C

def specificHeatResidualArray(params,C,P,T,R,M,I,fit_type):
	
	x = params['x'].value
	epsilon_2 = params['epsilon_2'].value
	A = params['A'].value
	B = params['B'].value
	Vratio = params['Vratio'].value
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	Tg_atm = params['Tg_atm'].value
	dP_dT_atm = params['dP_dT_atm'].value
	Mg=M[0]

	if 'Fit_below_Tg' in fit_type:
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}
		print 'A is:',A
		print 'B is:',B
		print '----------------------'

	if 'Fit_above_Tg' in fit_type:
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Mg':Mg,'epsilon_2':epsilon_2,'x':x,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}
		Rratio=Rratio_From_Own_Criteria_1(**kwargs)
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'Rratio':Rratio,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}
		print 'epsilon_2 is:',epsilon_2
		print 'Rratio is:',Rratio
		print 'x is:',x
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
# dPtilde_dT=dP_dT_atm/Pstar
# dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar
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
params_below_Tg.add_many(('x',			0.9,			    False,	0.0,	None,	None),
				(		'epsilon_2',	7000.0,				False,	0,		None,	None),
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

A=0.35758515	#PVME=0.35758515		#PMMA 140kilo=3.6870e-11				#PS=3.6870e-11
B=0.00280566	#PVME=0.00280566		#PMMA 140kilo=0.00390689				#PS=0.00372975

x= npy.linspace(0.10,0.60,20)
epsilon_2=npy.zeros(len(x))
Rratio=npy.zeros(len(x))
eps_help=npy.zeros(len(x)+1)
eps_help[0]=9000.0
eps_help_min=npy.zeros(len(x)+1)
eps_help_min[0]=eps_help[0]/2

for i in range(0,len(x)):

	print 'Now, Fitting Above Glass Transition'
	print 'Program is iterating for the cycle number = ',i+1,' with x= ', x[i]
	
	#Fitting above to the base curve above glass transition:
	params_above_Tg = Parameters()
	#The following code sets up the model's parameters. It includes both fitting parameters and parameters thabove will remain fixed
	#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
	#						(Name,			Value,		        Vary?,	Min,	Max,		Expr)
	params_above_Tg.add_many(('x',			x[i],			    False,	0,		None,		None),
					(		'epsilon_2',	eps_help[i],		True,	2000,	None,		None),
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
		
	if 'epsilon_2' in fit.params:

		epsilon_2[i] = fit.params['epsilon_2'].value
		x[i] = fit.params['x'].value
		Rratio[i]=bisect(OwnCriteria1,0.3,10,args=(epsilon_2[i],x[i],Vratio,Tg_atm,dP_dT_atm,Mg,Pstar,Tstar,Rstar))
		eps_help[i+1] = epsilon_2[i]
		eps_help_min[i+1] = epsilon_2[i]/2
		print 'epsilon_2 minimum value was:', eps_help_min[i]
		print 'Whereas epsilon_2 iterated value is:',epsilon_2[i]
		#kwargs = {'A':A,'B':B}

Rratio_min=min(Rratio)
index_min=npy.argmin(Rratio)
epsilon_2_min=epsilon_2[index_min]
x_min=x[index_min]

print Rratio_min
print epsilon_2_min
print x_min

print Rratio
print epsilon_2
# print eps_help
print x
