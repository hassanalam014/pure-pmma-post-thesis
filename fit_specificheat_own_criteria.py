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

'''
def OwnCriteria1_Deriv_Equation(Rratio,Tratio,Vratio,Tg_atm,dP_dT_atm,M,Pstar,Tstar,Rstar):  
	
	M=M[0]
	r = (Pstar*M)/(kB*Tstar*Rstar)

	T=Tg_atm
	dP_dT=dP_dT_atm
	P=0.101325
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar
	epsilon_2=kB*Tstarstar

	dPtilde_dTtilde=dP_dT*(Tstar/Pstar)

	F=(Rratio*exp(-Tratio**2/(Ttilde*Pratio)))/(1+Rratio*exp(-Tratio**2/(Ttilde*Pratio)))
	
	res=((((Tratio*Ttilde)/Pratio)*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((dPtilde_dTtilde)/(Rtilde))))/(Rtilde*(((Ttilde/Rtilde)*(((Rtilde)/(1-Rtilde))+(1/r)))-2)))+((((Tratio*Tstarstar)/(Pratio*T))**2)*F*(1-F))

	return res


def Rratio_From_Own_Criteria_1(M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	
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
	print epsilon_2
	# print Rratio

	Rratio = bisect(OwnCriteria1_Deriv_Equation,0.3,10,args=(Tratio,Vratio,Tg_atm,dP_dT_atm,M,Pstar,Tstar,Rstar))
	print 'Rratio is'
	print Rratio

	return Rratio
'''

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
		# print epsilon_2
		# print Rratio

	if 'Fit_below_Tg' in fit_type:
		C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
		C3=A+B*T
		C=C1+C3
		print A
		print B

	return	C


def specificHeatResidualArray(params,C,P,T,R,M,I,fit_type):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	epsilon_2 = params['epsilon_2'].value
	Rratio = params['Rratio'].value
	Vratio = params['Vratio'].value
	A = params['A'].value
	B = params['B'].value
	Tg_atm = params['Tg_atm'].value
	dP_dT_atm = params['dP_dT_atm'].value

	print epsilon_2
	print Rratio

	if 'Fit_below_Tg' in fit_type:
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'Rratio':Rratio,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}
	
	if 'Fit_above_Tg' in fit_type:
		#Falto: Rratio=Rratio_From_Own_Criteria_1(M,**kwargs)
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'Rratio':Rratio,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}

	# print epsilon_2
	
	residual=npy.zeros(len(C))

	for j in range(0,len(C)):
		C_calculated = specificHeat(P[j],T[j],R[j],M[j],fit_type,**kwargs)
		residual[j] = (C[j]-C_calculated)

	return residual

P = P_atm
T=Tg_atm
M=M_infinity
R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

Ptilde=P/Pstar
Ttilde=T/Tstar
Rtilde=R/Rstar
# vtilde=1/Rtilde	
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

#################################################################

print T0_below_Tg
print 'First Fitting Base Curve'

#Fitting Data to the base curve below glass transition:
params_below_Tg = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#						(Name,			Value,		        Vary?,	Min,	Max,	Expr)
params_below_Tg.add_many(('A',			0.1,			    True,	0.0,	None,	None),
				(		'B',			0.00258531,			True,	0,		None,	None),
				(		'Vratio',		1.0,				False,	0,		None,	None),
				(		'epsilon_2',	7000.0,				False,	0,		None,	None),
				(		'Rratio',		3.0,				False,	0,		None,	None),
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

A=3.6870e-11		#PMMA 140kilo=3.6870e-11				#PS=3.6870e-11
B=0.00551656		#PMMA 140kilo=0.00390689				#PS=0.00372975

print 'Now, Fitting Above Glass Transition'

#Fitting above to the base curve above glass transition:
params_above_Tg = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters thabove will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#						(Name,			Value,		        Vary?,	Min,	Max,		Expr)
params_above_Tg.add_many(('A',			A,				    False,	0,		None,		None),
				(		'B',			B,					False,	0,		None,		None),
				(		'Vratio',		1.0,				False,	0,		None,		None),
				(		'epsilon_2',	6251.31324,			True,	0,		10000.00,	None),
				(		'Rratio',		1.20630412,			True,	0,		None,		None),
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
	epsilon_2 = fit.params['epsilon_2'].value
	Rratio = fit.params['Rratio'].value
	Vratio = fit.params['Vratio'].value
	# kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'Rratio':Rratio,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}
	#Rratio=Rratio_From_Own_Criteria_1(M0_above_Tg,**kwargs)
	print epsilon_2
	print Rratio
'''


'''
ANSWER FOR PS

[[Variables]]
    A:          3.687e-11 (fixed)
    B:          0.00372975 (fixed)
    Vratio:     1 (fixed)
    epsilon_2:  6480.99556 +/- 87.7958505 (1.35%) (init = 7151)
    Rratio:     1.85399320 +/- 0.03753580 (2.02%) (init = 1.65)
    Pstar:      357 (fixed)
    Tstar:      735 (fixed)
    Rstar:      1.105 (fixed)
    Tg_atm:     374 (fixed)
    dP_dT_atm:  3.164557 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(epsilon_2, Rratio) = -0.999
6480.995558620286
1.8539932000439938

PMMA 140Kilo Fiting Answer:

[[Variables]]
    A:          3.687e-11 (fixed)
    B:          0.00390689 (fixed)
    Vratio:     1 (fixed)
    epsilon_2:  6251.31324 +/- 908.242883 (14.53%) (init = 6480.99)
    Rratio:     1.20630412 +/- 0.23442955 (19.43%) (init = 1.1)
    Pstar:      503 (fixed)
    Tstar:      696 (fixed)
    Rstar:      1.269 (fixed)
    Tg_atm:     378 (fixed)
    dP_dT_atm:  4.237288 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(epsilon_2, Rratio) = -0.997
6251.3132383031225
1.2063041202404259

ANSWER of PMMA 44kilo Fit:

[[Variables]]
    A:          3.687e-11 (fixed)
    B:          0.00551656 (fixed)
    Vratio:     1 (fixed)
    epsilon_2:  7550.38882 +/- 771.896319 (10.22%) (init = 6251.313)
    Rratio:     63.5105956 +/- 30.7724911 (48.45%) (init = 1.206304)
    Pstar:      503 (fixed)
    Tstar:      696 (fixed)
    Rstar:      1.269 (fixed)
    Tg_atm:     378 (fixed)
    dP_dT_atm:  4.237288 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(epsilon_2, Rratio) =  1.000
7550.388823290527
63.510595568883716

'''