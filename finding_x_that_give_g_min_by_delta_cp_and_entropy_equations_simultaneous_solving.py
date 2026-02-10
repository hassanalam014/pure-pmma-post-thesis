# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
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

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

### Programs My Theory

def glassTransitionCriteria(T,P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	# Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	# Pstarstar=Pratio*Pstar
	# Rstarstar=Rratio*Rstar

	epsilon_2=Tstarstar*kB

	# MY Theory:
	# Own_Criteria_1_incorrect=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
	# Own_Criteria_1_incorrect=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x)-((x/Vratio)*(ln(1+Rratio))))

	res=Own_Criteria_1

	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg
'''
def ResidualArray(params,P,Tg):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	M = params['M'].value
	epsilon_2 = params['epsilon_2'].value
	Rratio = params['Rratio'].value
	Vratio = params['Vratio'].value
	x = params['x'].value
	Tg_atm = params['Tg_atm'].value

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio,'x':x,'Tg_atm':Tg_atm}
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],**kwargs)
		residual[j] = abs((Tg[j]-Tg_calculated))
	
	return residual
'''

Program_Running_For=['PVME Casalini 60kilo']

Pick_List_Element = Program_Running_For[0]
Divide_List_Picked_Element = Pick_List_Element.split()

print(Divide_List_Picked_Element)

Polymer_Type=Divide_List_Picked_Element[0]
Reference=Divide_List_Picked_Element[1]
Polymer_Weight=Divide_List_Picked_Element[2]
# class Polymer_Type

kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

Abelow,Bbelow,Aabove,Babove,A,B,deltaCp,T0_excluding_Tg,M0_excluding_Tg,C0_excluding_Tg,P0_excluding_Tg,I0_excluding_Tg,Tg0_excluding_Tg,T0_above_Tg,M0_above_Tg,C0_above_Tg,P0_above_Tg,I0_above_Tg,Tg0_above_Tg,T0_at_Tg,M0_at_Tg,C0_at_Tg,P0_at_Tg,I0_at_Tg,Tg0_at_Tg,T0_below_Tg,M0_below_Tg,C0_below_Tg,P0_below_Tg,I0_below_Tg,Tg0_below_Tg,T0_complete_Tg,M0_complete_Tg,C0_complete_Tg,P0_complete_Tg,I0_complete_Tg,Tg0_complete_Tg=loadSpecificHeatExperimentalData(**kwargs)
Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

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

# print deltaCp
# deltaCp=0.4
#####################################################################################################
'''
#####################################################################################################
#Simultaneous Equation Solver, Finding g_min by changing 'x':

x= npy.linspace(0.20,0.34,50)
epsilon_2_array=npy.zeros(len(x))
Rratio_array=npy.zeros(len(x))
Vratio=1.0
mpmath.mp.dps = 15

for i in range(0,len(x)):
	
	print 'Program is iterating for the cycle number = ',i+1,' with x= ', x[i]

	Rratio = Symbol('Rratio')
	epsilon_2 = Symbol('epsilon_2')
	# Incorrect Equation below:
	# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))
	#All correct Equations below:
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x[i])-((x[i]/Vratio)*(ln(1+Rratio))))
	Cp_GlassRHS=((Pstar/(Rstar*Tstar))*(Vratio)*((epsilon_2/(kB*T))**2)*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))*(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))))-deltaCp
	# F=((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))
	# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((epsilon_2/(kB*T))*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))
	# Cp_GlassRHS=((Pstar/(Rstar*Tstar))*((((((Vratio*((epsilon_2)**2))*Rratio)/((kB*T)**2))*(exp(-(((Vratio*epsilon_2))/(kB*T)))))/((1+(Rratio*(exp(-(((Vratio*epsilon_2))/(kB*T))))))**2))))-deltaCp

	answer=nsolve((Own_Criteria_1, Cp_GlassRHS), (Rratio, epsilon_2), (1.5, 5000.0),verify=True)

	Rratio_array[i]=answer[0]
	epsilon_2_array[i]=answer[1]

	print Rratio_array[i]
	print epsilon_2_array[i]

Rratio=Rratio_array
epsilon_2=epsilon_2_array

print 'Thus, the answers are:'
Rratio_min=min(Rratio)
index_min=npy.argmin(Rratio)
epsilon_2_min=epsilon_2[index_min]
x_min=x[index_min]
x_min_upper=x[index_min+1]
x_min_lower=x[index_min-1]

print 'Vratio is:', Vratio
print 'Rratio_min is:',Rratio_min
print 'epsilon_2_min is:',epsilon_2_min
print 'x_min is:',x_min
print 'x_min_upper is:',x_min_upper
print 'x_min_lower is:',x_min_lower

# print Rratio
# print epsilon_2
# print x
######################################################################################################
'''
######################################################################################################
#Russell Trick Solve for 'x':
# x = Symbol('x')
# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
# answer= solve(Own_Criteria_1, x)
# print answer
######################################################################################################

P_line = npy.linspace(0.101325,200,15)
T_line = npy.zeros(len(P_line))
R_line=npy.zeros(len(P_line))

#Ideal Experimental Straight Line Data
for i in range(0,len(P_line)):
	T_line[i]=((P_line[i]-P)/dP_dT_atm)+T
	#R_line[i]=density(P_line[i],T_line[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

###############################################################################

###############################################################################
# Vratio,Rratio,epsilon_2,x=	1.0,	1.8324276399155552,	5387.40142223434,	0.2884848484848485
# Vratio,Rratio,epsilon_2,x= 	2.0,	4.911577067700816,	3207.618691759665,	0.28141414141414145
# Vratio,Rratio,epsilon_2,x= 	5.0,	25.989597422998806,	1762.5533180862994,	0.26653061224489794
# Vratio,Rratio,epsilon_2,x=	10.0,	143.97184062092578,	1175.499469035043,	0.25224489795918364
# Vratio,Rratio,epsilon_2,x= 	0.8,	1.3771433918540772,	6468.453995384623,	0.2889795918367347
# Vratio,Rratio,epsilon_2,x= 	0.5,	0.78128692515956,	9575.643341548659,	0.29306122448979594
# Vratio,Rratio,epsilon_2,x= 	0.3,	0.4384688373404174,	15087.456355816663,	0.2951020408163265
# Vratio,Rratio,epsilon_2,x= 	0.2,	0.2824873540080634,	22022.934776983886,	0.2951020408163265
# Vratio,Rratio,epsilon_2,x= 	0.1,	0.13642199049954273,42509.25205654484,	0.29714285714285715

Vratio=1.0
Rratio=1.311853510334076
epsilon_2=5105.083507998504
x=0.2742857142857143

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar


#Initializing the array of densities.
P = npy.linspace(0.101325,700,20)
R=npy.zeros(len(P))
Tg_calculated=npy.zeros(len(P))

for i in range(0,len(P)):

	Tg_calculated[i]=glassTemp(P[i],M=M,x=x,Rratio=Rratio,Tratio=Tratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	print P[i]

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

plt.plot(P,Tg_calculated,'k',color='g',lw=linewidth,ls='-',label='Fit')
# plt.plot(P_line,T_line,'k',color='r',lw=linewidth,ls='-',label='Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Exp. Data')

# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
# plt.axis([0,4,370,380])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()
