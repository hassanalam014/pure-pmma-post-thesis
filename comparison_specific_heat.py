# Date: April 2017
#
# Description: The purpose of this file is to plot Polystyrene (PS) density information based on experiment and theory for comparison.
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
# from all_p_params import *
# from loadSpecificHeatExperimentalData import *
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
from loadSpecificHeatExperimentalData import *
from Parameters_of_Different_Polymers import *
from matplotlib.ticker import FormatStrFormatter #To fix axis decimal places

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.0000000000000001,0.9999999999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def specificHeat_myTheory(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters and "three" flexibility parameters.
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	Pratio=Pstarstar/Pstar
	Tratio=Tstarstar/Tstar
	Rratio=Rstarstar/Rstar

	C_myTheory=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	
	return C_myTheory

def specificHeat_kier(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	C_kier=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
	
	return C_kier

def specificHeat_line_below_Tg(P,T,R,M,**kwargs):     #Cp/mass  and A, B
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	C_line_below_Tg=A+B*T
	
	return C_line_below_Tg

def specificHeat_line_above_Tg(P,T,R,M,**kwargs):     #Cp/mass  and A, B
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	C_line_above_Tg=A+B*T
	
	return C_line_above_Tg

def specificHeat_condoTheory(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters and "three" flexibility parameters.
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	Pratio=Pstarstar/Pstar
	Tratio=Tstarstar/Tstar
	Rratio=Rstarstar/Rstar

	epsilon_2=Tstarstar*kB
	F_condo=((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))

	C_condoTheory=(Pstar/(Rstar*Tstar))*((((r-2)/r)*((epsilon_2/(kB*T))**2)*(F_condo)*(1-F_condo)))

	return C_condoTheory

P0 = P_atm
M0=M_infinity
r = (Pstar*M0)/(kB*Tstar*Rstar)
Vratio=1.0
##########################################################################################################

##########################################################################################################

for_condoTheory=Polymer_Type+' '+Reference
for_myTheory=Polymer_Type+' '+Reference+' '+Polymer_Weight

condoTheory_List= 		['PMMA Grassia','PMMA Condo',	'PS Condo',		'PVAc Sandberg','PVAc Roland',	'PVME Roland',	'PC Zoller']
condoTheory_Rratio_List=		[2.0,			3.0,			3.0,			2.0,			3.0,			3.0,			2.0]
condoTheory_epsilon_2_List= 	[4714.115862,	7428.006111,	7160.480018,	4463.798377,	6029.978362,	4480.294232,	6247.794198]
index_condoTheory = condoTheory_List.index(for_condoTheory)
Rratio_condoTheory=	condoTheory_Rratio_List[index_condoTheory]
epsilon_2_condoTheory=condoTheory_epsilon_2_List[index_condoTheory]

myTheory_List=			['PMMA Grassia 44kilo',	'PMMA Grassia 140kilo','PMMA Condo 44kilo','PMMA Condo 140kilo','PS Condo 36kilo',	'PVAc Sandberg 00kilo','PVAc Sandberg 01kilo','PVAc Sandberg 189kilo','PVAc Roland 00kilo','PVAc Roland 01kilo','PVAc Roland 189kilo','PVME Roland 60kilo','PC Zoller 00kilo','PC Zoller 01kilo']
myTheory_Rratio_List=	[1.99735062,  			1.07693026,   	   		1.65676715,   	   	1.04721418,  	 	1.66780617,  		2.04023707, 		   	1.49366076, 		  1.64252435, 		  		2.13863669,   	   1.56116973, 	    	1.78836852, 		  	1.83250186,		   	0.83946858,	  		0.27061254]
myTheory_epsilon_2_List=[7854.51586017724,    	7065.33508891425,    	8085.95247387189, 	7517.92081962379, 	8104.91224307777, 	7042.06893815706,    	6698.88658094399,     6801.7047548495,      	6890.2761804729,   6537.88343666762,  	6688.19579102878,   	5365.73458238909,  	8154.30819514158,	7284.7806806012]
myTheory_x_List=		[0.271428571,  			0.271428571,		   	0.281632653, 	   	0.281632653, 		0.271428571,      	0.281632653, 		   	0.281632653	, 	      0.281632653,  		  	0.281632653,  	   0.281632653,			0.281632653, 		  	0.271428571, 	   	0.291836735, 	  	0.281632653]

index_myTheory = myTheory_List.index(for_myTheory)
Rratio_myTheory=	myTheory_Rratio_List[index_myTheory]
epsilon_2_myTheory=myTheory_epsilon_2_List[index_myTheory]
x_myTheory=	myTheory_x_List[index_myTheory]

Tstarstar_condoTheory=epsilon_2_condoTheory/kB
Tratio_condoTheory=Tstarstar_condoTheory/Tstar
Rstarstar_condoTheory=Rratio_condoTheory*Rstar
Pratio_condoTheory=Tratio_condoTheory/Vratio
Pstarstar_condoTheory=Pratio_condoTheory*Pstar

Tstarstar_myTheory=epsilon_2_myTheory/kB
Tratio_myTheory=Tstarstar_myTheory/Tstar
Rstarstar_myTheory=Rratio_myTheory*Rstar
Pratio_myTheory=Tratio_myTheory/Vratio
Pstarstar_myTheory=Pratio_myTheory*Pstar

#Initializing the array of densities.
T_max=T0_complete_Tg[-1]
T_min=T0_complete_Tg[0]
T0 = npy.linspace(T_min,T_max,100)

C_line_below_Tg=npy.zeros(len(T0))				#Cp = A + B*T
C_line_above_Tg=npy.zeros(len(T0))				#Cp = A + B*T
C_kier=npy.zeros(len(T0))				#Cp Kier
C_baseFit_below_Tg=npy.zeros(len(T0))			#Cp_kier + C_line  i.e. Cp Base Curve Fit
C_baseFit_above_Tg=npy.zeros(len(T0))			#Cp_kier + C_line  i.e. Cp Base Curve Fit
C_myTheory=npy.zeros(len(T0))			#Cp My Theory Only
C_condoTheory=npy.zeros(len(T0))		#Cp Condo Theory Only
C_total_myTheory=npy.zeros(len(T0))		#Total Cp from My Theory
C_total_condoTheory=npy.zeros(len(T0))	#Total Cp from Condo Theory
R0=npy.zeros(len(T0))

#Setting font size
axis_size = 20			#Size of axis labels
title_size = 20			#Size of title
size = 14				#Size of legender
label_size = 20			#Size of axis numbers
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'plot_specificHeat'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

for i in range(0,len(T0)):
	R0[i]=density(P0,T0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

for i in range(0,len(T0)):
	C_line_below_Tg[i] = specificHeat_line_below_Tg(P0,T0[i],R0[i],M0,A=A,B=B)
	C_line_above_Tg[i] = specificHeat_line_above_Tg(P0,T0[i],R0[i],M0,A=Aabove,B=Babove)
	C_kier[i] = specificHeat_kier(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	C_baseFit_below_Tg[i] = C_line_below_Tg[i]+C_kier[i]
	C_baseFit_above_Tg[i] = C_line_above_Tg[i]+C_kier[i]
	C_myTheory[i] = specificHeat_myTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar_myTheory,Tstarstar=Tstarstar_myTheory,Rstarstar=Rstarstar_myTheory)
	C_condoTheory[i] = specificHeat_condoTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar_condoTheory,Tstarstar=Tstarstar_condoTheory,Rstarstar=Rstarstar_condoTheory)
	C_total_myTheory[i] = C_myTheory[i] + C_kier[i] + C_line_below_Tg[i]
	C_total_condoTheory[i] = C_condoTheory[i] + C_kier[i] + C_line_below_Tg[i]

delete_first_n_elements=23
T0 = T0[delete_first_n_elements:]
C_total_myTheory = C_total_myTheory[delete_first_n_elements:]
C_total_condoTheory = C_total_condoTheory[delete_first_n_elements:]

delete_first_n_experiment_data=12
T0_complete_Tg=T0_complete_Tg[delete_first_n_experiment_data:]
C0_complete_Tg=C0_complete_Tg[delete_first_n_experiment_data:]
# Plotting Theoretical Curve For Difference Values of Pressure:

#General line properties.
linewidth = 1
markersize = 4

arrow_ls = 'dashdot'
show_arrows = True
#==================================================================================
figPUREPS=plt.figure(num=None, figsize=(8,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

# plt.axvline(x=Tg_atm,lw=0.5,color='k', linestyle='--',label='Glass Temperature')

# plt.plot(T0,C_line_below_Tg,'k',color='r',lw=linewidth,ls='-',label='C_line_below_Tg')
# plt.plot(T0,C_line_above_Tg,'k',color='r',lw=linewidth,ls='-',label='C_line_above_Tg')
# plt.plot(T0,C_kier,'k',color='b',lw=linewidth,ls='--',label='C_Kier theory')
# plt.plot(T0,C_baseFit_below_Tg,'k',color='g',lw=linewidth,ls='-.',label='Linear-Fit below Tg')
# plt.plot(T0,C_baseFit_above_Tg,'k',color='g',lw=linewidth,ls='-.',label='Linear-Fit above Tg')
# plt.plot(T0,C_myTheory,'k',color='c',lw=linewidth,ls=':',label='C_myTheory Only')
# plt.plot(T0,C_condoTheory,'k',color='m',lw=linewidth,ls=':',label='C_condoTheory Only')
plt.plot(T0,C_total_myTheory,'k',color='r',lw=linewidth,ls='-',label='Present Theory')
plt.plot(T0,C_total_condoTheory,'k',color='b',lw=linewidth,ls='-.',label='Condo et. al.')

plt.plot(T0_complete_Tg,C0_complete_Tg,'sk',ms=markersize,label='Experiment')
# plt.plot(T0_complete_Tg,C0_complete_Tg,'sk',ms=markersize,label='Experiment Data of {}'.format(Polymer_Type))
plt.xlabel('T(K)',fontsize=axis_size)
plt.ylabel(r'$\mathrm{C_P}$( $\mathrm{J/g.K}$ )',fontsize=axis_size)
# plt.axis([325,342,1.80,2.00])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(left=0.20,right=0.90,top=0.95,bottom=0.15)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))	#Number of decimal places is 2
# figPUREPS.savefig('./'+output_folder+r'\pure_PS_specificHeatFitted_36kilo'+img_extension,dpi=img_dpi)

plt.show()
