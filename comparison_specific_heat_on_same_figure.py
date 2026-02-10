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
	

	Rratio=Rstarstar/Rstar

	# Pratio=Pstarstar/Pstar
	# Tratio=Tstarstar/Tstar
	# C_myTheory_this_is_also_correct=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	
	epsilon_2=kB*Tstarstar
	Vratio=1.0

	C_myTheory=((Pstar/(Rstar*Tstar))*(Vratio)*((epsilon_2/(kB*T))**2)*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))*(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))))

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


# Program_Running_For=['PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PS Quach 36kilo',	'PVAc Sandberg 189kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo','LPP Passaglia 15kilo','LPP Hollander 15kilo','BPP Passaglia 15kilo','BPP Hollander 15kilo']
Program_Running_For=['PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PS Quach 36kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo']
Reference_of_Cp=['Martin','Agari','Richardson','Sasabe','Pyda','Wunderlich']

#Setting font size
axis_size = 7
title_size = 7
size = 6
label_size = 7
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
linewidth = 1.0
markersize = 3

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
fig=plt.figure(num=None, figsize=(5,7), dpi=img_dpi, facecolor='w', edgecolor='k')
# fig, ax = plt.subplots(nrows=4, ncols=2, sharex=False, sharey=False, figsize=(5, 7), dpi=img_dpi, facecolor='w', edgecolor='k')
fig.text(0.3, 0.04, 'T(K)', ha='center')														#Common x label
fig.text(0.04, 0.5,r'$\mathrm{C_P}$( $\mathrm{J/g.K}$ )', va='center', rotation='vertical')		#Common y label

ax = plt.axes()

for k in range(0,len(Program_Running_For)):

	Ref_Cp=Reference_of_Cp[k]
	Pick_List_Element = Program_Running_For[k]
	Divide_List_Picked_Element = Pick_List_Element.split()

	print(Divide_List_Picked_Element)

	Polymer_Type=Divide_List_Picked_Element[0]
	Reference=Divide_List_Picked_Element[1]
	Polymer_Weight=Divide_List_Picked_Element[2]
	# class Polymer_Type

	kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

	Abelow,Bbelow,Aabove,Babove,A,B,deltaCp,T0_excluding_Tg,M0_excluding_Tg,C0_excluding_Tg,P0_excluding_Tg,I0_excluding_Tg,Tg0_excluding_Tg,T0_above_Tg,M0_above_Tg,C0_above_Tg,P0_above_Tg,I0_above_Tg,Tg0_above_Tg,T0_at_Tg,M0_at_Tg,C0_at_Tg,P0_at_Tg,I0_at_Tg,Tg0_at_Tg,T0_below_Tg,M0_below_Tg,C0_below_Tg,P0_below_Tg,I0_below_Tg,Tg0_below_Tg,T0_complete_Tg,M0_complete_Tg,C0_complete_Tg,P0_complete_Tg,I0_complete_Tg,Tg0_complete_Tg=loadSpecificHeatExperimentalData(**kwargs)
	Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

	P0 = P_atm
	M0=M_infinity
	r = (Pstar*M0)/(kB*Tstar*Rstar)
	Vratio=1.0
	##########################################################################################################

	##########################################################################################################

	for_condoTheory=Polymer_Type+' '+Reference
	for_myTheory=Polymer_Type+' '+Reference+' '+Polymer_Weight

	condoTheory_List= 		['PMMA Grassia','PMMA Olabisi',	'PS Quach',		'PVAc Sandberg','PVAc Roland',	'PVME Casalini',	'PC Zoller',	'LPP Passaglia',	'LPP Hollander',	'BPP Passaglia',		'BPP Hollander']
	condoTheory_Rratio_List=		[2.0,			3.0,			3.0,			2.0,			3.0,			3.0,			2.0,			3.0,			3.0,				2.0,					3.0]
	condoTheory_epsilon_2_List= 	[4714.115862,	7428.006111,	7160.480018,	4463.798377,	6029.978362,	4480.294232,	6247.794198,	4168.372824,	4340.291764,		2647.071561,	4354.015059]
	index_condoTheory = condoTheory_List.index(for_condoTheory)
	Rratio_condoTheory=	condoTheory_Rratio_List[index_condoTheory]
	epsilon_2_condoTheory=condoTheory_epsilon_2_List[index_condoTheory]

	# myTheory_List_incorrect=			['PMMA Grassia 44kilo',	'PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PMMA Olabisi 140kilo','PS Quach 36kilo',	'PVAc Sandberg 00kilo','PVAc Sandberg 01kilo','PVAc Sandberg 189kilo','PVAc Roland 00kilo','PVAc Roland 01kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo','PC Zoller 01kilo',	'LPP Passaglia 15kilo','LPP Hollander 15kilo','BPP Passaglia 15kilo','BPP Hollander 15kilo']
	# myTheory_Rratio_List_incorrect=	[1.99735062,  			1.07693026,   	   		1.65676715,   	   	1.04721418,  	 	1.66780617,  		2.04023707, 		   	1.49366076, 		  1.64252435, 		  		2.13863669,   	   1.56116973, 	    	1.90502640242404, 		  	1.83250186,		   	0.83946858,	  		0.27061254,					2.46498986,			2.47102959,				2.15196345,			2.15851346]
	# myTheory_epsilon_2_List_incorrect=[7854.51586017724,    	7065.33508891425,    	8085.95247387189, 	7517.92081962379, 	8104.91224307777, 	7042.06893815706,    	6698.88658094399,     6801.7047548495,      	6890.2761804729,   6537.88343666762,  	6812.53215419597,   	5365.73458238909,  	8154.30819514158,	7284.7806806012,			5560.89299743581,	5744.54221056484,		5441.23463670387,	5615.10389303674]
	# myTheory_x_List_incorrect=		[0.271428571,  			0.271428571,		   	0.281632653, 	   	0.281632653, 		0.271428571,      	0.281632653, 		   	0.281632653	, 	      0.281632653,  		  	0.281632653,  	   0.281632653,			0.278383838383838, 		  	0.271428571, 	   	0.291836735, 	  	0.281632653,				0.278383838,		0.276363636,			0.274343434,		0.273333333]

	myTheory_List=			['PMMA Grassia 140kilo',	'PMMA Olabisi 44kilo',		'PS Quach 36kilo',		'PVAc Roland 189kilo',	'PVME Casalini 60kilo',		'PC Zoller 00kilo']
	myTheory_Rratio_List=	[1.0769079749464847,			1.6567603698535678,	 1.6674132169170892, 		1.905027683861342,		 1.8324276399155552,		 0.8392099749106786]
	myTheory_epsilon_2_List=[7093.578262388853,				8094.136627663586,	 	8012.605528610325, 		6814.75937995625, 			5387.40142223434,		 8273.34220202496]
	myTheory_x_List=		[0.2925252525252525,			0.32282828282828285,	 0.3107070707070707,	 0.3208080808080808,	 0.2884848484848485,		 0.3167676767676768]

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
	T_max=T0_above_Tg[-1]
	T_min=T0_above_Tg[0]
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

	# delete_first_n_elements=23
	# T0 = T0[delete_first_n_elements:]
	# C_total_myTheory = C_total_myTheory[delete_first_n_elements:]
	# C_total_condoTheory = C_total_condoTheory[delete_first_n_elements:]

	# delete_first_n_experiment_data=12
	# T0_complete_Tg=T0_complete_Tg[delete_first_n_experiment_data:]
	# C0_complete_Tg=C0_complete_Tg[delete_first_n_experiment_data:]
	# Plotting Theoretical Curve For Difference Values of Pressure:

	plt.subplot(3, 2, k+1)

	# plt.axvline(x=Tg_atm,lw=0.5,color='k', linestyle='--',label='Glass Temperature')

	# plt.plot(T0,C_line_below_Tg,'k',color='r',lw=linewidth,ls='-',label='C_line_below_Tg')
	# plt.plot(T0,C_line_above_Tg,'k',color='r',lw=linewidth,ls='-',label='C_line_above_Tg')
	# plt.plot(T0,C_kier,'k',color='b',lw=linewidth,ls='--',label='C_Kier theory')
	# plt.plot(T0,C_baseFit_below_Tg,'k',color='g',lw=linewidth,ls='-.',label='Linear-Fit below Tg')
	# plt.plot(T0,C_baseFit_above_Tg,'k',color='g',lw=linewidth,ls='-.',label='Linear-Fit above Tg')
	# plt.plot(T0,C_myTheory,'k',color='c',lw=linewidth,ls=':',label='C_myTheory Only')
	# plt.plot(T0,C_condoTheory,'k',color='m',lw=linewidth,ls=':',label='C_condoTheory Only')
	plt.plot(T0,C_total_myTheory,'k',color='k',lw=linewidth,ls='-',label='Present Theory')
	plt.plot(T0,C_total_condoTheory,'k',color='k',lw=linewidth,ls='--',label='Condo et al.')

	plt.plot(T0_above_Tg,C0_above_Tg,'sk',ms=markersize,label='Experiment')
	# plt.plot(T0_complete_Tg,C0_complete_Tg,'sk',ms=markersize,label='Experiment Data of {}'.format(Polymer_Type))

	# plt.axis([325,342,1.80,2.00])

	if k==0:
		# plt.xlabel('T(K)',fontsize=axis_size)
		# plt.ylabel(r'$\mathrm{C_P}$( $\mathrm{J/g.K}$ )',fontsize=axis_size)
		plt.legend(loc=4,fontsize=size,numpoints=1)

	# plt.yticks(rotation=45)
	plt.title(' {}) {} - {}'.format(chr(k+97),Polymer_Type,Ref_Cp), fontsize=title_size,x=0.0,y=0.84,horizontalalignment='left')

plt.subplots_adjust(left=0.15,right=0.95,top=0.95,bottom=0.10,wspace=0.30,hspace=0.25)
fig.savefig('./'+output_folder+r'\pure_Cp vs T_5'+img_extension,dpi=240)
plt.show()
