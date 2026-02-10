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
from loadSpecificHeatExperimentalData import *
from sympy import Symbol, nsolve
import sympy
import mpmath
from Parameters_of_Different_Polymers import *
from matplotlib.ticker import FormatStrFormatter #To fix axis decimal places


def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

#Condo Theroy Programs
'''
def TgCriteriaForRratio_condoTheory(Rratio,P,T,M,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
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

def Rrat_condoTheory(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Rratio = bisect(TgCriteriaForRratio_condoTheory,0.4,100,args=(P,T,M,Tratio,Vratio,Pstar,Tstar,Rstar))
	return Rratio

def TgCriteriaForEpsilon_condoTheory(epsilon_2,P,T,M,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
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

def eps_condoTheory(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	# epsilon_2 = bisect(TgCriteriaForEpsilon,1500,20000,args=(P,T,M,Rratio,Vratio,Pstar,Tstar,Rstar))
		
	for i in range(0,50000,500):
		
		epsilon_2=0.0
		try:
			epsilon_2 = bisect(TgCriteriaForEpsilon_condoTheory,i,i+500,args=(P,T,M,Rratio,Vratio,Pstar,Tstar,Rstar))
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
'''
def glassTransitionCriteria_condoTheory(T,P,M,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
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

def glassTemp_condoTheory(P,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria_condoTheory,100,10000,args=(P,M,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg

### Programs My Theory

def glassTransitionCriteria_myTheory(T,P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
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

def glassTemp_myTheory(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria_myTheory,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg

# Program_Running_For=['PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PS Quach 36kilo','PVAc Sandberg 189kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo','LPP Passaglia 15kilo','LPP Hollander 15kilo','BPP Passaglia 15kilo','BPP Hollander 15kilo']
Program_Running_For=['PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PS Quach 36kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo']

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
# fig, ax = plt.subplots(nrows=4, ncols=2, sharex=True, sharey=True, figsize=(5, 7), dpi=img_dpi, facecolor='w', edgecolor='k')
fig.text(0.3, 0.04, 'P(MPa)', ha='center')										#Common x label
fig.text(0.04, 0.5,r'$\mathrm{T_g} $(K)', va='center', rotation='vertical')		#Common y label

ax = plt.axes()

for k in range(0,len(Program_Running_For)):

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

	M=M_infinity
	Vratio=1.0
	dTg_dP_atm_extreme_ends = (T_upper-Tg_exp[0])/(P_upper-Pg_exp[0])         #Line passing through extreme ends of only linear range of exp data
	dP_dT_extreme_ends=1/dTg_dP_atm_extreme_ends

	P_line = npy.linspace(0.101325,P_upper,10)
	Tg_line = npy.zeros(len(P_line))

	#Ideal Experimental Straight Line Data for Fitting:
	for i in range(0,len(P_line)):
		Tg_line[i]=((P_line[i]-P_atm)/dP_dT_extreme_ends)+Tg_atm		#Extreme End Line Passing Line

	############################################################################################

	##########################################################################################

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

	Tstarstar_myTheory=epsilon_2_myTheory/kB
	Tratio_myTheory=Tstarstar_myTheory/Tstar

	P_max=Pg_exp[-1]

	#Initializing the array of densities.
	P=npy.linspace(0.101325,P_max,20)
	R=npy.zeros(len(P))
	Tg_calculated_condoTheory=npy.zeros(len(P))
	Tg_calculated_myTheory=npy.zeros(len(P))

	for i in range(0,len(P)):
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio_condoTheory,'Rratio':Rratio_condoTheory,'Vratio':Vratio}
		Tg_calculated_condoTheory[i] = glassTemp_condoTheory(P[i],M,**kwargs)
		Tg_calculated_myTheory[i]=glassTemp_myTheory(P[i],M=M,x=x_myTheory,Rratio=Rratio_myTheory,Tratio=Tratio_myTheory,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

		print P[i]

	plt.subplot(3, 2, k+1)


	plt.plot(P,Tg_calculated_myTheory,'k',color='k',lw=linewidth,ls='-',label='Present Theory')
	plt.plot(P,Tg_calculated_condoTheory,'k',color='k',lw=linewidth,ls='--',label='Condo et al.')
	# plt.plot(P_line,Tg_line,'k',color='r',lw=linewidth,ls='-',label='Ideal Straight Line')
	plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Experiment')
	# plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Experiment {}'.format(Polymer_Type))

	# plt.axhline(y=243.5,lw=0.5,color='k', linestyle='-.')
	# plt.axvline(x=0.101325,lw=0.5,color='k', linestyle='-.')
	# ax.text(0.3,0.9,'{}'.format(Polymer_Type),	horizontalalignment='left',transform=ax.transAxes)

	#plt.axis([300,500,0,1.5])
	if k==0:
		# plt.xlabel('P(MPa)',fontsize=axis_size)
		# plt.ylabel(r'$\mathrm{T_g} $(K)',fontsize=axis_size)
		plt.legend(loc=4,fontsize=size,numpoints=1)

	plt.title(' {}) {} - {}'.format(chr(k+97),Polymer_Type,Reference), fontsize=title_size,x=0.0,y=0.84,horizontalalignment='left')
	# ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))	#Number of decimal places is 2

plt.subplots_adjust(left=0.15,right=0.95,top=0.95,bottom=0.10,wspace=0.30,hspace=0.25)
fig.savefig('./'+output_folder+r'\pure_Tg vs P_5'+img_extension,dpi=240)
plt.show()
