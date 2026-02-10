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

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.0000000000000001,0.9999999999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

'''
def OwnCriteria1_Deriv_Equation(epsilon_2,Rratio,Vratio,Tg_atm,dP_dT_atm,M,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	T=Tg_atm
	dP_dT=dP_dT_atm
	P=0.101325
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
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

	# Tstarstar=epsilon_2/kB
	# Tratio=Tstarstar/Tstar
	division_gap=100
	for i in range(10,20000,division_gap):
		# print i
		try:
   			# print epsilon_2
			epsilon_2 = bisect(OwnCriteria1_Deriv_Equation,i,i+division_gap,args=(Rratio,Vratio,Tg_atm,dP_dT_atm,M,Pstar,Tstar,Rstar))
			print 'Your supposed range is sussessful'
			print epsilon_2
			# print 'Rratio is'
			# print Rratio
			pass
		
		except:
   			print 'Your supposed range does not work'
   			pass

	return epsilon_2

'''

def specificHeat(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters and "three" flexibility parameters.
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	# Pratio=Pstarstar/Pstar
	# Tratio=Tstarstar/Tstar
	# Rratio=Rstarstar/Rstar
		
	# C2=deltaCp
	# C1=C_kier
	C2=((Pstar/(Rstar*Tstar))*(Vratio)*((epsilon_2/(kB*T))**2)*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))*(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))))
	C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))

	# C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
	# Not Used: C2_with v equal v0=(Pstar/(Rstar*Tstar))*((((((Tratio**2)*Rratio)/(Ttilde**2))*(exp(-(Tratio/Ttilde))))/((1+(Rratio*(exp(-(Tratio/Ttilde)))))**2)))
	# C2=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	#Not used: C2=(Pstar/(Rstar*Tstar))*(((((   (Tratio**2)*r*Rratio/N   )/(Ttilde**2))*(exp(-(((r*Tratio)/N)/Ttilde))))/((1+(Rratio*(exp(-(((r*Tratio)/N)/Ttilde)))))**2)))

	# A=3.687e-11#0.4747231	#3.363225e-10	#3.363225e-10
	# B=0.00358531	#0.003023843	#0.004156399	#0.003992431# 0.00394548
	C3=A+B*T
	#Crotation=0.0								#Pstar/(Rstar*Tstar)
	#Tstarstarstar=700.0
	#Cvibration=2.0*(((Tstarstarstar/T)**2)*(exp(-(Tstarstarstar/T)))/((1-(exp(-(Tstarstarstar/T))))**2))#*(Pstar/(Rstar*Tstar))
	C=C1+C2+C3#+Crotation+Cvibration
	#CexcludingBending=C1+Crotation+Cvibration
	return C,C1,C2,C3

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

P0 = P_atm
# T0=Tg_atm
M0=M_infinity
# R0=density(P0,T0,M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M0)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

# Ptilde=P0/Pstar
# Ttilde=T0/Tstar
# Rtilde=R0/Rstar
# vtilde=1/Rtilde	
dPtilde_dT=dP_dT_atm/Pstar
dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar

##########################################################################################################

##########################################################################################################
Vratio,Rratio,epsilon_2,x=	1.0,	1.8324276399155552,	5387.40142223434,	0.2884848484848485
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
# x=0.2742857142857143

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Rstarstar=Rratio*Rstar
Pratio=Tratio/Vratio
Pstarstar=Pratio*Pstar

#Initializing the array of densities.
T0 = npy.linspace(100,600,100)
C0=npy.zeros(len(T0))		#Total Cp
C1=npy.zeros(len(T0))		#Cp Kier
C2=npy.zeros(len(T0))		#Cp at Glass Transition
C3=npy.zeros(len(T0))		#Cp = A + B*T
C4=npy.zeros(len(T0))		#Cp_kier + A + B*T  i.e. Cp excluding Glass Transition
R0=npy.zeros(len(T0))

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
output_folder = 'plot_specificHeat'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Defining linetype
# M179_line = '-'

for i in range(0,len(T0)):
	R0[i]=density(P0,T0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

for i in range(0,len(T0)):
	C0[i],C1[i],C2[i],C3[i] = specificHeat(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,epsilon_2=epsilon_2,Rratio=Rratio,Vratio=Vratio)


for i in range(0,len(T0)):
	C4[i]=C1[i]+C3[i]

# Plotting Theoretical Curve For Difference Values of Pressure:

'''
#==============================================================================================================
#Calculating Isobars.
#==============================================================================================================

press = ['0MPa']#,'15MPa','30MPa','45MPa','60MPa','125MPa','150MPa','175MPa','200MPa']
isobar= [0.101325]#,15,30,45,60,125,150,175,200]
#print press[3]
for i in range(0,len(press)):
	exec "P0_%s = npy.full((1,len(T0)),%s)" %(press[i],isobar[i])
	exec "C0_%s = npy.zeros(len(T0))" %(press[i])
	exec "C1_%s = npy.zeros(len(T0))" %(press[i])
	exec "C2_%s = npy.zeros(len(T0))" %(press[i])
	exec "C3_%s = npy.zeros(len(T0))" %(press[i])

	exec "P0_%s =P0_%s[0] " %(press[i],press[i])
	exec "R0_%s = calculateDensity(P0_%s[0],T0,M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (press[i],press[i])
	exec "R0_%s = R0_%s[1]" %(press[i],press[i])
	for j in range(0,len(T0)):
		exec "C0_%s[j],C1_%s[j],C2_%s[j],C3_%s[j] = specificHeat(P0_%s[j],T0[j],R0_%s[j],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)"  % (press[i],press[i],press[i],press[i],press[i],press[i])
	#exec "P%s_T = result[0]" % (press[i])
	#exec "vector_%s = findVectors(P%s_T,P0,T0_%s,R0_%s)" % (press[i],press[i],press[i],press[i])



#==============================================================================================================
#Calculating isomasss.
#==============================================================================================================

mass = ['179kgpermol','208kgpermol','248kgpermol','36kgpermol','268kgpermol']
isomass= [179000,208000,248000,36000,268000]
P0=0.0034
#print mass[3]
for i in range(0,len(mass)):
	exec "M0_%s = npy.full((1,len(T0)),%s)" %(mass[i],isomass[i])
	exec "C0_%s = npy.zeros(len(T0))" %(mass[i])
	exec "C1_%s = npy.zeros(len(T0))" %(mass[i])
	exec "C2_%s = npy.zeros(len(T0))" %(mass[i])
	exec "C3_%s = npy.zeros(len(T0))" %(mass[i])

	exec "M0_%s =M0_%s[0] " %(mass[i],mass[i])
	exec "R0_%s = calculateDensity(P0,T0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (mass[i],mass[i])
	exec "R0_%s = R0_%s[1]" %(mass[i],mass[i])
	for j in range(0,len(T0)):
		exec "C0_%s[j],C1_%s[j],C2_%s[j],C3_%s[j] = specificHeat(P0,T0[j],R0_%s[j],M0_%s[j],N=N,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)"  % (mass[i],mass[i],mass[i],mass[i],mass[i],mass[i])
	#exec "P%s_T = result[0]" % (mass[i])
	#exec "vector_%s = findVectors(P%s_T,P0,T0_%s,R0_%s)" % (mass[i],mass[i],mass[i],mass[i])
'''


# C0tilde=C0*Rstar*Tstar/Pstar
# T0tilde=T0/Tstar


# molarmass = ['179K','240K','36K']							#K=kilo; not Kelvin

# for i in range(0,len(molarmass)):
# 	exec "R0=calculateDensity(P0_%s[0],T0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (molarmass[i],molarmass[i])
# 	exec "R0=R0[1]"	#Because caclulateDensity returns nested list whose 2nd row is list of R0 values
# 	exec "result = calculateSpecificHeat(P0_%s[0],T0,R0,M0_%s[0],N=N,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)" % (molarmass[i],molarmass[i])
# 	exec "M%s_C = result[0]" % (molarmass[i])
# 	#exec "vector_%s = findVectors(T%s_P,R0,P0_%s,R0_%s)" % (temp[i],temp[i],temp[i],temp[i])

#General line properties.
linewidth = 1
markersize = 1

arrow_ls = 'dashdot'
show_arrows = True
#print M179K_C
#==================================================================================
#P versus R plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

plt.axvline(x=Tg_atm,lw=0.5,color='k', linestyle='-.')

plt.plot(T0,C0,'k',color='r',lw=linewidth,ls='-',label='C_Total')
plt.plot(T0,C1,'k',color='b',lw=linewidth,ls='--',label='C_Kier theory')
plt.plot(T0,C2,'k',color='g',lw=linewidth,ls='-.',label='C_Glass theory')
plt.plot(T0,C3,'k',color='c',lw=linewidth,ls=':',label='C_A+BT theory')
plt.plot(T0,C4,'k',color='m',lw=linewidth,ls=':',label='C_excluding_Tg theory')

# plt.plot(T0,C0_0MPa,'k',color='r',lw=linewidth,ls='-',label='C_0MPa')
# plt.plot(T0,C0_15MPa,'k',color='b',lw=linewidth,ls='-',label='C_15MPa')
# plt.plot(T0,C0_30MPa,'k',color='g',lw=linewidth,ls='-',label='C_30MPa')
# plt.plot(T0,C0_45MPa,'k',color='c',lw=linewidth,ls='-',label='C_45MPa')
# plt.plot(T0,C0_60MPa,'k',color='m',lw=linewidth,ls='-',label='C_60MPa')
'''
plt.plot(T0,C0_179kgpermol,'k',color='r',lw=linewidth,ls='-',label='C_179kgpermol')
plt.plot(T0,C0_208kgpermol,'k',color='b',lw=linewidth,ls='-',label='C_208kgpermol')
plt.plot(T0,C0_248kgpermol,'k',color='g',lw=linewidth,ls='-',label='C_248kgpermol')
plt.plot(T0,C0_36kgpermol,'k',color='c',lw=linewidth,ls='-',label='C_36kgpermol')
plt.plot(T0,C0_268kgpermol,'k',color='m',lw=linewidth,ls='-',label='C_268kgpermol')
'''
#plt.plot(T524K_P,R0,'y')
#plt.plot(T0_179K,C0_179K,'ok',ms=markersize)#,label='179kilo experiment')
#plt.plot(T0_240K,C0_240K,'^k',ms=markersize)#,label='240kilo experiment')
plt.plot(T0_complete_Tg,C0_complete_Tg,'sk',ms=markersize)#,label='36kilo experiment')
#plt.plot(P0_524K,R0_524K,'*y')
plt.xlabel('Temperature T (K)',fontsize=axis_size)
plt.ylabel(r'Specific Heat $c_P$ ($J/g.K$)',fontsize=axis_size)
# plt.axis([250,450,1.00,2.25])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)
# figPUREPS.savefig('./'+output_folder+r'\pure_PS_specificHeatFitted_36kilo'+img_extension,dpi=img_dpi)

plt.show()
