# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
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
import cmath

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

### Programs My Theory
'''
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

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))

	res=Own_Criteria_1

	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg
'''

def glassTransition_for_eps2(epsilon_2,P,T,R,M,x,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	Pratio=Tratio/Vratio

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
	# print Own_Criteria_1
	res=Own_Criteria_1

	return res

def eps_2(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
		
	for i in range(0,50000,500):
		
		epsilon_2=0.0
		try:
			epsilon_2 = bisect(glassTransition_for_eps2,i,i+500,args=(P,T,R,M,x,Rratio,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get epsilon_2")
			pass
		if epsilon_2!=0.0:
			print 'Hurry! epsilon_2_dependent is:', epsilon_2, 'Rratio_independent is:', Rratio, 'and x_independent=',x
			break
	
	if epsilon_2==0.0:
		print 'Program Failed to get value of epsilon_2 against Rratio_independent=', Rratio, 'and x_independent=',x
		# epsilon_2=50000

	return epsilon_2

def glassTransition_for_Rratio(Rratio,P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	Pratio=Tratio/Vratio

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))

	# print Own_Criteria_1

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
			print 'Hurry! Rratio_dependent is:', Rratio, 'against epsilon_2_independent =', epsilon_2, 'and x_independent=',x
			break
	
	if Rratio==0.0:
		print 'Program Failed to get value of Rratio agasint epsilon_2_independent=',epsilon_2, 'and x_independent=',x
		# Rratio=50000

	return Rratio

def glassTransition_for_x(x,P,T,R,M,epsilon_2,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	Pratio=Tratio/Vratio

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))

	# print Own_Criteria_1

	res=Own_Criteria_1

	return res

def x_atm(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	for i in range(0,11):
		
		lower_limit=i/10
		upper_limit=(i/10)+0.1
		x=0.0
		try:
			x = bisect(glassTransition_for_x,lower_limit,upper_limit,args=(P,T,R,M,epsilon_2,Rratio,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get Rratio")
			pass
		if x!=0.0:
			print 'Hurry! x_dependent is:', x, 'against epsilon_2_independent =', epsilon_2, 'and Rratio_independent =', Rratio
			break
	
	if x==0.0:
		print 'Program Failed to get value of x against epsilon_2_independent=',epsilon_2,'and Rratio_independent =', Rratio

	return x

'''
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
'''

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

Vratio=1.0
#####################################################################################################

######################################################################################################
number_of_point=5
epsilon_2_independent = npy.linspace(2000,10000,number_of_point)
Rratio_independent = npy.linspace(1.0,5,number_of_point)
x_independent = npy.linspace(0.2,0.6,number_of_point)

list_of_numbers = list(range(number_of_point))

for j in range(0,len(list_of_numbers)):
	exec "epsilon_2_dependent_constR_%s = npy.zeros(len(epsilon_2_independent))" %(list_of_numbers[j])
	exec "epsilon_2_dependent_constX_%s = npy.zeros(len(epsilon_2_independent))" %(list_of_numbers[j])
	exec "epsilon_2_isocurve_%s = npy.linspace(epsilon_2_independent[j],epsilon_2_independent[j],number_of_point)" %(list_of_numbers[j])
	exec "Rratio_dependent_constE_%s = npy.zeros(len(epsilon_2_independent))" %(list_of_numbers[j])
	exec "Rratio_dependent_constX_%s = npy.zeros(len(epsilon_2_independent))" %(list_of_numbers[j])
	exec "Rratio_isocurve_%s = npy.linspace(Rratio_independent[j],Rratio_independent[j],number_of_point)" %(list_of_numbers[j])
	exec "x_dependent_constR_%s = npy.zeros(len(epsilon_2_independent))" %(list_of_numbers[j])
	exec "x_dependent_constE_%s = npy.zeros(len(epsilon_2_independent))" %(list_of_numbers[j])
	exec "x_isocurve_%s = npy.linspace(x_independent[j],x_independent[j],number_of_point)" %(list_of_numbers[j])
	
	for i in range(0,len(epsilon_2_independent)):
		exec "epsilon_2_dependent_constR_%s[i]=eps_2(P,T,R,M,x=x_independent[i],Rratio=Rratio_independent[j],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" %(list_of_numbers[j])
		exec "epsilon_2_dependent_constX_%s[i]=eps_2(P,T,R,M,x=x_independent[j],Rratio=Rratio_independent[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" %(list_of_numbers[j])
		exec "Rratio_dependent_constE_%s[i]=Rrat(P,T,R,M,x=x_independent[i],epsilon_2=epsilon_2_independent[j],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" %(list_of_numbers[j])
		exec "Rratio_dependent_constX_%s[i]=Rrat(P,T,R,M,x=x_independent[j],epsilon_2=epsilon_2_independent[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" %(list_of_numbers[j])
		exec "x_dependent_constR_%s[i]=x_atm(P,T,R,M,epsilon_2=epsilon_2_independent[i],Rratio=Rratio_independent[j],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" %(list_of_numbers[j])
		exec "x_dependent_constE_%s[i]=x_atm(P,T,R,M,epsilon_2=epsilon_2_independent[j],Rratio=Rratio_independent[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" %(list_of_numbers[j])

###############################################################################

###############################################################################

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
markersize = 3

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

# plt.plot(epsilon_2_independent,Rratio_dependent,'k',color='g',lw=linewidth,ls='-',label='eps_indep.')

color=['r','b','g','m','y','c','k','r','b','g','m','y','c','k']

for j in range(0,len(list_of_numbers)):
	# exec "print epsilon_2_dependent_%s" %(list_of_numbers[j])
	# exec "print Rratio_isocurve_%s" %(list_of_numbers[j])

	# exec "plt.plot(epsilon_2_dependent_constR_%s,Rratio_isocurve_%s,'sk',color='%s',ms=markersize,label='eps_dep.')" %(list_of_numbers[j],list_of_numbers[j],color[j])
	# exec "plt.plot(epsilon_2_isocurve_%s,Rratio_dependent_constE_%s,'sk',color='%s',ms=markersize,label='eps_dep.')" %(list_of_numbers[j],list_of_numbers[j],color[j])
	exec "plt.plot(epsilon_2_dependent_constX_%s,Rratio_independent,'k',color='%s',ms=markersize,label='eps_dep.')" %(list_of_numbers[j],color[j])
	# exec "plt.plot(epsilon_2_independent,Rratio_dependent_constX_%s,'k',color='%s',ms=markersize,label='eps_dep.')" %(list_of_numbers[j],color[j])

plt.xlabel('ep_2',fontsize=axis_size)
plt.ylabel(r'Rratio',fontsize=axis_size)
plt.axis([3000,9000,0,5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()
