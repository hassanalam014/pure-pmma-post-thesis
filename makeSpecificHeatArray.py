# Date: April 2019
#
# Description: The purpose of this file is to make Excel PVT data files, in CSV format, based on SL-EOS.

#			   Wu et al, 2011. Journal of Physical and Chemical Reference Data (40).
#			   #Valid from the triple point (2.2Pa, 131.66K) to 550K and up to 50MPa.
#			   Not valid near the critical point (5.3368MPa, 400.378K).
#

import os,sys,math,matplotlib.pyplot as plt,numpy as npy
import csv
import calculatePureVariables
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
#from DME_EOS import *



#================================================================
# PARAMETERS
#================================================================

#Setting whether or not to plot results.
plot_result = True

#Setting the name of the output folder.
output_folder = 'Data'
#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Initializing the array of densities.
R0 = npy.linspace(0.02,0.65,50)

#Molar mass of PS:
M = 110000			#Toal Mass; Not Monomer Mass
Pstar = 421.762455
Tstar = 687.788143
Rstar = 1.11768655
[alpha,vh,epsilon]=calculatePureVariables.calculateMolecularParameters(Pstar,Tstar,Rstar,M)
T=303   #Kelvin
R=0.9   #g/cm^3
P=400   #MPa
P=calculatePureVariables.pressure(T,R,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
T=calculatePureVariables.temperature(P,R,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
R=calculatePureVariables.density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

# # Array of temperatures to create isotherms.
# Tlist = [423,433,443,453,463,473,483,503,523,543]

# #================================================================
# # CALCULATING THEORETICAL ISOTHERMS AND WRITING OUTPUT TO FILE
# #================================================================

# length = len(R0)+6

# for i in range(0,len(Tlist)):
	
# 	exec "T%s_P = calculatePressureDME(float(%s),R0)" % (Tlist[i],Tlist[i])
# 	exec "P0 = T%s_P" % (Tlist[i])
	
# 	output_array = npy.zeros((length,4))
	
# 	output_array[4][0] = len(R0)
	
# 	for j in range(0,len(R0)):
# 		output_array[j+6][0] = P0[j]
# 		output_array[j+6][1] = float(Tlist[i])
# 		output_array[j+6][2] = R0[j]
# 		output_array[j+6][3] = M
	
# 	exec "name = './Data/%sK_DME_PVT.csv'" % (Tlist[i])
# 	with open(name, 'wb') as f:
# 		writer = csv.writer(f)
# 		writer.writerows(output_array)

# #================================================================
# # PLOTTING RESULTS
# #================================================================

# if plot_result:
# 	plt.figure(num=None, figsize=(12, 10), dpi=80, facecolor='w', edgecolor='k')
# 	ax = plt.axes()
# 	plt.plot(T423_P,R0,'or')
# 	plt.plot(T463_P,R0,'sg')
# 	plt.plot(T483_P,R0,'^b')

# 	plt.xlabel('Pressure $P$ (MPa)',fontsize=14)
# 	plt.ylabel(r'Density $\rho$ (g/cm$^3$)',fontsize=14)
# 	plt.axis([0,85,0.0,0.6])

# 	plt.show()
