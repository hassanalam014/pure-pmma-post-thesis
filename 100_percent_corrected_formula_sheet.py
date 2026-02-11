R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M)/(kB*Tstar*Rstar)
#####################################################################################################
Condo Paper Equation Set
#####################################################################################################
F_condo=((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))
S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-epsilon_2/(kB*T)))/(1+(Rratio)*exp(-epsilon_2/(kB*T)))))-((((Rratio)*exp(-epsilon_2/(kB*T)))/(1+(Rratio)*exp(-epsilon_2/(kB*T))))*epsilon_2/(kB*T))))
S_condo_again=(Pstar/(Rstar*Tstar))*(-(((1-Rtilde)/Rtilde)*(ln(1-Rtilde)))-((ln(Rtilde))/r)+((ln(r))/r)-1-(((ln(2/(Rratio+2)))-1)/r)-(((r-2)/r)*((ln(1-((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))))-((((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T)))))
dS_dTg_condo=(1+(ln(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(ln(1-Rtilde)+Rtilde))-(((((Rratio)*exp(-epsilon_2/(kB*T)))/(1+((Rratio)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(((1-(((Rratio)*exp(-epsilon_2/(kB*T)))/(1+((Rratio)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)
dS_dTg_condo_again=(((r-2)/r)*((epsilon_2/(kB*T))**2)*(((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))*(1-((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T)))))))/T)*((Rtilde)**2)*(((Ttilde/Rtilde)*((1/r)+(Rtilde/(1-Rtilde))))-2))+((((ln(1-Rtilde))/(Rtilde))+1-(1/r))*((dPtilde_dT)+((1/Tstar)*((ln(1-Rtilde))+((1-(1/r))*Rtilde)))))
C_condoTheory_total_including_Kier_Part=(Pstar/(Rstar*Tstar))*(((((1-(1/r))+((ln(1-Rtilde))/(Rtilde)))**2)/((1/r)+(Rtilde/(1-Rtilde))-(2*Rtilde/Ttilde)))+(((r-2)/r)*((epsilon_2/(kB*T))**2)*(F_condo)*(1-F_condo)))
C_condoTheory_excluding_kier_part=(Pstar/(Rstar*Tstar))*((((r-2)/r)*((epsilon_2/(kB*T))**2)*(F_condo)*(1-F_condo)))
C_kier_part=(Pstar/(Rstar*Tstar))*((((1-(1/r))+((ln(1-Rtilde))/(Rtilde)))**2)/((1/r)+(Rtilde/(1-Rtilde))-(2*Rtilde/Ttilde)))

######################################################################################################
MY THEORY
######################################################################################################
Own_Criteria_1_wrong=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
Own_Criteria_1_wrong=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
dS_dT_Own_Criteria_1=(Vratio*Ttilde*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((1/Rtilde)*(dPtilde_dTtilde)))/(Rtilde*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))+(((Vratio*epsilon_2/(kB*T))**2)*((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))*(1-((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))))
#100% Correct Following:
Best: Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x[i])-((x[i]/Vratio)*(ln(1+Rratio))))
Best: Cp_GlassRHS=((Pstar/(Rstar*Tstar))*(Vratio)*((epsilon_2/(kB*T))**2)*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))*(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))))-deltaCp
Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((epsilon_2/(kB*T))*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))
Cp_GlassRHS=((Pstar/(Rstar*Tstar))*((((((Vratio*((epsilon_2)**2))*Rratio)/((kB*T)**2))*(exp(-(((Vratio*epsilon_2))/(kB*T)))))/((1+(Rratio*(exp(-(((Vratio*epsilon_2))/(kB*T))))))**2))))-deltaCp

# Following Specific Heat. Not 100% Sure But Very Very hopeful that following are also right:
Delta_Cp=((Pstar/(Rstar*Tstar))*((((((Vratio*((epsilon_2)**2))*Rratio)/((kB*T)**2))*(exp(-(((Vratio*epsilon_2))/(kB*T)))))/((1+(Rratio*(exp(-(((Vratio*epsilon_2))/(kB*T))))))**2))))
Cp_Kier=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
Cp_myTheory_excluding_kier_part=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
C_myTheory_total_including_kier_part=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((1/r)+(Rtilde/(1-Rtilde))))-2)))+(((Rratio*Tratio**3)/(Pratio*Ttilde**2))*((exp(-(Tratio**2)/(Pratio*Ttilde)))/((1+(Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))**2))))
