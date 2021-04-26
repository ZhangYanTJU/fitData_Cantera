import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import math


mech='mech.cti'
gas = ct.Solution(mech)
Tlow = 300
Tcommon = 1400
Thigh = 5000
deltaT = 100
T_low = np.arange(Tlow,Tcommon+deltaT,deltaT)
T_high = np.arange(Tcommon,Thigh+deltaT,deltaT)


species_name = "N2"
sp_index = gas.species_index(species_name)
species = gas.species(sp_index)
molecular_weight = gas.molecular_weights[sp_index]
nasa_coeffs = species.thermo.coeffs # orininal coeffs, size = 15, [0]=Tcommon, [1-7] is high coeffs, [8-14] is low coeffs

Tlist = np.concatenate([T_low, T_high])
cp = [species.thermo.cp(T)/molecular_weight for T in Tlist] # [J/(kg·K)]





a1_a5_org_low = nasa_coeffs[8:13]
a6_org_low = nasa_coeffs[13]
a1_a5_org_high = nasa_coeffs[1:6]
a6_org_high = nasa_coeffs[6]
                
cp_bad_low = [ np.sum(np.array([a1_a5_org_low[i]*T**(i) for i in range(5)])) for T in T_low] # cp/R
cp_bad_high = [ np.sum(np.array([a1_a5_org_high[i]*T**(i) for i in range(5)])) for T in T_high] # cp/R
R = ct.gas_constant/molecular_weight  # [J/(kmol·K)] / [kg/kmol] = [J/(kg·K)]
cp_bad = np.concatenate([cp_bad_low, cp_bad_high])*R


#enthalpy = [species.thermo.h(T)/molecular_weight for T in Tlist] # J/kg
#enthalpy_fit_low = [ a6_low + np.sum(np.array([a1_a5_low[i]*T**(i+1)/(i+1) for i in range(5)])) for T in T_low] # h/R
#enthalpy_fit_high = [ a6_high + np.sum(np.array([a1_a5_high[i]*T**(i+1)/(i+1) for i in range(5)])) for T in T_high] # h/R
#enthalpy_fit = np.concatenate([enthalpy_fit_low, enthalpy_fit_high])*R

fig = plt.figure()

ax1 = fig.add_subplot()
ax1.scatter(Tlist, cp, c='blue', label="cp_org, Tcommon=%d"%(species.thermo.coeffs[0]))
ax1.plot(Tlist, cp_bad, c='purple', label="cp_bad, Tcommon=%d"%(Tcommon))
ax1.set_xlabel("T [K]")
ax1.set_ylabel("cp [J/kg]")

#ax2 = plt.twinx()
#ax2.scatter(Tlist, enthalpy, c='black', label="h_org, Tcommon=%d"%(species.thermo.coeffs[0]))
#ax2.plot(Tlist, enthalpy_fit, c='green', label="h_fit, Tcommon=1000")
#ax2.set_ylabel("h [J/kg]")

fig.legend(loc="lower right", bbox_to_anchor=(1,0), bbox_transform=ax1.transAxes)
plt.savefig("%s.png"%(species_name), dpi=500, bbox_inches='tight', pad_inches=0.1)
plt.close()

