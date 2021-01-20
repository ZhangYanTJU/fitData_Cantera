import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import math

# reference about fitting thermo data: https://github.com/Upstream-Applied-Science/canteraJanaf
# reference about writing thermo file: https://github.com/jiweiqi/Cantera2Chemkin/blob/master/soln2ck.py

def build_nasa(nasa_coeffs, row):
            """
            Creates string of nasa polynomial coefficients
            :param nasa_coeffs
                cantera species thermo coefficients object
            :param row
                which row to write coefficients in
            """
            line_coeffs = ''
            lines = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14]]
            line_index = lines[row-2]
            for ix, c in enumerate(nasa_coeffs):
                if ix in line_index:
                    if c >= 0:
                        line_coeffs += ' '
                    line_coeffs += str('{:.8e}'.format(c))
            return line_coeffs


def write(gas, thermo_fileName_new, plot):
    """Function to write cantera gas object to inp file.
    :param gas:
        gas: Cantera gas object,
        thermo_fileName_new: Name of converted thermo file
        plot: whether plot or not
    """
    with open(thermo_fileName_new, 'w+') as f:
        #f.write('THERMO ALL' + '\n' +
        #        '   300.000  1000.000  5000.000' +'\n')

        Tlow = 300
        Tcommon = 1000
        Thigh = 5000
        deltaT = 100
        T_low = np.arange(Tlow,Tcommon,deltaT)
        T_high = np.arange(Tcommon,Thigh,deltaT)
        temp_range = str('{0:.3f}'.format(Tlow)) + '  ' + str('{0:.3f}'.format(Tcommon)) + '  ' + '{0:.3f}'.format(Thigh)


        f.write('THERMO ALL' + '\n' +
                '   ' + temp_range + '\n')

        #write data for each species in the Solution object
        for sp_index in range(len(gas.species_names)):
            species = gas.species(sp_index)
            species_name = gas.species_name(sp_index)
            molecular_weight = gas.molecular_weights[sp_index]
            nasa_coeffs = species.thermo.coeffs # orininal coeffs, size = 15, [0]=Tcommon, [1-7] is high coeffs, [8-14] is low coeffs
            temp_range = str('{0:.3f}'.format(Tlow)) + '  ' + str('{0:.3f}'.format(Thigh)) + '  ' + '{0:.3f}'.format(Tcommon)

            if (nasa_coeffs[0] != Tcommon):
                print("fitting species: %s"%(species_name), "from ",nasa_coeffs[0], "to ",Tcommon)
                cp_low = [species.thermo.cp(T)/molecular_weight for T in T_low] #  [J/(kmol·K)] / [kg/kmol] = [J/(kg·K)]
                cp_high = [species.thermo.cp(T)/molecular_weight for T in T_high] # J/kg

                cp_low_polynomials = np.polyfit(T_low, cp_low, 4)
                cp_high_polynomials = np.polyfit(T_high, cp_high, 4)

                R = ct.gas_constant/molecular_weight  # [J/(kmol·K)] / [kg/kmol] = [J/(kg·K)]

                #reverse order and divide by R (np.ployfit has calculated coeffs*R so we must divide by R to get just the coeffs)
                a1_a5_low = np.flip(cp_low_polynomials,0)/R
                a1_a5_high = np.flip(cp_high_polynomials,0)/R

                T_ref = 298.15
                h = species.thermo.h(T_ref)/molecular_weight # [J/kmol] / [kg/kmol] = [J/kg]
                s = species.thermo.s(T_ref)/molecular_weight # [J/kmol/K] / [kg/kmol] = [J/(kg·K)]

                a6_low = h/R - np.sum(np.array([a1_a5_low[i]*T_ref**(i+1)/(i+1) for i in range(5)]))
                a6_high = h/R - np.sum(np.array([a1_a5_high[i]*T_ref**(i+1)/(i+1) for i in range(5)]))

                a7_low = s/R - np.sum(np.array([a1_a5_low[i]*T_ref**(i)/(i) for i in range(1,5)]))-a1_a5_low[0]*math.log(T_ref)
                a7_high = s/R - np.sum(np.array([a1_a5_high[i]*T_ref**(i)/(i) for i in range(1,5)]))-a1_a5_high[0]*math.log(T_ref)

                nasa_coeffs = np.concatenate([[Tcommon], a1_a5_high, [a6_high], [a7_high], a1_a5_low, [a6_low], [a7_low]])


                if plot:
                    Tlist = np.concatenate([T_low, T_high])
                    cp = [species.thermo.cp(T)/molecular_weight for T in Tlist] # J/kg
                    cp_fit_low = [ np.sum(np.array([a1_a5_low[i]*T**(i) for i in range(5)])) for T in T_low]
                    cp_fit_high = [ np.sum(np.array([a1_a5_high[i]*T**(i) for i in range(5)])) for T in T_high]
                    cp_fit = np.concatenate([cp_fit_low, cp_fit_high])

                    a1_a5_org_low = species.thermo.coeffs[8:13]
                    a1_a5_org_high = species.thermo.coeffs[1:6]
                    cp_bad_low = [ np.sum(np.array([a1_a5_org_low[i]*T**(i) for i in range(5)])) for T in T_low]
                    cp_bad_high = [ np.sum(np.array([a1_a5_org_high[i]*T**(i) for i in range(5)])) for T in T_high]
                    cp_bad = np.concatenate([cp_bad_low, cp_bad_high])


                    enthalpy = [species.thermo.h(T)/molecular_weight for T in Tlist] # J/kg
                    enthalpy_fit_low = [ a6_low + np.sum(np.array([a1_a5_low[i]*T**(i+1)/(i+1) for i in range(5)])) for T in T_low]
                    enthalpy_fit_high = [ a6_high + np.sum(np.array([a1_a5_high[i]*T**(i+1)/(i+1) for i in range(5)])) for T in T_high]
                    enthalpy_fit = np.concatenate([enthalpy_fit_low, enthalpy_fit_high])

                    fig = plt.figure()

                    ax1 = fig.add_subplot()
                    ax1.scatter(Tlist, cp, c='blue', label="cp_org, Tcommon=%d"%(species.thermo.coeffs[0]))
                    ax1.plot(Tlist, cp_fit*R, c='red', label="cp_fit, Tcommon=1000")
                    ax1.plot(Tlist, cp_bad*R, c='purple', label="cp_bad")
                    ax1.set_xlabel("T [K]")
                    ax1.set_ylabel("cp [J/kg]")

                    ax2 = plt.twinx()
                    ax2.scatter(Tlist, enthalpy, c='black', label="h_org, Tcommon=%d"%(species.thermo.coeffs[0]))
                    ax2.plot(Tlist, enthalpy_fit*R, c='green', label="h_fit, Tcommon=1000")
                    ax2.set_ylabel("h [J/kg]")

                    fig.legend(loc="upper left", bbox_to_anchor=(0,1), bbox_transform=ax1.transAxes)
                    plt.savefig("%s.png"%(species_name), dpi=500, bbox_inches='tight', pad_inches=0.1)
                    plt.close()



            species_comp = ''
            for atom in species.composition:
                species_comp += '{:<3}'.format(atom)
                species_comp += '{:<2}'.format(str(int(species.composition[atom])))

            species_phase = 'G'

            line_1 = (
                    '{:<18}'.format(species_name) +
                    '{:<6}'.format('    ') +
                    '{:<20}'.format(species_comp) +
                    '{:<4}'.format(species_phase) +
                    '{:<31}'.format(temp_range) +
                    '{:<1}'.format('1') +
                    '\n')
            f.write(line_1)

            line_2_coeffs = build_nasa(nasa_coeffs, 2)
            line_2 = line_2_coeffs  + '    2\n'
            f.write(line_2)

            line_3_coeffs = build_nasa(nasa_coeffs, 3)
            line_3 = line_3_coeffs + '    3\n'
            f.write(line_3)

            line_4_coeffs = build_nasa(nasa_coeffs, 4)
            line_4 = line_4_coeffs + '                   4\n'
            f.write(line_4)

        f.write('END\n')



if __name__ == '__main__':
    mech='mech.cti'
    gas = ct.Solution(mech)
    thermo_fileName_new = 'thermoCorrected.dat'
    print(thermo_fileName_new)
    plot=True # True or False
    write(gas, thermo_fileName_new, plot)