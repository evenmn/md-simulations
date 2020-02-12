from lammps_simulator import AutoSim
from post_process import Log
from visualize_ovito import visualize
from pack_water import WaterPack

import os
import numpy as np
import matplotlib.pyplot as plt

number = 2000       # Number of water molecules
density = 0.9966    # Density of water at 25°C and 1atm

packer = WaterPack(number)
packer.den2len(density)
packer.packmol_gen(pbc=1.0)
packer.packmol_run()
packer.xyz2lmp("../data/water_lmps.data")

Z_Hs = [0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60]
thetas = np.arange(95, 106)

for Z_H in Z_Hs:
    for theta in thetas:

        path = "../data/ZH{}_theta{}/".format(Z_H, theta)
 
        try:
            # Create target Directory
            os.mkdir(path)
            print("Directory ", path, " created")
        except FileExistsError:
            print("Directory ", path, " already exists")

        Z_O = - 2 * Z_H
        params = {"HHH" : {"Zi" : Z_H, "Zj" : Z_H},
                  "OOO" : {"Zi" : Z_O, "Zj" : Z_O},
                  "OHH" : {"Zi" : Z_O, "Zj" : Z_H, "cos(theta)" : np.cos(theta)},
                  "HOO" : {"Zi" : Z_H, "Zj" : Z_O}}
            
        sim = AutoSim("water")
        sim.set_parameters(params)
        sim.generate_parameter_file(filename=path + "H2O.vashishta.Zh_{}".format(Z_H))
        sim.simulate("desk", lammps_exec="mpirun -n 4 lmp_mpi", extension=Z_H)

        logger = Log(path + "log.{}".format(Z_H))
        energy = logger.find("TotEng")
        enthal = logger.find("Enthalpy")
        temp = logger.find("Temp")
        time = logger.find("Time")
        pres = logger.find("Press")
        density = logger.find("Density")

        plt.figure()
        plt.plot(temp, energy)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Total energy [eV]")
        plt.savefig(path + "temp_eng_{}.png".format(Z_H))
        
        plt.figure()
        plt.plot(temp, enthal)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Enthalpy [eV]")
        plt.savefig(path + "temp_enth_{}.png".format(Z_H))
        
        plt.figure()
        plt.plot(time, temp)
        plt.xlabel("Time [ps]")
        plt.ylabel("Temperature [K]")
        plt.savefig(path + "time_temp_{}.png".format(Z_H))
        
        plt.figure()
        plt.plot(time, energy)
        plt.xlabel("Time [ps]")
        plt.ylabel("Total energy [eV]")
        plt.savefig(path + "time_eng_{}.png".format(Z_H))
        
        plt.figure()
        plt.plot(time, pres)
        plt.xlabel("Time [ps]")
        plt.ylabel("Pressure [Bar]")
        plt.savefig(path + "time_pres_{}.png".format(Z_H))
        
        plt.figure()
        plt.plot(time, density)
        plt.xlabel("Time [ps]")
        plt.ylabel("Density [g/cm³]")
        plt.savefig(path + "time_density_{}.png".format(Z_H))
        
        visualize(path + "minimize_300K.{}.data".format(Z_H),
                  path + "minimize_300K.{}.png".format(Z_H))
                  
        visualize(path + "water_300K.{}.data".format(Z_H),
                  path + "water_300K.{}.png".format(Z_H))
                  
        visualize(path + "water_301K.{}.data".format(Z_H),
                  path + "water_301K.{}.png".format(Z_H))
                  
        visualize(path + "vapor_450K.{}.data".format(Z_H),
                  path + "vapor_450K.{}.png".format(Z_H))
