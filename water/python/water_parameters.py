from lammps_simulator import AutoSim
from post_process import Log
from visualize_ovito import visualize
from pack_water import WaterPack

import os
import numpy as np
import matplotlib.pyplot as plt
plt.rc('figure', max_open_warning = 0)      # Avoiding RecursionError

number = 2000       # Number of water molecules
density = 0.9966    # Density of water at 25°C and 1atm
read_data = "../data/water_lmps.data"

# First pack data
packer = WaterPack(number)
packer.den2len(density)
packer.packmol_gen(pbc=1.0)
packer.packmol_run()
packer.xyz2lmp(read_data)

Z_Hs = [0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60]
thetas = np.arange(95, 106)
Bs = np.arange(20, 82, 2)

Z_Hs = [0.50]
thetas = [100]
Bs = [40]
Ds = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]

for Z_H in Z_Hs:
    for theta in thetas:
        for B in Bs:
            for D in Ds:
                path = "../data/ZH{}_theta{}_B{}_D{}/".format(Z_H, theta, B, D)
         
                try:
                    # Create target Directory
                    os.mkdir(path)
                    print("Directory ", path, " created")
                except FileExistsError:
                    print("Directory ", path, " already exists")

                Z_O = - 2 * Z_H
                params = {"HHH" : {"Zi" : Z_H, "Zj" : Z_H},
                          "OOO" : {"Zi" : Z_O, "Zj" : Z_O},
                          "OHH" : {"Zi" : Z_O, "Zj" : Z_H, "cos(theta)" : np.cos(np.deg2rad(theta)), "B" : B, "D" : D},
                          "HOO" : {"Zi" : Z_H, "Zj" : Z_O, "D" : D}}
                    
                sim = AutoSim("water")
                sim.set_parameters(params)
                sim.generate_parameter_file(filename=path + "H2O.vashishta")
                sim.simulate(read_data=read_data, lammps_exec="mpirun -n 4 lmp_mpi", path=path)

                logger = Log(path + "log.data", ignore_first=3)
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
                plt.savefig(path + "temp_eng.png")
                
                plt.figure()
                plt.plot(temp, enthal)
                plt.xlabel("Temperature [K]")
                plt.ylabel("Enthalpy [eV]")
                plt.savefig(path + "temp_enth.png")
                
                plt.figure()
                plt.plot(time, temp)
                plt.xlabel("Time [ps]")
                plt.ylabel("Temperature [K]")
                plt.savefig(path + "time_temp.png")
                
                plt.figure()
                plt.plot(time, energy)
                plt.xlabel("Time [ps]")
                plt.ylabel("Total energy [eV]")
                plt.savefig(path + "time_eng.png")
                
                plt.figure()
                plt.plot(time, pres)
                plt.xlabel("Time [ps]")
                plt.ylabel("Pressure [Bar]")
                plt.savefig(path + "time_pres.png")
                
                plt.figure()
                plt.plot(time, density)
                plt.xlabel("Time [ps]")
                plt.ylabel("Density [g/cm³]")
                plt.savefig(path + "time_density.png")
                
                visualize(path + "minimize_300K.data",
                          path + "minimize_300K.png")
                          
                visualize(path + "water_after_nvt.data",
                          path + "water_after_nvt.png")
                          
                visualize(path + "water_after_npt.data",
                          path + "water_after_npt.png")
                          
                visualize(path + "water_final.data",
                          path + "water_final.png")
                          
                visualize(path + "vapor_450K.data",
                          path + "vapor_450K.png")
