'''
This script was written to analyze the output files from LAMMPs.
This includes the dump files (containing information about all particles)
and log files (containing information about the whole system). The former
is treated by the Dump class, while the latter is treated by the Log class.

Author: Even Marius Nordhagen
'''

import numpy as np
import matplotlib.pyplot as plt
plt.style.use("bmh")
plt.rcParams["font.family"] = "Serif"   # Font

class Dump:
    ''' 
    Analyzing dump files, containing particle-specific information. 
    '''
    def __init__(self, filename, info_lines=9, particle_line=3):
        '''
        Initializing the class and loading the file.
        
        Arguments:
        ----------
        filename        {str}       : Which dump file to read
        info_lines      {int}       : Number of lines storing information
        particle_line   {int}       : The line that gives the number of particles
        '''
        print("Loading data. For large files, this might takes a while.")
        f = open(filename,'rt')
        self.particles = int(f.readlines()[particle_line])
        f.close()
        num_lines = sum(1 for line in open(filename))
        length = self.particles + info_lines
        self.steps = int(num_lines / length)
        data = []
        from tqdm import tqdm
        for t in tqdm(range(self.steps)):
            data.append(np.loadtxt(open(filename,'rt').readlines()[t * length + info_lines: (t+1) * length]))
        self.data = np.asarray(data)


    def plot_position_distribution(self, steps=[0], show=False, save=False):
        '''
        Plots a histogram of the radius of the particles. 
        
        Arguments:
        ----------
        steps           {list}      : The timesteps of interest
        show            {bool}      : Show plot yes/no (True/False). Default: False
        save            {bool}      : Save plot yes/no (True/False). Default: False
        '''
        for i in steps:
            positions = self.data[i][:, 2:5]
            radius = np.linalg.norm(positions, axis=1)
            plt.hist(radius, 100, density=True, facecolor='b', alpha=0.75)
            plt.xlabel('Radius')
            plt.ylabel('Density')
            plt.grid()
            if show: plt.show()
            if save: plt.savefig('../fig/position_distribution_{}.png'.format(t))


    def plot_velocity_distribution(self, steps=[0], show=False, save=False):
        '''
        Plots a histogram of the speed of all particles at selected timesteps.
        
        Arguments:
        ----------
        steps           {list(int)} : The timesteps of interest
        show            {bool}      : Show plot yes/no (True/False). Default: False
        save            {bool}      : Save plot yes/no (True/False). Default: False
        '''
        for i in steps:
            velocity = self.data[i][:, 5:8]
            speed = np.linalg.norm(velocity, axis=1)
            plt.hist(speed, 100, density=True, facecolor='b', alpha=0.75)
            plt.xlabel('Speed')
            plt.ylabel('Density')
            plt.grid()
            if show: plt.show()
            if save: plt.savefig('../fig/velocity_distribution_{}.png'.format(t))
            
    def plot_diffusion(self, show=False, save=False):
        '''
        Plot the diffusion and estimate the diffusion constant.
        '''
        pos = self.data[:,:,2:5]
        initial_pos = pos[0]
        diff = pos - initial_pos[np.newaxis]
        res = np.einsum('ijk,ijk->i',diff,diff)
        plt.plot(res)
        if save: plt.savefig("../fig/diffusion.png")
        if show: plt.show()
        
        t = np.linspace(0, 500, len(res))
        D = np.divide(res, t) / 6
        print(D)
        
    def plot_radial(self, show=False, save=False, L=3):
        '''
        Plot the radial distribution function as a function of relative distance.
        '''
        pos = self.data[:,:,2:5]
        from scipy.spatial import distance_matrix
        d = distance_matrix(pos, pos)
        print(d)
        d_ = d.flatten()
        bins = np.linspace(0, L/2, 1000)
        inds = np.digitize(d.flatten(), bins)
        plt.plot(inds)
        if show: plt.show()
       
        
class Log:
    ''' Analyzing log files, containing system information. '''
    def __init__(self, filename, ignore_first=0):
        '''Initialize class by reading log file.
        
        Arguments:
        ----------
        filename            {str}   : String with file to load.
        ignore_first        {int}   : Ignore equilibriation fixes.
        '''
        self.read_log_file(filename, ignore_first)

    def read_log_file(self, filename, ignore_first):
        ''' Reading log file by going through file line after line. 
        
        Arguments:
        ----------
        filename            {str}   : String with file to load.
        ignore_first        {int}   : Ignore equilibriation fixes.
        '''
        f = open(filename, "r")

        self.timestep = 0.005        # Default
        self.mass = 1                # Default
        self.lst = []

        read = False                # True if the line should be read
        for i, line in enumerate(f.readlines()):
            # Search for variables
            if line.startswith("Step"):
                self.variables = line.split()
                num_variables = len(self.variables)
                self.lst.append([[] for _ in range(num_variables)])
                read = True
            elif line.startswith("Loop time of"):
                read = False
            elif read:
                #print(line)
                strings = line.split()
                for j, string in enumerate(strings):
                    self.lst[-1][j].append(float(line.split()[j]))
            # Search for timestep
            elif line.startswith("timestep"):
                self.timestep = line.split()[1]
            # Search for mass
            elif line.startswith("mass"):
                self.mass = float(line.split()[1])
        f.close()
        self.array = np.hstack(self.lst[ignore_first:])

    def categorize(self):
        '''
        Trying to sort data into categories automatically.
        '''
        for i, variable in enumerate(self.variables):
            if variable == "Step":
                self.times = np.asarray(self.array[i]) * self.timestep
            elif variable == "Temp":
                self.temp = self.array[i]
            elif variable == "Press":
                self.pres = self.array[i]
            elif variable == "E_pair":
                self.E_int = self.array[i]
            elif variable == "E_mol":
                self.E_mol = self.array[i]
            elif variable == "TotEng":
                self.E_tot = self.array[i]
            elif variable == "Volume":
                self.vol = self.array[i]
            elif variable == "PotEng":
                self.pot = self.array[i]
            elif variable == "KinEng":
                self.kin = self.array[i]
                
    def find(self, key):
        '''
        Search for a category (Step, Temp, Press etc...). If the 
        keyword exists, if returns the associated array containing
        the quantity as a function of timesteps.
        
        Arguments:
        ----------
        key         {str}   : String containing keyword.
        '''
        array = None
        for i, variable in enumerate(self.variables):
            if variable == key:
                array = self.array[i]
        if array is None:
            raise KeyError("No category named {} found.".format(key))
        else:
            return np.array(array)
            
    def step2time(self, steps):
        '''
        Converting an array of steps to actual times in 
        SI units, based on the timestep found in the file.
        '''
        self.time = np.asarray(steps) * self.timestep
        return self.time

    def convert2si(self):
        from lj_units import Convert_LJ
        convert = Convert_LJ(self.mass, 1, 1)

        self.times = convert.time(self.times)
        self.temp = convert.temp(self.temp)
        self.pres = convert.pressure(self.pres)
        self.E_int = convert.energy(self.E_int)
        self.E_mol = convert.energy(self.E_mol)
        self.E_tot = convert.energy(self.E_tot)
        self.pot = convert.energy(self.pot)
        self.kin = convert.energy(self.kin)
        self.vol = convert.volume(self.vol)

    def plot_temperature(self, show=False, save=False):
        plt.plot(self.times, self.temp)
        plt.xlabel("Time [s]")
        plt.ylabel("Temperature [K]")
        if save: plt.savefig("../fig/temp_{}.png".format(self.timestep))
        if show: plt.show()
        
    def plot_pressure(self, show=False, save=False):
        plt.plot(self.times, self.pres)
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [m$^{-3}$]")
        if save: plt.savefig("../fig/pressure_{}.png".format(self.timestep))
        if show: plt.show()

    def plot_energy(self, show=False, save=False, total=True, kin=False, pot=False):
        plt.plot(self.times, self.E_tot, label="Total")
        if kin: plt.plot(self.times, self.kin, label="Kinetic")
        if pot: plt.plot(self.times, self.pot, label="Potential")
        plt.xlabel("Time [s]")
        plt.ylabel("Energy [J]")
        if save: plt.savefig("../fig/energy_{}.png".format(self.timestep))
        if show: plt.show()
        
    def plot_volume(self, show=False, save=False):
        plt.plot(self.times, self.vol)
        plt.xlabel("Time [s]")
        plt.ylabel("Volume [m$^3$]")
        if save: plt.savefig("../fig/volume_{}.png".format(self.timestep))
        if show: plt.show()
        
    def plot_pressure(self, show=False, save=False):
        plt.plot(self.temp, self.pres)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Pressure [m$^{-3}$]")
        if save: plt.savefig("../fig/temp_pres_{}.png".format(self.timestep))
        if show: plt.show()
