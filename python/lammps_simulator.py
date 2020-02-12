'''
Simulation interface for LAMMPS using the Vashishta potential.

Prerequisites:
- re
- shutil
- os / subprocess
'''

class AutoSim:
    def __init__(self, substance):
        '''
        Initialize class
        
        Arguments:
        ----------
        substance   {str}   :   which potential to use. 
                                Candidates: 'SiO2', 'H2O', 'SiO2H2O'.
        '''
        from default_parameters import get_parameters
        self.substance = substance.lower()
        self.parameters, self.masses = get_parameters(self.substance)
            
    def set_parameters(self, parameters):
        '''
        This function overwrites the default parameters.
        
        Arguments:
        ----------
        params      {dct}   :   Nested dictionary with new parameters. Has to
                                be in the form of: 
                                {"comb1": {"param1": value1, "param2": value2, ...}, 
                                 "comb2": {"param1": value1, "param2": value2, ...},
                                 ...}.
        '''
        for comb, params in parameters.items():
            for param, value in params.items():
                self.parameters[comb][param] = value
        
    def ordered_parameter_string(self, params, param_suffices, list, string):
        '''
        Returning an ordered list of all the parameter values.
        
        Arguments:
        ----------
        params          {dct}       :   Dictionary with all parameters.
        param_suffices  {list(str)} :   List of all parameter suffices.
        list            {list}      :   Initial list to append parameter to.
        string          {str}       :   Initial string that will be extended.
        '''
        for suffix in param_suffices:
            list.append(params[suffix])
        for param in list:
            string += str(param) + 2 * " "
        return string + "\n"
        
    def append_type_to_file(self, name, params, filename):
        '''
        Append the actual parameter values to the parameter file.
        
        Arguments:
        ----------
        name            {str}       :   Name of interaction combo (e.g. "SiSiSi")
        params          {dct}       :   Dictonary with all parameters
        filename        {str}       :   Filename 
        '''
        # Split name
        from re import findall
        prefix_list = findall('[A-Z][^A-Z]*', name)
        params_line1 = ["H", "eta", "Zi", "Zj", "lambda1", "D", "lambda4"] # correctly ordered
        params_line2 = ["W", "rc", "B", "gamma", "r0", "C", "cos(theta)"]  # correctly ordered
        string_line1 = self.ordered_parameter_string(params, params_line1, prefix_list, "")
        string_line2 = self.ordered_parameter_string(params, params_line2, [], (len(name) + 6) * " ")
            
        with open(filename, 'a') as file:
            file.write("\n")
            file.write(string_line1)
            file.write(string_line2)
        
        
        
    def generate_parameter_file(self, filename="../data/dest.vashishta", 
                                      header_filename="../data/header.vashishta"):
        '''
        Generates input parameter file for the potential. The default
        parameters are the ones specified in Wang et al., so parameters
        that are not specified will fall back on these default parameters.
        
        Arguments:
        ----------
        filename        {str}   :   Filename
        header_filename {str}   :   Header file name
        '''
        # Add header to file
        from shutil import copyfile
        self.filename = filename
        copyfile(header_filename, filename)
        
        # Add parameters to file
        for name, params in self.parameters.items():
            self.append_type_to_file(name, params, filename)
        
            
    def modify_shell(self, read_data, input_script, extension):
        '''
        Modify shell
        '''
        element_string = ""
        masses = []
        for key, value in self.masses.items():
            element_string += key + " "
            masses.append(value)
        
        self.input_script = input_script
        f = open("../lammps/shell.in", "r")
        contents = f.readlines()
        f.close()
        
        if self.substance == "silica" or self.substance == "sio2":
            contents.insert(6, "variable extension equal {}\n".format(extension))
            contents.insert(7, "\n")
            contents.insert(8, "read_data ../data/silica_slab.data\n")
            contents.insert(9, "pair_style vashishta\n")
            contents.insert(10, "pair_coeff * * {} {}\n".format(self.filename, element_string))
            for i in range(len(masses)):
                contents.insert(11+i, "mass" + 12 * " " + str(i+1) + " " + str(masses[i]) + "\n")

            f = open(input_script, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
            
        elif self.substance == "water" or self.substance == "h2o":
            contents.insert(6, "variable extension equal {}\n".format(extension))
            contents.insert(7, "\n")
            contents.insert(8, "read_data ../data/water_lmps.data\n")
            contents.insert(9, "pair_style vashishta\n")
            contents.insert(10, "pair_coeff * * {} {}\n".format(self.filename, element_string))
            for i in range(len(masses)):
                contents.insert(11+i, "mass" + 12 * " " + str(i+1) + " " + str(masses[i]) + "\n")

            f = open(input_script, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
            
    def call_lammps(self, lammps_exec):
        '''
        Call LAMMPS.
        '''
        way = "os"
        if way == "os":
            from os import system
            call_string = "{} -in {}".format(lammps_exec, self.input_script)
            system(call_string)
        elif way == "subprocess":
            from subprocess import call
            call([lammps_exec, "-in", self.input_script])

        
    def simulate(self, read_data="../data/water_lmps.data",
                       lammps_exec="mpirun -n 4 lmp_mpi", 
                       input_script="../lammps/script.in", 
                       extension="0"):
        '''
        Run LAMMPs simulation with the parameters. 
        
        Arguments:
        ----------
        where_to_simulate   :   specify where to simulate. Cluster should
                                be an option.
        '''
    
        self.modify_shell(read_data, input_script, extension)
        self.call_lammps(lammps_exec)
        return None
        
    def estimate_boiling_temperature(self):
        '''
        Reading LAMMPs log file and finds the boiling point.
        '''
        boiling_temperature=100
        return boiling_temperature
        
    def estimate_boiling_enthalpy(self):
        '''
        Reading LAMMPs log file and finds the enthalpy.
        '''
        enthalpy = 10
        return enthalpy
        
    def grid_search(self, parameter_span):
        '''
        Do a grid search with several parameters.
        '''
        return None
        
if __name__ == "__main__":
    params = {"OOSi" : {"B" : 2.3, "H" : 700}}
    
    obj = AutoSim("water")
    #obj.set_parameters(params)
    obj.generate_parameter_file(filename="../data/H2O.vashishta")
    obj.simulate("desk")
    
    
    
