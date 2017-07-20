# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:45:49 2017

@author: cssilva
"""

# ATMOS 1 - plotting the spectral features of individual functionals within molecules

import math
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum
import cPickle as pickle
import re

boltzmann = 1.38064852 * (10**(-23))
light_speed = 2.998 * 10**8
plank = 6.626070040 * (10**(-34))
cm1_joules = 5.03445 * 10**(22)

temperature = 300.0


# intensity_centre = 1.0
# freq_centre = 0.0
# bcon = float(10**(45)*(plank/(8*math.pi*math.pi*light_speed*5)))
# jmax = int(0.5 * (np.sqrt(((2*boltzmann*temperature)/(bcon*10**(-24))))) - 0.5)
# #
# for j in range(0, jmax):
#     dcon = (bcon * 10 ** (-3)) / (j + 1)
#     dcon_plus = (bcon * 10 ** (-3)) / (j + 2)
#     spacing = 2*bcon - ((4 * dcon_plus)*((j + 2)**3)) + ((4 * dcon)*((j + 1)**3))
#     intensity_j = intensity_centre * 5 * ((2 * j) + 1) * (10 ** (-2)) * np.e**(-(((bcon*10**(-24)) * (2 * j) * ((2*j) + 1)) / (boltzmann * temperature)))
#     position_j_pbranch = freq_centre - spacing
#     position_j_qbranch = freq_centre + spacing
#     print 'pbranch', j, position_j_pbranch, intensity_j
#     print 'qbranch', j, position_j_qbranch, intensity_j
# #




class Molecule:
    def __init__(self, code):
        self.code = code
        self.functionals = []

    def addFunctional(self, functional, number):
        self.functionals.append((functional, number))

    def contains(self, element_name):
        return self.code.contains(elements[element_name])

    def line_shapes(self):
        lines = []

        for functional_tuple in self.functionals:
            functional = functional_dictionary[functional_tuple[0]]
            for symmetry in functional.symmetries:
                for property in symmetry.properties:
                    x = np.linspace(property.low, property.high)
                    y = functional.line_function(x, property.frequency_average(), property.intensity.value)
                    lines.append((x, y))

        return lines

    def branches(self):
        branches = []

        for functional_tuple in self.functionals:
            functional = functional_dictionary[functional_tuple[0]]
            for symmetry in functional.symmetries:
                for property in symmetry.properties:
                    branches.append(self.prBranches(property))

        return branches

    def atom_count(self):
        atoms = len(re.sub(r"[^A-Z]+",'',self.code))

        return atoms

    def prBranches(self, property):
        pr_branch_x = []
        pr_branch_y = []

        bcon = float(plank / (8 * math.pi * math.pi * light_speed * self.atom_count() * 10**(-44)))
        jmax = int(np.sqrt((boltzmann * temperature)/(2 * plank * light_speed * bcon)) - 0.5)

        for j in range(0, jmax):
            dcon = (bcon * 10 ** (-3)) / (j + 1)
            dcon_plus = (bcon * 10 ** (-3)) / (j + 2)
            spacing = 2 * bcon - ((4 * dcon_plus) * ((j + 2) ** 3)) + ((4 * dcon) * ((j + 1) ** 3))
            intensity_j = property.intensity.value * ((2* j) + 1) * np.e**(-((plank * light_speed * bcon * j * (j + 1))/(boltzmann * temperature)))

            position_j_pbranch = property.frequency_average() - spacing
            position_j_rbranch = property.frequency_average() + spacing
            pr_branch_x.append(position_j_pbranch)
            pr_branch_y.append(intensity_j)
            pr_branch_x.append(position_j_rbranch)
            pr_branch_y.append(intensity_j)

        return (pr_branch_x, pr_branch_y)

    def average_points(self):
        points = []

        for functional_tuple in self.functionals:
            functional_code = functional_tuple[0]
            if functional_code in functional_dictionary.keys():
                functional = functional_dictionary[functional_code]
                for symmetry in functional.symmetries:
                    for property in symmetry.properties:
                        points.append((property.frequency_average(), property.intensity.value))
                    

        return points
    
    def high_and_low_frequencies(self):
        frequencies = []

        for functional_tuple in self.functionals:
            functional_code = functional_tuple[0]
            if functional_code in functional_dictionary.keys():
                functional = functional_dictionary[functional_code]
                for symmetry in functional.symmetries:
                    for property in symmetry.properties:
                        frequencies.append((property.low, property.high, property.intensity.value))
                    

        return frequencies
    

class Functional:
    def __init__(self, code, a = 1, b = 1):
        self.code = code
        self.a = a
        self.b = b
        self.symmetries = []

    def addSymmetry(self, symmetry):
        symmetry.functional = self
        self.symmetries.append(symmetry)

    def line_function(self, x, translateX, scaleY):
        print("Calculating graph for functional '" + self.code + "': using default f()")

        return (1/(1+pow((x - translateX),2))) * scaleY



# Specify a subclass of functional that has a different graphing function
class ExpFunctional(Functional):
    def line_function(self, x, translateX, scaleY):
        print("Calculating graph for functional '" + self.code + "': using exp f()")

        return (self.a * np.exp(-pow((x - translateX), 2))) * scaleY

class Symmetry:
    def __init__(self, type):
        self.type = type
        self.properties = []

    def addProperty(self, property):
        property.symmetry = symmetry
        self.properties.append(property)

class Property:
    def __init__(self, low, high, intensity):
        self.low = low
        self.high = high
        self.intensity = intensity

    def frequency_average(self):
        return np.mean([self.high,self.low])
#        return self.low + ((self.high - self.low) / 2)

class Intensity(Enum):
    weak = 1
    medium = 2
    strong = 3

    @classmethod
    def fromString(self, str):
        if str == 'strong':
            return Intensity.strong
        elif str == 'medium':
            return Intensity.medium
        elif str == 'weak':
            return Intensity.weak

# Load Functionals

# Example data
# COC C-O-C sbend 2500 2720 weak
# COC C-O-C abend 2800 2920 strong

functional_dictionary = {}
functional_data = open('func_temp', "r")

for line in functional_data:
    columns = line.split()

    code = columns[0]
    name = columns[1]
    symmetry_name = columns[2]
    low = float(columns[3])
    high = float(columns[4])
    intensity = Intensity.fromString(columns[5])

    if not functional_dictionary.has_key(code):
        # Extra logic for determining which type of functional we are dealing with.
        if "CF" in code:
            f = ExpFunctional(code, 2) # in this case we are setting 'a' for this functional's graphing function to 2.
        else:
            f = Functional(code)

        functional_dictionary[code] = f

    symmetry = Symmetry(symmetry_name)
    property = Property(low, high, intensity)

    functional_dictionary[code].addSymmetry(symmetry)
    symmetry.addProperty(property)

# Load Molecules

molecules = {}
#molecule_dictionary = pickle.load(open("dict_sorted_results_func_intra_test2_numbers.p", "rb"))
molecule_dictionary = pickle.load(open("dict_sorted_results_func_intra_table_part.p", "rb"))
#print 'test molecule', molecule_dictionary['[H]OC([H])(C)C']
print 'Molecule dictionary sample', molecule_dictionary.items()[:5]

print 'Number of molecules', len(molecule_dictionary.items())



for molecule_code, molecule_functionals in molecule_dictionary.iteritems():
    molecule = Molecule(molecule_code)

    for functional_tuple in molecule_functionals:
        if '#' not in functional_tuple[0]:
            functional_code = functional_tuple[0]
            functional_incidence = functional_tuple[1]
            molecule.addFunctional(functional_code, functional_incidence)

    molecules[molecule_code] = molecule


# Finds all the strong features whose average frequency is within a window
window_low = 450
window_high = 950
window_molecules = []
strong_window_molecules = []
for molecule_code in molecules:
    molecule = molecules[molecule_code]

#    for frequency in molecule.high_and_low_frequencies():
#            print frequency[0]

    
    window_filtered_list = list(filter(lambda x: window_low < x[0] < window_high, molecule.average_points()))
#    low_filtered_list = list(filter(lambda x: 600 < x[0] < 800, molecule.high_and_low_frequencies()[0]))
    
    for point in window_filtered_list:
  #      print point
        window_molecules.append(molecule.code)
        if point[1] == 3:
            strong_window_molecules.append(molecule_code)
#            print point
               
#    if len(window_filtered_list) > 0:
#        window_molecules.append(molecule.code)
#        print points
#        strong_filtered_list = list(filter(lambda y: y[1] == 3, points))
#        if len(strong_filtered_list) > 0:
#            strong_window_molecules.append(molecule_code)
        
#        print 'Molecule exits in window', molecule.code, (filtered_list)

print 'Window: ', window_low, ' to ', window_high
print 'Number of molecules in window', len(window_molecules)

print 'Number of strong molecules in window', len(strong_window_molecules)
print 'List of Molecules', strong_window_molecules


example = molecules['C(Cl)(CC(F)I)']
#print len(example.functionals)



# print example.atom_count()

# Plot points
xs, ys = zip(*example.average_points())
#print 'points', example.points()
#plt.plot(xs, ys, linestyle='None', marker='o', color='black', linewidth=2)

window = range(3250, 3450)
filtered_list = list(filter(lambda x: 1600 < x[0] < 1900, example.average_points()))
#print 'filter:',(filtered_list)


markerline, stemlines, baseline = plt.stem(xs, ys, '-')
plt.setp(baseline, 'color', 'r', 'linewidth', 1)

#Plot branches
#for line in example.branches():
#     x, y = line
#     markerline, stemlines, baseline = plt.stem(x, y, '-')
#     plt.setp(baseline, color='r', linewidth=1, marker='None')

#plt.xlim(1459,1459.5)
plt.show()

# Plot lines
#for line in example.lines():
#    x, y = line

#    plt.plot(x, y)


