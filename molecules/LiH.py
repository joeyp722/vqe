from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
import numpy as np
from scipy.io import savemat
import os
import math

# Parameters
name = 'LiH'
distance_interval=0.05
min_distance=0.10
max_distance=4.50

occupied_indices=[0]
active_indices=[1, 2, 3]

states = np.array([0, 1, 2, 0, 1, 2])
spins = np.array([0, 0, 0, 1, 1, 1])

run_mp2=True
run_cisd=True
run_ccsd=True
run_fci=True

path=''


# Create target Directory if don't exist

directory_name=name

if not os.path.exists(directory_name):
    os.mkdir(directory_name)
    print("Directory " , directory_name ,  " Created ")
else:    
    print("Directory " , directory_name ,  " already exists")



# Calculate used variables

max_distance_index=math.ceil((max_distance-min_distance+distance_interval)/distance_interval)
matrix_size=2*len(active_indices)

for distance_index in range(1,max_distance_index+1):

    distance=(distance_index-1)*distance_interval+min_distance

    geometry = [['Li', [0, 0, 0]], ['H', [0, 0, distance]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0

    molecule = MolecularData(geometry, basis, multiplicity, charge)
    molecule = run_psi4(molecule, run_mp2=run_mp2, run_cisd=run_cisd, run_ccsd=run_ccsd, run_fci=run_fci)
    core_electron_matrix, oneelectronintegrals, twoelectronintegrals =molecule.get_active_space_integrals(occupied_indices=occupied_indices, active_indices=active_indices)
 
    one_electron_matrix = np.zeros((matrix_size, matrix_size))
    two_electron_matrix = np.zeros((matrix_size, matrix_size, matrix_size, matrix_size))

    for i in range(matrix_size):
        for j in range(matrix_size):
            if spins[i] == spins[j]:
                one_electron_matrix[i, j] = oneelectronintegrals[states[i], states[j]]

    for i in range(matrix_size):
        for j in range(matrix_size):
            for k in range(matrix_size):
                for l in range(matrix_size):
                    if spins[i] == spins[l] and spins[j] == spins[k]:
                        two_electron_matrix[i, j, k, l] = twoelectronintegrals[states[i], states[j], states[k], states[l]]

    savemat(directory_name+'/core_electron_matrix' + str(distance_index) + '.mat', mdict={'matrix': core_electron_matrix})
    savemat(directory_name+'/one_electron_matrix'+str(distance_index)+'.mat', mdict={'matrix': one_electron_matrix})
    savemat(directory_name+'/two_electron_matrix' + str(distance_index) + '.mat', mdict={'matrix': two_electron_matrix})

    print(distance)
