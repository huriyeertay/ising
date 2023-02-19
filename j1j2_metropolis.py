#%%
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from icecream import ic
import os

'''
This file is created by Huriye Ertay for the purpose of studying Ising Lattice with respect to changing T and J models
'''

#%%

#Creating classes to allow easy manipulation of file
class IsingLattice:

    def __init__(self, lattice_size, j1, j2):
        self.lattice_size = lattice_size
        self.num_sites = lattice_size*lattice_size
        self.j1 = j1
        self.j2 = j2

        lattice_state = np.random.randint(2, size=(lattice_size,lattice_size))
        self.lattice_state = np.where(lattice_state==0,-1,1)

    def plot_state(self):
        plt.imshow(self.lattice_state)
        plt.axis("off")

    def flip_spin(self,i,j):
        self.lattice_state[i,j]*= -1

    def print_info(self):
        print("Lattice Size = ",self.lattice_size)
        print("j1 = ",self.j1)
        print("j2 = ",self.j2)
        print("Total number of spins = ",self.num_sites)

    def spin_energy(self,i,j):
        # Define nearest neighbour interaction
        l=self.lattice_size
        spin_ij = self.lattice_state[i,j]
        nearest_up = self.lattice_state[(i-1)%l,j]
        nearest_down = self.lattice_state[(i+1)%l,j]
        nearest_left = self.lattice_state[i,(j-1)%l]
        nearest_right = self.lattice_state[i,(j+1)%l]

        sum_nearest_neighbour = nearest_down + nearest_left + nearest_right + nearest_up

        # Define second nearest neighbour interaction
        right_up = self.lattice_state[(i-1)%l,(j+1)%l]
        right_down = self.lattice_state[(i+1)%l,(j-1)%l]
        left_up = self.lattice_state[(i-1)%l,(j-1)%l]
        left_down = self.lattice_state[(i+1)%l,(j+1)%l]

        sum_second_nearest_neighbour = right_down + right_up + left_down + left_up
        # Define total Ising energy of spin E = j1 * sum(no of nearest neighbour) + j2 * sum(no of second nearest neigbour)
        if self.j2 == 0:
            spin_energy = (-self.j1 * sum_nearest_neighbour * spin_ij)
        else:
            spin_energy = (-self.j1 * sum_nearest_neighbour * spin_ij) + (-self.j2 * sum_second_nearest_neighbour * spin_ij)

        return spin_energy
    
    def total_energy(self):
        #total energy of all spins, but we divide it to no of spins at the end to normalize the energy, we also divide by 2 two to prevent repetition
        total_energy = 0
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                total_energy += self.spin_energy(i,j)
        
        return total_energy/2.0/self.num_sites
    
    def spin_mag(self):
    #Magnetization is 1/number of spins * total of spins (include nominal value of state (-1,1))
        return self.lattice_state.sum() / self.num_sites


#Metropolis Algorithm 
def sweep_lattice(lattice:IsingLattice,T):
    '''
    Scan the lattice once ly, flip the spin if Metropolis criterion holds: if energy<=0 or generated random number<np.exp(-energy_change/temperature)
    '''

    flip_count = 0
    no_flip_count = 0

    for x in range(lattice.num_sites):
        i, j = np.random.randint(lattice.lattice_size, size=2)
        energy_initial = lattice.spin_energy(i,j)
        lattice.flip_spin(i,j)
        energy_flipped = lattice.spin_energy(i,j)
        exchange_energy = energy_flipped - energy_initial
        lattice.flip_spin(i,j)
        if exchange_energy <= 0 or \
            np.random.random()<=np.exp(-exchange_energy/T):
            lattice.flip_spin(i,j)

#Thermalize divides the T into many pieces for sweeping the lattice by decreasing the rate of change of T
def thermalize(lattice:IsingLattice,sweep_num,T_hot,T_cold):
    for T in np.linspace(T_hot, T_cold, sweep_num):
        sweep_lattice(lattice, T)

#Gathering and Saving Data Necessary
def mc_simulation(lattice:IsingLattice, T:np.array, num_sweeps, save_data = True):
    energy = []
    magnetization = []
    lattice_states = []
    file_name = f'L_{lattice.lattice_size}_J1_{lattice.j1:.1f}_J2_{lattice.j2:.3f}'
    os.makedirs(file_name, exist_ok=True)

    thermalize(lattice, 1000, 100, T[0])
    # Get data for T[0]
    sweep_lattice(lattice, T[0])
    energy.append(lattice.total_energy())
    magnetization.append(lattice.spin_mag())
    lattice_states.append(lattice.lattice_state)
   
    if save_data:
        plt.imshow(lattice.lattice_state)
        plt.axis('off')
        plt.savefig(file_name + "/" + file_name + f'_T_{T[0]:.3f}.png')
        plt.close()

    for i in tqdm(range(T.size-1)):
        thermalize(lattice,1000, T[i], T[i+1])
        # Now we want to collect data
        sweep_lattice(lattice, T[i+1])
        energy.append(lattice.total_energy())
        magnetization.append(lattice.spin_mag())
        lattice_states.append(lattice.lattice_state)
        if save_data:
            plt.imshow(lattice.lattice_state)
            plt.axis('off')
            plt.savefig(file_name + "/" +file_name + f'_T_{T[i+1]:.3f}.png')
            plt.close()

    energy = np.asarray(energy)
    magnetization = np.asarray(magnetization)
    lattice_states = np.asarray(lattice_states)
    print('Simulation done.')
    if save_data:
        print('Saving data...')
        np.savetxt(file_name + "/" +file_name + '_energy.txt', energy, delimiter=',')
        np.savetxt(file_name + "/" +file_name + '_magnetization.txt', magnetization, delimiter=',')
    return energy, magnetization, lattice_states
#%%
#Set Lattice Parameters
test_lattice_2 = IsingLattice(100, 1, 0)
T = np.arange(4,0.95,-0.05)

#initiate sequence
ene, mag, lat = mc_simulation(test_lattice_2,T,1000, save_data=True)
#%%
plt.scatter(T,ene)
plt.scatter(T,np.abs(mag))
