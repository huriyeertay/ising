# Ising

### Data Generation : 

    - [X] Create IsingLattice Class
    - [X] Implement Metropolis Algorithm
    - [X] Implement data generation structure for many samples at once
    - [] Organize data generation directory tree appropriately
    - [] Generate data for 22 $J_2$ values given by reference.

### Data Analysis : 
    - [] Make a function to read energy and magnetization, given $J_2$ from the csv files
    - [] Plot average energy graph for each $J_2$ data
    - [] Plot average absolute magnetization for each $J_2$ value
    - [] Plot the variances of energy and magnetizations
    - [] Get the maximum values of the variances. They are the desired $T_c$ values/
    - [] From the $T_c$ values, plot $(J_2,T_c)$ graph i.e. phase diagram
to run exec :

```shell
python j1j2_metropolis.py 
```
