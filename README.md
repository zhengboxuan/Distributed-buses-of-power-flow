# Distributed-buses-of-power-flow
    This project modifies two standard optimal power flow algorithms to solve the distributed variation of the optimal power flow
    The algorithms are implemented in a Matlab function.
    The main function is 'Distribute.m', the basic function is 'solveNR_Distribute.m' and 'solveFDPF_Distribute.m'.

## Setup
    Clone the code and put them into the folder, create a new matlab to use the function of 'distribute'
    
## Package
    Note that the all the code shown above requires installation of the Matpower package.
    
## Distribute Slack Variants of Power Flow Algorithms
    The power flow problem discussed in class uses a single slack bus with specified voltage magnitude and phase angle (Vslack = 1 and θslack = 0°) and unspecified active and reactive power injections (Pslack and Qslack). A more general power flow formulation called a distributed slack can better model the actual behavior of typical power systems.
    
## Parameters
### Distribute(mpc,method,part) solves the Distributed Slack variation of the OPF, using either the Fast-Decoupled Power Flow or the Newton-Raphson algorithm. 
Inputs: `mpc` - MATPOWER case <br>
        `method` - algorithm used to solve OPF <br>
                 `1`: Fast-Decoupled Power Flow <br>
                 `2`: Newton-Raphson <br>
        `part` - method used to compute participation factors <br>
                 `1`: participation factors based on the cost function. <br>
                 `2`: participation factors based on the generation capacity of the buses. <br>
Output: `res` - result structure with fields `V` (voltage phasors), `success` (1 if convergence achieved and 0 if failed), `et` <br> (elapsed  time), and `nither` (number of iterations).
### solveNR_Distribute (mpc, method) solve the power flow problem basing on the NR method.
Inputs: `mpc` - MATPOWER case <br>
        `method` - method used to compute participation factors<br>
                 `1`: participation factors based on the cost function.<br>
                 `2`: participation factors based on the generation capacity of the buses.<br>
Output: `V` (Amplitude of voltage), `success` (1 if convergence achieved and 0 if failed), `et` (elapsed time), and `nither` <br>
(number of iterations).
### solveFDPF_Distribute (mpc, method) solve the power flow problem basing on the FDPF method.
Inputs: `mpc` - MATPOWER case <br>
        `method` - method used to compute participation factors <br>
                 `1`: participation factors based on the cost function. <br>
                 `2`: participation factors based on the generation capacity of the buses.<br>
Output: `V` (Amplitude of voltage), `success` (1 if convergence achieved and 0 if failed), `et` (elapsed time), and `nither` <br> (number of iterations).<br>

## performance
The code run over many MATPOWER cases to test the overall performance. The x-axis measured elapsed time before the code converged and the y-axis measured the fraction of all cases that converged at a given time.<br>
![performance of NR](https://github.com/zhengboxuan/Distributed-buses-of-power-flow/blob/master/NR.png)<br>
![performance of FDPF](https://github.com/zhengboxuan/Distributed-buses-of-power-flow/blob/master/FDPF.png)

