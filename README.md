# One-dimensional colliding turbidity currents

## Simulation
In `solver1d.c` we solve the shallow-water equations using a finite-volume scheme with a Harten-Lax-Leer (HLL) flux. A fifth-order weighted essentially non-oscillatory (WENO) method is used to construct the spatial stencil for the finite-volume scheme and a third-order (3-stage) explicit strong stability preserving Runge-Kutta (SSPERK) scheme is used for time integration. The data will be stored in plain text files to be post-processed in python.

A C compiler is needed to compile the program. One example to compile, using `gcc`, is as follows:
```bash
gcc -O3 -o sim.out solver1d.c
```
If there is any error compiling, first try removing the -O3 flag. 

After compiling, run the executable as follows: 
```bash
./sim.out 1.0 1.0 0.02
```

The code is for scientific research and intended to run different test cases with varying initial height, initial concentration, and settling speeds. Therefore, the exectuable takes three arguments `h2init`, `c2init`, `U_s`, respectively. To lower the initial height and concentration of the current on the right to 0.7 and 0.9 and reduce the settling speed to 0.01, run the executable as follows: 
```bash
./sim.out 0.7 0.9 0.01
```

## Post-processing
In `PostProc.py` we post-process the data for visualization and scientific research.

### Installation
Install the required packages using pip:

```bash
pip install numpy matplotlib celluloid 
```

### Usage

Run the post-process code as follows in Python's interactive mode. 
```bash
python -i PostProc.py
```

Then, the `LoadSim` class can be called to load the specific simulation. 
```python
sim=LoadSim(1.0,1.0,0.02)
```
To visualize the data either plot a specific variable using the `plot_times` method, the first argument is the variable (options are `'h'`,`'u'`,`'c1'`,`'c2'`) and the second is a list of times, or generate an MP4 with the `makeMP4` method. 
```python
sim.makeMP4()
sim.plot_times('h',[0,2,4,6,8])
```

To see all the relevant information of the simulation, run the `sim_info` method. 
```python
sim.sim_info()
```
## Customizing the simulation.

The file `solver1d.c` is modular and can be adapted to run different simulations of a system of 4 conservation laws. Here are the main places to edit the code for a different simulation.  
- **Parameters:** At the beginning of the file there are several parameters that can be changed, `N` (number of cells), `T` (final time), etc...
- **Initial conditions:** Initial conditions are defined in the `initialize` function.
- **Boundary conditions:** How the stencil behaves at the boundary is defined in the `BC` function. 
- **Flux:** For different conservation laws the flux  
- **Size of system:** This is written for a system with 4-depenednet variables. For systems with additional variables a new variable will need to be defined and added throughout the entire code. 

