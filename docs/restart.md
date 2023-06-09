(tutorial:restart)=

# Restarting a simulation

In OasisMove, you can restart a simulation that was previously stopped or ended prematurely. The restarted simulation
will load solutions, the mesh, and parameters representing the last state of your simulation, which are stored in
the `Checkpoint` folder. The files will be stored to the `Checkpoint` folder during simulation based on the value of
the `checkpoint` argument. For instance, if `checkpoint=100`, the current velocity and pressure solution, mesh, and
parameters will be stored every *100 time steps*, overwriting the previous files. The default value is `checkpoint=10`,
which can be changed either through the command line as an argument, or hard-coded into the problem parameters in the
problem file.

## Restarting a simulation in OasisMove

To demonstrate restarting a simulation, we will consider
the [`DrivenCavity`](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/DrivenCavity.py)
problem. For this problem, the default values are $T=10$, a time step of $\Delta t=0.005$, and checkpointing every 500
time step (`checkpoint=500`). Assume we have run the problem using the following command:

``` console
$ oasism NSfracStep solver=IPCS_ABCN problem=DrivenCavity 
```

After 500 time steps, or when the simulation is finished, there will be located a `Checkpoint` folder within
the `results_driven_cavity` folder with the checkpointed solutions, mesh, and parameters. We now consider two scenarios
for restarting a simulation.

### Restarting a prematurely ended simulation

If the simulation ended prematurely for any reason, and given that the simulation lasted for more than 500 time steps,
we can restart the simulation at the latest checkpoint. Assuming we have fixed whatever caused the solution to end
prematurely, we may restart the simulation by adding the  `restart_folder` parameter to the run command, which points to
the `Checkpoint` folder:

``` console
$ oasism NSfracStep solver=IPCS_ABCN problem=DrivenCavity restart_folder=results_driven_cavity/data/1/Checkpoint 
```

The simulation should now restart at the latest checkpoint, and continue as normal until $T=10.

### Restarting a finished simulation

If the simulation was successful, but we want to simulate the problem for longer, we can restart the simulation at
the latest checkpoint, often representing the final solution. In this case, we will have to supply the `restart_folder` parameter to the run command, pointing
to the `Checkpoint` folder, and a new end-time ($T$) for the simulation. With a default value of $T=10$, we can restart
the `DrivenCavity` problem and continue running it until $T=20$ with the following command:

``` console
$ oasism NSfracStep solver=IPCS_ABCN problem=DrivenCavity T=20 restart_folder=results_driven_cavity/data/1/Checkpoint  
```

The simulation should now restart at $T=10$, and continue running until $T=20$.

