# Parameters 

Every problem will contain parameters which area both specific to the problem, but also general to Oasis.
These can be altered by adding the as command line arguments. 
As an exaple, we will consider the simplest rigid flow problem, the default `DrivenCavity` problem. 

The parameters are 
in 



```python
NS_parameters.update(
    # Fluid parameters
    nu=0.001,
    # Simulation parameters
    T=1.0,
    dt=0.005,
    folder="results_driven_cavity",
    # Oasis parameters
    testing=False,
    max_iter=2,
    dynamic_mesh=False,
    save_solution_frequency=1,
    checkpoint=500,
    print_intermediate_info=100,
    velocity_degree=1,
    pressure_degree=1,
    use_krylov_solvers=True,
    max_error=1e-8)
```


``` console
$ oasism NSfracStepMove 
```



```{bibliography}
:filter: docname in docnames
```
