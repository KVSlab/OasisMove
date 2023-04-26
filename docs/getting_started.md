(gs:gs)=
# Getting started

Running simulations in OasisMove is performed using a solver based on an incremental pressure correction scheme, as
presented in {cite}`simo1994unconditional`. The algorithm is based on a fractional step method, and can be summarized in
the following pseudocode:

```python
while t < T:
    t += dt
    solve mesh equation
    update mesh coordinates
    for i in range(max_inner_iters):
        solve tentative velocity
        solve pressure
    update velocity
```

In OasisMove, the solver has been extended for transient flows in moving domains, and the main implemented is located in
the [NSfracStepMove.py](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/NSfracStepMove.py) Python file.
Problems are solved by running the ``NSfracStepMove.py`` script, or by using the `oasism` executable. To run the default
problem `DrivenCavity.py`, you can run the following command:

``` console
$ oasism NSfracStepMove 
```

To specify a specific problem, you can run the following command:

``` console
$ oasism NSfracStepMove problem=DrivenCavity 
```

```{note}
When loading a problem file, OasisMove will start by looking inside the `problems/NSfracStep` folder, and then
look in the current working directory. If the problem you have specified after `problem=` is not located in either locations, an error will occur.
```

```{bibliography}
:filter: docname in docnames
```
