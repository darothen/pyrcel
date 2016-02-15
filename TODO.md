# Master TODO list

1. Currently uses virtual temperature ($T_v = T(1 + 0.61w)$) when approximating density; should rather use the density temperature for consistency.
    - updates in thermo.py, parcel_aux.pyx, and fortran code 

2. Flesh out basic plotting routines for quick visualization.

3. Add a simple IO package
    - save state for an initial simulation, and its end point
    - read in that state to initialize a new model
    - better flexibility between pandas / xray as basic storage options

4. Re-structure integration logic
    - model should accept two timesteps - one of the numerical integration, one for the output
        + Actually not a problem for the variable-step size solvers; e.g. for CVODE, dt will be the output points
    - the integration should proceed piecewise between output timesteps
        + Already does for Assimulo; 1-minute chunks reporting at the desired timestep
        
5. Add activation diagnostics to integration loop

