# Examples

For all simulations:

```text
Legend:
✅ : converged
👍 : converging
🥵 : diverging
😱 : found NaN or Inf
❋  : non-scaled max(R)
?  : no info abut convergence
```

## Heat: Arpaci Nonlinear 1d

Arpaci's Example 3-8 on page 130 (variable conductivity)

* Arpaci V. S. (1966) Conduction Heat Transfer, Addison-Wesley, 551p

### Test goal

This tests verifies the nonlinear solver for the diffusion equation
with a variable conductivity coefficient.

### Mesh

```text
o-----------------------------------------------------------o
|    |    |    |    |    |    |    |    |    |    .....     | h = 1
o-----------------------------------------------------------o
                     <-  L = 10 ->
```

### Initial conditions

Temperature T = 0 at all points

### Boundary conditions

Temperature T = 0 on right side @ x = L

### Configuration and parameters

* Steady simulation
* Source = 5
* Variable conductivity (k = (1 + β T) kᵣ I) with kᵣ = 2

### Note

The temperature at the right T = 0 (T_inf) must be zero in order to
result in k(T_inf) = kᵣ as required by the analytical solution.

### Simulation file

[heat_arpaci_nonlinear_1d.rs]

### Results

```text
                                                  _   
timestep             t            Δt  iter    max(R)  
       1   1.000000e-1   1.000000e-1     .        .  
       .             .             .     1    2.50e0❋ 
       .             .             .     2    1.28e0? 
       .             .             .     3   9.36e-2👍
       .             .             .     4   1.25e-3👍
       .             .             .     5   1.42e-7👍
       .             .             .     6  1.35e-14✅

T(0) = 87.08286933869708  (87.08286933869707)
```

## Heat: Bhatti 1d5 Convection

## Heat: Bhatti 6d22 Convection

## Heat: Lewis Transient 1d

## Heat: Mathematica Axisymmetric Nafems

## Heat: Mathematica Axisymmetric Simple

## Rod: Bhatti 1d4 Truss

## Solid Bhatti 1d6 Plane Stress

## Solid Felippa Thick Cylinder Axisymmetric

## Solid Smith 5d11 Qua4 Plane Strain_uy

## Solid Smith 5d15 Qua8 Plane Strain

## Solid Smith 5d17 Qua4 Axisymmetric

## Solid Smith 5d24 Hex20 3D

## Solid Smith 5d27 Qua9 Plane Strain

## Solid Smith 5d2 Tri3 Plane Strain

## Solid Smith 5d30 Tet4 3D

## Solid Smith 5d7 Tri15 Plane Strain
