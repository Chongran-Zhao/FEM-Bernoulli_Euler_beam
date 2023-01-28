# Intro
A CPP code developed to solve 1D Bernoulli Euler beam problem in FEM using continuous Hermite-shape function.

### Phsical problem
$u_{,xxxx} = l(x)$ on $\mathbf{\Omega} = (0,L)$ (transverse equilibrium)

$u(L) = 0$ (zero transverse displacement)

$u_{,x}(L) = 0$ (zero slope)

$EI u_{,xx}(0) = M$ (prescribed moment)

$EI u_{,xxx}(0) = Q$ (prescribed shear)

Notes: $E$ is Young's modulus and $I$ is the moment of inertia, both of which are assumed to be constant.

### Exact solution
We can solve exact solution of this 4th order ODE with four boundary conditions

### Code details
`void Gauss(int N, double a, double * qp, double * wq)` assigns values to `qp` and `wq`.

`double Hermitebasis(const double &x1, const double &x2, int i, int der, double &x)` returns specific value of Hermite-shape function.

Hermite-shape function has $C^1$ continuity on the element nodes.

`uu = solver.solve(F)` is filled with $[u_{1},u_{1,x},u_{2},u_{2,x},...]^T$.

### Run the code
1. Run `mkdir build` under `FEM-Bernoulli_Euler_beam/` to create a new directory.
2. Run `cd build` to step in the new directory.
3. Run `cmake ..` to create `makefile` using `CMakeLists.txt` from upper directory.
4. Run `make` to create an executable file named `Bernoulli_Euler_beam`.
5. Run `./Bernoulli_Euler_beam`.
