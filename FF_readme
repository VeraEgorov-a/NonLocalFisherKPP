Front-Fixing Method for Nonlocal Fisher-KPP Problem

This MATLAB function implements the **Front-Fixing method** for solving a Fisher-KPP equation with nonlocal diffusion and free boundaries. It utilizes a two-step transformation combined with explicit finite difference methods (FDM).

#### Overview:
The method employs:
1. A **front-fixing transformation** to map the original domain with free boundaries onto a fixed domain.
2. Explicit **finite difference schemes** for solving the transformed equations.
3. **Simpson’s rule** for numerical integration of nonlocal terms.

#### Inputs:
- **scheme**: Numerical scheme indicator (0 = Lax, 1 = Forward-Time Backward-Space (FTB-FS), 2 = Forward-Time Centered-Space (FTCS)).
- **M**: Number of spatial grid points.
- **T**: Final simulation time.
- **mu**: Parameter controlling the boundary dynamics.
- **h0**: Initial length of the domain.
- **u0**: Function handle for the initial population distribution \( u_0(x) \).
- **alpha2**: Diffusion coefficient.
- **f**: Function handle for the growth term \( f(u) \).
- **J**: Kernel function for nonlocal diffusion.
- **K**: Cumulative distribution function of the kernel.
- **dt** *(optional)*: Time step size (default is computed based on stability criteria).

#### Outputs:
- **X**: Evolving spatial grid at each time step.
- **W**: Solution matrix, where each row corresponds to the population density at a specific time.
- **H**: Upper boundary values \( h(t) \).
- **G**: Lower boundary values \( g(t) \).
- **t**: Time vector.

#### Core Computations:
1. **Grid Setup**:
   - Spatial step size \( h = \frac{1}{M} \) and time step size \( k \) are determined.
   - Initial conditions for \( x \), \( w \), and domain boundaries are established.

2. **Numerical Integration**:
   - Simpson’s rule is used for evaluating integrals involving \( w \) and kernel functions \( J \) and \( K \).

3. **Finite Difference Schemes**:
   - Three explicit schemes are supported:
     - **Lax-Wendroff** (default for high accuracy).
     - **Forward-Time Backward-Space (FTB-FS)**.
     - **Forward-Time Centered-Space (FTCS)**.
   - The schemes update \( w \) at each time step based on advection and diffusion contributions.

4. **Boundary Updates**:
   - Stefan-like conditions are applied to update \( g(t) \) and \( h(t) \) based on population fluxes.

5. **Stability Considerations**:
   - Courant-Friedrichs-Lewy (CFL) condition is checked to ensure numerical stability.

6. **Performance Logging**:
   - Execution time and step sizes are displayed for monitoring.

#### Example Usage:
```matlab
% Define parameters and functions
M = 100; T = 50; mu = 1.5; h0 = 0.7; alpha2 = 2;
u0 = @(x) max(h0^2 - x.^2, 0); % Initial population density
f = @(u) 2 * u .* (1 - u);     % Logistic growth
J = @(x) 1 ./ (pi * (1 + x.^2)); % Cauchy kernel
K = @(x) atan(x) / pi + 0.5;    % Cumulative of J

% Run the front-fixing method
[X, W, H, G, t] = FF(0, M, T, mu, h0, u0, alpha2, f, J, K);

% Plot the results
plot(t, H, '-r', t, G, '-b');
legend('h(t)', 'g(t)');
xlabel('Time');
ylabel('Boundaries');
```


#### Notes:
- Ensure that the kernel \( J(x) \) satisfies normalization (\( \int J(x) dx = 1 \)) for consistent results.
- Adjust \( M \) and \( dt \) for higher resolution and stability.

