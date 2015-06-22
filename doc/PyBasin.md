PyBasin
=======


## Decompaction

Porosity is calculated assuming an exponential decay of porosity with depth [@Athy1930;Fowler2008]:

\begin{equation}
	n = n_0 e^{-c z}
\end{equation}


The thickness of the pore space $b_w$ is related to porosity $n$ by [@Allen2005]:

\begin{equation}
	b_w = \int _{z1}^{z2} n_0 e^{-c z}
\end{equation}

which after integration gives:

\begin{equation}
	b_w = \frac{n_0}{c} \left(e^{-cz_1} - e^{-cz_2} \right)
\end{equation}


The total thickness is equal to the thickness of the pore space plus the thickness of the matrix: $b_t = b_w + b_m$

The initial, decompacted thickness can be calculated by setting $z_1 = 0$ and $z_2 = b_0$, where $b_0$ is the initial thickness:

\begin{equation}
	b_0 = b_m + \frac{n_0}{c} \left(1 - e^{-cb_0} \right)
\end{equation}

This equation cannot be solved analytically. We use iteration to find the value of $b_0$:

\begin{equation}
	b_0^{i+1} = b_m + \frac{n_0}{c} \left(1 - e^{-cb_0^{i}} \right)
\end{equation}

where $^i$ denotes a step in the iteration. In general 10 to 20 iterations are required to find a value within 0.01 m of the 'real' decompacted thickness.


## Decompaction with exhumation phases

Things get significantly more complicated if we included exhumation phases in the simulated burial history. 

One solution is to first calculate the decompacted thickness assuming maximum present-day burial. Then go back to each of the exhumation phases and calculate the estimated burial depth of each unit. If the burial depth exceeds present-day burial depths then the thickness of the pore space $b_w$ and the matrix thickness $b_m$ need to be recalculated. To keep things simple we assume no thickness change during exhumation, only during burial.

In this case the pore thickness can be calculated as:

\begin{equation}
	b_w = \frac{n_0}{c} \left(e^{-cz_{1,m}} - e^{-cz_{2,m}} \right)
\end{equation}

where $z_{1,m}$ and $z_{1,m}$ are the top and the base of the unit at maximum burial depth. Because we assume no thickness changes during uplift, the thickness at maximum burial is the same as the present-day thickness ($b_p$): 

\begin{equation}
	z_{2,m} = z_{1,m} + b_p
\end{equation}

However, we do not know in advance if the burial depth of a unit was deeper than the present-day burial or not.

The workflow is to calculate the initial and matrix thickness first on the basis of the present-day burial depth. Subsequently we go back to the exhumation phase and add the eroded section on top. Then we stack all the other units below the eroded unit and calculate their porosity and thickness. ...

For the first unit the thickness is simply the eroded thickness:

For the subsequent units, denoted by $_j$:

\begin{equation}
	b_w,j = \frac{n_0}{c} \left(e^{-cz_{b,j-1}} - e^{-cz_{b,j}} \right)
\end{equation}

and the total thickness is again the sum of the matrix and the pore water thickness:

\begin{equation}
    b_j = b_w,j + b_m,j
\end{equation}

and the base depth of unit $j$ is equal to the top plus the thickness:

\begin{equation}
    z_j = z_{j-1} + b_j
\end{equation}

again this leaves us with an equation that cannot be solved analytically:

\begin{equation}
	b_j = b_m + \frac{n_0}{c} \left(e^{-cz_{b,j-1}} - e^{-c (z_{b,j-1}+b_j)} \right)
\end{equation}



# Heat flow

Heat flow is governed by the following equation

\begin{equation}
    \rho c \frac{\partial T}{\partial t} = \nabla q + Q
\end{equation}
where $\rho$ is density (kg m^-3^), *c* is heat capacity (), *T* is temperature (K), *t* is time (sec), *q* is heat flux () and Q is heat production (W m^-3^).

We discretize the drivative of temperature in time as 
\begin{equation}
    \frac{\partial T}{\partial t} = \frac{T^{t+\Delta t} - T^t}{\Delta t}
\end{equation}

heat flow is defined as
\begin{equation}
    q = K \nabla T
\end{equation}
where *K* is thermal conductivity (W m^-1^ K^-1^).

For 1D domains, the heat flow equation can be discretized spatially as
\begin{equation}
    \nabla q_i = \frac{q_{i+1/2} - q_{i-1/2}}{\Delta x_i}
\end{equation}

\begin{equation}
    q_{i+1/2} = K_{i+1/2} \frac{T_{i+1} - T_i}{\Delta x_{i+1/2}}
\end{equation}

\begin{equation}
    q_{i-1/2} = K_{i-1/2} \frac{T_{i} - T_{i-1}}{\Delta x_{i-1/2}}
\end{equation}

And writing everything using x coordinates:

\begin{equation} 
    \Delta x_{i+1/2} = x_{i+1} - x_i
\end{equation}

\begin{equation}
    q_{i+1/2} = K_{i+1/2} \frac{T_{i+1} - T_i}{x_{i+1} - x_i}
\end{equation}

\begin{equation}
    q_{i-1/2} = K_{i-1/2} \frac{T_{i} - T_{i-1}}{x_i - x_{i-1}}
\end{equation}

therefore

\begin{equation}
    \nabla q_i =  \frac{2}{x_{i+1} - x_{i-1}} \left( K_{i+1/2} \frac{T_{i+1} - T_i}{x_{i+1} - x_i} - K_{i-1/2} \frac{T_{i} - T_{i-1}}{x_i - x_{i-1}} \right)
\end{equation}

## Boundaries

For the left-hand boundary, with a fixed flux $q_l$ from the left:

\begin{equation}
    \nabla q_0 = \dfrac{q_{1/2} - q_l}{\Delta x_0}
\end{equation}

\begin{equation}
    \nabla q_0 = \dfrac{q_{1/2} - q_l}{\dfrac{x_1 - x_0}{2}}
\end{equation}

\begin{equation}
    \nabla q_0 = \dfrac{2}{x_1 - x_0} \left( q_{1/2} - q_l \right)
\end{equation}

\begin{equation}
    \nabla q_0 = \dfrac{2}{x_1 - x_0} \left( \left(K_{1/2} \frac{T_1 - T_0}{x_1 - x_0} \right) - q_l \right)
\end{equation}

For the right-hand boundary at the last node $n$, with a fixed flux $q_r$ from the right:

\begin{equation}
    \nabla q_n = \dfrac{q_r - q_{n-1/2}}{\Delta x_0}
\end{equation}

\begin{equation}
    \nabla q_n = \dfrac{2}{x_n - x_{n-1}} \left(q_r - q_{n-1/2} \right)
\end{equation}

\begin{equation}
    q_{n-1/2} = K_{n-1/2} \frac{T_{n} - T_{n-1}}{x_n - x_{n-1}}
\end{equation}

\begin{equation}
    \nabla q_n = \dfrac{2}{x_n - x_{n-1}} \left(q_r - \left( K_{n-1/2} \frac{T_{n} - T_{n-1}}{x_n - x_{n-1}} \right) \right)
\end{equation}


## Constant $K$ and $\Delta x$

To keep things simple we first try to solve this for a simple case where thermal conductivity ($K$) and grid size ($\Delta x$) are constant.

in the simple case with constant $\Delta x$ and constant $K$ we get:

\begin{equation}
    \nabla q =  \frac{2}{\Delta x} \left( K \frac{T_{i+1} - T_i}{\Delta x} - K \frac{T_{i} - T_{i-1}}{\Delta x} \right)
\end{equation}

\begin{equation}
    \nabla q =  \frac{2K}{\Delta x ^2} \left( T_{i+1} + T_{i-1} - 2 T_i \right)
\end{equation}

and for the full heat flow equation

\begin{equation}
    \rho c \frac{T_i^{t+1} - T_i^t}{\Delta T} = \frac{2K}{\Delta x ^2} \left( T_{i+1} + T_{i-1} - 2 T_i \right) + Q
\end{equation}

implicit FD:
\begin{equation}
    \rho c \frac{T_i^{t+1} - T_i^t}{\Delta T} = \frac{2K}{\Delta x ^2} \left( T_{i+1}^{t+1} + T_{i-1}^{t+1} - 2 T_i^{t+1} \right) + Q
\end{equation}
where T^t^ is known and T^t+1^ is unknown

reorganizing:
\begin{equation}
    T_i^{t+1} - T_i^t = s \left( T_{i+1}^{t+1} + T_{i-1}^{t+1} - 2 T_i^{t+1} \right) + \frac{Q \Delta t}{\rho c}
\end{equation}

where
\begin{equation}
    s = \frac{2K \Delta t}{\rho c \Delta x^2}
\end{equation}

and putting all unknowns ($T_{...}^{t+1}$) on the left hand side
\begin{equation}
    s \left( T_{i+1}^{t+1} + T_{i-1}^{t+1} - 2 T_i^{t+1} \right) - T_i^{t+1}  =  - \frac{Q \Delta t}{\rho c} - T_i^t
\end{equation}

\begin{equation}
    -s T_{i-1}^{t+1} + (1+ 2s) T_i^{t+1} - s T_{i+1}^{t+1}  =  T_i^t + \frac{Q \Delta t}{\rho c} 
\label{eq:HF_FD_simple}
\end{equation}


## Matrices

We now have discretized the heat flow equation for a uniform grid. We have one equation for each grid cell, and one unknown for each grid cell: the temperature at the next timestep $T_i^{t+1}$. Recall from your math class that if you have an equal number of unkowns and equations one can a set of equations to find the unknowns. 

The most efficient way to solve a system of equations is to use matrices. First we convert equation \ref{eq:HF_FD_simple} to matrix notation:

\begin{equation}
    A \vec{x} = \vec{b}
	\label{eq:matrix}
\end{equation}

where $A$ is a matrix and $\vec{x}$ and $\vec{v}$ are vectors. We can fit equation \ref{eq:HF_FD_simple} in this equation by moving all the coefficients to matrix $A$, the unknowns ($T_0^{t+1}$, $T_1^{t+1}$, $T_2^{t+1}$, etc..) to vector $\vec{x}$ and the known variables ($T_i^t + \frac{Q \Delta t}{\rho c} $) to vector $\vec{b}$.

For an example grid containing 6 cells, the matrix $A$ takes the shape of:

\begin{equation}
    A = 
    \begin{bmatrix} 
        1  & 0 & 0  & 0 & 0 & 0 \\
        -s & (1+2s) & -s & 0 & 0 & 0 \\
        0  & -s & (1+2s) & -s & 0 & 0 \\   
        0  & 0 & -s & (1+2s) & -s & 0 \\   
        0  & 0  & 0 & -s & (1+2s) & -s \\   
        0  & 0 & 0  & 0 & -s & (1+s) \\   
    \end{bmatrix}
\end{equation}

Note that the top and bottom line are different from what you would expect from eq. \ref{eq:HF_FD_simple}. The reason for this is that we have to account for boundary conditions, more on this in the next section. 

The vector $\vec{b}$ is:
\begin{equation}
    \vec{b} = 
    \begin{bmatrix} 
        T_s \\
        T_1^t + \frac{Q \Delta t}{\rho c} \\
        T_2^t + \frac{Q \Delta t}{\rho c} \\
        T_3^t + \frac{Q \Delta t}{\rho c} \\
        T_4^t + \frac{Q \Delta t}{\rho c} \\
        T_5^t + \frac{(Q +w/\Delta x) \Delta t}{\rho c} \\
    \end{bmatrix}
\end{equation}

and the vector containing the unknowns is:

\begin{equation}
    \vec{x} = 
    \begin{bmatrix} 
        T_0^{t+1} \\
        T_1^{t+1} \\
        T_2^{t+1} \\
        T_3^{t+1} \\
        T_4^{t+1} \\
        T_5^{t+1} \\
    \end{bmatrix}
\end{equation}


## Solving system of equations

One way to solve a system of equations is to calculate the inverse matrix. If we multiply the left and right hand side of equation \ref{eq:matrix} with the inverse matrix we get:

\begin{equation}
    A^{-1} A x = A^{-1} b
	\label{eq:inverse}
\end{equation}

multiplying a matrix with its inverse results in the identity matrix $A^{-1} A = I$:

\begin{equation}
I = 
\begin{bmatrix} 
    1  & 0 & 0  & 0 & 0 & 0 \\
    0  & 1 & 0  & 0 & 0 & 0 \\
    0  & 0 & 1  & 0 & 0 & 0 \\
    0  & 0 & 0  & 1 & 0 & 0 \\
    0  & 0 & 0  & 0 & 1 & 0 \\
    0  & 0 & 0  & 0 & 0 & 1 \\
\end{bmatrix}
\end{equation}

and therefore we can simplify equation \ref{eq:inverse}:

\begin{equation}
    x = A^{-1} b
\end{equation}

more on how to find the inverse matrix.... with python ...


## Checking solutions with analytical equations

We check whether our numerical model is correct by comparing our solution to an analytical solution for the cooling of an intrusive in the subsurface. The solution for temperature change of an initially perturbed temperature field is [@Carslaw1959]:

\begin{equation}
    T(z,t) = T_b + \frac{T_i - T_b}{2} \left( erf \left(\frac{L-z}{2 \sqrt{\kappa t}} \right) + erf \left(\frac{L+z}{2 \sqrt{\kappa t}} \right) \right)
\end{equation}

where *T~b~* is the background temperature, *T~i~* is the temperature of the intrusive, *L* is the length of the intrusive.

Following [@Ehlers2005] we use this equation to simulate cooling of an intrusive in the subsurface. The numerical and analytical solutions for cooling are shown in Figure 1. The Intrusive body has an initial temperature of 700 degrees C, and stretches from 0 to 500 m distance. The background temperature is 50 degrees C. The solutions match to within ... degree C and show the gradual decrease in temperatures. 

![Comparison of numerical and analytical solution for the cooling of an intrusive in the subsurface. The analytical and numerical solution overlap, the maximum difference is 1 degree C. $\Delta$x = 5 m, $\Delta$t = 10.0 years, $\kappa$ = 1.01 $\times$ 10^-6^ m^2^ s^-1^ (=32 km^2^ Ma^-1^).](fig/intrusive_test_heat_flow.pdf)


## Variable $K$ and $\Delta x$

Now we have successfully solved heat flow for a simple case with constant thermal conductivity, lets make things more complex and realistic by introducing non-constant grid spacing and thermal conductivity.

Recall that the discretized equation for $\nabla q$ is:

\begin{equation}
	\nabla q =  \frac{2}{x_{i+1} - x_{i-1}} \left( K_{i+1/2} \frac{T_{i+1} - T_i}{x_{i+1} - x_i} - K_{i-1/2} \frac{T_{i} - T_{i-1}}{x_i - x_{i-1}} \right)
\end{equation}

if we define 
\begin{equation}
	s = \frac{2}{x_{i+1} - x_{i-1}}
\end{equation}

\begin{equation}
	t = \frac{K_{i+1/2}}{x_{i+1} - x_i}
\end{equation}

\begin{equation}
	u = \frac{K_{i-1/2}}{x_{i} - x_{i-1}}
\end{equation}

we end up with a substituted equation:

\begin{equation}
	\nabla q = su T_{i-1}^{t+1} - (su + st) T_i^{t+1} + st T_{i+1}^{t+1}
\end{equation}

adding this to the full heat flow equation yields:

\begin{equation}
	\rho_i c_i \frac{T_i^{t+1} - T_i^t}{\Delta t} = su T_{i-1}^{t+1} - (su + st) T_i^{t+1} + st T_{i+1}^{t+1} + Q
\end{equation}

if we define $v$ as:
\begin{equation}
	v = \frac{\Delta t}{\rho_i c_i}
\end{equation}

we can rewrite this as

\begin{equation}
	T_i^{t+1} - T_i^t = suv T_{i-1}^{t+1} - (suv + stv) T_i^{t+1} + stv T_{i+1}^{t+1} + Qv
\end{equation}

which after rearranging gives:

\begin{equation}
	T_i^t + Q v = -suv T_{i-1}^{t+1} + (1 + (suv + stv)) T_i^{t+1} - stv T_{i+1}^{t+1}
\end{equation}

### Lower fixed flux boundary:

For the lower node $n$ with a specified flux $q_r$ the flux is given by:

\begin{equation}
    \nabla q_n = \dfrac{2}{x_n - x_{n-1}} \left(q_r - \left( K_{n-1/2} \frac{T_{n} - T_{n-1}}{x_n - x_{n-1}} \right) \right)
\end{equation}

and therefore the full equation for the last node is:

\begin{equation}
    \rho_n c_n \frac{T_n^{t+1} - T_n^t}{\Delta t} = \dfrac{2}{x_n - x_{n-1}} \left(q_r - \left( K_{n-1/2} \frac{T_{n} - T_{n-1}}{x_n - x_{n-1}} \right) \right) + Q
\end{equation}

if we define:

\begin{equation}
	s_n = \dfrac{2}{x_n - x_{n-1}}
\end{equation}

\begin{equation}
	u_n =  \dfrac{K_{n-1/2} }{x_n - x_{n-1}}
\end{equation}

\begin{equation}
    T_n^{t+1} - T_n^t = s_n v \left(q_r - \left( u_n (T_{n} - T_{n-1}) \right) \right) + Q v
\end{equation}

\begin{equation}
    T_n^{t+1} - T_n^t = s_n v q_r - \left( s_n v u_n (T_{n} - T_{n-1}) \right) + Q v
\end{equation}

\begin{equation}
    T_n^{t+1} - T_n^t = s_n v q_r - s_n v u_n T_{n}^{t+1} + s_n v u_n T_{n-1}^{t+1} + Q v
\end{equation}

\begin{equation}
     - T_n^t - Q_v - s_n v q_r =  -(1+ s_n v u_n T_{n}^{t+1}) + s_n v u_n T_{n-1}^{t+1}
\end{equation}

\begin{equation}
     T_n^t + Q_v + s_n v q_r =  - s_n u_n v T_{n-1}^{t+1} + (1+ s_n u_n  v T_{n}^{t+1}) 
\end{equation}

and therefore the last term in the matrix is 

### Final matrix


casting this equation into matrix form we end up with a slightly more complex matrix $A$:

\begin{equation}
    A = 
    \begin{bmatrix} 
        1  & 0 & 0  & 0 & 0 & 0 \\
        -suv & (1+stv+suv) & -stv & 0 & 0 & 0 \\
        0  & -suv & (1+stv+suv) & -stv & 0 & 0 \\   
        0  & 0 & -suv & (1+stv+suv) & -stv & 0 \\   
        0  & 0  & 0 & -suv & (1+stv+suv) & -stv \\   
        0  & 0 & 0  & 0 & -suv & (1+suv) \\   
    \end{bmatrix}
\end{equation}

The vector $\vec{b}$ is:
\begin{equation}
    \vec{b} = 
    \begin{bmatrix} 
        T_s \\
        T_1^t + Q v \\
        T_2^t + Q v \\
        T_3^t + Q v \\
        T_4^t + Q v \\
        ... \\
			  T_n^t + Q v + s v q_r\\
    \end{bmatrix}
\end{equation}

\pagebreak

# References