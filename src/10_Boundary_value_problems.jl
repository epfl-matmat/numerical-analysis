### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 2cff355a-0339-11ef-085c-a5f9a6f5549e
begin
	using Plots
	using PlutoUI
	using PlutoTeachingTools
	using LaTeXStrings
	using LinearAlgebra
	using HypertextLiteral
end

# ╔═╡ b72b45ad-6191-40cb-9e9f-950bf1bfe212
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/10_Boundary_value_problems.pdf)
"""

# ╔═╡ 206ae56c-fcfa-4d6f-93e4-30f03dee8f90
TableOfContents()

# ╔═╡ 8c818dad-6789-4771-9eab-4e4ba78f5612
md"""
# Boundary value problems

The goal of boundary value problems is to reconstruct the function values $u(x, t)$
of a function $u$ on an entire domain $x \in [a, b]$ 
and for all times $t$ knowing only two things
- Some information about $u$ at the boundary $x=a$ and $x=b$
  for the initial time $t=0$. For example we may know
  the values $u(a, 0)$ and $u(b, 0)$.
- The relationship between $u$ and its derivative at each point in space-time $(x, t)$.

Classic examples of boundary value problems include
- the [heat equation](https://en.wikipedia.org/wiki/Heat_equation),
  where $u(x, t)$ is the temperature at time $t$ and position $x$
  of a material --- and the flux of heat at $(x, t)$ itself depends on the temperature at this point in space and time, or
- [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation)
  where $u(x, t)$ is the electrostatic potential at position $x$ and time $t$
  corresponding to a given time-dependent density of electronic charges --- where the charge itself determines the electrostatic potential acting at each point $(x, t)$, or
- the [diffusion equation](https://en.wikipedia.org/wiki/Diffusion_equation)
  where $u(x, t)$ is the concentration of some solute at time $t$
  at position $x$ in the solution.

We will first consider the heat equation as one example.
"""

# ╔═╡ cf1055db-a773-4dc2-bd66-0ae840c656a4
md"""
## Stationary heat equation in one dimension
"""

# ╔═╡ 93b28014-c767-4f40-b27c-4646f63faa89
RobustLocalResource("https://raw.githubusercontent.com/epfl-matmat/numerical-analysis/ab92c724d96689f4ebbc5059ca9bac3baa6ac61b/src/img/heat-equation.svg", "img/heat-equation.svg", :width => 500)

# ╔═╡ f1b662de-b4cf-4a73-9c00-580a2ea8a428
md"""
We consider a metal bar extending over the interval $[a, b]$. For a given moment in time $t$ we denote the **temperature** at a point $x \in [a, b]$ by $u(x, t)$. The **heat flux** we fruther denote by $\phi(x, t)$ and allow for an **external source of heat** $f(x)$, e.g. by one more more heat guns pointed at the metal bar.
"""

# ╔═╡ b5c3ed4b-3c4d-40d5-9171-94187416b363
md"""
The **rate of change of the internal heat**
$\frac{\partial Q}{\partial t}$
at point $x$ of the rod
is proportional to the change in temperature $\frac{\partial u}{\partial t}$,
that is
```math
\frac{\partial Q}{\partial t} = c ρ \, \frac{\partial u}{\partial t}(x, t)
```
where proportionality constants are the heat capacity $c$ and the mass density $ρ$.

Due to **conservation of energy** the change of internal heat of a slice $[x, x+δx]$ is equal to the incoming heat $f(x)\,δx + \phi(x, t)$ minus the outgoing heat $\phi(x + δx, t)$, such that we obtain
```math
\frac{\partial Q}{\partial t} δx = \phi(x, t) + f(x)\, δx - \phi(x + δx, t).
```
Dividing by $δx$ and **taking the limit** $δx\to 0$ thus yields
```math
\frac{\partial Q}{\partial t} = - \frac {\partial \phi}{\partial x}(x, t) + f(x).
```
"""

# ╔═╡ 226a30be-f721-4f3d-8d49-abddca22f661
md"""
Combining with **Fourier's law**
```math
\phi(x, t) = -k \frac{\partial u}{\partial x} (x, t)
```
where $k$ is the thermal conductivity, we obtain the equation
```math
c ρ \, \frac{\partial u}{\partial t}(x, t) = k \frac{\partial^2 u}{\partial x^2}(x, t) + f(x) \qquad x \in (a, b), \ t > 0
```
which is the **heat equation** in one dimension,
a **partial differential equation**.
The heat equation describes the **evolution of the temperature in a system**
(here the metal bar) and has as its unknown the temperature distribution $u(x, t)$ at any location at any point in time.
"""

# ╔═╡ 59f11469-c86e-4dc4-9eed-039101aaa058
md"""
To simplify matters in this lecture 
**we restrict** ourselves to the case of the **stationary boundary value problems**, i.e. where all time derivatives vanish. In the case of the heat equation we would thus assume the **temperature $u(x, t)$** to be in a state where it **no longer changes with time**, i.e. $\frac{\partial u}{\partial t} = 0$.
Dropping thus the zero terms and the dependency on $t$ we obtain the **stationary heat equation**
```math
-k\, u''(x) = f(x) \qquad x \in (a, b).
```
where we used the common notation $u'' = \frac{\partial^2 u}{\partial x^2}$ for spatial derivatives.
"""

# ╔═╡ 3e10cf8e-d5aa-4b3e-a7be-12ccdc2f3cf7
md"""
To fully specify the problem and find a unique solution we still **need to indicate what happens at the extremities of the bar**, i.e. at $u(a)$ and at $u(b)$.
Since these two points sit at the boundary of our computational domain $[a, b]$ we usually call these conditions on $u$ **boundary conditions**. 

One simple type of boundary conditions is to just set these extremal points to some values, i.e. enforce $u(a) = γ_a$ and $u(b) = γ_b$. Physically this implies a condition where the bar is in contact with heat baths of constant temperature, such that the temperature values at the boundary are known. We obtain the **heat equation using Dirichlet boundary conditions**:
```math
\left\{
\begin{aligned}
-k \, u''(x) &=  f(x) \qquad x \in (a, b).\\
u(a) &= γ_a\\
u(b) &= γ_b
\end{aligned}
\right.
```

Note that other types of boundary conditions are possible, see for example [chapter 10.1](https://fncbook.com/tpbvp/) of Driscoll, Brown: Fundamentals of Numerical Computation for a more general treatment.

For completeness we now provide the general definition of a boundary value problem (BVP) before specialising to the case we will focus on in this course.
"""

# ╔═╡ 23f9e709-1f77-4a96-bb09-2d358cdecc76
md"""
!!! info "Definition: Boundary value problems (general definition)"
	Given functions $ϕ$, $g_{a}$ and $g_{b}$ as well as a problem domain $[a, b]$
	a boundary-value problem is the set of equations
	```math
	\left\{
	\begin{aligned}
	u''(x) &=  ϕ(x, u, u'') \qquad x \in (a, b),\\
	g_{a}(u(a), u'(a)) &= 0, \\
	g_{b}(u(b), u'(b)) &= 0
	\end{aligned}
	\right.
	```
	to be solved for the unknown function $u : [a, b] → \mathbb{R}$.

	Notice that $ϕ$ thus takes as input a point $x$, the *entire* function $u$ and the *entire* function $u'$.

	The functions $g_{a}$ and $g_{b}$ take two values, namely the function value and the derivative value at the boundary points $a$, respectively $b$. These functions are called **boundary conditions**.

To simplify matters for this course we will make two restrictions to this setting:

- We will only consider **linear boundary value problems**, that is problems where the the functions $ϕ$, $g_{a}$ and $g_{b}$ are linear in the arguments involving $u$ and $u'$. Specifically this means that $\phi(x,u,u') = -p(x)\,u'(x) - q(x)\,u(x) + r(x)$, where $p$, $q$ and $r$ are known functions from $[a, b]$ to $\mathbb{R}$.

- We will only consider **Dirichlet boundary conditions**, that is cases where $g_{a}(u(a), u'(a)) = u(a) - γ_a$ and $g_{b}(u(b), u'(b)) = u(b) - γ_b$ where $γ_a, γ_b \in \mathbb{R}$.
"""

# ╔═╡ 6df8394f-6457-4db0-86e5-7a9dd7d65c9a
md"""
!!! info "Definition: Linear Dirichlet boundary value problem (our setting)"
	Given a domain $[a, b]$ and functions
	$p : [a, b] \to \mathbb{R}$, 
	$q : [a, b] \to \mathbb{R}$, 
	$r : [a, b] \to \mathbb{R}$
	as well as boundary values $γ_a, γ_b \in \mathbb{R}$
	find the function $u : [a, b] \to \mathbb{R}$ such that
	```math
	\tag{1}
	\left\{
	\begin{aligned}
	u''(x) + p(x) \, u'(x) + q(x) \, u(x) &= r(x) \qquad x \in (a, b),\\
	u(a) &= γ_a \\
	u(b) &= γ_b
	\end{aligned}
	\right.
	```
"""

# ╔═╡ b4b33be9-bb9a-4fe9-af48-8883c2fbc7e9
md"""
In line with the discussion at the beginning of this notebook, we recover the 
stationary heat equation from equation (1) by setting $p(x) = 0$, $q(x) = 0$ and $r(x) = -f(x) / k$, resulting in the problem
```math
\tag{2}
\left\{
\begin{aligned}
u''(x) &= -\frac{1}{k} f(x) \qquad x \in (a, b).\\
u(a) &= γ_a, \, u(b) = γ_b
\end{aligned}
\right.
```
"""

# ╔═╡ b00a37d1-37ad-41b5-ba88-d6582ef680cb
# TODO("Sketch also how the Poisson equation and the diffusion equation arise.")

# ╔═╡ be847fe0-7447-493b-9087-f4a1d961edba
md"""
## Differentiation matrices

Having introduced the problem setting of boundary value problems (1), we will now build a numerical scheme for solving such problems in practice.

Similar to our discussion in the previous chapters on 
[Interpolation](https://teaching.matmat.org/numerical-analysis/07_Interpolation.html),
[Numerical integration](https://teaching.matmat.org/numerical-analysis/08_Numerical_integration.html)
or [Numerical differentiation](https://teaching.matmat.org/numerical-analysis/09_Numerical_differentiation.html) 
we first **perform a discretisation**: Instead of solving the problem on all of the continous interval $[a, b]$ we discretise the interval into $N$ subintervals. Taking $N$ pieces of equal length $h = \frac{b-a}{N}$ leads to the nodes
```math
	\tag{3}
	x_i = a + i \, h, \qquad \forall i = 0, \ldots, N, \qquad h = \frac{b-a}{N}.
```
Note, that similar to previous chapters our node indexing starts at $0$.

Since the functions $p$, $q$ and $r$ are given to us, we can easily compute discretised representations
```math
p_i = p(x_i), \quad q_i = q(x_i), \quad r_i = r(x_i), \qquad \forall i = 0, \ldots, N,
```
which we can denote as the vectors, for example
```math
\mathbf{r} = \begin{pmatrix} r_0 \\ r_1 \\ \vdots \\ r_{N} \end{pmatrix}
 = \begin{pmatrix} r(x_0) \\ r(x_1) \\ \vdots \\ r(x_{N}) \end{pmatrix}\in \mathbb{R}^{n+1}.
```
In analogy, **solving a boundary value problem (1) in the discretised setting**
boils down to **finding a vector ${\mathbf u} \in \mathbb{R}^{n+1}$** such that its **entries ${u}_i$ are approximations to $u(x_i)$** for $i=0, \ldots, N$, i.e. the exact solution $u$ evaluated at the nodes $x_i$.
Before we discuss in details how to do this (see the section on the [collocation method](#Collocation-method) below), we first need to understand how one can express the derivative $u'(x_i)$ and $u''(x_i)$ in terms of the values $u(x_i)$.
"""

# ╔═╡ a184ac2c-97ae-4221-a64e-215e9cadeb80
md"""
To answer this point, let us consider a general function $f$ for which we know the values $f(x_i)$ on the nodal points for $i=0,\ldots,N$. Our goal is to find a vector  $\mathbf{g}$ such that $g_i ≈ f'(x_i)$, i.e. we have a **representation of the derivative of $f$ on the nodal points**.
From [chapter 09 on Numerical differentiation](https://teaching.matmat.org/numerical-analysis/09_Numerical_differentiation.html) we already know the *forward* finite difference formula, which suggests
```math
g_i = D_h^+ f(x_i) = \frac{f_{i+1} - f_i}{h} \qquad \text{for $i=0, \ldots, N-1$}.
```
Notably we cannot apply this formula to $i=N$ because then the formula would refer to the value $f_{N+1}$, which we do not know.
However, in this setting the *backwards* finite-difference formula applies:
```math
g_N = D_h^- f(x_N) = \frac{f_{N} - f_{N-1}}{h}
```
Collecting the function values of $f$ in a vector
```math
\mathbf{f} = \begin{pmatrix} f(x_0) \\ f(x_1) \\ \vdots \\ f(x_N) \end{pmatrix}\in \mathbb{R}^{N+1}
```
leads to the vector equation
```math
\begin{pmatrix} f'(x_0) \\ f'(x_1) \\ \vdots \\ f'(x_N) \end{pmatrix}
\approx \mathbf{g} = \widetilde{\mathbf{D}}_{x} \mathbf{f},
\qquad
\widetilde{\mathbf{D}}_{x} = \frac{1}{h} \begin{pmatrix}
-1 & 1 \\
& -1 & 1 \\
&& \ddots & \ddots \\
&&& -1 & 1 \\
&&& -1 & 1
\end{pmatrix}.
```
Notice that for the square matrix $\widetilde{\mathbf{D}}_{x} \in \mathbb{R}^{(N+1)\times (N+1)}$ all elements we do not show are zero.

Here, the matrix $\widetilde{\mathbf{D}}_{x} \in \mathbb{R}^{(N+1)\times (N+1)}$ provides a discrete approximation of the derivative operator $\frac{d}{dx}$: it maps a vector $\mathbf{f}$ (representing a discretised form of the function $f$) to a vector $\mathbf{g}$, which represents a discretised form of the derivative $f'$.
"""

# ╔═╡ 3be4502b-95af-4ca8-8fd4-4c8cb659a858
md"""
In the construction of the matrix $\widetilde{\mathbf{D}}_{x}$ we used only first-order finite-difference formulas. As a result $\widetilde{\mathbf{D}}_{x}$ only provides a first-order approximation to the derivative operator.

A more accurate construction is to employ a second-order finite difference formula wherever possible and fall back to the second-order forward / second-order backward formulas on the boundary points. This leads to

!!! info "Definition: Second-order differentiation matrix for first derivative"
	A second-order differentiation matrix for the **first** derivative is
	```math
	\mathbf{D}_x = \frac{1}{h} \begin{pmatrix}
	-\frac{3}{2} & 2 & -\frac12 \\
	-\frac12 & 0 & \frac12 \\
	& \ddots& \ddots & \ddots \\
	&&-\frac12 & 0 & \frac12 \\
	&& \frac12 & -2 & \frac32
	\end{pmatrix} \in \mathbb{R}^{(N+1)\times(N+1)}
	```
	where the nodes are given in equation (3).
"""

# ╔═╡ f44f3230-0062-49db-a488-f19910b58792
md"""
For solving boundary value problems (1) we not only need to be able to numerically compute first derivatives, but we also need to be able to compute second derivatives as well, that is we seek a differentiation matrix $\mathbf{D}_{xx}$ such that
```math
\begin{pmatrix}
	f''(x_0) \\ f''(x_1) \\ \vdots \\ f''(x_N)
\end{pmatrix}
≈
\mathbf{D}_{xx}
\begin{pmatrix}
	f(x_0) \\ f(x_1) \\ \vdots \\ f(x_N).
\end{pmatrix}
```
Using a second-order central finite difference scheme where possible and else a second-order one-sided scheme gives:

!!! info "Definition: Second-order differentiation matrix for second derivative"
	A second-order differentiation matrix for the **second** derivative is
	```math
	\mathbf{D}_{xx} = \frac{1}{h^2} \begin{pmatrix}
	2 & -5 & 4 & -1 \\
	1 & -2 & 1 \\
	& 1 & -2 & 1 \\
	&& \ddots& \ddots & \ddots \\
	&&&1 & -2 & 1 \\
	&& -1 & 4 & -5 & 2
	\end{pmatrix} \in \mathbb{R}^{(N+1)\times(N+1)}
	```
	where the nodes are given in equation (3).

An implementation returning the second-order matrices $\mathbf D_x$ and $\mathbf D_{xx}$ for approximating $\frac{d}{dx}$ and $\frac{d^2}{dx^2}$ is given below:
"""

# ╔═╡ 393313f4-f2b0-4f2e-b175-ae24c0514e58
function diffmats(N, a, b)
	# N:  Number of subintervals for discretisation, i.e. there will be N+1 nodes 
	# a:  Starting point of the computational domain [a, b]
	# b:  End point of the computational domain

	h = (b-a) / N
	x = [a + i*h for i in 0:N]  # Nodes

	# Fill Dₓ with central finite difference formula (correct for most rows)
	dx_sub   = fill(-1/2 / h, N)  # Subdiagonal   of Dₓ
	dx_super = fill( 1/2 / h, N)  # Superdiagonal of Dₓ
	Dₓ = diagm(-1 => dx_sub, 1 => dx_super)

	# Adjust first and last row of Dₓ
	Dₓ[1, 1:3]         = [-3/2,  2, -1/2] / h  # Fix first row (forward  finite diff)
	Dₓ[end, end-2:end] = [ 1/2, -2,  3/2] / h  # Fix last  row (backward finite diff)

	# Fill Dₓₓ with central finite differences (correct for most rows)
	dxx_main = fill(-2 / h^2, N+1)  # Main diagonal  of Dₓₓ
	dxx_side = fill( 1 / h^2, N)    # Side diagonals of Dₓₓ
	Dₓₓ = diagm(-1 => dxx_side, 0 => dxx_main, 1 => dxx_side)

	# Adjust first and last row of Dₓₓ
	Dₓₓ[1,   1:4]       = [2, -5, 4, -1] / h^2  # First row (forward  finite diff)
	Dₓₓ[end, end-3:end] = [-1, 4, -5, 2] / h^2  # Last row  (backward finite diff)

	(; x, Dₓ, Dₓₓ, h)
end

# ╔═╡ e2ed03ad-f720-4e29-afdc-1c1c6bc233c0
md"""
!!! warning "Example: Using differentiation matrices"
	In this example we want to test the first-order and second-order differentiation matrices $D_x$ and $D_{xx}$. We consider the following function
	```math
	f_\text{exp} = x + e^{\sin(4x)}
	```
	on the interval $[-1, 1]$.
	The exact first and second derivatives are
	```math
	\begin{aligned}
	f_\text{exp}'(x)  &= 1 + 4\cos(4x) e^{\sin(4x)} \\
	f_\text{exp}''(x) &= 4 (4\cos^2(4x) - 4 \sin(4x)) e^{\sin(4x)}
	\end{aligned}
	```
	We start by discretising on $20$ equispaced nodes (i.e. $N=19$ subintervals):
"""

# ╔═╡ 93bf2a7a-a926-4e0d-b127-889a0b82dcd2
begin
	fexp(x)     = x + exp( sin(4x) )
	dx_fexp(x)  = 1 + 4cos(4x) * exp( sin(4x) )  # Reference derivative
	dx2_fexp(x) = 4 * ( 4cos(4x)^2 - 4sin(4x) ) * exp( sin(4x) )  # and 2nd deriv.
end;

# ╔═╡ 94b97d6c-6c45-4887-b15e-5731ae5bbe9c
t, Dₓ, Dₓₓ = diffmats(19, -1, 1)

# ╔═╡ 80df1d2f-24f9-4503-9075-eb5a37a486bd
md"and we evaluate $f_\text{exp}$ on the nodes:"

# ╔═╡ cdc1fc28-714b-45e1-ad51-31c6d70dcce3
fexp_t = fexp.(t)

# ╔═╡ 22876147-cac6-4c56-921a-8d86e605bdeb
md"Next we obtain approximate derivatives using the matrix-vector multiplication with the differentiation matrices:"

# ╔═╡ b0c8cf85-b034-41c6-a808-6562ace45768
begin
	dx_fexp_t  = Dₓ  * fexp_t
	dx2_fexp_t = Dₓₓ * fexp_t
end;

# ╔═╡ 0e66bac1-a66a-47e9-aa22-b1f95ce2e45d
md"Comparing against the reference values shows a rather poor accuracy:"

# ╔═╡ 796c699d-f6d3-40b3-aa48-fff0d627bb51
let
	p = plot(dx_fexp, xlims=[-1.05, 1.05], xaxis=L"x", yaxis=L"f'(x)", label="reference", lw=2, legend=:topleft)
	p = scatter!(p, t, dx_fexp_t, label="approximation")

	q = plot(dx2_fexp, xlims=[-1.05, 1.05], xaxis=L"x", yaxis=L"f''(x)", label="reference", lw=2, legend=false)
	q = scatter!(q, t, dx2_fexp_t, label="approximation")

	plot(p, q)
end

# ╔═╡ 78cf3006-47c8-48c2-be3a-2a518ff4731d
md"""
Since we took a second-order finite difference formula, we would expect a second-order convergence. Plotting the maximal error of the approximated derivative values `dx_fexp_t` respectively `dx2_fexp_t`, that is
```julia
maximum(abs.(dx_fexp.(t)  - dx_fexp_t))
```
and
```julia
maximum(abs.(dx2_fexp.(t) - dx2_fexp_t))
```
against the subinterval size $h$ on a log-log plot confirms this:
"""

# ╔═╡ ba98bfcb-a716-49e2-bc0b-7294bdaff012
let
	Ns = round.(Int, 2.0 .^ (4:0.5:11))
	hs = zeros(length(Ns))
	error_df  = zeros(length(Ns))
	error_df2 = zeros(length(Ns))
	for k in 1:length(Ns)
		t, Dₓ, Dₓₓ, h = diffmats(Ns[k], -1, 1)		
		fexp_t = fexp.(t)
		error_df[k]  = maximum(abs.(dx_fexp.(t)  - Dₓ  * fexp_t))
		error_df2[k] = maximum(abs.(dx2_fexp.(t) - Dₓₓ * fexp_t))
		hs[k] = h
	end

	p = plot(; xaxis=:log, xlabel=L"h", yaxis=:log, ylabel="maximal error", xflip=true)
	plot!(p, hs, error_df,  mark=:o, lw=2, label=L"f'")
	plot!(p, hs, error_df2, mark=:o, lw=2, label=L"f''")

	# Guide line
	plot!(p, hs, 30 * hs .^ 2, lw=2, ls=:dash, c=:black, label=L"h^2")
	
	xticks!(p, 10 .^ (-1:-0.5:-3))
	yticks!(p, 10 .^ (1:-1.0:-4))

	p
end

# ╔═╡ 8d1cf059-05b1-44da-a51c-0157bb06efe9
md"""
## Collocation method

We will now construct a numerical method based on the finite differences differentiation matrices for solving the linear boundary value problem
```math
\left\{
\begin{aligned}
u''(x) + p(x) \, u'(x) + q(x) \, u(x) &= r(x) \qquad x \in (a, b),\\
u(a) &= γ_a \\
u(b) &= γ_b
\end{aligned}
\right.
```
with equispaced nodes
```math
x_i = a + i \, h, \qquad \forall i = 0, \ldots, N, \qquad h = \frac{b-a}{N}.
```
Rather than solving for the function $u : [a,b] \to \mathbb{R}$
we will instead seek the vector $\tilde{\mathbf u} \in \mathbb{R}^{n+1}$
of its approximate function values at the nodes $x_i$:
```math
\tilde{\mathbf{u}} = \begin{pmatrix} \tilde{u}_0 \\ \tilde{u}_1 \\ \vdots \\ \tilde{u}_{N} \end{pmatrix}
\textcolor{red}{\approx}
\begin{pmatrix} u(x_0) \\ u(x_1) \\ \vdots \\ u (x_{N}) \end{pmatrix}.
```
The key idea of **collocation** is to use the BVP equation (1) evaluated *at the same nodal points* $x_i$ to impose additional constraints on the solution vector. Concretely we use the differentiation matrices to approximate derivatives of the solution, i.e.
```math
\tag{4}
\begin{pmatrix} u'(x_0) \\ u'(x_1) \\ \vdots \\ u'(x_{N}) \end{pmatrix}
\textcolor{red}{\approx}
\tilde{\mathbf{u}}' := \mathbf{D}_x \tilde{\mathbf{u}}
\quad \text{and}\quad
\begin{pmatrix} u''(x_0) \\ u''(x_1) \\ \vdots \\ u''(x_{N}) \end{pmatrix}
\textcolor{red}{\approx}
\tilde{\mathbf{u}}'' := \mathbf{D}_{xx} \tilde{\mathbf{u}}
```
"""

# ╔═╡ 4f7bd62c-792c-4abb-ac93-28bb0169e739
md"""
With these definitions the discrete form of the first line of equation (1)
therefore
becomes
```math
\tilde{\mathbf{u}}'' + \mathbf{P} \tilde{\mathbf{u}}' + \mathbf{Q} \tilde{\mathbf{u}} = \mathbf{r}
```
where
```math
\mathbf{P} = \begin{pmatrix} p(x_0) \\ &\ddots \\ &&p(x_N) \end{pmatrix},
\quad
\mathbf{Q} = \begin{pmatrix} q(x_0) \\ &\ddots \\ &&q(x_N) \end{pmatrix},
\quad
\mathbf{r} = \begin{pmatrix} r(x_0) \\\vdots \\ r(x_{N}) \end{pmatrix}.
```
Using equation (4) we can rewrite this as
```math
\tag{5}
\mathbf{L} \tilde{\mathbf{u}} = \mathbf{r}
\quad \text{where} \quad \mathbf{L} = \mathbf{D}_{xx} + \mathbf{P} \mathbf{D}_x + \mathbf{Q}
```
where we note $\mathbf{L} \in \mathbb{R}^{(N+1)\times(N+1)}$,
i.e. this is a linear system with $N+1$ equations and $N+1$ unknowns.
Notice further that since $\mathbf{D}_{xx}$, $\mathbf{D}_{x}$, $\mathbf{P}$, $\mathbf{Q}$ are all sparse, banded matrices also $\mathbf{L}$ is sparse and banded.
"""

# ╔═╡ 41326216-32f2-412f-b1a4-45c287f3b640
md"""
The final step is to notice that the boundary conditions (last two lines of equation (1)) take the form of two additional conditions
```math
u_0 = γ_a, \qquad u_N = γ_b,
```
which are so far not yet taken into account in equation (5).
The easiest way to achieve this is to replace the first and last equation in the linear system of equation (5) exactly by these two equations given by the boundary conditions.
To make this clear, let us first write equation (5) more explicitly as
```math
\begin{pmatrix}
L_{00} & L_{01} & \cdots & L_{0N} \\
L_{10} & L_{11} & \cdots & L_{1N} \\
\vdots & \vdots & \ddots & \vdots \\
L_{N-1,0} & L_{N-1,1} & \cdots & L_{N-1,N} \\
L_{N0} & L_{N1} & \cdots & L_{NN}
\end{pmatrix}
\begin{pmatrix} \tilde{u}_0 \\ \tilde{u}_1 \\ \vdots \\ \tilde{u}_{N-1} \\ \tilde{u}_N \end{pmatrix}
=
\begin{pmatrix} r_0 \\ r_1 \\ \vdots \\ r_{N-1} \\ r_N \end{pmatrix}.
```
Now we replace the first equation by $u_0 = γ_a$
and the last equation by $u_N = γ_b$, this yields
```math
\begin{pmatrix}
\textcolor{red}{1} & \textcolor{red}{0} & \textcolor{red}{\cdots} & \textcolor{red}{0} \\
L_{10} & L_{11} & \cdots & L_{1N} \\
\vdots & \vdots & \ddots & \vdots \\
L_{N-1,0} & L_{N-1,1} & \cdots & L_{N-1,N} \\
\textcolor{red}{0} & \textcolor{red}{0} &  \textcolor{red}{\cdots}& \textcolor{red}{1}
\end{pmatrix}
\begin{pmatrix} \tilde{u}_0 \\ \tilde{u}_1 \\ \vdots \\ \tilde{u}_{N-1} \\ \tilde{u}_N \end{pmatrix}
=
\begin{pmatrix} \textcolor{red}{γ_a} \\ r_1 \\ \vdots \\ r_{N-1} \\ \textcolor{red}{γ_b} \end{pmatrix}.
```
Finally we make the notation more compact by introducing the $(N-1) \times (N+1)$ matrix
```math
\mathbf{E} = \left(\begin{array}{c|ccc|c}
0     & 1& & &  0 \\
\vdots&  &\ddots& &  \vdots \\ 
0     &  & & 1&  0
\end{array} \right)
```
which is constructed such that the product $\mathbf{E}\mathbf{L}$ of this matrix with an $(N+1) \times (N+1)$ matrix $\mathbf{L}$ deletes the first and last rows of $\mathbf{L}$. Similarly $\mathbf{E} \mathbf{r}$ deletes the first and last rows of a matrix $\mathbf{r}$. Let us demonstrate this on a $4\times4$ example:
"""

# ╔═╡ 51552b34-01b4-4ef0-866c-6ff6c33a3462
E_2by4 = [0 1 0 0;
		  0 0 1 0]

# ╔═╡ 5d361452-2eef-4503-a80f-151277d3b855
M = [00 01 02 03;
	 10 11 12 13;
	 20 21 22 23;
     30 31 32 33]

# ╔═╡ 8ddf3426-1d5d-42ee-b44a-3c4ea1ed5f4c
E_2by4 * M

# ╔═╡ b2418c32-854f-469c-ac9d-4fb7b93d5b79
md"""
Using $\mathbf{E}$ we can write the linear system more compactly as
```math
\tag{6}
\underbrace{
\begin{pmatrix}
\textcolor{red}{1} & \textcolor{red}{0} & \textcolor{red}{\cdots} & \textcolor{red}{0} & \textcolor{red}{0} \\
\hline
&&\mathbf{E}\mathbf{L} \\
\hline
\textcolor{red}{0} & \textcolor{red}{0} &  \textcolor{red}{\cdots}&  \textcolor{red}{0} & \textcolor{red}{1}
\end{pmatrix}}_{=\mathbf{A}}
\tilde{\mathbf{u}}
=
\underbrace{\begin{pmatrix} \textcolor{red}{γ_a} \\\hline \mathbf{E}\mathbf{r} \\\hline \textcolor{red}{γ_b} \end{pmatrix}}_{=\mathbf{b}}.
```
This is an $(N+1) \times (N+1)$ linear system $\mathbf{A} \mathbf{u} = \mathbf{b}$ with unknown $\mathbf{u}$ and --- as we will observe in more detail below --- a sparse, banded matrix $\mathbf{A}$.
Problem (6) can therefore be **efficiently solved** using [direct methods based on (sparse) LU factorisation (chapter 5)](https://teaching.matmat.org/numerical-analysis/05_Direct_methods.html) or an [iterative approaches (chapter 6)](https://teaching.matmat.org/numerical-analysis/06_Iterative_methods.html).
"""

# ╔═╡ f1991ea0-9270-4c27-8032-4f4db84330e4
md"""
!!! warning "Example: Setting up a collocation linear system"
	As an example we discretise the boundary value problem
	```math
	u'' = u - x u', \quad u(0) = -2, \quad u(2) = 3.
	```

	- **Step 1: Selection of nodes:**
	  We choose an equispaced grid with $N=4$ subintervals. As a result we have
	  the following discretisation nodes
	  ```math
	  h = \frac{2-0}{4} = \frac12 \quad \Rightarrow \quad x_0 = 0, x_1 = \frac12, x_2 = 1, x_3=\frac32, x_4 = 2
	  ```

	- **Step 2: Setup of operator $L$:**
	  Bringing the differential equation into the form of equation (1),
	  i.e. $u'' + x u' - u = 0$ we notice $p(x) = x$, $q(x) = -1$ and $r(x) = 0$.
	  Therefore
	  ```math
	  \mathbf{P} = \begin{pmatrix}0 \\&1/2\\&&1\\&&&3/2\\&&&&2\end{pmatrix},
	  \quad
	  \mathbf{Q} = -\begin{pmatrix}1 \\&1\\&&1\\&&&1\\&&&&1\end{pmatrix},
	  \quad
	  \mathbf{r} = \begin{pmatrix}0\\0\\0\\0\\0\end{pmatrix}
	  ```
	  Using the second-order differentiation matrices $\mathbf{D}ₓ$ and $\mathbf{D}ₓₓ$ we therefore obtain:
	  ```math
	  \begin{aligned}
	  \mathbf{L} &= \mathbf{D}ₓₓ + \mathbf{P} \mathbf{D}ₓ + \mathbf{Q} \\
	  &= \frac{1}{(1/2)^2} \begin{pmatrix}
	  2 & -5 & 4 & -1 \\
	  1 & -2 & 1 \\
	  & 1 & -2 & 1 \\
	  && 1 & -2 & 1 \\
	  & -1 & 4 & -5 & 2
	  \end{pmatrix}
	  + \mathbf{P} \frac{1}{1/2} \begin{pmatrix}
	  -\frac{3}{2} & 2 & -\frac12 \\
	  -\frac12 & 0 & \frac12 \\
	  & -\frac12 & 0 & \frac12 \\
	  &&-\frac12 & 0 & \frac12 \\
	  && \frac12 & -2 & \frac32
	  \end{pmatrix}
	  + \mathbf{Q} \\
	  &= \begin{pmatrix}
	  7   & -20 & 16  & -4 & 0 \\
	  7/2 & -9  & 9/2 &  0 & 0 \\
	  0   & 3   & -9  & 5  & 0 \\
	  0   & 0   & 5/2 & -9 & 11/2 \\
	  0   & -4  & 18  & -28 & 13
	  \end{pmatrix}
	  \end{aligned}
	  ```

	- **Step 3: Application of boundary conditions:**
	  To apply the boundary conditions we change the first and last lines
	  of the discretised differential equation $\mathbf{L} \mathbf{u} = \mathbf{r}$
	  arriving at the linear system $\mathbf{A} \mathbf{u} = \mathbf{b}$ with
	  ```math
	  \mathbf{A} = \begin{pmatrix}
	  1   & 0 & 0  & 0 & 0 \\
	  7/2 & -9  & 9/2 &  0 & 0 \\
	  0   & 3   & -9  & 5  & 0 \\
	  0   & 0   & 5/2 & -9 & 11/2 \\
	  0   & 0  & 0  & 0 & 1
	  \end{pmatrix} \qquad
	  \mathbf{b} = \begin{pmatrix} -2 \\ 0 \\ 0 \\ 0 \\ 3 \end{pmatrix}
	  ```

	- **Step 4: Solution of linear system (6):**
	  The collocated system $\mathbf{A} \mathbf{u} = \mathbf{b}$ we can now
	  solve using LU factorisation:
"""

# ╔═╡ 27338f87-6f0a-4577-8cec-8eaba0c9d411
let
	A = [1     0  0     0  0;
		 7/2  -9  9/2   0  0;
		 0     3  -9    5  0;
		 0     0  5/2  -9  11/2;
		 0     0  0     0    1]
	b = [-2;
		  0;
		  0;
		  0;
		  3]
	A \ b
end

# ╔═╡ b1d08025-3557-4045-9646-a42faae8b91c
md"""
Finally we present a generic implementation of the collocation method for solving linear boundary value problems based on the `diffmats` function and LU factorisation to solve the arising linear system:
"""

# ╔═╡ df13e0b3-ff18-4dc7-8d46-0a2f786cdbeb
function fd_bvp(N, a, b, p, q, r, γa, γb)
	# Solves the differential equation
	#   u''(x) + p(x) u'(x) + q(x) u(x) = r(x)
	# in the domain [a, b] with Dirichlet boundary conditions u(a) = γa, u(b) = γb
	#
	# N:  Number of subintervals for discretisation, i.e. there will be N+1 nodes 
	# a:  Starting point of the computational domain [a, b]
	# b:  End point of the computational domain
	# p, q, r:  Julia functions to define the unknowns in the equation above
	# γa, γb:   Dirichlet boundary values

	(; x, Dₓ, Dₓₓ, h) = diffmats(N, a, b)

	# Assemble discretised differential equation operator L
	P = Diagonal(p.(x))
	Q = Diagonal(q.(x))
	L = Dₓₓ + P * Dₓ + Q

	# Apply bounary conditions
	A = zeros(N+1, N+1)
	A[2:end-1, :] = L[2:end-1, :]  # Set, but exclude first and last row of L 
	A[1, 1] = A[end, end] = 1      # Set upper left and lower-right corners

	b = zeros(N+1)
	b[2:end-1] = r.(x[2:end-1])  # Set, but exclude first and last entry
	b[1]   = γa  # Set first entry
	b[end] = γb  # Set last entry
	
	# Solve the linear system
	u = A \ b

	(; x, u, h, A, b, Dₓ, Dₓₓ)
end

# ╔═╡ 79d958dd-7b0e-428f-ba17-eb9983857504
md"""
To test this implementation we consider the following heat equation
```math
\left\{
\begin{aligned}
- u''(x) &= \sin(x) \qquad x \in (0, 2π)\\
u(0) &= 1, \quad u(2π) = 2,
\end{aligned}
\right.
```
which is to be solved in the interval $[0, 2π]$.

Rewriting this equation as $u''(x) = -\sin(x)$ and comparing with the
interface of the function `fd_bvp` we notice the following parameters:
"""

# ╔═╡ 3bf1ee87-d42c-429d-8a07-f080f369f61f
begin
	p(x) = 0.0
	q(x) = 0.0
	r(x) = -sin(x)

	a = 0.0  # Interval begin
	b = 2π   # Interval end
	
	γa = 1   # Left  boundary value
	γb = 2   # Right boundary value
end;

# ╔═╡ 9f87dc91-c763-46b8-9277-2cc5f9a92cc0
md"""
In this case the exact solution can be found analytically as
"""

# ╔═╡ 785b3c50-b008-4ec3-943e-d2e874942f7f
u_exact(x) = sin(x) + x/2π + 1

# ╔═╡ fadf5e83-d5b4-48b9-911b-7f1309f02234
md"We solve this using finite differences"

# ╔═╡ f32bc6ac-f84d-453d-8cab-c06491410c24
md"and plot the result:"

# ╔═╡ 0bb8239c-22c0-4420-8ae7-59c5c00436a2
md"`N = ` $(@bind N Slider(6:1:30; default=6, show_value=true))"

# ╔═╡ 452d1b69-42db-4ab0-889b-bc48ae0ec365
res_fd = fd_bvp(N, a, b, p, q, r, γa, γb);

# ╔═╡ 474dee17-627a-49b1-84e8-d692346d8aa0
let
	plot(u_exact; xlims=(a-0.1, b+0.1), legend=:topright,
	     label="reference", lw=2, xlabel=L"x", ylabel=L"u(t)",
	     title=L"-\frac{\partial^2 u}{\partial x^2} = \sin(x)")

	plot!(res_fd.x, res_fd.u; label="Finite differences N=$N", mark=:o, lw=2, ls=:dash)
end

# ╔═╡ 69ca61cc-5482-45c1-99fc-0ef780771759
md"`h = ` $(round((b-a)/N; digits=2))"

# ╔═╡ 814ac12b-5ed7-4512-b053-fbe729b90ce6
md"""
We observe that as the number of points `N` is increased,
the solution visually becomes more and more accurate.
The output of our numerical procedure are 
the computed points $u_0, u_1, \ldots u_N$.
These points approximate the true function values 
$u(x_0), u(x_1), \ldots, u(x_N)$ of $u$ evaluated at the nodal points.
A measure for the accuracy of our obtained solution is thus
the maximal deviation any of these computed points has from
from the true $u(x_i)$.
This measure is called the **global error** defined by
```math
|e_i| = |u(x_i) - u_i| \qquad \text{for $i = 1, \ldots, N$}.
```

Since we used a quadratic difference formulas for building
the differentiation matrices $Dₓ$ and $Dₓₓ$ we obtain
**quadratic convergence** in this global error with increasing $N$:
"""

# ╔═╡ 7a647ffe-000e-4bb7-a348-033cd75b1369
let
	Ns = [ round(Int, 5 * 10^k) for k in 0:0.25:2.5 ]
	errors = Float64[]
	for N in Ns
		res_fd = fd_bvp(N, a, b, p, q, r, γa, γb)
		x = res_fd.x
		u = res_fd.u
		error = [ u_exact(x[j]) - u[j] for j in 1:length(x)]
		push!(errors, maximum(error))		
	end
	
	p = plot(; title="Convergence of Dirichlet finite differences",
			   xaxis=:log, yaxis=:log,
	           xlabel=L"N", ylabel="Largest global error",
	           legend=:bottomleft,
	           xlims=(1, 5e3), ylims=(1e-7, 10))
	plot!(p, Ns, errors; lw=2, mark=:o, c=1, label="Error")
	plot!(p, Ns, 0.02(first(Ns) ./ Ns).^2, ls=:dash, lw=2, label=L"N^{-2} = O(h^2)", c=1)

	xticks!(p, 10.0 .^ (0:1:3))
	yticks!(p, 10.0 .^ (0:-2:-6))
	p
end

# ╔═╡ 52961854-695e-44e7-b8e2-20ed7ff07cda
md"""
In line with this discussion we define convergence for numerical schemes for boundary value problems as follows:

!!! info "Definition: Convergence order for boundary value problems"
	A numerical scheme approximating a boundary value problem (9)
	**converges with order $p$** if there exists a constant $α>0$
	such that
	```math
	\max_{j=0,\ldots,N+1} |u(x_j) - u_j| ≤ α \, h^p
	```
	provided that the solution $u$ is sufficiently regular.
"""

# ╔═╡ 359e61a3-d7fc-402f-ad6c-24b4301885d6
md"""
Since we used a quadratic difference formulas for building
the differentiation matrices $Dₓ$ and $Dₓₓ$ we obtain indeed quadratic convergence $p=2$, as is formally summarised in the following theorem.
"""

# ╔═╡ f05ac977-45b4-4bd1-baf0-8a2d50525ec9
md"""
!!! info "Theorem 1: Global error of finite differences"
	Under the assumption that the exact solution $u$ of the Dirichlet boundary value problem (1) is **four times differentiable**,
	the finite-difference scheme of equations (5) and (6)
	involving the **second-order differentiation matrices**
	$\mathbf{D}_{x}$ and $\mathbf{D}_{xx}$
	converges as 
	```math
	\max_{j=1,\ldots,N} |u(x_j) - u_j| ≤ α \, \|u''''\|_\infty \, h^2
	```
	where $α > 0$ is a constant that independent of $u$ and $h$.
	It thus achieves **quadratic convergence**.
"""

# ╔═╡ b3a1222b-52cc-49ff-af74-cc472bd40999
md"""
The proof of this result goes beyond the scope of this course and can for example be found in chapter 7.4 of *Introduction to Numerical Analysis* by Stoer and Bulirsch. The key observation to take a way from this developments is:

!!! info "Observation: Convergence order of BVP scheme and FD formula"
	Provided our solution function $u$ can be infinitely often differentiated continously, than we can **obtain $p$-th order convergence** (as $h\to0$ or $N\to\infty$) of the solution vector  $\tilde{\mathbf u}$ obtained in (6)
	provided that we consistently use **$p$-th order finite-difference formulas** to build *both* the differentiation matrices $Dₓ$ and $Dₓₓ$ in (5).
"""

# ╔═╡ d9fb073d-042c-4724-834f-bcdb0cc4aec6
md"""
### Numerical stability

We want to push the convergence plot from above a little further. For this we employ an implementation of `fd_bvp` that directly sets up the matrix `A` as a sparse matrix. Unfold the cell below to see the implementation.
"""

# ╔═╡ 75311119-2d3a-4ca4-bddb-5a595ffb09a4
function fd_bvp_fast(N, a, b, p, q, r, γa, γb)
	@assert iszero(p(randn()))  # Assume p(x) = 0.0 for now

	# Faster version of fd_bvp by building a sparse matrix A
	h = (b-a) / N
	x = [a + i*h for i in 0:N]  # Nodes

	# Build Dₓₓ, but without first and last row
	dxx_main  = [0.0; fill(-2.0, N-1); 0.0] ./ h^2
	dxx_sub   = [fill(1.0, N-1); 0.0] ./ h^2
	dxx_super = [0.0; fill(1.0, N-1)] ./ h^2
	Dₓₓ = Tridiagonal(dxx_sub, dxx_main, dxx_super)

	# Assemble discretised differential equation operator
	# (but without first and last row)
	Q = Diagonal([0; q.(x[2:end-1]); 0])
	A = Dₓₓ + Q

	# Apply boundary conditions
	A[1, 1] = A[end, end] = 1.0

	b = [γa; r.(x[2:end-1]); γb]
	
	# Solve the linear system
	u = A \ b

	(; x, u, h, A, b, Dₓ, Dₓₓ)
end

# ╔═╡ b482d8e2-901a-4d8b-85af-e0fbc5c46d4f
fd_bvp_fast(N, a, b, p, q, r, γa, γb).A

# ╔═╡ 2fd3f387-b950-462f-a06a-eabd458fe10e
let
	Ns = [ round(Int, 5 * 10^k) for k in 2:0.1:5.5 ]
	errors = Float64[]
	for N in Ns
		res_fd = fd_bvp_fast(N, a, b, p, q, r, γa, γb)
		x = res_fd.x
		u = res_fd.u
		error = [ u_exact(x[j]) - u[j] for j in 1:length(x)]
		push!(errors, maximum(error))		
	end
	
	p = plot(; title="Convergence of Dirichlet finite differences",
			   xaxis=:log, yaxis=:log,
	           xlabel=L"N", ylabel="Largest global error",
	           legend=:bottomleft,
	           xlims=(100, 3e6), ylims=(1e-13, 10))
	plot!(p, Ns, errors; lw=2, mark=:o, c=1, label="Error")
	plot!(Ns, 0.002(first(Ns) ./ Ns).^2, ls=:dash, lw=2, label=L"O(N^{-2})", c=1)

	xticks!(p, 10.0 .^ (2:1:6))
	yticks!(p, 10.0 .^ (0:-2:-12))
	p
end

# ╔═╡ 8039250a-c79a-4009-9a53-81307b8e0ace
md"""
While initially the convergence thus nicely follows the expected convergence curve, **for larger $N$ the convergence degrades and the error starts increasing again**.

Similar to our discussion on numerical stability
in the [chapter on numerical differentiation](https://teaching.matmat.org/numerical-analysis/09_Numerical_differentiation.html)
this error plot is the result of a balance between two error contributions:
- The **discretisation error** due to the choice of $N$, where as $N$ gets larger
  this error **decreases** as $O(N^{-2})$.
- The **error due to finite floating-point precision**, which turns out to increase
  as $N$ increases.

With respect to the second error contribution
the underlying problem is related to solving the linear system $\textbf{A} \textbf{u} = \textbf{b}$ in Algorithm 1. 
As it turns out the **condition number of the system matrix $\mathbf{A}$ grows** as $O(N^2) = O(h^{-2})$:
"""

# ╔═╡ f3819d92-c61a-49b2-a77b-9837cbed0cae
let
	Ns = [5, 10, 50, 80, 200, 500, 1000]
	condition_numbers = [cond(fd_bvp_fast(N, a, b, p, q, r, γa, γb).A) for N in Ns]

	p = plot(Ns, condition_numbers; mark=:o, yaxis=:log, xaxis=:log, lw=2, xlabel=L"N", ylabel=L"κ(A)", legend=:topleft, label=L"Condition number of $A$", ylims=(1, 10^13), xlims=(1, 3e6), yticks=(10.0 .^ (0:2:12)), xticks=(10.0 .^ (0:1:6)))
	plot!(p, Ns -> Ns^2, ls=:dash, lw=2, label=L"O(N^2) = O(h^{-2})")
end

# ╔═╡ 30b8bac2-7c95-4fee-bbde-b04a36258857
md"""
As a result **the finer we make the discretisation**, i.e. the **larger $N$**,
the harder it becomes to solve the system $\mathbf{A} \mathbf{u} = \mathbf{b}$,
such that the **floating-point error in computing the solution $\mathbf{u}$ itself** increases quadratically.

To make this clear let's consider $N = 10^5$.
At this stage the condition number is about $10^{10}$,
such that the typical $10^{-16}$ error we make when representing
any floating-point operation (such as the computation of $\mathbf{b}$)
gets amplified to a relative error in the solution $\mathbf{u}$
of around $10^{-6}$: at most around $6$ digits of the solution
can be known exactly.

As a result around $N = 10^5$ the floating-point error becomes
the leading error contribution,
such that increasing $N$ beyond $10^5$ thus causes the total
observed error of the numerical procedure to increase again.
"""

# ╔═╡ 858a60ba-8d0c-43d8-908d-0bb413c79285
TODO("Needs an adaptive pass from here.")

# ╔═╡ 23307b35-c76e-487d-906a-18f65b691d79
md"""
## Galerkin methods 

Let us return to the general heat equation using Dirichlet boundary conditions:

```math
\tag{17}
\left\{
\begin{aligned}
-k \, \frac{\partial^2 u}{\partial x^2}(x) &=  f(x) \qquad x \in (0, L).\\
u(0) &= b_0, \quad u(L) = b_L
\end{aligned}
\right.
```

We saw that **with finite differences** the increasing condition number of the system matrix $\mathbf{A}$ as we increase $N$ makes it **impossible to obtain the solution to arbitrary accuracy**.
In the above example we were unable to obtain a solution
to higher accuracy than $10^{-10}$, irrespective of what value for $N$ we chose.
In this section we will develop an alternative approach to solve boundary value problems, which is more general.
"""

# ╔═╡ 7ae614fb-f2e0-46af-a557-6a378f6553e8
md"""
Let us assume we are given some smooth function $\psi(x)$,
the so-called **test function**,
for which we will additionally assume that it **vanishes at the boundary**,
i.e. $\psi(0) = \psi(L) = 0$.
If we multiply the first line of (17) by this function and integrate over the full computational domain $[0, L]$, this yields:
```math
\tag{18}
\begin{aligned}
\int_{0}^L f(x) \psi(x)\, dx
&= \int_0^L -k\, u''(x) \psi(x)\, dx\\
&\stackrel{(\ast)}{=} [-k\, u'(x) \psi(x)]_0^L + \int_0^L k u'(x) \psi'(x)\, dx\\
&= \int_0^L k u'(x) \psi'(x)\, dx
\end{aligned}
```
where we used partial integration in step $(\ast)$ and the property  $\psi(0) = \psi(L) = 0$ in the final step.
"""

# ╔═╡ 962ca311-c2ed-46eb-b043-827ecfb64ac3
md"""
Let us define the **set of all test functions** as
```math
V_0 = \Big\{ \psi : [0, L] \to \mathbb{R} \ \Big|\ \text{$\psi$ smooth and $\psi(0) = \psi(L) = 0$} \Big\}.
```
If $u : [0, L] \to \mathbb{R}$ is a solution to (17) than our discussion implies that (18) is satisfied for all functions from $V_0$:
```math
\tag{19}
\left\{
\begin{aligned}
\int_0^L k\, u'(x)\, \psi'(x)\, dx &= \int_{0}^L f(x)\, \psi(x)\, dx
\qquad \forall\,\psi \in V_0\\
u(0) &= b_0, \quad u(L) = b_L
\end{aligned}
\right.
```
Equation (19) is known as the **weak form** of the 1D heat equation equation (17)
and solving it is interesting in its own right:
"""

# ╔═╡ 05999b01-4a41-480a-9f3a-2e04d97748b3
md"""
!!! info "Definition: Weak form and weak solution"
	If a function $u : [0, L] \to \mathbb{R}$ satisfies
	```math
	\left\{
	\begin{aligned}
	\int_0^L k\, u'(x)\, \psi'(x)\, dx &= \int_{0}^L f(x)\, \psi(x)\, dx
	\qquad \forall\,\psi \in V_0\\
	u(0) &= b_0, \quad u(L) = b_L,
	\end{aligned}
	\right.
	```
	for all choices of the test function $\psi$ we call $u$ a **weak solution**
	of the boundary value problem (17).
"""

# ╔═╡ 7a86c1db-e5be-4310-951a-4de51bec4e98
md"""
The original problem (17) is additionally called the **strong form** of the heat equation BVP and its solution $u$ the **strong solution**.

!!! info "Observation: Strong and weak solutions"
	If $u$ is a solution to the strong form (17) of a BVP,
	than it is also a solution to the weak form (19).
	However, **not every weak solution is also a strong solution**.
"""

# ╔═╡ 79816267-531b-4b99-bdea-1194574dd160
md"""
The weak form (19) looks unusual at first sight
and additionally one might wonder why it is useful,
since it can produce solutions, which do not satisfy
the original strong problem (17).

But surprisingly for many problems from physics or engineering
the situation turns out to be the reverse:
it turns out that the somewhat "looser" form provided
by **the weak formulation is actually more physical**. 
In some cases --- such as the atomistic modelling of materials ---
using the strong form  can lead to unphysical artifacts,
which are in fact avoided when using the weak form.

While we do not have time to discuss this further in this course,
the interested reader is referred to the master course [MATH-500: Error control in scientific modelling](https://teaching.matmat.org/error-control/)
where some of this is discussed.
"""

# ╔═╡ bbb839d2-599f-40e2-a1f7-2419d07c96a8
md"""
### Discretisation of the weak form

In order to solve (19) our goal is to rewrite the problem in the form of linear algebra --- vectors and matrices. For simplicity we will only consider the case $b_0=b_L=0$ in this section, i.e. we restrict ourselves to the **heat equation with a zero Dirichlet boundary**
```math
\tag{20}
\left\{
\begin{aligned}
-k \, \frac{\partial^2 u}{\partial x^2}(x) &=  f(x) \qquad x \in (0, L).\\
u(0) &= u(L) = 0
\end{aligned}
\right.
```
and **corresponding weak form**
```math
\tag{21}
\left\{
\begin{aligned}
\int_0^L k\, u'(x)\, \psi'(x)\, dx &= \int_{0}^L f(x)\, \psi(x)\, dx
\qquad \forall\,\psi \in V_0\\
u(0) &= u(L) = 0.
\end{aligned}
\right.
```
"""

# ╔═╡ 43ade623-1b58-4035-a1e5-803911063833
md"""
One difficulty in this equation is that the condition $\forall \psi \in V_0$ corresponds to an infinite number of constraints to satisfy (as the set $V_0$ has infinitely many members). One idea is to approximate this condition by satisfying only finitely many constraints.
"""

# ╔═╡ e6f3450d-6fb7-4c59-a209-ae3587a43bab
md"""
To do so we assume that **$\psi$ can be written as a linear combination**
```math
\tag{22}
\psi(x) = \sum_{i=1}^m ξ_i\, \varphi_i(x),
```
where $φ_1, φ_2, \ldots, φ_m$ are a selection of
$m$ linearly independent functions from $V_0$,
i.e. smooth functions which each satisfy $φ_i(0) = φ_i(L) = 0$.
Inserting (22) into the first line of (21) leads to
```math
\begin{aligned}
0 &= \int_0^L k\,u'(x)\,\psi'(x) - f(x)\psi(x)\,dx \\
  &= \sum_{i=1}^m ξ_i \left[ \int_0^L k\,u'(x)\,\varphi_i'(x)  - f(x)\,\varphi_i(x) \,dx\right]
\end{aligned}
```
which should be true **independent of the values of $ξ_i$**.
As a result the above condition can be achieved if and only if
the term in the **square bracket is zero all our functions $\varphi_i$**,
that is if
```math
\tag{23}
\forall i = 1, \ldots, m: \qquad
\int_0^L k\,u'(x)\,\varphi_i'(x)\,dx = \int_{0}^L f(x)\,\varphi_i(x)\,dx.
```
"""

# ╔═╡ c7610b00-7277-46c3-881d-85d417db2e9e
md"""
Notice, that the condition (23) is independent of the actual
values taken by the coefficients $\{ξ_i\}_{i=1}^m$.
Provided that (23) holds we thus do not need to worry about
determining the values of $\{ξ_i\}_{i=1}^m$.
In fact in all following development these coefficients
will play no further role.

With this development we can replace the infinite number of constraints encoded
by the $\forall \psi \in V_0$ in equation (21) by the approximate **version (23)**,
which **only involves $m < \infty$ conditions** to satisfy.
"""

# ╔═╡ 5e7387cc-f9eb-484e-ade7-c844bddeb716
md"""
Considering $u(x)$ we make an additional approximation,
namely that --- similar to $\psi$ --- this function can **also be
approximated by a linear combination of finitely many functions**.
Here, for simplicity we employ the same selection of functions,
that we used for $\psi$ leading to an ansatz
```math
\tag{24}
u(x) = \sum_{j=1}^m c_j φ_j(x).
```
Inserting (24) into (23) we obtain
```math
\forall i = 1, \ldots, m: \qquad
\int_0^L k\,\sum_{j=1}^m c_j φ_j'(x)\,φ_i'(x)\,dx = \int_{0}^L f(x) \, \varphi_i(x)\,dx
```
or arranged differently
```math
\tag{25}
\forall i = 1, \ldots, m: \qquad
\sum_{j=1}^m c_j \int_0^L k\,φ_j'(x)\,φ_i'(x)\,dx = \int_{0}^L f(x)\,\varphi_i(x)\, dx
```
"""

# ╔═╡ e32d65f1-9bf2-4043-811c-d47801d0f85e
md"""
Since $φ_j(0) = φ_j(L) = 0$ an immediate consequence is that
such a $u$ satisfies the boundary condition $u(0) = u(L) = 0$
independent of the choice of the coefficients $c_j$.
Solving equation (25) therefore automatically ensures that the second line of the weak problem (21) is satisfied.
Additionally this equation satisfies (23)
--- our approximation to the first line of (21).

In summary a **solution to equation (25)** thus **satisfies** (our approximation to)
**both conditions required to be a solution
to the weak problem** (21).
Equation (25) is thus given a special name:
"""

# ╔═╡ ef1ad23b-89bc-4b75-9560-4fa4fb96c802
md"""
!!! info "Definition: Galerkin conditions"
	Given a basis $φ_1, φ_2, \ldots, φ_m$
	of $m$ linearly independent functions
	from the test function set $V_0$,
	the linear combination $u = \sum_{j=1}^m c_j φ_j$
	is an approximate solution to the weak problem (21)
	if the coefficients $c_1, c_2, \ldots, c_m$
	satisfy the **Galerkin conditions**
	```math
	\tag{25}
	\forall i = 1, \ldots, m: \qquad
	\sum_{j=1}^m c_j \int_0^L k\,φ_j'(x)\,φ_i'(x)\,dx = \int_{0}^L f(x)\,\varphi_i(x)\, dx
	```

	Collecting the unknown coefficients into a vector $\mathbf{c} = (c_1, c_2, \ldots, c_m)^T \in \mathbb{R}^m$ and introducing
	the matrix $\mathbf{A} \in \mathbb{R}^{m\times m}$
	and the vector $\mathbf{f} \in \mathbb{R}^m$ with elements
	```math
	\begin{aligned}
	A_{ij} &= \int_0^L k\,φ_i'(x)\,φ_j'(x)\,dx, &
	f_i    &= \int_0^L f(x)\,φ_i(x)\,dx
	\end{aligned}
	```
	we can also write the Galerkin conditions more conveniently as the linear system
	```math
	\tag{26}
	\mathbf{A} \mathbf{c} = \mathbf{f},
	```
	which one needs to solve for $\mathbf{c}$.
"""

# ╔═╡ fcd48a38-eca2-44d9-a32e-9dd65c58fb89
md"""
### Example: Sine basis

Let us consider the Dirichlet heat equation
```math
\left\{
\begin{aligned}
- \frac{\partial^2 u}{\partial x^2}(x) &= x \qquad x \in (0, π)\\
u(0) &= u(L) = 0,
\end{aligned}
\right.
```
i.e. where $k=1$, $L=π$ and $f(x) = x$.
"""

# ╔═╡ 4d7ae336-be84-4960-835f-aa81ab91c4e2
md"""
For the basis functions $\varphi_j$ we select the sine basis
```math
\varphi_j(x) = \sin(j\,x) \qquad \text{for $j = 1, \ldots, m$}.
```
It is easy to see that all of them satisfy $\varphi_j(0) = \varphi_j(2π) = 0$.
With these basis functions
we can obtain the entries of $\mathbf{A}$ and $\mathbf{f}$
in equation (26) as
```math
\begin{aligned}
	A_{ij}
	&= \int_0^{π} \frac{d}{dx}\big[ \sin(i\, x)\big]\, \frac{d}{dx}\big[ \sin(j\, x)\big]\, dx \\
	&= i\,j \int_0^{π} \cos(i\,x) \cos(j\,x)\,dx = \frac{π\,i^2}{2}\, δ_{ij} = \left\{\begin{array}{ll} \frac{π\,i^2}{2} &i=j \\ 0 & i\neq j\end{array}\right.\\[1em]
	f_i &= \int_0^{π} x \sin(i\,x) \,dx = \frac{(-1)^{i+1}\, π}{i}
\end{aligned}
```
"""

# ╔═╡ 27a67eb3-76e5-4e23-a848-7a9e58485721
md"""
For example for $m=5$ we find the matrix and vector
```math
\mathbf{A} = \frac{π}{2} \begin{pmatrix}
1 \\
& 4 \\
&&9\\
&&&16\\
&&&&25
\end{pmatrix}
\qquad
\textbf{f} = π \begin{pmatrix}
	1 \\ -\frac12 \\ \frac13 \\ -\frac14 \\ \frac15
\end{pmatrix}.
```
"""

# ╔═╡ 52716cfc-6dea-44ad-aa96-842332fe607b
md"""
Since $\mathbf{A}$ is diagonal we can directly compute
the coefficient vector $\mathbf{c}$ as 
```math
\mathbf{c} = 
\begin{pmatrix}
f_1 / A_{11} \\
f_2 / A_{22} \\
f_3 / A_{33} \\
f_4 / A_{44} \\
f_5 / A_{55}
\end{pmatrix}
= \begin{pmatrix}
2 \\ -\frac28 \\ \frac2{27} \\ -\frac2{64} \\ \frac2{125}
\end{pmatrix}
```
and the approximate solution becomes
```math
u(x) = 2\sin(x) -\frac28 \sin(2x) + \frac2{27} \sin(3x) -\frac2{64} \sin(4x) + \frac2{125} \sin(5x).
```
"""

# ╔═╡ dd44672f-b793-449a-8239-b5b2e1c5ef35
md"""
Generalising to arbitrary $m$ we find the solution coefficients
```math
c_i = \frac{f_i}{A_{ii}} = \frac{2\, (-1)^{i+1}}{i^3} 
\qquad \forall \, i = 1, \ldots, m
```
leading to the solution function
```math
\tag{27}
u(x) = \sum_{i=1}^m \frac{2\, (-1)^{i+1}}{i^3} \sin(i\, x),
```
which can be implemented as:
"""

# ╔═╡ a3705c5f-b833-4594-8ad3-c505f08660e7
function sine_solution(m, x)
	# Compute the function value of the numerical solution
	# for given m and x by evaluating (27)
	
	result = 0.0
	for i in 1:m
		result += 2 * (-1)^(i+1) / i^3 * sin(i * x)
	end
	result
end

# ╔═╡ d0a2133c-47d0-40d5-8508-18bb046d92d1
md"""
Choosing different values for $m$ we obtain:
- Show reference: $(@bind show_solution CheckBox(default=true))
- Show $m=10$:    $(@bind show_10 CheckBox(default=true))
"""

# ╔═╡ aa589aa5-1a5c-4782-af66-2509f21ac891
let
	x = range(0, π; length=1000)
	reference(x) = sine_solution(300, x)

	p = plot(title="Solution")
	if show_solution
		plot!(p, x, reference.(x); lw=2.5, label="Reference (m = 300)")
	end
	q = plot(; yaxis=:log, ylims=(1e-6, 1e-1), xlabel=L"x", title="Error")
	c = 1
	for m in [3, 5, 10]
		m == 10 && !show_10 && continue
		
		c += 1
		plot!(p, x, sine_solution.(m, x); lw=1.5, ls=:dash, label="Solution m=$m", c)
		
		error = [abs(sine_solution(m, xᵢ) - reference(xᵢ)) for xᵢ in x] .+ 1e-16
		plot!(q, x, error; label="", c, ls=:dash, lw=1.5)
	end

	plot(p, q; layout=grid(2, 1; heights=[0.75, 0.25]))
end

# ╔═╡ 4fee47fc-077d-4921-8572-7820ab84988e
md"""We see that already for 5 $\sin$ basis functions it becomes visually very hard to see the difference to the reference solution with $300$ basis functions.

Numerically we observe a **quadratic convergence** in a log-log plot:
"""

# ╔═╡ 1d55c98e-4ed9-4ac4-a847-b18432483472
let
	x = range(0, π; length=1000)
	reference(x) = sine_solution(800, x)

	ms = [5, 10, 20, 40, 50, 60, 80, 100, 300]
	errors = Float64[]
	for m in ms
		error = [abs(sine_solution(m, xᵢ) - reference(xᵢ)) for xᵢ in x]
		push!(errors, maximum(error))
	end
	
	p = plot(ms, errors, yaxis=:log, xaxis=:log, lw=2, mark=:o, title=L"Convergence in $m$", xlabel=L"m", ylabel="maximal error", label="Error sin Galerkin")
	plot!(p, ms, 1 ./ms.^2, ls=:dash, lw=2, label=L"O(m^2)")
	
	ylims!(p, (1e-6, 1e-1))
	xticks!(p, 10.0 .^ (0:0.5:3))
	yticks!(p, 10.0 .^ (0:-1:-6))

	p
end

# ╔═╡ 5f23627d-1d8e-4705-a9b8-155bb745cca9
md"""
## Comparison and conclusion
"""

# ╔═╡ 27fb5fde-566d-481b-b0d6-41c78c198eac
TODO("Some conclusion chapter to provide direct comparison between finite differences and Galerkin methods and from that some conclusion why Galerkin methods can be better (numerical stability, less unknowns for a targeted accuracy etc.)")

# ╔═╡ e88a64b8-1259-4320-9914-664bf30122a5
TODO("To introduce FEM we need the concept of a mass matrix in the Galerkin development above as the basis functions are not orthonormal.

	 See chapter 10.5 of Driscoll Brown
	 ")

# ╔═╡ 83fdd745-33e1-4e77-80a2-9008b0eb59b8
md"""
### Optional: Finite elements

A widely employed set of basis functions for Galerkin approximations
are the hat functions $φ_i = H_i$,
which we already discussed in the chapter
on [Interpolation (chapter 7)](https://teaching.matmat.org/numerical-analysis/07_Interpolation.html).
Recall, that given a set of nodes $x_0 < x_1 < \cdots < x_{n}$
the hat functions are defined as
```math
H_i(x) = \left\{ \begin{array}{ll}
\frac{x - x_{i-1}}{x_i - x_{i-1}} & \text{if $i>0$ and $x\in [x_{i-1}, x_i]$}\\
\frac{x_{i+1} - x}{x_{i+1} - x_{i}} & \text{if $i< n$ and $x\in [x_{i}, x_{i+1}]$}\\
0 & \text{otherwise}
\end{array}\right.
```
for $i = 0, \ldots, n$.
In this discussion we will again only consider equispaced nodes
for simplicity, i.e. we take the case with
$0 = x_0 < x_1 < \cdots < x_{n} = L$ and where $x_{i+1} - x_i = \frac{L}{n}$.
Defining $h = \frac{L}{n}$ the hat functions can thus equally be expressed as
```math
\tag{28}
H_i(x) = \left\{ \begin{array}{ll}
\frac{x - x_{i-1}}{\textcolor{red}{h}} & \text{if $i>0$ and $x\in [x_{i-1}, x_i]$}\\
\frac{x_{i+1} - x}{\textcolor{red}{h}} & \text{if $i<n$ and $x\in [x_{i}, x_{i+1}]$}\\
0 & \text{otherwise}
\end{array}\right.
```
for $i = 0, \ldots, n$.
"""

# ╔═╡ 57d511cc-aaae-42b3-80ec-7caf7d6f68ac
md"""
Due to the cardinality property of the hat functions $φ_i = H_{i}(x)$,
i.e. $H_i(x_j) = δ_{ij}$, the expansion (24) of the unknown solution function
can be simplified to
```math
\tag{29}
u(x) = \sum_{j=1}^{m} c_j φ_j(x) = \sum_{i=\textcolor{red}{1}}^{\textcolor{red}{n-1}} u_i\, H_i(x),
```
where $u_i = u(x_i)$, i.e. the numerical solution at the nodal points.
Note that the last sum in (29) is deliberately truncated to omit
the hat functions $H_0$ and $H_n$.
The nodal points of these two functions are $x_0 = 0$ and $x_n = L$,
respectively,
where these functions take by definition the value $1$.
As a result these functions are unable to satisfy
the condition $0 = φ_j(0) = φ_j(L)$ underlying $V_0$
and are thus excluded.
"""

# ╔═╡ 40c8b1a5-b61c-4af2-8b4b-f8b2b5e182dc
md"""
The importance of the hat functions for Galerkin methods stems from the fact that each hat function $H_i(x)$ is only non-zero in the interval $[x_{i-1}, x_{i+1}]$.
As a result when evaluating the Galerkin conditions (26) one never has to integrate over the entire domain $[0, L]$, but only a much smaller subset of this interval where the involved hat functions are non-zero
--- a noteworthy reduction of the computational effort.
"""

# ╔═╡ 4de7a48b-e825-4115-9269-edcbcaace19a
md"""
As an example consider the evaluation of the elements of the matrix $\mathbf{A} \in \mathbb{R}^{m\times m} = \mathbb{R}^{(n-1) \times (n-1)}$.
First note that the derivatives $H'_{i}$ are non-zero only where $H_{i}$
itself is non-zero, i.e. $[x_{i-1}, x_{i+1}]$.
Therefore the product $H'_{i}\, H'_{j}$
is only non-zero in the intersection
of $[x_{i-1}, x_{i+1}]$ and $[x_{j-1}, x_{j+1}]$.
"""

# ╔═╡ e72ea328-9ab3-4fc1-ad80-59091074c6f5
md"""
This intersection is empty, thus the matrix element $A_{ij}$ zero,
whenever $|i - j| > 2$, i.e. when the two indices are more than two
apart from each other.
For $|i - j| = 2$ the intersection contains only a single point,
which similarly implies $A_{ij} = 0$.
As a consequence
```math
A_{i,i+d} = \int_0^L k\,H_{i}'(x)\,H_{i+d}'(x)\,dx = 0
\qquad \text{if $d ≥ 2$}.
```
where $i = 1, \ldots, n-3$ and similarly by symmetry
```math
A_{i+d,i} = \int_0^L k\,H_{i+d}'(x)\,H_{i}'(x)\,dx = 0
\qquad \text{if $d ≥ 2$}.
```
again for $i = 1, \ldots, n-3$

"""

# ╔═╡ 153850e9-8b5c-45f2-a27f-0cf0b7389b21
md"""
Now we consider the remaining three cases, namely
```math
\begin{aligned}
A_{ii} &= \int_0^L k\,H_{i}'(x)\,H_{i}'(x)\,dx\\
	    &= \int_{x_{i-1}}^{x_{i}} k\,  H_{i}'(x)\,H_{i}'(x)\,dx
			+ \int_{x_{i}}^{x_{i+1}} k\,  H_{i}'(x)\,H_{i}'(x)\,dx \\
	    &= \int_{x_{i-1}}^{x_{i}} k\cdot  \frac{1}{h} \cdot\frac{1}{h}\,dx
			+ \int_{x_{i}}^{x_{i+1}} k\cdot \frac{-1}{h}\cdot\frac{-1}{h}\,dx \\ 
		&= \frac{k}{h^2} \left[ \int_{x_{i-1}}^{x_{i}} 1 \, dx
			+ \int_{x_{i}}^{x_{i+1}} 1 \, dx
		\right] \\
		&= \frac{2k}{h}
\end{aligned}
```
for $i=1,\ldots,n-1$ as well as
"""

# ╔═╡ 64977b68-8822-43e2-b2de-06b83b8a64e8
md"""
```math
\begin{aligned}
A_{i+1,i} =
A_{i,i+1} &= \int_0^L k\,H_{i}'(x)\,H_{i+1}'(x)\,dx\\
	    &= \int_{x_{i}}^{x_{i+1}} k\, H_{i}'(x)\,H_{i+1}'(x)\,dx\\
		&= \int_{x_{i}}^{x_{i+1}} k\cdot \frac{-1}{h} \cdot \,\frac{1}{h}\,dx\\
		&= -\frac{k}{h}.
\end{aligned}
```
for $i=1,\ldots,n-2$. We recover a tridiagonal matrix
```math
\tag{30}
\mathbf{A} = \frac{k}{h} \begin{pmatrix}
2  & -1 \\
-1 & 2 & -1\\
   & -1 & 2 & \ddots \\
   &  &  \ddots & \ddots & -1 \\
&&& -1 & 2
\end{pmatrix} \in \mathbb{R}^{(n-1)\times(n-1)}
```
"""

# ╔═╡ a30354e0-6254-46fd-a10b-f28d1c9dd6d4
TODO("Properly work the mass matrix stuff into this.")

# ╔═╡ 8af2b992-5072-47e0-85fc-dd7d332d505f
md"""
Mass matrix
```math
\begin{aligned}
M_{ii} &= \int_0^L H_{i}(x)\,H_{i}(x)\,dx\\
	    &= \int_{x_{i-1}}^{x_{i}} H_{i}(x)\,H_{i}(x)\,dx
			+ \int_{x_{i}}^{x_{i+1}} H_{i}(x)\,H_{i}(x)\,dx \\
	    &= \int_{x_{i-1}}^{x_{i}} \frac{(x-x_{i-1}) (x-x_{i-1})}{h^2} \,dx
			+ \int_{x_{i}}^{x_{i+1}} \frac{(x_{i+1}-x) (x_{i+1}-x)}{h^2} \,dx \\ 
		&= \frac{(x_i-x_{i-1})^3}{3h^2} + \frac{(x_{i+1}-x_i)^3}{3h^2}\\
		&= \frac{h^3}{3h^2} + \frac{h^3}{3h^2}\\
		&= \frac{2}{3}h
\end{aligned}
```
and
```math
\begin{aligned}
M_{i,i+1} &= \int_0^L H_{i}(x)\,H_{i+1}(x)\,dx
	    = \int_{x_{i}}^{x_{i+1}} H_{i}(x)\,H_{i+1}(x)\,dx \\
	    &= \int_{x_{i}}^{x_{i+1}} \frac{(x_{i+1}-x)(x-x_{i+1})}{h^2}\,dx\\ 
		&= - \frac{(x_{i+1} - x_i)^3}{3h^2}\\
		&= - \frac{h}{3}
\end{aligned}
```

"""

# ╔═╡ fade56d1-7270-4070-baf8-9b32c07a0b9d
md"""
Finally we consider the elements of the vector $\mathbf{f}$.
Since again $H_i$ is only non-zero in $[i-1, i+1]$,
we need to perform the two integrals
```math
f_i  = \int_{0}^{L} f(x)\,H_i(x)\,dx
= \int_{x_{i-1}}^{x_i} f(x)\,H_i(x)\,dx + \int_{x_i}^{x_{i+1}} f(x)\,H_i(x)\,dx
```
for $i = 1, \ldots, n-1$.
One can show that to a very good approximation one can replace
$f(x)$ by the average value of the function over the respective integrals, i.e.
```math
\begin{aligned}
\int_{x_{i-1}}^{x_i} f(x)\,H_i(x)\,dx &≈  \frac{f(x_{i-1}) + f(x_{i})}{2} \int_{x_{i-1}}^{x_i} H_i(x)\,dx = \frac{f(x_{i-1}) + f(x_{i})}{2} \cdot \frac{h}{2} \\
\int_{x_{i}}^{x_{i+1}} f(x)\,H_i(x)\,dx &≈  \frac{f(x_i) + f(x_{i+1})}{2}\int_{x_{i}}^{x_{i+1}}  H_i(x)\,dx = \frac{f(x_i) + f(x_{i+1})}{2} \cdot \frac{h}{2}.
\end{aligned}
```
In particular this approximation converges quadratically as $h\to0$,
which turns out to be exactly the order of approximation of the finite element method itself. Therefore putting in additional efforts to compute these integrals more accurately would not even improve the accuracy of our final result.
"""

# ╔═╡ 5a91ecbb-7ad9-4253-8dfb-f9afe2664f9b
md"""
Putting these developments tothere we obtain for $i=1,\ldots,n-1$ that
```math
\tag{31}
\begin{aligned}
f_i &= \frac{f(x_i) + f(x_{i-1})}{2} \frac{h}{2   }
+ \frac{f(x_{i+1}) + f(x_{i})}{2} \frac{h}{2}  \\
&= \frac{h}{4} \Big( f(x_{i-1}) + 2f(x_i) + f(x_{i+1}) \Big).
\end{aligned}
```
"""

# ╔═╡ 387244c2-2824-44cd-b413-45be7e447902
md"""
By solving the linear system
```math
(\mathbf{A} + \mathbf{M})\, \mathbf{u} = \mathbf{f}
```
we then obtain the coefficients $\mathbf{u} = (u_1, u_2, \ldots, u_{n-1})^T$
in the expansion (29).
Notably, due to the cardinality property of the hat functions $u_i = u(x_i)$, i.e. these coefficients are equal to evaluating our numerical approximation $u(x)$
at the nodal points $\{x_i\}_{i=1}^{n-1}$.


An implementation of this approch is given below:
"""

# ╔═╡ 563e2c1f-3c51-4be4-9693-5a2cfb8e6eda
function heat_equation_1d_fem(f, k, L, n)
	# f:  Function describing the external heat source
	# k:  thermal conductivity
	# L:  Length of the metal rod
	# N:  Number of nodal points for the finite element scheme
	
	h = L / n                 # Step size
	x = [i * h for i in 0:n]  # Nodal points (x₀, x₁, ..., xₙ)
	x_inner = x[2:n]          # Inner nodal points (x₁, ..., xₙ₋₁)
	

	# Build A and f: Note that these are a (n-1) × (n-1) matrix
	# respectively a vector of length (n-1)
	A = k/h * SymTridiagonal(2ones(n-1), -ones(n-2))
	M = h/3 * SymTridiagonal(2ones(n-1), -ones(n-2))
	f = h/4 * [f(x[i-1]) + 2f(x[i]) + f(x[i+1])
	           for i in 2:n]   # Skip first/last nodal point, i.e. x[1] = x₀

	u = (A + M) \ f
	(; x=x_inner, u)
end

# ╔═╡ d1115511-a798-47e7-aeea-fc87ad77d2f7
md"""
Let us return to the example we considered with the sine basis, i.e.
```math
\left\{
\begin{aligned}
- \frac{\partial^2 u}{\partial x^2}(x) &= x \qquad x \in (0, π)\\
u(0) &= u(L) = 0,
\end{aligned}
\right.
```
i.e. the Dirichlet heat equation (20) with $k=1$, $L=π$ and $f(x) = x$.

As a reference solution we consider the solution using the sine basis with $m=300$, which we know to be very accurate for this problem:
"""

# ╔═╡ 72b37ada-df37-4bfa-a2da-ae56026625dd
reference(x) = sine_solution(300, x)

# ╔═╡ 4042ff5d-3da6-418f-816b-905318e10f39
md"""
- `n = ` $(@bind n Slider(4:2:50; show_value=true, default=6))
"""

# ╔═╡ 1c8966d7-8a12-4b25-9c51-93dea0faeb14
let
	# Parameters of the problem
	f(x) = x
	k = 1
	L = π

	# Plot the reference
	p = plot(reference; lw=2.5, label="Reference", title="Solution", xlims=(-0.05, L + 0.05))

	# Compute FEM solution and plot it
	result = heat_equation_1d_fem(f, k, L, n)
	plot!(p, result.x, result.u; lw=1.5, ls=:dash, label="Solution n=$n", mark=:o)

	scatter!([0, L], [0, 0], label="Boundary values", mark=:o, c=3)
end

# ╔═╡ 25b07bb7-3ed3-4c7e-bbd6-d74c6a094127
let
	f(x) = x
	k = 1
	L = π

	ns = [5, 10, 20, 40, 50, 60, 80, 100, 200, 300, 500, 700]
	errors = Float64[]
	for n in ns
		result = heat_equation_1d_fem(f, k, L, n)
		error = abs.(result.u - reference.(result.x))
		push!(errors, maximum(error))
	end
	
	p = plot(ns, errors, yaxis=:log, xaxis=:log, lw=2, mark=:o, title=L"Convergence in $n$", xlabel=L"n", ylabel="maximal error", label="Error")
	plot!(p, ns, 1 ./ns.^2, ls=:dash, lw=2, label=L"O(n^2)")

	xticks!(p, 10.0 .^ (0:0.5:3))
	yticks!(p, 10.0 .^ (0:-1:-6))
end

# ╔═╡ 8fc6ab34-4c99-4764-a58e-8ad300f3a6ff
md"""
The finite element method is one of the most widely employed approaches in science and engineering to solve differential equations. With this short instroduction we only scratched the surface.

More information on the approach can be found for example in
[chapter 10.6](https://tobydriscoll.net/fnc-julia/bvp/galerkin.html) of Driscoll, Brown: *Fundamentals of Numerical Computation*.
"""

# ╔═╡ 9143e06a-9f27-44d1-ae2c-c50ce8dfdd30
md"Show outline $(@bind show_outline CheckBox(default=true))"

# ╔═╡ 3348115a-2fd4-4935-9dba-b85d9e0c95f5
if show_outline
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 400)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~1.0.0"
LaTeXStrings = "~1.4.0"
Plots = "~1.41.6"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "3a4816e2df6a778f82d0d0d5844c73951d782f6d"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cb7fe11da6adb8683cbacf8aa9b5237941e3a75"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.5+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "cac41ca6b2d399adfc95e51240566f8a60a80806"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.0+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "9e0fb9e54594c47f278d75063980e43066e26e20"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "67c6f1f085cb2671c93fe34244c9cccde30f7a26"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.5.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58e5ed5e386e156bd93e86b305ebd21ac63d2d04"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "86f5831495301b2a1387476cb30f86af7ab99194"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.0"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "850b06095ee71f0135d644ffd8a52850699581ed"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─b72b45ad-6191-40cb-9e9f-950bf1bfe212
# ╠═2cff355a-0339-11ef-085c-a5f9a6f5549e
# ╟─206ae56c-fcfa-4d6f-93e4-30f03dee8f90
# ╟─8c818dad-6789-4771-9eab-4e4ba78f5612
# ╟─cf1055db-a773-4dc2-bd66-0ae840c656a4
# ╟─93b28014-c767-4f40-b27c-4646f63faa89
# ╟─f1b662de-b4cf-4a73-9c00-580a2ea8a428
# ╟─b5c3ed4b-3c4d-40d5-9171-94187416b363
# ╟─226a30be-f721-4f3d-8d49-abddca22f661
# ╟─59f11469-c86e-4dc4-9eed-039101aaa058
# ╟─3e10cf8e-d5aa-4b3e-a7be-12ccdc2f3cf7
# ╟─23f9e709-1f77-4a96-bb09-2d358cdecc76
# ╟─6df8394f-6457-4db0-86e5-7a9dd7d65c9a
# ╟─b4b33be9-bb9a-4fe9-af48-8883c2fbc7e9
# ╠═b00a37d1-37ad-41b5-ba88-d6582ef680cb
# ╟─be847fe0-7447-493b-9087-f4a1d961edba
# ╟─a184ac2c-97ae-4221-a64e-215e9cadeb80
# ╟─3be4502b-95af-4ca8-8fd4-4c8cb659a858
# ╟─f44f3230-0062-49db-a488-f19910b58792
# ╠═393313f4-f2b0-4f2e-b175-ae24c0514e58
# ╟─e2ed03ad-f720-4e29-afdc-1c1c6bc233c0
# ╠═93bf2a7a-a926-4e0d-b127-889a0b82dcd2
# ╠═94b97d6c-6c45-4887-b15e-5731ae5bbe9c
# ╟─80df1d2f-24f9-4503-9075-eb5a37a486bd
# ╠═cdc1fc28-714b-45e1-ad51-31c6d70dcce3
# ╟─22876147-cac6-4c56-921a-8d86e605bdeb
# ╠═b0c8cf85-b034-41c6-a808-6562ace45768
# ╟─0e66bac1-a66a-47e9-aa22-b1f95ce2e45d
# ╠═796c699d-f6d3-40b3-aa48-fff0d627bb51
# ╟─78cf3006-47c8-48c2-be3a-2a518ff4731d
# ╠═ba98bfcb-a716-49e2-bc0b-7294bdaff012
# ╟─8d1cf059-05b1-44da-a51c-0157bb06efe9
# ╟─4f7bd62c-792c-4abb-ac93-28bb0169e739
# ╟─41326216-32f2-412f-b1a4-45c287f3b640
# ╠═51552b34-01b4-4ef0-866c-6ff6c33a3462
# ╠═5d361452-2eef-4503-a80f-151277d3b855
# ╠═8ddf3426-1d5d-42ee-b44a-3c4ea1ed5f4c
# ╟─b2418c32-854f-469c-ac9d-4fb7b93d5b79
# ╟─f1991ea0-9270-4c27-8032-4f4db84330e4
# ╠═27338f87-6f0a-4577-8cec-8eaba0c9d411
# ╟─b1d08025-3557-4045-9646-a42faae8b91c
# ╠═df13e0b3-ff18-4dc7-8d46-0a2f786cdbeb
# ╟─79d958dd-7b0e-428f-ba17-eb9983857504
# ╠═3bf1ee87-d42c-429d-8a07-f080f369f61f
# ╟─9f87dc91-c763-46b8-9277-2cc5f9a92cc0
# ╠═785b3c50-b008-4ec3-943e-d2e874942f7f
# ╟─fadf5e83-d5b4-48b9-911b-7f1309f02234
# ╠═452d1b69-42db-4ab0-889b-bc48ae0ec365
# ╟─f32bc6ac-f84d-453d-8cab-c06491410c24
# ╠═474dee17-627a-49b1-84e8-d692346d8aa0
# ╟─0bb8239c-22c0-4420-8ae7-59c5c00436a2
# ╟─69ca61cc-5482-45c1-99fc-0ef780771759
# ╟─814ac12b-5ed7-4512-b053-fbe729b90ce6
# ╠═7a647ffe-000e-4bb7-a348-033cd75b1369
# ╟─52961854-695e-44e7-b8e2-20ed7ff07cda
# ╟─359e61a3-d7fc-402f-ad6c-24b4301885d6
# ╟─f05ac977-45b4-4bd1-baf0-8a2d50525ec9
# ╟─b3a1222b-52cc-49ff-af74-cc472bd40999
# ╟─d9fb073d-042c-4724-834f-bcdb0cc4aec6
# ╟─75311119-2d3a-4ca4-bddb-5a595ffb09a4
# ╠═b482d8e2-901a-4d8b-85af-e0fbc5c46d4f
# ╠═2fd3f387-b950-462f-a06a-eabd458fe10e
# ╟─8039250a-c79a-4009-9a53-81307b8e0ace
# ╟─f3819d92-c61a-49b2-a77b-9837cbed0cae
# ╟─30b8bac2-7c95-4fee-bbde-b04a36258857
# ╠═858a60ba-8d0c-43d8-908d-0bb413c79285
# ╟─23307b35-c76e-487d-906a-18f65b691d79
# ╟─7ae614fb-f2e0-46af-a557-6a378f6553e8
# ╟─962ca311-c2ed-46eb-b043-827ecfb64ac3
# ╟─05999b01-4a41-480a-9f3a-2e04d97748b3
# ╟─7a86c1db-e5be-4310-951a-4de51bec4e98
# ╟─79816267-531b-4b99-bdea-1194574dd160
# ╟─bbb839d2-599f-40e2-a1f7-2419d07c96a8
# ╟─43ade623-1b58-4035-a1e5-803911063833
# ╟─e6f3450d-6fb7-4c59-a209-ae3587a43bab
# ╟─c7610b00-7277-46c3-881d-85d417db2e9e
# ╟─5e7387cc-f9eb-484e-ade7-c844bddeb716
# ╟─e32d65f1-9bf2-4043-811c-d47801d0f85e
# ╟─ef1ad23b-89bc-4b75-9560-4fa4fb96c802
# ╟─fcd48a38-eca2-44d9-a32e-9dd65c58fb89
# ╟─4d7ae336-be84-4960-835f-aa81ab91c4e2
# ╟─27a67eb3-76e5-4e23-a848-7a9e58485721
# ╟─52716cfc-6dea-44ad-aa96-842332fe607b
# ╟─dd44672f-b793-449a-8239-b5b2e1c5ef35
# ╠═a3705c5f-b833-4594-8ad3-c505f08660e7
# ╟─d0a2133c-47d0-40d5-8508-18bb046d92d1
# ╟─aa589aa5-1a5c-4782-af66-2509f21ac891
# ╟─4fee47fc-077d-4921-8572-7820ab84988e
# ╠═1d55c98e-4ed9-4ac4-a847-b18432483472
# ╟─5f23627d-1d8e-4705-a9b8-155bb745cca9
# ╠═27fb5fde-566d-481b-b0d6-41c78c198eac
# ╠═e88a64b8-1259-4320-9914-664bf30122a5
# ╟─83fdd745-33e1-4e77-80a2-9008b0eb59b8
# ╟─57d511cc-aaae-42b3-80ec-7caf7d6f68ac
# ╟─40c8b1a5-b61c-4af2-8b4b-f8b2b5e182dc
# ╟─4de7a48b-e825-4115-9269-edcbcaace19a
# ╟─e72ea328-9ab3-4fc1-ad80-59091074c6f5
# ╟─153850e9-8b5c-45f2-a27f-0cf0b7389b21
# ╟─64977b68-8822-43e2-b2de-06b83b8a64e8
# ╠═a30354e0-6254-46fd-a10b-f28d1c9dd6d4
# ╟─8af2b992-5072-47e0-85fc-dd7d332d505f
# ╟─fade56d1-7270-4070-baf8-9b32c07a0b9d
# ╟─5a91ecbb-7ad9-4253-8dfb-f9afe2664f9b
# ╟─387244c2-2824-44cd-b413-45be7e447902
# ╠═563e2c1f-3c51-4be4-9693-5a2cfb8e6eda
# ╟─d1115511-a798-47e7-aeea-fc87ad77d2f7
# ╠═72b37ada-df37-4bfa-a2da-ae56026625dd
# ╟─4042ff5d-3da6-418f-816b-905318e10f39
# ╠═1c8966d7-8a12-4b25-9c51-93dea0faeb14
# ╠═25b07bb7-3ed3-4c7e-bbd6-d74c6a094127
# ╟─8fc6ab34-4c99-4764-a58e-8ad300f3a6ff
# ╟─9143e06a-9f27-44d1-ae2c-c50ce8dfdd30
# ╟─3348115a-2fd4-4935-9dba-b85d9e0c95f5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
