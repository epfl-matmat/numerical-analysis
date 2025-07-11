### A Pluto.jl notebook ###
# v0.20.6

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
- The behaviour of (some) derivatives of $u$ and
- Some information about $u$ at the boundary $x=a$ and $x=b$
  for the initial time $t=0$. For example we may know
  the values $u(a, 0)$ and $u(b, 0)$.

Classic examples of boundary value problems include
- the [heat equation](https://en.wikipedia.org/wiki/Heat_equation),
  where $u(x, t)$ is the temperature at time $t$ and position $x$
  of a material, or
- [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation)
  where $u(x, t)$ is the electrostatic potential at position $x$ and time $t$
  corresponding to a given time-dependent density of electronic charges, or
- the [diffusion equation](https://en.wikipedia.org/wiki/Diffusion_equation)
  where $u(x, t)$ is the concentration of some solute at time $t$
  at position $x$ in the solution.

In this lecture, to simplify matters and to make this more concrete,
we will only consider one example, namely the **heat equation**.
"""

# ╔═╡ cf1055db-a773-4dc2-bd66-0ae840c656a4
md"""
## Stationary heat equation in one dimension
"""

# ╔═╡ 93b28014-c767-4f40-b27c-4646f63faa89
RobustLocalResource("https://raw.githubusercontent.com/epfl-matmat/numerical-analysis/30a76271c3fa3be9d24896ecd2450d18fa505df3/src/img/heat-equation.svg", "img/heat-equation.svg", :width => 500)

# ╔═╡ f1b662de-b4cf-4a73-9c00-580a2ea8a428
md"""
We consider a metal bar extending over the interval $[0, L]$. For a given moment in time $t$ we denote the **temperature** at a point $x \in [0, L]$ by $u(x, t)$. The **heat flux** we fruther denote by $\phi(x, t)$ and allow for an **external source of heat** $f(x)$, e.g. by one more more heat guns pointed at the metal bar.
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
c ρ \, \frac{\partial u}{\partial t}(x, t) = k \frac{\partial^2 u}{\partial x^2}(x, t) + f(x) \qquad x \in (0, L), \ t > 0
```
which is the **heat equation** in one dimension,
a **partial differential equation**.
The heat equation describes the **evolution of the temperature in a system**
(here the metal bar) and has as its unknown the temperature distribution $u(x, t)$ at any location at any point in time.
"""

# ╔═╡ 59f11469-c86e-4dc4-9eed-039101aaa058
md"""
To simplify matters in this lecture 
**we restrict** ourselves to the case of the **stationary heat equation**, i.e. the case where we assume the **temperature $u(x, t)$** to be in a state where it **no longer changes with time**, i.e. $\frac{\partial u}{\partial t} = 0$.
Dropping thus the zero terms and the dependency on $t$ we obtain the **stationary heat equation**
```math
-k\, \frac{\partial^2 u}{\partial x^2}(x) = f(x) \qquad x \in (0, L).
```
"""

# ╔═╡ 3e10cf8e-d5aa-4b3e-a7be-12ccdc2f3cf7
md"""
To fully specify the problem and find a unique solution we still **need to indicate what happens at the extremities of the bar**, i.e. at $u(0)$ and at $u(L)$.
Since these two points sit at the boundary of our computational domain $[0, L]$ we usually call these conditions on $u$ **boundary conditions**. 
"""

# ╔═╡ 7fd851e6-3180-4008-a4c0-0e08edae9954
md"""
We sketch the two most common cases:

- Consider the case where the bar is in contact with heat baths of constant temperature. A temperature $b_0$ for the left-hand side bath and $b_L$ for the right-hand side bath. In this case we know $u(0) = b_0$ and $u(L) = b_L$. We obtain the **heat equation using Dirichlet boundary conditions**:
  ```math
  \tag{1}
  \left\{
  \begin{aligned}
  -k \, \frac{\partial^2 u}{\partial x^2}(x) &=  f(x) \qquad x \in (0, L).\\
  u(0) &= b_0\\
  u(L) &= b_L
  \end{aligned}
  \right.
  ```
- The bar may also be insulating at the boundary, such that no heat may be allowed to flow in or out of the extremities. In this case the boundary conditions are $0 = -k \frac{\partial u}{\partial x}(0)$ and $0 = -k \frac{\partial u}{\partial x}(L)$ (zero heat flux).
  We obtain a **heat equation with Neumann boundary conditions**, where the flux at the boundary is controlled:
  ```math
  \tag{2}
  \left\{
  \begin{aligned}
  -k \, \frac{\partial^2 u}{\partial x^2}(x) &=  f(x) \qquad x \in (0, L).\\
  \frac{\partial u}{\partial x}(0) &= - \frac{\phi_0}{k} \\
  \frac{\partial u}{\partial x}(L) &= - \frac{\phi_L}{k}
  \end{aligned}
  \right.
  ```
  where in the insulating case $\phi_0 = \phi_L = 0$. We will not consider Neumann problems any further here.

Since for such kind of problems knowing the boundary is imperative to obtain a unique solution one refers to such problems as **boundary value problems**.
"""

# ╔═╡ 52c7ce42-152d-40fd-a910-78f755fcae47
md"""
## Finite difference approximation

We focus on the case of a Dirichlet boundary (1) with $k = 1$.
Our goal is thus to find a function $u : [0, L] \to \mathbb{R}$ with
```math
\left\{
\begin{aligned}
- \frac{\partial^2 u}{\partial x^2}(x) &=  f(x) \qquad x \in (0, L)\\
u(0) &= b_0, \quad u(L) = b_L,
\end{aligned}
\right.
```
were $b_0, b_L \in \mathbb{R}$.

TODO IVP is now after BCP, adjust reference accordingly

Similar to our approach when [solving initial value problems (chapter 12)](https://teaching.matmat.org/numerical-analysis/12_Initial_value_problems.html)
we **divide the full interval $[0, L]$ into $N+1$ subintervals** $[x_j, x_{j+1}]$
of uniform size $h$, i.e. 
```math
x_j = j\, h \quad j = 0, \ldots, N+1, \qquad h = \frac{L}{N+1}.
```
Our goal is thus to find approximate points $u_j$ such that $u_j ≈ u(x_j)$ at the nodes $x_j$.
"""

# ╔═╡ 82788dfd-3462-4f8e-b0c8-9e196dac23a9
md"""
Due to the Dirichlet boundary conditions $u(0) = b_0$ and $u(L) = b_L$.
To satisfy these we neccessarily need $u_0 = b_0$ and $u_{N+1} = b_L$.
The unknowns of our problems are therefore those
$u_j$ with $1 ≤ j ≤ N$
--- the **internal nodes** of the interval.
"""

# ╔═╡ d43ecff3-89a3-4edd-95c2-7262e317ce29
md"""
These internal nodes $u(x_j)$ need to satisfy
```math
- \frac{\partial^2 u}{\partial x^2}(x_j) = f(x_j) \qquad \forall\, 1 ≤ j ≤ N.
```
As the derivatives of $u$ are unknown to us we employ a
**[central finite-difference formula](https://teaching.matmat.org/numerical-analysis/09_Numerical_differentiation.html)**
to replace this derivative by the approximation
```math
\tag{3}
f(x_j) = -\frac{\partial^2 u}{\partial x^2}(x_j) ≈ -\frac{u(x_{j-1}) - 2u(x_j) + u(x_{j+1})}{h^2} \qquad \forall\, 1 ≤ j ≤ N
```
or in terms of $u_j$:
```math
\tag{4}
\frac{-u_{j-1} + 2u_j - u_{j+1}}{h^2} = f(x_j) \qquad \forall\, 1 ≤ j ≤ N.
```
"""

# ╔═╡ 1fb53091-89c8-4f70-ab4b-ca2371b830b2
md"""
Due to our **imposed boundary conditions** $u_0 = b_0$ and $u_{N+1} = b_L$,
such that we can simplify the equations for $j = 1$ and $j=N$ as follows:
```math
\begin{aligned}
j&=1: & \frac{-b_0 + 2u_1 - u_2}{h^2} &= f(x_1) \quad &\Leftrightarrow \quad \frac{2u_1 - u_2}{h^2} &= f(x_1) + \frac{b_0}{h^2}\\
j&=N: & \frac{-u_{N-1} + 2u_N - b_L}{h^2} &= f(x_N) \quad &\Leftrightarrow \quad \frac{-u_{N-1} +2u_N}{h^2} &= f(x_N) + \frac{b_L}{h^2}\\
\end{aligned}
```
Collecting everything we obtain the system of equations
```math
\tag{5}
\left\{
\begin{aligned}
\frac{2u_1 - u_2}{h^2} &= f(x_1) + \frac{b_0}{h^2} && j=1\\
\frac{-u_{j-1} + 2u_j - u_{j+1}}{h^2} &= f(x_j)  && \forall\, 2 ≤ j ≤ N-1 \\
\frac{-u_{N-1} +2u_N}{h^2} &= f(x_N) + \frac{b_L}{h^2} && j=N
\end{aligned}
\right.
```
"""

# ╔═╡ 129d2e86-49e6-4c5f-aca1-0cb45d923294
Foldable("Optional: Does the problem (5) always have a solution ?",
md"""
An important question in the mathematical literature for problems such as (5)
is whether this problem is **well-posed**, that is whether a unique solution to such problems even exists at all.

For the boundary value problems we consider in this chapter it turns out that such a solution always exists, provided that the incoming heat $f$ and the boundary values $b_0$ and $b_L$ are finite. More precisely we have the following result:

!!! info "Theorem 1"
	For every $\mathbf{f} = (f(x_1), f(x_2), \ldots f(x_n))^T \in \mathbb{R}^N$
	and all $b_0, b_L \in \mathbb{R}$, the system of equations shown in (5)
	admits a unique solution
	$\mathbf{u} = (u_1, \ldots, u_N)^T \in \mathbb{R}^N$ with
	```math
	\tag{6}
	\max_{j=1,\ldots,N} |u_j| ≤ \frac18 \max_{j=1,\ldots,N} |f(x_j)| + \max(|b_0|, |b_L|).
	```

Colloquially speaking this Theorem ensures that the solution $\mathbf{u}$ cannot take too large values and moreover it tells us that the solution values are always controlled by the norm of the vector $\mathbf{f}$ and the boundary values $b_0$ and $b_L$.
""")

# ╔═╡ cbcd1b1e-7755-4177-b52a-f361ee025e01
md"""
Finally, to write problem (5) more compactly we introduce
a **vector of all unknowns** $\mathbf{u} = (u_1, \ldots, u_N)^T \in \mathbb{R}^N$,  define the **system matrix** $\mathbf{A} \in \mathbb{R}^{N\times N}$
as well as the **right-hand side** $\mathbf{b} \in \mathbb{R}^N$ with
```math
\tag{7}
\mathbf{A} = \frac{1}{h^2} \begin{pmatrix}
2  & -1 \\
-1 & 2 & -1\\
   & -1 & 2 & \ddots \\
   &  &  \ddots & \ddots & -1 \\
&&& -1 & 2
\end{pmatrix}
\qquad\qquad
\mathbf{b} = \begin{pmatrix}
f(x_1) + \frac{b_0}{h^2} \\ f(x_2) \\ \vdots \\ f(x_{N-1}) \\ f(x_N) + \frac{b_L}{h^2}
\end{pmatrix}.
```
With these objects the problem can be written as a linear system
```math
\tag{8}
\mathbf{A} \mathbf{u} = \mathbf{b}.
```
which is to be solved for the unknows $\mathbf{u}$.
"""

# ╔═╡ c21502ce-777f-491a-a536-ff499fc172fc
md"""
We notice that $\mathbf{A}$ is **symmetric and tridiagonal**. Additionally one can show $\mathbf{A}$ to be **positive definite**.
Problem (8) can therefore be **efficiently solved** using [direct methods based on (sparse) LU factorisation (chapter 5)](https://teaching.matmat.org/numerical-analysis/05_Direct_methods.html) or an [iterative approaches (chapter 6)](https://teaching.matmat.org/numerical-analysis/06_Iterative_methods.html), e.g. the conjugate gradient method.
"""

# ╔═╡ c2bb42b3-4fee-4ad4-84c0-06f58c7f7665
md"""
We summarise:

!!! info "Algorithm 1: Finite differences for Dirichlet bounary value problems"
	Given a boundary value problem with Dirichlet boundary conditions
	```math
	\tag{9}
	\left\{
	\begin{aligned}
	- \frac{\partial^2 u}{\partial x^2}(x) &=  f(x) \qquad x \in (0, L)\\
	u(0) &= b_0, \quad u(L) = b_L,
	\end{aligned}
	\right.
	```
	and nodes $x_j = j\, h$ with $j = 0, \ldots, N$ and $h = \frac{L}{N+1}$
	solve the linear system $\mathbf{A} \mathbf{u} = \mathbf{b}$
	with $\mathbf{A}$ and $\mathbf{b}$ given in (7)
	for $\mathbf{u} = (u_1, \ldots, u_N)^T$ and $u_0 = b_0$, $u_{N+1} = b_L$.

	The values $u_0, u_1, \ldots, u_{N+1}$ approximate the solution $u(x)$
	at the nodal points $\{x_j\}_{j=0}^{N+1}$.
"""

# ╔═╡ e38a1d5b-a35c-4942-aeb2-0e91668ac8f3
md"""
One implementation of Algorithm 1, which solves the linear system using LU factorisation is:
"""

# ╔═╡ 3bc63106-3d1c-4b4e-a6a6-682e09475cb4
function fd_dirichlet(f, L, b₀, bₗ, N)
	# f:  Function describing the external heat source
	# L:  Length of the metal rod
	# b₀: Left-hand side boundary value
	# bₗ: Right-hand side boundary value
	# N:  Number of nodal points for the finite-differences scheme
	
	h = L / (N+1)             # Step size
	x = [i * h for i in 1:N]  # Nodal points

	# Build system matrix A
	diagonal    = 2ones(N)   / h^2
	offdiagonal = -ones(N-1) / h^2
	A = SymTridiagonal(diagonal, offdiagonal)

	# Build right-hand side
	b = [f(xⱼ) for xⱼ in x]  # Evaluate function f at nodal points

	# Set terms due to boundary conditions
	b[1] = f(x[1]) + b₀ / h^2
	b[N] = f(x[N]) + bₗ / h^2

	# Solve problem using LU factorisation
	u = A \ b

	(; u, x, h, A, b)  # Return intermediates and results
end

# ╔═╡ 4c9e3775-bb3e-4ce9-b18e-251aac342b5d
md"""
Let us consider the example
```math
\left\{
\begin{aligned}
- \frac{\partial^2 u}{\partial x^2}(x) &= \sin(x) \qquad x \in (0, 2π)\\
u(0) &= 1, \quad u(L) = 2,
\end{aligned}
\right.
```
i.e. the case of 
"""

# ╔═╡ 9f4df9d4-e0cb-4706-8158-0da7a6505442
begin
	b₀ = 1
	bₗ = 2
	L = 2π
	f(x) = sin(x)
end;

# ╔═╡ 9f87dc91-c763-46b8-9277-2cc5f9a92cc0
md"""
In this case the exact solution can be found analytically as
"""

# ╔═╡ 785b3c50-b008-4ec3-943e-d2e874942f7f
u_exact(x) = sin(x) + x/2π + 1

# ╔═╡ e49d8cb7-dcf6-4449-9ecc-ab49f3c1bf28
md"""
Let's solve it with finite differences
"""

# ╔═╡ 0dbcfaa4-3e0e-470b-8412-4eb0a471bd40
md"and plot the result:"

# ╔═╡ 0bb8239c-22c0-4420-8ae7-59c5c00436a2
md"`N = ` $(@bind N Slider(10:5:100; default=10, show_value=true))"

# ╔═╡ b0342965-0c56-4ae8-a092-9f4ac565d249
res_fd = fd_dirichlet(f, L, b₀, bₗ, N);

# ╔═╡ 474dee17-627a-49b1-84e8-d692346d8aa0
let
	plot(u_exact; xlims=(-0.1, L+0.1), legend=:topright,
	     label="reference", lw=2, xlabel=L"x", ylabel=L"u(t)",
	     title=L"-\frac{\partial^2 u}{\partial x^2} = \sin(x)")

	plot!(res_fd.x, res_fd.u; label="Finite differences N=$N", mark=:o, lw=2, ls=:dash)

	scatter!([0, L], [b₀, bₗ], label="Boundary values", mark=:o, c=3)
end

# ╔═╡ 814ac12b-5ed7-4512-b053-fbe729b90ce6
md"""
Note that the plot of the finite-difference solution does not include the nodal point at the boundaries ($x = 0$ and $x = 2π$ --- shown in green in the plot).


We observe that as the number of points `N` is increased,
the solution visually becomes more and more accurate.
The output of our numerical procedure are 
the computed points $u_1, u_2, \ldots u_N$.
These points are meant to approximate the true function values 
$u(x_1), u(x_2), \ldots, u(x_N)$ of $u$ evaluated at the nodal points.
A measure for the accuracy of our obtained solution is thus
the maximal deviation any of these computed points has from
from the true $u(x_j)$.
This measure is called the **global error** defined by
```math
|e_j| = |u(x_j) - u_j| \qquad \text{for $j = 1, \ldots, N$}.
```

Numerically, by observing how this global error changes
as $N$ increases, we observe **quadratic convergence**:
"""

# ╔═╡ 7a647ffe-000e-4bb7-a348-033cd75b1369
let
	Ns = [ round(Int, 5 * 10^k) for k in 0:0.5:4 ]
	errors = Float64[]
	for N in Ns
		res_fd = fd_dirichlet(f, L, b₀, bₗ, N)
		x = res_fd.x
		u = res_fd.u
		error = [ u_exact(x[j]) - u[j] for j in 1:length(x)]
		push!(errors, maximum(error))		
	end
	
	p = plot(; title="Convergence of Dirichlet finite differences",
			   xaxis=:log, yaxis=:log,
	           xlabel=L"N", ylabel="Largest global error",
	           legend=:bottomleft,
	           xlims=(1, 3e5), ylims=(1e-13, 10))
	plot!(p, Ns, errors; lw=2, mark=:o, c=1, label="Error")
	plot!(p, Ns, 0.02(first(Ns) ./ Ns).^2, ls=:dash, lw=2, label=L"O(n^{-2})", c=1)

	xticks!(p, 10.0 .^ (0:1:5))
	yticks!(p, 10.0 .^ (0:-2:-12))
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

	For the **central finite difference scheme** discussed here $p=2$.
"""

# ╔═╡ d9fb073d-042c-4724-834f-bcdb0cc4aec6
md"""
### Numerical stability

Let us push the convergence plot from above a little further:
"""

# ╔═╡ 2fd3f387-b950-462f-a06a-eabd458fe10e
let
	Ns = [ round(Int, 5 * 10^k) for k in 2:0.1:5.2 ]
	errors = Float64[]
	for N in Ns
		res_fd = fd_dirichlet(f, L, b₀, bₗ, N)
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
	condition_numbers = [cond(fd_dirichlet(f, L, b₀, bₗ, N).A) for N in Ns]

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

# ╔═╡ e35398a8-d8e5-4b19-b5a9-2c7ab6b497e9
md"""
### Optional: Error analysis of the discretisation error

Let us come back to the observed quadratic convergence
in the regime of small $N$ where the discretisation error dominates.
Ignoring thus the floating-point error we conclude that the global error 
```math
|e_j| = |u(x_j) - u_j| \qquad \text{for $j = 1, \ldots, N$}.
```
should scale as $\max_{j=1,\ldots,N} |u(x_j) - u_j| ≤ α\, h^2$ for some constant $α$.
In this section we will make this more quantitative.
"""

# ╔═╡ e70af7c2-6df6-4724-8a82-7f3fc5c1fadf
md"""
First, since the value $u_j$ at each nodal point is determined by solving (4),
i.e. by replacing the *exact* function values $u(x_{j-1})$, $u(x_j)$ and $u_{j+1}$
by the approximations $u_{j-1}$, $u_{j}$ and $u_{j+1}$,
there are in fact two contributions to the global error:
1. The error due to employing the finite difference formula in (3)
   instead of the exact partial derivative $\frac{\partial^2 u}{\partial x^2}$.
2. The propagation of error from the neighbours $u_{j-1}$ / $u_{j+1}$ to $u_j$ it self and vice versa.
"""

# ╔═╡ 7fc0e2c7-fdf5-458a-9c54-63f906cb63ad
md"""
We start by understanding the first contribution,
which is the **local error** due to employing the finite difference formula (3).

!!! info "Definition: Local truncation error"
	Given the exact solution $u(x_j)$ of problem (9) evaluated at the nodal points $x_j = j\,h$ for $j = 0, \ldots, N+1$ of the finite difference scheme of Algorithm 1, we define the **local truncation error** at the node $x_j$ as
	```math
	\tag{10}
	\tau^h_j = \frac{-u(x_{j-1}) + 2u(x_j) - u(x_{j+1})}{h^2} - f(x_j)
	\qquad j = 1, \ldots, N.
	```
	That is the local truncation error is the residual of the numerical scheme (4)
	when we replace the approximated solution $u_j$ by the exact solution $u(x_j)$.
"""

# ╔═╡ 1f3b6898-3ff1-40b5-bcc2-97db023fc0e3
md"""
In our case one can show for the local truncation error:

!!! info "Lemma 2: Bound on the local truncation error"
	If the solution $u$ is four times differentiable the local truncation error (10) satisfies
	```math
	\tag{11}
	\max_{j=1,\ldots,N} |τ^h_j| ≤ \frac{h^2}{12} \max_{x\in[0,L]}|u''''(x)|
	= \frac{h^2}{12} \|u''''\|_\infty.
	```
"""

# ╔═╡ 5fa3718d-d58c-492f-9e04-1bca93d8b9fd
md"""
> **Proof sketch:**
> Since $u$ is the exact solution to (9) we have that $f(x_j) = -\frac{\partial^2 u}{\partial x^2}(x_j)$ for all $j = 1, \ldots, N$. Therefore
> ```math
> \begin{aligned}
> \tau^h_j &= \frac{-u(x_{j-1}) + 2u(x_j) - u(x_{j+1})}{h^2} - f(x_j)\\
> &= \frac{-u(x_{j-1}) + 2u(x_j) - u(x_{j+1})}{h^2} + \frac{\partial^2 u}{\partial x^2}(x_j) \\
> &= -\frac{u(x_{j} - h) - 2u(x_j) + u(x_{j}+h)}{h^2} + u''(x_j)
> \end{aligned}
> ```
> Considering a Taylor expansion to fourth order of $u(x_j \pm h)$
> around $x_j$ than leads to the result.
"""

# ╔═╡ b811d679-e80b-41e1-8e85-dd81d01f00a7
md"""
Next we tackle the second aspect,
namely how this local truncation error propagates to give global error.

By rearranging the definition of the local truncation error (10),
we find that the exact solution satisfies the discrete problem
```math
\tag{12}
\frac{-u(x_{j-1}) + 2u(x_j) - u(x_{j+1})}{h^2} = f(x_j) + \tau^h_j 	\qquad j = 1, \ldots, N.
```
In contrast according to (4) the approximate solution satisfies
```math
\tag{4}
\frac{-u_{j-1} + 2u_j - u_{j+1}}{h^2} = f(x_j) 	\qquad j = 1, \ldots, N.
```
Subtracting (4) from (12) thus leads to
```math
\frac{-e_{j-1} + 2e_j - e_{j+1}}{h^2}
= \tau^h_j 	\qquad j = 1, \ldots, N,
```
where we used the definition of the error $e_j = u(x_j) - u_j$.
At the boundary, where we set $u_0 = b_0 = u(0)$ and $u_{N+1} = b_L = u(L)$
we have $e_0 = e_{N+1} = 0$, such that the error satisfies
```math
\tag{13}
\left\{
\begin{aligned}
\frac{2e_1 - e_2}{h^2} &= τ^h_1  && j=1\\
\frac{-e_{j-1} + 2e_j - e_{j+1}}{h^2} &= τ^h_j  && \forall\, 2 ≤ j ≤ N-1 \\
\frac{-e_{N-1} +2e_N}{h^2} &= τ^h_N && j=N,
\end{aligned}
\right.
```
which is again a discretised Dirichlet bounary value problem (5) with the conditions that $f(x_j) = τ^h_j$ and $b_0=b_L=0$.

We can therefore apply equation (14) of Theorem 1 
(see the folded section *Optional: Does the problem (5) always have a solution ?*)
to problem (13) and conclude
```math
\tag{15}
\max_{j=1,\ldots,N} |e_j| ≤ \frac{1}{8} \max_{j=1,\ldots,N} |τ^h_j|,
```
which relates the local truncation error to the error at each nodal point.
"""

# ╔═╡ f472502a-678f-479a-ab18-9ba0af7cde2f
md"""
Combining this result with Lemma 2 yields

!!! info "Theorem 3: Global error of finite differences"
	Under the assumption that the exact solution $u$ of the Dirichlet boundary value problem (9) is four times differentiable,
	the finite differences scheme of Algorithm 1 converges as
	```math
	\tag{16}
	\max_{j=1,\ldots,N} |u(x_j) - u_j| ≤ α \, h^2
	```
	where $α = \frac{1}{96} \|u''''\|_\infty$.
	It thus achieves **quadratic convergence**.
"""

# ╔═╡ b3d4248c-6f83-4594-a41a-d6c3c7e8d405
TODO("Show application of finite difference method to a second example Boundary Value Problem.

	 A good example could be the Allan-Cahn equation (Demo 10.5.5 in Driscoll Brown)
	 
")

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

# ╔═╡ 3348115a-2fd4-4935-9dba-b85d9e0c95f5
let
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
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.4"
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "386992f2ee9496aea372b754edcbbad6c546e5b7"

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
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "062c5e1a5bf6ada13db96a4ae4749a4c2234f521"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.9"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

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
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

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
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

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
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "0ff136326605f8e06e9bcf085a356ab312eef18a"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.13"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "9cb62849057df859575fc1dda1e91b82f8609709"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.13+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

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
git-tree-sha1 = "f93655dc73d7a0b4a368e3c0bce296ae035ad76e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.16"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "1d4015b1eb6dc3be7e6c400fbd8042fe825a6bac"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.10"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "872cd273cb995ed873c58f196659e32f11f31543"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.44"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd10d2cc78d34c0e2a3a36420ab607b611debfbb"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.7"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

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
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "688d6d9e098109051ae33d126fcfc88c4ce4a021"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

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
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a9697f1d06cc3eb3fb3ad49cc67f2cfabaac31ea"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.16+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3b31172c032a1def20c98dae3f2cdc9d10e3b561"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
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
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "41c9a70abc1ff7296873adc5d768bff33a481652"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.12"

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

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "8252b5de1f81dc103eb0293523ddf917695adea1"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "d3de2694b52a01ce61a036f18ea9c0f61c4a9230"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.62"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "5cf59106f9b47014c58c5053a1ce09c0a2e0333c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.7.3"

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

    [deps.Revise.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

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
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

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
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

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
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

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

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

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
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "85c7811eddec9e7f22615371c3cc81a504c508ee"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+2"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5db3e9d307d32baba7067b13fc7b5aa6edd4a19a"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.36.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

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
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+3"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "807c226eaf3651e7b2c468f687ac788291f9a89b"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.3+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d7155fea91a4123ef59f42c4afb5ab3b4ca95058"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+3"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "6fcc21d5aea1a0b7cce6cab3e62246abd1949b86"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.0+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "984b313b049c89739075b8e2a94407076de17449"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.2+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a1a7eaf6c3b5b05cb903e35e8372049b107ac729"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.5+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "b6f664b7b2f6a39689d822a6300b14df4668f0f4"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.4+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a490c6212a0e90d2d55111ac956f7c4fa9c277a6"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+1"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "dbc53e4cf7701c6c7047c51e17d6e64df55dca94"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+1"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "ab2221d309eda71020cdda67a973aa582aa85d69"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+1"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6e50f145003024df4f5cb96c7fce79466741d601"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.56.3+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0ba42241cb6809f1a278d0bcb976e0483c3f1f2d"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+1"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "068dfe202b0a05b8332f1e8e6b4080684b9c7700"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.47+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "63406453ed9b33a0df95d570816d5366c92b7809"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+2"
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
# ╟─7fd851e6-3180-4008-a4c0-0e08edae9954
# ╟─52c7ce42-152d-40fd-a910-78f755fcae47
# ╟─82788dfd-3462-4f8e-b0c8-9e196dac23a9
# ╟─d43ecff3-89a3-4edd-95c2-7262e317ce29
# ╟─1fb53091-89c8-4f70-ab4b-ca2371b830b2
# ╟─129d2e86-49e6-4c5f-aca1-0cb45d923294
# ╟─cbcd1b1e-7755-4177-b52a-f361ee025e01
# ╟─c21502ce-777f-491a-a536-ff499fc172fc
# ╟─c2bb42b3-4fee-4ad4-84c0-06f58c7f7665
# ╟─e38a1d5b-a35c-4942-aeb2-0e91668ac8f3
# ╠═3bc63106-3d1c-4b4e-a6a6-682e09475cb4
# ╟─4c9e3775-bb3e-4ce9-b18e-251aac342b5d
# ╠═9f4df9d4-e0cb-4706-8158-0da7a6505442
# ╟─9f87dc91-c763-46b8-9277-2cc5f9a92cc0
# ╠═785b3c50-b008-4ec3-943e-d2e874942f7f
# ╟─e49d8cb7-dcf6-4449-9ecc-ab49f3c1bf28
# ╠═b0342965-0c56-4ae8-a092-9f4ac565d249
# ╟─0dbcfaa4-3e0e-470b-8412-4eb0a471bd40
# ╠═474dee17-627a-49b1-84e8-d692346d8aa0
# ╟─0bb8239c-22c0-4420-8ae7-59c5c00436a2
# ╟─814ac12b-5ed7-4512-b053-fbe729b90ce6
# ╠═7a647ffe-000e-4bb7-a348-033cd75b1369
# ╟─52961854-695e-44e7-b8e2-20ed7ff07cda
# ╟─d9fb073d-042c-4724-834f-bcdb0cc4aec6
# ╟─2fd3f387-b950-462f-a06a-eabd458fe10e
# ╟─8039250a-c79a-4009-9a53-81307b8e0ace
# ╟─f3819d92-c61a-49b2-a77b-9837cbed0cae
# ╟─30b8bac2-7c95-4fee-bbde-b04a36258857
# ╟─e35398a8-d8e5-4b19-b5a9-2c7ab6b497e9
# ╟─e70af7c2-6df6-4724-8a82-7f3fc5c1fadf
# ╟─7fc0e2c7-fdf5-458a-9c54-63f906cb63ad
# ╟─1f3b6898-3ff1-40b5-bcc2-97db023fc0e3
# ╟─5fa3718d-d58c-492f-9e04-1bca93d8b9fd
# ╟─b811d679-e80b-41e1-8e85-dd81d01f00a7
# ╟─f472502a-678f-479a-ab18-9ba0af7cde2f
# ╠═b3d4248c-6f83-4594-a41a-d6c3c7e8d405
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
# ╟─3348115a-2fd4-4935-9dba-b85d9e0c95f5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
