### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 8197b6f0-00b7-11ef-2142-4b7cbbaefd90
begin
	using Plots
	using PlutoUI
	using PlutoTeachingTools
	using LaTeXStrings
	using DifferentialEquations
	using ForwardDiff
	using HypertextLiteral
	using LinearAlgebra
end

# ╔═╡ ba9b6172-0234-442c-baaa-876b12f689bd
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/11_Initial_value_problems.pdf)
"""

# ╔═╡ d8406b01-e36f-4953-a5af-cd563005c2a1
TableOfContents()

# ╔═╡ 692aa6e2-21a9-4dad-a719-9d8ad88bd467
md"""
# Initial value problems
"""

# ╔═╡ ebabe949-5492-4b8b-959d-05c5728da043
md"""
In computational science we are often faced with quantitites that change continuously in space or time. For example the temperature profile of a hot body or the population of an animal species over time. Such problems are often modelled by differential equations. If there is only a single indepentent variable $t$ we call the model an **ordinary differential equation**. The usual setup is that for initial value $t = 0$ we have some knowledge about our problem, e.g. by performing a measurement, and the key question is thus how the situation evolves for $t > 0$. 
"""

# ╔═╡ 05ef8174-e9a6-4280-8640-08d74635fba2
md"""
!!! warning "Example: Ducks on the Leman"
	Suppose we want to model the population of ducks on the Leman, i.e. let $u(t)$ denote the number of ducks on the lake at time $t$. To simplify the mdoelling we  allow for "fractional ducks", i.e. we allow $u$ to be any real number.

	First we assume a constant growth rate $C$ for the number of ducks at any time $t$, i.e. that the number of ducks born minus the number of ducks deceased in one unit of time is proportional to $u(t)$ --- the current population of ducks. This leads to the ordinary differential equation (ODE)
	```math
		\left\{\begin{aligned}
		\frac{d u(t)}{d t} &= C\, u(t) &&t > 0\\
		u(0) = u_0
		\end{aligned}\right.
	```
	where $u_0$ is the initial number of ducks we observed at $t = 0$.
    In comparison to the general definition shown above
    we thus have $f(t, u) = C\, u$.

	The solution of this linear equation is
	```math
	u(t) = u_0 \, e^{Ct},
	```
	i.e. exponential growth.

	Clearly as the size of the Leman is finite, this is not a reasonable model
	as both food and space on the lake is limited, which should cap the growth.
	To improve the model we assume that the death rate itself is proportional
	to the size of the population, i.e. $d\,u$ for a constant $d > 0$.
	Keeping the idea of a constant birth rate $b > 0$ one thus obtains
	an improved model
	```math
		\left\{\begin{aligned}
		\frac{d u(t)}{d t} &= b\, u(t) - d \, \left(u(t) \right)^2 &&t > 0\\
		u(0) = u_0
		\end{aligned}\right.
	```
	i.e. we have the case of $f(u, t) = b\, u - d \, u^2$
	when considering the above definition.

	This is the **logistic equation**, which in fact has multiple solutions.
	The solution relevant for population models has the form
	```math
	u(t) = \frac{b/d}{1 + \left( \frac{b}{d\, u_0} -1 \right) e^{-b t}}
	```
	As a plot shows, this solution varies smoothly from some initial population
	$u_0$ to a final population $b/d$:
"""

# ╔═╡ 31332dd6-1ef2-4181-a638-35980f470552
let
	u₀ = 1
	b = 4
	d = 0.5
	u(t) = b/d / (1 + (b/(d*u₀) - 1) * exp(-b*t))

	plot(u; xlims=(0, 3), label="b=4, d=0.5, u₀=1", xlabel=L"t", ylabel=L"u", title="Duck population", lw=2)
end

# ╔═╡ 1e65223e-9a12-4a4c-8b5e-31bfe4f813ec
md"""This problem is an example for the a class of problems one calls **initial value problem**, because based on some initial knowledge at $t=0$ one wants to know how a quantity (e.g. here the population) evolves.
"""

# ╔═╡ 64fe575e-0d47-4949-a9e8-2056ddee45df
md"""
!!! info "Definition: Initial-value problem"
	A scalar, first-order initial value problem (IVP) can be formulated as
	```math
	\tag{1}
	\left\{
	\begin{aligned}
	\frac{d u(t)}{d t} &= f(t, u(t)) && a ≤ t ≤ b \\
	u(a) &= u_0,
	\end{aligned}\right.
	```
	We call $t$ the **indepentent variable**, $u$ the **dependent variable**
	and $u_0$ the initial conditions.

	A **solution** of an initial-value problem is a differentiable, continuous
	function $u : \mathbb{R} \to \mathbb{R}$
	which makes both $u'(t) = f(t, u(t))$ (for all $a ≤ t ≤ b$)
	and $u(a) = u_0$ true equations.

	For the specific case where $f(t, u) = g(t) + u h(t)$ for two functions $g$ and $h$ the differential equation (1) is called **linear**.

Often (but not always) $t$ plays the role of time and (1) thus models the time-dependence of a quantity $u$.
"""

# ╔═╡ 1dbe5d72-17d9-4f20-b365-ad913cd607c3
md"""
## Numerical solution in DifferentialEquations.jl
For simple examples like the Duck problem above
analytic solutions can still be found with a little practice.
However, initial value problems quickly become more involved.
For example consider the problem
```math
\tag{2}
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= \sin\left[(u+t)^2\right] && 0 ≤ t ≤ 4 \\
u(0) &= -1,
\end{aligned}\right.
```
for which an analytical solution is not so easy to obtain.
"""

# ╔═╡ 6c2f1cd1-e79a-4790-9871-0a8825b5c02e
md"""
Before we introduce a few key ideas how to solve such problems,
we first need a reference technique to compare against.
Here the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) Julia package provides a range of production-grade numerical solvers for ODE problems.

We first translate our problem (2) into Julia code.
"""

# ╔═╡ 8f159939-4bf4-4151-bd3e-a84c5e1493db
begin
	f(u, t) = sin((t + u)^2)  # defines du/dt
	u₀ = -1.0                 # initial value
	tstart = 0.0              # Start time
	tend   = 4.0              # Final time
end;

# ╔═╡ 9b71f2b6-34bb-416a-947d-6e009f565bf4
md"""This we can now solve using `DifferentialEquations`, the details of which is beyond the scope of this course and hidden in the function `solve_reference`.
"""

# ╔═╡ 3f59bd74-1bb0-41df-bea2-bf5990a2b89a
function solve_reference(f, u₀, tstart, tend; kwargs...)
	# Construct ODEProblem object
	ivp = ODEProblem((u, p, t) -> f(u, t), u₀, (tstart, tend))

	# Solve problem using the Tsit5() algortihm, a pretty good default.
	solve(ivp, Tsit5(); kwargs...)
end

# ╔═╡ 9cb5fd8e-d923-434e-b543-a0c987a5b5a7
sol = solve_reference(f, u₀, tstart, tend);

# ╔═╡ 5d2886f2-87a5-43ae-9064-448478f957ca
md"The obtained solution can easily be visualised:"

# ╔═╡ 82321d1e-50b9-4362-afc5-d256a506bed7
plot(sol, label="solution", lw=2, xlabel="t", ylabel=L"u(t)",
	 title=L"\frac{du}{dt} = \sin((t+u)^2)")

# ╔═╡ 6b4ba69d-0a70-4c66-b4bc-63b717348471
md"""
## Existence and uniqueness of solutions

Let us remark that the **existence** of a solution to problem (1) is not guaranteed for all $t$ with $a ≤ t ≤ b$. For example consider the problem
```math
\frac{d u(t)}{dt} = \left(u(t)\right)^2,\, 0≤t≤2, \qquad u(0) = 1,
```
which has solution $u(t) = \frac{1}{1-t}$. This solution only exists for $t<1$.
When attempting a numerical solution beyond $t=1$ of such a problem we are faced with problems:
"""

# ╔═╡ 7a0fa3cd-a270-4947-b7ba-d30110139795
let
	f(u, t) = u^2
	u₀      = 1.0
	tstart  = 0
	tend    = 10
	
	sol = solve_reference(f, u₀, tstart, tend)
	plot(sol; lw=2, label="", xlabel=L"t", yaxis=:log, ylabel=L"u(t)", xlims=(0, 2))
end

# ╔═╡ c94d264b-ad1b-47fd-8b58-fd434cfe0f11
md"""
A second issue is that the solution may not necessarily be **unique**. Consider for example the equation
```math
\frac{d u(t)}{dt} = \sqrt{u(t)}, \, t>0, \qquad u(0) = 0,
```
which has the two solutions $u_1(t) = 0$ and $u_2(t) = \frac14 t^2$.
"""

# ╔═╡ 3dc8b41f-3ba0-4c0a-b017-731fd0958920
md"""
Both cases cause additional difficulties to numerical techniques, which we do not want to discuss here. For the scope of the course **we will assume that all IVP problems we consider, have a unique solution.**

If you are curious, the precise conditions for existence and uniqueness are given below:
"""

# ╔═╡ bc840e8d-77d9-4a4a-a908-3e9b4a8d253b
Foldable("Optional: Theorem governing Existence and uniquens of first-order ODEs",
md"""
!!! info "Theorem 1: Existence and uniquens of first-order ODEs"
	Given $f : \mathbb{R} \times \mathbb{R} \to \mathbb{R}$ a continuous function
	with respect to both arguments and **Lipschitz-continuous** with respect to its second argument, i.e. there exists a constant $L>0$ such that
	```math
	| f(t, x) - f(t, y)| \leq L |x-y|, \qquad \forall x,y \in \mathbb{R}, \, \forall t \text{ with } a ≤ t ≤ b
	```
	Then the ODE problem (1) has a unique solution $u(t)$ defined for all
	$a ≤ t ≤ b$ and $u$ is moreover continuously differentiable.
""")

# ╔═╡ 66affc7b-f0b0-48f4-b8bd-a1ff16d0357e
md"""
## Forward Euler
"""

# ╔═╡ f60ecd21-fda2-4925-ba86-ac5ad5b7c7d0
TODO(md"Mention the general idea of solving the problem: Knowing an approximation to the solution at some time $t_n$ we want to understand how we can find the solution at some later time $t_{n+1}$.")

# ╔═╡ 41cc0b76-5e30-41dd-a648-5491e1b4bd22
md"""

Let's assume we want to solve the IVP (1) in the interval $[a, b]$, i.e.
```math
\tag{3}
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= f(t, u(t)) && a ≤ t ≤ b \\
u(a) &= u_0,
\end{aligned}\right.
```
A natural idea is to divide this interval into $N$ subintervals $[t_n, t_{n+1}]$ of uniform length $h$, i.e.
```math
\tag{4}
t_n = a + n\, h \quad n = 0, \ldots, N  \qquad h = \frac{b-a}{N},
```
where $h$ is also called the **stepsize**.
At time $t_n$ we thus need to satisfy
```math
\frac{d\, u(t_n)}{dt} = f(t_n, u(t_n)),
```
where both $u(t_n)$ are unknown as well as the derivative on the LHS are unknown.
To reduce this to a single unknown we replace the derivative by one of the finite differences formulas discussed in the previous chapter on numerical differentiation. 

The simplest approach is the **forward Euler method** in which one approximates
$\frac{d\, u(t_n)}{d t}$ by forward finite differences, i.e.
```math
\tag{5}
\frac{d\, u(t_n)}{d t} ≈ D^+_h u(t_n) =\frac{1}{h} \Big( u(t_n+h) - u(t_n) \Big)
= \frac{1}{h} \Big( u(t_{n+1}) - u(t_n) \Big).
```
If by $u^{(n)}$ we denote an approximation to $u(t_n)$ and similarly for $u(t_{n+1})$ we obtain from (3):
```math
\tag{6}
\left\{
\begin{aligned}
\frac{u^{(n+1)} - u^{(n)}}{h} &= f(t_n, u^{(n)}), \qquad\forall n = 0, \ldots, N-1\\
u^{(0)} = u(a) &= u_0
\end{aligned}\right.
```
This defines a numerical scheme:

!!! info "Algorithm 1: Forward Euler method"
	Given an initial value problem $\frac{d u}{d t} = f(u, t)$, $u(0) = u_0$
	with $t \in [a, b]$ and nodes $t_n = a+n\, h$, $h = \frac{b-a}{N}$
	iteratively compute the sequence
	```math
	\tag{7}
	u^{(n+1)} = u^{(n)} + h\, f(t_n, u^{(n)})\qquad \forall n = 0, \ldots, N-1.
	```
	Then $u^{(n)}$ is *approximately* the solution $u$ at $t = t_n$.

Notice, that the replacement of the derivative $\frac{d\, u(t_n)}{d t}$ in (6) by the  forward finite difference approximation implies that the computed points $(u^{(0)}, u^{(1)}, \ldots, u^{(N)})$ are **not equal to** the function values of the exact solution $u$.
"""

# ╔═╡ c641cdf2-aabf-45e7-bb12-83d73a93b5b7
md"""
A basic implementation of Euler's method is shown below.
It expects the number of time intervals $N$
as well as the IVP --- in form of an `ODEProblem`.
All other problem parameters ($a$, $b$, $f$, ...) are then extract from
the latter object.
The output is the vector of nodes $t_n$ and the vector of approximate solution
values $u^{(n)} ≈ u(t_n)$.

"""

# ╔═╡ 9b0581c3-6fa8-44d8-8ebb-d2fda2841230
function forward_euler(f, u₀, a, b, N)
	# f:  Derivative function
	# u0: Initial value
	# a:  Start of the time interval
	# b:  End of the time interval

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration
	u[1] = u₀
    for i in 1:N
		u[i+1] = u[i] + h * f(u[i], t[i])
    end

    (; t, u)
end

# ╔═╡ f0164280-14fd-48a5-9ae8-8b31f4580da8
md"""
Let us compare this approach to the test problem (2)
```math
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= \sin\left[(u+t)^2\right] && 0 ≤ t ≤ 4 \\
u(a) &= 1,
\end{aligned}\right.
```
which we solved previously using `DifferentialEquations`. The variables `f`, `u₀`, `tstart` and `tend` already contains a setup of this initial value problem
"""

# ╔═╡ 6052f3cf-8a1e-4755-8724-a3f885702d62
f

# ╔═╡ 90a2176f-7986-4140-8e34-555993070324
u₀

# ╔═╡ 15e135e3-e95f-4494-8a20-2c571c4019f2
(tstart, tend)

# ╔═╡ 14edb601-ab10-4228-81a1-d27bb16ea202
md"""
and we can directly proceed to solving it:
"""

# ╔═╡ fb38ada2-26ce-42e9-bfdc-e4c677f8c501
md"and plot the result:"

# ╔═╡ 98f64da0-0d01-4c7e-8782-3d97dcf4c950
md"- Number of subintervals: `Neuler = ` $(@bind Neuler Slider(5:1:50; show_value=true, default=20)) "

# ╔═╡ 50007ba8-8160-48c1-b5c1-27a9332c1d58
res_euler = forward_euler(f, u₀, tstart, tend, Neuler);

# ╔═╡ 65807ddc-0c8b-4b88-bbc2-4ff3c7a1d2f8
let
	res_euler = forward_euler(f, u₀, tstart, tend, Neuler)

	plot(sol, label="reference", lw=2, xlabel="t", ylabel=L"u(t)",
	 title=L"\frac{du}{dt} = \sin((t+u)^2)", ylims=(-2, 0))

	plot!(res_euler.t, res_euler.u; label="Euler N=$Neuler", mark=:o, lw=2, ls=:dash)
end

# ╔═╡ ecb868c3-44a2-4f95-bafe-b3324770ee6b
md"""
We observe that the approximated solution becomes more and more accurate as we increase the number of sub-intervals $N$ and as the step size $h$ decreases.

But we also observe that for too small values of $N$ (e.g. $N=9$) that the Euler solution notably deviates from the reference.
"""

# ╔═╡ 667fbe5a-dbee-43a3-bb93-0160b61cf31e
md"""
Forward Euler is just one example of a broader class of numerical algorithms for initial value problems:

!!! info "Definition: One-step methods for initial value problems"
	A **one-step method** to solve (1) is a method of the form
	```math
	\tag{8}
	u^{(n+1)} = u^{(n)} + h\, \phi_f (u^{(n)}, t_n, h)
	\qquad \text{for $n=0, 1, \ldots, N-1$},
	```
	where $\phi_f(u, t, h)$ is a real-valued function
	and $u^{(0)} = u₀$.

!!! warning "Examples"
	- For **Forward Euler** we have $\phi_f(u^{(n)}, t_n, h) = f(t_n, u^{(n)})$
	- More involved methods exist. For example the [**midpoint method**](https://en.wikipedia.org/wiki/Midpoint_method), which we will discuss below, 	employs $\phi_f(u^{(n)}, t_n, h) = f\left(t_n + \frac{h}{2}, u^{(n)} + \frac{h}{2} f(t_n, u^{(n)})\right)$.
"""

# ╔═╡ 46f8284f-7751-4c1e-a617-386d86e7f4d0
md"""
## Error analysis

As discussed above applying Forward Euler to an initial value problem (1) with $N$ subintervals leads to a sequence of computed points $u^{(0)}, u^{(1)}, \ldots, u^{(N)}$, which approximate $u(t_0), u(t_1), \ldots, u(t_N)$ --- the solution $u$ evaluated at the nodes $\{t_i\}_{i=0}^N$. In this section our goal is to quantify the error
```math
\tag{9}
|e^{(n)}| = |u(t_n) - u^{(n)}| \qquad \text{for $n = 0, 1, \ldots, N$}
```
between the computed points and the exact values of the solution.

A notable feature of Forward Euler (Algorithm 1) is that the expression for computing $u^{(n+1)}$ results from two approximations:
1. Instead of evaluating the exact derivative $\frac{d\, u(t_n)}{d t}$ we employ a finite differences formula in (5).
1. Instead of using the exact value $u(t_n)$ to compute $u^{(n+1)}$ in (6) we employ the *approximation* $u^{(n)}$.
Both approximations contribute to the total error (9).
"""

# ╔═╡ 9385fd22-95f4-42b6-83ca-d4abbfe83a30
md"""
### Local truncation error
We first want to understand the error due to employing the finite difference formula
instead of the exact derivative.
To isolate this contribution let us assume that we do  have access to the *exact* values $u(t_n)$.
We than follow a Forward Euler step (7), just replacing the approximation $u^{(n)}$ by the exact $u(t_n)$.
This leads to the quantity 
```math
\tilde{u}^{(n+1)} = u(t_n) + h f(t_n, u(t_n)),
```
which thus isolates the error due to employing the finite difference formula.
We investigate the error of this quantity using a Taylor series
```math
\begin{aligned}
u(t_{n+1}) - \tilde{u}^{(n+1)}
&= u(t_n + h) - \Big[ u(t_n) + h\, \underbrace{f(t_n, u(t_n))}_{=u'(t_n)} \Big] \\
&= u(t_n) + h u'(t_n) + \frac12 h^2 u''(t_n) + O(h^3) - u(t_n) - h u'(t_n)\\
&= \frac12 h^2 u''(t_n) + O(h^3),
\end{aligned}
```
where we used the fact that $u$ satisfies the differential equation.
This expression only describes the error of a single step,
which advances us by a stepsize $h$ in going from $a$ to $b$.
To make the scaling of errors comparable as $h\to0$
it is usually more convenient ot investigate its size relative to the
chosen stepsize $h$, i.e.
 ```math
\tag{10}
\frac{u(t_{n+1}) - \tilde{u}^{(n+1)}}{h} = \frac{u(t_{n+1}) - u(t_n)}{h} - f(t_n, u(t_n)) = \frac12 h\, u''(t_n) + O(h^2).
```
This last expression is the so-called **local truncation error** of forward Euler.

Generalising towards the general class of one-step methods we define:
"""

# ╔═╡ 266c5a9a-8cd5-4dc0-a7de-492d64306aab
md"""
!!! info "Definition: Local truncation error"
	The **local truncation error** of a one-step method (8) is the quantity
	```math
	\tag{11}
	\tau^{(n)}_h = \frac{u(t_{n+1}) - u(t_n)}{h} - \phi_f\big(u(t_n), t_n, h\big)
	```

There are two ways to interpret the local truncation error $\tau^{(n)}_h$:
- As we discussed above (10) with respect to the Forward Euler method,
  it measures the **error of the numerical scheme** going from step $n$
  to $n+1$, provided that we start from the exact $u(t_n)$.
- From definition (9) we can also view  $τ^{(n)}_h$ as the (relative)
  **residual of the numerical method** when replacing the numerical solution 
  $u^{(n)}$ by the exact solution $u(t_n)$.
"""

# ╔═╡ d7150f04-e3e8-41f7-a173-9d458f72cdde
md"""
### Global error
As the name suggests the local truncation error $\tau^{(n)}$ only makes a statement about the error committed in a single isolated time step. However, what we were after was estimating the *global* error
```math
|e^{(n)}| = |u(t_n) - u^{(n)}| \qquad \text{for $n = 0, 1, \ldots, N$}.
```

Naively one might think that simply adding all local error contributions provides an estimate for this global error. However, this neglects something important. Namely, due to the algorithmic structure of one-step methods (8), where one obtains $u^{(n+1)}$ by adding to the approximation $u^{(n)}$, the error of step $n$ *propagates* to step $n+1$. In particular for small values of $N$ this can cause the error $|u(t_n) - u^{(n)}|$ in later time steps to grow signifficantly over the earlier ones. For example consider:
"""

# ╔═╡ 6d73f60b-0ccf-492f-9700-09d0bbf15124
let
	res_euler = forward_euler(f, u₀, tstart, tend, 7)
	plot(sol, label="reference", lw=2, xlabel="t", ylabel=L"u(t)",
	 title=L"\frac{du}{dt} = \sin((t+u)^2)", ylims=(-2, 1))
	plot!(res_euler.t, res_euler.u; label=L"Euler $N=7$", mark=:o, lw=2, ls=:dash)
end

# ╔═╡ 80748b7b-2546-4c4c-b60d-fdf7df98c8a7
md"""
where from time step 5 the solution is qualitatively wrong.

We therefore need to account for the propagation of error from one step to the next on top of the local error measured by $τ^{(n)}_h$. This is the subject of the following theorem:

!!! info "Theorem 2"
	Given a one-step method (8) with local truncation error satisfying
	```math
	\tag{12}
	|τ^{(n)}_h| ≤ C \, h^p
	```
	with a constant $C > 0$ and $p \in \mathbb{N}$ as well as[^1]
	```math
	\tag{13}
	\frac{\partial \phi_f(u, t, h)}{\partial u} ≤ L \qquad \text{for all $t \in [a, b]$ and all $h>0$}.
	```
	Then the global error satisfies
	for all $n = 0, 1, \ldots, N$
	```math
	\tag{14}
	|e^{(n)}| = |u(t_n) - u^{(n)}| ≤ \frac{C h^p}{L} \left(e^{L (t_n-a)} -1 \right)
	```
	as $h\to 0$,
	that is the numerical method (8) is **convergent with order $p$**.
"""

# ╔═╡ 5ac4eca2-0372-4b2d-add6-2347ce9461a3
details("Show the proof",
md"""
> **Proof:**
> Based on $e^{(n)} = u(t_n) - u^{(n)}$ we obtain
> ```math
> \begin{aligned}
> e^{(n+1)} - e^{(n)} &= [u(t_{n+1}) - u^{(n+1)}] - [u(t_n) - u^{(n)}] \\
> &= [u(t_{n+1}) - u(t_n)] - [u^{(n+1)} - u^{(n)}] \\
> &\stackrel{(8)}{=} [u(t_{n+1}) - u(t_n)] - h\, \phi_f (u^{(n)}, t_n, h).
> \end{aligned}
> ```
> By adding and subtracting $h\,\phi_f \big(u(t_n), t_n, h\big)$ we develop this to
> ```math
> \begin{aligned}
> e^{(n+1)} = e^{(n)} &+ \overbrace{\left[u(t_{n+1}) - u(t_n) - h\,\phi_f \big(u(t_n), t_n, h\big)\right]}^{=h\, τ^{(n)}_h}\\
> &+ h \left[ \phi_f \big(u(t_n), t_n, h\big) - \phi_f \big(u^{(n)}, t_n, h\big) \right]
> \end{aligned}
> ```
> Employing the triangle inequality yields
> ```math
> |e^{(n+1)}| ≤ |e^{(n)}| + h |τ^{(n)}_h| + h \left|\phi_f \big(u(t_n), t_n, h\big) - \phi_f \big(u^{(n)}, t_n, h\big)\right|
> ```
> Concerning the last term, the Fundamental Theorem of Calculus implies
> ```math
> \begin{aligned}
> \left|\phi_f \big(u(t_n), t_n, h\big) - \phi_f \big(u^{(n)}, t_n, h\big)\right|
> &= \left| \int_{u^{(n)}}^{u(t_n)} \frac{\partial \phi_f (u, t_n, h)}{\partial u} du \right|\\
> &≤ \int_{u^{(n)}}^{u(t_n)} \left|\frac{\partial \phi_f (u, t_n, h)}{\partial u}\right| du \\
> &\stackrel{(13)}{≤} L \left| u(t_n) - u^{(n)} \right|\\ &= L |e^{(n)}|.
> \end{aligned}
> ```
> Therefore using (12)
> ```math
> \begin{aligned}
> |e^{(n+1)}| ≤ (1 + h\,L)\, \left|e^{(n)}\right| + h |τ^{(n)}_h| \stackrel{(12)}{≤}  C\, h^{p+1} + (1 + h\,L)\,  \left|e^{(n)}\right|
> \end{aligned}
> ```
> with which we develop
> ```math
> \begin{aligned}
> |e^{(n+1)}| &≤ C\,h^{p+1} + (1 + h\,L)\,  \left|e^{(n)}\right| \\
> &≤ C\,h^{p+1} + (1 + h\,L) \left( C\,h^{p+1} + (1 + h\,L)\,\left|e^{(n-1)}\right|  \right) \\
> &= C\,h^{p+1} + (1 + h\,L) \, C\,h^{p+1} + (1+h\,L)^2\, \left|e^{(n-1)}\right| \\
> &\ \vdots\\
> &≤ \sum_{i=0}^n (1+h\,L)^i\, C\,h^{p+1} + (1+h\,L)^{n+1} |e_0| \\
> &= \frac{1-(1+h\,L)^{n+1}}{1 - (1+h\,L)} C\,h^{p+1} + (1+h\,L)^{n+1} |e_0|.
> \end{aligned}
> ```
> where in the last step we used the identity
> ```math
> \sum_{i=0}^n r^i = \frac{1-r^{n+1}}{1-r} \qquad \forall r \neq 1
> ```
> Now we notice that $e_0 = u(a) - u_0 = 0$, such that the expression simplifies to
> ```math
> |e^{(n+1)}| ≤ \frac{1-\left(1+h\,L\right)^{n+1}}{-h\,L} C\,h^{p+1} 
> = \frac{C\,h^{p}}{L} \left(\left(1+h\,L\right)^{n+1} -1 \right)
> ```
> Finally notice that $1 + x ≤ e^x$ for all $x ≥ 0$ and therefore that
> $(1+h\,L)^{n+1} \leq e^{(n+1)\,hL} = e^{(t_{n+1}-a) L}$. Replacing $n+1$ by $n$
> we obtain
> ```math
> |e^{(n)}| ≤ \frac{C\,h^{p}}{L} \left( e^{(t_{n}-a) L} - 1\right)
> ```
> which concludes the proof.

""")

# ╔═╡ 74a4bb06-d784-423c-9658-3cb9e4f545ae
md"""
[^1]: This result can also be prooven with minor modifications
	  under the conditions of Theorem 1, i.e. that
      $|\phi_f(u, t, h) - \phi_f(v, t, h)| \leq L |u-v|$ for all $u,v \in \mathbb{R}$,
      all $h>0$ and $a ≤ t ≤ b$.
"""

# ╔═╡ 285cbf5b-7097-47c0-950c-0c84d69a65bd
md"""
**Remark:**
One could paraphrase Theorem 2 by saying if the local truncation error $τ^{(n)}_h$ converges with order $p$, then the one-step methods also converges globally with order $p$. However, when thinking about the theorem in this way one should not forget the prefactor $(e^{L (t_n-a)} -1)$ in (14), which grows *exponentially in time*. Thus if we use a one-step method even for a higher order method our results may well deterioriate if $b \gg a$ ! This point we will pick up in the section on *Stability and implicit methods* below.

!!! warning "Example: Convergence of Forward Euler"
	As discussed in (10) for Forward Euler has local truncation error
	$\tau^{(n)}_h = \frac12 h\, u''(t_n) + O(h^2)$, such that 
	```math
	|\tau^{(n)}_h| ≤ \frac12 \max_{t\in[t_n,t_{n+1}]} |u''(t)| \, h
	```
	and using Theorem 2 we conclude that **Forward Euler converges linearly**.

	Let us verify this graphically on the test problem (2)
	```math
	\left\{
	\begin{aligned}
	\frac{d u(t)}{d t} &= \sin\left[(u+t)^2\right] && 0 ≤ t ≤ 4 \\
	u(0) &= -1,
	\end{aligned}\right.
	```
	which is contained in the variable `ivp` from before.
"""

# ╔═╡ 20096950-17ef-4d60-8059-dc7f00265349
let
	# Obtain extremely tight solution using DifferentialEquations.jl
	# Note, that ivp still contains the ODEProblem (2)
	u_exact = solve_reference(f, u₀, tstart, tend; reltol=1e-14, abstol=1e-14)

	Ns = [ round(Int, 5 * 10^k) for k in 0:0.5:3 ]
	errors = Float64[]
	for N in Ns
		res_euler = forward_euler(f, u₀, tstart, tend, N)
		u = res_euler.u
		t = res_euler.t

		# Compute global errors in each node tₙ
		global_error = [abs(u_exact(t[n]) - u[n]) for n in 1:length(u)]

		# Push the largest error to the array for plotting
		push!(errors, maximum(global_error))
	end
	
	plot(Ns, errors; mark=:o, label="Errors", lw=2, xaxis=:log, yaxis=:log,
		xlabel=L"N", ylabel="Largest global error",
		title="Convergence of Forward Euler")

	# Add line for perfect 1st order convergence:
	plot!(Ns, 0.5 ./ Ns, ls=:dash, lw=2, label=L"O(n^{-1})")
end

# ╔═╡ 38a5d61d-aff7-44dd-90e2-e690cd2643f7
md"""
## Runge-Kutta methods

Similar to the previous chapters we now ask the question: How can we improve the convergence of a numerical method for an IVP ? The answer to this question leads to the major and most widely employed types of methods for solving initial-value problems, the **Runge-Kutta** (RK) **methods**.

We start with a second-order RK method, namely the midpoint method.
"""

# ╔═╡ 04b3c441-a527-4b15-99e1-042479c4fdf2
md"""
### Optional: Midpoint method
The goal of the midpoint method is to construct a 
second-order method for solving the IVP (1)
```math
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= f\Big(t, u(t)\Big) && a ≤ t ≤ b \\
u(a) &= u_0.
\end{aligned}\right.
```
We stay within the framework of one-step methods, i.e. we iterate
for $n=0, 1, \ldots, N-1$
```math
u^{(n+1)} = u^{(n)} + h\, \phi_f (u^{(n)}, t_n, h).
```
The ideal method would obtain the exact $u(t_{n+1})$ from the exact $u(t_{n})$,
implying a local truncation error $τ^{(n)}_h = 0$. To build such a method
we expand
```math
\tag{15}
\begin{aligned}
u(t_{n+1}) = u(t_{n} + h) 
&= u(t_n) + h\,u'(t_n)  \hspace{27pt}+ \frac12 h^2 u''(t_n) + O(h^3)\\
&= u(t_n) + h\,f\Big(t_n, u(t_n)\Big) + \frac12 h^2 u''(t_n) + O(h^3) .
\end{aligned}
```
where we used the definition of the IVP to replace $u'$.
If in this equation we keep only the first two terms we recover
the Forward Euler method (7). To obtain more accuracy
we therefore need to compute or approximate $u''(t_n)$ as well.
Using the chain rule we note
```math
\frac{d^2 u}{d t^2} = \frac{d f}{d t}
\stackrel{(*)}{=} \frac{\partial f}{\partial t}
+ \frac{\partial f}{\partial u}
\frac{\partial u}{\partial t}
= \frac{\partial f}{\partial t}
+ \frac{\partial f}{\partial u} f,
```
where notably step $(*)$ is neccessary because both $f$ and its argument $u$ depend on $t$. Inserting this expression into (15) we obtain
```math
\tag{16}
\begin{aligned}
u(t_{n+1}) = u(t_{n} + h) 
&= u(t_n) + h\,\left[f\Big(t_n, u(t_n)\Big) + \frac h2 
\frac{\partial f}{\partial t}\Big(t_n, u(t_n)\Big)
\right.\\&\hspace{70pt}\left.
+ \textcolor{blue}{\frac h2} \frac{\partial f}{\partial u}\Big(t_n, u(t_n)\Big)
\textcolor{blue}{f\Big(t_n, u(t_n)\Big)}
\right] + O(h^3).
\end{aligned}
```
Since computing and implementing the partial derivatives $\frac{\partial f}{\partial t}$ and $\frac{\partial f}{\partial u}$ can in general become ticky,
we will also approximate these further.
Expanding $f$ in a multi-dimensional Taylor series we observe
```math
\begin{aligned}
f\Big(t_n + ζ, u(t_n) + \textcolor{blue}{ξ} \Big)
= f\Big(t_n, u(t_n) \Big) &+ ζ\, \frac{\partial f}{\partial t}\Big(t_n, u(t_n)\Big)
+ \textcolor{blue}{ξ}\, \frac{\partial f}{\partial u}\Big(t_n, u(t_n)\Big)\\
&+ O(ζ^2 + |ζ\textcolor{blue}{ξ}| + \textcolor{blue}{ξ}^2).
\end{aligned}
```
With the choice $ζ = \frac h2$ and $\textcolor{blue}{ξ = \frac h2 f\big(t_n, u(t_n)\big)}$
this expressions becomes equal to lowest order with the term
in the square brackets of (16). This leads to
```math
\begin{aligned}
u(t_{n+1}) = u(t_{n} + h) 
&= u(t_n) + h\,f\Big(t_n + \frac h2, u(t_n) + \textcolor{blue}{\frac h2 f\big(t_n, u(t_n)\big)} \Big)\\
&\hspace{35pt}+ O(h^3 + h ζ^2 + h |ζξ| + h ξ^2).
\end{aligned}
```
"""

# ╔═╡ 4d20e43b-ba3d-4eab-adfc-f3a7c6feeb0e
md"""
This is the **midpoint method** we mentioned earlier:

!!! info "Algorithm 2: Midpoint method"
	Given an initial value problem $\frac{d u}{d t} = f(u, t)$, $u(0) = u_0$
	with $t \in [a, b]$ and nodes $t_n = a+n\, h$, $h = \frac{b-a}{N}$
	iteratively compute the sequence
	```math
	u^{(n+1)} = u^{(n)} + h\, \phi_f(u^{(n)}, t_n, h)\qquad \forall n = 0, \ldots, N-1.
	```
	with $\phi_f(u^{(n)}, t_n, h) = f\left(t_n + \frac{h}{2}, \textcolor{green}{u^{(n)} + \frac{h}{2} f(t_n, u^{(n)}})\right)$

From our derivation the local trucnation error for the midpoint method can be identified as
```math
\tau^{(n)}_h = \frac{u(t_{n+1}) - u(t_n)}{h} - \phi_f\big(u(t_n), t_n, h\big)
= \frac{1}{h} O(h^3 + h ζ^2 + h |ζξ| + h ξ^2) = O(h^2).
```
From Theorem 2 we observe that the **midpoint method is** indeed **a second-order method**.
"""

# ╔═╡ 3430cee3-d5ac-4f7d-99ce-713008cd3553
md"""
Runge-Kutta methods, such as the midpoint methods, are what one calls **multi-stage methods**. To make this term clear, let us rewrite equation (7) in two stages
by isolating the green part:
```math
\begin{aligned}
	\textcolor{green}{u^{(n+\frac12)}} &\,\textcolor{green}{:= u^{(n)} + \frac h2 f(t_n, u^{(n)})} \\
	u^{(n+1)} &= u^{(n)} + h\, f\big(t_n + \frac h2, \textcolor{green}{u^{(n+\frac12)}}\big)
\end{aligned}
```
Written like this we notice that the **first stage** (in green)
performs a Forward Euler step with **half the stepsize**,
namely from time $t_n$ to $t_n + \frac h2$.
The **second stage** then performs an Euler-style update over the whole
timestep, but employing the slope $u^{(n+\frac12)}$
from the first stage in the computation of $f$.

This interpretation also explains why the method is called **midpoint method**.

Employing again the `ODESystem` to supply the setup of our problem,
we obtain an implementation of the midpoint method as:
"""

# ╔═╡ 76978387-3b7b-4644-a27e-246f84186eb9
function midpoint(f, u₀, a, b, N)
	# f:  Derivative function
	# u0: Initial value
	# a:  Start of the time interval
	# b:  End of the time interval

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration ... this is what changes over forward_euler
	u[1] = u₀
    for i in 1:N
		uhalf  = u[i] + h/2 * f(u[i],  t[i]      )
		u[i+1] = u[i] + h   * f(uhalf, t[i] + h/2)
    end

    (; t, u)
end

# ╔═╡ 87345743-d3a5-4ec4-bf66-15e6209a729a
md"""In comparison with Forward Euler we notice this method to be clearly more accurate for the case of rather small number of timesteps:"""

# ╔═╡ 2c3ad400-acdd-41bb-ac40-16cdc4b6fb07
let
	N = 9
	res_euler    = forward_euler(f, u₀, tstart, tend, N)
	res_midpoint = midpoint(f, u₀, tstart, tend, N)
	plot(sol, label="reference", lw=2, xlabel="t", ylabel=L"u(t)",
	 title=L"\frac{du}{dt} = \sin((t+u)^2)", ylims=(-2, 0.5))
	plot!(res_euler.t, res_euler.u; label=L"Euler $N=9$", mark=:o, lw=2, ls=:dash)
	plot!(res_midpoint.t, res_midpoint.u; label=L"Midpoint $N=9$", mark=:o, lw=2, ls=:dash)
end

# ╔═╡ 262a2a67-3fd3-497c-b56c-71fe26abffc1
md"""
### Optional: Higher-order Runge-Kutta methods

The ideas we used for constructing the midpoint method point to a clear path how one could construct even higher-order methods: All we need to do is to introduce further intermediate stages, i.e. half, third or other intermediate timesteps. In this way we can generate additional equations and effectively match the higher-order derivatives in the Taylor expansion (15), which will in turn reduce the local truncation error $τ^{(n)}_h$. The methods generated in this form are called **Runge-Kutta methods**.

The algebra required to work out the details grows in complexity, so we will not attempt to do this here and only present a general overview. Constructing an $s$-stage Runge-Kutta methods leads to the set of equations
```math
\begin{aligned}
v_1 &= h\,f(t_n, \hphantom{+ c_1\,h,\ }u_n) \\
v_2 &= h\,f(t_n + c_1\,h,\ u_n + a_{11}\, v_1) \\
v_3 &= h\,f(t_n + c_2\,h,\ u_n + a_{21}\, v_1 + a_{22}\, v_2) \\
&\vdots\\
v_s &= h\,f\left(t_n + c_{s-1}\,h,\ u_n + \sum_{i=1}^{s-1} a_{s-1,i}\,v_i\right)\\
u^{(n+1)} &= u^{(n)} + \sum_{i=1}^{s} b_i\, v_i.
\end{aligned}
```
where specifying both the number of stages $s$ as well as the constants $a_{ij}$, $b_i$ and $c_i$ uniquely determines an RK method.
Both one-step methods we discussed so far actually match this framework:
- **Forward Euler:** $s=1$, $b_1 = 1$ and no $a_{ij}$ and no $c_i$.
- **Midpoint method:** $s=2$, $b_1 = 0$, $b_2=1$, $c_1 = \frac12$ and $a_{11} = \frac12$.
"""

# ╔═╡ 850928f3-ab31-42c6-a51a-0dca28175bdf
md"""
Additionally we want to specify one additionall **fourth-order RK** method,
often called **RK4**, which is the most commonly used IVP approach:

!!! info "Algorithm 3: Fourth-order Runge-Kutta method (RK4)"
	Given an initial value problem $\frac{d u}{d t} = f(u, t)$, $u(0) = u_0$
	with $t \in [a, b]$ and nodes $t_n = a+n\, h$, $h = \frac{b-a}{N}$
	iteratively compute the sequence
	```math
	\begin{aligned}
	v_1 &= h\,f(t_n\hphantom{\,\, + \frac h2}, u^{(n)}) \\
	v_2 &= h\,f(t_n + \frac h2, u^{(n)} + \frac{v_1}{2}) \\
	v_3 &= h\,f(t_n + \frac h2, u^{(n)} + \frac{v_2}{2}) \\
	v_4 &= h\,f(t_n + h, \ \, u^{(n)} + v_3) \\[0.5em]
	u^{(n+1)} &= u^{(n)} + \frac{v_1}{6} + \frac{v_2}{3} + \frac{v_3}{3} + \frac{v_4}{6}
	\end{aligned}
	```

Let us mention in passing, that our reference method `Tsit5()` is also based on a  Runge-Kutta scheme.

An implementation of RK4 is given by:
"""

# ╔═╡ 5f7562f6-4329-4c1e-b5d4-74c9b1980346
function rk4(f, u₀, a, b, N)
	# f:  Derivative function
	# u0: Initial value
	# a:  Start of the time interval
	# b:  End of the time interval

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration ... this is what changes over forward_euler
	u[1] = u₀
    for i in 1:N
		v₁ = h * f(u[i],        t[i]      )
		v₂ = h * f(u[i] + v₁/2, t[i] + h/2)
		v₃ = h * f(u[i] + v₂/2, t[i] + h/2)
		v₄ = h * f(u[i] + v₃,   t[i] + h  )

		u[i+1] = u[i] + (v₁/6 + v₂/3 + v₃/3 + v₄/6)
    end

    (; t, u)
end

# ╔═╡ f0f756fc-f7fe-459c-8419-4e60c741fce5
md"""
Let us compare all methods we saw in this chapter:
- `N = ` $(@bind N Slider(5:1:25; show_value=true, default=7))
"""

# ╔═╡ f9e8c1ec-3bc6-40d3-b7e3-f465bd53ecab
let
	res_euler    = forward_euler(f, u₀, tstart, tend, N)
	res_midpoint = midpoint(f, u₀, tstart, tend, N)
	res_rk4      = rk4(f, u₀, tstart, tend, N)
	plot(sol, label="reference", lw=2, ylabel=L"u(t)", xlabel="",
	 title=L"\frac{du}{dt} = \sin((t+u)^2)", ylims=(-2, 0.5), titlefontsize=12)
	plot!(res_euler.t, res_euler.u; label="Euler", mark=:o, lw=2, ls=:dash, markersize=3)
	plot!(res_midpoint.t, res_midpoint.u; label="Midpoint", mark=:o, lw=2, ls=:dash, markersize=3)
	p = plot!(res_rk4.t, res_rk4.u; label="RK4", mark=:o, lw=2, ls=:dash, markersize=3)

	q = plot(; yaxis=:log, ylims=(1e-6, 1), legend=false, ylabel=L"$|u - u_{\textrm{ref}}|$", xlabel=L"t", title="Error", titlefontsize=10)
	plot!(q, res_euler.t, abs.(res_euler.u - sol.(res_euler.t) .+ 1e-16);
	      label="Euler", mark=:o, lw=2, markersize=3, ls=:dash, c=2)
	plot!(q, res_midpoint.t, abs.(res_midpoint.u - sol.(res_midpoint.t) .+ 1e-16);
	      label="Midpoint", mark=:o, lw=2, markersize=3, ls=:dash, c=3)
	plot!(q, res_rk4.t, abs.(res_rk4.u - sol.(res_rk4.t) .+ 1e-16);
	      label="RK4", mark=:o, lw=2, markersize=3, ls=:dash, c=4)

	plot(p, q; layout=grid(2, 1; heights=[0.75, 0.25]))
end

# ╔═╡ bdb783bc-1648-43f0-b706-462d2a5d6dcf
md"""
RK4 converges faster than the other two methods. As the name suggests it is indeed a fourth-order method:
"""

# ╔═╡ 99922f12-34b6-4a9d-b3af-2b61708725fa
begin
	# Obtain extremely tight solution using DifferentialEquations.jl
	# Note, that ivp still contains the ODEProblem (2)
	u_exact = solve_reference(f, u₀, tstart, tend; reltol=1e-14, abstol=1e-14)

	# Function to compute the global error
	function global_error(t, u, u_exact)
		errors = [abs(u_exact(t[n]) - u[n]) for n in 1:length(u)]
		maximum(errors)
	end

	Ns = [ round(Int, 5 * 10^k) for k in 0:0.5:3 ]
	errors_euler    = Float64[]
	errors_midpoint = Float64[]
	errors_rk4 = Float64[]
	for N in Ns
		result_euler = forward_euler(f, u₀, tstart, tend, N)
		push!(errors_euler, global_error(result_euler.t, result_euler.u, u_exact))
		
		result_midpoint = midpoint(f, u₀, tstart, tend, N)
		push!(errors_midpoint, global_error(result_midpoint.t, result_midpoint.u, u_exact))
		
		res_rk4 = rk4(f, u₀, tstart, tend, N)
		push!(errors_rk4, global_error(res_rk4.t, res_rk4.u, u_exact))
	end

	p = plot(; title="Convergence of RK methods", xaxis=:log, yaxis=:log,
	           xlabel=L"n", ylabel="Largest global error", legend=:bottomleft,
	           xlims=(1, 3e4), ylims=(1e-13, 10))
	plot!(p, Ns, errors_euler;    lw=2, mark=:o, c=1, label="Euler")
	plot!(p, Ns, errors_midpoint; lw=2, mark=:o, c=2, label="Midpoint")
	plot!(p, Ns, errors_rk4;      lw=2, mark=:o, c=3, label="RK4")

	# Add line for perfect 1st order convergence:
	plot!(Ns, 0.08(first(Ns) ./ Ns).^1, ls=:dash, lw=2, label=L"O(n^{-1})", c=1)
	plot!(Ns, 0.08(first(Ns) ./ Ns).^2, ls=:dash, lw=2, label=L"O(n^{-2})", c=2)
	plot!(Ns, 0.08(first(Ns) ./ Ns).^4, ls=:dash, lw=2, label=L"O(n^{-4})", c=3)
end

# ╔═╡ 0bd09a85-2fe9-491d-b6f9-f2d74872e195
md"""
On the other hand the four stages of RK4 also have a **disadvantage**: For each timestep four times as many evaluations of the function $f$ are required than Forward Euler. Since this is typically the cost-dominating step we might ask if RK4 should always be employed or if other schemes with less stages (like Midpoint or Euler) are sometimes more appropriate.

To answer this question we take our previous plot and multiply the x-values by the number of stages in order to track the error versus the number of evaluations of the function $f$. With this we obtain:
"""

# ╔═╡ f5c8f93c-c09a-4c4f-b6c4-18b12ff56878
let
	p = plot(; title="Convergence of RK methods", xaxis=:log, yaxis=:log,
	           xlabel="number of function evaluations",
	           ylabel="Largest global error", legend=:bottomleft,
			   xlims=(1, 3e4), ylims=(1e-13, 10))

	# Euler has only a single stage
	plot!(p, 1Ns, errors_euler;    lw=2, mark=:o, c=1, label="Euler")

	# Midpoint has two stages
	plot!(p, 2Ns, errors_midpoint; lw=2, mark=:o, c=2, label="Midpoint")

	# RK4 has four stages
	plot!(p, 4Ns, errors_rk4;      lw=2, mark=:o, c=3, label="RK4")
end

# ╔═╡ 591d6a3a-fc71-4bde-b499-db10637e0945
md"""
We observe that even though RK4 needs four $f$ evaluations per time step it still provides the best accuracy versus $f$-evaluations ratio for almost the entire range of evaluations we considered.
"""

# ╔═╡ b6022950-3711-470c-b7fe-ae02232a0749
md"This function rapidly decays to zero as $t \to \infty$:"

# ╔═╡ a09b8d7b-9687-487d-940d-93c6957b027c
md"""Applying our numerical methods from the previous sections, we would expect to recover this behaviour. So let's check this.

We first set up the `ODEProblem`:"""

# ╔═╡ d57e6a62-cc98-4dbc-b471-b21448c0343e
md"""... and try it on our previously implemented methods:"""

# ╔═╡ 1f1af0d4-3c5d-4e38-9eda-c307100bd61f
md"""
- `C = `: $(@bind C Slider(10 .^ (-4:0.2:1); show_value=true, default=10^0.8))
- `Ndecay = `: $(@bind Ndecay Slider(15:5:50; show_value=true, default=15))
"""

# ╔═╡ 85b121ef-9f4c-483e-a04b-a49c1c6c4624
md"""
## Stability and implicit methods

Let us consider the innocent looking initial value problem
```math
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= -C && t > 0 \\
u(0) &= 10.
\end{aligned}\right.
```
This model governs for example the decay of a radioactive species with initial concentration $10$ over time. The rate of decay is $C$ ≈  $(round(C; sigdigits=3)).

By simple integration we find the exact solution of this problem as:
"""

# ╔═╡ 5348acce-d58b-4b33-97ad-44f44717e1fe
u(t) = 10 * exp(-C * t)

# ╔═╡ f2bde123-95b7-4370-9914-577985c4efeb
plot(u; xlims=(0, 10), lw=2, title="du/dx = $(round(-C; sigdigits=3))", titlefontsize=12, label="Reference", xlabel=L"t", ylabel=L"u(t)")

# ╔═╡ af8f49e2-a876-4bb9-8036-46c1f510f369
begin
	dcy_f(u, t) = -C * u
	dcy_u₀      = 10.0
	dcy_tstart  =  0.0
	dcy_tend    = 10.0
end;

# ╔═╡ e4fa9b2d-8a81-4439-84ff-8529a7421ae9
let
	h = (dcy_tend - dcy_tstart) / Ndecay  # Stepsize of our numerical methods
	
	soln_euler    = forward_euler(dcy_f, dcy_u₀, dcy_tstart, dcy_tend, Ndecay)
	soln_midpoint = midpoint(dcy_f, dcy_u₀, dcy_tstart, dcy_tend, Ndecay)
	soln_rk4      = rk4(dcy_f, dcy_u₀, dcy_tstart, dcy_tend, Ndecay)

	p = plot(u; lw=2, title="du/dx = $(round(-C; sigdigits=3))", titlefontsize=12, label="Reference", xlabel=L"t", ylabel=L"u(t)", xlims=(-0.1, 10.1))
	
	if h ≥ 2/C  # See discussion below to understand
	            # why this formula appears here
		ylims!(p, (-1_000, 10_000))
	end

	plot!(p, soln_euler.t, soln_euler.u;
		  mark=:o, c=2, lw=2, markersize=3, ls=:dash, label="Forward Euler")
	plot!(p, soln_midpoint.t, soln_midpoint.u;
		  mark=:o, c=3, lw=2, markersize=3, ls=:dash, label="Midpoint")
	plot!(p, soln_rk4.t, soln_rk4.u;
		  mark=:o, c=4, lw=2, markersize=3, ls=:dash, label="RK4")
end

# ╔═╡ 60a34642-3824-4a3f-a11c-e76039c91e4a
md"""We notice for too few nodal points or too large decay constants $C$ our numerical methods **all fail** to recover the decay behaviour.
In other words while the exact solution satisfies $\lim_{t\to\infty} u(t) = 0$,
the numerical methods are **all qualitatively wrong**:
instead of the reproducing the correct long-time limit numerically,
i.e. $\lim_{n\to\infty} u^{(n)} = 0$,
the numerical solution actually grows over time.

Recall that in Theorem 2 we found that the global error satisfies
```math
|e^{(n)}| ≤ \frac{C h^p}{L} \left(e^{L (t_n-a)} -1 \right)
= \frac{C h^p}{L} \left(e^{L\,hn} -1 \right).
```
In the limit of infinite time,
i.e. $n\to \infty$, the upper bound on the RHS thus grows exponentially.
This bound is in general pessimistic, i.e. not all numerical methods may actually reproduce the exponentially increasing error behaviour, but it clearly shows that unless one demands additional properties from a method for IVPs, one may fail to reproduce a limiting behaviour like $\lim_{t\to\infty} u(t) = 0$.

Without going into further details the property we ask for is called **stability**:

!!! info "Definition: Absolute stability"
	Given an initial value problem (1) with exact solution $u(t)$ satisfying
	$\lim_{t\to\infty} u(t) = u_\ast$ for any initial condition $u_0$. In other words we have a problem where the solution converges to $u_\ast$ independent of the initial value $u_0$.

	Then a numerical method is called **absolutely stable** for a fixed stepsize $h$ if
	```math
	\lim_{n\to\infty} u^{(n)} = u_\ast \qquad \text{for any } u_0.
	```

	A method is further called **unconditionally absolutely stable** if it is stable for all $h > 0$.

In the case of our decay problem, one can show for example that the forward Euler method is absolutely stable if
```math
h < \frac{2}{C}
```
--- which explains the condition we chose to switch the axis limits in the plot above.

In particular we also observe that none of our methods is unconditionally absolutely stable.
"""

# ╔═╡ 2263c105-7d59-43a9-9196-b4762a53b1b0
md"""
### Backward Euler

The construction and implementation of unconditionally absolutely stable methods is in general more involved than the methods we considered so far. We will therefore only construct a single approach here.

For this we return to the derivation of the Forward Euler method. Recall that in order to solve the IVP (1)
```math
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= f(t, u(t)) && a ≤ t ≤ b \\
u(a) &= u_0,
\end{aligned}\right.
```
we introduced a time discretisation into $N$ subintervals $[t_n, t_{n+1}]$,
such that in each interval we needed to solve the problem
```math
\tag{20}
\frac{d\, u(t_n)}{dt} = f(t_n, u(t_n)).
```
Instead of employing the forward finite differences formula,
leading to the Forward Euler method (7),
we now instead employ the backward finite differences, leading to
```math
\tag{21}
\frac{d\, u(t_{n+1})}{d t} ≈ D^-_h u(t_{n+1}) = \frac{1}{h} \Big( u(t_{n+1}) - u(t_{n+1}-h) \Big)
= \frac{1}{h} \Big( u(t_{n+1}) - u(t_n) \Big).
```
Following the same idea as before,
i.e. to employ for $u(t_{n+1})$ the approximation $u^{(n+1)}$
of the $n+1$-st time step and instead of $u(t_n)$ employ $u^{(n)}$
we obtain from (20) and (21):
```math
\tag{22}
\frac{u^{(n+1)} - u^{(n)}}{h} = f(t_{n+1}, u^{(n+1)}), \qquad\forall n = 0, \ldots, N-1,
```
which is the **Backwards Euler method**.
Notice, that in contrast to the case of Forward Euler we cannot immediately deduce an explicit update formula like (7) from this expression, since the dependency on $u^{(n+1)}$ is both on the LHS as well as in $f$. This is why such methods are called **implicit methods** as $u^{(n+1)}$ is only implicitly defined.

To obtain $u^{(n+1)}$ from $u^{(n)}$ one needs to solve (22) iteratively:
For this we introduce the map
```math
\tag{23}
x = g^{(n)}(x) = u^{(n)} + h \, f(t_{n+1}, x)
```
and notice that its fixed-point $x_\ast = g^{(n)}(x_\ast)$ is exactly $u^{(n+1)}$.
For each time step $n$ we thus need to solve a fixed-point problem using one of the methods we described in chapter 4.

In summary our algorithm becomes:

!!! info "Algorithm 4: Backward Euler method"
	Given an initial value problem $\frac{d u}{d t} = f(u, t)$, $u(0) = u_0$
	with $t \in [a, b]$ and nodes $t_n = a+n\, h$, $h = \frac{b-a}{N}$
	we iterate for $n = 1, \ldots N$:
	- Find a fixed-point $x_\ast$ of the map $g^{(n)}(x) = u^{(n)} + h \, f(t_{n+1}, x)$, e.g. using the Newton method.
	- Set $u^{(n+1)} = x_\ast$
"""

# ╔═╡ 5075f5e4-af28-4a93-b8e6-d71f96cf2b89
md"""An implementation of Backward Euler employing Newton's method to solve the fixed-point problem is given below. The implementation of Newton's method is repeated in the  appendix."""

# ╔═╡ 37bb847f-aa6b-4efe-8538-e580159d2892
md"""
- `Nbw = ` $(@bind Nbw Slider(15:5:40; show_value=true, default=15) )
"""

# ╔═╡ e031fb26-dc5e-44e9-ba4c-a642c62b7930
md"""
While Forward Euler stays absolutely stable only for large values of `Nbw` (respectively small values of $h$), Backward Euler stays stable no matter what value of `Nbw` is chosen.

Let us conclude by mentioning that similar to Forward Euler, **Backward Euler** is also only a **first order method**. Higher-order implicit methods can also be constructed. E.g. the Runge-Kutta family of methods can be extended to the implicit setting as well. An example is the [Crank–Nicolson method](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method), a second-order implicit Runge-Kutta method.

In standard libraries, such as DifferentialEquatios.jl a [zoo of ODE methods](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) is typically available.
"""

# ╔═╡ 35b8a591-62ea-4836-aa02-e1aaa584621d
md"""
## Optional: Higher-order differential equations

Consider the setting of a mass $m$ hanging on a spring from the ceiling:
"""

# ╔═╡ 6b5ac847-ea65-4d44-9379-5cc755908521
TODO("Image")

# ╔═╡ d8e1eb69-9f8f-47c3-a3f8-bf7f725959de
md"""
In its rest position show above, gravity is exactly cancelled by the upward force excerted by the spring, such that we can ignore gravity for our discussion. Now displacing the mass by a length $x$ leads to a restoring force $- k x$ where $k$ is the spring constant of the spring. By Newton's law this accelerates the mass by $- k/m x$. Employing $m=1$ for simplicity this leads to the **second-order ODE**
```math
\left\{
\begin{aligned}
\frac{d^2 x(t)}{d t^2} &= -k\, x(t) && t > 0 \\
x(0) &= 0,
\end{aligned}\right.
```
which is usually referred to as the **simple harmonic oscillator model**,
abbreviated HO.
To solve this problem numerically we first introduce the velocity function $v(t)$ and rewrite it as a system of first-order ODEs:
```math
\tag{17}
\left\{
\begin{aligned}
\frac{d\, x(t)}{d t} &= v(t) && t > 0 \\
\frac{d\, v(t)}{d t} &= -k\, x(t) && t > 0 \\
x(0) &= 0\\
v(0) &= 0.5,
\end{aligned}\right.
```
where we choose to start at the rest position, but with an initial velocity of $0.5$.
Collecting the position and velocity into a single vector-valued variable, i.e.
```math
\textbf{u}(t) = \begin{pmatrix} u_1(t) \\ u_2(t) \end{pmatrix} = \begin{pmatrix} x(t) \\ v(t) \end{pmatrix},
\qquad
\frac{d \textbf{u}(t)}{d t} = \begin{pmatrix} \frac{d\, u_1(t)}{dt} \\ \frac{d\, u_2(t)}{dt} \end{pmatrix} = \begin{pmatrix} \frac{d\, x(t)}{dt} \\ \frac{d\, v(t)}{dt} \end{pmatrix},
\qquad
\textbf{u}_0 = \begin{pmatrix} 0 \\ 0.5 \end{pmatrix}
```
and introducing the function $\textbf{f} : \mathbb{R}^2 \times \mathbb{R} \to \mathbb{R}^2$ with
```math
\textbf{f}(\textbf{u}, t) = \begin{pmatrix} u_2(t) \\ -k\, u_1(t) \end{pmatrix}
```
we can rewrite this as
```math
\tag{18}
\left\{
\begin{aligned}
\frac{d\, \textbf{u}(t)}{d t} &= \mathbf{f}(\mathbf{u}, t) && t > 0\\
\textbf{u}(0) &= \textbf{u}_0.
\end{aligned}\right.
```
Notice, that this has exactly the same structure as (1),
just all quantities are vector valued.

All algorithms we discussed in this chapter so far can be extended to this setting. For example Algorithm 1 (Forward Euler) simply becomes:
"""

# ╔═╡ d09ee1f0-c9ca-4888-9573-6ca0c105dccd
md"""
!!! info "Algorithm 1a: Forward Euler for vector-valued problems"
	Given an initial value problem $\frac{d \mathbf{u}}{d t} = \textbf{f}(t, \mathbf{u})$, $\mathbf{u}(0) = \mathbf{u}_0$
	with $t \in [a, b]$ and nodes $t_n = a+n\, h$, $h = \frac{b-a}{N}$
	iteratively compute the sequence
	```math
	\tag{19}
	\textbf{u}^{(n+1)} = \textbf{u}^{(n)} + h\, \textbf{f}(t_n, \textbf{u}^{(n)})\qquad \forall n = 0, \ldots, N-1.
	```

which is identical to Algorithm 1, just with all quantities replaced by their vector-valued analogues. In fact even the implementations of `forward_euler` can just be used without changing a single line of code as we will demonstrate now.

First we setup the HO model function and parameters:
"""

# ╔═╡ f35f9be4-74ae-4a05-ab6a-55aa799d98a2
k = 2  # Force constant

# ╔═╡ 01e9b24f-7819-4e72-9463-2a6aca835d22
begin
	# To setup f we assume u is a vector of size 2 and return a vector of size 2
	ho_f(u, t) = [     u[2];
				  -k * u[1]]

	# As the initial value again we supply a vector of size 2
	ho_u₀ = [0.0;
		     1.0]

	# and we are interested in the behaviour from 0 to tend (defined below)
	ho_tstart =  0.0
	ho_tend   = 10.0
end

# ╔═╡ 3834116a-be19-4690-a4be-542d7824f410
md"""Then running the dynamics just takes a call to `forward_euler` as before:"""

# ╔═╡ 612e93db-322b-4b7c-a13e-cfd18058178a
ho_euler = forward_euler(ho_f, ho_u₀, ho_tstart, ho_tend, 50);

# ╔═╡ 94a39f52-0653-4f01-a08a-ee9bc887da7d
let
	x = [u[1] for u in ho_euler.u]  # Extract the particle position
	plot(ho_euler.t, x, mark=:o, c=2, lw=2, xlabel=L"t", ylabel=L"x(t)", label="", title=L"\frac{d^2 x}{d t^2} = -k x", titlefontsize=12, markersize=3, ls=:dash)
end

# ╔═╡ 823685fd-221a-40c3-855a-20825fda7bf3
md"""
In fact the `midpoint` and `rk4` implementations we provide are similarly generic with respect to being used with scalar-valued or vector-valued $u$. Applying them all we obtain:
"""

# ╔═╡ e195e04c-fd23-4f92-8241-235d69404d8e
md"""
- `Nho = ` $(@bind Nho Slider(5:5:100; show_value=true, default=50) )
"""

# ╔═╡ 828d31d2-6816-47fa-8789-fe370e1deb06
let
	# Find a reference solution using Tsit5():
	reference = solve_reference(ho_f, ho_u₀, ho_tstart, ho_tend;
								abstol=1e-14, reltol=1e-14);
	
	# Solve the problem using all our implemented RK methods
	soln_euler    = forward_euler(ho_f, ho_u₀, ho_tstart, ho_tend, Nho)
	soln_midpoint = midpoint(ho_f, ho_u₀, ho_tstart, ho_tend, Nho)
	soln_rk4      = rk4(ho_f, ho_u₀, ho_tstart, ho_tend, Nho)
	
	# Function to compute the pointwise error in the positions
	function compute_error(solution)
		t = solution.t
		u = solution.u
		[abs(reference(t[n])[1] - u[n][1]) .+ 1e-16 for n in 1:length(u)]
	end
	
	# HO plot x versus t
	xlims = (ho_tstart, ho_tend)
	p = plot(t -> reference(t)[1];  # plot only position
	         lw=2, label="reference", ylabel=L"x(t)", xlims,
		     title=L"\frac{d^2 x}{d t^2} = -k x", titlefontsize=12, ylims=(-3, 3),
			 legend=:bottomleft)
	plot!(p, soln_euler.t, [u[1] for u in soln_euler.u];
	      label="Forward Euler", lw=2, mark=:o, ls=:dash, c=2, markersize=3)
	plot!(p, soln_midpoint.t, [u[1] for u in soln_midpoint.u];
		  label="Midpoint", lw=2, mark=:o, ls=:dash, c=3, markersize=3)
	plot!(p, soln_rk4.t, [u[1] for u in soln_rk4.u];
		  label="RK4", lw=2, mark=:o, ls=:dash, c=4, markersize=3)

	# Total energy plots
	q = plot(; legend=false, xlabel=L"t", title="Error", titlefontsize=10,
	           yaxis=:log, ylims=(1e-5, 10), xlims)
	
	plot!(q, soln_euler.t, compute_error(soln_euler); label="Euler", lw=2, ls=:dash, c=2)
	plot!(q, soln_midpoint.t, compute_error(soln_midpoint); label="Midpoint", lw=2, ls=:dash, c=3)
	plot!(q, soln_rk4.t, compute_error(soln_rk4); label="RK4", lw=2, ls=:dash, c=4)

	plot(p, q; layout=grid(2, 1; heights=[0.75, 0.25]))
end

# ╔═╡ 366477ed-5532-44fa-97e0-fc0ab564fabe
md"""## Appendix"""

# ╔═╡ d3855794-f57d-4fc5-bf97-63ce44582139
function newton(f, df, xstart; maxiter=40, tol=1e-6)
	# f:  Function of which we seek the roots
	# df: Function, which evaluates its derivatives
	# xstart: Start of the iterations
	# maxiter: Maximal number of iterations
	# tol: Convergence tolerance

	history_x = [float(xstart)]
	history_r = empty(history_x)

	r = Inf  # Dummy to enter the while loop
	k = 0

	# Keep running the loop when the residual norm is beyond the tolerance
	# and we have not yet reached maxiter
	while norm(r) ≥ tol && k < maxiter
		k = k + 1
		
		# Pick most recent entry from history_x (i.e. current iterate)
		x = last(history_x)

		# Evaluate function, gradient and residual
		r = - f(x) / df(x)
		
		push!(history_r, r)      # Push residual and
		push!(history_x, x + r)  # next iterate to history
	end

	(; root=last(history_x), n_iter=k, history_x, history_r)
end

# ╔═╡ 7867834d-e5d5-4ae5-946b-8ea6ec6a91c9
function fixed_point_newton(g, dg, xstart; maxiter=40, tol=1e-6)
	f(x)  = g(x)  - x
	df(x) = dg(x) - 1
	(; root, n_iter, history_x, history_r) = newton(f, df, xstart; maxiter, tol)
	(; fixed_point=root, n_iter, history_x, history_r)
end

# ╔═╡ 58f37900-0e6a-4d00-adb3-ec82cbfbcf4c
function backward_euler(f, u₀, a, b, N; tol=1e-10)
	# f:  Derivative function
	# u₀: Initial value
	# a:  Start of the time interval
	# b:  End of the time interval
	
	
	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration
	u[1] = u₀
    for i in 1:N
		g(x)   = u[i] + h * f(x, t[i+1])
		dg(x)  = ForwardDiff.derivative(g, x)
		u[i+1] = fixed_point_newton(g, dg, u[i]; tol).fixed_point
    end

    (; t, u)
end

# ╔═╡ cceb1700-55fc-4e5a-a2f2-81137fad5c1f
let
	soln_euler    = forward_euler(dcy_f, dcy_u₀, dcy_tstart, dcy_tend, Nbw)
	soln_bw_euler = backward_euler(dcy_f, dcy_u₀, dcy_tstart, dcy_tend, Nbw)

	p = plot(u; lw=2, title="du/dx = $(round(-C; sigdigits=3))", titlefontsize=12, label="Reference", xlabel=L"t", ylabel=L"u(t)", xlims=(-0.1, 10.1), ylims=(-1, 25))

	plot!(p, soln_euler.t, soln_euler.u; mark=:o, c=2, lw=2, markersize=3, ls=:dash, label="Forward Euler")
	plot!(p, soln_bw_euler.t, soln_bw_euler.u; mark=:o, c=3, lw=2, markersize=3, ls=:dash, label="Backward Euler")
end

# ╔═╡ 403993d0-1613-4ac1-a527-9362993e28e6
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 500)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DifferentialEquations = "~7.16.1"
ForwardDiff = "~0.10.38"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.12"
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.62"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "073d0fedb2bdd79c9801f78082bc44acadfdc9ec"

[[deps.ADTypes]]
git-tree-sha1 = "e2478490447631aedba0823d4d7a80b2cc8cdb32"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.14.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.AlmostBlockDiagonals]]
deps = ["ConcreteStructs"]
git-tree-sha1 = "743abe5e5fe8cff96dad4123f263c0d8eee281c0"
uuid = "a95523ee-d6da-40b5-98cc-27bc505739d5"
version = "0.1.10"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "4e25216b8fea1908a0ce0f5d87368587899f75be"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "614c6aba1d562191d9832df2af24f594aa7ebf61"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.9.3"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "ArrayInterface", "BoundaryValueDiffEqAscher", "BoundaryValueDiffEqCore", "BoundaryValueDiffEqFIRK", "BoundaryValueDiffEqMIRK", "BoundaryValueDiffEqMIRKN", "BoundaryValueDiffEqShooting", "DiffEqBase", "FastClosures", "ForwardDiff", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "e3829b5aa0cb49348956c81b927b5edf64cdf6bf"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.16.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.BoundaryValueDiffEqAscher]]
deps = ["ADTypes", "AlmostBlockDiagonals", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a3ed69c1c0249a53622bd4435384c4e76ac547d9"
uuid = "7227322d-7511-4e07-9247-ad6ff830280e"
version = "1.5.0"

[[deps.BoundaryValueDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolveFirstOrder", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings"]
git-tree-sha1 = "832ade257129d0c222a53b66e2d7e6f5d937ae34"
uuid = "56b672f2-a5fe-4263-ab2d-da677488eb3a"
version = "1.8.0"

[[deps.BoundaryValueDiffEqFIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a92feb2cbb12c6c9adc4d3c4e7427709e9477540"
uuid = "85d9eb09-370e-4000-bb32-543851f73618"
version = "1.6.0"

[[deps.BoundaryValueDiffEqMIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "4cd74dc128326804f780ad6e18ec4886279293de"
uuid = "1a22d4ce-7765-49ea-b6f2-13c8438986a6"
version = "1.6.0"

[[deps.BoundaryValueDiffEqMIRKN]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "0db565e02c9784e254325b616a8dd6c0dfec7403"
uuid = "9255f1d6-53bf-473e-b6bd-23f1ff009da4"
version = "1.5.0"

[[deps.BoundaryValueDiffEqShooting]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "7429a95010c57e67bd10e52dd3f276db4d2abdeb"
uuid = "ed55bfe0-3725-4db6-871e-a1dc9f42a757"
version = "1.6.0"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "bfdafcc043eb34fe21a2dae769734fd918546d6b"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.1.3"
weakdeps = ["ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

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

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqRosenbrock", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack", "SymbolicIndexingInterface"]
git-tree-sha1 = "7123a01ba4ec2d4058bd14478afd5318c49ea6c1"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.52.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "TruncatedStacktraces"]
git-tree-sha1 = "e384a2cf3bb402e6dc66b1503ade22c7c1471c4d"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.167.2"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "Functors", "LinearAlgebra", "Markdown", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "10481b5d8e046df5280cf081117c8058a7ce7376"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.4.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "516d553f5deee7c55b2945b5edf05b6542837887"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.24.1"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "afdc7dfee475828b4f0286d63ffe66b97d7a3fa7"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.16.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "e41b6696c84291c4ad15f5f6eaf071b4dfbfda06"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.51"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "0b4190661e8a4e51a842070e7dd4fae440ddb7f4"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.118"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EnzymeCore]]
git-tree-sha1 = "0cdb7af5c39e92d78a0ee8d0a447d32f7593137e"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.8"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

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

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cae251c76f353e32d32d76fae2fea655eab652af"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.27.0"
weakdeps = ["StaticArrays"]

    [deps.ExponentialUtilities.extensions]
    ExponentialUtilitiesStaticArraysExt = "StaticArrays"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

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

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "3f03d94c71126b6cfe20d3cbcc41c5cd27e1c419"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.4"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "fd923962364b645f3719855c88f7074413a6ad92"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.0.2"

[[deps.FastPower]]
git-tree-sha1 = "df32f07f373f06260cd6af5371385b5ef85dd762"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.2"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "f089ab1f834470c525562030c8cfde4025d5e915"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.27.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

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

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

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

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "3169fd3440a02f35e549728b0890904cfd4ae58a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.1"

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

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

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

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "f2bdec5b4580414aee3178c8caa6e46c344c0bbc"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.14.3"
weakdeps = ["FastBroadcast"]

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b29d37ce30fa401a4563b18880ab91f979a29734"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.10"

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

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "866ce84b15e54d758c11946aacd4e5df0e60b7a3"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.1"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

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

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "2bbbdcce6d80a4aed929365d0d97b15b264bb9e7"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.7.2"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"

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

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "16a726dba99685d9e94c8d0a8f655383121fc608"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.0.1"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

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

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "453de0fc2be3d11b9b93ca4d0fddd91196dcf1ed"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.5"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "b14c7be6046e7d48e9063a0053f95ee0fc954176"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.9.1"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "7ae7322d658544bd8f6b24a1a0374a6b4ac1fc7e"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.5.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "f8ece81557f7e42879f017fa089e5283988a5f67"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.5.2"
weakdeps = ["BandedMatrices", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearSolve", "SparseArrays", "SparseMatrixColorings"]

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseDiffEqBaseExt = "DiffEqBase"
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "53e9df7c663c5b9ee5443ce4591f145143095c23"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.3.1"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "61341153cec9ab307b6bae09f2661ebbd34cc1f9"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.2.1"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "DiffEqBase", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "7e94679f6989398b90a4d665f6b0d16c0263a451"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.1.1"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

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
version = "0.8.1+4"

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

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "31b3b1b8e83ef9f1d50d74f1dd5f19a37a304a1f"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.12.0"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "97037e44313e33cd29e8b08e2ec82dd157f866ae"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.93.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "82f78099ecf4e0fa53545811318520d87e7fe0b8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.2.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "6c8114a2f5a7649367241a99724cc7e79f7d2d40"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.4.0"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "ba84fa52a477a537213b7d0a6a83ba9b2f3aaa00"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.22.0"
weakdeps = ["EnzymeCore"]

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport"]
git-tree-sha1 = "835c06684b6ff1b8904ceae4d18cc8fe45b9a7cc"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.3.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "dafdc71c19e05744876302babdc91149829c2e0a"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.6.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "8d2ab84d7fabdfde995e5f567361f238069497f5"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.4.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "80a636aac325c546b04e3bf20f0c80eaa0173dd4"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.5.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "7d2c82c13a634f7400a3f398d33f1354ab38a090"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.10.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a7cc74d3433db98e59dc3d58bc28174c6c290adf"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.1.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "925a91583d1ab84f1f0fea121be1abf1179c5926"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.1.1"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "103e017ff186ac39d731904045781c9bacfca2b0"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.1.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "095bab73a3ff185e9ef971fc42ecc93c7824e589"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.3.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "0f81a77ede3da0dc714ea61e81c76b25db4ab87a"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.1.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "d4bb32e09d6b68ce2eb45fb81001eab46f60717a"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.2.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "52ec7081e65291fa5c19749312df0818db2fa1bc"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.3.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "dbc4ca97a9051595d86a0d0fc31ded2de9cdfd7e"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.6.0"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "ef44754f10e0dfb9bb55ded382afed44cd94ab57"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.1.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "StaticArrays"]
git-tree-sha1 = "f74b27b8b811a83d77a9cad6293e793ab0804cdc"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.3.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport"]
git-tree-sha1 = "da525d277962a1b76102c79f30cb0c31e13fe5b9"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.1.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "332f9d17d0229218f66a73492162267359ba85e9"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.1.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "41c09d9c20877546490f907d8dffdd52690dd65f"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.1.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "baa4a9b4380b2fb65f1e2b4ec01d3bd019a6dcea"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.9.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "b3a7e3a2f355d837c823b435630f035aef446b45"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.3.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "b97b9691437063ae7d9b1ef130e8b0d81415116f"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.2.1"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "111c23b68ad644b47e38242af920d5805c7bedb1"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.3.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "1b0d894c880e25f7d0b022d7257638cf8ce5b311"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.1.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "a13d59a2d6cfb6a3332a7782638ca6e1cb6ca688"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.3.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "96552f7d4619fabab4038a29ed37dd55e9eb513a"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.1.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "81d7841e73e385b9925d5c8e4427f2adcdda55db"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.1.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "0e1340b5d98971513bddaa6bbed470670cebbbfe"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.34"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3b31172c032a1def20c98dae3f2cdc9d10e3b561"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

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

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6d38fea02d983051776a856b7df75b30cf9a3c1f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.16"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "4406f9a118bfcf362290d755fcb46c0c4894beae"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.26"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

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

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

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

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "112c876cee36a5784df19098b55db2b238afc36a"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.31.2"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

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

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "5cf59106f9b47014c58c5053a1ce09c0a2e0333c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.7.3"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "c4ce89e19f2a220e34c0f377fb41516b642ec899"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.83.1"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "15634a7c06849c6871a3a2391346d48ef3ba0fbe"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.2"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "1c4b7f6c3e14e6de0af66e66b86d525cae10ecb4"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.13"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
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

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "a65e81385d9c45c6abe49d7676d71b07b5b6bbc1"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.2.1"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveDiffEqBaseExt = "DiffEqBase"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

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

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "15dd194e46a5e74b6f7f361e9eb3ed869617ccd5"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.16"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerDataInterpolationsExt = "DataInterpolations"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "59bad850b1fc622051bf80a2be86c95b487e0243"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.24.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "0582fd1410a01a667a2a2a79cdc98a7c478d11d8"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.18"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

    [deps.StaticArrayInterface.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "35b09e80be285516e52c9054792c884b9216ae3c"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.4.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NonlinearSolveBase", "Reexport", "SciMLBase"]
git-tree-sha1 = "66a028f9a2bb44d0f6de0814a2b9840af548143a"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.5.0"

[[deps.StochasticDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "e48c760e73ebf9c0726892625a072cca3b96ea96"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.75.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "c135b599cec3558be36eaf86ab1ce7e259ef9534"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.27.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrettyTables", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "7530e17b6ac652b009966f8ad53371a4ffd273f2"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.39"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

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

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "18ad3613e129312fe67789a71720c3747e598a61"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.3"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f57facfd1be61c42321765d3551b3df50f7e09f6"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.28"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

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
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
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

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

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
# ╟─ba9b6172-0234-442c-baaa-876b12f689bd
# ╠═8197b6f0-00b7-11ef-2142-4b7cbbaefd90
# ╟─d8406b01-e36f-4953-a5af-cd563005c2a1
# ╟─692aa6e2-21a9-4dad-a719-9d8ad88bd467
# ╟─ebabe949-5492-4b8b-959d-05c5728da043
# ╟─05ef8174-e9a6-4280-8640-08d74635fba2
# ╠═31332dd6-1ef2-4181-a638-35980f470552
# ╟─1e65223e-9a12-4a4c-8b5e-31bfe4f813ec
# ╟─64fe575e-0d47-4949-a9e8-2056ddee45df
# ╟─1dbe5d72-17d9-4f20-b365-ad913cd607c3
# ╟─6c2f1cd1-e79a-4790-9871-0a8825b5c02e
# ╠═8f159939-4bf4-4151-bd3e-a84c5e1493db
# ╟─9b71f2b6-34bb-416a-947d-6e009f565bf4
# ╟─3f59bd74-1bb0-41df-bea2-bf5990a2b89a
# ╠═9cb5fd8e-d923-434e-b543-a0c987a5b5a7
# ╟─5d2886f2-87a5-43ae-9064-448478f957ca
# ╠═82321d1e-50b9-4362-afc5-d256a506bed7
# ╟─6b4ba69d-0a70-4c66-b4bc-63b717348471
# ╠═7a0fa3cd-a270-4947-b7ba-d30110139795
# ╟─c94d264b-ad1b-47fd-8b58-fd434cfe0f11
# ╟─3dc8b41f-3ba0-4c0a-b017-731fd0958920
# ╟─bc840e8d-77d9-4a4a-a908-3e9b4a8d253b
# ╟─66affc7b-f0b0-48f4-b8bd-a1ff16d0357e
# ╠═f60ecd21-fda2-4925-ba86-ac5ad5b7c7d0
# ╟─41cc0b76-5e30-41dd-a648-5491e1b4bd22
# ╟─c641cdf2-aabf-45e7-bb12-83d73a93b5b7
# ╠═9b0581c3-6fa8-44d8-8ebb-d2fda2841230
# ╟─f0164280-14fd-48a5-9ae8-8b31f4580da8
# ╠═6052f3cf-8a1e-4755-8724-a3f885702d62
# ╠═90a2176f-7986-4140-8e34-555993070324
# ╠═15e135e3-e95f-4494-8a20-2c571c4019f2
# ╟─14edb601-ab10-4228-81a1-d27bb16ea202
# ╠═50007ba8-8160-48c1-b5c1-27a9332c1d58
# ╟─fb38ada2-26ce-42e9-bfdc-e4c677f8c501
# ╟─65807ddc-0c8b-4b88-bbc2-4ff3c7a1d2f8
# ╟─98f64da0-0d01-4c7e-8782-3d97dcf4c950
# ╟─ecb868c3-44a2-4f95-bafe-b3324770ee6b
# ╟─667fbe5a-dbee-43a3-bb93-0160b61cf31e
# ╟─46f8284f-7751-4c1e-a617-386d86e7f4d0
# ╟─9385fd22-95f4-42b6-83ca-d4abbfe83a30
# ╟─266c5a9a-8cd5-4dc0-a7de-492d64306aab
# ╟─d7150f04-e3e8-41f7-a173-9d458f72cdde
# ╟─6d73f60b-0ccf-492f-9700-09d0bbf15124
# ╟─80748b7b-2546-4c4c-b60d-fdf7df98c8a7
# ╟─5ac4eca2-0372-4b2d-add6-2347ce9461a3
# ╟─74a4bb06-d784-423c-9658-3cb9e4f545ae
# ╟─285cbf5b-7097-47c0-950c-0c84d69a65bd
# ╟─20096950-17ef-4d60-8059-dc7f00265349
# ╟─38a5d61d-aff7-44dd-90e2-e690cd2643f7
# ╟─04b3c441-a527-4b15-99e1-042479c4fdf2
# ╟─4d20e43b-ba3d-4eab-adfc-f3a7c6feeb0e
# ╟─3430cee3-d5ac-4f7d-99ce-713008cd3553
# ╠═76978387-3b7b-4644-a27e-246f84186eb9
# ╟─87345743-d3a5-4ec4-bf66-15e6209a729a
# ╟─2c3ad400-acdd-41bb-ac40-16cdc4b6fb07
# ╟─262a2a67-3fd3-497c-b56c-71fe26abffc1
# ╟─850928f3-ab31-42c6-a51a-0dca28175bdf
# ╠═5f7562f6-4329-4c1e-b5d4-74c9b1980346
# ╟─f0f756fc-f7fe-459c-8419-4e60c741fce5
# ╟─f9e8c1ec-3bc6-40d3-b7e3-f465bd53ecab
# ╟─bdb783bc-1648-43f0-b706-462d2a5d6dcf
# ╟─99922f12-34b6-4a9d-b3af-2b61708725fa
# ╟─0bd09a85-2fe9-491d-b6f9-f2d74872e195
# ╟─f5c8f93c-c09a-4c4f-b6c4-18b12ff56878
# ╟─591d6a3a-fc71-4bde-b499-db10637e0945
# ╟─85b121ef-9f4c-483e-a04b-a49c1c6c4624
# ╠═5348acce-d58b-4b33-97ad-44f44717e1fe
# ╟─b6022950-3711-470c-b7fe-ae02232a0749
# ╠═f2bde123-95b7-4370-9914-577985c4efeb
# ╟─a09b8d7b-9687-487d-940d-93c6957b027c
# ╠═af8f49e2-a876-4bb9-8036-46c1f510f369
# ╟─d57e6a62-cc98-4dbc-b471-b21448c0343e
# ╟─1f1af0d4-3c5d-4e38-9eda-c307100bd61f
# ╠═e4fa9b2d-8a81-4439-84ff-8529a7421ae9
# ╟─60a34642-3824-4a3f-a11c-e76039c91e4a
# ╟─2263c105-7d59-43a9-9196-b4762a53b1b0
# ╟─5075f5e4-af28-4a93-b8e6-d71f96cf2b89
# ╠═58f37900-0e6a-4d00-adb3-ec82cbfbcf4c
# ╟─37bb847f-aa6b-4efe-8538-e580159d2892
# ╠═cceb1700-55fc-4e5a-a2f2-81137fad5c1f
# ╟─e031fb26-dc5e-44e9-ba4c-a642c62b7930
# ╟─35b8a591-62ea-4836-aa02-e1aaa584621d
# ╠═6b5ac847-ea65-4d44-9379-5cc755908521
# ╟─d8e1eb69-9f8f-47c3-a3f8-bf7f725959de
# ╟─d09ee1f0-c9ca-4888-9573-6ca0c105dccd
# ╠═f35f9be4-74ae-4a05-ab6a-55aa799d98a2
# ╠═01e9b24f-7819-4e72-9463-2a6aca835d22
# ╟─3834116a-be19-4690-a4be-542d7824f410
# ╠═612e93db-322b-4b7c-a13e-cfd18058178a
# ╠═94a39f52-0653-4f01-a08a-ee9bc887da7d
# ╟─823685fd-221a-40c3-855a-20825fda7bf3
# ╟─e195e04c-fd23-4f92-8241-235d69404d8e
# ╠═828d31d2-6816-47fa-8789-fe370e1deb06
# ╟─366477ed-5532-44fa-97e0-fc0ab564fabe
# ╠═7867834d-e5d5-4ae5-946b-8ea6ec6a91c9
# ╠═d3855794-f57d-4fc5-bf97-63ce44582139
# ╠═403993d0-1613-4ac1-a527-9362993e28e6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
