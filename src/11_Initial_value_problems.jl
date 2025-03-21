### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
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

# ╔═╡ e9ef08a6-e823-4f2b-9fcc-e3646171cefb
TODO(md"Start with example and introduce general definition thereafter.")

# ╔═╡ 6e9b0961-04b7-436b-90f9-9804e06d10f8
md"""
In computational science we are often faced with quantitites that change continuously in space or time. For example the temperature profile of a hot body or the population of an animal species over time. Such problems are often modelled by differential equations. If there is only a single indepentent variable $t$ we call the model an **ordinary differential equation**. The usual setup is that for initial value $t = 0$ we have some knowledge about our problem, e.g. by performing a measurement, and the key question is thus how the situation evolves for $t > 0$. Mathematically we call this an **initial value problem**, for example:

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

Often (but not always) $t$ plays the role of time and (1) thus models the time-dependence of a quantity $u$. Let us consider an example:
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

# ╔═╡ d811efdc-c133-4146-a507-0786cbde68ad
TODO(md"Avoid introducing DifferentialEquations.jl and the ODEProblem object.")

# ╔═╡ c2359767-5256-4403-b6fe-0ded5ba41efa
md"""
## Numerical solution in DifferentialEquations.jl
For such simple examples analytic solutions can still be found with a little practice. However, for more involved cases analytical solutions can be more tricky or impossible and we will discuss numerical techniques to solve such equations. 
Before we discuss a few simple methods in detail, we first consider the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) Julia package, which provides a range of production-grade numerical solvers for ODE problems and will serve us as the reference method in our discussion.

As an example we consider the more involved initial value problem
```math
\tag{2}
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= \sin\left[(u+t)^2\right] && 0 ≤ t ≤ 4 \\
u(0) &= -1,
\end{aligned}\right.
```
for which an analytical solution is not so easy to obtain.
To solve it using DifferentialEquations
we first construct an `ODEProblem` object,
which is essentially a Julia representation of (2):
"""

# ╔═╡ e3c3abaf-0a17-4d84-945d-a1bf9bdfc0de
begin
	f = (u, p, t) -> sin((t + u)^2)  # defines du/dt, must include p argument
	u₀ = -1.0                        # initial value
	tspan = (0.0, 4.0)               # t interval 
	ivp = ODEProblem(f, u₀, tspan)
end

# ╔═╡ 49d045f8-a510-4f7b-a8c9-21ed3b4e0e6c
md"""
A good default solver for such problems is the `Tsit5` method, which can be employed as follows:
"""

# ╔═╡ 0f37dc8a-fb6b-424c-84aa-0e46ef89437b
sol = solve(ivp, Tsit5());

# ╔═╡ 5d2886f2-87a5-43ae-9064-448478f957ca
md"The obtained solution can easily be visualised:"

# ╔═╡ 82321d1e-50b9-4362-afc5-d256a506bed7
plot(sol, label="solution", lw=2, xlabel="t", ylabel=L"u(t)",
	 title=L"\frac{du}{dt} = \sin((t+u)^2)")

# ╔═╡ 6b4ba69d-0a70-4c66-b4bc-63b717348471
md"""
## Optional: Existence and uniqueness of solutions

Before we discuss some numerical techniques for solving initial value problems
in detail, a quick reminder that in general such an IVP can by itself not be solvable or admit more than one solution.

First consider the case where the **existence** of a solution to problem (1) is not guaranteed for all $t$ with $a ≤ t ≤ b$. For example consider the problem
```math
\frac{d u(t)}{dt} = \left(u(t)\right)^2,\, 0≤t≤2, \qquad u(0) = 1,
```
which has solution $u(t) = \frac{1}{1-t}$. This solution only exists for $t<1$.
When attempting a numerical solution beyond $t=1$ of such a problem we are faced with problems:
"""

# ╔═╡ 7a0fa3cd-a270-4947-b7ba-d30110139795
let
	f  = (u, p, t) -> u^2
	u₀ = 1.0
	tspan = (0, 10)
	
	ivp = ODEProblem(f, u₀, tspan)
	sol = solve(ivp,Tsit5())
	
	plot(sol; lw=2, label="", xlabel=L"t", yaxis=:log, ylabel=L"u(t)", xlims=(0, 2))
end

# ╔═╡ 0025be9a-6320-415a-b2ac-70f24bb0d12b
md"""
A second issue is that the solution may not necessarily be **unique**. Consider for example the equation
```math
\frac{d u(t)}{dt} = \sqrt{u(t)}, \, t>0, \qquad u(0) = 0,
```
which has the two solutions $u_1(t) = 0$ and $u_2(t) = \frac14 t^2$. 

Both cases cause additional difficulties to numerical techniques, which we do not want to discuss here. Instead we thus recall the following result ensuring the existence of a unique solution to (1):

!!! info "Theorem 1: Existence and uniquens of first-order ODEs"
	Given $f : \mathbb{R} \times \mathbb{R} \to \mathbb{R}$ a continuous function
	with respect to both arguments and **Lipschitz-continuous** with respect to its second argument, i.e. there exists a constant $L>0$ such that
	```math
	| f(t, x) - f(t, y)| \leq L |x-y|, \qquad \forall x,y \in \mathbb{R}, \, \forall t \text{ with } a ≤ t ≤ b
	```
	Then the ODE problem (1) has a unique solution $u(t)$ defined for all
	$a ≤ t ≤ b$ and $u$ is moreover continuously differentiable.

For the rest of this chaper we will always assume this condition to be satisfied.
"""

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
function forward_euler(ivp, N)
	a, b = ivp.tspan   # Extract time interval [a, b]
	u₀   = ivp.u0      # Initial condition
	f    = ivp.f       # Derivatve function f
	p    = ivp.p       # Further parameters (not used in our implementation)

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration
	u[1] = u₀
    for i in 1:N
		u[i+1] = u[i] + h * f(u[i], p, t[i])
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
which we solved previously using `DifferentialEquations`. The `ivp` object already contains a setup of this initial value problem
"""

# ╔═╡ 5374c9db-6b90-4775-a3f0-4c226111d7d5
ivp

# ╔═╡ 14edb601-ab10-4228-81a1-d27bb16ea202
md"""
and we can directly proceed to solving it:
"""

# ╔═╡ fb38ada2-26ce-42e9-bfdc-e4c677f8c501
md"and plot the result:"

# ╔═╡ 98f64da0-0d01-4c7e-8782-3d97dcf4c950
md"- Number of subintervals: `Neuler = ` $(@bind Neuler Slider(5:1:50; show_value=true, default=20)) "

# ╔═╡ 50007ba8-8160-48c1-b5c1-27a9332c1d58
res_euler = forward_euler(ivp, Neuler);

# ╔═╡ 65807ddc-0c8b-4b88-bbc2-4ff3c7a1d2f8
let
	res_euler = forward_euler(ivp, Neuler)

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
	res_euler = forward_euler(ivp, 7)
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
	u_exact = solve(ivp, Tsit5(); reltol=1e-14, abstol=1e-14)

	Ns = [ round(Int, 5 * 10^k) for k in 0:0.5:3 ]
	errors = Float64[]
	for N in Ns
		res_euler = forward_euler(ivp, N)
		u = res_euler.u
		t = res_euler.t

		# Compute global errors in each node tₙ
		global_error = [abs(u_exact(t[n]) - u[n]) for n in 1:length(u)]

		# Push the largest error to the array for plotting
		push!(errors, maximum(global_error))
	end
	
	plot(Ns, errors; mark=:o, label="Errors", lw=2, xaxis=:log, yaxis=:log,
		xlabel=L"N
		", ylabel="Largest global error",
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
function midpoint(ivp, N)
	a, b = ivp.tspan   # Extract time interval [a, b]
	u₀   = ivp.u0      # Initial condition
	f    = ivp.f       # Derivatve function f
	p    = ivp.p       # Further parameters (not used in our implementation)

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration ... this is what changes over forward_euler
	u[1] = u₀
    for i in 1:N
		uhalf  = u[i] + h/2 * f(u[i],  p, t[i]      )
		u[i+1] = u[i] + h   * f(uhalf, p, t[i] + h/2)
    end

    (; t, u)
end

# ╔═╡ 87345743-d3a5-4ec4-bf66-15e6209a729a
md"""In comparison with Forward Euler we notice this method to be clearly more accurate for the case of rather small number of timesteps:"""

# ╔═╡ 2c3ad400-acdd-41bb-ac40-16cdc4b6fb07
let
	N = 9
	res_euler    = forward_euler(ivp, N)
	res_midpoint = midpoint(ivp, N)
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
function rk4(ivp, N)
	a, b = ivp.tspan   # Extract time interval [a, b]
	u₀   = ivp.u0      # Initial condition
	f    = ivp.f       # Derivatve function f
	p    = ivp.p       # Further parameters (not used in our implementation)

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration ... this is what changes over forward_euler
	u[1] = u₀
    for i in 1:N
		v₁ = h * f(u[i],        p, t[i]      )
		v₂ = h * f(u[i] + v₁/2, p, t[i] + h/2)
		v₃ = h * f(u[i] + v₂/2, p, t[i] + h/2)
		v₄ = h * f(u[i] + v₃,   p, t[i] + h  )

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
	res_euler    = forward_euler(ivp, N)
	res_midpoint = midpoint(ivp, N)
	res_rk4      = rk4(ivp, N)
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
	u_exact = solve(ivp, Tsit5(); reltol=1e-14, abstol=1e-14)

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
		result_euler = forward_euler(ivp, N)
		push!(errors_euler, global_error(result_euler.t, result_euler.u, u_exact))
		
		result_midpoint = midpoint(ivp, N)
		push!(errors_midpoint, global_error(result_midpoint.t, result_midpoint.u, u_exact))
		
		res_rk4 = rk4(ivp, N)
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
- `logC = `: $(@bind logC Slider(-4:0.2:1; show_value=true, default=0.8))
- `Ndecay = `: $(@bind Ndecay Slider(15:5:50; show_value=true, default=15))
"""

# ╔═╡ 85b121ef-9f4c-483e-a04b-a49c1c6c4624
md"""
## Stability and implicit methods

Let us consider the innocent looking initial value problem
```math
\left\{
\begin{aligned}
\frac{d u(t)}{d t} &= -10^\textrm{logC} && t > 0 \\
u(0) &= 10.
\end{aligned}\right.
```
This model governs for example the decay of a radioactive species with initial concentration $10$ over time. The rate of decay is $10^\textrm{logC}$ ≈  $(round(10^logC; sigdigits=3)).

By simple integration we find the exact solution of this problem as:
"""

# ╔═╡ 5348acce-d58b-4b33-97ad-44f44717e1fe
u(t) = 10 * exp(-10^logC * t)

# ╔═╡ f2bde123-95b7-4370-9914-577985c4efeb
plot(u; xlims=(0, 10), lw=2, title="du/dx = $(round(-10^logC; sigdigits=3))", titlefontsize=12, label="Reference", xlabel=L"t", ylabel=L"u(t)")

# ╔═╡ af8f49e2-a876-4bb9-8036-46c1f510f369
decay = let
	f  = (u, p, t) -> -10^logC * u
	u₀ = 10.0
	tspan = (0.0, 10.0)
	decay = ODEProblem(f, u₀, tspan)
end

# ╔═╡ e4fa9b2d-8a81-4439-84ff-8529a7421ae9
let
	h = 10 / Ndecay  # Stepsize of our numerical methods
	
	soln_euler    = forward_euler(decay, Ndecay);
	soln_midpoint = midpoint(decay, Ndecay);
	soln_rk4      = rk4(decay, Ndecay);

	p = plot(u; lw=2, title="du/dx = $(round(-10^logC; sigdigits=3))", titlefontsize=12, label="Reference", xlabel=L"t", ylabel=L"u(t)", xlims=(-0.1, 10.1))

	if h ≥ 2/10^logC  # See discussion below to understand
	                  # why this formula appears here
		ylims!(p, (-1_000, 10_000))
	end

	plot!(p, soln_euler.t, soln_euler.u; mark=:o, c=2, lw=2, markersize=3, ls=:dash, label="Forward Euler")
	plot!(p, soln_midpoint.t, soln_midpoint.u; mark=:o, c=3, lw=2, markersize=3, ls=:dash, label="Midpoint")
	plot!(p, soln_rk4.t, soln_rk4.u; mark=:o, c=4, lw=2, markersize=3, ls=:dash, label="RK4")
end

# ╔═╡ 60a34642-3824-4a3f-a11c-e76039c91e4a
md"""We notice for too few nodal points or too large decay constants $10^\textrm{logC}$ our numerical methods **all fail** to recover the decay behaviour.
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
	\lim_{n\to\infty} u^{(n)} = u_\ast \qquad \forall u_0.
	```

	A method is further called **unconditionally absolutely stable** if it is stable for all $h > 0$.

In the case of our decay problem, one can show for example that the forward Euler method is absolutely stable if
```math
h < \frac{2}{10^{\textrm{logC}}}
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

First we setup the HO model as an `ODEProblem`:
"""

# ╔═╡ f35f9be4-74ae-4a05-ab6a-55aa799d98a2
k = 2  # Force constant

# ╔═╡ 01e9b24f-7819-4e72-9463-2a6aca835d22
ho = let
	# To setup f we assume u is a vector of size 2 and return a vector of size 2
	f = (u, p, t) -> [     u[2];
				      -k * u[1]]

	# As the initial value again we supply a vector of size 2
	u₀ = [0.0;
		  1.0]

	# And we are interested in the behaviour from 0 to tend (defined below)
	tspan = (0.0, 10.0)

	ho = ODEProblem(f, u₀, tspan)
end

# ╔═╡ 3834116a-be19-4690-a4be-542d7824f410
md"""Then running the dynamics just takes a call to `forward_euler` as before:"""

# ╔═╡ 612e93db-322b-4b7c-a13e-cfd18058178a
ho_euler = forward_euler(ho, 50);

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
	reference = solve(ho, Tsit5(); abstol=1e-14, reltol=1e-14);
	
	# Solve the problem using all our implemented RK methods
	soln_euler    = forward_euler(ho, Nho)
	soln_midpoint = midpoint(ho, Nho)
	soln_rk4      = rk4(ho, Nho)
	
	# Function to compute the pointwise error in the positions
	function compute_error(solution)
		t = solution.t
		u = solution.u
		[abs(reference(t[n])[1] - u[n][1]) .+ 1e-16 for n in 1:length(u)]
	end
	
	# HO plot x versus t
	p = plot(t -> reference(t)[1];  # plot only position
	         lw=2, label="reference", ylabel=L"x(t)", xlims=ho.tspan,
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
	           yaxis=:log, ylims=(1e-5, 10), xlims=ho.tspan)
	
	plot!(q, soln_euler.t, compute_error(soln_euler); label="Euler", lw=2, ls=:dash, c=2)
	plot!(q, soln_midpoint.t, compute_error(soln_midpoint); label="Midpoint", lw=2, ls=:dash, c=3)
	plot!(q, soln_rk4.t, compute_error(soln_rk4); label="RK4", lw=2, ls=:dash, c=4)

	plot(p, q; layout=grid(2, 1; heights=[0.75, 0.25]))
end

# ╔═╡ 366477ed-5532-44fa-97e0-fc0ab564fabe
md"""## Appendix"""

# ╔═╡ 448d1e06-c4b2-4954-978e-00c2ed274ae0
function fixed_point_iterations(g, xstart; tol=1e-6, maxiter=100)
	# g:      Fixed-point function
	# xstart: Initial guess
	# tol:    Tolerance

	history_x = [float(xstart)]
	history_r = empty(history_x)
	
	r⁽ᵏ⁾ = Inf  # For initial pass in while loop
	k  = 0
	while k < maxiter && abs(r⁽ᵏ⁾) ≥ tol
		x⁽ᵏ⁾ = last(history_x)  # Take out the most recent point from the history
		x⁽ᵏ⁺¹⁾ = g(x⁽ᵏ⁾)
		r⁽ᵏ⁾   = x⁽ᵏ⁺¹⁾ - x⁽ᵏ⁾
		push!(history_r, r⁽ᵏ⁾)

		k = k + 1  # Update k
		push!(history_x, x⁽ᵏ⁺¹⁾)  # Push next point to the history
	end

	# Return results as a named tuple
	(; fixed_point=last(history_x), residual=r⁽ᵏ⁾, n_iter=k, history_x, history_r)
end

# ╔═╡ be538898-0493-45dc-b2ac-5d079981ef3b
function newton(f, df, xstart; maxiter=100, tol=1e-6)
	# f:  Function of which we seek the roots
	# df: Function, which evaluates its derivatives
	# xstart: Start of the iterations
	# maxiter: Maximal number of iterations
	# tol: Convergence tolerance

	# Define the fixed-point function g_Newton using f and df
	g_Newton(x) = x - f(x) / df(x)

	# Solve for its fixed point:
	fixed_point_iterations(g_Newton, xstart; tol, maxiter)
end

# ╔═╡ 7867834d-e5d5-4ae5-946b-8ea6ec6a91c9
function fixed_point_newton(g, dg, xstart; maxiter=40, tol=1e-6)
	f(x)  = g(x)  - x
	df(x) = dg(x) - 1
	newton(f, df, xstart; maxiter, tol)
end

# ╔═╡ 58f37900-0e6a-4d00-adb3-ec82cbfbcf4c
function backward_euler(ivp, N; tol=1e-10)
	a, b = ivp.tspan   # Extract time interval [a, b]
	u₀   = ivp.u0      # Initial condition
	f    = ivp.f       # Derivatve function f
	p    = ivp.p       # Further parameters (not used in our implementation)

	# Discrete time nodes to consider:
    h = (b - a) / N
    t = [a + i * h for i in 0:N]

	# Setup output vector with all elements initialised to u₀
	u = [copy(u₀) for i in 0:N]

	# Time integration
	u[1] = u₀
    for i in 1:N
		g(x)   = u[i] + h * f(x, p, t[i+1])
		dg(x)  = ForwardDiff.derivative(g, x)
		u[i+1] = fixed_point_newton(g, dg, u[i]; tol).fixed_point
    end

    (; t, u)
end

# ╔═╡ cceb1700-55fc-4e5a-a2f2-81137fad5c1f
let
	soln_euler    = forward_euler(decay,  Nbw)
	soln_bw_euler = backward_euler(decay, Nbw)

	p = plot(u; lw=2, title="du/dx = $(round(-10^logC; sigdigits=3))", titlefontsize=12, label="Reference", xlabel=L"t", ylabel=L"u(t)", xlims=(-0.1, 10.1), ylims=(-1, 25))

	plot!(p, soln_euler.t, soln_euler.u; mark=:o, c=2, lw=2, markersize=3, ls=:dash, label="Forward Euler")
	plot!(p, soln_bw_euler.t, soln_bw_euler.u; mark=:o, c=3, lw=2, markersize=3, ls=:dash, label="Backward Euler")
end

# ╔═╡ 403993d0-1613-4ac1-a527-9362993e28e6
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 570)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DifferentialEquations = "~7.13.0"
ForwardDiff = "~0.10.36"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Plots = "~1.40.4"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "f21ad0f661f44b93652d1dfece925c44e4b08838"

[[deps.ADTypes]]
git-tree-sha1 = "016833eb52ba2d6bea9fcb50ca295980e728ee24"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "c0d491ef0b135fd7d63cbc6404286bc633329425"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.36"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "133a240faec6e074e07c31ee75619c90544179cf"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.10.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "33207a8be6267bc389d0701e97a9bce6a4de68eb"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.9.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "30b7ea34abc4fe816eb1a5f434a43da804836163"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.7.0"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "005b55fa2eebaa4d7bf3cfb8097807f47116175f"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.7.1"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "585a387a490f1c4bd88be67eea15b93da5e85db7"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "4b270d6465eb21ae89b732182c20dc165f8bf9f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.25.0"

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
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
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
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
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

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "5959ae76ebd198f70e9af81153644543da0cfaf2"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.47.3"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces"]
git-tree-sha1 = "d520d3007de793f4fca16c77a25a9774ebe4ad6d"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.150.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NonlinearSolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "c959cfd2657d16beada157a74d52269e8556500e"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "3.6.2"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "65cbbe1450ced323b4b17228ccd96349d96795a7"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.21.0"

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
git-tree-sha1 = "81042254a307980b8ab5b67033aca26c2e157ebb"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.13.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "22c595ca4146c07b16bcf9c8bea86f731f7109d2"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.108"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "18394bc78ac2814ff38fe5e0c9dc2cd171e2810c"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.7.2"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "9dc913faf8552fd09b92a0d7fcc25f1d5609d795"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.1"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "a6e756a880fc419c8b41592010aebe6a5ce09136"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.8"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f4102aab9c7df8691ed09f9c42e34f5ab5458ab9"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "2de436b72c3422940cbe1367611d137008af7ec3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.23.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

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
deps = ["LinearAlgebra"]
git-tree-sha1 = "d3e63d9fa13f8eaa2f06f64949e2afc593ff52c2"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.10"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ddda044ca260ee324c5fc07edb6d7cf3f0b9c350"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "278e5e0f820178e8a26df3184fcb2280717c79b1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.5+0"

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
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "4f2b57488ac7ee16124396de4f2bbdd51b2602ad"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.11.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

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
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be50fe8df3acbffa0274a744f1a99d29c45a57f4"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "e7cbed5032c4c397a6ac23d1493f3289e01231c4"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.14"
weakdeps = ["Dates"]

    [deps.InverseFunctions.extensions]
    DatesExt = "Dates"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e9648d90370e2d0317f9518c9c6e0841db54a90b"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.31"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "ed08d89318be7d625613f3c435d1f6678fba4850"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.11.1"
weakdeps = ["FastBroadcast"]

    [deps.JumpProcesses.extensions]
    JumpProcessFastBroadcastExt = "FastBroadcast"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "267dad6b4b7b5d529c76d40ff48d33f7e94cb834"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays"]
git-tree-sha1 = "35079a6a869eecace778bcda8641f9a54ca3a828"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "1.10.0"
weakdeps = ["StaticArrays"]

    [deps.LazyArrays.extensions]
    LazyArraysStaticArraysExt = "StaticArrays"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "efd815eaa56c0ffdf86581df5aaefb7e901323a0"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.30.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = ["Enzyme", "EnzymeCore"]
    LinearSolveFastAlmostBandedMatricesExt = ["FastAlmostBandedMatrices"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8f6786d8b2b3248d79db3ad359ce95382d5a6df8"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.170"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "c6a36b22d2cca0e1a903f00f600991f97bf5f426"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.6"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "80b2833b56d466b3858d565adcd16a4a05f2089b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "6731e0574fa5ee21c02733e397beb133df90de35"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "2.2.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "b1f2f92feb0bc201e91c155ef575bcc7d9cc3526"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.2"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

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

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "dc0d78eeed89323526203b8a11a4fa6cdbe25cd6"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.11.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSymbolicsExt = "Symbolics"
    NonlinearSolveZygoteExt = "Zygote"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.OffsetArrays]]
git-tree-sha1 = "e64b4f5ea6b7389f6f046d13d4896a8f9c1ba71e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "4cf03bfe9c6159f66b57cda85f169cd0eff0818d"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.76.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

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

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

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
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "b3e2bae88cf07baf0a051fe09666b8ef97aefe93"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.14"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "a660e9daab5db07adf3dedfe09b435cc530d855e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.21"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

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

[[deps.PtrArrays]]
git-tree-sha1 = "077664975d750757f30e739c870fbbdc01db7913"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.1.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "758bc86b90e9fee2edc4af2a750b0d3f2d5c02c5"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.19.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

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
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "12aa2d7593df490c407a3bbd8b86b8b515017f3e"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.14"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

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

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "265f1a7a804d8093fa0b17e33e45373a77e56ca5"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.38.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools", "Setfield", "SparseArrays", "StaticArraysCore"]
git-tree-sha1 = "10499f619ef6e890f3f4a38914481cc868689cd5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.8"

[[deps.SciMLStructures]]
git-tree-sha1 = "d778a74df2f64059c38453b34abad1953b2b8722"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.2.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "c020028bb22a2f23cbd88cb92cf47cbb8c98513f"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.8.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolvePolyesterForwardDiffExt = "PolyesterForwardDiff"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveStaticArraysExt = "StaticArrays"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

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

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "cce98ad7c896e52bb0eded174f02fc2a29c38477"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.18.0"

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

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d2fdac9ff3906e27f7a618d47b676941baa6c80c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.10"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "9ae599cd7529cfce7fea36cf00a62cfc56f0f37c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.4"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "1158cfdf0da5b0eacdfcfba7c16b174a37bdf6c7"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.2.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "97e5d0b7e5ec2e68eec6626af97c59e9f6b6c3d0"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.65.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "25349bf8f63aa36acbff5e3550a86e9f5b0ef682"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.6"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "e15f5a73f0d14b9079b807a9d1dac13e4302e997"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.24.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "ba4d38faeb62de7ef47155ed321dce40a549c305"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.2+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "b479c7a16803f08779ac5b7f9844a42621baeeda"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.21"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

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

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
git-tree-sha1 = "5d54d076465da49d6746c647022f3b3674e64156"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.8"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "66c68a20907800c0b7c04ff8a6164115e8747de2"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.0"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "dd260903fdabea27d9b6021689b3cd5401a57748"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.20.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "6129a4faf6242e7c3581116fbe3270f3ab17c90d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.67"

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
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "532e22cf7be8462035d092ff21fada7527e2c488"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.6+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

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
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

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
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─ba9b6172-0234-442c-baaa-876b12f689bd
# ╠═8197b6f0-00b7-11ef-2142-4b7cbbaefd90
# ╟─d8406b01-e36f-4953-a5af-cd563005c2a1
# ╟─692aa6e2-21a9-4dad-a719-9d8ad88bd467
# ╠═e9ef08a6-e823-4f2b-9fcc-e3646171cefb
# ╟─6e9b0961-04b7-436b-90f9-9804e06d10f8
# ╟─05ef8174-e9a6-4280-8640-08d74635fba2
# ╠═31332dd6-1ef2-4181-a638-35980f470552
# ╠═d811efdc-c133-4146-a507-0786cbde68ad
# ╟─c2359767-5256-4403-b6fe-0ded5ba41efa
# ╠═e3c3abaf-0a17-4d84-945d-a1bf9bdfc0de
# ╟─49d045f8-a510-4f7b-a8c9-21ed3b4e0e6c
# ╠═0f37dc8a-fb6b-424c-84aa-0e46ef89437b
# ╟─5d2886f2-87a5-43ae-9064-448478f957ca
# ╠═82321d1e-50b9-4362-afc5-d256a506bed7
# ╟─6b4ba69d-0a70-4c66-b4bc-63b717348471
# ╠═7a0fa3cd-a270-4947-b7ba-d30110139795
# ╟─0025be9a-6320-415a-b2ac-70f24bb0d12b
# ╟─66affc7b-f0b0-48f4-b8bd-a1ff16d0357e
# ╠═f60ecd21-fda2-4925-ba86-ac5ad5b7c7d0
# ╟─41cc0b76-5e30-41dd-a648-5491e1b4bd22
# ╟─c641cdf2-aabf-45e7-bb12-83d73a93b5b7
# ╠═9b0581c3-6fa8-44d8-8ebb-d2fda2841230
# ╟─f0164280-14fd-48a5-9ae8-8b31f4580da8
# ╠═5374c9db-6b90-4775-a3f0-4c226111d7d5
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
# ╠═99922f12-34b6-4a9d-b3af-2b61708725fa
# ╟─0bd09a85-2fe9-491d-b6f9-f2d74872e195
# ╠═f5c8f93c-c09a-4c4f-b6c4-18b12ff56878
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
# ╠═be538898-0493-45dc-b2ac-5d079981ef3b
# ╠═448d1e06-c4b2-4954-978e-00c2ed274ae0
# ╟─403993d0-1613-4ac1-a527-9362993e28e6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
