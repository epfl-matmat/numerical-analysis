### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 62ff56f9-caec-4352-b547-973ce2bfcb8f
begin
	using Plots
	using PlutoUI
	using PlutoTeachingTools
	using LaTeXStrings
	using Printf
	using ForwardDiff
	using LinearAlgebra
	using HypertextLiteral
end

# ╔═╡ 0441b7ed-b62f-417c-a4f6-33051faac08e
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/04_Nonlinear_equations.pdf)
"""

# ╔═╡ b258697b-4c4f-4e47-8472-73b72178c108
TableOfContents()

# ╔═╡ 6755e1bd-d225-46e6-b80f-fd7287c9dcef
md"""
# Root finding and fixed-point problems

We saw in the introduction that problems where one wants to find a **fixed point** or a **root** (**zero**) of a function arise naturally in scientific questions.
Moreover if the underlying function is more complicated than just a simple polynomial, it becomes **quickly challenging to find such a fixed point**
using pen and paper.
In this chapter we develop numerical methods for solving such problems.

We start by defining these problems mathematically.
In the past you certainly met **root finding problems**
where one seeks to find a scalar $x$ such that a function $f : \mathbb{R} \to \mathbb{R}$ is zero, i.e. $f(x) = 0$.

The vector-valued analogue to this problem is
that one has a family of functions $f_1, f_2, \ldots, f_n$
each depending on multiple variables $x_1, x_2, \ldots, x_n$
and one seeks to simultaniously satisfy the system of equations
```math
\left\{ \begin{aligned}
f_1(x_1, \ldots, x_n) &= 0\\
f_2(x_1, \ldots, x_n) &= 0\\
\vdots\\
f_n(x_1, \ldots, x_n) &= 0
\end{aligned} \right.
```
Introducing the compact vector notation
```math
\textbf{x} = \left( \begin{array}{c} x_1 \\ \vdots \\ x_n \end{array}\right),
\qquad
\textbf{f}\left(\textbf{x}\right)
= \left( \begin{array}{c}
f_1(x_1, \ldots, x_n)\\
f_2(x_1, \ldots, x_n)\\
\vdots\\
f_n(x_1, \ldots, x_n)
\end{array}\right),
```
we can define the vector version of the root-finding problem:
"""

# ╔═╡ b03b35dd-13c3-4684-88ef-b74c91e4a8dc
md"""
!!! info "Definition: Rootfinding problem"
    Given a continuous vector-valued function
    $f : \mathbb{R}^n \to \mathbb{R}^n$
    find a vector $\textbf{x}_\ast \in \mathbb{R}^n$
    such that $f(\textbf{x}_\ast) = \textbf{0} \in \mathbb{R}^n$.
    Such a value $\mathbf{x}_\ast$ is called a **root** of ${f}$.
"""

# ╔═╡ 4e3ba905-31dd-4d29-9e45-4a77228efea0
md"""
!!! warning "Example: Intersecting circle and parabola"
    Consider the problem of the plot below, where we seek to find the
	marked intersection points of a circle and a parabola.
	These intersection points need to satisfy both the
	equation for a parabola ($x_2 = 2x_1^2$)
	and the equation for the circle with radius $1$
	(i.e. $x_1^2 + x_2^2 = 1$).
	Rearringing the equations to
	```math
	\begin{aligned}
	0 &= 1 - x_1^2 - x_2^2\\
	0 &= 2x_1^2 - x_2
	\end{aligned}
	```
	thus defines the root finding problem
	${f}(\textbf{x}) = \textbf{0}$
	with ${f} : \mathbb{R}^2 \to \mathbb{R}^2$ with
	```math
	\left\{ \begin{aligned}
	f_1(x_1, x_2) &= 1 - x_1^2 - x_2^2 \\
	f_2(x_1, x_2) &= 2x_1^2 - x_2\\
	\end{aligned} \right. .
	```

	In this example we can still find the solution analytically.
	E.g. insert the parabola equation $x_2 = 2x_1^2$ into
	the circle equation $x_1^2 + x_2^2 = 1$
	to yield $x_2^2 + \frac12 x_2 -1 = 0$
	which can be solved for $x_2$ to yield
	```math
	\begin{aligned}
	x_{2,\ast} &= \frac{ \sqrt{17} - 1 }{4} \approx 0.781 &
	x_{1,\ast} &= \sqrt{ \frac{1}{2} x_{2,\ast}  } \approx 0.624
	\end{aligned}
	```
	Not the nicest values to compute by hand ...
"""

# ╔═╡ 3e12bce6-701d-4ed1-baf3-13fe888ede07
begin
	x2_ast = (sqrt(17)-1)/4
	x1_ast = sqrt( x2_ast/2 )
end;

# ╔═╡ beba7886-fb36-45b1-b416-689733aecf4b
let
	xs = range(-1+1e-14, 1-1e-14; length=1000)
	
	p = plot(xs, @. sqrt(1 - xs^2); lw=2, c=1, aspect_ratio=1, label=L"Unit circle $x_1^2 + x_2^2 = 1$", xlabel=L"x_1", ylabel=L"x_2", xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
	plot!(p, xs, @. -sqrt(1 - xs^2); lw=2, c=1, label="")
	plot!(p, x -> 2x^2, lw=2, c=2, label=L"Parabola $x_2 = 2x_1^2$")
	scatter!(p, [-x1_ast, x1_ast], [x2_ast, x2_ast], label="Intersections")
end

# ╔═╡ 126314f4-8488-4338-a419-fc87e8965ee3
md"""
An equivalent and sometimes more intuitive way to think about solving equations ${f}(\mathbf{x}) = \mathbf{0}$ is to recast them as a fixed point-problem:

!!! info "Definition: Fixed-point problem"
    Given a continous function ${g} : \mathbb{R}^n \to \mathbb{R}^n$ a point
    $\mathbf{x}_\ast$, such that ${g}(\mathbf{x}_\ast) = \mathbf{x}_\ast$,
    is called a **fixed point** of $g$.
    The problem of finding a fixed point of ${g}$
	is equivalent to finding a root
    of the function ${f}(\mathbf{x}) = {g}(\mathbf{x}) - \mathbf{x}$.

Such transformations can also be reversed, i.e. given an $\mathbf{f}$
for root finding, we could also seek a fixed-point of $\mathbf{g}(\mathbf{x}) = \mathbf{x} + \mathbf{f}(\mathbf{x})$.
Moreover the transformations are not unique.
Indeed, for any $c\neq 0$ finding a fixed-point $\mathbf{x}_\ast$ of ${g}(\mathbf{x}) = \mathbf{x} + c \, {f}(\mathbf{x})$
implies ${f}(\mathbf{x}_\ast) = \mathbf{0}$.
"""

# ╔═╡ 65e1a532-64c3-427b-a5ec-3fcf0f7eccec
md"""
!!! warning "Example: Intersecting circle and parabola (continued)"
	We continue our example of a circle intersecting a parabola.
	Recall that the two equations to satisfy are
	$x_2 = 2x_1^2$ (parabola) and $x_1^2 + x_2^2 = 1$ (circle).
	We restrict ourselves to the upper-right quadrant, i.e. we will assume $x_1 \geq 0$ and $x_2 \geq 0$.

	- Rearranging the circle equation to $x_1 = \sqrt{1-x_2}$
	  thus leads to the fixed-point problem
	  ```math
	  \mathbf{x} = {g}_D(\mathbf{x})
	  \qquad \text{with}\qquad
	  {g}_D(\mathbf{x})
	  = \begin{pmatrix} \sqrt{1-x_2^2} \\ 2x_1^2 \end{pmatrix}.
	  ```
	- Rearranging instead the parabola equation to $x_1 = \sqrt{\frac{x_2}2}$
	  and the circle equation to $x_2 = \sqrt{1-x_1}$ yields the problem
	  ```math
	  \mathbf{x} = {g}_C(\mathbf{x})
	  \qquad \text{with}\qquad
	  {g}_C(\mathbf{x})
	  = \begin{pmatrix} \sqrt{\frac{x_2}2} \\ \sqrt{1-x_1^2} \end{pmatrix}.
	  ```
"""

# ╔═╡ c454473a-a5ef-4b08-af45-da4a7acbb1d0
md"""
!!! exercise "Exercise"
	Verify that the fixed-point problems $\mathbf{x} = {g}_D(\mathbf{x})$
	and $\mathbf{x} = {g}_C(\mathbf{x})$
	are equivalent to the root-finding problem $\mathbf{f}(\mathbf{x}) = \mathbf{0}$.
"""

# ╔═╡ d3ba1da7-41eb-4afd-a949-19123b7a6724
md"""
## Fixed-point iterations

In the [introduction](https://teaching.matmat.org/numerical-analysis/01_Introduction.html) when discussing the diode model,
we already mentioned fixed-point iterations as one way of solving the fixed-point problem.

!!! info "Algorihm 1: Fixed-point iteration"
    Given a fixed-point map ${g}$ and initial guess $\mathbf{x}^{(0)}$, iterate
    ```math
    \mathbf{x}^{(k+1)} = {g}(\mathbf{x}^{(k)})  \qquad{\text{for $k = 1, 2, \ldots$}}.
    ```
    until a **stopping criterion** is reached, see text below.
"""

# ╔═╡ e6c13ceb-fb13-4159-b0bc-a854406d9f21
md"""
If $\lim_{k\to\infty} \mathbf{x}^{(k)} = \mathbf{x}_\ast$ for the sequence $\mathbf{x}^{(k)}$ generated by Algorithm 1,
then ${g}(\mathbf{x}_\ast) = \mathbf{x}_\ast$,
i.e. $\mathbf{x}_\ast$ is a fixed point of ${g}$.
A typical stopping criterion is to check if the vector norm
between subsequent iterates, i.e. $\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\|$
drops below a certain tolerance.
Here $\|x\|$ denotes the Euclidean norm
$\|x\| = \sqrt{x_1^2 + x_2^2 + \cdots + x_n^2}$.
We will discuss in the [Convergence analysis section](#Stopping-criteria-and-residual) why this is a good choice.
"""

# ╔═╡ 8b9bae6b-02ae-455e-a6f0-a10ca6db694e
md"""
If we are faced with a fixed point problem (finding $\mathbf{x}_\ast$ such that ${g}(\mathbf{x}_\ast) = \mathbf{x}_\ast$) then Algorithm 1 can be directly applied. However, to **apply** the **fixed-point method** to a **root-finding problem** (seek $\mathbf{x}_\ast$ s.t. ${f}(\mathbf{x}_\ast) = \mathbf{0}$)
we first need to **rewrite the non-linear equation** ${f}(\mathbf{x}) = \mathbf{0}$ into a fixed-point problem, i.e. identify a suitable ${g}$, such that if $\mathbf{x}_\ast$ is a root of $\mathbf{f}$, then
```math
{f}(\mathbf{x}_\ast) = 0 \qquad \Longleftrightarrow \qquad \mathbf{x}_\ast = {g}(\mathbf{x}_\ast).
```
On $\mathbf{g}$ we then apply fixed-point iteration.

We saw one example how to achieve this rewriting in the discussion
of the "Intersecting circle and parabola" example in the [section above](#Root-finding-and-fixed-point-problems), where we first defined a root-finding problem and then two equivalent fixed-point problems for the same task.
"""

# ╔═╡ 9f80abeb-432c-4b68-a570-9e32e2510fe0
md"""
We already noted in the [introduction](https://teaching.matmat.org/numerical-analysis/01_Introduction.html) that these **fixed-point iterations do not always converge**.

Let us investigate this closer for the circle intersection problem of the [section above](http://localhost:1234/edit?id=5189964e-fa08-11f0-88ac-47e68bd3c143#Root-finding-and-fixed-point-problems). We noted that the problem can be formulated as two fixed-point problems 
$\mathbf{x} = {g}_D(\mathbf{x})$
and $\mathbf{x} = {g}_C(\mathbf{x})$
involving the functions
```math
\begin{aligned}
{g}_D(\mathbf{x})
	  &= \begin{pmatrix} \sqrt{1-x_2^2} \\ 2x_1^2 \end{pmatrix} &
 {g}_C(\mathbf{x})
	&= \begin{pmatrix} \sqrt{\frac{x_2}2} \\ \sqrt{1-x_1^2} \end{pmatrix}
\end{aligned}
```
which we code up as follows:
"""

# ╔═╡ cabca1d2-e287-4d61-9d03-16108e10ac9c
gD(x) = [sqrt(1 - x[2]^2);
	              2x[1]^2]

# ╔═╡ c8cb64e6-a12a-4308-89bd-73bf178c4986
gC(x) = [    sqrt(x[2]/2);
		 sqrt(1 - x[1]^2)]

# ╔═╡ b4a0754f-dbb3-4de1-ba94-35832ece1048
md"""They take as input and provide as output a two-dimensional vector."""

# ╔═╡ aec943b1-5865-4306-a780-fd730afd208a
md"""
We visualise the iterations for `gD` and `gC` by visualising the sequence of points $(x_1, x_2)$ they generate. In the graphics below both the selection of points, the starting point as well as the number of iterations can be selected.
"""

# ╔═╡ 1b35b5da-7e1d-42d9-a189-c20d6957e939
md"""
- Fixed-point function `gfix = ` $(@bind gfix Select([gD, gC])) 
- Starting point `xinit = ` $(@bind xinit Select([[0.7, 0.8], [0.2, 0.8], [0.0, 0.4]]))
- Number of iterations `n_iter = ` $(@bind n_iter Slider(1:10; show_value=true, default=2))
"""

# ╔═╡ 08578e28-28c0-4cb1-9934-734034dec251
let
	p = plot(x ->  sqrt(1 - x^2); lw=2, c=1, aspect_ratio=1, label=L"Unit circle $x_1^2 + x_2^2 = 1$", xlabel=L"x_1", ylabel=L"x_2", xlims=(-0.2, 1.2), ylims=(-0.2, 1.2), legend=:bottomright)
	plot!(p, x -> -sqrt(1 - x^2); lw=2, c=1, label="")
	plot!(p, x -> 2x^2, lw=2, c=2, label=L"Parabola $x_2 = 2x_1^2$")
	scatter!(p, [x1_ast], [x2_ast], label="Fixed point", c=:black, ms=3)

	
	scatter!(p, [xinit[1]], [xinit[2]]; color=:grey, mark=:o, ms=5, label="Starting point")
	x = xinit
	for k = 1:n_iter
		xnext = gfix(x)

		plot!(p, [x[1], xnext[1]], [x[2], xnext[2]];
			  arrow=true, color=3, label="", lw=2)

		if xnext[1] > 1 || xnext[2] > 1
			break
		end
		
		x = xnext
	end

	p
end

# ╔═╡ 2f670892-be3e-43f3-b49f-803d8ddf756e
md"""
We notice:
- While $g_C$ converges (i.e. the arrows get closer and closer to the fixed-point shown in black), $g_D$ diverges.
- While the sequence of points in both cases "circles around" the fixed point. For the converging function $g_C$ the circle spirals inwards, while for $g_D$ the circle spirals outwards.
- Intuitively the fixed-point method therefore converges if it makes the distance between the iteration points and the fixed-point smaller.

We will make this observation more precise in the section on [Convergence analysis](#Convergence-analysis).
"""

# ╔═╡ b385d67b-c7e0-4a5e-87e1-b11f00ba638a
md"""
### Optional: Rewriting root-finding and fixed-point problems
"""

# ╔═╡ c6d53a9e-0b6d-4ad3-a7c8-c34070f75ff6
md"""
We consider solving non-linear equation $f(x) = 0$ with
```math
\tag{1}
f(x) = x + \log(x + 1) - 2
```
"""

# ╔═╡ 1319313d-e7da-4ff5-9725-548df0e7f699
md"""
In this one can construct four equivalent fixed-point equations $x = g_i(x)$ with $i=1,2,3,4$, namely
"""

# ╔═╡ 858510eb-7768-4521-ba20-af5b1aa33e08
begin
	g₁(x) = x - (1/2) * (log(x + 1) + x - 2)
	g₂(x) = 2 - log(x + 1)
	g₃(x) = exp(2-x) - 1
	g₄(x) = (1/2) * x * (log(x+1) + x)
end;

# ╔═╡ 7bdb39f4-4fab-43a1-a9d9-12e198cef107
md"""
We show how to construct $g_3$ and $g_4$ explicitly.

First $g_3$. At convergence we have that
```math
\begin{aligned}
&& 0 &= f(x) = x + \log(x + 1) - 2 && \text{(add $2-x$)} \\
\Leftrightarrow && 2-x &= \log(x+1) && \text{(take $\exp$)}\\
\Leftrightarrow && \exp(2-x) &= x+1 && \text{(subtract $1$)}\\
\Leftrightarrow && \exp(2-x) -1 &= x
\end{aligned}
```
This is now a fixed-point problem in $x$ where we seek a fixed point
of the function $g_3(x) = \exp(2-x) -1$,
which is just the left-hand side of the expression.

Now $g_4$. Again starting from
 ```math
\begin{aligned}
&& 0 &= f(x) = x + \log(x + 1) - 2 && \text{(add $2$)} \\
\Leftrightarrow && 2 &= x + \log(x+1) && \text{(divide by $2$)}\\
\Leftrightarrow && 1 &= \frac12 \left(x + \log(x+1) \right) && \text{(multiply by $x$)}\\
\Leftrightarrow && x &= \underbrace{\frac x 2 \left(x + \log(x+1)\right)}_{g_4(x)} 
\end{aligned}
```
where again we get a problem $x = g_4(x)$, i.e. finding a fixed point of $g_4$.
"""

# ╔═╡ 5e21e4c5-5538-40cf-a955-837809f7c3c3
md"""
### Convergence analysis

#### Scalar functions

To understand the aforementioned behaviour mathematically,
we first stay within the setting of a **scalar fixed-point problem**.
That is we consider the setting finding a fixed-point $x_\ast = g(x_\ast)$ for a differentiable function $g : \mathbb{R} \to \mathbb{R}$.

To understand convergence more rigorously we study the beviour of the error sequence
```math
e^{(k)} = x^{(k)} - x_\ast
```
The goal is to find an expression that relates $e^{(k+1)}$ with $e^{(k)}$ and in this way *recursively* to $e^{(0)}$, the error of the initial guess.
The hope is this enables to find a condition for $g$ that tells us when the fixed point iterations converge (i.e. the error goes to $0$ as $k\to\infty$) or not
(the error stays non-zeroas $k\to\infty$).
"""

# ╔═╡ 05989b7d-e694-4860-9a8e-6e6cb52aea8b
md"""
First we note $x^{(k)} = x_\ast + e^{(k)}$ and we develop the 
Taylor expansion of $g$ around $x_\ast$:
```math
\tag{3}
\begin{align*}
g(x^{(k)}) &= \underbrace{g(x_\ast)}_{=x_\ast} + g'(x_\ast) (x^{(k)}-x_\ast) + R(x^{(k)}) \\
\end{align*}
```
where 
```math
R(x) = O(|x - x_\ast|^2) \qquad \text{as $x \to x_\ast$},
```
which means that there exist positive constants $M, δ > 0$ such that
```math
|R(x)| ≤ M |x - x_\ast|^2 \qquad \forall 0 < |x-x_\ast| < δ,
```
i.e. that $R(x)$ stays bounded by the next term of the Taylor expansion
(up to a multiplicative constant).
The notation $O(|x - x_\ast|^2)$ is thus a
mathematically precise way of saying
that there are more terms that we don't show
but their order is at most $|x - x_\ast|^2$.
See also the discussion in [Revision and preliminaries](https://teaching.matmat.org/numerical-analysis/03_Preliminaries.html) on Taylor approximations.
"""

# ╔═╡ 01db98ec-daf2-4779-9f31-c3271039f44c
md"""
Using (3) and the key fixed-point iterations equation, $x^{(k+1)} = g(x^{(k)})$
we obtain
```math
\begin{aligned}
e^{(k+1)}
&= x^{(k+1)} - x_\ast = g(x^{(k)}) - x_\ast \\
&\stackrel{(3)}{=} x_\ast + g'(x_\ast) \underbrace{(x^{(k)}-x_\ast)}_{=e^{(k)}} + O(|e^{(k)}|^2) - x_\ast \\
&=  g'(x_\ast)\, e^{(k)} + O(|e^{(k)}|^2).
\end{aligned}
```
Taking moduli on both sides:
```math
|e^{(k+1)}| = |g'(x_\ast)| \ |e^{(k)}| + O(|e^{(k)}|^2) 
```

We employ this relation now in a recursive argument.
Assume we choose a good initial guess,
then $x^{(0)}$ is close enough to $x_\ast$, such that $O(|e^{(0)}|^2)$
is neglibile compared to $|g'(x_\ast)| \ |e^{(0)}|$.
Similarly, provided that the iteration makes progress,
 $O(|e^{(1)}|^2)$ is in turn
smaller than $|g'(x_\ast)| \ |e^{(1)}|$ and so forth.
Therefore
```math
\begin{aligned}
|e^{(k+1)}| &= |g'(x_\ast)| \ |e^{(k)}| + O(\text{small}) \\
&= |g'(x_\ast)|^2 \ |e^{(k-1)}| + O(\text{small}) \\
&= \ldots \\
&= |g'(x_\ast)|^{k+1} \ |e^{(0)}| + O(\text{small})
\end{aligned}
```
In other words as $k \to \infty$, i.e. the iteration progresses,
$|e^{(k+1)}|$ approaches zero
if $|g'(x_\ast)| < 1$.
"""

# ╔═╡ 9176b666-41f7-436e-b5ad-61b196a8b35b
md"""
!!! note "Theorem 1 (scalar version)"
    Let $g : \mathbb{R} \to \mathbb{R}$
	be a once differentiable function
    and $x_\ast \in \mathbb{R}$ be a fixed point of $g$.
    If $|g'(x_\ast)| < 1$, then
    there exists an $ε > 0$,
    such that for all $x^{(0)} \in [x_\ast - ε, x_\ast + ε]$
    the fixed-point iterations
	$x^{(k+1)} = g(x^{(k)})$ converge to $x_\ast$,
    i.e. $\lim_{k\to\infty} x^{(k)} = x_\ast$

We will generalise this theorem to the vector case in the following secition.
"""

# ╔═╡ 511fcf1a-5ea8-4827-9d24-6a43e7dcccd6
md"""
#### Higher dimensions

We now consider the generalisation of the above argument to the vector setting,
i.e. finding a fixed-point $\mathbf{x}_\ast = {g}(\mathbf{x}_\ast) \in \mathbb{R}^n$ of a function ${g} : \mathbb{R}\to\mathbb{R}$.
To make a similar argument to the scalar case, we need to consider again the Talyor expansion of ${g}(\mathbf{x}^{(k)}) = {g}(\mathbf{x}_\ast + \mathbf{e}^{(k)})$ around $\mathbf{x}_\ast$, where as before $\mathbf{e}^{(k)} = \mathbf{x}^{(k)} - \mathbf{x}_\ast$.

First recall that for a vector-valued function the Taylor expansion to first order reads
```math
\tag{4}
{g}(\textbf{y})
= {g}(\textbf{x}) + \textbf{J}_{g}(\textbf{x}) \, (\textbf{y} - \textbf{x})
+ O\Big(( \textbf{y} - \textbf{x} )^2\Big).
```
In this the **Jacobian matrix**
$\textbf{J}_{g}(\textbf{x}) \in \mathbb{R}^{n\times n}$
is the collection of all partial derivatives of ${g}$ *evaluated at $\mathbf{x}$*, i.e.
```math
\textbf{J}_{g} = \left(\begin{array}{cccc}
\frac{\partial g_1}{\partial x_1} &
\frac{\partial g_1}{\partial x_2} &
\ldots &
\frac{\partial g_1}{\partial x_n} \\

\frac{\partial g_2}{\partial x_1} &
\frac{\partial g_2}{\partial x_2} &
\ldots &
\frac{\partial g_2}{\partial x_n} \\

\vdots & & \ddots & \vdots\\

\frac{\partial g_n}{\partial x_1} &
\frac{\partial g_n}{\partial x_2} &
\ldots &
\frac{\partial g_n}{\partial x_n}
\end{array}\right).
```
See also the discussion on multi-dimensional Talyor approximations in [Revision and preliminaries](https://teaching.matmat.org/numerical-analysis/03_Preliminaries.html).

The Jacobian very much plays the role of a generalised derivative
of a multidimensional function ${g}$. This is why we will sometimes also
borrow the notation $g'$ from the scalar case to denote the Jacobian.

Also not that (just like any derivative) it is a function of an independent
variable (here $\textbf{x}_\ast$).
"""

# ╔═╡ c7068f2e-7313-44bc-85a4-785f4d4adc60
md"""
!!! warning "Intersecting circle and parabola (continued)"
	As an example we compute the Jacobian of the two fixed-point maps
	that arose in the study of the curve intersection problem between circle
	and parabola, which we considered above, namely
	```math
	\begin{aligned}
	{g}_D(\mathbf{x})
		  &= \begin{pmatrix} \sqrt{1-x_2^2} \\ 2x_1^2 \end{pmatrix} &
	 {g}_C(\mathbf{x})
		&= \begin{pmatrix} \sqrt{\frac{x_2}2} \\ \sqrt{1-x_1^2} \end{pmatrix}.
	\end{aligned}
	```
	For the first we find
	```math
	\mathbf{J}_{{g}_D}(\mathbf{x}) = \begin{pmatrix}
	\frac{\partial}{\partial x_1} \sqrt{1-x_2^2} &
	\frac{\partial}{\partial x_2} \sqrt{1-x_2^2} \\
	\frac{\partial}{\partial x_1} 2x_1^2 &
	\frac{\partial}{\partial x_2} 2x_1^2
	\end{pmatrix}
	= \begin{pmatrix}
	0 &
	\frac{-x_2}{\sqrt{1-x_2^2}} \\
	4 x_1 &
	0
	\end{pmatrix}
	```
	and for the second
	```math
	\mathbf{J}_{{g}_C}(\mathbf{x}) = \begin{pmatrix}
	\frac{\partial}{\partial x_1} \sqrt{\frac{x_2}2} &
	\frac{\partial}{\partial x_2} \sqrt{\frac{x_2}2} \\
	\frac{\partial}{\partial x_1} \sqrt{1-x_1^2} &
	\frac{\partial}{\partial x_2} \sqrt{1-x_1^2}
	\end{pmatrix}
	= \begin{pmatrix}
	0 &
	\frac{1}{4\sqrt{\frac{x_2}2}} \\
	\frac{-x_1}{\sqrt{1-x_1^2}} &
	0
	\end{pmatrix}
	```
"""

# ╔═╡ 8a8ac29f-f664-4c62-bda7-eeb5214b2bc8
md"""
Now, using the Taylor expansion of (4) with $\mathbf{y} = \mathbf{x}_\ast + \mathbf{e}^{(k)}$
and $\mathbf{x} = \mathbf{x}_\ast$
we obtain analogously to the scalar case for the **error in the fixed-point iterations**
```math
\tag{5}
\begin{aligned}
\mathbf{e}^{(k+1)}
&= \mathbf{g}(\mathbf{x}_\ast + \mathbf{e}^{(k)}) - \mathbf{x}_\ast \\
&\stackrel{(4)}{=} {g}(\mathbf{x}_\ast) + \mathbf{J}_{g}(\mathbf{x}_\ast) \, \mathbf{e}^{(k)} + O\Big((\mathbf{e}^{(k)})^2\Big) - \mathbf{x}_\ast \\
&=  \mathbf{J}_{g}(\mathbf{x}_\ast)\, \mathbf{e}^{(k)} + O\Big((\mathbf{e}^{(k)})^2\Big).
\end{aligned}
```
i.e. we see the Jacobian indeed play the role of the derivative of $\mathbf{g}$.
"""

# ╔═╡ d4e3f192-583e-4d3f-af7d-0767d2c9c0a1
md"""
The above expression again relates the error in step $k+1$ to the error in step $k$.  Following the analogy from the scalar case we now want to understand how the size of this error varies between iterations. We need one further ingredient, namely the generalisation of the modulus of the gradient $|g'(x_\ast)|$ to higher dimensions, where notably $\mathbf{J}_{g}$ is now an $n\times n$ matrix.

We need the

!!! info "Definition: Matrix norm"
    Given $\mathbf{M} \in \mathbb{R}^{m \times n}$ a real matrix (not necessarily square), we define the matrix norm of $\mathbf{M}$ as
	```math
		\|\mathbf{M}\| = \max_{\stackrel{\mathbf{v} \in \mathbb{R}^n}{\mathbf{v}\neq\mathbf{0}}} \frac{\|\mathbf{M}\mathbf{v}\|}{\|\mathbf{v}\|}.
	```

The definition of the matrix norm implies in particular that

!!! info ""
	```math
	\|\mathbf{M}\mathbf{v}\| \leq \|\mathbf{M}\| \, \|\mathbf{v}\|
	\quad \forall \mathbf{v} \in \mathbb{R}^n.
	```
"""

# ╔═╡ 9c9719e3-ec6c-4bdc-b05b-ab4bd4119cb9
md"""
We now take vector norms on either side of (5) and make use of this last inequality to obtain

```math
\|\mathbf{e}^{(k+1)}\| \leq \left\| \mathbf{J}_{g}(\mathbf{x}_\ast) \right\| \, \left\| \mathbf{e}^{(k)} \right\| + O(\|\mathbf{e}^{(k)}\|^2 ).
```

Under the assumption that our initial guess $\mathbf{x}^{(0)}$ is sufficiently close
to $\mathbf{x}_\ast$ we can again follow a recursive argument to obtain
```math
\begin{aligned}
\|\mathbf{e}^{(k+1)}\| &=  \left\| \mathbf{J}_{g}(\mathbf{x}_\ast) \right\|^{k+1} \ \|\mathbf{e}^{(0)}\| + O(\text{small})
\end{aligned}
```
where $O(\text{small})$ is a small term that we do not make more precise for simplicity.

We are again faced with the conclusion that  as $k \to \infty$, i.e. the iteration progresses, that the error norm $\|\mathbf{e}^{(k+1)}\|$ approaches zero
if $\left\| \mathbf{J}_{g}(\mathbf{x}_\ast) \right\| <  1$.

The following theorem summarises our argument
"""

# ╔═╡ fe6a0ff6-f70d-404a-8fd5-8a9185da2ee7
md"""
!!! note "Theorem 1"
    Let ${g} : \mathbb{R}^n \to \mathbb{R}^n$
	be a once differentiable function
    and $\mathbf{x}_\ast \in \mathbb{R}^n$ be a fixed point of $g$.
    If $\left\| \mathbf{J}_{g}(\mathbf{x}_\ast) \right\| <  1$, then
    there exists an $ε > 0$, such that for all $x^{(0)}$
	with $\|\mathbf{x}^{(0)} - \mathbf{x}_\ast\| < ε$
    - the fixed-point iterations
      $\mathbf{x}^{(k+1)} = {g}(\mathbf{x}^{(k)})$ converge
	  to $\mathbf{x}_\ast$,
      i.e. $\lim_{k\to\infty} \mathbf{x}^{(k)} = \mathbf{x}_\ast$
    - Moreover the **convergence rate**
      *([formal definition below](#Convergence-order))*
	  is given by
      ```math
      \lim_{k\to\infty} \frac{\|\mathbf{x}^{(k+1)} - \mathbf{x}_\ast\|}{\|\mathbf{x}^{(k)} - \mathbf{x}_\ast\|} = \| \mathbf{J}_{g}(\mathbf{x}_\ast) \|,
      ```
       i.e. the smaller the norm of the Jacobian, the faster the convergence.
"""

# ╔═╡ 33c6b6ed-267e-4062-8390-fa02b8624545
md"""
We notice the central role the Jacobian matrix norm plays for both determining if the iterations converge and at which rate. It will therefore be important for us to compute such matrix norms. We will provide more details on this point in the final sections of [Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/05_Direct_methods.html). For now we
only note the key result that

!!! info ""
	For any matrix  $\mathbf{M} \in \mathbb{R}^{m \times n}$
	```math
	\|\mathbf{M}\| = \sqrt{ \lambda_\text{max}(\mathbf{M}^T\mathbf{M}) }.
	```
	where $\lambda_\text{max}(\mathbf{M}^T \mathbf{M})$ is the *largest* eigenvalue
	of the matrix $\mathbf{M}^T \mathbf{M}$.
"""

# ╔═╡ 29931ae9-bcb7-4ec0-b397-a89c491d950e
md"""
!!! warning "Intersecting circle and parabola (continued)"
	Based on Theorem 1 we now want to understand why searching the intersection
	of circle and parabola using fixed-point iterations on the map ${g}_C$
	is successful while the iterations using ${g}_D$ fails.
	We already computed the respective Jacobians as
	```math
	\begin{aligned}
	\mathbf{J}_{{g}_C}(\mathbf{x})
	&= \begin{pmatrix}
	0 &
	\frac{1}{4\sqrt{\frac{x_2}2}} \\
	\frac{-x_1}{\sqrt{1-x_1^2}} &
	0
	\end{pmatrix} &
	\mathbf{J}_{{g}_D}(\mathbf{x}) &= \begin{pmatrix}
	0 &
	\frac{-x_2}{\sqrt{1-x_2^2}} \\
	4 x_1 &
	0
	\end{pmatrix}
	\end{aligned}
	```

	According to Theorem 1 to understand the convergence properties
	we need to investigate 
	$\left\| \mathbf{J}_{g_C}(\mathbf{x}_\ast) \right\|$
	and
	$\left\| \mathbf{J}_{g_D}(\mathbf{x}_\ast) \right\|$
	and for this compute respectively the largest eigenvalues of
	```math
	\begin{align*}
	\mathbf{J}_{g_C}(\mathbf{x}_\ast)^T \mathbf{J}_{g_C}(\mathbf{x}_\ast)
	&= \begin{pmatrix}
	\frac{1}{8x_{2,\ast}} & 0 \\
	0 & \frac{x_{1,\ast}^2}{1-x_{1,\ast}^2}
	\end{pmatrix}
	& \mathbf{J}_{g_D}(\mathbf{x}_\ast)^T  \mathbf{J}_{g_D}(\mathbf{x}_\ast)
	&=\begin{pmatrix}
	\frac{x_{2,\ast}^2}{1-x_{2,\ast}^2} & 0\\
	0 & 16 x_{1,\ast}^2
	\end{pmatrix}
	\end{align*}
	```
	Since these are diagonal matrices their eigenvalues can be directly read off.
	In particular we compute numerically:
"""

# ╔═╡ 5409350b-4845-46e5-ac3a-f7342fc28d3d
JgC = [0                           1 / (4sqrt(x2_ast/2));
	   -x1_ast / sqrt(1-x1_ast^2)                      0]

# ╔═╡ 9d7ca2e1-41b0-4651-a734-31bf883cee37
JgC*JgC'

# ╔═╡ fffd0bd2-7b66-4a21-8311-2965e974d88b
md"""
!!! warning ""
	Therefore $\left\| \mathbf{J}_{g_C}(\mathbf{x}_\ast) \right\| \approx \sqrt{0.640} \approx 0.400 < 1$
"""

# ╔═╡ 50639d02-55d5-4fcb-8335-13fd7f6b7624
JgD = [0         -x2_ast / sqrt(1-x2_ast^2);
	   4x1_ast                            0]

# ╔═╡ da60beec-74e7-4b3f-aa09-27b806054896
JgD*JgD'

# ╔═╡ ec8b8273-3724-4ea9-91d3-259390abc55d
md"""
!!! warning ""
	and $\left\| \mathbf{J}_{g_D}(\mathbf{x}_\ast) \right\| \approx \sqrt{6.25} \approx 2.50 > 1$.

	Provided we start sufficiently close to $\mathbf{x}_\ast$
	fixed-point iterations for $g_C$ therefore converge
	while fixed-point iterations for $g_D$ diverge.
"""

# ╔═╡ 5d7b3d35-3456-48df-ad22-0ffccaa2f529
md"""
### Stopping criteria and residual

Let us come back to the question at which point to stop the fixed-point iteration algorithm. Let $\epsilon$ denote the tolerance to which we want to determine $\mathbf{x}_\ast$, i.e. we would like to stop the iterations as soon as the error is smaller, $\|\mathbf{x}^{(k)} - \mathbf{x}_\ast\| < \epsilon$.
Since $\mathbf{x}_\ast$ is not known, this expression cannot be computed during the iterations. We thus seek an alternative approach.

In a given step $\mathbf{x}^{(k)}$ we have in general not yet achieved our goal,
i.e. ${g}(\mathbf{x}^{(k)}) \neq \mathbf{x}^{(k)}$.
An idea is thus to consider exactly the descrepancy 
```math
\mathbf{r}^{(k)} = {g}(\mathbf{x}^{(k)}) - \mathbf{x}^{(k)},
```
the so-called **residual**. A natural stopping criterion is thus
"""

# ╔═╡ 284f3fa4-ce24-4b99-8bb7-6f74a4589550
md"""
!!! note "Algorithm: Fixed-point iteration stopping criterion"
    ```math
    \|\mathbf{r}^{(k)}\| = \|{g}(\mathbf{x}^{(k)}) - \mathbf{x}^{(k)}\| < \epsilon.
    ```

Employing this stopping criteria, the algorithm to find a fixed point of the function $\mathbf{g}$ becomes
"""

# ╔═╡ 899698d1-9dee-4f82-9171-f1b49aefcabe
function fixed_point_iterations_simple(g, xstart; tol=1e-6)
	# g:      Fixed-point function
	# xstart: Initial guess
	# tol:    Tolerance

	rᵏ = Inf
	xᵏ = xstart
	k  = 0
	while norm(rᵏ) ≥ tol
		xᵏ⁺¹ = g(xᵏ)
		rᵏ   = xᵏ⁺¹ - xᵏ

		k  = k + 1  # Update k
		xᵏ = xᵏ⁺¹   # Update xᵏ accordingly
	end

	# Return results as a named tuple
	(; fixed_point=xᵏ, residual=rᵏ, n_iter=k)
end

# ╔═╡ 072cdb08-b279-4583-b98b-4a28433b5b8d
md"""We apply this function:"""

# ╔═╡ 3c0d4752-de1a-4819-a63f-aa49c94ccf82
res_simple = fixed_point_iterations_simple(gC, [0.4, 0.3]; tol=1e-14)

# ╔═╡ b5ec49d7-1876-404e-b6ed-03217612b5e4
md"Since `res_simple` is now a named tuple, we can access its individual fields to e.g. get the fixed point"

# ╔═╡ 58cf700e-8bb4-4370-ae55-875d551a3308
res_simple.fixed_point

# ╔═╡ 0eb10646-4830-4d97-b574-54438a20a3d9
md"or the number of iterations required"

# ╔═╡ b422909a-2f84-4815-b2f3-9902482806e6
res_simple.n_iter

# ╔═╡ c1b1c4ff-e2c9-49fb-8794-b9ee11c5d543
md"""In practice it is often useful to **also include a cap on the maximal number of iterations** (for cases where the algorithm does not converge) and to **record a history** of the visited points, which is useful for later analysis."""

# ╔═╡ 7cb34ba9-bcfc-4596-8eb3-49f4b99eab58
function fixed_point_iterations(g, xstart; tol=1e-6, maxiter=100)
	# g:      Fixed-point function
	# xstart: Initial guess
	# tol:    Tolerance

	history_x = [xstart]
	history_r = empty(history_x)

	rᵏ = Inf     # For initial pass in while loop
	xᵏ = xstart  # Starting point of iterations
	k  = 0
	while k < maxiter && abs(rᵏ) ≥ tol
		xᵏ⁺¹ = g(xᵏ)
		rᵏ   = xᵏ⁺¹ - xᵏ
		push!(history_r, rᵏ)

		k  = k + 1  # Update k
		xᵏ = xᵏ⁺¹   # Update xᵏ accordingly
		push!(history_x, xᵏ)  # Push next point to the history
	end

	# Return results as a named tuple
	(; fixed_point=xᵏ, residual=rᵏ, n_iter=k, history_x, history_r)
end

# ╔═╡ 54f314c1-d1cf-41f1-96e5-5aca90d82b95
fixed_point_iterations_simple(gC, [0.4, 0.3]; tol=1e-14)

# ╔═╡ 73f0c043-24a6-401f-8695-7ac2ca154a7c
md"""
#### Additional remarks on the residual
- The residual is in general only an **error indicator**.
  This means that **there is no guarantee**
  that $\|g(\mathbf{x}^{(k)}) - \mathbf{x}^{(k)}\| < \epsilon$
  **always implies** $\|\mathbf{x}^{(k)} - \mathbf{x}_\ast\| < \epsilon$.
- To explain this, let us consider the **scalar case** $g : \mathbb{R} \to \mathbb{R}$.
  In this setting we can derive the **residual-error relationship** *(see derivation below)*
  ```math
  \tag{5}
  |x^{(k)} - x_\ast| = \frac{1}{|1 - g'(\xi^{(k)})|} |r^{(k)}|.
  ```
  for some $\xi^{(k)} \in [x_\ast, x^{(k)}]$.
  Note, that this is just a **conceptional expression** as determining
  $\xi^{(k)}$ is in general *as hard* as finding $x_\ast$.
  But it will be useful in some theoretical arguments.
- For converging fiterations $x^{(k)} \to x_\ast$ as $k \to \infty$.
  Therefore the interval $[x_\ast, x^{(k)}]$ gets smaller and smaller,
  such that necessarily $\xi^{(k)} \to x_\ast$
  and $g'(\xi^{(k)}) \to g'(x_\ast)$ as $k \to \infty$.
  We note it is the **gradient at the fixed point**,
  $g'(x_\ast)$, which determines **how reliable our error indicator** is.
- In particular if $|g'(x_\ast)|$ **is close to $1$**, then the denominator
  of the prefactor may blow up and the **residual criterion** $|r^{(k)}| < \epsilon$
  may well **stop the iterations too early**.
  That is the **actual error** $|x^{(k)} - x_\ast|$ may still be
  **way larger than the residual** $|r^{(k)}|$ and thus way larger than our
  desired accuracy $\epsilon$.
- In contracst if $|g'(x_\ast)| = 0$ than $|r^{(k)}| < \epsilon$
  is an excellent stopping criterion as $|x^{(k)} - x_\ast| = |r^{(k)}|$ as $k\to\infty$.
- For the **higher-dimensional case** (see detailed discussion below) the conclusion is similarly that **the residual is a reliable error indicator** if **no eigenvalue of the Jacobian** $\mathbf{J}_g$ **is close to 1**. See below for a more detailed discussion.
"""

# ╔═╡ 62ec67d7-10d4-4624-8dd7-d4dd0507d6c2
Foldable("Derivation of the residual-error relationship for the scalar case",
md"""
For simplicity we perform the derivation of the residual-error relationship
for a scalar fixed-point problem $g : \mathbb{R} \to \mathbb{R}$.
		 
Employing the [Lagrange form of the reminder](https://teaching.matmat.org/numerical-analysis/03_Preliminaries.html) of a 0-th order Taylor expansion of $g$ around $x_\ast$ we obtain
```math
g(x^{(k)}) = g({x}_\ast) + g'(\xi^{(k)}) \ (x^{(k)} - x_\ast)
```
for some $\xi^{(k)} \in [x_\ast, x^{(k)}]$.
This expression is also known as the [mean value theorem](https://en.wikipedia.org/wiki/Mean_value_theorem) (MVT).

Using the MVT we can show
```math
\begin{aligned}
x^{(k)} - x_\ast
&= x^{(k)} \underbrace{- g(x^{(k)}) + g(x^{(k)})}_{=0} - x_\ast \\
&= x^{(k)} - g(x^{(k)}) + g(x^{(k)}) - g(x_\ast) \\
&\stackrel{(MVT)}{=} x^{(k)} - g(x^{(k)}) + g'(\xi^{(k)}) (x^{(k)} - x_\ast)
\end{aligned}
```
where $\xi^{(k)} \in [x_\ast, x^{(k)}]$.
Rearranging and adding moduli we arrive at the cited **residual-error relationship**
```math
|x^{(k)} - x_\ast| = \frac{1}{|1 - g'(\xi^{(k)})|} |r^{(k)}|.
```
In particular if $g'(\xi^{(k)}) \simeq 1$ the **true error
can thus be considerably larger than the residual** !
"""
)

# ╔═╡ de918a0d-05ee-40b3-8cbd-f45274c37ed3
Foldable("In-depth discussiong of the higher-dimensional case",
md"""
We now consider a function $g : \mathbb{R}^n \to \mathbb{R}^n$. We first derive the residual-error relationship for the higher-dimensional case.
		 
For this consider the expansion of the function $t \mapsto g\big(\mathbf{x}_\ast + t (\mathbf{x}^{(k)} - \mathbf{x}_\ast)\big)$ around $t = 0$ to 0-th order. By the chain rule and again invoking the Lagrange form of the remainder we obtain
```math
g\big(\mathbf{x}_\ast + t (\mathbf{x}^{(k)} - \mathbf{x}_\ast)\big)
= g(\mathbf{x}_\ast) + \underbrace{\mathbf{J}_g\Big(\mathbf{x}_\ast + \tau (\mathbf{x}^{(k)} - \mathbf{x}_\ast)\Big)}_{=\widetilde{\mathbf{J}}^{(k)}} \ (\mathbf{x} - \mathbf{x}_\ast) \ t
```
where $\tau \in [0, 1]$ and we defined $\widetilde{\mathbf{J}}^{(k)}$. Evaluating this at $t=1$ therefore:
```math
g(\mathbf{x}^{(k)})
= g(\mathbf{x}_\ast) + \widetilde{\mathbf{J}}^{(k)} \ (\mathbf{x}^{(k)} - \mathbf{x}_\ast)
```
Following the same line of argument as in the scalar case we show
```math
\mathbf{x}^{(k)} - \mathbf{x}_\ast = \mathbf{x}^{(k)} - g(\mathbf{x}^{(k)}) + \widetilde{\mathbf{J}}^{(k)}(\mathbf{x}^{(k)} - \mathbf{x}_\ast)
```
or (assuming $\mathbf{I} - \widetilde{\mathbf{J}}^{(k)}$ to be invertible)
```math
\mathbf{x}^{(k)} - \mathbf{x}_\ast = \big(\mathbf{I} - \widetilde{\mathbf{J}}^{(k)}\big)^{-1} \mathbf{r}^{(k)}
```	 
where $\mathbf{I}$ is the identity matrix. Taking norms and using the matrix norm inequality we obtain
```math
\|\mathbf{x}^{(k)} - \mathbf{x}_\ast\| \leq \left\|\big(\mathbf{I} - \widetilde{\mathbf{J}}^{(k)}\big)^{-1}\right\| \ \|\mathbf{r}^{(k)}\|.
```

We notice
- As $k\to \infty$ then $\widetilde{\mathbf{J}}^{(k)} \to \mathbf{J}_g(\mathbf{x}_\ast)$.
- If $\mathbf{J}_g(\mathbf{x}_\ast)$ has eigenvalues near $1$ than $\mathbf{I} - \mathbf{J}_g(\mathbf{x}_\ast)$ has eigenvalues near zero. As a result $\left(\mathbf{I} - \mathbf{J}_g(\mathbf{x}_\ast)\right)^{-1}$ has very large eigenvalues (for more details see [Spectral transformations](https://teaching.matmat.org/numerical-analysis/11_Eigenvalue_problems.html#Spectral-transformations)) and thus the $\left\|\big(\mathbf{I} - \widetilde{\mathbf{J}}^{(k)}\big)^{-1}\right\|$ prefactor becomes large and a small residual norm again does not imply a small error.	 
""")

# ╔═╡ a377563c-63da-4a50-886b-002644085eab
md"""
!!! note "General principle: Residual"
    Note, that the **residual is a general terminology**,
    which is not only applied to such an error indicator in the context
    of fixed-point iterations, but used in general for iterative procedures.
    The idea is that the **residual provides
    the discrepancy from having fully solved the problem**
    and is thus a natural error indicator.
    The **functional form**, however, is **different for each type of iterative
    procedure** as we will see.
"""

# ╔═╡ 36860565-f44b-452b-a97f-e2ae356319cb
TODO("Vectorify from here")

# ╔═╡ 91060f66-9540-43dd-bffa-3ecc022efbdc
md"""
### Convergence order

When performing numerical methods one is usually not only interested whether an iteration converges, but also how quickly, i.e. how the error approaches zero.

!!! info "Definition: Convergence order and rate"
    A sequence $x^{(k)}$, which converges to $x_\ast$, is said to have
    **convergence order $q$** and **convergence rate $C$**
    when there exists a $q > 1$ and $C > 0$,
    such that
    ```math
    \tag{7}
    \lim_{k \to \infty} \frac{|x^{(k+1)} - x_\ast|}{|x^{(k)} - x_\ast|^q} = C
    ```
    If $q = 1$ (convergence order 1) we additionally require $0 < C < 1$.
"""

# ╔═╡ deff7837-e80e-4fcc-8eb2-baafc047a992
md"""
Some remarks:
  - Convergence order $q = 1$ is also called **linear convergence**,
    any convergence $q > 1$ is usually called **superlinear convergence**.
    In particular $q = 2$ is called **quadratic convergence**
  - These names become more apparent if one considers a logarithmic scale.
    Suppose for simplicity that $q=1$ and all ratios in (7) are equal to $C$
    (perfect linear convergence), then
    $|x^{(k+1)} - x_\ast| = α C^k$
    where $α$ is some constant. Taking the logs we get
    ```math
    \log |x^{(k)} - x_\ast| = k \log(C) + \log(α)
    ```
    which is a straight line.
"""

# ╔═╡ 7f2c43c0-a87e-44bb-ac52-77dc70302cd2
md"""
#### Visual inspection: Log of error

The last remark provides an idea how to **visually inspect the convergence order**,
namely by plotting the error $|x^{(k)} - x_\ast|$ on a logscale. For our fixed point iterations,
the (hopefully) converging sequence is exactly generated by the relationship $x^{(k+1)} = g(x^{(k)})$.

So let us inspect the convergence of $g_1$ and $g_2$ graphically:
"""

# ╔═╡ b2e96aa1-69a6-4d07-89b5-15ecc9a41398
let
	fp = 1.2079400315693  # Approximate fixed point
	p = plot(; yaxis=:log, xlabel="k", ylabel=L"|x^{(k)} - x_\ast|")

	results_g₁ = fixed_point_iterations(g₁, 4.0; maxiter=15)
	errors_g₁ = [abs(xn - fp) for xn in results_g₁.history_x]
	plot!(p, errors_g₁, label=L"g_1", mark=:x, lw=2)
	
	results_g₂ = fixed_point_iterations(g₂, 4.0, maxiter=15)
	errors_g₂ = [abs(xn - fp) for xn in results_g₂.history_x]
	plot!(p, errors_g₂, label=L"g_2", mark=:x, lw=2)

	yticks!(p, 10.0 .^ (-6:0))
end

# ╔═╡ 86f4c162-41bd-4e20-9e3f-681ff36bac75
md"""
So clearly both $g_1$ and $g_2$ converge linearly,
but $g_1$ has a larger convergence rate.
"""

# ╔═╡ dc136c3d-4534-4d66-86ca-36114bb825bb
md"""
#### Visual inspection: Residual ratio

One caveat with this analysis is that we cheated a little by assuming that we already *know* the solution. An alternative approach is to **build upon our 
residual-error relationship** (6), i.e.
```math
|x^{(k)} - x_\ast| = \frac{1}{|1 - g'(\xi^{(k)})|} r^{(k)}.
```
and investigate the limit of the **residual ratio**
```math
\lim_{k\to\infty}
\frac{|r^{(k+1)}|}{\big| r^{(k)} \big|^q}
= \lim_{k\to\infty} \frac{|1 - g'(\xi^{(k+1)})|}{|1 - g'(\xi^{(k)})|^q}
\frac{|x^{(k+1)} - x_\ast|}{|x^{(k)} - x_\ast|^q}
\stackrel{(7)}{=} \frac{1}{|1 - g'(x_\ast)|^{q-1}} C .
```
In other words if the residual ratios approach a constant for a chosen
$q$, then we have $q$-th order convergence.

In particular for **linear order** ($q=1$) we have
```math
\lim_{k\to\infty} \frac{|r^{(k+1)}|}{\big| r^{(k)} \big|} = C,
```
i.e. that the **residual ratio should approach a constant** as $k\to\infty$,
which is the **convergence rate**.

This is a condition we can check also
*without* knowing the solution:
"""

# ╔═╡ d22a7de4-67ef-4664-951b-8fb46116a7bc
let
	p = plot(xlabel="k", ylabel=L"|r^{(k+1)}| / |r^{(k)}|")

	results_g₁   = fixed_point_iterations(g₁, 4.0; maxiter=15)
	residuals_g₁ = results_g₁.history_r
    ratios_g₁  = [abs(residuals_g₁[i+1] / residuals_g₁[i])
	              for i in 1:length(residuals_g₁)-1]
	plot!(p, ratios_g₁, label=L"g_1", mark=:x, lw=2)

	results_g₂ = fixed_point_iterations(g₂, 4.0, maxiter=15)
	residuals_g₂ = results_g₂.history_r
    ratios_g₂  = [abs(residuals_g₂[i+1] / residuals_g₂[i])
	              for i in 1:length(residuals_g₂)-1]
	plot!(p, ratios_g₂, label=L"g_2", mark=:x, lw=2)
end

# ╔═╡ 314f733c-c3c2-4679-bb7b-c94b96b54961
md"""
Clearly in both cases these ratios become approximately constant as $k$ gets larger.

!!! info "Observations"
    - For **linear convergence**, the **error reduces** in each iteration by a **constant factor**
    - Looking at the **error norm** or **residual ratio**
      on a **`log`-scale** as the iteration proceeds
      is often extremely insightful to understand the convergence behaviour
      (and debug implementation bugs)!
"""

# ╔═╡ 68f63b13-1326-4cb6-8db1-3b043b2ad78e
md"""
## Optional: Bisection method

Fixed-point iterations are a very useful tool to solve non-linear equations.
However, their condition for convergence, namely $g'(x_\ast)$ is very hard
to verify *a priori*, i.e. before even attempting a numerical solution
of the problem we have at hand.

We will now develop a simple algorithm for root finding,
which has the appealing feature, that under an easily verifyable
condition, it is guaranteed to converge.

Going back to our diode circuit problem, we find that the non-linear, scalar fixed-point problem $g_\text{exp}(v_D) = v_D$ can be written as the root-finding problem $f(x) = 0$ for the function
```math
f(x) =  R \, i_0 \left( e^{x/v_0} - 1 \right) + x - V.
```
"""

# ╔═╡ b96fd4b2-13d1-4e42-ac6c-074d595f4750
md"We select the parameters"

# ╔═╡ 7839caac-64e9-443d-bd49-03930fbe7aba
begin
	i0 = 1.0
	v0 = 0.1
	R = 1.0
	V = 1.0
end;

# ╔═╡ 19327a6c-c552-4115-8d69-1a08491d3bc4
md"and first plot $f$ to get an idea of the function:"

# ╔═╡ 0a9b876e-31de-4199-9804-7837ff9fe8a1
f(x) = R*i0*(exp(x/v0) - 1) + x - V

# ╔═╡ ccc8345d-9db8-443e-837f-af631aa9f41c
plot(f, ylims=(-3, 1), xlims=(-1, 0.5), label=L"f", lw=2)

# ╔═╡ 71c04ce7-bb99-4094-a04c-f53a15606cca
md"""
This function has clearly a pronounced root near zero, where it changes its sign.
This observation is put on more rigorous footing by the following

!!! info "Theorem 2"
    Let $f$ be a continuous function defined on $[a, b]$
    with $f(a) f(b) < 0$ (i.e. the function takes different signs at the boundary
    of the interval). Then there exists a $x_\ast \in [a, b]$ such that $f(x_\ast) = 0$.

This is a direct consequence of the [intermediate value theorem](https://en.wikipedia.org/wiki/Intermediate_value_theorem),
which we now want to exploit to find a root numerically.

Suppose $f$ on $[a, b]$ satisfies the conditions of the theorem $f(a) f(b) < 0$
and let us consider the midpoint of the interval $x_m = \frac{a+b}{2}$.
There are three possible cases:
- If $f(x_m) f(a) < 0$, then there is a root in $[a, x_m]$
- If $f(x_m) f(b) < 0$, then root is in $[x_m, b]$
- If $f(x_m) = 0$, we found a root.

By simply repeating this procedure (now on the *smaller*
interval $[a, x_m]$ or $[x_m, b]$)
we obtain the **bisection method**:
"""

# ╔═╡ b89d2d83-7fb8-44da-b6c8-a4781b16a384
function bisection_method(f, a, b; tol=1e-6)
	@assert a ≤ b
	@assert f(a) * f(b) < 0  # Otherwise the assumptions are not true

	# Initialise
	k  = 0
	xᵏ = (a + b) / 2
	
	history_x = Float64[]  # Empty Array, but only for Float64 numbers
	while abs(b - a) / 2 ≥ tol
		k = k + 1
		if f(xᵏ) * f(a) < 0
			b = xᵏ  # New interval [a, xᵏ]
		else
			a = xᵏ  # New interval [xᵏ, b]
		end

		xᵏ = (a + b) / 2
		push!(history_x, xᵏ)
	end

	(; root=xᵏ, history_x, n_iter=k)
end

# ╔═╡ 54ff2cfe-6d22-4451-a496-f8f87e8bb1d7
md"Its convergence plot again suggests a linear convergence:"

# ╔═╡ 17e5586b-6c5f-4c3a-abba-5dc9ed8865f5
let
	# First get an almost exact root
	reference = bisection_method(f, -0.5, 0.5; tol=1e-12)

	# Now run again to plot convergence
	result = bisection_method(f, -0.5, 0.5)
	errors  = [abs(xn - reference.root) for xn in result.history_x]
	
	plot(errors, label="bisection error", mark=:x, lw=2, yaxis=:log)
end

# ╔═╡ 1bbe92e6-3a06-4c0f-bc7d-e2117da4f6c1
md"""
In our suggested implementation we choose to stop
either when the length of the interval $[a, b]$ drops below
the desired tolerance or when a maximal number of iterations is reached.
However, for the bisection method we could even determine *a priori*
how many iterations to perform as we will now demonstrate.

The bisection method starts from an interval $I^{(0)} = [a, b]$,
which has width $|I^{(0)}| = b - a$.
In each iteration we split the interval into two parts of equal size,
therefore the interval in iteration $k$ has size
$|I^{(k)}| = \frac{b - a}{2^k}$.
By construction the root $x_\ast \in I^{(k)}$.
Our estimate for the root is always the midpoint of $I^{(k)}$.
Therefore the error in the $k$-th step is bounded by
```math
\tag{7}
|x^{(k)} - x_\ast| \leq \frac12 |I^{(k)}| = \left(\frac{1}{2}\right)^{k+1} (b-a)
```
Suppose now we want to achieve an error less than $\epsilon$,
i.e. $|x^{(k)} - x_\ast| \leq ϵ$. According to (7) it is sufficient to achieve
```math
\tag{8}
\left(\frac{1}{2}\right)^{k+1} (b-a) \leq ϵ
```
in order to have an error below the threshold $\epsilon$.
Rearringing (8) leads to
```math
k > \underbrace{\log_2 \left( \frac{b - a}{ϵ} \right) - 1}_{=K},
```
which thus provides a lower bound on the number of iterations
we need to achieve convergence to $\epsilon$.
Note that $K$ can be computed *a priori* before even starting
the bisection algorithm.
Moreover since (7) and (8) are *guaranteed* bounds,
iterating for at least $K$ iterations guarantees
that an error below $ϵ$ is obtained.
In comparison to our residual-based stopping criterion
for the fixed point iterations, this is a much stronger result.

From our analysis we can characterise the bisection method as follows:
- If we can find an interval $[a, b]$ wherea given function $f$ changes sign, then the bisection method almost certainly converges to a root.[^2]
- We can precisely control the error up to the point where we know *a priori* how many iterations are needed.
- However, the algorithm cannot be employed if such an interval $[a, b]$ cannot be found.
"""

# ╔═╡ ac693757-bb8e-467e-b19d-c15ed7fd244a
md"""
!!! exercise "Exercise"
    Proove the linear convergence of the bisection method. What is the convergence rate $C$ ?
"""

# ╔═╡ 46620c64-a43c-4d08-bb3e-5fda3c8ff57c
md"""
[^2]: Our analysis does not include the effect of finite-precision floating point arithmetic, which in theory and in practice can inhibit convergence for some tricky cases.
"""

# ╔═╡ d365a1b6-ca5d-4719-99af-ce16b078b6b2
md"""
## Newton's method

### Achieving higher-order convergence

So far we we only constructed methods with linear convergence order.
In this subsection we first want to understand what is required
to achieve quadratic or higher-order convergence
and then use this to construct a second-order method.

Let us return to the fixed-point iterations
$x^{(k+1)} = g(x^{(k)})$ and revisit
the **Taylor expansion** (3) of the **fixed-point map** $g$.
Continuing the expansion to **second order**
we notice for the error
```math
\tag{9}
\begin{aligned}
x^{(k+1)} - x_\ast &= g(x^{(k)}) - x_\ast\\
&= g'(x_\ast) \, (x^{(k)} - x_\ast)
+ \frac12 g''(x_\ast) \, (x^{(k)} - x_\ast)^2
+ O(|x^{(k)} - x_\ast|^3).
\end{aligned}
```
**Assume** now that $g'(x_\ast) = 0$, such that
```math
x^{(k+1)} - x_\ast = \frac12 g''(x_\ast) \, (x^{(k)} - x_\ast)^2
+ O(|x^{(k)} - x_\ast|^3).
```
By neglecting the small terms and rearranging we observe
```math
\frac{x^{(k+1)} - x_\ast}{(x^{(k)} - x_\ast)^2} \simeq \frac12 g''(x_\ast)
```
when $x^{(k)}$ is close to $x_\ast$.
Comparing with the condition of order-$q$ convergence,
i.e. that the limit
```math
    \lim_{k \to \infty} \frac{|x^{(k+1)} - x_\ast|}{|x^{(k)} - x_\ast|^q} = C
```
is a constant, we thus would **expect such a fixed-point method**
with $g'(x_\ast) = 0$ to **give quadratic convergence**.
More generally if $g'(x_\ast) = g''(x_\ast) = \cdots = g^{(q-1)}(x_\ast) = 0$
and $g^{(q)}(x_\ast) \neq 0$ we obtain a method of order $q$.
We summarise in a theorem:
"""

# ╔═╡ ac5a79d7-07e0-499e-9dab-519ae0b5f638
md"""
!!! info "Theorem 3"
    Let $g : [a, b] \to \mathbb{R}$ be a $p$ times
    continuously differentiable fixed-point map
    with fixed point $x_\ast \in [a, b]$.
    If
    ```math
    g'(x_\ast) = g''(x_\ast) = \cdots = g^{(q-1)}(x_\ast) = 0
    \quad \text{and} \quad g^{(q)}(x_\ast) \neq 0
    ```
    then there exists a $δ>0$ such that for all
    starting points $x^{(0)} \in [x_\ast - δ, x_\ast + δ]$
    the following holds:
    * The **fixed-point iterations** $x^{(k+1)} = g(x^{(k)})$
      converge to $x_\ast$ with **convergence order $q$**.
    * The convergence rate is
      ```math
      \lim_{k\to\infty} \frac{|x^{(k+1)} - x_\ast|}{|x^{(k)} - x_\ast|^q}
      = \frac{1}{q!} \left|g^{(q)}(x_\ast)\right|
      ```

Recall the residual-error relationship (6)
```math
|x^{(k)} - x_\ast| = \frac{1}{|1 - g'(\xi^{(k)})|} r^{(k)}.
```
A corollary of our arguments is that for superlinear methods
we have that $g'(x_\ast) = 0$, such that for $x^{(k)}$ close to
$x_\ast$ (and thus $\xi^{(k)} \simeq x_\ast$) we have that
```math
|x^{(k)} - x_\ast| \simeq r^{(k)}.
```
As a result for superlinear methods a residual-based stopping criterion
becomes extremely reliable.
"""

# ╔═╡ 92f0e4d8-46b1-45ae-9c0a-1c3485713c94
md"""
### Construction of Newton's method

Newton's method and its variants are the **most common approaches**
to **solve non-linear equations** $f(x) = 0$.

To develop these methods assume that the fixed-point is $x_\ast$
and we are given a point $x$, which is close to $x_\ast$.
We consider a Taylor expansion of $f$ around $x$:
```math
0 = f(x_\ast) = f(x) + f'(x) (x_\ast - x) + O\big( (x_\ast - x)^2 \big)
```
Where we made the condition $f(x_\ast) = 0$ has been made explicict.
Now assume $f'(x) \neq 0$ to rearrange this to
```math
0 = \frac{f(x)}{f'(x)} + (x_\ast - x) + O\big( (x - x_\ast)^2 \big).
```
If $x$ is close to $x_\ast$, thus $x - x_\ast$ is small,
then the last $O \big((x - x_\ast)^2\big)$ term is even smaller.
Neglecting it we can further develop this to an approximation for $x_\ast$ as
```math
0 \simeq \frac{f(x)}{f'(x)} + x_\ast - x
\qquad
\Longrightarrow
\qquad
x_\ast \simeq x - \frac{f(x)}{f'(x)}.
```
If we denote the RHS by $g_\text{Newton}(x) = x - \frac{f(x)}{f'(x)}$,
then a root of $f$ is a fixed point of $g_\text{Newton}$
and vice versa.
Performing fixed-point iterations on $g_\text{Newton}$
is the idea of Newton's method:

!!! info "Algorithm 1: Newton's method (fixed-point formulation)"
    Given a  $C^1$ function $f$ and an initial guess $x^{(0)}$
    perform fixed-point iterations
    ```math
    x^{(k+1)} = g_\text{Newton}(x^{(k)})
    ```
    on the map
    ```math
    \tag{10}
    g_\text{Newton}(x) = x - \frac{f(x)}{f'(x)}.
    ```
"""

# ╔═╡ 7930a403-bd4d-4f64-a79b-54cff28b869b
md"""
### Graphical interpretation
- See the chapter on [Newton's method](https://computationalthinking.mit.edu/Fall23/images_abstractions/newton_method/)
  in the MIT computational thinking class.

- See [chapter 4.3](https://tobydriscoll.net/fnc-julia/nonlineqn/newton.html)
  of Driscoll, Brown: *Fundamentals of Numerical Computation*.
"""

# ╔═╡ 06d37eb4-b3e6-4219-82bd-fc9d779885c9
md"""
### Convergence analysis

Our goal is to apply Theorem 3 in order to obtain both the result
that Newton's method converges as well as an understanding of its
convergence order. We thus study the
derivatives of $g_\text{Newton}$ at the fixed point $x_\ast$.
We obtain
```math
\begin{aligned}
g_\text{Newton}'(x)
&= 1 - \frac{f'(x)\,f'(x) - f''(x)\,f(x)}{\big(f'(x)\big)^2}
= \frac{f''(x)\,f(x)}{\big(f'(x)\big)^2}\\
g_\text{Newton}''(x)
&= \frac{
\big(f'(x)\big)^3\,f''(x)
+ f(x)\, \big(f'(x)\big)^2\, f'''(x)
- 2 f(x) \, \big(f''(x)\big)^2 \, f'(x)
}{\big(f'(x)\big)^4}
\end{aligned}
```
such that under the assumption that $f'(x_\ast) \neq 0$
and $f''(x_\ast) \neq 0$ we obtain
```math
\begin{aligned}
g_\text{Newton}'(x_\ast) &= 0\\
g_\text{Newton}''(x_\ast) &= \frac{f''(x_\ast)}{f'(x_\ast)}  \neq 0
\end{aligned}
```
where we used $f(x_\ast) = 0$.
We summarise

!!! info "Theorem 4: Convergence of Newton's method" 
    Let $f$ be a twice differentiable ($C^2$) function and $x_\ast$ a root of $f$.
    If $f'(x_\ast) \neq 0$ and $f''(x_\ast) \neq 0$, then
    Newton's method converges quadratically for every $x^{(0)}$
    sufficiently close to $x_\ast$. The rate is
    ```math
    \lim_{k\to\infty} \frac{|x^{(k+1)} - x_\ast|}{|x^{(k)} - x_\ast|^2}
    = \frac12 \left|\frac{f''(x_\ast)}{f'(x_\ast)}\right|
    ```
Some remarks:
- Theorem 4 only makes a *local* convergence statement,
  i.e. it requires the initial value $x^{(0)}$ to be close enough
  to $x_\ast$.
- If $f'(x_\ast) = 0$ we can show that Newton's method is only
  of first order.
"""

# ╔═╡ 0d4425d5-964a-4ce5-aab6-0b388b523ee3
md"""
### Implementation

Since in our construction Newton's method (Algorithm 1) is obtained in form of the fixed-point map $g_\text{Newton}$ the implementation is straightforward by employing the `fixed_point_iterations` function we already implemented above:
"""

# ╔═╡ 88cf8bbf-55cd-45f0-9ab4-ae2bdf6c554a
TODO(md"""Code up this next algorithm live""")

# ╔═╡ 929010af-2af6-4b63-bfb4-7987f1b790ca
function newton_fp(f, df, xstart; maxiter=40, tol=1e-6)
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

# ╔═╡ 2b04681b-6704-4172-9bf1-279051fc1d78
md"""
In this setting where $g_\text{Newton}$ exhibits quadratic convergence,
the **residual-based stopping criterion** of `fixed_point_iterations` is actually
**extremely reliable**, see the discussion after Theorem 3.
"""

# ╔═╡ 0d61a532-e185-49bc-baf9-c22b8443a979
begin
	text_algo2_newton = md"""
!!! info "Algorithm 2: Newton's method (conventional)"
	Given a once differentiable function $f : \mathbb{R} \to \mathbb{R}$,
	a starting value $x^{(0)}$ and a convergence tolerance $\epsilon$,
    perform for $k = 1, 2, 3, \ldots$:
    1. Compute the residual $r^{(k)} = - \frac{f(x^{(k)})}{f'(x^{(k)})}$
    2. Update ${x}^{(k+1)} = {x}^{(k)} + {r}^{(k)}$
    Loop 1. and 2. until $|{r}^{(k)}| < \epsilon$.
"""
	
md"""
	More conventionally one "inlines" the function $g_\text{Newton}$ into the fixed point iterations and expresses the problem as $text_algo2_newton
	A corresponding implementation of Algorithm 2 is:
	"""
end

# ╔═╡ d8091ef9-7350-45e9-8381-55c9297429cd
function newton(f, df, xstart; maxiter=40, tol=1e-6)
	# f:  Function of which we seek the roots
	# df: Function, which evaluates its derivatives
	# xstart: Start of the iterations
	# maxiter: Maximal number of iterations
	# tol: Convergence tolerance

	history_x = [float(xstart)]
	history_r = empty(history_x)

	r = Inf     # Dummy to enter the while loop
	x = xstart  # Initial iterate
	k = 0

	# Keep running the loop when the residual norm is beyond the tolerance
	# and we have not yet reached maxiter
	while norm(r) ≥ tol && k < maxiter
		k = k + 1

		# Evaluate function, gradient and residual
		r = - f(x) / df(x)

		# Evaluate next iterate
		x = x + r
		
		push!(history_r, r)  # Push residual and
		push!(history_x, x)  # next iterate to history
	end

	(; root=x, n_iter=k, history_x, history_r)
end

# ╔═╡ ecd71729-ea89-40f6-86c8-598e42cd787d
md"""
To compare Algorithm 2 to Algorithm 1 note that steps 1 and 2 jointly apply the function $g_\text{Newton}$ and that the residual is defined as
```math
r^{(k)} = g_\text{Newton}(x^{(k)}) - x^{(k)} =  x^{(k)} - \frac{f(x^{(k)})}{f'(x^{(k)})} - x^{(k)} = - \frac{f(x^{(k)})}{f'(x^{(k)})}
```
"""

# ╔═╡ baabfa2d-04b8-4305-9a73-2008cd018a09
md"""
To see how this method performs we compare against bisection
and the plain fixed-point iterations in $g_\text{log}$ we saw earlier.
"""

# ╔═╡ e50fb224-63c1-4885-8165-21491110b804
md"""
```math
g_\text{log}(x) = v_0 \log\left( \frac{V - x}{R \, i_0} + 1 \right)
```
Let's code up this function in Julia:
"""

# ╔═╡ 03d48dc9-a7a6-480b-a755-2c4e049d3f7b
function glog(vD)
	v0 * log((V - vD) / (R * i0) + 1)
end

# ╔═╡ 2e1751b5-a2cd-4a39-a7a8-4601f4448ff4
let
	# First get an almost exact root
	reference = bisection_method(f, -0.5, 0.5; tol=1e-14)

	p = plot(yaxis=:log, xlims=(0, 20), ylims=(1e-12, Inf),
	         xlabel="k", ylabel=L"|x^{(k)} - x_\ast|")

	# Run bisection
	result = bisection_method(f, -0.5, 0.5; tol=1e-12)
	errors = [abs(xn - reference.root) for xn in result.history_x]
	plot!(p, errors, label="Bisection", mark=:x, lw=2)

	# Run fixed-point on glog
	result = fixed_point_iterations(glog, 0.0; tol=1e-12)
	errors  = [abs(xn - reference.root) for xn in result.history_x]
	plot!(p, errors, label="Fixed point", mark=:x, lw=2)

	# For Newton we need the derivative of f.
	# An easy way to obtain this derivative is to use algorithmic differentiation:
	df(x) = ForwardDiff.derivative(f, x)

	# With this we run Newton
	result = newton(f, df, 0.0; tol=1e-12)
	errors = [abs(xn - reference.root) for xn in result.history_x]
	plot!(p, errors, label="Newton", mark=:x, lw=2)
	yticks!(p, 10.0 .^ (-12:-1))

	p
end

# ╔═╡ 28088d7f-cbcb-4c17-95bd-c3790da34ae2
md"""
On the log-scale of the plot the quadratic convergence behaviour of Newton's method is clearly visible.

Let us investigate a little the stability of this algorithm, especially with respect
to the requirement to choose a sufficiently close initial guess:
"""

# ╔═╡ 82286823-1240-41b5-92d9-9b3e091129e5
md"
- `x0 = ` $(@bind x0 Slider(-1:0.1:3; default=2.5, show_value=true))
"

# ╔═╡ e1a1e4a7-21b6-4a88-92f9-a4d3e0dacbb9
let
	f(x)  = 0.5x^3 - x^2 - 2x
	df(x) = 1.5x^2 - 2x - 2
	result = newton(f, df, x0; tol=0)

	p = plot(f, xlims=(-4, 5), label=L"f", lw=2, legend=:topleft)
	hline!(p, [0.0], color=:grey, lw=1, label="")

	colors = palette(:algae, 6; rev=true)
	for (i, x) in enumerate(result.history_x)
		i > 6 && break
		vline!(p, [x], color=colors[i], lw=2, ls=:dash, label="Iteration $(i-1)")
	end
	scatter!(p, [result.root], [0.0];
	         mark=:o, color=colors[end], label="Root")
	scatter!(p, [x0], [0.0], mark=:o, color=colors[1], label="")

	p
end

# ╔═╡ e5725d2b-527f-4a50-9a2a-89d0514ae6bb
md"""
Finally, let us investigate what quadratic convergence
means in terms of our error estimate,
the residuals of the Newton fixed-point map. Since Newton converges
fast and thus very quickly maxes out the roughly 16 digits of precision
in standard `Float64` numbers, we employ Julia's arbitrary precision
`BigFloat` type to see more closely what is going on:
"""

# ╔═╡ b5a13fb0-6d04-411f-ab9e-5ed61f2451e5
let
	# We want to solve e^x * x = 2, which has a solution near 1.
	f(x)  = x*exp(x) - 2
	df(x) = (x+1) * exp(x)

	xstart = BigFloat(1.0)
	result = newton(f, df, xstart; tol=1e-60)

	for (i, r) in enumerate(result.history_r)
		@printf "%3i %e\n" i r
	end
end

# ╔═╡ 9be996e7-9a33-4399-9043-28d1f93cd890
md"We observe that from the starting point only $7$ iterations are required to get the result accurate to 61 digits."

# ╔═╡ efcd0f6a-4c48-49c8-9c7a-f49421170a33
md"""
!!! info "Obervations"
    - For quadratic convergence the error roughly squares in each iteration
    - Newton only converges well if the initial guess is chosen sufficiently
      close to the fixed point.
	- Unlike the bisection algorithm the convergence behaviour of Newton
	  is thus sometimes less reliable. In particular if $f$ has **multiple roots** it is
	  **not guaranteed**, that **Newton converges to the *closest* root**.
	  A good graphical representation of this phaenomenon
	  are [Newton fractals](https://en.wikipedia.org/wiki/Newton_fractal).
"""

# ╔═╡ 04cbf165-4786-4980-8055-5efca14dfa3f
md"""
## Lessons to learn from fixed-point iterations

In this chapter we discussed fixed point iterations as well as Newton's method (Algorithm 2) as two examples for iterative algorithms.

!!! info "General form of iterative algorithms"
	More generally an iterative algorithm has the form:
	1. **Initialisation:** Choose a starting point $x^{(0)} = x_\text{start}$.
	2. **Iteration:** For $k = 1, 2, 3, ...$ we perform the same
	  iterative procedure, advancing $x^{(k-1)}$ into $x^{(k)}$.
	3. **Check for convergence:** Once $x^{(k)}$ is similar enough to $x^{(k-1)}$ we consider the iteration converged and exit the iterations.

Writing the second step more mathematically we can consider it as the application of a function $g$, i.e. $x^{(k)} = g(x^{(k-1)})$. In this formulation step 3 thus does nothing else than checking whether the iterates $x^{(k)}$ no longer change. Or, put in other words, if we have found a **fixed point** of $g$.
"""

# ╔═╡ a0390836-e86b-46fb-bb42-fe6e5b09bd64
md"""
!!! info "Observation: Iterative algorithms are fixed-point problems"
	By identifying the iteration step $x^{(k)} = g(x^{(k-1)})$ of *any*
	iterative algorithm with a function $g$ we can view this algorithm
	as a fixed-point problem, where at convergence a fixed-point of $g$
	is found. 

	As a result **any technique we discussed** for understanding *when
	fixed-point iterations converge* and at *which convergence rate*
	can be used **in general** to analyse the convergence of
	*any iterative procedure*. 

We will consider this aspect further,
for example in [Iterative methods for linear systems](https://teaching.matmat.org/numerical-analysis/06_Iterative_methods.html).
"""

# ╔═╡ bdff9554-58b6-466e-9c93-6b1367262b50
md"""
## Optional: Secant method

See [chapter 4.4](https://tobydriscoll.net/fnc-julia/nonlineqn/secant.html)
of Driscoll, Brown: *Fundamentals of Numerical Computation*.

"""

# ╔═╡ 6c823b3d-2332-4299-8fef-bc392d13da94
md"""
## Non-linear equation systems
"""

# ╔═╡ 74ff4f4f-66a9-4de7-826b-4cc66d445a73
md"""
!!! warning "Running example: Definition"
	As the running example in this section we will consider the problem 
	$\textbf{f}(\textbf{x}) = \textbf{0}$
	with $\textbf{f} : \mathbb{R}^3 \to \mathbb{R}^3$ given by
	```math
	\left\{ \begin{aligned}
	f_1(x_1, x_2, x_3) &= -x_1 \cos(x_2) - 1\\
	f_2(x_1, x_2, x_3) &= x_1 \, x_2 + x_3\\
	f_3(x_1, x_2, x_3) &= e^{-x_3} \, \sin(x_1 + x_2) + x_1^2 - x_2^2
	\end{aligned} \right. .
	```
"""

# ╔═╡ 0449d301-0896-4c2e-952c-0d380e60c0b9
md"""
This example very much illustrates the previous point of the increased complexity:
Even guessing an approximate solution is not obvious,
so we proceed to develop a numerical technique.
Taking inspiration from Newton's method in 1D
we proceed to solve such equation systems **by linearisation**.
That is to say we start by developing $\textbf{f}$ to first order
around some initial point $\textbf{x}$,
which we assume to be close enough to the solution $\textbf{x}_\ast$
```math
\tag{11}
\textbf{0} = \textbf{f}(\textbf{x}_\ast)
= \textbf{f}(\textbf{x}) + \textbf{J}_\textbf{f}(\textbf{x}) \, (\textbf{x}_\ast - \textbf{x})
+ O(\| \textbf{x}_\ast - \textbf{x} \|^2).
```
In this the **Jacobian matrix**
$\textbf{J}_\textbf{f}(\textbf{x}) \in \mathbb{R}^{n\times n}$
is the collection of all partial derivatives of $\textbf{f}$, i.e.
```math
\textbf{J}_\textbf{f} = \left(\begin{array}{cccc}
\frac{\partial f_1}{\partial x_1} &
\frac{\partial f_1}{\partial x_2} &
\ldots &
\frac{\partial f_1}{\partial x_n} \\

\frac{\partial f_2}{\partial x_1} &
\frac{\partial f_2}{\partial x_2} &
\ldots &
\frac{\partial f_2}{\partial x_n} \\

\vdots & & \ddots & \vdots\\

\frac{\partial f_n}{\partial x_1} &
\frac{\partial f_n}{\partial x_2} &
\ldots &
\frac{\partial f_n}{\partial x_n}
\end{array}\right).
```
See also the discussion on multi-dimensional Talyor approximations in [Revision and preliminaries](https://teaching.matmat.org/numerical-analysis/03_Preliminaries.html).


The Jacobian very much plays the role of a generalised derivative
of a multidimensional function $\textbf{f}$.
Also not that (just like any derivative) it is a function of the independent
variable $\textbf{x}$.

!!! warning "Running example: Computing the Jacobian"
	For the running example identified above we can compute the Jacobian as
	```math
	\textbf{J}_\textbf{f}(\textbf{x})
	= \left(\begin{array}{ccc}
	-\cos(x_2) & x_1 \sin(x_2) & 0 \\
	x_2 & x_1 & 1 \\
	e^{-x_3} \, \cos(x_1 + x_2) + 2x_1 &
	e^{-x_3} \, \cos(x_1 + x_2) - 2x_2 &
	- e^{-x_3} \, \sin(x_1 + x_2)
    \end{array}\right)
    ```
"""

# ╔═╡ 4cf5f022-4047-49f1-8fc5-86eda6d61e30
md"""
In expansion (11) the terms $\textbf{f}(\textbf{x}) + \textbf{J}_\textbf{f}(\textbf{x}) \, (\textbf{x}_\ast - \textbf{x})$
represent the linear part of $\textbf{f}$ around $\textbf{x}_\ast$.
Assuming that these dominate, i.e. that the remaining term $O(\| \textbf{x}_\ast - \textbf{x} \|^2)$ is indeed small and can be neglected, we obtain:
```math
\textbf{0} \simeq \textbf{f}(\textbf{x}) + \textbf{J}_\textbf{f}(\textbf{x}) \, (\textbf{x}_\ast - \textbf{x})
```
In the same spirit as in the 1D Newton case we want to employ this relation
in an iterative scheme, where in iteration $k$ we have $\textbf{x}^{(k)}$
and want to compute an improved iterate $\textbf{x}^{(k+1)}$.
Inserting
$\textbf{x}^{(k)}$ as $\textbf{x}$ 
and
$\textbf{x}^{(k+1)}$ as $\textbf{x}_\ast$ 
as in the 1D case, we obtain
```math
\textbf{0} =
\textbf{f}(\textbf{x}^{(k)}) + \textbf{J}_\textbf{f}(\textbf{x}^{(k)}) \, (\textbf{x}^{(k+1)} - \textbf{x}^{(k)})
```
or
```math
\tag{12}
\Big( \textbf{J}_\textbf{f}(\textbf{x}^{(k)}) \Big) \, (\textbf{x}^{(k+1)} - \textbf{x}^{(k)}) = - \textbf{f}(\textbf{x}^{(k)}).
```
Assuming that $\det \textbf{J}_\textbf{f}(\textbf{x}^{(k)}) \neq 0$,
i.e. that the Jacobian is non-singular,
this linear system can be solved and thus $\textbf{x}^{(k+1)}$ computed.
"""

# ╔═╡ 7e3558c8-b16d-47a4-87ea-d2af5614ea2e
md"""
Note, that in accordance with the 1D case the **residual in this
multi-dimensional version** is $\textbf{r}^{(k)} = \textbf{x}^{(k+1)} - \textbf{x}^{(k)}$,
i.e. exactly the solution to the linear system in (12).
Similar to the 1D case we will thus employ the stopping criterion
$\|\textbf{r}^{(k)}\| < \epsilon$,
where $\|x\|$ denotes the Euclidean norm
$\|x\| = \sqrt{x_1^2 + x_2^2 + \cdots + x_n^2}$.
All combined we obtain the algorithm
"""

# ╔═╡ fa140f9b-af22-4060-9a66-027189f80dbc
md"""
!!! info "Algorithm 3: Multidimensional Newton's method"
    Given a once differentiable function $\textbf{f} : \mathbb{R}^n \to \mathbb{R}^n$,
	a starting value $\textbf{x}^{(0)}$ and a convergence tolerance $\epsilon$,
    perform for $k = 1, 2, 3, \ldots$:
    1. Compute
       the right-hand side $\textbf{y}^{(k)} = \textbf{f}(\textbf{x}^{(k)})$
       and Jacobian
       $\textbf{A}^{(k)} = \textbf{J}_\textbf{f}(\textbf{x}^{(k)})$.
    2. **Newton step:** Solve the linear system $\textbf{A}^{(k)} \textbf{r}^{(k)} = - \textbf{y}^{(k)}$ for $\textbf{r}^{(k)}$.
    3. Update $\textbf{x}^{(k+1)} = \textbf{x}^{(k)} + \textbf{r}^{(k)}$
    Loop 1. to 3. until $\|\textbf{r}^{(k)}\| < \epsilon$.
"""

# ╔═╡ fc9d11e1-b3b6-4d74-9cbe-1faf2d8f3abc
md"""
An implementation of this algorithm is given below:
"""

# ╔═╡ 4e746757-26a0-45c4-a190-04dd92b3cab0
function newtonsys(f, jac, xstart; maxiter=40, tol=1e-8)
	history_x = [float(xstart)]
	history_r = empty(history_x)

	r = Inf     # Dummy to enter the while loop
	x = xstart  # Initial iterate
	k = 0
	while norm(r) ≥ tol && k < maxiter
		k = k + 1
		
		y = f(x)      # Function value
		A = jac(x)    # Jacobian
		r = -(A \ y)  # Newton step
		x = x + r     # Form next iterate
		
		push!(history_r, r)  # Push newton step and
		push!(history_x, x)  # next iterate to history
	end

	(; root=x, n_iter=k, history_x, history_r)
end

# ╔═╡ f1850ccb-e682-40db-91f1-ca2c4f60b3f3
md"""
Note that he linear system $\textbf{A}^{(k)} \textbf{r}^{(k)} = - \textbf{y}^{(k)}$ is solved in Julia using the backslash operator `\`, which employs a numerically more stable algorithm than explicitly computing the inverse `inv(A)` and then applying this to `y`.
We will discuss these methods in
[Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/05_Direct_methods.html).
"""

# ╔═╡ 702ffb33-7fbe-4673-aed7-d985a76b455a
Foldable("Remark: Connection to conventional 1D Newton algorithm (Algorithm 2)",
md"""
Algorithm 3 (Multidimensional Newton) is just a multi-dimensional generalisation of Algorithm 2 (conventional 1D Newton's method). 

This can be easily seen by simplifying Algorithm 3 for the special case of a 1D problem, that is a function $f : \mathbb{R} \to \mathbb{R} : x_1 \mapsto f(x_1) $. In that case we can identify the Jacobian
(which formally is a $1 \times 1$ matrix) by the single element
$\frac{\partial f}{\partial x_1} = f'(x_1)$. With this in mind the three steps of Algorithm 3 become: For $k = 1, 2, 3, \ldots$
1. Compute the right-hand side $y^{(k)} = f(x^{(k)})$ and $A^{(k)} = J_f(x^{(k)}) = f'(x^{(k)})$.
2. Solve $A^{(k)} r^{(k)} = -y^{(k)}$ for $r^{(k)}$. Since all entities are real numbers this can be done by simple division, i.e. $r^{(k)} = - \frac{y^{(k)}}{A^{(k)}} = -\frac{f(x^{(k)})}{f'(x^{(k)})}$.
3. Update $x^{(k+1)} = x^{(k)} + r^{(k)} = x^{(k)} - \frac{f(x^{(k)})}{f'(x^{(k)})}$.

For reference recall Algorithm 2 was $text_algo2_newton

We note that the algorithms are idential, in essence Algorithm 3 only uses two steps for what Algorithm 2 expressed as step "1.".
""")

# ╔═╡ 6192f3a8-b9ff-4187-9561-00c2638a483a
md"""
Let's apply `newtonsys` (Algorithm 2) to our running example. First we implement the functions computing $\textbf{f}(\textbf{x})$ and $\textbf{J}_\textbf{f}(\textbf{x})$ for a given $\textbf{x}$.
"""

# ╔═╡ 1307929a-f151-4868-a43a-85488599f0df
begin
	func(x) = [
		-x[1] * cos(x[2]) - 1,
		x[1] * x[2] + x[3],
		exp(-x[3]) * sin(x[1] + x[2]) + x[1]^2 - x[2]^2
	]
	jac_func(x) = [
		-cos(x[2])       x[1]*sin(x[2])    0;
		x[2]             x[1]              1;
		exp(-x[3])*cos(x[1]+x[2]) + 2x[1]  exp(-x[3])*cos(x[1]+x[2]) - 2x[2] exp(-x[3])*sin(x[1]+x[2]) 
	]
end;

# ╔═╡ 8919edbb-5ea7-4956-ad75-3eb68a47408c
md"""
Since we want to estimate the convergence order we again
run the Newton solver using arbitrary precision floating-point numbers
by using Julia's `BigFloat` number type:
"""

# ╔═╡ b6510e80-6b64-421e-b8aa-bfa803d56e81
res = newtonsys(func, jac_func, BigFloat.([1.5, -1.5, 5]), tol=1e-50);

# ╔═╡ 42216744-eaed-4bb5-b882-dc4589d8e63d
md"""
Plotting the residual norm (our estimate of the error) in a log-plot gives a strong indication this is again quadratic convergence:
"""

# ╔═╡ 9bc67152-9220-4905-9d2b-ec3f37b43707
plot(norm.(res.history_r); yaxis=:log, label="", xlabel=L"k", ylabel=L"\Vert  e^{(k)} \Vert \simeq \Vert r^{(k)} \Vert")

# ╔═╡ b6edac36-ee05-4efa-8d75-7840e9e185e4
md"""
Using the residual norms stored in the Newton result, we can now also
look at the ratios
```math
\frac{\| \textbf{x}^{(k+1)} - \textbf{x}^{(k)} \| }
{\| \textbf{x}^{(k)} - \textbf{x}^{(k-1)} \|^q }
= \frac{\| \textbf{r}^{(k)} \| }{\| \textbf{r}^{(k-1)} \|^q }
```
of two consecutive increments.
Recall that for a $q$-th order convergence these should converge to a constant.
"""

# ╔═╡ 95120495-15c7-4567-84f1-cb6f535c3885
let
	for q in (1, 2, 3)
		println("# Checking order q=$q")
		for k in 2:length(res.history_r)
			ratio = norm(res.history_r[k]) / norm(res.history_r[k-1])^q
			@printf "%i  %.5f\n" k ratio
		end
		println()
	end
end

# ╔═╡ 6686d3bc-abc8-4a3a-a3cd-f6cac50fdf34
md"""
As can be see the most constant is the sequence corresponding to $q=2$,
such that we conclude that the method converges quadratically.
"""

# ╔═╡ 6d25666f-68c2-4a15-bf01-732fa2ebd4e3
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 530)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
ForwardDiff = "~1.3.1"
HypertextLiteral = "~1.0.0"
LaTeXStrings = "~1.4.0"
Plots = "~1.41.4"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "3211ff2f9b6aeb94c22c3c1a9129467822a9c79b"

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
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

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
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

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
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "b2977f86ed76484de6f29d5b36f2fa686f085487"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.1"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "f2fd0ae89599c11473fbd846ea1b9ed8c24613fb"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.21"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "55965bc822094d9e5d314d2845fb77f1b74fc8f7"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.21+0"

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
git-tree-sha1 = "6b4d2dc81736fe3980ff0e8879a9fc7c33c44ddf"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.2+0"

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
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

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
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

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
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

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
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

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
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

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

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39a11854f0cba27aa41efaedf43c77c5daa6be51"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.0+0"

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
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

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
git-tree-sha1 = "063ef757a1e0e15af77bbe92be92da672793fd4e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.4"

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
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

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
git-tree-sha1 = "34f7e5d2861083ec7596af8b8c092531facf2192"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "8f528b0851b5b7025032818eb5abbeb8a736f853"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+2"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
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
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

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
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

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
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

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
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

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
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

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
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

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
# ╟─0441b7ed-b62f-417c-a4f6-33051faac08e
# ╠═62ff56f9-caec-4352-b547-973ce2bfcb8f
# ╟─b258697b-4c4f-4e47-8472-73b72178c108
# ╟─6755e1bd-d225-46e6-b80f-fd7287c9dcef
# ╟─b03b35dd-13c3-4684-88ef-b74c91e4a8dc
# ╟─4e3ba905-31dd-4d29-9e45-4a77228efea0
# ╠═3e12bce6-701d-4ed1-baf3-13fe888ede07
# ╟─beba7886-fb36-45b1-b416-689733aecf4b
# ╟─126314f4-8488-4338-a419-fc87e8965ee3
# ╟─65e1a532-64c3-427b-a5ec-3fcf0f7eccec
# ╟─c454473a-a5ef-4b08-af45-da4a7acbb1d0
# ╟─d3ba1da7-41eb-4afd-a949-19123b7a6724
# ╟─e6c13ceb-fb13-4159-b0bc-a854406d9f21
# ╟─8b9bae6b-02ae-455e-a6f0-a10ca6db694e
# ╟─9f80abeb-432c-4b68-a570-9e32e2510fe0
# ╠═cabca1d2-e287-4d61-9d03-16108e10ac9c
# ╠═c8cb64e6-a12a-4308-89bd-73bf178c4986
# ╟─b4a0754f-dbb3-4de1-ba94-35832ece1048
# ╟─aec943b1-5865-4306-a780-fd730afd208a
# ╟─1b35b5da-7e1d-42d9-a189-c20d6957e939
# ╟─08578e28-28c0-4cb1-9934-734034dec251
# ╟─2f670892-be3e-43f3-b49f-803d8ddf756e
# ╟─b385d67b-c7e0-4a5e-87e1-b11f00ba638a
# ╟─c6d53a9e-0b6d-4ad3-a7c8-c34070f75ff6
# ╟─1319313d-e7da-4ff5-9725-548df0e7f699
# ╠═858510eb-7768-4521-ba20-af5b1aa33e08
# ╟─7bdb39f4-4fab-43a1-a9d9-12e198cef107
# ╟─5e21e4c5-5538-40cf-a955-837809f7c3c3
# ╟─05989b7d-e694-4860-9a8e-6e6cb52aea8b
# ╟─01db98ec-daf2-4779-9f31-c3271039f44c
# ╟─9176b666-41f7-436e-b5ad-61b196a8b35b
# ╟─511fcf1a-5ea8-4827-9d24-6a43e7dcccd6
# ╟─c7068f2e-7313-44bc-85a4-785f4d4adc60
# ╟─8a8ac29f-f664-4c62-bda7-eeb5214b2bc8
# ╟─d4e3f192-583e-4d3f-af7d-0767d2c9c0a1
# ╟─9c9719e3-ec6c-4bdc-b05b-ab4bd4119cb9
# ╟─fe6a0ff6-f70d-404a-8fd5-8a9185da2ee7
# ╟─33c6b6ed-267e-4062-8390-fa02b8624545
# ╟─29931ae9-bcb7-4ec0-b397-a89c491d950e
# ╠═5409350b-4845-46e5-ac3a-f7342fc28d3d
# ╠═9d7ca2e1-41b0-4651-a734-31bf883cee37
# ╟─fffd0bd2-7b66-4a21-8311-2965e974d88b
# ╠═50639d02-55d5-4fcb-8335-13fd7f6b7624
# ╠═da60beec-74e7-4b3f-aa09-27b806054896
# ╟─ec8b8273-3724-4ea9-91d3-259390abc55d
# ╟─5d7b3d35-3456-48df-ad22-0ffccaa2f529
# ╟─284f3fa4-ce24-4b99-8bb7-6f74a4589550
# ╠═899698d1-9dee-4f82-9171-f1b49aefcabe
# ╟─072cdb08-b279-4583-b98b-4a28433b5b8d
# ╠═3c0d4752-de1a-4819-a63f-aa49c94ccf82
# ╟─b5ec49d7-1876-404e-b6ed-03217612b5e4
# ╠═58cf700e-8bb4-4370-ae55-875d551a3308
# ╟─0eb10646-4830-4d97-b574-54438a20a3d9
# ╠═b422909a-2f84-4815-b2f3-9902482806e6
# ╟─c1b1c4ff-e2c9-49fb-8794-b9ee11c5d543
# ╠═7cb34ba9-bcfc-4596-8eb3-49f4b99eab58
# ╠═54f314c1-d1cf-41f1-96e5-5aca90d82b95
# ╟─73f0c043-24a6-401f-8695-7ac2ca154a7c
# ╟─62ec67d7-10d4-4624-8dd7-d4dd0507d6c2
# ╟─de918a0d-05ee-40b3-8cbd-f45274c37ed3
# ╟─a377563c-63da-4a50-886b-002644085eab
# ╠═36860565-f44b-452b-a97f-e2ae356319cb
# ╟─91060f66-9540-43dd-bffa-3ecc022efbdc
# ╟─deff7837-e80e-4fcc-8eb2-baafc047a992
# ╟─7f2c43c0-a87e-44bb-ac52-77dc70302cd2
# ╠═b2e96aa1-69a6-4d07-89b5-15ecc9a41398
# ╟─86f4c162-41bd-4e20-9e3f-681ff36bac75
# ╟─dc136c3d-4534-4d66-86ca-36114bb825bb
# ╠═d22a7de4-67ef-4664-951b-8fb46116a7bc
# ╟─314f733c-c3c2-4679-bb7b-c94b96b54961
# ╟─68f63b13-1326-4cb6-8db1-3b043b2ad78e
# ╟─b96fd4b2-13d1-4e42-ac6c-074d595f4750
# ╠═7839caac-64e9-443d-bd49-03930fbe7aba
# ╟─19327a6c-c552-4115-8d69-1a08491d3bc4
# ╠═0a9b876e-31de-4199-9804-7837ff9fe8a1
# ╠═ccc8345d-9db8-443e-837f-af631aa9f41c
# ╟─71c04ce7-bb99-4094-a04c-f53a15606cca
# ╠═b89d2d83-7fb8-44da-b6c8-a4781b16a384
# ╟─54ff2cfe-6d22-4451-a496-f8f87e8bb1d7
# ╠═17e5586b-6c5f-4c3a-abba-5dc9ed8865f5
# ╟─1bbe92e6-3a06-4c0f-bc7d-e2117da4f6c1
# ╟─ac693757-bb8e-467e-b19d-c15ed7fd244a
# ╟─46620c64-a43c-4d08-bb3e-5fda3c8ff57c
# ╟─d365a1b6-ca5d-4719-99af-ce16b078b6b2
# ╟─ac5a79d7-07e0-499e-9dab-519ae0b5f638
# ╟─92f0e4d8-46b1-45ae-9c0a-1c3485713c94
# ╟─7930a403-bd4d-4f64-a79b-54cff28b869b
# ╟─06d37eb4-b3e6-4219-82bd-fc9d779885c9
# ╟─0d4425d5-964a-4ce5-aab6-0b388b523ee3
# ╠═88cf8bbf-55cd-45f0-9ab4-ae2bdf6c554a
# ╠═929010af-2af6-4b63-bfb4-7987f1b790ca
# ╟─2b04681b-6704-4172-9bf1-279051fc1d78
# ╟─0d61a532-e185-49bc-baf9-c22b8443a979
# ╠═d8091ef9-7350-45e9-8381-55c9297429cd
# ╟─ecd71729-ea89-40f6-86c8-598e42cd787d
# ╟─baabfa2d-04b8-4305-9a73-2008cd018a09
# ╟─e50fb224-63c1-4885-8165-21491110b804
# ╠═03d48dc9-a7a6-480b-a755-2c4e049d3f7b
# ╠═2e1751b5-a2cd-4a39-a7a8-4601f4448ff4
# ╟─28088d7f-cbcb-4c17-95bd-c3790da34ae2
# ╟─82286823-1240-41b5-92d9-9b3e091129e5
# ╟─e1a1e4a7-21b6-4a88-92f9-a4d3e0dacbb9
# ╟─e5725d2b-527f-4a50-9a2a-89d0514ae6bb
# ╠═b5a13fb0-6d04-411f-ab9e-5ed61f2451e5
# ╟─9be996e7-9a33-4399-9043-28d1f93cd890
# ╟─efcd0f6a-4c48-49c8-9c7a-f49421170a33
# ╟─04cbf165-4786-4980-8055-5efca14dfa3f
# ╟─a0390836-e86b-46fb-bb42-fe6e5b09bd64
# ╟─bdff9554-58b6-466e-9c93-6b1367262b50
# ╟─6c823b3d-2332-4299-8fef-bc392d13da94
# ╟─74ff4f4f-66a9-4de7-826b-4cc66d445a73
# ╟─0449d301-0896-4c2e-952c-0d380e60c0b9
# ╟─4cf5f022-4047-49f1-8fc5-86eda6d61e30
# ╟─7e3558c8-b16d-47a4-87ea-d2af5614ea2e
# ╟─fa140f9b-af22-4060-9a66-027189f80dbc
# ╟─fc9d11e1-b3b6-4d74-9cbe-1faf2d8f3abc
# ╠═4e746757-26a0-45c4-a190-04dd92b3cab0
# ╟─f1850ccb-e682-40db-91f1-ca2c4f60b3f3
# ╟─702ffb33-7fbe-4673-aed7-d985a76b455a
# ╟─6192f3a8-b9ff-4187-9561-00c2638a483a
# ╠═1307929a-f151-4868-a43a-85488599f0df
# ╟─8919edbb-5ea7-4956-ad75-3eb68a47408c
# ╠═b6510e80-6b64-421e-b8aa-bfa803d56e81
# ╟─42216744-eaed-4bb5-b882-dc4589d8e63d
# ╠═9bc67152-9220-4905-9d2b-ec3f37b43707
# ╟─b6edac36-ee05-4efa-8d75-7840e9e185e4
# ╠═95120495-15c7-4567-84f1-cb6f535c3885
# ╟─6686d3bc-abc8-4a3a-a3cd-f6cac50fdf34
# ╟─6d25666f-68c2-4a15-bf01-732fa2ebd4e3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
