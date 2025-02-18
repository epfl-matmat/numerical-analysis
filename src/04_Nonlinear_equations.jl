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

# ╔═╡ b258697b-4c4f-4e47-8472-73b72178c108
TableOfContents()

# ╔═╡ 74f7ec18-1aa1-4de7-b19f-a95ad11c8cf6
md"""
# Root finding and fixed-point problems

We saw in the introduction that problems where one wants to find the **root** or **zero** of a function $f$ arise rather naturally in scientific questions. Moreover as soon as $f$ is more complicated than just a simple polynomial, it becomes **quickly challenging to find its roots** using pen and paper.
We thus have to develop numerical methods for solving such problems.

Let us first define the problem formally:

!!! info "Definition: Rootfinding problem"
    Given a continous scalar function $f : [a, b] \to \mathbb{R}$
    find a value $x_\ast \in [a, b]$ such that
    ```math
    f(x_\ast) = 0.
    ```
    Such a value $x_\ast$ is called a **root** of $f$.
"""

# ╔═╡ 126314f4-8488-4338-a419-fc87e8965ee3
md"""
In fact an equivalent and sometimes more intuitive way to think about solving equations $f(x) = 0$ is to recast them as a fixed point-problem:

!!! info "Definition: Fixed-point problem"
    Given a continous function $g : [a, b] \to \mathbb{R}$ a point
    $x_\ast$, such that $g(x_\ast) = x_\ast$, is called a **fixed point** of $g$.
    The problem of finding a fixed point of $g$ is equivalent to finding a root
    of the function $f(x) = g(x) - x$.

Such transformations can also be reversed, i.e. given an $f$
for root finding, we could also seek a fixed-point of $g(x) = x + f(x)$.
Moreover the transformations are not unique.
Indeed, for any $c\neq 0$ finding a fixed-point $x_\ast$ of $g(x) = x + c f(x)$
implies $f(x_\ast) = 0$.

To start of our investigation we will explore different ways of numerically solving our circuit problem.
"""

# ╔═╡ fe49686f-d5b2-4601-88b5-88c3d8d5fe5b
md"""
## Revisiting the diode model

Recall the circuit diagram of the diode model
"""

# ╔═╡ 4170d04a-063f-4f79-bd57-67e5cfd2af63
RobustLocalResource("https://raw.githubusercontent.com/epfl-matmat/numerical-analysis/09d365542432dc652a9ed1ce5c5f54075590aebd/notes/img/circuit.png", "img/circuit.png")

# ╔═╡ 941fa7db-20ee-4900-a193-80ac0b45df5b
begin
	i0 = 1.0
	v0 = 0.1
	R = 1.0
	V = 1.0
end;

# ╔═╡ 6bc129ec-2485-480e-b8c8-3e167562581e
md"""
where we had the relationships
```math
\tag{1}
\begin{aligned}
V &= v_D + v_R \\
v_R &= R\, i \\
i &= i_0 (e^{v_D / v_0} - 1) \qquad \text{(Shockley diode model)}
\end{aligned}
```
and wanted to solve for the diode voltage $v_D$. We consider two iterative methods to find it, motivated from the physical setting.

**Method 1**

- In many physical settings one often already has a good general idea about what $v_D$ should be like. In this case a reasonable starting point would be to assume that the diode is open and thus the voltage across the diode is zero, thus we set $v_D^{(0)} = 0$.
- Given the voltage $v_D^{(0)} = 0$ across the diode, we find the voltage across the resistor as $v_R^{(0)} = V - v_D^{(0)}$ and thus the current
  ```math
   i^{(0)} = v_R^{(0)} / R = \frac{V - v_D^{(0)}}{R}
  ```
- Since the current in all elements of the circuit is the same,
  we can estimate the voltage across the diode by reversing the Shockley diode model:
  ```math
  v^{(1)}_D = v_0 \log\left(\frac{V - v_D}{R\,i_0} + 1\right)
  ```
  We hope for this new estimation $v^{(1)}_D$ to be better than the
  initial one $v^{(0)}_D$.
- Finally we close the loop and repeat the process from top to bottom using $v^{(1)}_D$ instead of $v^{(0)}_D$, thus obtaining $v^{(2)}_D$, etc.
  The resulting sequence $v^{(0)}_D, v^{(1)}_D, v^{(2)}_D$ we hope to converge
  to the true diode voltage, that is
  ```math
  \lim_{k\to\infty} v^{(k)}_D = v_D  \qquad \text{(hopefully)}
  ```

Mathematically, the voltage estimate obtained in the $k$-th iteration is given by
```math
v^{(k+1)}_D = g_\text{log}(v^{(k)}_D) \qquad \text{for $k = 1, 2, \ldots$},
```
where 
```math
\tag{2}
g_\text{log}(x) = v_0 \log\left( \frac{V - x}{R \, i_0} + 1 \right).
```
If the method converges to a value $v_D$, then $v_D = g_\text{log}(v_D)$, i.e. we found a fixed point of $g_\text{log}$.

Let's see if this actually works for our example. We select the parameters:
  * `i0 =` $(i0)
  * `v0 =` $(v0)
  * `R  =` $(R)
  * `V  =` $(V)

and define the fixed-point map
"""

# ╔═╡ 99be8b53-7c56-4a9b-8fb2-cb5be7a6c8a2
function glog(vD)
	v0 * log((V - vD) / (R * i0) + 1)
end

# ╔═╡ 657e28a9-0e6c-4b99-9122-50cae4b5670a
md"""
finally we iterate 10 steps, recording along the way the produced iterates
"""

# ╔═╡ 4b337c8b-d7d5-4bfc-b8a6-b07075da6ade
let
	vD = 0.1
	for i in 1:10
		vD = glog(vD) # Call fixed-point map
		@printf "iter %2i:  vD = %.15f\n" i vD  # Print formatted data
	end
end

# ╔═╡ 7b834e4b-68d5-432d-a7b8-c38fa8303c18
md"""
Hooray! After 10 steps about 12 significant figures have stabilised and we thus expect the error to be smaller than $10^{-12}$.

But setting up such iterations is in fact not unique:
"""

# ╔═╡ bd73ab5b-ccd5-4e47-b5ad-2459ced9e914
md"""
**Method 2**

- Assume again a $v_D^{(0)}$ is provided by physical intuition.
- Using the Shockley relation we obtain the current through the diode as $i^{(0)} = i_0 \left( e^{v_D^{(0)}/v_0} - 1 \right)$, which again equals the current through the resistor.
- The voltage across the resistor is thus $v_R^{(0)} = R\, i^{(0)}$, from which
  we can compute an updated voltage estimation across the diode as $v^{(1)}_D = V - v_R^{(0)}$

Again closing the loop and iterating multiple times, we obtain at the $k$-th iteration
```math
v^{(k+1)}_D = g_\text{exp}(v^{(k)}_D) \qquad \text{for $k = 1, 2, \ldots$},
```
where
```math
\tag{3}
g_\text{exp}(x) = V - R \, i_0 \left( e^{x/v_0} - 1 \right).
```

Let's test this method as well:
"""

# ╔═╡ 1f0c7a58-174d-46a3-9dad-b4ee3536f6c1
function gexp(vD)
	V - R * i0 * (exp(vD/v0)-1) 
end

# ╔═╡ b016f09e-972d-4d22-bee3-546febaf342c
let
	vD = 0.1
	for i in 1:10
		vD = gexp(vD) # Call fixed-point map
		@printf "iter %2i:  vD = %.15f\n" i vD  # Print formatted data
	end
end

# ╔═╡ 1b6aa472-a5f3-4953-89d0-a7c579ce33d6
md"""
Even though the method seems equally plausible as method 1 on first sight, it clearly performs worse and does not converge at all. For finding the voltage $v_D$ across the diode it is thus not useful.

Let us now formalise these methods mathematically in order to analyse them more carefully.
"""

# ╔═╡ d3ba1da7-41eb-4afd-a949-19123b7a6724
md"""
## Fixed-point iterations

As discussed in the previously the equations of the diode model (1) lead to the non-linear problem
```math
\tag{4}
f(v_D) = R \, i_0 \left( e^{v_D/v_0} - 1 \right) + v_D - V,
```
where its root is the desired diode voltage. In both **Method 1** and **Method 2** we effectively rewrote this equation into a different fixed-point problem:
```math
\begin{aligned}
g_\text{log}(v_D) &= v_D & g_\text{log}(x) &= v_0 \log\left( \frac{V - x}{R \, i_0} + 1 \right)\\
g_\text{exp}(v_D) &= v_D & g_\text{exp}(x) &= V - R \, i_0 \left( e^{x/v_0} - 1 \right)
\end{aligned}
```
To solve these fixed-point problems we then applied

!!! info "Algorihm: Fixed-point iteration"
    Given a fixed-point map $g$ and initial guess $x^{(0)}$, iterate
    ```math
    x^{(k+1)} = g(x^{(k)})  \qquad{\text{for $k = 1, 2, \ldots$}}.
    ```
    until a stopping criterion is reached.
"""

# ╔═╡ 673ef83d-0242-4e24-a21f-85a56c2b3ba1
md"""
If $\lim_{k\to\infty} x^{(k)} = x_\ast$ for the sequence $x^{(k)}$ generated by Algorithm 1, then $g(x_\ast) = x_\ast$, i.e. $x_\ast$ is a fixed point of $g$.
At this point we do not yet specify what is a good stopping criterion. We will return to this point in the *Convergence analysis* section.

If we are faced with a fixed point problem (finding $x_\ast$ such that $g(x_\ast) = x_\ast$) then Algorithm 1 can be directly applied. However, to **apply** the **fixed-point method** to a **root-finding problem** (seek $x_\ast$ s.t. $f(x_\ast) = 0$)
we first need to **rewrite the non-linear equation** $f(x) = 0$ into a fixed-point problem, i.e. identify a suitable $g$, such that if $x_\ast$ is a root of $f$, then
```math
f(x_\ast) = 0 \qquad \Longleftrightarrow \qquad x_\ast = g(x_\ast).
```
On $g$ we then apply fixed-point iteration.

We will see an example of this rewriting in the next section.
"""

# ╔═╡ d54f9a44-d41d-4d98-85ad-821a2a06dadc
md"""
### Visualising fixed-point iterations

Before we proceed to a closer mathematical analysis of fixed-point methods and ultimately to understand the difference between employing $g_\text{log}$ and $g_\text{exp}$, let us first consider a simpler case, where it is easier to obtain a visual understanding.

We consider solving non-linear equation $f(x) = 0$ with
```math
f(x) = x + \log(x + 1) - 2
```
"""

# ╔═╡ 1319313d-e7da-4ff5-9725-548df0e7f699
md"""
In this one can easily construct four equivalent fixed-point equations $x = g_i(x)$ with $i=1,2,3,4$, namely
"""

# ╔═╡ 858510eb-7768-4521-ba20-af5b1aa33e08
begin
	g₁(x) = x - (1/2) * (log(x + 1) + x - 2)
	g₂(x) = 2 - log(x + 1)
	g₃(x) = exp(2-x) - 1
	g₄(x) = (1/2) * x * (log(x+1) + x)
end;

# ╔═╡ 4715b982-6ba9-415d-b33b-6a2f145904d9
md"""
!!! warning "Example: Rewriting f(x)=0 as a fixed point problem"
    We are given the problem to solve $f(x) = 0$ with $f(x) = x + \log(x + 1) - 2$.
    We show how to construct $g_3$ and $g_4$.

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

# ╔═╡ e4e37eb1-9013-4dc2-887f-fcafb51aac41
md"""
To visualise the fixed-point iterations, we choose yet another understanding
of fixed point problems:

!!! info "Observation: Fixed point problems are about curve intersections"
    We can understand a fixed-point problem, where we seek
    an $x_\ast$ such that $g(x_\ast) = x_\ast$ as the problem of finding the
    **intersection of the curve** $y = g(x)$ with the **line** $y = x$.
    In other words we want to find a point $(x, y)$ such that
    $y = g(x)$ and at the same time $x = y$.

The following plot uses this observation to visualise fixed-point iterations. Use the slider and the drop-down menu to switch between the fixed-point functions and the starting point.
"""

# ╔═╡ 424bb138-36b5-46a2-8af3-1eb8afc1c2cb
md"""
- `xstart = ` $(@bind xstart Slider(0.0:0.1:5.0; default=4.0, show_value=true))
- `g = ` $(@bind g Select([g₁, g₂, g₃, g₄]))
- Show gradient: $(@bind show_gradient CheckBox())
"""

# ╔═╡ ed74f330-28f0-4550-bc3d-dbf73f220bf6
let
	# Plot function g and identity
	p = plot(g, label=L"y = g(x)";
	         aspect_ratio=1, lw=2, xlims=(0, 5),
	         ylims=(0, 5), xlabel=L"x", ylabel=L"y", color=1)
	plot!(p, x -> x, label=L"y = x", lw=2, color=2)

	# Iteration 1
	# We evaluate g at x, which gives us the initial point (x, y)
	x = xstart
	vline!(p, [xstart], label="xstart", color=:grey, lw=2, ls=:dash)
	y = g(x)
	scatter!(p, [xstart], [y], label="", color=:grey, mark=:o, ms=5)
	
	# Iteration 2
	# Using the requirement x = y, we will next update x to be y.
	# This moves the point *horizontally*
	plot!(p, [x, y], [y, y]; arrow=true, color=3, label="", lw=2)
	x = y
	# Next we re-evaluate y = g(x) at the *new* x, this updates y
	# therefore we move *vertically*
	y = g(x)
	plot!(p, [x, x], [x, y]; arrow=true, color=4, label="", lw=2)

	# Now just repeat until 5 iterations are reached.
	for k = 3:5
		plot!([x, y], [y, y]; color=3, label="", lw=2)
		x = y
		y = g(x)
		plot!(p, [x, x], [x, y]; color=4, label="", lw=2)
	end

	# Plot guiding lines through the fixed point
	# (useful for discussion below)
	fp = 1.2079400315693  # Approximate fixed point
	if show_gradient
		# compute gradient at the fixed point
		grad = ForwardDiff.derivative(g, fp)

		plot!(p, [0, fp, 0], [0, fp, 2fp], fill=(0, :orange), α=0.3, label="", lw=0)
		plot!(p, [5, fp, 5], [2fp-5, fp, 5], fill=(0, :orange), α=0.3, label="",lw=0)
		plot!(p, x -> fp - (x-fp), label=L"slope $-1$", lw=2, color=2, ls=:dash)
		plot!(p, x -> fp + grad * (x-fp), label=L"gradient of $g$", lw=2, color=1, ls=:dash)
	end
	
	p
end

# ╔═╡ b4100dc2-8e77-4455-a5d8-0714145b1979
md"""
We notice that that $g_1$ and $g_2$ converge, while $g_3$ and $g_4$ do not.

Their distinctive feature becomes clear when we also show the gradient of $g$ at the fixed point (click the checkbox above).

Convergen occurrs when the slope is between $-1$ and $1$, i.e. the graph sits between the regions spanned by the orange lines of slope $1$ and $-1$ through the fixed point.
"""

# ╔═╡ 49c6c848-dbc1-4da7-8c0a-fe90ab07ef93
md"""
### Convergence analysis

We consider the setting finding a fixed-point $x_\ast = g(x_\ast)$ for a differentiable function $g$. Our hypothesis is that the fixed-point method converges if $g'(x_\ast) \in (-1, 1)$, while it diverges when $|g'(x_\ast)| > 1$.

To establish this result more rigorously we study the beviour of the error sequence
```math
e^{(k)} = x^{(k)} - x_\ast
```
The goal is to find an expression that relates $e^{(k+1)}$ with $e^{(k)}$ and in this way *recursively* to $e^{(0)}$, the error of the initial guess.

First we note $x^{(k)} = x_\ast + e^{(k)}$ and we develop the 
Taylor expansion of $g$ around $x_\ast$:
```math
\tag{5}
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

Using (5) and the key fixed-point iterations equation, $x^{(k+1)} = g(x^{(k)})$
we obtain
```math
\begin{aligned}
e^{(k+1)}
&= x^{(k+1)} - x_\ast = g(x^{(k)}) - x_\ast \\
&\stackrel{(5)}{=} x_\ast + g'(x_\ast) \underbrace{(x^{(k)}-x_\ast)}_{=e^{(k)}} + O(|e^{(k)}|^2) - x_\ast \\
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
if $|g'(x_\ast)| < 1$, exactly as we concluded from the plot.
Our argument proves the following

!!! note "Theorem 1"
    Let $g : [a, b] \to \mathbb{R}$ be a function of class $C^1$[^1]
    and $x_\ast \in [a, b]$ be a fixed point of $g$.
    If $|g'(x_\ast)| < 1$, then
    there exists an $ε > 0$,
    such that for all $x^{(0)} \in [x_\ast - ε, x_\ast + ε]$
    - the fixed-point iterations
      $x^{(k+1)} = g(x^{(k)})$ converge to $x_\ast$,
      i.e. $\lim_{k\to\infty} x^{(k)} = x_\ast$
    - Moreover the **convergence rate** *(formal definition below)* is given by
      ```math
      \lim_{k\to\infty} \frac{|x^{(k+1)} - x_\ast|}{|x^{(k)} - x_\ast|} = |g'(x_\ast)|,
      ```
       i.e. the smaller the gradient, the faster the convergence.

[^1]: differentiable with continuous first derivative.
"""

# ╔═╡ 5949bb96-3f2c-4d31-a7df-b5cc681914e1
md"""
!!! example "Exercise"
    Verify the theorem for the fixed-point problems we considered so far,
    i.e. show that
    - The map $g_\text{log}$ has a gradient modulus less than $1$ at the fixed point while $g_\text{exp}$ has one larger $1$.
    - similarly verify analytically the convergence of the fixed-point iterations
       of $g_1$ and $g_2$
       and the divergence of $g_3$ and $g_4$.
"""

# ╔═╡ 5d7b3d35-3456-48df-ad22-0ffccaa2f529
md"""
### Stopping criteria and residual

Let us come back to the question at which point to stop the fixed-point iteration algorithm. Let $\epsilon$ denote the tolerance to which we want to determine $x_\ast$, i.e. we would like to stop the iterations as soon as the error is smaller, $|x^{(k)} - x_\ast| < \epsilon$.

Since $x_\ast$ is not known, this expression cannot be exactly computed during the iterations. We thus need to seek an alternative approach. In a given step $x^{(k)}$ during the iteration we have likely not yet achieved our goal,
i.e. $g(x^{(k)}) \neq x^{(k)}$.
A natural idea is thus to consider exactly the descrepancy 
```math
r^{(k)} = g(x^{(k)}) - x^{(k)},
```
the so-called **residual**. A natural stopping criterion is thus
"""

# ╔═╡ 284f3fa4-ce24-4b99-8bb7-6f74a4589550
md"""
!!! note "Algorithm: Fixed-point iteration stopping criterion"
    ```math
    |r^{(k)}| = |g(x^{(k)}) - x^{(k)}| < \epsilon.
    ```

Employing this stopping criteria, the algorithm to find a fixed point of the function $g$ becomes
"""

# ╔═╡ 899698d1-9dee-4f82-9171-f1b49aefcabe
function fixed_point_iterations_simple(g, xstart; tol=1e-6)
	# g:      Fixed-point function
	# xstart: Initial guess
	# tol:    Tolerance

	r⁽ᵏ⁾ = Inf
	x⁽ᵏ⁾ = xstart
	k  = 0
	while abs(r⁽ᵏ⁾) ≥ tol
		x⁽ᵏ⁺¹⁾ = g(x⁽ᵏ⁾)
		r⁽ᵏ⁾   = x⁽ᵏ⁺¹⁾ - x⁽ᵏ⁾

		k    = k + 1  # Update k
		x⁽ᵏ⁾ = x⁽ᵏ⁺¹⁾ # Update x⁽ᵏ⁾ accordingly
	end

	# Return results as a named tuple
	(; fixed_point=x⁽ᵏ⁾, residual=r⁽ᵏ⁾, n_iter=k)
end

# ╔═╡ 072cdb08-b279-4583-b98b-4a28433b5b8d
md"""We apply this function:"""

# ╔═╡ 3c0d4752-de1a-4819-a63f-aa49c94ccf82
res_simple = fixed_point_iterations_simple(g₁, 4.0; tol=1e-14)

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

# ╔═╡ 54f314c1-d1cf-41f1-96e5-5aca90d82b95
fixed_point_iterations(g₁, 4.0; tol=1e-14).fixed_point

# ╔═╡ 0adb95c9-d86f-4ad0-b4c6-ebe616587168
md"""
**Additional remarks on the residual**
- The residual is in general only an **error indicator**.
  This means that **there is no guarantee**
  that $|g(x^{(k)}) - x^{(k)}| < \epsilon$
  **always implies** $|x^{(k)} - x_\ast| < \epsilon$.

- In fact using the [mean value theorem](https://en.wikipedia.org/wiki/Mean_value_theorem) (MVT) we can show
  ```math
  \tag{6}
  \begin{aligned}
  x^{(k)} - x_\ast
  &= x^{(k)} - g(x^{(k)}) + g(x^{(k)}) - x_\ast \\
  &= x^{(k)} - g(x^{(k)}) + g(x^{(k)}) - g(x_\ast) \\
  &\stackrel{(MVT)}{=} x^{(k)} - g(x^{(k)}) + g'(\xi^{(k)}) (x^{(k)} - x_\ast)
  \end{aligned}
  ```
  where $\xi^{(k)} \in [x_\ast, x^{(k)}]$,
  thus we obtain the **residual-error relationship**
  ```math
  |x^{(k)} - x_\ast| = \frac{1}{|1 - g'(\xi^{(k)})|} r^{(k)}.
  ```
  In particular if $g'(\xi^{(k)}) \simeq 1$ the **true error
  can thus be considerably larger than the residual** !

- Since a convering fixed-point method implies that
  $g'(\xi^{(k)}) \to g'(x_\ast)$ as $k \to \infty$,
  we note that it is again the gradient at the fixed point,
  $g'(x_\ast)$, which determines how reliable our error indicator is.
  Notably, if $g'(x_\ast) \simeq 1$ we might stop the iterations
  *well before* achieving $|x^{(k)} - x_\ast| < \epsilon$.

- Note, that the **residual is a general terminology**,
  which is not only applied to such an error indicator in the context
  of fixed-point iterations, but used in general for iterative procedures.
  The idea is that the **residual provides
  the discrepancy from having fully solved the problem**
  and is thus a natural error indicator.
  The functional form, however, is different for each type of iterative
  procedure as we will see.


"""

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
    If $q = 1$ (convergence order 1) we additionally require $C < 1$.
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
	p = plot(; yaxis=:log)

	results_g₁ = fixed_point_iterations(g₁, 4.0; maxiter=15)
	errors_g₁ = [abs(xn - fp) for xn in results_g₁.history_x]
	plot!(p, errors_g₁, label=L"g_1", mark=:x, lw=2)
	
	results_g₂ = fixed_point_iterations(g₂, 4.0, maxiter=15)
	errors_g₂ = [abs(xn - fp) for xn in results_g₂.history_x]
	plot!(p, errors_g₂, label=L"g_2", mark=:x, lw=2)
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
\stackrel{(7)}{=} \frac{1}{|1 - g'(x_\ast)|^{q-1}} C
```
In other words if the residual ratios approach a constant for a chosen
$q$, then we have $q$-th order convergence.
In particular for **linear order** ($q=1$) we need
$\frac{|r^{(k+1)}|}{\big| r^{(k)} \big|}$ to approach a constant as $k\to\infty$.

This is a condition we can check also
*without* knowing the solution:
"""

# ╔═╡ d22a7de4-67ef-4664-951b-8fb46116a7bc
let
	p = plot(xlabel="k")

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

For this let us revisit our diode circuit problem. the non-linear equation
to solve was
```math
f(x) =  R \, i_0 \left( e^{x/v_0} - 1 \right) + x - V.
```
We will now directly attempt to find a root $f(x) = 0$ of this equation.
"""

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
	k = 0
	x⁽ᵏ⁾ = (a + b) / 2
	
	history_x = Float64[]  # Empty Array, but only for Float64 numbers
	while abs(b - a) / 2 ≥ tol
		k = k + 1
		if f(x⁽ᵏ⁾) * f(a) < 0
			b = x⁽ᵏ⁾  # New interval [a, x⁽ᵏ⁾]
		else
			a = x⁽ᵏ⁾  # New interval [x⁽ᵏ⁾, b]
		end

		x⁽ᵏ⁾ = (a + b) / 2
		push!(history_x, x⁽ᵏ⁾)
	end

	(; root=x⁽ᵏ⁾, history_x, n_iter=k)
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
the **Taylor expansion** (5) of the **fixed-point map** $g$.
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
0 = f(x_\ast) = f(x) + f'(x) (x - x_\ast) + O\big( (x - x_\ast)^2 \big)
```
Where we made the condition $f(x_\ast) = 0$ has been made explicict.
Now assume $f'(x) \neq 0$ to rearrange this to
```math
0 = \frac{f(x)}{f'(x)} (x - x_\ast) + O\big( (x - x_\ast)^2 \big).
```
If $x$ is close to $x_\ast$, thus $x - x_\ast$ is small,
then the last $O \big((x - x_\ast)^2\big)$ term is even smaller.
Neglecting it we can further develop this to an approximation for $x_\ast$ as
```math
x_\ast \simeq x - \frac{f(x)}{f'(x)}
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
    Let $f$ be a $C^2$ function and $x_\ast$ a root of $f$.
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
md"""
More conventionally one "inlines" the function $g_\text{Newton}$ into the fixed point iterations and expresses the problem as

!!! info "Algorithm 2: Newton's method (conventional)"
	Given a once differentiable function $f : \mathbb{R} \to \mathbb{R}$,
	a starting value $x^{(0)}$ and a convergence tolerance $\epsilon$,
    perform for $k = 1, 2, 3, \ldots$:
    1. Compute the residual $r^{(k)} = - \frac{f(x^{(k)})}{f'(x^{(k)})}$
    2. Update ${x}^{(k+1)} = {x}^{(k)} + {r}^{(k)}$
    Loop 1. to 3. until $|{r}^{(k)}| < \epsilon$.

A corresponding implementation of Algorithm 2 is:
"""

# ╔═╡ d8091ef9-7350-45e9-8381-55c9297429cd
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

# ╔═╡ 2e1751b5-a2cd-4a39-a7a8-4601f4448ff4
let
	# First get an almost exact root
	reference = bisection_method(f, -0.5, 0.5; tol=1e-14)

	p = plot(yaxis=:log, xlims=(0, 20), ylims=(1e-12, Inf))

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
for example in [Iterative methods for linear systems](https://teaching.matmat.org/numerical-analysis/07_Iterative_methods.html).
"""

# ╔═╡ bdff9554-58b6-466e-9c93-6b1367262b50
md"""
## Optional: Secant method

See [chapter 4.4](https://tobydriscoll.net/fnc-julia/nonlineqn/secant.html)
of Driscoll, Brown: *Fundamentals of Numerical Computation*.

"""

# ╔═╡ 9eb4a16f-899c-4b9b-be88-ff11f6c293f8
md"""
## Non-linear equation systems

Many applications are characterised by by more than one degree of freedom
and more than a single equation to satisfy. Here, we will generalise our discussion
and consider a system of equations

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
we can define the multi-dimensional version of the root-finding problem:

!!! info "Definition: Multidimensional root-finding problem"
    Given a continuous vector-valued function
    $\textbf{f} : \mathbb{R}^n \to \mathbb{R}^n$
    find a vector $\textbf{x}_\ast \in \mathbb{R}^n$
    such that $\textbf{f}(\textbf{x}_\ast) = \textbf{0}$.

Solving such nonlinear multi-dimensional equation systems
is much more involved. Even establishing basic mathematical
properties, such as the existance or uniqueness of solutions
is typically quite difficult, let alone solving such equation
systems analytically.

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

	r = Inf  # Dummy to enter the while loop
	k = 0
	while norm(r) ≥ tol && k < maxiter
		k = k + 1
		
		x = last(history_x)
		y = f(x)      # Function value
		A = jac(x)    # Jacobian
		r = -(A \ y)  # Newton step
		
		push!(history_r, r)      # Push newton step and
		push!(history_x, x + r)  # next iterate to history
	end

	(; root=last(history_x), n_iter=k, history_x, history_r)
end

# ╔═╡ f1850ccb-e682-40db-91f1-ca2c4f60b3f3
md"""
Note that he linear system $\textbf{A}^{(k)} \textbf{r}^{(k)} = - \textbf{y}^{(k)}$ is solved in Julia using the backslash operator `\`, which employs a numerically more stable algorithm than explicitly computing the inverse `inv(A)` and then applying this to `y`.
We will discuss these methods in
[Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.html).
"""

# ╔═╡ 702ffb33-7fbe-4673-aed7-d985a76b455a
Foldable("Remark: Connection to 1D Newton algorithm (Algorithm 1)",
md"""
Even though Algorithm 2 (Multidimensional Newton) and Algorithm 1 (1D Newton's method) are presented in different formalisms, Algoritm 2 is just a multi-dimensional generalisation of Algorithm 1.

This can be easily seen by simplifying Algorithm 2 for the special case of a 1D problem, that is a function $f : \mathbb{R} \to \mathbb{R} : x_1 \mapsto f(x_1) $. In that case we can identify the Jacobian
(which formally is a $1 \times 1$ matrix) by the single element
$\frac{\partial f}{\partial x_1} = f'(x_1)$. With this in mind the three steps of Algorithm 2 become: For $k = 1, 2, 3, \ldots$
1. Compute the right-hand side $y^{(k)} = f(x^{(k)})$ and $A^{(k)} = J_f(x^{(k)}) = f'(x^{(k)})$.
2. Solve $A^{(k)} r^{(k)} = -y^{(k)}$ for $r^{(k)}$. Since all entities are real numbers this can be done by simple division, i.e. $r^{(k)} = - \frac{y^{(k)}}{A^{(k)}} = -\frac{f(x^{(k)})}{f'(x^{(k)})}$.
3. Update $x^{(k+1)} = x^{(k)} + r^{(k)} = x^{(k)} - \frac{f(x^{(k)})}{f'(x^{(k)})}$.

Comparing with Algorithm 1 we note that in step 3 we indeed recover $x^{(k+1)} = g_\text{Newton}(x^{(k)})$ with $g_\text{Newton}(x) = x - \frac{f(x)}{f'(x)}$ as before in (10).
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
ForwardDiff = "~0.10.36"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.9"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "48728d14071e58de9f09326dc3a99600cb737f60"

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
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "545a177179195e442472a1c4dc86982aa7a1bef0"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.7"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "26ec26c98ae1453c692efded2b17e15125a5bea1"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.28.0"

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

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

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
git-tree-sha1 = "f36e5e8fdffcb5646ea5da81495a5a7566005127"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.3"

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
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

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
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

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
git-tree-sha1 = "21fac3c77d7b5a9fc03b0ec503aa1a6392c34d2b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.15.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "786e968a8d2fb167f2e4880baba62e0e26bd8e4e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "846f7026a9decf3679419122b49f8a1fdb48d2d5"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.16+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "9bf00ba4c45867c86251a7fd4cb646dcbeb41bf0"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.12"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "36d5430819123553bf31dfdceb3653ca7d9e62d7"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.12+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "01979f9b37367603e2848ea225918a3b3861b606"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "c67b33b085f6e2faf8bf79a61962e7339a81129c"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.15"

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
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "71b48d857e86bf7a1838c4736545699974ce79a2"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.9"

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
git-tree-sha1 = "a729439c18f7112cbbd9fcdc1771ecc7f071df6a"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.39"

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
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

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

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "df37206100d39f79b3376afb6b9cee4970041c61"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.51.1+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "89211ea35d9df5831fca5d33552c02bd33878419"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e888ad02ce716b319e6bdb985d2ef300e7089889"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.3+0"

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
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

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
git-tree-sha1 = "cc0a5deefdb12ab3a096f00a6d42133af4560d71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.2"

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
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+3"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

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
git-tree-sha1 = "ed6834e95bd326c52d5675b4181386dfbe885afb"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.55.5+0"

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
git-tree-sha1 = "dae01f8c2e069a683d3a6e17bbae5070ab94786f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.9"

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
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

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
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "9bb80533cb9769933954ea4ffbecb3025a783198"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.7.2"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

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
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

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
git-tree-sha1 = "a2fccc6559132927d4c5dc183e3e01048c6dcbd6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56c6604ec8b2d82cc4cfe01aa03b00426aac7e1f"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.4+1"

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
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+3"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e9216fdcd8514b7072b43653874fd688e4c6c003"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.12+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "807c226eaf3651e7b2c468f687ac788291f9a89b"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.3+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "89799ae67c17caa5b3b5a19b8469eeee474377db"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.5+0"

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

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c57201109a9e4c0585b208bb408bc41d205ac4e9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.2+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "1a74296303b6524a0472a8cb12d3d87a78eb3612"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+3"

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
git-tree-sha1 = "6dba04dbfb72ae3ebe5418ba33d087ba8aa8cb00"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.1+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "622cf78670d067c738667aaa96c553430b65e269"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+0"

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
git-tree-sha1 = "055a96774f383318750a1a5e10fd4151f04c29c5"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.46+0"

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
# ╠═62ff56f9-caec-4352-b547-973ce2bfcb8f
# ╟─b258697b-4c4f-4e47-8472-73b72178c108
# ╟─74f7ec18-1aa1-4de7-b19f-a95ad11c8cf6
# ╟─126314f4-8488-4338-a419-fc87e8965ee3
# ╟─fe49686f-d5b2-4601-88b5-88c3d8d5fe5b
# ╟─4170d04a-063f-4f79-bd57-67e5cfd2af63
# ╟─6bc129ec-2485-480e-b8c8-3e167562581e
# ╟─941fa7db-20ee-4900-a193-80ac0b45df5b
# ╠═99be8b53-7c56-4a9b-8fb2-cb5be7a6c8a2
# ╟─657e28a9-0e6c-4b99-9122-50cae4b5670a
# ╠═4b337c8b-d7d5-4bfc-b8a6-b07075da6ade
# ╟─7b834e4b-68d5-432d-a7b8-c38fa8303c18
# ╟─bd73ab5b-ccd5-4e47-b5ad-2459ced9e914
# ╠═1f0c7a58-174d-46a3-9dad-b4ee3536f6c1
# ╠═b016f09e-972d-4d22-bee3-546febaf342c
# ╟─1b6aa472-a5f3-4953-89d0-a7c579ce33d6
# ╟─d3ba1da7-41eb-4afd-a949-19123b7a6724
# ╟─673ef83d-0242-4e24-a21f-85a56c2b3ba1
# ╟─d54f9a44-d41d-4d98-85ad-821a2a06dadc
# ╟─1319313d-e7da-4ff5-9725-548df0e7f699
# ╠═858510eb-7768-4521-ba20-af5b1aa33e08
# ╟─4715b982-6ba9-415d-b33b-6a2f145904d9
# ╟─e4e37eb1-9013-4dc2-887f-fcafb51aac41
# ╟─424bb138-36b5-46a2-8af3-1eb8afc1c2cb
# ╠═ed74f330-28f0-4550-bc3d-dbf73f220bf6
# ╟─b4100dc2-8e77-4455-a5d8-0714145b1979
# ╟─49c6c848-dbc1-4da7-8c0a-fe90ab07ef93
# ╟─5949bb96-3f2c-4d31-a7df-b5cc681914e1
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
# ╟─0adb95c9-d86f-4ad0-b4c6-ebe616587168
# ╟─91060f66-9540-43dd-bffa-3ecc022efbdc
# ╟─deff7837-e80e-4fcc-8eb2-baafc047a992
# ╟─7f2c43c0-a87e-44bb-ac52-77dc70302cd2
# ╠═b2e96aa1-69a6-4d07-89b5-15ecc9a41398
# ╟─86f4c162-41bd-4e20-9e3f-681ff36bac75
# ╟─dc136c3d-4534-4d66-86ca-36114bb825bb
# ╠═d22a7de4-67ef-4664-951b-8fb46116a7bc
# ╟─314f733c-c3c2-4679-bb7b-c94b96b54961
# ╟─68f63b13-1326-4cb6-8db1-3b043b2ad78e
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
# ╠═929010af-2af6-4b63-bfb4-7987f1b790ca
# ╟─2b04681b-6704-4172-9bf1-279051fc1d78
# ╟─0d61a532-e185-49bc-baf9-c22b8443a979
# ╠═d8091ef9-7350-45e9-8381-55c9297429cd
# ╟─ecd71729-ea89-40f6-86c8-598e42cd787d
# ╟─baabfa2d-04b8-4305-9a73-2008cd018a09
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
# ╟─9eb4a16f-899c-4b9b-be88-ff11f6c293f8
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
