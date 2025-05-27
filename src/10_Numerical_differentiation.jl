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

# ╔═╡ 9688d440-d231-11ee-0510-8db284ccd0c9
begin
	using Plots
	using PlutoUI
	using PlutoTeachingTools
	using Printf
	using LaTeXStrings
	using HypertextLiteral
	using Symbolics
	using ForwardDiff
end

# ╔═╡ 4103c9d2-ef89-4c65-be3f-3dab59d1cc47
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/10_Numerical_differentiation.pdf)
"""

# ╔═╡ e9151d3f-8d28-4e9b-add8-43c713f6f068
TableOfContents()

# ╔═╡ bc9f0276-6b7b-4b19-8ef0-6371da964a9c
md"""
# Numerical differentiation

Taking derivatives by hand can be a error-prone and time-consuming task.
Recall the innocent-looking function
"""

# ╔═╡ 7b74f71f-a5fa-47ba-b218-c1798fa540f7
v(t) = 64t * (1 − t) * (1 − 2t)^2 * (1 − 8t + 8t^2)^2

# ╔═╡ 80df99a7-3e31-4418-9314-399116831e72
let
	p = plot(v;
			 ylims=(-10, 2),
			 xlims=(-0.5, 1.5),
			 linewidth=2,
		     label="",
			 xlabel="t",
		     ylabel="v")
end

# ╔═╡ 87753dbe-f075-435b-9ac5-75bdcd6019c0
md"""
we already considered in the [Introduction](https://teaching.matmat.org/numerical-analysis/01_Introduction.html). Recall that it's derivative was a rather lengthy and technical expression, which one would like to avoid computing by hand:
"""

# ╔═╡ 0f0c0f37-f751-4950-9a8e-8dadb038325e
let
	@variables t          # Define a variable
	dt = Differential(t)  # Differentiating wrt. t
	dv_dt = dt( v(t) )    # Compute dv / dt
	dv_dt = expand_derivatives(dv_dt)
end

# ╔═╡ 657a065f-738f-4512-9a48-84fc9dae76da
md"""
In this chapter we will consider numerical techniques for computing such derivatives.
"""

# ╔═╡ 50cb2913-618f-4fd8-93a2-dd4f4411bbe1
md"""
## First order derivatives

Given a regular function $f : [a, b] \to \mathbb{R}$ our goal is to **approximate numerically its derivative** $f'$ in a point $x \in [a, b]$. We will allow ourselves **only to perform $n$ pointwise evaluations** $f(t_i)$ with $t_i \in [a, b]$ and $i = 1, \ldots, n$ and **take linear combinations of these results**. That is we work towards approximations of the form
```math
    f'(x) ≈ \sum_{i=1}^m α_i f(x_i),
```
where $x_i$ are points around $x$ and $α_i$ some coefficients,
with details how to choose these points and coefficients to be specified.
"""

# ╔═╡ 1df757af-8d84-429b-8005-b17e9814e48c
md"""
A first idea to achieve such a **numerical differentiation formula** takes us back to the definition of the derivative $f'(x)$ as the limit $h\to 0$ of the slope of secants over an intervall $[x, x+h]$, i.e.
```math
f'(x) = \lim_{h\to 0} \frac{f(x+h) - f(x)}{h}.
```
A natural idea is to not fully take the limit, i.e. to take a small $h > 0$
and approximate
```math
f'(x) ≈ \frac{f(x+h) - f(x)}{h}.
```
This immediately leads to the **forward finite difference formula**
```math
\tag{1}
D^+_h f(x) = \frac{1}{h} \Big( f(x+h) - f(x) \Big).
```
Note that in this expression
we interpret $D^+_h f(x)$ in the same way as $\frac{d}{dx} f(x)$,
i.e. as an operation applied to the function $f$ at the point $x$.
"""

# ╔═╡ 0d8522c0-8824-47d5-97d1-d2472e9c1adf
md"""
Without any surprise this can indeed be employed to approximate derivatives.
We will consider the function $f(x) = \sin(e^{x+1})$, or in code
"""

# ╔═╡ 40634b5e-d328-4da0-8274-48f166ffdbba
f(x) = sin(exp(x + 1))

# ╔═╡ 965625d6-b9da-4279-ad6d-54b2b62f1247
md"""
which has a slightly more tractable analytical derivative
compared to $v(t)$, so makes it easier for us to compare results.

We note $f'(x) = e^{x+1} \cos\left(e^{x + 1}\right)$, i.e. $f'(0) = e \cos(e)$ or
"""

# ╔═╡ 0c3831e6-0353-4dbc-927b-f55521ccd8e4
exact_value = ℯ * cos(ℯ)

# ╔═╡ 3b0f6640-5389-40f3-a16b-03dbaa14e3c6
md"""
Simply evaluating (1) with $h = 10^{-4}$ and $x = 0$ yields
"""

# ╔═╡ 832c0172-94b0-4298-92fe-071e71f1a1c9
let
	x = 0
	h = 1e-4
	D⁺ₕ = 1/h * (f(x + h) - f(x))
	error = abs(exact_value - D⁺ₕ)
	
	(D⁺ₕ, error)
end

# ╔═╡ 1151c2f5-7dd1-427f-8442-0332efd3f997
md"""
which is already a pretty good approximation to $f'(0)$. For a smaller $h = 10^{-6}$ the result is even better
"""

# ╔═╡ 16df9d30-879e-4f63-822d-34d57fc28433
let
	x = 0
	h = 1e-6
	D⁺ₕ = 1/h * (f(x + h) - f(x))
	error = abs(exact_value - D⁺ₕ)
		
	(D⁺ₕ, error)
end

# ╔═╡ 05022c02-5324-43ec-b86f-6dfa338c7fc4
md"Continuing this further we observe **linear convergence**."

# ╔═╡ 977c6fe3-8c0b-4dfc-b354-8fdab0c76884
md"""
Let us investigate the convergence behaviour more closely.
We consider a Talyor series of $f$ around $x$:
```math
f(x+h) = f(x) + f'(x)h + \frac12 f''(\xi) h^2 \qquad \text{for $\xi \in (x, x+h)$}
```
Considering (1) we thus obtain
```math
D^+_h f(x) = \frac{1}{h} \Big( f(x+h) - f(x) \Big)
= f'(x) + \frac12 f''(\xi) h
```
or
```math
\left| f'(x) - D^+_h f(x)\right| \leq \frac h2 \max_{x\in[a,b]} |f''(x)|
= \frac h2 \|f''\|_\infty
```
We just proved the linear convergence of $D^+_h$:
"""

# ╔═╡ 544bfd07-d57c-45f3-9a1a-37ea0bb7c0e3
md"""
!!! info "Theorem 1: Convergence of forward finite differences"
	Given $f : [a, b] \to \mathbb{R}$ a twice differentiable function, then
	the forward finite difference formula (1) converges linearly
    ```math
	\left| f'(x) - D^+_h f(x)\right| \leq α\, h
	```
	with constant $α = \frac12 \|f''\|_\infty$.
"""

# ╔═╡ 10295701-a493-43bc-ae0c-977b27151902
md"""
In building such finite difference formulas we are not restricted to two nodes and evaluations of $f$. In general one can imagine to approximate the derivative of $f$ at $x$ by taking $n$ points to the right of $x$, i.e. $x + h, x + 2h, \ldots, x + nh$ and $m$ points to its left, $x - mh, \ldots, x-2h, x-h$:
"""

# ╔═╡ 1d853f67-a95b-4ddf-b31d-c00f9fb92787
RobustLocalResource("https://raw.githubusercontent.com/epfl-matmat/numerical-analysis/3720244719e1b0a9299a1675315cc8c5a0fe3e5d/notes/img/finite-difference-nodes.svg", "img/finite-difference-nodes.svg", :width => 600)

# ╔═╡ 84595c02-e2bc-4ef7-8a24-eb501c4ef0e9
md"""
A general definition is:

!!! info "Definition: Finite differences formula"
    A **finite differences formula** for the $n+m+1$ equally equispaced **nodes**
	$x + i h$, $i = -m, \ldots, -1, 0, 1, \ldots, n$
	around $x \in [a, b]$ is an expression of the form
    ```math
	\tag{2}
    D_h \, f(x) = \frac{1}{h}\, \sum_{i=-m}^n w_i f(x + i h),
	```
	where the **coefficients** $w_{-m}, \ldots, w_0, \ldots w_n$
	do not depend on $h$.

	Such a formula is called **consistent** if
	for a regular function $f : [a, b] \to \mathbb{R}$
	it approximates its first derivative for $h\to 0$, i.e.
	```math
	f'(x) = \lim_{h\to0} D_h \, f(x).
	```
"""

# ╔═╡ 12c00f60-63d6-4be3-9e0b-45934e29e9ca
md"""
The forward differences formula (1) only employed the two nodes $x$ and $x + h$, i.e. $m=0$ and $n=1$. Two other formulas with only two nodes are frequently employed, namely:
- **Backward finite differences:** Instead of constructing the line through $x$ and $x+h$, one can also construct the line through $x-h$ and $x$, which leads to
  ```math
  \tag{3}
  D^-_h f(x) = \frac{1}{h} \Big(f(x) - f(x - h)\Big)
  ```
  i.e. a two-point formula with $m=-1$ and $n=0$.
- **Central finite differences:**
  ```math
  \tag{4}
  D^c_h f(x) = \frac{1}{2h} \Big(f(x+h) - f(x - h)\Big)
  ```
  One can think of this formula as the symmetrised form of the
  forward and backward formulas (1) and (3),
  obtained by adding half of (1) to half of (3).
"""

# ╔═╡ 4e996740-ac6f-4c83-87ed-f0e09aa45e89
md"""
Visualisation of the derivatives obtained by these methods as slopes:
- Exact derivative: $(@bind show_exact CheckBox(default=false))
- Forward finite differences: $(@bind show_forward CheckBox(default=true))
- Backward finite differences: $(@bind show_backward CheckBox(default=true))
- Central finite differences: $(@bind show_central CheckBox(default=true))
"""

# ╔═╡ 2ea81c7f-fd18-48c5-9676-2e460c66e35f
let
	x₀ = 0.5
	h  = 0.1

	fdash      = exp(x₀+1) * cos(exp(x₀+1))
	df_left    = 1/h  * (f(x₀ + h) - f(x₀))
	df_right   = 1/h  * (f(x₀) - f(x₀ - h))
	df_central = 1/2h * (f(x₀+h) - f(x₀-h))
	
	p = plot(f, xlims=(0.25, 0.75), lw=3, label="f", ylims=(-1.2, -0.7), legend=:bottomleft)
	if show_exact
	plot!(p, x -> f(x₀) + fdash * (x - x₀);
		  label="Exact derivative", lw=1.5, ls=:dash, c=2)
	end
		
	if show_forward
	plot!(p, x -> f(x₀) + df_left * (x - x₀);
		  label="Forward finite differences", lw=1.5, ls=:dash, c=4)
	end
	if show_backward
	plot!(p, x -> f(x₀) + df_right * (x - x₀);
		  label="Backward finite differences", lw=1.5, ls=:dash, c=3)
	end
	if show_central
	plot!(p, x -> f(x₀) + df_central * (x - x₀) + 0.091;
		  label="Central finite differences", lw=1.5, ls=:dash, c=5)
	end
	
	scatter!(p, [x₀],   [f(x₀)],   label="x", c=2)
	scatter!(p, [x₀-h], [f(x₀-h)], label="x-h", c=3)
	scatter!(p, [x₀+h], [f(x₀+h)], label="x+h", c=4)
end

# ╔═╡ 6c316ead-2b3e-47e3-b625-2b33ff410d89
md"""Let us consider the convergence of these three variants for computing the derivative of $f(x) = e^{\sin(x)}$ at $x=0$. First we construct the range of $h$ values:
"""

# ╔═╡ d9649df4-8507-4668-9c21-ab8946de21fc
hs = [10^e for e in 0:-0.5:-7];  # Range of h parameter

# ╔═╡ 2f5714fe-187e-4984-889b-fe3855ba390e
md"Then we compute the finite-difference approximations of $f'$ at $x$:"

# ╔═╡ 9076a083-e044-416b-a183-dc6fa906514a
begin
	x = 0
	deriv_forward  = [1/h    * (f(x+h) - f(x))   for h in hs]
	deriv_backward = [1/h    * (f(x)   - f(x-h)) for h in hs]
	deriv_center   = [1/(2h) * (f(x+h) - f(x-h)) for h in hs]
end;

# ╔═╡ 4d265f09-143c-4cd9-ac5c-aea0aad85bea
let
	p = plot(; yaxis=:log, xaxis=:log, xflip=true, ylims=(1e-9, 10),
		   xlabel=L"h", ylabel="absolute error")
	hs = [10^e for e in 0:-0.5:-7];
	deriv_forward = [1/h * (f(x + h) - f(x)) for h in hs]
	error_forward  = abs.(deriv_forward  .- exact_value)
	plot!(p, hs, error_forward;  label="forward",  lw=2, mark=:x)
end

# ╔═╡ 45126823-9ba8-4f89-a9c5-81789244ed99
md"Next we compute the errors agains. the reference $f'(x) = 1$
and plot the convergence in a log-log plot:"

# ╔═╡ e922280b-b0d2-4159-9597-0e1e1b71569d
let
	p = plot(; yaxis=:log, xaxis=:log, xflip=true, ylims=(1e-9, 10),
		   title="Convergence of finite-difference formulas", xlabel=L"h",
		   ylabel="absolute error")
	
	error_forward  = abs.(deriv_forward  .- exact_value)
	error_backward = abs.(deriv_backward .- exact_value)
	error_center   = abs.(deriv_center   .- exact_value)

	plot!(p, hs, error_forward;  label="forward",  lw=2, mark=:x)
	plot!(p, hs, error_backward; label="backward", lw=2, mark=:x)
	plot!(p, hs, error_center;   label="center",   lw=2, mark=:x)

	# Lines for perfect 1st and 2nd order convergence.
	plot!(p, hs, hs;       ls=:dash, label=L"O(h)",   lw=2)
	plot!(p, hs, 0.5hs.^2; ls=:dash, label=L"O(h^2)", lw=2)
end

# ╔═╡ dd5220d0-0290-490a-8309-c8b66164d332
md"""
We observe that both forward and backward finite differences converge at about the same rate. Indeed, one easily proves that the **backward finite differences** formula is also of **first order**.
"""

# ╔═╡ 7d48ab14-1638-4c03-90bc-188068183391
md"""
However, the center finite differences formula converges much faster. We again analyse using Taylor's formula:
```math
\begin{aligned}
D_h^c f(x) &= \frac{1}{2h} \left( f(x+h) - f(x-h) \right) \\
&= \frac{1}{2h} \left[
f(x) + h f'(x) + \frac{h^2}2 f''(x) + \frac{h^3}{6} f'''(\xi_1) \right.\\
&\hspace{4em} \left. - \left(
f(x) - h f'(x) + \frac{h^2}2 f''(x) - \frac{h^3}{6} f'''(\xi_2)
\right) \right] \\
&= \frac{1}{2h} \left[2\, h f'(x) + \frac{2h^3}{6} \left(f'''(\xi_1) + f'''(\xi_2)\right) \right]\\
&= f'(x) + \frac{h^2}{6} \left(f'''(\xi_1) + f'''(\xi_2)\right)\\
&≤ f'(x) + \frac{2 h^2}{6} \max_{x\in [a, b]} |f'''(x)|
\end{aligned}
```
where $\xi_1 \in (x, x+h)$ and $\xi_2 \in (x-h, x)$.
"""

# ╔═╡ 94e8d3cb-f664-4f26-b234-e955ab3f61cb
md"""
Therefore central finite differences is of second order:

!!! info "Theorem 2: Convergence of central finite differences"
	Given $f : [a, b] \to \mathbb{R}$ a three times differentiable function, then
	the central finite difference formula (3) converges **quadratically**
    ```math
	\left| f'(x) - D^c_h f(x)\right| \leq α \, h^2
	```
	with constant $α = \frac13 \|f'''\|_\infty = \frac13 \max_{x\in [a, b]} |f'''(x)|$.

"""

# ╔═╡ d28173db-0465-418f-bbd1-92b8fc64540c
md"""
In light of this discussion let us formalise the definition of convergence order for the context of finite differences formulas:

!!! info "Definition: Convergence order of finite differences"
	A finite differences formula $D_h \, f$ of the form (2)
	with equally spaced nodes of separation $h$ is of **order $p$**
	if a constant $α > 0$ indepentent of $h$
	(but possibly dependent on $f$) exists, such that
	```math
	\left| f'(x) - D_h f(x) \right| \leq α \, h^p
	```
	as long as the function $f$ is sufficiently regular.
"""

# ╔═╡ c6a94a08-3c8d-43e2-85e1-dc15d4f8b071
md"""
## Numerical stability

One of the key principles behind all finite difference formulas we considered was that they approximate the derivative, in the sense that $f'(x) = \lim_{h\to0} D_h \, f(x)$. As a result we would expect results to become more and more accurate as $h$ gets smaller.
"""

# ╔═╡ 0d89e57d-a8e8-405f-8962-e6377ab87842
md"""
We stick to our example function
```math
\begin{aligned}
f(x)  &= \sin(e^{x+1}) \\
f'(x) &= e^{x+1} \cos\left(e^{x + 1}\right)
\end{aligned}
```
where we consider the function on the interval $[a, b] = [-1, 1]$.
In particular we evaluate the derivative at $0$,
such that $f'(0) = e \, \cos(e)$.

Using the forward finite differences formula we compute
"""

# ╔═╡ 56636368-9ff6-4b45-a3be-8b1e05f6ee20
begin
	X = 0.0
	derivative_at_X = ℯ * cos(ℯ)
	
	all_h = [10^i for i in -1:-1.0:-14]
	all_error = Float64[]  # Error of D⁺ₕ
	for h in all_h
		D⁺ₕ = 1/h * (f(X + h) - f(X))
		error = abs(D⁺ₕ - derivative_at_X)
		@printf "h=%.1e   D⁺ₕ=%.12f   error=%.6e\n" h D⁺ₕ error
		push!(all_error, error)
	end	
end

# ╔═╡ 1e426a63-f3ea-46e7-b89d-1de6b2a0c421
md"""Plotting this error graphically yields:"""

# ╔═╡ e165dccd-c909-428f-95d6-9942fc8adce2
begin
	# Plot log-log plot of error versus step size:
	p = plot(all_h, all_error; xaxis=:log, yaxis=:log, mark=:x, lw=2, label="error", xflip=true, xlabel=L"h", ylims=(1e-9, 2), xlims=(1e-14, 1))

	# Plot guide line for linear convergence:
	plot!(p, all_h, 0.5all_h, label=L"O(h)", lw=1.5, ls=:dash)

	# Make y-axis labels denser
	yticks!(p, 10.0 .^ (-8:2:0))
end

# ╔═╡ 4bc09748-a008-4e62-b3f0-cc401b27af03
md"""
We notice that the derivative formula gives a good approximation to the exact derivative $f'(0) = e\,\cos(e)$ for $h = 10^{-8}$. However, if the node distance $h$ is **further descreased** the **approximation deteriorates**.
This is a result of the **round-off error in the finite-precision floating-point arithmetic** of the computer as we will discuss now in more detail.
"""

# ╔═╡ 77d5c198-ee33-40a4-a6e0-057c2e9c8281
md"""
We take another look at the forward finite difference formula
```math
D^+_h f(x) = \frac{f(x+h) - f(x)}{h}.
```
As $h$ gets smaller computing the difference $f(x+h) - f(x)$ becomes problematic.
Both the values of $f(x+h)$ and $f(x)$ can **only be computed to finite precision** in the computer's arithmetic. 
As these **values become more and more similar** (after all we take $h$ smaller and smaller) and this makes the **difference $f(x+h) - f(x)$
less and less precise**.
For example let us assume $16$ digits are accurate for both $f(x+h)$ and $f(x)$
and that we have chosen $h$ so small, that the first $10$ digits of both $f(x+h)$ and $f(x)$ agree.
Than taking the difference  $f(x+h) - f(x)$ will effectively nullify these first $10$ digits,
leaving us with only $6$ correct digits in the final answer.
We conclude: As **$h$ gets smaller**, the **difference $f(x+h) - f(x)$ has less and less correct digits**, in turn making the **numerical derivative** $D^+_h f(x)$
**less and less accurate**.
"""

# ╔═╡ 724f7887-9d5b-4b9c-b2fd-70e2a7062782
md"""
This rationalises our observations in the above plot,
but it does not yet tell us how we should choose the best $h$.
To answer this question we consider a more **quantitative mathematical analysis**.
When evaluate the function $f(x)$ using a numerical procedure
we always suffer from a small round-off error.
Instead of evaluating $f(x_1)$ the **computer thus actually evaluates**
```math
\tag{5}
\widetilde{f}(x_1) = f(x_1) (1 + ϵ_1),
```
where the round-off error $ϵ_1$ is on the order of $10^{-16}$.
Note, that $|ϵ_1|$ is the relative error of $\widetilde{f}$
```math
\frac{|\widetilde{f}(x_1) - f(x_1)|}{|{f}(x_1)|} = \frac{|ϵ_1| \, |f(x_1)|}{|f(x_1)|} = |ϵ_1|.
```
In general $ϵ_1$ will depend on $x_1$,
i.e. evaluating at a different point $x_2$ will lead to
a slightly different error $ϵ_2$.
However, standard floating-point furthermore come with the guarantee that
```math
|\epsilon_i| ≤ \epsilon_M
```
For double-precision floating-point numbers $\epsilon_M$ has the value
"""

# ╔═╡ 062f5835-42f5-4c57-a17b-ec3e617b542d
eps(Float64)

# ╔═╡ 6f4da5cf-9d05-42ec-8c27-1ae20632e94b
Foldable("Optional: Effect of additional error contributions when evaluating f",
md"""Note that on top of round-off errors there may be additional error contributions. A typical case is when evaluating $f(\mathbf x)$ involves itself an iterative procedure. For example if $f(\mathbf{x}) = \mathbf{z}$ where $\mathbf{z}$ is the solution to a linear systems $\mathbf A \mathbf z = \mathbf x$ computed using a conjugate gradient algorithm. In this case the iterative procedure is rarely solved perfectly, but stopped once a relative residual norm has been reached, which causes additional errors, effectively leading to an increase of the $\epsilon_M$ to be employed when working out the $h_\text{opt}$.
""")

# ╔═╡ 4f1ea545-4d37-4e78-b5de-6c8071e09d67
md"""
Employing this error model (5) within the
evaluation of the **first-order finite-difference formula** (1)
the computer will thus actually compute
```math
\begin{aligned}
\widetilde{D}^+_h f(x)
&= \frac{\widetilde{f}(x+h) - \widetilde{f}(x)}{h} \\
&= \frac{f(x+h) (1+ϵ_1) - f(x) (1+ϵ_2)}{h} \\
&= \frac{f(x+h) - f(x)}{h} + \frac{ϵ_1}{h}f(x+h) - \frac{ϵ_2}{h} f(x) \\
&= f'(x) + \frac{f''(\xi)}{2} h + \frac{ϵ_1}{h}f(x+h) - \frac{ϵ_2}{h} f(x),
\end{aligned}
```
where $|ϵ_1| ≤ \epsilon_M$, $|ϵ_2| ≤ \epsilon_M$
and $\xi \in (x, x+h)$. Collecting everything we obtain the error
of the computed finite-difference approximation to $f'(x)$ as
```math
\tag{6}
\begin{aligned}
|f'(x) - \widetilde{D}^+_h f(x)|
&= \left|\frac{f''(\xi)}{2} h + \frac{ϵ_1}{h}f(x+h) - \frac{ϵ_2}{h} f(x)\right| \\
&≤ \frac{h}{2} |f''(\xi)| + \frac{|ϵ_1|}{h} |f(x+h)| + \frac{|ϵ_2|}{h} |f(x)| \\
&≤ \frac{h}{2} \max_{x \in [a, b]} |f''(x)| + 2 \, \max_{x \in [a, b]} |f(x)| \, \frac{\epsilon_M}{h} \\
&= \underbrace{\frac{h}{2} \|f''\|_\infty}_{\text{trunc. error}}
+ \underbrace{\frac{2 \epsilon_M}{h} \|f\|_\infty}_{\text{round-off error}}
\end{aligned}
```
"""

# ╔═╡ 99974a63-9ba9-43c1-95ac-b02d83363e74
md"""
We notice there are **two error terms** in (6).
One is **proportional to $h$** (finite differences **truncation error**)
and one is **proportional to $1/h$**
(due to **round-off error**).
As $h\to0$ the first term thus derceases as the finite-difference approximation gets better. However, if $h$ is taken too small the second error growing as $1/h$ will dominate and our approximation will be bad.
"""

# ╔═╡ 896ee051-09d0-4849-8173-6f894ce66357
md"""
To obtain which $h$ gives the best sweet spot
we want to **balance both errors**. Utilising equation (6)
the total error of the approximation is bounded by
```math
\text{error}(h) = \frac{h}{2} \|f''\|_\infty + \frac{2 \epsilon_M}{h} \|f\|_\infty.
```
This function has a **single minimum**, which we can compute by
```math
0 = \frac{d\,\text{error}}{dh} \qquad \Rightarrow \qquad \frac{\|f''\|_\infty}{2} - \frac{2\, \epsilon_M \, \|f\|_\infty}{h^2} = 0
```
The **optimal node distance** $h_\text{opt}$, which **minimises the error**,
is thus
```math
\tag{7}
h_\text{opt} = \sqrt{\frac{4 \|f\|_\infty}{\|f''\|_\infty}} \, \sqrt{\epsilon_M}.
```
In our example we have $\|f\|_\infty = 1$ and $\|f''\|_\infty ≈ |f''(1)| ≈ e^4 ≈ 55$ and therefore $h_\text{opt}$ is of order $\sqrt{10^{-16}} = 10^{-8}$,
which we also observed numerically.
"""

# ╔═╡ 3cda6ecd-a0a2-452b-8cb6-b36432c0bdda
let
	# Plot log-log plot of error versus step size:
	p = plot(all_h, all_error; xaxis=:log, yaxis=:log, mark=:x, lw=2, label="error", xflip=true, xlabel=L"h", ylims=(1e-12, 2), xlims=(1e-14, 1))

	# Plot guide line for linear convergence:
	plot!(p, all_h, 0.5all_h, label=L"O(h)", lw=1.5, ls=:dash)

	plot!(p, all_h, 1e-16 ./ all_h, label=L"O(1/h)", lw=1.5, ls=:dash)


	# Make y-axis labels denser
	yticks!(p, 10.0 .^ (-12:2:0))
end

# ╔═╡ 31644da8-5501-477b-a65d-54f38c649718
md"""
**Optimal $h$ for higher-order formulas:**
Note that in (6) the $h$ dependence of the **first error term** (finite difference truncation error) depends on the **order of the finite difference formula**.

For an order $p$ method the error will thus have the form
```math
\text{error}(h) = C_1 h^p + C_2 \frac{\epsilon_M}{h}
```
with appropriate constants $C_1$ and $C_2$.
By a similar argument to minimise this error wrt. $h$
one can show that the optimal value of $h$ is
on the order of $\sqrt[p+1]{\epsilon_M}$     .
We summarise:

!!! info "Observation: Optimal nodal spacing h for finite differences"
    When computing a numerical derivative of $f$ using a **finite-difference
	method of order $p$**
	the optimal spacing of nodes satisfies roughly
	```math
	h_\text{opt} \approx \sqrt[p+1]{\epsilon_M}.
	```
	where $\epsilon_M$ is the maximal relative error in the evaluation of $f$.
"""

# ╔═╡ 86d719b4-494b-489d-9c6a-40c8d2477bad
md"""
To illustrate this graphically we apply three finite-difference (FD) formulas

- **1st order:**  $\displaystyle D_1\, g(x_0) = \frac{g(x_0+h) - g(x_0)}{h} \qquad \qquad$ *(forward finite differences)*
- **2nd order:**  $\displaystyle D_2\, g(x_0) = \frac{g(x_0+h) - g(x_0-h)}{2h} \qquad$ *(central finite differences)*
- **4th order:**  $\displaystyle D_4\, g(x_0) = \frac{g(x_0-2h) -8\, g(x_0-h) + 8\, g(x_0+h) - g(x_0+2h)}{12h}$

to the selected function to compute the derivative at at $x_0 = 0.2$.

- `g = ` $(@bind g Select([(x -> exp(-1.3x)) => "h(x) = exp(-1.3 x)", v => "v(t) = 64t * (1 − t) * (1 − 2t)^2 * (1 − 8t + 8t^2)^2", f => "f(x) = sin( exp(x+1) )"]))
"""

# ╔═╡ 512f4377-c80a-47da-a195-b63a746fcf75
md"We **compute a reference** using `ForwardDiff`, a more precise algorithmic technique to compute derivatives (outside of the scope of this course):"

# ╔═╡ 1469326c-92b1-4fad-a54f-e6397596c39f
x₀ = 0.2

# ╔═╡ 2db471c0-ebfb-4e30-9a95-f281d9467a4f
let
	p = plot(g, xlims=(-1, 1), label="Selected function g", lw=2, ylims=(-2, 2))
	scatter!(p, [x₀], [g(x₀)], label="x₀")
end

# ╔═╡ c7e410cc-aeba-4630-9ebb-87af7ef093a5
reference = ForwardDiff.derivative(g, x₀)

# ╔═╡ 6b30e9a0-b220-49c9-aa31-76c8d24d93b2
md"... and apply the three formulas:"

# ╔═╡ b0c56709-a77e-4417-bd2b-bf4543cb2a13
let
	error_fd1 = Float64[]
	error_fd2 = Float64[]
	error_fd4 = Float64[]

	all_h = [10^i for i in -1:-1.0:-14]	
	for h in all_h
		diff_fd1 = 1/h   * (g(x₀ + h) - g(x₀))
		diff_fd2 = 1/2h  * (g(x₀ + h) - g(x₀ - h))
		diff_fd4 = 1/12h * (g(x₀ -2h) - 8g(x₀-h) + 8g(x₀+h) - g(x₀+2h))

		push!(error_fd1, abs.(diff_fd1 - reference))
		push!(error_fd2, abs.(diff_fd2 - reference))
		push!(error_fd4, abs.(diff_fd4 - reference))
	end

	p = plot(; xaxis=:log, yaxis=:log, mark=:x, lw=2, xflip=true, xlabel=L"h", legend=:bottomright)
	plot!(p, all_h, error_fd1; label="1st order FD", mark=:x, lw=2)
	plot!(p, all_h, error_fd2; label="2nd order FD", mark=:x, lw=2,)
	plot!(p, all_h, error_fd4; label="4th order FD", mark=:x, lw=2)

	vline!(p, [sqrt(eps())], label=L"\sqrt{ε_M}", c=1, ls=:dash)
	vline!(p, [cbrt(eps())], label=L"ε_M^{1/3}",  c=2, ls=:dash)
	vline!(p, [eps()^(1/5)], label=L"ε_M^{1/5}",  c=3, ls=:dash)

	yticks!(p, 10. .^ (-14:2:0))
end

# ╔═╡ 2b940279-24e5-495b-af0c-b0b9e547f33a
md"""
**As $h$ shrinks** the **errors** are **initially** dominated by the **truncation error**,
which dicreases most rapidly for the 4th order formula.
However, the increasing **round-off error eventually dominates** the behaviour as the truncation error continues to decrease.

**As the order increases** the **crossover point** between truncation error and round-off error **moves further to the left** and further down. Thus **higher-order methods are generally more accurate** as the numerically problematic small $h$ values can be avoided.

The optimal value for $h$ does indeed scale with $\sqrt[p+1]{\epsilon_M}$ .
However, as a visual inspection of the error plot shows, **this formula does
only provide a rough orientation**: Sometimes the $h$ with lowest error can
in fact be smaller or larger.
In practice finding the best $h$ can still be rather challenging.
"""

# ╔═╡ 83be5e79-e71c-4f32-95e9-854f2fdc024d
md"""
## Construction of finite difference formulas

As we saw above, the delicate balance between truncation error and
round-off error implies that --- even when choosing the best $h$ ---
the **accuracy of low-order finite-difference formulas
can remain limited**.
In particular for accuracies beyond $10^{-6}$
second and higher-order formulas are almost always needed
and in fact a wide range of such formulas are available in the literature.

We will not discuss all details in this lecture and only sketch **two ideas how other finite difference formulas can be obtained** in practice.
For further discussion and a table of common
finite difference formulas
see for example [chapter 5.4](https://tobydriscoll.net/fnc-julia/localapprox/finitediffs.html#arbitrary-nodes) Driscoll, Brown: Fundamentals of Numerical Computation.
"""

# ╔═╡ cb27d648-dc8b-4840-925a-b122ce5d4fd5
md"""
### Determination of finite differences coefficients
The definition (2) of a finite difference formula already
introduced the general expression
```math
D_h \, f(x) = \frac{1}{h}\, \sum_{i=-m}^n w_i f(x + i h).
```
with nodes
"""

# ╔═╡ 8ab17ea3-a8e4-4471-9d54-9798a4b01a26
RobustLocalResource("https://raw.githubusercontent.com/epfl-matmat/numerical-analysis/3720244719e1b0a9299a1675315cc8c5a0fe3e5d/notes/img/finite-difference-nodes.svg", "img/finite-difference-nodes.svg", :width => 600)

# ╔═╡ a16bd249-755b-4b6a-b197-c299d1fd959d
md"""
The integers $m$ and $n$ determine the limits of the sum and thus the number of nodes at which one needs to evaluate the function.

Since **function evaluation is usually the expensive step**, 
and such a formula needs around $n + m + 1$ function evaluations,
the values of $m$ and $n$ should remain small.
Thus the values for $m$ and $n$ often need to be set due to practical limitations, such as the available computational time or the structure of the computational problem.

Assuming that $n$ and $m$ are given the
**main unknown** are thus the **coefficients $w_i$**.
Using a Taylor expansions of $f$
these can be obtained such that $D_h\,f(x)$ matches
$f'(x)$ as close as possible.
We consider an example:
"""

# ╔═╡ 75f3aa38-57b5-498f-b80e-97e28859316f
md"""
!!! warning "Example: m=n=1"
	We consider the case $m=n=1$, i.e. we want to derive the finite-differences formula using the three nodal points $x-h$, $x$ and $x+h$. For this setting the general formula (2) becomes
	```math
	\tag{8}
	{D}_h \, f(x) = \frac{1}{h}\left[
		w_{-1} f(x-h) + w_0 f(x) + w_1 f(x+h)
	\right].
	```
	Expanding $f(x-h)$ and $f(x+h)$ using Taylor series we have
	```math
	\begin{aligned}
	{D}_h \, f(x) &= \frac{w_{-1}}{h}\left[
		f(x) - f'(x) h + \frac{f''(x)}{2} h^2 + O(h^3)
	\right] \\
	&+ \frac{w_0}{h} f(x) \\
	&+ \frac{w_1}{h}\left[
		f(x) + f'(x) h + \frac{f''(x)}{2} h^2 + O(h^3)
	\right].
	\end{aligned}
	```
	Collecting coefficients we obtain
	```math
	{D}_h \, f(x) = \underbrace{\frac{w_{-1} + w_0 + w_1}{h}}_{=0}\, f(x)
		+ \underbrace{\left(w_1 - w_{-1}\right)}_{=1} \, f'(x)
		+ \underbrace{\frac{h\, (w_{-1} + w_1)}{2}}_{=0} \,f''(x)
		+ O(h^2),
	```
	which we want to be as close as possible to $f'(x)$.
	On the coefficients this imposes the conditions
	```math
	\left\{
	\begin{aligned}
	w_{-1} + w_0 + w_1 &= 0\\
	 - w_{-1} \phantom{+w_0} +w_{1} &= 1\\
	w_{-1} \phantom{+w_0} + w_1  &= 0
	\end{aligned}
	\right.
	\qquad \Leftrightarrow \qquad
	\underbrace{\begin{pmatrix} 1&1&1\\ -1&0&1 \\1&0&1 \end{pmatrix}}_{=\textbf{A}}
	\begin{pmatrix} w_{-1}\\ w_0 \\ w_1 \end{pmatrix}
	=
	\underbrace{\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}}_{=\textbf{b}}
	```
	The solution to this system can be easily computed as
"""

# ╔═╡ 478ad0b4-160c-4134-9ad7-3f93d754b131
let
	A = [ 1 1 1;
         -1 0 1;
          1 0 1]
	b = [0;
		 1;
		 0]
	w = A \ b
end

# ╔═╡ 99ff07c9-f2f0-474a-8c33-9a26fad3ee92
md"""
!!! warning ""
	i.e. $w_{-1} = -\frac12$, $w_0 = 0$, $w_1 = \frac12$,
	which we insert into (8) to obtain
	```math
	\tilde{D}_h \, f(x) = \frac{1}{h}\left(\frac{-1}2 f(x-h) + \frac12 f(x+h)\right)
	= \frac{1}{2h} \Big(f(x+h) - f(x - h)\Big).
	```
	We thus **recover the central finite differences formula** (4) from earlier.
"""

# ╔═╡ 37d7e297-82b4-4c15-b590-34c47cc41bb9
md"""
### Using interpolating polynomials

Recall our visualisation of the various finite differences formulas:
"""

# ╔═╡ 933aada3-991a-4a5c-ba6c-d2b85e5b62e1
let
	x₀ = 0.5
	h  = 0.1

	fdash      = exp(x₀+1) * cos(exp(x₀+1))
	df_left    = 1/h  * (f(x₀ + h) - f(x₀))
	df_right   = 1/h  * (f(x₀) - f(x₀ - h))
	df_central = 1/2h * (f(x₀+h) - f(x₀-h))
	
	p = plot(f, xlims=(0.25, 0.75), lw=3, label="f", ylims=(-1.2, -0.7), legend=:bottomleft)		
	plot!(p, x -> f(x₀) + df_left * (x - x₀);
		  label="Forward finite differences", lw=1.5, ls=:dash, c=4)
	plot!(p, x -> f(x₀) + df_right * (x - x₀);
		  label="Backward finite differences", lw=1.5, ls=:dash, c=3)

	plot!(p, x -> f(x₀) + df_central * (x - x₀) + 0.091;
		  label="Central finite differences", lw=1.5, ls=:dash, c=5)

	
	scatter!(p, [x₀],   [f(x₀)],   label="x", c=2)
	scatter!(p, [x₀-h], [f(x₀-h)], label="x-h", c=3)
	scatter!(p, [x₀+h], [f(x₀+h)], label="x+h", c=4)
end

# ╔═╡ e1d3d69d-380c-4dac-ab45-5e3517d542a6
md"""
In this image the value of the finite-difference derivative was the slope of the plotted lines. Most notably these lines are interpolating lines between two of the points $\big(x-h, f(x-h)\big)$, $\big(x, f(x)\big)$ or $\big(x+h, f(x+h)\big)$.

For example the forward finite differences formula (1)
```math
D^+_h f(x) = \frac{1}{h} \Big( f(x+h) - f(x) \Big) = \frac{f(x + h) - f(x)}{h}.
```
can be interpreted as the slope of a line going through
$\big(x, f(x)\big)$ and $\big(x+h, f(x+h)\big)$.
To put in another way we can interpret this formula
as the result of a two-step procedure:
1. **Interpolate a polynomial** through $\big(x, f(x)\big)$ and $\big(x+h, f(x+h)\big)$, leading to the linear polynomial interpolation $p_1$.
2. **Take the derivative of this interpolation** $p_1$ at $x$.
"""

# ╔═╡ 91bc8975-072b-4c98-9e81-311a155acab5
md"""
Indeed, using Lagrange polynomials we easily construct
the interpolating polynomial $p_1$ as
```math
\begin{aligned}
p_1(\xi) &= f(x) \frac{\xi - (x + h)}{x - (x + h)} + f(x + h) \frac{\xi - x}{x + h - x}\\
&= -f(x) \frac{\xi - (x + h)}{h} + f(x + h) \frac{\xi - x}{h}\\
&= \frac{f(x+h) - f(x)}{h} \xi + \frac{x f(x) + hf(x) -x f(x+h)}{h}
\end{aligned},
```
which has derivative
```math
p_1'(x) = \frac{f(x+h) - f(x)}{h}
```
as required.
"""

# ╔═╡ a3ddd54a-0171-4e7f-aa65-11224858f8be
md"""
This leads to a **natural generalisation** to obtain finite-difference formulas:
Interpolate a polynomial through $n+m+1$ nodes $\big((x+ih), f(x+ih)\big)$
for $i = -m, \ldots, n$,
leading to the $(n+m)$-th degree polynomial $p_{n+m}$.
Then take its derivative to obtain the finite differences formula as
```math
D_h\, f(x) = p'_{n+m}(x).
```
"""

# ╔═╡ ca3564d3-059e-40ce-bbcd-a04594e8ac24
md"""
!!! warning "Example n=m=1"
	We apply the procedure again to the case $m=n=1$ with the three nodal points $x-h$, $x$ and $x+h$. Using a Lagrange basis we find the second-degree interpolating
	polynomial $p_2$ through the points $\big(x-h, f(x-h)\big)$, $\big(x + 0\cdot h, f(x + 0\cdot h)\big)$ and $\big(x+h, f(x+h)\big)$ as
	```math
	\tag{9}
	\begin{aligned}
	p_2(\xi) &= f(\textcolor{red}{x-h}) \frac{\big(\xi-\textcolor{blue}{x}\big)\big(\xi - \textcolor{green}{(x+h)}\big)}{\big(\textcolor{red}{(x-h)} - \textcolor{blue}{x}\big) \big(\textcolor{red}{(x-h)} - \textcolor{green}{(x+h)}\big)} \\
			 &\hspace{3em}+ f(\textcolor{blue}{x})   \frac{\big(\xi - \textcolor{red}{(x-h)}\big) \big(\xi - \textcolor{green}{(x+h)}\big)}{\big(x - \textcolor{red}{(x-h)})(\textcolor{blue}{x} - \textcolor{green}{(x+h)}\big)} \\
			 &\hspace{3em}+ f(\textcolor{green}{x+h}) \frac{\big(\xi - \textcolor{red}{(x-h)}\big)\big(\xi - \textcolor{blue}{x}\big)}{\big(\textcolor{green}{(x+h)} - \textcolor{red}{(x-h)}\big) \big(\textcolor{green}{(x+h)} - \textcolor{blue}{x}\big)} \\
			&= f(x-h) \frac{\big(\xi-x\big)\big(\xi - x-h\big)}{2h^2} 
			 + f(x)   \frac{\big(\xi - x+h\big) \big(\xi - x-h\big)}{-h^2} \\
			 &\hspace{3em}+ f(x+h) \frac{\big(\xi - x+h\big)\big(\xi - x\big)}{2h^2} .
	\end{aligned}
	```
	Its derivative is
	```math
	\begin{aligned}
	p_2'(\xi) &= f(x-h)\frac{\xi - x + \xi - x - h}{2h^2}
	+ f(x) \frac{\xi - x+h + \xi - x-h}{-h^2}\\
	&+ f(x+h) \frac{\xi - x+h + \xi - x}{2h^2}
	\end{aligned}
	```
	such that at $\xi = x$ we get the finite-difference formula
	```math
	p_2'(x) = \frac{-f(x-h)}{2h} + \frac{f(x+h)}{2h},
	```
	which coincides again with the central finite difference formula (4).
"""

# ╔═╡ 56368f98-6f51-49be-8333-6fa1b9c64cf1
md"""
## Computing higher-order derivatives

The methods sketched in the aforementioned section also allow us to build finite-difference formulas to build higher-order derivatives.

For example, based on the interpolated polynomial (9) on the nodal points $x-h$, $x$ and $x+h$, we obtain a formula for **approximating the second derivatives** as
```math
\tag{10}
D^2_h f(x) = p_2''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2}.
```
This formula is also of **second order** as can be checked using a Taylor series:
!!! exercise
    Prove that $D_h^2f(x)$ approximates the second derivative of $f$ to second order.
"""

# ╔═╡ 2adfac60-ec91-4c70-ae0c-f460f47d37dc
md"""
!!! danger "Potential confusion: n-th derivative versus order n"
    Do not confuse the **derivation order** (how many times we differentiate) and the **approximation order** (the leading power of $h$ in the approximation error).

    In this case, we are approximating a second derivative, so the derivation order is 2, and it turns out that this formula has an approximation error in $\mathcal O (h^2)$ so the approximation order is also 2.
"""

# ╔═╡ 9c8ad9b0-50eb-478e-b63c-26a868765230
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 310)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
ForwardDiff = "~1.0.1"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.1"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.57"
Symbolics = "~6.38.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "9722ec798b822e3ffa8d656da2f354c706fd6116"

[[deps.ADTypes]]
git-tree-sha1 = "e2478490447631aedba0823d4d7a80b2cc8cdb32"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.14.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

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

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

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

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

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

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "9a3ae38b460449cc9e7dd0cfb059c76028724627"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.1"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

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
git-tree-sha1 = "910febccb28d493032495b7009dce7d7f7aee554"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.0.1"
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

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

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

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "453de0fc2be3d11b9b93ca4d0fddd91196dcf1ed"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.5"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "a414039192a155fb38c4599a60110f0018c6ec82"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.16.0"
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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "0e1340b5d98971513bddaa6bbed470670cebbbfe"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.34"

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
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

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
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

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

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrettyTables", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "7530e17b6ac652b009966f8ad53371a4ffd273f2"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.39"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper", "WeakValueDicts"]
git-tree-sha1 = "ca5929df933a8b7272bc7f01dcd14b2d976c56e2"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.25.1"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "e46dbf646bc3944c22a37745361c2e0a94f81d66"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.38.0"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

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

[[deps.TaskLocalValues]]
git-tree-sha1 = "d155450e6dff2a8bc2fcb81dcb194bd98b0aeb46"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.2"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

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
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

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

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

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
# ╟─4103c9d2-ef89-4c65-be3f-3dab59d1cc47
# ╠═9688d440-d231-11ee-0510-8db284ccd0c9
# ╟─e9151d3f-8d28-4e9b-add8-43c713f6f068
# ╟─bc9f0276-6b7b-4b19-8ef0-6371da964a9c
# ╠═7b74f71f-a5fa-47ba-b218-c1798fa540f7
# ╟─80df99a7-3e31-4418-9314-399116831e72
# ╟─87753dbe-f075-435b-9ac5-75bdcd6019c0
# ╟─0f0c0f37-f751-4950-9a8e-8dadb038325e
# ╟─657a065f-738f-4512-9a48-84fc9dae76da
# ╟─50cb2913-618f-4fd8-93a2-dd4f4411bbe1
# ╟─1df757af-8d84-429b-8005-b17e9814e48c
# ╟─0d8522c0-8824-47d5-97d1-d2472e9c1adf
# ╠═40634b5e-d328-4da0-8274-48f166ffdbba
# ╟─965625d6-b9da-4279-ad6d-54b2b62f1247
# ╠═0c3831e6-0353-4dbc-927b-f55521ccd8e4
# ╟─3b0f6640-5389-40f3-a16b-03dbaa14e3c6
# ╠═832c0172-94b0-4298-92fe-071e71f1a1c9
# ╟─1151c2f5-7dd1-427f-8442-0332efd3f997
# ╠═16df9d30-879e-4f63-822d-34d57fc28433
# ╟─05022c02-5324-43ec-b86f-6dfa338c7fc4
# ╟─4d265f09-143c-4cd9-ac5c-aea0aad85bea
# ╟─977c6fe3-8c0b-4dfc-b354-8fdab0c76884
# ╟─544bfd07-d57c-45f3-9a1a-37ea0bb7c0e3
# ╟─10295701-a493-43bc-ae0c-977b27151902
# ╟─1d853f67-a95b-4ddf-b31d-c00f9fb92787
# ╟─84595c02-e2bc-4ef7-8a24-eb501c4ef0e9
# ╟─12c00f60-63d6-4be3-9e0b-45934e29e9ca
# ╟─4e996740-ac6f-4c83-87ed-f0e09aa45e89
# ╟─2ea81c7f-fd18-48c5-9676-2e460c66e35f
# ╟─6c316ead-2b3e-47e3-b625-2b33ff410d89
# ╠═d9649df4-8507-4668-9c21-ab8946de21fc
# ╟─2f5714fe-187e-4984-889b-fe3855ba390e
# ╠═9076a083-e044-416b-a183-dc6fa906514a
# ╟─45126823-9ba8-4f89-a9c5-81789244ed99
# ╠═e922280b-b0d2-4159-9597-0e1e1b71569d
# ╟─dd5220d0-0290-490a-8309-c8b66164d332
# ╟─7d48ab14-1638-4c03-90bc-188068183391
# ╟─94e8d3cb-f664-4f26-b234-e955ab3f61cb
# ╟─d28173db-0465-418f-bbd1-92b8fc64540c
# ╟─c6a94a08-3c8d-43e2-85e1-dc15d4f8b071
# ╟─0d89e57d-a8e8-405f-8962-e6377ab87842
# ╠═56636368-9ff6-4b45-a3be-8b1e05f6ee20
# ╟─1e426a63-f3ea-46e7-b89d-1de6b2a0c421
# ╟─e165dccd-c909-428f-95d6-9942fc8adce2
# ╟─4bc09748-a008-4e62-b3f0-cc401b27af03
# ╟─77d5c198-ee33-40a4-a6e0-057c2e9c8281
# ╟─724f7887-9d5b-4b9c-b2fd-70e2a7062782
# ╠═062f5835-42f5-4c57-a17b-ec3e617b542d
# ╟─6f4da5cf-9d05-42ec-8c27-1ae20632e94b
# ╟─4f1ea545-4d37-4e78-b5de-6c8071e09d67
# ╟─99974a63-9ba9-43c1-95ac-b02d83363e74
# ╟─896ee051-09d0-4849-8173-6f894ce66357
# ╟─3cda6ecd-a0a2-452b-8cb6-b36432c0bdda
# ╟─31644da8-5501-477b-a65d-54f38c649718
# ╟─86d719b4-494b-489d-9c6a-40c8d2477bad
# ╟─2db471c0-ebfb-4e30-9a95-f281d9467a4f
# ╟─512f4377-c80a-47da-a195-b63a746fcf75
# ╠═1469326c-92b1-4fad-a54f-e6397596c39f
# ╠═c7e410cc-aeba-4630-9ebb-87af7ef093a5
# ╟─6b30e9a0-b220-49c9-aa31-76c8d24d93b2
# ╟─b0c56709-a77e-4417-bd2b-bf4543cb2a13
# ╟─2b940279-24e5-495b-af0c-b0b9e547f33a
# ╟─83be5e79-e71c-4f32-95e9-854f2fdc024d
# ╟─cb27d648-dc8b-4840-925a-b122ce5d4fd5
# ╟─8ab17ea3-a8e4-4471-9d54-9798a4b01a26
# ╟─a16bd249-755b-4b6a-b197-c299d1fd959d
# ╟─75f3aa38-57b5-498f-b80e-97e28859316f
# ╠═478ad0b4-160c-4134-9ad7-3f93d754b131
# ╟─99ff07c9-f2f0-474a-8c33-9a26fad3ee92
# ╟─37d7e297-82b4-4c15-b590-34c47cc41bb9
# ╟─933aada3-991a-4a5c-ba6c-d2b85e5b62e1
# ╟─e1d3d69d-380c-4dac-ab45-5e3517d542a6
# ╟─91bc8975-072b-4c98-9e81-311a155acab5
# ╟─a3ddd54a-0171-4e7f-aa65-11224858f8be
# ╟─ca3564d3-059e-40ce-bbcd-a04594e8ac24
# ╟─56368f98-6f51-49be-8333-6fa1b9c64cf1
# ╟─2adfac60-ec91-4c70-ae0c-f460f47d37dc
# ╟─9c8ad9b0-50eb-478e-b63c-26a868765230
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
