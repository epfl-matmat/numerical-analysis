### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 9688d440-d231-11ee-0510-8db284ccd0c9
begin
	using Plots
	using PlutoUI
	using PlutoTeachingTools
	using Printf
	using LaTeXStrings
	using HypertextLiteral
end

# ╔═╡ 4103c9d2-ef89-4c65-be3f-3dab59d1cc47
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/10_Numerical_differentiation.pdf)
"""

# ╔═╡ e9151d3f-8d28-4e9b-add8-43c713f6f068
TableOfContents()

# ╔═╡ d4e538c0-5529-4ec3-bc11-ec797bda8cfc
md"""
# Numerical differentiation

In the last chapter we discussed numerical techniques for integration based on quadrature formulas. In this chapter we will consider numerical techniques for performing the inverse operation, namely differentiation.
"""

# ╔═╡ fe24a732-672f-470d-9174-315df4e4a7ca
md"""
## First order derivatives

Given a regular function $f : [a, b] \to \mathbb{R}$ our goal is to approximate numerically its derivative $f'$ in a point $x \in [a, b]$. Similar to the case of numerical integration we will allow ourselves only to perform $n$ pointwise evaluations $f(t_i)$ with $t_i \in [a, b]$ and $i = 1, \ldots, n$ and take linear combinations of these results.

A first idea goes back to the definition of the derivative $f'(x)$ as the limit $h\to 0$ of the slope of secants over an intervall $[x, x+h]$, i.e.
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

Without any surprise this can indeed be employed to approximate derivatives.
Let us consider the function $f(x) = \sin\left(e^{x + 1}\right)$ or
"""

# ╔═╡ 40634b5e-d328-4da0-8274-48f166ffdbba
f(x) = sin(exp(x + 1))

# ╔═╡ 965625d6-b9da-4279-ad6d-54b2b62f1247
md"""
which has analytical derivative $f'(x) = e^{x+1} \cos\left(e^{x + 1}\right)$, i.e. $f'(0) = e \cos(e)$ or
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

# ╔═╡ 1d818a50-65a1-4788-985b-c7a1c9e9e5b0
md"""
This experiment already points to **linear convergence** (Why ?),
so let us investigate the convergence behaviour more closely.
We consider a Talyor serios of $f$ around $x$:
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

!!! info "Theorem 1: Convergence of forward finite differences"
	Given $f : [a, b] \to \mathbb{R}$ a $C^2$ function, then
	the forward finite difference formula (1) converges linearly
    ```math
	\left| f'(x) - D^+_h f(x)\right| \leq C h
	```
	with constant $C = \frac12 \|f''\|_\infty$.
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

# ╔═╡ 63f156df-afa5-4f1d-9091-dc5d4e8b1ce2
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

# ╔═╡ 0e0936aa-d759-4fe7-8375-89b9f7ea37b4
TODO("Some image explaining these formulas using the slope of lines between points would be good.")

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

# ╔═╡ 48442160-36e1-4604-9a7a-c8bc6014494a
md"""
We observe that both forward and backward finite differences converge at about the same rate. Indeed, one easily proves that the **backward finite differences** formula is also of **first order**.

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
where $\xi_1 \in (x, x+h)$ and $\xi_2 \in (x-h, x)$. Therefore central finite differences is of second order:

!!! info "Theorem 2: Convergence of central finite differences"
	Given $f : [a, b] \to \mathbb{R}$ a $C^3$ function, then
	the central finite difference formula (3) converges **quadratically**
    ```math
	\left| f'(x) - D^c_h f(x)\right| \leq C h^2
	```
	with constant $C = \frac16 \|f'''\|_\infty = \frac16 \max_{x\in [a, b]} |f'''(x)|$.

"""

# ╔═╡ d28173db-0465-418f-bbd1-92b8fc64540c
md"""
In light of this discussion let us formalise the definition of convergence order for the context of finite differences formulas:

!!! info "Definition: Convergence order of finite differences"
	A finite differences formula $D_h \, f$ of the form (2)
	with equally spaced nodes of separation $h$ is of **order $p$**
	if a constant $C$ indepentent of $h$
	(but possibly dependent on $f$) exists, such that
	```math
	\left| f'(x) - D_h f(x) \right| \leq C h^p
	```
	as long as the function $f$ is sufficiently regular.
"""

# ╔═╡ 83be5e79-e71c-4f32-95e9-854f2fdc024d
md"""
## Construction of finite difference formulas

Above we introduced three finite difference formulas, which have been motivated directly from the definition of of the derivative itself and all employed only two points. The construction of more involved formulas is feasible and in fact frequently employed to obtain an even higher convergence order.

We will not discuss all details in this lecture and only sketch two ideas how other finite difference formulas can be obtained in practice. For further discussion and a table of common
finite difference formulas
see for example [chapter 5.4](https://tobydriscoll.net/fnc-julia/localapprox/finitediffs.html#arbitrary-nodes) Driscoll, Brown: Fundamentals of Numerical Computation.
"""

# ╔═╡ dfd068d8-d8e6-4b66-a51a-7da122a95362
md"""
### Determination of finite differences coefficients
The definition (2) of a finite difference formula already
introduced the general expression
```math
\tilde{D}_h \, f(x) = \frac{1}{h}\, \sum_{i=-m}^n w_i f(x + i h).
```
The integers $m$ and $n$ determine the limits of the sum and thus the number of nodes at which one needs to evaluate the function.

Since function evaluation is usually the expensive step, the values for $m$ and $n$ are usually set based on practical limitations (such as the available computational time or the structure of the computational problem). The key unknown is thus to determine the coefficients $w_i$, which can be obtained by matching $\tilde{D}_h\,f(x)$ as closely as possible to $f'(x)$ using Taylor expansions of $f$. We consider an example:

!!! warning "Example: m=n=1"
	We consider the case $m=n=1$, i.e. we want to derive the finite-differences formula using the three nodal points $x-h$, $x$ and $x+h$. For this setting the general formula (2) becomes
	```math
	\tag{5}
	\tilde{D}_h \, f(x) = \frac{1}{h}\left[
		w_{-1} f(x-h) + w_0 f(x) + w_1 f(x+h)
	\right].
	```
	Expanding $f(x-h)$ and $f(x+h)$ using Taylor series we have
	```math
	\begin{aligned}
	\tilde{D}_h \, f(x) &= \frac{w_{-1}}{h}\left[
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
	\tilde{D}_h \, f(x) = \underbrace{\frac{w_{-1} + w_0 + w_1}{h}}_{=0}\, f(x)
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
	which we insert into (5) to obtain
	```math
	\tilde{D}_h \, f(x) = \frac{1}{h}\left(\frac{-1}2 f(x-h) + \frac12 f(x+h)\right)
	= \frac{1}{2h} \Big(f(x+h) - f(x - h)\Big).
	```
	We thus recover the central finite differences formula (4) from earlier.
"""

# ╔═╡ 85e2c3d5-3e63-460a-83ae-6c4bce19d9d1
md"""
### Using interpolating polynomials

Similar to the construction of quadrature formulas,
we can also obtain finite difference formulas based on interpolating polynomials.

To rationalise this idea recall the forward finite differences formula (1)
```math
D^+_h f(x) = \frac{1}{h} \Big( f(x+h) - f(x) \Big) = \frac{f(x + h) - f(x)}{h}.
```
We already noted that this equation can be interpreted
as the slope of a line going through the points
$\big(x, f(x)\big)$ and $\big(x+h, f(x+h)\big)$.
Therefore it is equal to the derivative $p_1'(x)$
of a (linear) polynomial interpolating
$\big(x, f(x)\big)$ and $\big(x+h, f(x+h)\big)$
as can be easily verified.
Indeed, using Lagrange polynomials we easily construct
the interpolating polynomial as
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

This procedure can be generalised:
Given the $n+m+1$ nodes $x + ih$ for $i = -m, \ldots, n$
of a finite differences formula,
we first construct the $(n+m)$-th degree polynomial $p_{n+m}$
interpolating the data $\big((x+ih), f(x+ih)\big)$. Then we build the finite differences formula as
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
	\tag{6}
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

For example, based on the interpolated polynomial (6) on the nodal points $x-h$, $x$ and $x+h$, we obtain a formula for **approximating the second derivatives** as
```math
\tag{7}
D^2_h f(x) = p_2''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2}.
```
This formula is also of **second order** as can be checked easily using a Taylor series.
"""

# ╔═╡ 19d8b680-6992-4c9f-aaf4-c665600a27ce
md"""
## Numerical stability

One of the key principles behind all finite difference formulas we considered was that they approximate the derivative, in the sense that $f'(x) = \lim_{h\to0} D_h \, f(x)$. As a result we would expect results to become more and more accurate as $h$ gets smaller.

Let's check this with the function
"""

# ╔═╡ ca5b9f56-9494-4720-a861-9c3627761146
g(x) = exp(-1.3x)

# ╔═╡ 3dbd632b-7330-464e-b3c2-db138daad238
md"""
which at $x = 0$ clearly has the derivative $g'(0) = -1.3$.
Using the forward finite differences formula we compute
"""

# ╔═╡ 56636368-9ff6-4b45-a3be-8b1e05f6ee20
begin
	X = 0.0
	reference = -1.3  # g'(X) = g'(0)
	
	all_h = [10^i for i in -1:-1.0:-14]
	all_error = Float64[]  # Error of D⁺ₕ
	for h in all_h
		D⁺ₕ = 1/h * (g(X + h) - g(X))
		error = abs(D⁺ₕ - reference)
		@printf "h=%.1e   D⁺ₕ=%.12f   error=%.6e\n" h D⁺ₕ error
		push!(all_error, error)
	end	
end

# ╔═╡ 1e426a63-f3ea-46e7-b89d-1de6b2a0c421
md"""Plotting this error graphically yields:"""

# ╔═╡ e165dccd-c909-428f-95d6-9942fc8adce2
begin
	plot(all_h, all_error; xaxis=:log, yaxis=:log, mark=:x, lw=2, label="error", xflip=true, xlabel=L"h", ylims=(1e-12, 2), xlims=(1e-14, 1))
	plot!(all_h, 0.2all_h, label=L"O(h)", lw=1.5, ls=:dash)
end

# ╔═╡ 6988bfd1-6664-47e9-874c-69437924bec0
md"""
We notice that the derivative formula gives a good approximation to the exact derivative $g'(0) = -1.3$ for $h = 10^{-8}$. However, if the node distance $h$ is further descreased the approximation deteriorates. This is a result of the round-off error in the finite-precision floating-point arithmetic of the computer. We take another look at the forward finite difference formula
```math
D^+_h f(x) = \frac{f(x+h) - f(x)}{h}.
```
In particular computing the difference  $f(x+h) - f(x)$ becomes problematic as $h$ gets smaller, since this implies that the values of $f(x+h)$ and $f(x)$ get more and more similar. As a result the difference $f(x+h) - f(x)$ (and thus the overall numerical derivative $D^+_h f(x)$) become less and less precise.

Let us quantify this effect more precisely. When computing the function $f(x)$ the employed algorithm will always suffer from a small round-off error[^1].
Instead of evaluating $f(x)$ the computer thus actually evaluates
```math
\hat{f}(x) = f(x) (1 + \eta),
```
where the round-off error $\eta$ is on the order of $10^{-16}$.
Note that in standard floating-point formats the round-off
error is always a *relative error*.
The precise value of $\eta$ in general depends on $x$,
but we further have $|\eta| ≤ \epsilon_M$
where for standard double-precision
floating-point numbers $\epsilon_M$ is about $2 \cdot 10^{-16}$.
"""

# ╔═╡ 4ac1b1e3-e570-4c00-93e3-e25fa14cba85
md"""
When evaluating the finite-difference
formula the computer will thus actually evaluate
```math
\begin{aligned}
\widehat{D}^+_h f(x)
&= \frac{\hat{f}(x+h) - \hat{f}(x)}{h} \\
&= \frac{f(x+h) (1+\eta_1) - f(x) (1+\eta_2)}{h} \\
&= \frac{f(x+h) - f(x)}{h} + \frac{\eta_1}{h}f(x+h) - \frac{\eta_2}{h} f(x) \\
&= f'(x) + \frac{f''(\xi)}{2} h + \frac{\eta_1}{h}f(x+h) - \frac{\eta_2}{h} f(x),
\end{aligned}
```
where $|\eta_1| ≤ \epsilon_M$, $|\eta_2| ≤ \epsilon_M$
and $\xi \in (x, x+h)$. Collecting everything we obtain the error
of the computed finite-difference approximation to $f'(x)$ as
```math
\tag{8}
\begin{aligned}
|f'(x) - \widehat{D}^+_h f(x)|
&≤ \frac{h}{2} \max_{x \in [a, b]} |f''(x)| + 2 \, \max_{x \in [a, b]} |f(x)| \, \frac{\epsilon_M}{h} \\
&= \underbrace{\frac{h}{2} \|f''\|_\infty}_{\text{trunc. error}}
+ \underbrace{\frac{2 \epsilon_M}{h} \|f\|_\infty}_{\text{round-off error}}
\end{aligned}
```
We notice there are two error terms in (8).
One is proportional to $h$ (finite differences truncation error)
and one is proportional to $1/h$
(due to round-off errors).
As $h\to0$ the first term thus derceases as the finite-difference approximation gets better. However, if $h$ is taken too small the second error growing as $1/h$ will dominate and our approximation will be bad.

A natural question is to ask is which $h$ gives the best sweet spot
and balances both errors. Utilising equation (8)
the total error of the approximation is bounded by
```math
ε(h) = \frac{h}{2} \|f''\|_\infty + \frac{2 \epsilon_M}{h} \|f\|_\infty.
```
This function has a single minimum, which we can compute by
```math
0 = \frac{dε}{dh} \qquad \Rightarrow \qquad \frac{\|f''\|_\infty}{2} - \frac{2\, \epsilon_M \, \|f\|_\infty}{h^2} = 0
```
The optimal node distance $h_\text{opt}$, which minimises the error, is thus
```math
\tag{9}
h_\text{opt} = \sqrt{\frac{4 \|f\|_\infty}{\|f''\|_\infty}} \, \sqrt{\epsilon_M}.
```
In our example we have $\|f\|_\infty ≈ |f(0)| = 1$ and $\|f''\|_\infty ≈ |f''(0)| = 1.69$ and therefore $h_\text{opt}$ is of order $\sqrt{10^{-16}} = 10^{-8}$,
which we also observed numerically.

[^1]: Note that on top of round-off errors there may be additional error contributions. A typical case is when evaluating $f(\mathbf x)$ involves itself an iterative procedure, e.g. if $f(\mathbf{x}) = \mathbf{z}$ and $\mathbf{z}$ is the solution to a linear systems $\mathbf A \mathbf z = \mathbf x$ computed using a conjugate gradient algorithm. In this case the iterative procedure is rarely solved perfectly, but stopped once a relative residual norm has been reached, which causes additional error contributions.
"""

# ╔═╡ 7f04041d-1077-4d8a-91ed-a794379e55de
TODO("Show the theoretical error behaviour in the plot")

# ╔═╡ 31644da8-5501-477b-a65d-54f38c649718
md"""
Note that in (8) the $h$ dependence of the first error term (finite difference truncation error) depends on the order of the finite difference formula.
For an order $p$ method the error will thus have the form
```math
ε(h) = C_1 h^p + C_2 {\epsilon_M}{h}
```
with appropriate constants $C_1$ and $C_2$,
which leads to an optimal value of order $\sqrt[p+1]{\epsilon_M}$  .
We summarise:

!!! info "Observation: Optimal nodal spacing h for finite differences"
    For computing an approximate derivative with a finite-difference
	method of order $p$ in the presence of round-off errors in the
	function evaluations $f$, the optimal spacing of nodes satisfies
	roughly
	```math
	h_\text{opt} \approx \sqrt[p+1]{\epsilon_M}.
	```
"""

# ╔═╡ 86d719b4-494b-489d-9c6a-40c8d2477bad
md"""
To illustrate this graphically we continue our example of computing the derivative
of $g(x) = e^{-1.3 x}$ at $x=0$,
which has the reference result
"""

# ╔═╡ 08574e94-69f5-4a6e-8ffd-82d8e08be082
reference

# ╔═╡ 2b940279-24e5-495b-af0c-b0b9e547f33a
md"""
We consider three finite difference (FD) formulas, namely

- **1st order:**  $\displaystyle D_1\, g(x) = D_h^+\, g(x) = \frac{g(x+h) - g(x)}{h} \qquad$ (forward finite differences)
- **2nd order:**  $\displaystyle D_2\, g(x) = D_h^c\, g(x) = \frac{g(x+h) - g(x-h)}{2h} \qquad$ (central finite differences)
- **4th order:**  $\displaystyle D_4\, g(x) = \frac{g(x-2h) -8\, g(x-h) + 8\, g(x+h) - g(x+2h)}{12h}$

We see that as $h$ shrinks the errors are initially dominated by the truncation error,
which dicreases most rapidly for the 4th order formula. However, the increasing round-off error eventually dominates the behaviour as the truncation error continues to decrease. As the order increases the crossover point between truncation error and round-off error moves further to the left and further down. Thus higher-order methods are generally more accurate as the numerically problematic small $h$ values can be avoided.
"""

# ╔═╡ b0c56709-a77e-4417-bd2b-bf4543cb2a13
let
	error_fd1 = Float64[]
	error_fd2 = Float64[]
	error_fd4 = Float64[]
	
	for h in all_h
		diff_fd1 = 1/h   * (g(X + h) - g(X))
		diff_fd2 = 1/2h  * (g(X + h) - g(X - h))
		diff_fd4 = 1/12h * (g(X -2h) - 8g(X-h) + 8g(X+h) - g(X+2h))

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
end

# ╔═╡ 9c8ad9b0-50eb-478e-b63c-26a868765230
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 350)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Plots = "~1.40.1"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.57"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "4a6beb87deb138ad5ddbc7a1b079f1152373f7c8"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

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
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

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
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "d2c021fbdde94f6cdaa799639adfeeaa17fd67f5"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.13.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "9c4708e3ed2b799e6124b5673a712dda0b596a9b"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "ac67408d9ddf207de5cfa9a97e114352430f01ed"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.16"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

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
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "3458564589be207fa6a77dbbf8b97674c9836aab"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "77f81da2964cc9fa7c0127f941e8bce37f7f1d70"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.2+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ac7b73d562b8f4287c3b67b4c66a5395a19c1ae8"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "7b762d81887160169ddfc93a47e5fd7a6a3e78ef"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.29"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

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

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "20ce1091ba18bcdae71ad9b71ee2367796ba6c48"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.4"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

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
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

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
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

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
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

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
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "c4fa93d7d66acad8f6f4ff439576da9d2e890ee0"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.1"

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
git-tree-sha1 = "89f57f710cc121a7f32473791af3d6beefc59051"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.14"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "a6783c887ca59ce7e97ed630b74ca1f10aefb74d"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.57"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "3fe4e5b9cdbb9bbc851c57b149e516acc07f8f72"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.13"

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

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
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
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

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

[[deps.TranscodingStreams]]
git-tree-sha1 = "54194d92959d8ebaa8e26227dbe3cdefcdcd594f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.3"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

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
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

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
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

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
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

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
git-tree-sha1 = "873b4f805771d3e4bafe63af759a26ea8ca84d14"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.42+0"

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
# ╟─4103c9d2-ef89-4c65-be3f-3dab59d1cc47
# ╠═9688d440-d231-11ee-0510-8db284ccd0c9
# ╟─e9151d3f-8d28-4e9b-add8-43c713f6f068
# ╟─d4e538c0-5529-4ec3-bc11-ec797bda8cfc
# ╟─fe24a732-672f-470d-9174-315df4e4a7ca
# ╠═40634b5e-d328-4da0-8274-48f166ffdbba
# ╟─965625d6-b9da-4279-ad6d-54b2b62f1247
# ╠═0c3831e6-0353-4dbc-927b-f55521ccd8e4
# ╟─3b0f6640-5389-40f3-a16b-03dbaa14e3c6
# ╠═832c0172-94b0-4298-92fe-071e71f1a1c9
# ╟─1151c2f5-7dd1-427f-8442-0332efd3f997
# ╠═16df9d30-879e-4f63-822d-34d57fc28433
# ╟─1d818a50-65a1-4788-985b-c7a1c9e9e5b0
# ╟─10295701-a493-43bc-ae0c-977b27151902
# ╟─1d853f67-a95b-4ddf-b31d-c00f9fb92787
# ╟─84595c02-e2bc-4ef7-8a24-eb501c4ef0e9
# ╟─63f156df-afa5-4f1d-9091-dc5d4e8b1ce2
# ╠═0e0936aa-d759-4fe7-8375-89b9f7ea37b4
# ╟─6c316ead-2b3e-47e3-b625-2b33ff410d89
# ╠═d9649df4-8507-4668-9c21-ab8946de21fc
# ╟─2f5714fe-187e-4984-889b-fe3855ba390e
# ╠═9076a083-e044-416b-a183-dc6fa906514a
# ╟─45126823-9ba8-4f89-a9c5-81789244ed99
# ╠═e922280b-b0d2-4159-9597-0e1e1b71569d
# ╟─48442160-36e1-4604-9a7a-c8bc6014494a
# ╟─d28173db-0465-418f-bbd1-92b8fc64540c
# ╟─83be5e79-e71c-4f32-95e9-854f2fdc024d
# ╟─dfd068d8-d8e6-4b66-a51a-7da122a95362
# ╠═478ad0b4-160c-4134-9ad7-3f93d754b131
# ╟─99ff07c9-f2f0-474a-8c33-9a26fad3ee92
# ╟─85e2c3d5-3e63-460a-83ae-6c4bce19d9d1
# ╟─ca3564d3-059e-40ce-bbcd-a04594e8ac24
# ╟─56368f98-6f51-49be-8333-6fa1b9c64cf1
# ╟─19d8b680-6992-4c9f-aaf4-c665600a27ce
# ╠═ca5b9f56-9494-4720-a861-9c3627761146
# ╟─3dbd632b-7330-464e-b3c2-db138daad238
# ╠═56636368-9ff6-4b45-a3be-8b1e05f6ee20
# ╟─1e426a63-f3ea-46e7-b89d-1de6b2a0c421
# ╠═e165dccd-c909-428f-95d6-9942fc8adce2
# ╟─6988bfd1-6664-47e9-874c-69437924bec0
# ╟─4ac1b1e3-e570-4c00-93e3-e25fa14cba85
# ╠═7f04041d-1077-4d8a-91ed-a794379e55de
# ╟─31644da8-5501-477b-a65d-54f38c649718
# ╟─86d719b4-494b-489d-9c6a-40c8d2477bad
# ╠═08574e94-69f5-4a6e-8ffd-82d8e08be082
# ╟─2b940279-24e5-495b-af0c-b0b9e547f33a
# ╠═b0c56709-a77e-4417-bd2b-bf4543cb2a13
# ╟─9c8ad9b0-50eb-478e-b63c-26a868765230
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
