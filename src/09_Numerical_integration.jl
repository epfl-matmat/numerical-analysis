### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ dc68915e-d230-11ee-09c8-4fce8a858118
begin
	using Plots
	using PlutoUI
	using PlutoTeachingTools
	using LaTeXStrings
	using QuadGK
	using Printf
	using HypertextLiteral
end

# ╔═╡ d34833b7-f375-40f7-a7a6-ab925d736320
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/09_Numerical_integration.pdf)
"""

# ╔═╡ 47114a9b-0e74-4e48-bb39-b49f526f1e9b
TableOfContents()

# ╔═╡ 85c29f83-4617-41b7-90e4-f863c3f09b1d
md"""
# Numerical integration

Integration is an important operation in engineering and the physical sciences, but evaluating integrals analytically can be quite challenging. For such cases employing *numerical* integration techniques is a good alternative, which sometimes turns out to be no less accurate.

Let's start with an easy analytic case. The integral $\int_0^1 e^x dx$ can be computed rather elegantly by analytical means, since the anti-derivative of $e^x$ is again $e^x$. Therefore
```math
\int_0^1 e^x dx = [e^x]_0^1 = e^1 - 1 \simeq 1.718281828459045.
```
"""

# ╔═╡ 10e9343f-4f0b-4d95-9c5a-3e7a6f6b9794
exact = exp(1) - 1

# ╔═╡ 9a300e03-8ec2-41ed-85c2-b7c23c140010
md"""
To perform the integration numerically, we can employ the Julia package `QuadGK`, which implements an all-purpose numerical integration routine:
"""

# ╔═╡ 53347e5e-00b1-4317-8912-331e90c3aa07
let
	Q, error_estimate = quadgk(x -> exp(x), 0, 1)
	Q
end

# ╔═╡ 30b6c88e-ebf8-49b4-acbb-9dacdb650b28
md"""
As we observe, the difference to the analytical result is tiny and on the order of the floating-point precision.

However, the numerical approach is much more flexible. For example $e^{\sin(x)}$ has no useful anti-derivative. Still the numerical approach works just as well to compute $\int_0^1 e^{\sin(x)} dx$:
"""

# ╔═╡ f95ba764-ac44-45e0-a714-3fffe2066c7f
let
	Q, error_estimate = quadgk(x -> exp(sin(x)), 0, 1)
	Q
end

# ╔═╡ bfccad43-b336-4271-83eb-40af1ed36849
md"""
What is remarkable when looking at the graphs of these two functions is that they are very similar ... and for one case the area under the curve boils down to simple basic calculus and the other is impenetrable analytically. However, from a numerical standpoint they are basically the same problem.
"""

# ╔═╡ 81c6042e-20a1-4aa7-b7ff-4a369126c22e
begin
	p = plot(exp, 0, 1, fill=0, fillalpha=0.3, xlabel=L"x", ylabel=L"e^x", ylims=(0, 2.7), label="")
	q = plot(x -> exp(sin(x)), 0, 1, fill=0, fillalpha=0.3, xlabel=L"x", ylabel=L"e^{\sin(x)}", ylims=(0, 2.7), label="")

	plot(p, q, layout=(2, 1))
end

# ╔═╡ 899f6285-331e-4cc5-8b0a-0e68c67082d2
md"""
## Trapezoidal rule
"""

# ╔═╡ 86ff3f84-f234-4de1-9d78-54a289c9319e
md"""
The task of **numerical integration** is to **approximate an integral** $\int_a^b f(x) \, dx$. We want to achieve this by **sampling the function** at $n$ carefully selected points $t_i$, $i=0, \ldots n$, followed by taking linear combinations of the results. We thus work towards approximations of the form
```math
\int_a^b f(x) \, dx ≈ \sum_{i=1}^m α_i f(t_i).
```
In general the $t_i$ can be distributed arbitrarily in the interval $[a, b]$.
However, for simplicity we will assume **equally spaced nodes** $t_i$
using the definition
```math
\tag{1}
t_i = a + i \, h, \qquad h = \frac{b-a}{n} \qquad \text{where $i = 0, 1, \ldots, n$}.
```
"""

# ╔═╡ d060f303-2dc6-4322-96eb-5a8deb3b7846
md"""
A first idea goes back to **polynomial interpolation**:
As **polynomials are easy to integrate analytically**,
we could just **fit a** $n$-th degree **polynomial**
through our $n+1$ nodes $t_0, t_1, \ldots, t_n$
and **then integrate that** instead of $f$ itself.

Since the integration of the polynomial is essentially exact,
the error of such a scheme is **dominated by the error of the polynomial
interpolation**.
Recall the [chapter on Interpolation](https://teaching.matmat.org/numerical-analysis/05_Interpolation.html), where we noted polynomials
through equispaced nodes to become numerically unstable and
possibly inaccurate for large $n$ due to Runge's phaenomenon.

Therefore we will pursue **piecewise linear polynomial interpolation**
instead of fitting an $n$-th degree polynomial.

The following graphics illustrates the idea:
"""

# ╔═╡ 2da93494-ceff-40ec-90d3-3b214b134b5e
let
	f(x) = 1 + (x + 1.5) * (x + 1.3) * (x - 1.5) * (x - 3)
	
	n = 6

	a = -2
	b = 1
	h = (b-a) / n
	ts = [a + i * h for i in 0:n]

	p = plot(f, -2, 1, lw=3, showaxis=false, grid=false, label=L"f(x)")

	it = 4
	plot!(p, [ts[it], ts[it+1], ts[it+1], ts[it]],
		     [0, 0, f(ts[it+1]), f(ts[it])], fill=0, c=cgrad(:grays)[0.8], label="")
	
	plot!(p, ts, f.(ts), c=:black, lw=3, mark=:o, label="piecewise linear")

	for i in 0:n
		plot!(p, [a + i * h, a + i * h], [0, f(a + i * h)];
			  ls=:dash, c=:black, label="", lw=1.5)
	end
	annotate!(p, [(a,   -0.5, L"a = t_0"),
				  (a+h, -0.5, L"t_1"),
				  (a+n*h, -0.5, L"t_n = b")])

	ylims!(p, (-0.9, 10.5))
	p
end

# ╔═╡ 03290f7e-46fd-4511-aa00-673a75de180f
md"""
To approximate integral $I = \int_a^b f(x) dx$, that is the area under the blue curve,
we evaluate the function $f$ at $n+1$ equispaced nodes
leading to data points $(t_i, f(t_i))$.
From these we construct a piecewise linear polynomial interpolation $p_{1,h}$ (black line).
Integrating $p_{1,h}$ on $[t_0, t_n]$ then yields an approximation to $I$:
```math
\begin{aligned}
I = \int_a^b f(x)\, dx &\approx \int_a^b p_{1,h}(x)\, dx = \int_a^b \sum_{i=0}^n f(t_i) H_i(x)\, dx
= \sum_{i=0}^n f(t_i) \int_a^b H_i(x)\, dx.
\end{aligned}
```
As a reminder the $H_i$ are the hat functions defined as
```math
    H_i(x) = \left\{ \begin{array}{ll}
    \frac{x - x_{i-1}}{x_i - x_{i-1}} & \text{if $i>1$ and $x\in [x_{i-1}, x_i]$}\\
    \frac{x_{i+1} - x}{x_{i+1} - x_{i}} & \text{if $i\leq n$ and $x\in [x_{i}, x_{i+1}]$}\\
    0 & \text{otherwise}
    \end{array}\right.
```
"""

# ╔═╡ 12235876-6af1-4091-b0ed-7e029b70992d
md"""
Recall that $t_{i+1} - t_{i} = h$, that is to say the distance between two nodes is always $h$.
With this and using the formulas for computing the areas of triangles
the integrals of the hat functions can be evaluated as
```math
\int_a^b H_i(x)\, dx 
= \left\{ \begin{array}{ll}
h & \text{if } i = 1, \ldots, n-1 \\
\frac h 2 & \text{if } i = 0 \text{ or } i = n
\end{array}\right.
```
leading to the final formula
```math
\begin{aligned}
I = \int_a^b f(x)\, dx \approx T_n(f) &:= \frac h2 f(t_0) + h f(t_1) + \cdots + h f(t_{n-1}) + \frac h2 f(t_n)\\
&= \frac h2 f(t_0) + \sum_{i=1}^{n-1} h\, f(t_i) + \frac h2 f(t_n).
\end{aligned}
```
This formula is called the **trapezoidal rule**.
This name stems from a more geometric way of constructing the trapezoidal formula:
as can be captured in the visualisation above, one can think of the formula as  summing the areas of the trapezoids defined by the quadrature nodes (see shaded grey).

An implementation of the trapezoidal rule is:
"""

# ╔═╡ fb68914d-fc61-4160-abe3-7f311412de64
function trapezoid(f, a, b, n)
	# f:  Function
	# [a, b]: Interval to integrate over
	# n:  Number of pieces to break the interval into
	h = (b - a) / n
	t = range(a, b, length=n+1)
	y = [f(tₙ) for tₙ in t]
	integral = h * (0.5y[1] + sum(y[2:n]) + 0.5y[n+1])
	(; integral, h)
end

# ╔═╡ a1d83cb2-6e0d-4a53-a11f-60dc020249d4
md"""
Recall that in Theorem 4 of [chapter 05 (Interpolation)](https://teaching.matmat.org/numerical-analysis/05_Interpolation.html) we found that
the piecewise polynomial interpolation shows quadratic convergence
```math
\|f - p_{1,h}\|_\infty \leq α h^2 \| f'' \|_\infty,
```
where $α > 0$.
With this in mind we can bound the error of the trapezoidal rule as
```math
\begin{aligned}
I - T_n(f) &= \int_a^b f(x)\,dx - \int_a^b p_{1,h}(x)\,dx\\
&= \int_a^b [f(x) - p_{1,h}(x)] \, dx\\
&\leq \|f-p_{1,h}\|_\infty\, \int_a^b\, dx\\
&\leq \underbrace{(b-a)\, α\, \|f''\|_\infty}_{=\widetilde{α}} \, h^2 = \widetilde{α} \, h^2 = O(h^2).
\end{aligned}
```
We thus expect **quadratic convergence** with $h$.

Let us confirm this numerically on the simple integral from earlier, i.e.
"""

# ╔═╡ 9c2919cb-8580-422e-82f3-cd70b9e15112
f(x) = exp(x)

# ╔═╡ 3c37067a-941f-4ec2-863d-741e0e76c694
md"""
and 
```math
I = \int_0^1 e^x\,dx = e - 1
```
which is approximately $(round(exact; digits=7)).

We consider a sequence of results where we double the number of integration points:
"""

# ╔═╡ f2381a2a-f940-4593-91e3-a37aef8127a0
ns = [2^i for i in 1:11]

# ╔═╡ 734ab4de-5340-4f4f-adfa-8dd0e8001e6d
md"""The corresponding approximate integrals using the trapezoidal formula $T_n(f)$ and quadrature node distances are:
"""

# ╔═╡ bc197c36-9294-49ce-af73-514c35860d81
begin
	Tfs    = Float64[]
	hs     = Float64[]
	errors = Float64[]
	for n in ns
		res = trapezoid(f, 0, 1, n)
		Tf   = res.integral  # Trapezoidal approximation
		h     = res.h        # Quadrature node distances
		error = abs(Tf - exact)
		@printf "n =%4d   T(f) = %.10f   error = %.10f\n" n Tf error
		
		push!(Tfs, Tf)
		push!(hs, h)
		push!(errors, error)
	end
end

# ╔═╡ 18506fc5-9019-4076-b77e-1c000d01b946
md"""
We notice that the error decreases roughly by a factor $4$ when we double the number of quadrature nodes (i.e. half the distance $h$ between the quadrature nodes). This again confirms the idea of a second-order convergence.

This becomes even clearer if we plot the error versus $h$ in a log-log plot along with a line of slope $2$:
"""

# ╔═╡ f189cf70-d72e-483c-af61-a3346ddd201d
let
	# Log-log plot where we flip the x-axis
	p = plot(hs, errors, yaxis=:log, xflip=true, xaxis=:log,
		     xlabel=L"h", ylabel="error", label="Trapezoidal", lw=2)

	# Plot guiding line with slope 2 (quadratic convergence)
	plot!(p, hs, hs.^2, ls=:dash, c=:black, label=L"$O(h^2)$; slope $2$")

	xticks!(p, 10.0 .^ (0:-0.5:-3))
	yticks!(p, 10.0 .^ (0:-1:-7))
end

# ╔═╡ 5f0178ad-c05b-4ac4-a008-b35d105cc7b2
md"""
The trapezoidal rule is just one representative of the wide class of numerical integration formulas. A general definition is:

!!! info "Definition: Numerical integration formula"
    A **numerical integration formula** for the $N+1$ equispaced **quadrature nodes**
	$t_i$, $i=0, \ldots, N$ is the set of **weights** $w_0, w_1, \ldots w_N$,
    such that for an integrand $f : [a, b] \to \mathbb{R}$
    ```math
	\tag{2}
    \int_a^b f(x) dx \approx Q_a^b(f) = \frac{b-a}{N}\, \sum_{i=0}^N w_i f(t_i).
	```
	The weights $w_i$ are independent of $f$.

	An older and still frequently used name for numerical integration is **quadrature**.

By comparing with our derivation above, we realise:

!!! info "Definition: Trapezoid formula"
	The trapezoid formula is the numerical integration of the form (2) with $n+1$ nodes (i.e. $N = n$) and weights
	```math
	w_i = \left\{ \begin{array}{ll} 1 & \text{if } 0 < i < n \\ \frac12& \text{if }i = 0 \text{ or } i = n \end{array}\right.
	```
"""

# ╔═╡ c31c7012-6986-441f-ae99-5e2bb2b469e5
md"""
## Simpson's rule
"""

# ╔═╡ 4fdd346f-a279-4686-b775-675cfe9f9f21
TODO("Show an illustrative drawing as well")

# ╔═╡ becdcc4e-ecff-46ad-8f3a-92d117e374f7
md"""
Considering the construction of the trapezoidal rule
we may easily wonder: why stop at using only linear polynomials
to proximate $f$ within each interval $[t_i, t_{i+1}]$ ?

Indeed, Simpson's formula takes the idea one step further and
constructs a quadratic polynomial within each interval $[t_i, t_{i+1}]$
by evaluate $f$ on $t_i$, $t_{i+1}$
as well as the midpoint $m_i = \frac{t_i + t_{i+1}}{2}$.
This overall leads to a **piecewise quadratic interpolant**
 $p_{2,h}$ of the integrand $f$ on $[a, b]$, which we again integrate exactly.
This is a little harder to compute and will be done as an exercise. The resulting formula is Simpson's formula
```math
\tag{3}
\begin{aligned}
\int_a^b f(x)\, dx &\approx  S_{2n}(f)  = \int_a^b p_{2,h}(x)\, dx \\
&= \sum_{i=1}^n \frac h6 \big( f(t_{i-1}) + 4f(m_{i-1}) + f(t_i) \big)\\
&= \frac h6 f(t_0) + \frac h3 \sum_{i=1}^{n-1} f(t_i) + \frac {2h}3 \sum_{i=0}^{n-1} f(m_i) + \frac h6 f(t_n).
\end{aligned}
```
"""

# ╔═╡ 5ee0ec01-54c3-48d8-8ba8-4460144002dd
md"""
While a little harder to see, this formula can also be brought into the form of (2):
it employs **$2n + 1$ equispaced nodes**
--- namely the collection of both the $t_i$ for $i=0, \ldots, n$ *and* the $m_i$ for $i=0,\ldots n-1$.
Therefore $N = 2n$ in (2) leading to a **nodal distance** of $\frac{b-a}{2n} = \frac{h}{2}$, where we used that $h = t_{i+1} - t_i = \frac{b-a}{n}$
"""

# ╔═╡ abcafa59-e8a1-4438-9ff6-3e8fc9fbd28d
md"""
!!! exercise
    Derive Simpson's rule, i.e. show that
    ```math
    \int_a^b p_{2,h}(x)\, dx = \sum_{i=1}^n \frac h6 \big( f(t_{i-1}) + 4f(m_{i-1}) + f(t_i) \big)
    ```
"""

# ╔═╡ 277ec8d1-949b-457a-8c1a-12d357a76efc
md"""
A Julia implementation of Simpson's rule is given below:
"""

# ╔═╡ 538b816c-5cc3-4e5c-b259-33825bdc39c3
function simpson(f, a, b, n)
	# f:  Function
	# [a, b]: Interval to integrate over
	# n:  Number of pieces to break the interval into
	h = (b - a) / n
	t = range(a, b, length=n+1)        # Subinterval boundaries
	m = range(a+h/2, b-h/2, length=n)  # Subinterval midpoints

	# Evaluate f and compute approximation
	ft = f.(t)
	fm = f.(m)
	integral = h * (ft[1]/6 + sum(ft[2:n])/3 + 2sum(fm)/3 + ft[n+1]/6)

	(; integral, h=h/2)  # Note h/2 since the actual nodal distance is half of h
end

# ╔═╡ a255138c-74bd-4aad-b49a-0a16746f3bda
md"""
Comparing the error of the trapezoidal and Simpson's quadratures against the exact integral of $\int e^x dx$ we obtain a much faster convergence for Simpson's rule, numerically looking like an $O(h^4)$ convergence.
"""

# ╔═╡ df5373cc-b997-4a80-906a-9bde0689c84b
p_convergence = let
	# Trapezoidal rule
	Tfs     = [trapezoid(f, 0, 1, n).integral for n in 1:3:800]
	hs_trap = [trapezoid(f, 0, 1, n).h        for n in 1:3:800]
	p = plot(hs_trap, abs.(Tfs .- exact);
	         yaxis=:log, xflip=true, xaxis=:log,
	         xlabel=L"h", ylabel="error", label="Trapezoidal", lw=2)

	# Simpson's rule
	Sfs      = [simpson(f, 0, 1, n).integral for n in 1:2:400]
	hs_simps = [simpson(f, 0, 1, n).h        for n in 1:2:400]
	plot!(p, hs_simps, abs.(Sfs .- exact), label="Simpson", lw=2)

	# Guiding lines
	plot!(p, hs_trap, hs_trap.^2, ls=:dash, c=1, label=L"$O(h^2)$; slope $2$")
	plot!(p, hs_trap, hs_trap.^4/1000, ls=:dash, c=2, label=L"$O(h^4)$; slope $4$")

	xticks!(p, 10.0 .^ (0:-0.5:-3))
	yticks!(p, 10.0 .^ (0:-2:-14))
end

# ╔═╡ 3aa18f97-3bbb-4fa3-9314-70900d0832c8
md"""
In line with the definition of algebraic convergence (and convergence order) for other approximation techniques we define the accuracy of a quadrature formula as:

!!! info "Definition: Convergence order of numerical integration"
	A numerical integration formula $Q_a^b(f)$ of the form (2)
	with equally spaced quadrature notes of separation $h$ is of **order $p$**
	if a constant $α>0$ indepentent of $h$
	(but possibly dependent on $f$) exists, such that
	```math
	\left| \int_a^b f(x)\,dx - Q_a^b(f) \right| \leq α\, h^p
	```
	as long as the function $f$ is sufficiently regular.
"""

# ╔═╡ c7f38ce7-41f2-4725-b5e7-68f20a8df7a5
md"""
We notice that our numerical investigation suggests:

!!! info "Convergence order of common numerical integration techniques"
	- **Trapezoidal rule:** Convergence order $p=2$
	- **Simpon's rule:** Convergence order $p=4$
"""

# ╔═╡ adf895be-b9d1-40e1-9c2b-146b30b996be
md"""
## Error analysis

In this lecture we only consider so-called **composite quadrature formulas**,
i.e. formulas which satisfy
```math
Q_a^b(f) = \sum_{i=0}^{N-1} Q_{t_i}^{t_{i+1}}(f)
```
where $Q_{t_i}^{t_{i+1}}(f)$ is the quadrature formula applied to the subinterval $[t_{i-1}, t_i]$, i.e. the application of the quadrature formula to $\int_{t_{i-1}}^{t_i} f(x) dx$ with no no subdivision of the integration interval.
As a consequence we can decompose the error
```math
\tag{4}
\int_a^b f(x)\,dx - Q_a^b(f)
= \sum_{i=0}^{N-1} \int_{t_i}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f)
```
into error contributions from each of the intervals $[t_{i-1}, t_i]$.

Assume for simplicity that the function $f$ is smooth and we can thus
build a Taylor expansion
```math
f(x) = \sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m_i) \, (x-m_i)^k
```
around the midpoint $m_i = \frac{t_{i+1} + t_i}{2}$ of the interval $[t_i, t_{i+1}]$.
Based on this we can deduce a series for the exact integral:
```math
\begin{aligned}
\int_{t_i}^{t_{i+1}} f(x) dx
&= \int_{t_i}^{t_{i+1}} \left[\sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m_i) (x-m_i)^k \right]dx
\\
&= \sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m_i) \int_{t_i}^{t_{i+1}} (x-m_i)^k dx \\
&= \sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m_i) \int_{t_i}^{t_{i+1}} q_k(x) \,dx,
\end{aligned}
```
where we defined $q_k(x) = (x-m_i)^k$.
"""

# ╔═╡ 4144570b-f8d9-49c8-af6b-732966864755
md"""
If we apply a quadrature formula (2) to $f$ we can follow through with
similar steps.
Here note that for $Q_{t_i}^{t_{i+1}}$ the number of nodes
is exactly two --- the beginning of the interval and the end,
such that $N=1$ and the term $\frac{b-a}{N} = \frac{t_{i+1} - t_i}{1} = h$.
We obtain:
```math
\begin{aligned}
Q_{t_i}^{t_{i+1}}(f) &= h\, \sum_{j=i}^{i+1} w_j f(t_j) \\
&= h\, \sum_{j=i}^{i+1} w_j  \left[\sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m_j) \, (t_j-m_j)^k \right]\\
&= h\, \sum_{j=i}^{i+1} w_j  \left[\sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m_j) \, q_k(t_j) \right]\\
&= \sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m) \left[
h\, \sum_{j=i}^{i+1} w_j \, q_k(t_j)
\right] \\
&= \sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m) Q_{t_i}^{t_{i+1}}(q_k)
\end{aligned}.
```
The difference between these expressions is exactly the error
contribution from the interval $[t_{i}, t_{i+1}]$, namely
```math
\tag{5}
\begin{aligned}
\int_{t_i}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f)
&= \sum_{k=0}^\infty \frac{1}{k!} f^{(k)}(m) \left[ \int_{t_i}^{t_{i+1}} q_k(x) \,dx - Q_{t_i}^{t_{i+1}}(q_k) \right].
\end{aligned}
```
"""

# ╔═╡ 7812f9d7-cda7-4d05-a800-3f7bea2e0e8c
md"""
The error of the integration formula can thus be completely understood
by studying the error $\int_{t_i}^{t_{i+1}} q_k(x) - Q_{t_i}^{t_{i+1}}(q_k)$.
Since $q_k$ is just a polynomial of degree $k$
we can thus understand the accuracy of quadrature formulas
*for arbitrary functions* by studying the accuracy of quadrature formulas
for *polynomials*,
which is considerably simpler.

One property of quadrature formulas is their **degree of exactness**:

!!! info "Definition: Degree of exactness"
	A numerical integration formula $Q_a^b(f) = \frac{b-a}{N}\, \sum_{i=1}^n w_i f(t_i)$ as given in (2)
	has a **degree of exactness $r$** if it integrates all
	monomials $x^s$ with $0 \leq s \leq r$ ($s$ integer) exactly, i.e. if
	```math
	\tag{6}
		Q_a^b(x^s) \textcolor{red}{=} \int_a^b x^s \,dx \qquad \forall s \in \mathbb{N} \text{ with } 0 \leq s \leq r,
	```
	but not for $s = r+1$.
"""

# ╔═╡ bc2043be-41e0-4083-9f8b-82b3ce6a13af
md"""
Note that the polynomial
```math
q_k(x) = ( x - m_i )^{k} = x^k + \left(\begin{smallmatrix}k\\1\end{smallmatrix}\right)\, x^{k-1} m_i + \left(\begin{smallmatrix}k\\2\end{smallmatrix}\right)\, x^{k-2} m_i^2
+ \cdots + \left(\begin{smallmatrix}k\\k-1\end{smallmatrix}\right)\, x \, m_i^{k-1}
+ m_i^k
```
only features monomials $x^s$ with $0 \leq s \leq k$.
Therefore a formula with degree of exactness $r$ will have
$\int_{t_i}^{t_{i+1}} q_k(x) \,dx - Q_{t_i}^{t_{i+1}}(q_k) = 0$ for $k \leq r$.
In (5) the first non-zero error term is thus
```math
\begin{aligned}
\left|\int_{t_i}^{t_{i+1}} q_{r+1}(x) \,dx - Q_{t_i}^{t_{i+1}}(q_{r+1})\right|
&\stackrel{(\ast)}{=} \left|\int_{t_i}^{t_{i+1}} x^{r+1} \,dx - Q_{t_i}^{t_{i+1}}(x^{r+1})\right| \\
&\stackrel{(\S)}{\leq}
\widetilde{α}_i h^{r+2}
\end{aligned}
```
where in $(\ast)$ all powers in $x$ less than $r+1$ drop again because of $Q$'s degree of exactness and in $(\S)$ we skipped a few non-trivial steps,
which are optional and will be presented below.
This is also the leading-order error term, such that
```math
\left|\int_{t_i}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f)\right| ≤ \widetilde{α}_i h^{r+2}
```
"""

# ╔═╡ 9cb44731-9d50-4718-9d4d-8fffca3a387f
md"""
The error in each of the the $N$ subintervals thus converges with $(r+2)$-th order,
such that combining with (4) and using the triangle inequality
we obtain the total error as
```math
\begin{aligned}
\left|\int_a^b f(x)\,dx - Q_a^b(f)\right|
&\leq \sum_{i=0}^{N-1} \left|\int_{t_i}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f)\right|\\
&\leq h^{r+2} \underbrace{\sum_{i=0}^{N-1} \widetilde{α}_i}_\text{$N$ terms}\\
&\leq h^{r+1} \frac{b-a}{N} N \, \max_i \widetilde{α}_i \\
&= C \, h^{r+1}
\end{aligned}
```
where $α = (b-a)\, \max_i \widetilde{α}_i$.

We notice:
"""

# ╔═╡ 4a048166-315d-4bf6-b099-866c9e3f8813
md"""
!!! info "Theorem 1 (Convergence order, simple version)"
	A numerical integration formula
	$Q_a^b(f) = \frac{b-a}{N}\, \sum_{i=1}^n w_i f(t_i)$ as given in (2)
	with **degree of exactness $r$ 
	is of order $r+1$**
	as long as the function $f$ is at least $r+1$ times differentiable.
"""

# ╔═╡ 9bad55de-c2d7-44f9-9478-7b66751bfbfe
Foldable("Optional: More details on Theorem 1",
md"""
!!! info "Theorem 1"
	Let $Q_a^b(f)$ be a numerical integration formula of form (2)
	with $N+1$ quadrature nodes of equal spacing $h = (b-a) / N$,
	which approximates $I = \int_a^b f(x)\,dx$.
	If $Q(f)$ has a degree of exactness $r$ then it is of order $r+1$
	as long as the function is $(r+1)$-times 
	continuously differentiable on the interval $[a, b]$.
	In particular we have
	```math
	\left| \int_a^b f(x)\,dx - Q(f) \right| \leq α \, \|f^{(r+1)}\|_\infty \, h^{r+1}
	```
	where the constant $α = \frac{b-a}{2^r \, (r+1)!}$ depends neither on $f$ nor on $h$.

		 
**Proof:**
As discussed in this lecture we consider only composite formulas, such that
```math
Q_a^b(f) = \sum_{i=0}^{N-1} Q_{t_i}^{t_{i+1}}(f)
```
and the error decomposition (4)
```math
\int_a^b f(x)\,dx - Q_a^b(f)
= \sum_{i=0}^{N-1} \int_{t_i}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f)
```
hold.

We first the error in each interval $[t_{i-1}, t_i]$
for which we perform a Taylor expansion of $f$ around the midpoint
$m_i = \frac{t_i + t_{i+1}}{2}$:
```math
\begin{aligned}
f(x) &= \underbrace{f(m_i) + f'(m_i) (x - m_i) + \cdots + 
\frac{f^{(r)}(m_i)}{r!} (x - m_i)^r}_{T^r_f(x)}\\
&+
\underbrace{\frac{f^{(r+1)}(\xi_i)}{(r+1)!} (x-m_i)^{r+1}}_{R^r_f(x)}
\end{aligned}
```
where $\xi_i \in (m_i, x)$ and we denoted by $T^r_f$
the $r$-th degree Taylor polynomial
and by $R^r_f$ the $r+1$-st degree Lagrange remainder term.

Since the quadrature formula has degree of exactness $r$ we have that
```math
Q_{t_i}^{t_{i+1}}(T_f^r) = \int_{t_{i}}^{t_{i+1}} T_f^r(x)\,dx.
```
Therefore
```math
\tag{P}
\begin{aligned}
\left| \int_{t_{i}}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f) \right|
&= 
\left| \int_{t_{i}}^{t_{i+1}} R^r_f(x)\,dx - Q_{t_i}^{t_{i+1}}(R^r_f) \right|\\
&\leq
\left| \int_{t_{i}}^{t_{i+1}} R^r_f(x)\,dx\right| + \left| Q_{t_i}^{t_{i+1}}(R^r_f) \right|\\
&\leq
\int_{t_{i}}^{t_{i+1}} \left| R^r_f(x)\right|\,dx + Q_{t_i}^{t_{i+1}}\left(\left| R^r_f \right|\right).
\end{aligned}
```
Now note that
```math
| R^r_f(x)| \leq 
\frac{\max_{x \in [t_{i}, t_{i+1}]} |f^{(r+1)}(x)|}{(r+1)!}
(h/2)^{r+1}
```
such that we can develop equation (P) further to 
```math
\begin{aligned}
\left| \int_{t_{i}}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f) \right| &\leq
\frac{(h/2)^{r+1}}{(r+1)!} \max_{x \in [t_{i}, t_{i+1}]} |f^{(r+1)}(x)|
\left( \int_{t_{i}}^{t_{i+1}} 1 \, dx + Q_{t_i}^{t_{i+1}}(1) \right) \\
&\leq
\frac{h^{r+1}}{2^{r+1} (r+1)!} \max_{x \in [t_{i}, t_{i+1}]} |f^{(r+1)}(x)|
\, 2h \\
&=
\frac{h^{r+1}}{2^r (r+1)!} \max_{x \in [t_{i}, t_{i+1}]} |f^{(r+1)}(x)| \, h
\end{aligned}
```
Therefore for a quadrature formula over the entire domain:
```math
\begin{aligned}
\left| \int_a^b f(x)\,dx - Q_a^b(f) \right|
&= \left| \sum_{i=0}^{N-1} \int_{t_{i}}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f) \right|\\
&\leq \sum_{i=0}^{N-1} \left| \int_{t_{i}}^{t_{i+1}} f(x)\,dx - Q_{t_i}^{t_{i+1}}(f) \right|\\
&\leq \frac{h^{r+1}}{2^r (r+1)!} \max_{x \in [a, b]} |f^{(r+1)}(x)| \, \sum_{i=0}^{N-1} h \\
&=\frac{h^{r+1}}{2^r (r+1)!} \left\|f^{(r+1)}\right\|_\infty \, (b-a),
\end{aligned}
```
which proves the result.

""")

# ╔═╡ a4cf020a-bb95-4fed-ad50-11ebae1934f5
md"""
In short by **determining the degree of exactness** of a quadrature formula
we thus **automatically obtain its convergence order**.

In fact in agreement with our expected convergence orders (2 and 4)
one can show that:

!!! info "Degree of exactness of common numerical integration techniques"
	- **Trapezoidal rule:** Degree of exactness $r=1$
	- **Simpson's rule:** Degree of exactness $r=3$

We will show this now explicitly for the case of Simpson's rule.
"""

# ╔═╡ 1c5167d5-0475-4522-bc5f-28df2f806d6b
md"""
!!! warning "Example: Simpson's formula has degree of exactness r=3."
	We will confirm that Simpson's formula has degree of exactness $r=3$.

	Recall that the idea of Simpson's formula was to split the integral
	into $n$ equally sized subintervals $[t_{i-1}, t_i]$ as
	```math
	\int_a^b f(x)\,dx = \sum_{i=1}^{n} \int_{t_{i-1}}^{t_i} f(x) \, dx
	```
	and then on such a subinterval approximate
	```math
	\tag{7}
	\int_{t_{i-1}}^{t_i} f(x) \, dx \approx
	\frac h6 \big( f(t_{i-1}) + 4f(m_{i-1}) + f(t_i) \big).
	```
	If the approximation in (7) is an equality on all subintervals
	for the monomials with degree $\leq 3$, then Simpson's formula
	is an equality for these monomials on the full interval $[a, b]$ as well.
	Inserting $f(x) = x^s$ into (7) we notice that our task is to show
	```math
	\int_{t_{i-1}}^{t_i} x^s \, dx \textcolor{blue}{\,\mathbf{=}\,} 
	\frac h6 \big( t^s_{i-1} + 4m^s_{i-1} + t^s_i \big)
	\qquad \text{for $0 \leq s \leq 3$},
	```
	but that this does not hold for $s=4$.

	An additional trick we can employ to simplify our calculations
	is to assume that the integration interval is
	$[t_{i-1}, t_i] = [-1, 1]$ with a length of
	$h=2$, $t_{i-1} = -1$, $t_i = 1$, $m_{i-1} = 0$.
	While perhaps suprising this does not actually
	change the generality of our result  [^1].
    Inserting the values for $h$, $t_{i-1}$, $t_i$ and $m_i$ into the above expression
    we need to show that $0 \leq s \leq 3$
	```math
	\int_{-1}^{1} x^s \, dx \textcolor{blue}{\,\mathbf{=}\,} 
	\frac 26 \big( (-1)^s + 4 \cdot 0^s + (+1)^s \big)
    = \frac{(-1)^s + 4\cdot0^s + (+1)^s}{3},
	```
	but that this equality fails for $s\geq 4$.

	We have
	```math
	\begin{aligned}
    \text{$s=0$:} && \int_{-1}^1 1\, dx &= 2 = \frac{(-1)^0 + 4\cdot0^0 + 1^0}{3} && \checkmark\\
	\text{$s=1$:} && \int_{-1}^1 x \, dx &= 0 = \frac{-1 + 0 + 1}{3} && \checkmark \\
	\text{$s=2$:} && \int_{-1}^1 x^2 \, dx &= \frac23 = \frac{(-1)^2 + 4 \cdot 0^2 + 1^2}{3}  && \checkmark \\
	\text{$s=3$:} && \int_{-1}^1 x^3 \, dx &= 0 =
	\frac{(-1)^3 + 4 \cdot 0^3 + 1^3}{3}   && \checkmark 
	\end{aligned}
	```
	but in contrast
	```math
	\begin{aligned}
	\text{$s=4$:} && \int_{-1}^1 x^4 \, dx = \frac25 &\textcolor{red}{\neq}
	\frac{(-1)^4 + 4 \cdot 0^4 + 1^4}{3} = \frac23
	\end{aligned}
	```
	Therefore Simpson's formula has **degree of exactness $r = 3$**.

	[^1]: The reason is that by an appropriate change of variables
		we can always rescale the integration coordinate,
		such that integration does run over $[-1, 1]$.
		Note that this in principle does generate some extra terms
		on top of $x^s$, but these are always of a lower power in $x$ than $x^s$ itself.
		So if we proceed as shown here, namely to consider integrating monomials with increasing powers of $x$, then all terms generated by this change of variables
		are already known to be integrated exactly at this stage.
"""

# ╔═╡ 145edec6-a9c9-4d36-99fb-cb4c914ef603
md"""## Extrapolation techniques

In numerical integration the **computationally expensive step** is usually the **evaluation of the function $f$** at the points of the quadrature nodes.
Generally higher-order numerical integration formulas are able to achieve higher accuracy with the same or less function evaluations.
We will discuss one systematic technique to obtain higher-order quadrature formulas
called **extrapolation**.
"""

# ╔═╡ 1080859f-c505-46ea-9d21-1a336325d8da
TODO("Make the general treatment optional and focus mostly directly applying this to the Trapezoidal rule.")

# ╔═╡ 635ae6b3-1770-431b-a799-0eefe20dfd0d
TODO("What is confusing here is that before we did expansions (I - T_n(f) = stuff, but now we do I = T_n(f) + stuff. )")

# ╔═╡ 9587ec30-95b1-4bcd-b8b3-18c0ce089e3c
md"""
Suppose a quantity $A_0$ is approximated by an algorithm $A_h$ with error expansion
```math
\tag{8}
A_0 = A_h + c_1 h + c_2 h^2 + c_3 h^3 + \cdots = A_h + c_1 h + c_2 h^2 + O(h^3).
```
Crucially it is not necessary to know the coefficients $c_i$, the only requirement is for them to be indepentent of $h$.
If we now look at the better approximation $A_{h/2}$ with half the node spacing then
```math
\tag{9}
A_0 = A_{h/2} + \frac{c_1}2 h + \frac{c_2}4 h^2 + O(h^3)
```
Forming a **linear combination** $2 \cdot \text{(9)} - \text{(8)}$
**between both methods** we obtain
```math
A_0 = \underbrace{2A_{h/2} - A_h}_{=\widetilde{A}_{h/2}} - \frac{c_2}{2} h^2 + O(h^3)
```
which **defines a new approximation algorithm**
$\widetilde{A}_{h/2} = 2A_{h/2} - A_h$.
Note that while $A_h$ is a first-order algorithm, $\tilde{A}_{h/2}$
is a **second-order algorithm**.

This idea to cancel the leading-order error term by taking a linear combination of two formulas of simple and half node spacing is termed **Richardson extrapolation**.
Importantly in contrast to simply reducing the node spacing
this scheme is able to **increase the convergence order**.
"""

# ╔═╡ 15b22584-5f21-4f00-ac04-0643fb1dfd56
md"""
Let us **apply this to the trapezoidal formula** for approximating the integral $I = \int_a^b f(x)\, dx$. We use $n+1$ quadrature nodes of equal separation $h = (b-a)/n$.
As we have discussed above the trapezoidal formula is of order $2$,
so the leading-order error term is $h^2$.
However, in this fortunate case one can even show
(using the [Euler–Maclaurin formula](https://en.wikipedia.org/wiki/Euler%E2%80%93Maclaurin_formula))
that the odd powers of $h$ are missing, i.e.
```math
\tag{10}
I = T_n(f) + c_2 h^2 + c_4 h^4 + O(h^6),
```
where
```math
c_2 = - \frac{1}{12} \left(f'(b) - f'(a)\right)
\qquad
c_4 = \frac{1}{740} \left(f'''(b) - f'''(a)\right).
```
While it could happen that $f'(b) = f'(a)$
we will make the general assumption that $c_2 \neq 0$.
"""

# ╔═╡ 70eddba0-82fd-4b5e-8d52-a961c03db5a6
md"""
For convenience of notation we rewrite (10)
in terms of the number of subintervals. Using $n = O(h^{-1})$ we obtain
```math
\tag{11}
I = T_n(f) + c_2\, n^{-2} + c_4\, n^{-4} + O(n^{-6}).
```
The Trapezoidal rule with twice the number of subintervals (half the spacing) similarly reads
```math
\tag{12}
I = T_{2n}(f) + \frac14 c_2\,n^{-2} + \frac{1}{16} c_4\, n^{-4} + O(n^{-6}).
```
Using an appropriate linear combination we can again cancel out the second-order term. Specifically, define
```math
\tag{13}
S_{2n}(f) = \frac{1}{3} \Big( 4\,T_{2n}(f) - T_n(f) \Big)
```
and consider $\frac43 \cdot \text{(12)} - \frac13 \text{(11)}$:
```math
\tag{14}
I = S_{2n}(f) - \frac{1}{4} c_4 \, n^{-4} + O(n^{-6}).
```
The quadrature formula $S_{2n}(f)$ is thus a 4-th order method.
While not immediately clear (14) turns out to be
identical to Simpson's formula (3), see the details at the end of this section.
In other words we can view **Simpson's formula** as
the numerical integration formula obtained either
- by integrating a piecewise quadratic approximation of the integrand or
- by applying Richardson extrapolation to the trapezoidal rule.
"""

# ╔═╡ 3be40834-84c6-45ed-9b20-0b054011dfc3
md"""
We formulate a general procedure for  Richardson extrapolation
in the context of numerical integration:

!!! info "Observation: General procedure for Richardson extrapolation"
	Given an error expansion of a quadrature $Q_n$ with $n$ subintervals
	```math
	\tag{15}
	I = Q_n + c_p n^{-p} + O(n^{-q}),
	```
	where $p$ is integer and $q > p$.
	Then a a **quadrature of at least order $q$**
	can be found by the linear combination
	```math
	\tag{16}
	\widetilde{Q}_{2n} = \frac{2^p Q_{2n} - Q_n}{2^p - 1}.
	```
This can be verified as follows:
using twice the number of nodes in (15) we obtain
```math
\tag{17}
I = Q_{2n} + 2^{-p} c_p n^{-p} + O(n^{-q}).
```
Then considering the linear combination
$\frac{2^p}{2^p-1} \text{(17)} - \frac{1}{2^p-1}\text{(15)}$
between the equations of single and
double the number of subintervals we obtain
```math
I = \widetilde{Q}_{2n} + O(n^{-q}),
```
i.e. that $\widetilde{Q}_{2n}$ is of order at least $q$.
"""

# ╔═╡ 7165d25c-7cd3-43f6-b697-623648d43258
md"""
Note that Simpson's formula (13) turns out to be another error expansion of the form (8), so **we can apply Richardson extrapolate another time**,
this time using Simpson's rule with $2n$ and $4n$ quadrature subintervals. 
From (16) we note the appropriate linear combination to be
```math
\tag{18}
R_{4n}(f) = \frac{1}{15}\Big( 16 S_{4n}(f) - S_{2n}(f) \Big),
```
which is a **6-th order numerical integration formula**.

This process can of course be repeated using $8n$, $16n$, $\ldots$ subintervals,
leading to higher and higher quadrature orders.
One usually calls this **Romberg integration**,
which we will not present in full generality here.
"""

# ╔═╡ b265c24f-9b37-48ed-9241-1a95b1902048
details("Optional: Details on the equality of the two forms of Simpson's formula",
md"""
We now want to understand why (13) is equal to the form (3) of Simpson's formula we obtained earlier. First recall the definition
```math
\begin{aligned}
T_n(f) &= \frac h2 f(t_0) + \sum_{i=1}^{n-1} h\, f(t_i) + \frac h2 f(t_n) \\
       &= \frac h2 f(a)   + \sum_{i=1}^{n-1} h\, f(a + i\, h) + \frac h2 f(b)
\end{aligned}
```
of the Trapezoid formula, where we used $h = \frac{b-a}{n}$ and $t_i = a + i \, h$ for $i = 0, \ldots, n$.
Similarly we have
```math
\begin{aligned}
T_{2n}(f) &= \frac{\tilde{h}}{2} f(a) + \sum_{j=1}^{2n-1} \tilde{h}\, f(a + j\, \tilde{h}) + \frac{\tilde{h}}{2} f(b) \\
&= \frac{\tilde{h}}{2} f(a) + \sum_{i=0}^{n-1} \tilde{h}\, f(a + (2i+1)\, \tilde{h})
+ \sum_{i=1}^{n-1} \tilde{h}\, f(a + 2i\, \tilde{h})
+ \frac{\tilde{h}}{2} f(b)
\end{aligned}
```
where $\tilde{h} = \frac{h}{2} = \frac{b-a}{2n}$
and in the last line we have split the sum
over $j$ into two sums, one for only the odd $j$ and one for only the even $j$.

First we note the even $j$ to be identical to the quadrature nodes
of $T_n$, i.e. $f(a + 2i\, \tilde{h}) = f(a + i\, h) = f(t_i)$
for $i = 0, \ldots, n$.
The remaining odd $j$ we identify to be equal to the midpoints
of our previous construction of Simpson's formula
```math
f(a + (2i+1)\, \tilde{h}) = f\left(\frac{a + (i+1)\,h}{2} + \frac{a + i\,h}{2}\right)
= f\left( \frac{t_i + t_{i+1}}{2} \right) = f(m_i)
```
for $i = 0, \ldots, n-1$.
With this we rewrite
```math
\tag{19}
T_{2n}(f) = \frac{h}{4} f(t_0) + \frac{h}{2} \sum_{i=0}^{n-1} f(m_i) + \frac{h}{2} \sum_{i=1}^{n-1}  f(t_i) + \frac{h}{4} f(t_n),
```
which we insert into (13) to obtain
```math
\begin{align*}
S_{2n}(f) &= \frac{1}{3} \Big( 4\,T_{2n}(f) - T_n(f) \Big) \\
&= 
\frac{h}{3} f(t_0) + \frac{2h}{3} \sum_{i=0}^{n-1} f(m_i) + \frac{2h}{3} \sum_{i=1}^{n-1} f(t_i) + \frac{h}{3} f(t_n)\\
&- \frac h6 f(t_0) \hspace{7.6em} - \frac{h}{3}\sum_{i=1}^{n-1} f(t_i) - \frac h6 f(t_n)\\
&= \frac h6 f(t_0) + \frac {2h}3 \sum_{i=0}^{n-1} f(m_i) + \frac h3 \sum_{i=1}^{n-1} f(t_i) + \frac h6 f(t_n).
\end{align*}
```
This is again the expression (3) we had previously for Simpson's formula.
""")

# ╔═╡ 5534dfa0-f5fb-47d7-9b43-b0f97adf76be
md"""
### Node doubling

If we consider the integration formula $R_{4n}$ we notice that
$R_{4n}$ depends on $S_{4n}$ and $S_{2n}$, which in turn depend on
$T_n$, $T_{2n}$ and $T_{4n}$.
Continuing along the **Romberg integration procedure**
we would **recursively perform more and more levels of extrapolation**.
At the next level we also require $T_{8n}$, then $T_{16n}$ and so on,
such that **in each step** the **number of nodes** to consider roughly **doubles**.

This **node doubling is particularly advantageous in practice** and much
preferred over other schemes.
The reason is illustrated in the figure below,
which divides the node spacing by half from one level to the next.
In each layer new nodes are only introduces at the midpoints (shown in red).
Notably, about half of the total number of nodes in each layer is red.
Therefore moving from $h$ to $h/2$ (respectively from $n+1$ to $2n+1$ nodes)
**only about half of the nodes are new** and need to be evaluated.
On the others the function values are already known (as they were needed at the prevous level) and can be re-used without extra cost.
"""

# ╔═╡ 57e9518b-b8d8-4d4e-911c-9034afd807f6
RobustLocalResource("https://github.com/epfl-matmat/numerical-analysis/blob/30a76271c3fa3be9d24896ecd2450d18fa505df3/src/img/node-doubling.svg", "img/node-doubling.svg", :width => 600)

# ╔═╡ ea69e253-ff02-45b8-a760-4467260f36b4
TODO("Using the drawing above this can be better explained, in particular T_{2n} easily re-written from T_n if one considers halfing the interval size ... perhaps some extra labels need to be introduced in the drawing above")

# ╔═╡ d5885f27-7171-4a78-a32f-ead4a4791b0b
md"""
To make this explicit we use a formulation for the trapezoidal rule, which is developed in the part *Details on the equality of the two forms of Simpson's formula*, equation (19).
As it states the trapezoidal rule with $2n$ subintervals can be expressed
in terms of quantities of the trapezoidal rule with $n$ subintervals,
that is in terms of the quadrature nodes $t_i = a + i \, h$ (for $i = 0, \ldots, n$),
the nodal spacing $h = (b-a) / n$ and the midpoints $m_i = \frac{t_i + t_{i+1}}{2}$:
```math
T_{2n}(f) = \frac{h}{4} f(t_0) + \frac{h}{2} \sum_{i=0}^{n-1} f(m_i) + \frac{h}{2} \sum_{i=1}^{n-1}  f(t_i) + \frac{h}{4} f(t_n).
```
Noting $m_i = a + (2i+1)\, \frac{h}{2}$
for $i=0, \ldots, n-1$
and using the definition of $T_n(f)$ this can be reformulated as
```math
\tag{20}
T_{2n}(f) = \frac12 T_{n}(f) + \frac{h}{2} \sum_{i=0}^{n-1} f\left(a + (2i+1)\, \frac{h}{2}\right).
```
Since $a + (2i+1)\, \frac{h}{2} = a + (2i+1)\, \frac{b-a}{2n}$
are just the *odd* quadrature nodes of the $T_{2n}(f)$
quadrature formula, we see that only these $n$ odd nodes
need to be evaluated additionally when
considering $T_{2n}(f)$ instead of $T_n(f)$.
"""

# ╔═╡ 18d1f510-0395-4ce9-8ae8-087bd11468a5
md"""
Let's see this in an example. We want to compute
```math
V = \int_0^2 x^2 e^{-2x} \, dx
```
using extrapolation. First we use `quadgk` to get an accurate value:
"""

# ╔═╡ 64bf6442-2b6d-4052-8756-6190f7570038
begin
	g(x) = x^2 * exp(-2x)
	a = 0
	b = 2
	V, _ = quadgk(g, a, b, atol=1e-14, rtol=1e-14)
end;

# ╔═╡ 54c4c55b-b235-4215-a6ce-7df732826a56
V

# ╔═╡ 0753e01f-4f35-40dc-a094-b989c573d20f
md"Based on the trapezoidal rule on $n = 20$ subintervals we get a first estimate:"

# ╔═╡ 5b70fa56-ad42-4efb-956f-8aa4d38c2152
T20 = let
	n = 20
	h = (b - a) / n
	t = h * (0:n)
	y = g.(t)
	T20 = h*y[1]/2 + h*sum(y[2:n]) + h*y[n+1]/2
end

# ╔═╡ c9f99d3f-a5bb-4fb9-b0b7-51260ca5ffa9
md"Now we double to $n = 40$, but we only need to evaluate the odd nodes:"

# ╔═╡ 992d3f7e-8de6-4346-900f-04b8ac0d8b9d
T40 = let
	n = 40
	h = (b - a) / n
	t = h * (0:n)
	y_odd = g.(t[2:2:n])
	T40 = T20/2 + h*sum(y_odd)
end

# ╔═╡ dbccd71f-c6bd-404a-a6e3-3b90d688f98a
md"We repeat a second time and double $n$ again:"

# ╔═╡ a436e787-4b2f-4ee2-b138-5c65f6a54854
T80 = let
	n = 80
	h = (b - a) / n
	t = h * (0:n)
	y_odd = g.(t[2:2:n])
	T80 = T40/2 + h*sum(y_odd)
end

# ╔═╡ 73fdd6be-567c-4523-ac15-23c210a81794
md"Using (13) we can perform a first level of extrapolation and get the two Simpson values $S_{40}(f)$ and $S_{80}(f)$:"

# ╔═╡ 94bc43ca-06f6-4957-87ee-26380db97800
S40 = (4T40 - T20) / 3

# ╔═╡ e51b187e-cb35-4d7b-b2c4-1aaaecdd1787
S80 = (4T80 - T40) / 3

# ╔═╡ 4456d0c6-5db8-4167-8a6f-e163934491a4
md"Finally, we perform one more level of extrapolation to get the sxth-order accurate result $R_{80}(f)$:"

# ╔═╡ 865e86a3-1d65-4e54-be16-988d2a0125a9
R80 = (16S80 - S40) / 15

# ╔═╡ cf6e7935-268e-458a-990a-39432ab3e9f3
md"We compute all errors to 10 digits:"

# ╔═╡ cc54215f-3af2-48c0-83ec-04cc786fd2ce
begin
	eT = round.(V .- [T20 T40 T80]; digits=10)
	eS = round.(V .- [S40 S80];     digits=10)
	eR = round.(V .- [R80];         digits=10)
end;

# ╔═╡ 88c53514-4d5a-4328-a8eb-967d8f031b37
md"""
and summarise them in a table along with the number of function evaluations:

  number of evals |   order 2    | order 4      | order 6
  --------------- | ------------ | ------------ | ---------------
  21              |    $(eT[1])  |              | 
  41              |    $(eT[2])  | $(eS[1])     |
  81              |    $(eT[3])  | $(eS[2])     | $(eR[1])

We notice that the 6th order **result obtained using Richardson extrapolation**
is about **twice as accurate** as the order 2 result,
even though it uses the **same number of function evaluations**.

Since the cost of function evaluation is usually dominating
in numerical integration,
the take-away is that one should **never just use low-order quadratures**,
but always **employ extrapolation techniques.**
"""

# ╔═╡ 6c21de1a-d9fa-479d-b47c-b6304bf87d16
md"""
## Optional: *A posteriori* error estimation

If we want to numerically evaluate an integral $I= \int_a^b f(x)\, dx$
using a chosen quadrature formula $Q_n$, a natural question
to ask is how many quadrature nodes $n$ are necessary
to obtain an approximation of the the integral $I$
within a predefined error tolerance $ϵ$.

In the previous sections we discussed,
that Richardson extrapolation provides a receipe to obtain
a higher-order quadrature $\widetilde{Q}_{2n}(f)$
from the numerical integration values $Q_n(f)$ and $Q_{2n}(f)$.
As a result $\widetilde{Q}_{2n}(f)$ is in general
more accurate than $Q_{2n}(f)$, which motivates to employ the difference
```math
η_{2n} = |\widetilde{Q}_{2n}(f) - Q_{2n}(f)|
```
as an estimate of the error committed by the formula $Q_{2n}(f)$.
Indeed assuming an error expansion (15),
i.e. $I = Q_n(f) + c_p n^{-p} + O(n^{-q})$,
and and a corresponding Richardson extrapolation (16)
one verifies
```math
\begin{aligned}
|I - Q_{2n}(f)| &\leq |I - \widetilde{Q}_{2n}(f)| + |\widetilde{Q}_{2n}(f) - Q_{2n}(f)|
= \underbrace{η_{2n}}_{=O(n^{-p})} + \underbrace{|I - \widetilde{Q}_{2n}(f)|}_{=O(n^{-q})}\\
&= η_{2n} + O(n^{-q})
\end{aligned}
```
Therefore a simple adaptive strategy is to start with
$2n$ nodes by computing and $Q_n(f)$ and $Q_{2n}(f)$.
Then check whether $η_{2n} ≤ ϵ$ is satisfied
and if this is not the case keep doubling the number of nodes until it is.
Sticking to a single layer of extrapolation for simplicity
we obtain the following algorithm:
"""

# ╔═╡ 622b4050-b6ca-46a2-aed3-de1acb4b7748
md"""
!!! info "Algorithm: Adaptive quadrature based on Richardson extrapolation"
	Given the problem to compute $\int_a^b f(x)\,dx$,
	quadrature formula $Q_n(f)$ an initial number of subintervals $n_0$
	and a requested tolerance $ϵ$, compute:
	- Initialise $n = n_0$, $h = (b-a) / n_0$
	- Initial estimate of integral $I_n = Q_n(f)$
	- While $η_n > ϵ$ iterate:
	  * Update $h = \frac{b-a}{2n}$
	  *   $I_{2n} = Q_{2n}(f)$  $\qquad$ (using updated node distance $h$)
	  * Extrapolate $\tilde{I}_{2n} = \frac{2^p I_{2n} - I_n}{2^p - 1}$ using (14)
	  *   $η_{2n} = |\tilde{I}_{2n} - I_{2n}|$
	  * Update $n = 2n$

Applied to the trapezoidal formula this is implemented as follows:
"""

# ╔═╡ bce0faed-d46b-4bc1-872d-62a494387b58
function trapezoid_adaptive(f, a, b; n₀=2, tol=1e-12)
	# f:  Function
	# [a, b]: Interval to integrate over
	# n₀: Initial number of subintervals (i.e. n₀+1 quadrature nodes)
	# tol: Desired tolerance
	n = n₀
	h = (b - a) / n₀  # Node separation
	t = (0:n) .* h   # Quadrature nodes
	y = f.(t)        # Evaluate function on all nodes
	Tₙ = h * (0.5y[1] + sum(y[2:n]) + 0.5y[n+1])  # Trapezoidal formula with n₀ nodes

	extrapolated = 0.0
	ηₙ = 10tol
	while ηₙ > tol
		h = (b - a) / 2n  # New node separation
		t = (0:2n) .* h   # New set of quadrature nodes
		y_odd = f.(t[2:2:2n])  # Evaluate only at odd nodes

		# Use node doubling formula (20) to evaluate Trapezoidal formula with 2n nodes.
		T₂ₙ = Tₙ/2 + h * sum(y_odd)

		# Perform extrapolation: Note that p = 2 for the Trapezoidal formula
		extrapolated = (4T₂ₙ - Tₙ) / 3

		n  = 2n
		Tₙ = T₂ₙ
		ηₙ = abs(extrapolated - T₂ₙ)
	end
	
	(; integral=extrapolated, h, n)
end

# ╔═╡ b7c4c712-c1f6-48c7-98e8-f7d21a04c043
md"We test this works as expected using our example from the node doubling section:"

# ╔═╡ 42f7241a-7bea-4071-b45a-0f3b5c283863
(; integral, h, n) = trapezoid_adaptive(g, 0, 2; tol=1e-5)

# ╔═╡ 850c093f-ac5c-4175-b7d6-b8713718ecd5
md"The error is indeed below the requested tolerance:"

# ╔═╡ 571c5934-d5b4-4e7b-84f9-72437e8c794f
abs(integral - V)

# ╔═╡ b41adad4-2cf5-44ca-a1b8-d70ab1ab97ef
md"""
The next step to improve this approach would be to additionally employ
Romberg integration, i.e. to recursive extrapolate not only
from $T_{2n}(f)$ & $T_{4n}(f)$ to $S_{4n}(f)$,
but to also employ $T_n(f)$ and evaluate
the 6-th order formula $R_{4n}(f)$
like we did in the numerical example of the node doubling section.
However, these ideas are out of scope for us.
"""

# ╔═╡ 7d6e26e5-96ee-447d-86e4-27588726e429
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 320)
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
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.1"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.57"
QuadGK = "~2.11.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "0b137bedfe8a154be757099565fbd14062a17cca"

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
# ╟─d34833b7-f375-40f7-a7a6-ab925d736320
# ╠═dc68915e-d230-11ee-09c8-4fce8a858118
# ╟─47114a9b-0e74-4e48-bb39-b49f526f1e9b
# ╟─85c29f83-4617-41b7-90e4-f863c3f09b1d
# ╠═10e9343f-4f0b-4d95-9c5a-3e7a6f6b9794
# ╟─9a300e03-8ec2-41ed-85c2-b7c23c140010
# ╠═53347e5e-00b1-4317-8912-331e90c3aa07
# ╟─30b6c88e-ebf8-49b4-acbb-9dacdb650b28
# ╠═f95ba764-ac44-45e0-a714-3fffe2066c7f
# ╟─bfccad43-b336-4271-83eb-40af1ed36849
# ╠═81c6042e-20a1-4aa7-b7ff-4a369126c22e
# ╟─899f6285-331e-4cc5-8b0a-0e68c67082d2
# ╟─86ff3f84-f234-4de1-9d78-54a289c9319e
# ╟─d060f303-2dc6-4322-96eb-5a8deb3b7846
# ╟─2da93494-ceff-40ec-90d3-3b214b134b5e
# ╟─03290f7e-46fd-4511-aa00-673a75de180f
# ╟─12235876-6af1-4091-b0ed-7e029b70992d
# ╠═fb68914d-fc61-4160-abe3-7f311412de64
# ╟─a1d83cb2-6e0d-4a53-a11f-60dc020249d4
# ╠═9c2919cb-8580-422e-82f3-cd70b9e15112
# ╟─3c37067a-941f-4ec2-863d-741e0e76c694
# ╠═f2381a2a-f940-4593-91e3-a37aef8127a0
# ╟─734ab4de-5340-4f4f-adfa-8dd0e8001e6d
# ╠═bc197c36-9294-49ce-af73-514c35860d81
# ╟─18506fc5-9019-4076-b77e-1c000d01b946
# ╟─f189cf70-d72e-483c-af61-a3346ddd201d
# ╟─5f0178ad-c05b-4ac4-a008-b35d105cc7b2
# ╟─c31c7012-6986-441f-ae99-5e2bb2b469e5
# ╠═4fdd346f-a279-4686-b775-675cfe9f9f21
# ╟─becdcc4e-ecff-46ad-8f3a-92d117e374f7
# ╟─5ee0ec01-54c3-48d8-8ba8-4460144002dd
# ╟─abcafa59-e8a1-4438-9ff6-3e8fc9fbd28d
# ╟─277ec8d1-949b-457a-8c1a-12d357a76efc
# ╠═538b816c-5cc3-4e5c-b259-33825bdc39c3
# ╟─a255138c-74bd-4aad-b49a-0a16746f3bda
# ╠═df5373cc-b997-4a80-906a-9bde0689c84b
# ╟─3aa18f97-3bbb-4fa3-9314-70900d0832c8
# ╟─c7f38ce7-41f2-4725-b5e7-68f20a8df7a5
# ╟─adf895be-b9d1-40e1-9c2b-146b30b996be
# ╟─4144570b-f8d9-49c8-af6b-732966864755
# ╟─7812f9d7-cda7-4d05-a800-3f7bea2e0e8c
# ╟─bc2043be-41e0-4083-9f8b-82b3ce6a13af
# ╟─9cb44731-9d50-4718-9d4d-8fffca3a387f
# ╟─4a048166-315d-4bf6-b099-866c9e3f8813
# ╟─9bad55de-c2d7-44f9-9478-7b66751bfbfe
# ╟─a4cf020a-bb95-4fed-ad50-11ebae1934f5
# ╟─1c5167d5-0475-4522-bc5f-28df2f806d6b
# ╟─145edec6-a9c9-4d36-99fb-cb4c914ef603
# ╠═1080859f-c505-46ea-9d21-1a336325d8da
# ╠═635ae6b3-1770-431b-a799-0eefe20dfd0d
# ╟─9587ec30-95b1-4bcd-b8b3-18c0ce089e3c
# ╟─15b22584-5f21-4f00-ac04-0643fb1dfd56
# ╟─70eddba0-82fd-4b5e-8d52-a961c03db5a6
# ╟─3be40834-84c6-45ed-9b20-0b054011dfc3
# ╟─7165d25c-7cd3-43f6-b697-623648d43258
# ╟─b265c24f-9b37-48ed-9241-1a95b1902048
# ╟─5534dfa0-f5fb-47d7-9b43-b0f97adf76be
# ╟─57e9518b-b8d8-4d4e-911c-9034afd807f6
# ╠═ea69e253-ff02-45b8-a760-4467260f36b4
# ╟─d5885f27-7171-4a78-a32f-ead4a4791b0b
# ╟─18d1f510-0395-4ce9-8ae8-087bd11468a5
# ╠═64bf6442-2b6d-4052-8756-6190f7570038
# ╠═54c4c55b-b235-4215-a6ce-7df732826a56
# ╟─0753e01f-4f35-40dc-a094-b989c573d20f
# ╠═5b70fa56-ad42-4efb-956f-8aa4d38c2152
# ╟─c9f99d3f-a5bb-4fb9-b0b7-51260ca5ffa9
# ╠═992d3f7e-8de6-4346-900f-04b8ac0d8b9d
# ╟─dbccd71f-c6bd-404a-a6e3-3b90d688f98a
# ╠═a436e787-4b2f-4ee2-b138-5c65f6a54854
# ╟─73fdd6be-567c-4523-ac15-23c210a81794
# ╠═94bc43ca-06f6-4957-87ee-26380db97800
# ╠═e51b187e-cb35-4d7b-b2c4-1aaaecdd1787
# ╟─4456d0c6-5db8-4167-8a6f-e163934491a4
# ╠═865e86a3-1d65-4e54-be16-988d2a0125a9
# ╟─cf6e7935-268e-458a-990a-39432ab3e9f3
# ╠═cc54215f-3af2-48c0-83ec-04cc786fd2ce
# ╟─88c53514-4d5a-4328-a8eb-967d8f031b37
# ╟─6c21de1a-d9fa-479d-b47c-b6304bf87d16
# ╟─622b4050-b6ca-46a2-aed3-de1acb4b7748
# ╠═bce0faed-d46b-4bc1-872d-62a494387b58
# ╟─b7c4c712-c1f6-48c7-98e8-f7d21a04c043
# ╠═42f7241a-7bea-4071-b45a-0f3b5c283863
# ╟─850c093f-ac5c-4175-b7d6-b8713718ecd5
# ╠═571c5934-d5b4-4e7b-84f9-72437e8c794f
# ╟─b41adad4-2cf5-44ca-a1b8-d70ab1ab97ef
# ╟─7d6e26e5-96ee-447d-86e4-27588726e429
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
