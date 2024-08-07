### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 679517fb-6090-4259-9b93-571ded17ab88
begin
	using Plots
	using Polynomials
	using PlutoUI
	using PlutoTeachingTools
	using LaTeXStrings
	using LinearAlgebra
	using HypertextLiteral
end

# ╔═╡ 61e5ef66-a213-4b23-9406-9cc63a58104c
TableOfContents()

# ╔═╡ 8c00d7cc-c050-11ee-23f8-e90afbb1a5e9
md"""
# Interpolation

In this chapter we will return to the first problem we already briefly discussed in the introductory lecture, namely:

!!! info "Definition: Interpolation problem"
    Suppose we are given data $(x_i, y_i)$ with $i = 1, 2, 3, \ldots n$,
    where the $x_i$ are all distinct. Find a function $f$
    such that $f(x_i) = y_i$.

The function $f$ is usually called the **interpolant** or **interpolating function**. This immediately leads to the question how to choose $f$.
In practice there is often not the *one* correct choice,
much rather the answer depends on the properties of our data
and the properties we want the final $f$ to have.
E.g. frequently we want $f$ to be continuous or differentiable
or easy to integrate, because we want to do some post-processing on the obtained function.
Moreover in practical applications the quality of the fit usually depends
on the function $f$ to agree with the behaviour of the data itself.
For example, if the data is periodic (like an oscillatory system)
we should choose an $f$ to be periodic as well.

In general our ansatz will be to take $f$ from a family of suitable functions
(e.g. polynomials, trigonometric functions etc),
which form a vector space.
If $φ_1, φ_2, \ldots φ_n$ is a basis for for this vector space,
then we can write $f$ as the linear combination
```math
\tag{1}
f(x) = \sum_{j=1}^n c_j φ_j(x)
```
of $n$ basis functions. Employing further the condition
that $f$ needs to pass through the data points, i.e. $f(x_i) = y_i$ for all $i$,
then we get $n$ equations
```math
f(x_i) = \sum_{j=1}^n c_j φ_j(x_i) = y_i \qquad i = 1, \ldots n
```
which is a system of $n$ equations with $n$ unknowns.
In matrix form we can write
```math
\tag{2}
\textbf{A} \textbf{c} = \textbf{y}
```
where
```math
\textbf{A} = \left(\begin{array}{ccc}
φ_1(x_1) & \ldots & φ_n(x_1) \\
φ_1(x_2) & \ldots & φ_n(x_2) \\
\vdots \\
φ_1(x_n) & \ldots & φ_n(x_n) \\
\end{array}\right),
\qquad 
\textbf{c} = \left(\begin{array}{c} c_1\\c_2\\\vdots\\c_n\end{array}\right),
\qquad
\textbf{y} = \left(\begin{array}{c} y_1\\y_2\\\vdots\\y_n\end{array}\right)
```

*Note:* The number of basis functions and the number of data points does not need to agree. We will consider such more general regression problems later.
"""

# ╔═╡ eb43ef87-343e-4a9b-9a0f-3a18364d6b4b
md"""
## Polynomial interpolation

Let us come back to our earlier problem of finding an interpolating function $f$
for the given temperature-dependent reaction rates.

Note, that compared to the pathological example of the introduction, we start with a slightly simpler case:
"""

# ╔═╡ 9e933ebe-1805-46fa-9e0a-a64c03a71c5c
data = """
# Temperature(K)  Rate(1/s)
  250.0           1.65657
  260.0           1.70327
  270.0           1.74472
  280.0           1.78110
  290.0           1.81259
  300.0           1.83940
  310.0           1.86171
  320.0           1.87971
  330.0           1.89358
  340.0           1.90352
  350.0           1.90968
""";

# ╔═╡ f9f4ddec-fd3d-4624-ace4-f316b9cfa1e4
begin
	lines = split(data, "\n")
	temperature = [parse(Float64, split(line)[1]) for line in lines[2:end-1]]
	rate        = [parse(Float64, split(line)[2]) for line in lines[2:end-1]]
end;

# ╔═╡ 8be732e2-927b-4c56-8cab-9079c7598e86
scatter(temperature, rate; label="data", xlabel="temperature", ylabel="rate")

# ╔═╡ b73ea6d2-5751-4926-a28d-418deb60da2a
md"""
### Monomial basis

Let us first consider the case of **polynomial interpolation**:
given $n+1$ data points $(x_1, y_1)$ up to $(x_{n+1}, y_{n+1})$
we want to find a $n$-th degree polynomial
```math
\tag{3}
p_n(x) = \sum_{j=0}^n c_j x^j,
```
which is an interpolating function,
i.e. satisfies $p_n(x_i) = y_i$
for all $i = 1, \ldots, n + 1$.
A result from algebra ensures that such a polynomial of degree $n$,
which goes through all $n+1$ data points can always be found.
Moreover this polynomical (thus its coefficients $c_j$)
is uniquely defined by the data. 

To find this polynomial a
natural idea is thus to employ the monomials $1, x, x^2, \ldots x^n$
directly as the basis for our interpolation:

!!! info "Definition: Monomial basis"
    The basis of $n+1$ monomials are the functions
    $φ_1(x) = 1, φ_2(x) = x, \ldots, φ_{n+1}(x) = x^n$.

Employing these in equations (1) and (2) to perform our interpolation
leads to the linear system in the $n+1$ unknowns $c_0, c_1, \ldots, c_n$
```math
\tag{4}
\underbrace{
\left(\begin{array}{ccc}
1 & x_1 & \ldots & x_1^n \\
1 & x_2 & \ldots & x_2^n \\
\vdots \\
1 &x_n & \ldots & x_n^n \\
\end{array}\right)}_{= \mathbf{V}}
\,
\underbrace{
\left(\begin{array}{c} c_1\\c_2\\\vdots\\c_n\end{array}\right)
}_\textbf{c}
= 
\underbrace{
\left(\begin{array}{c} y_1\\y_2\\\vdots\\y_n\end{array}\right)
}_\textbf{y},
```
where the $(n+1)\times(n+1)$ matrix $\mathbf{V}$ is called the **Vandermonde matrix**. Assuming the data points $(x_i, y_i)$ to be distinct, we know that only one interpolating polynomial exists. The linear system (4) has therefore exactly one solution and the matrix $V$ is thus always invertible, $\det(\textbf{V}) \neq 0$.
"""

# ╔═╡ 2a088353-0b78-45a8-b354-6ced35b1c7f9
md"First we build the Vandermonde matrix:"

# ╔═╡ 4c1412f2-c1ac-4f95-afec-2c3ef1b09dd2
md"Then we solve the linear system to find the polynomial coefficients:"

# ╔═╡ 317058ec-7690-4eea-8fed-b9cd8190cc7b
md"... and plot the result:"

# ╔═╡ 2580774f-b4be-4665-9123-bc7c22210010
md"""
Use the Slider to regulate how many of the data points are included in the interpolation.
- `n_data_monomial =` $(@bind n_data_monomial Slider(1:11; show_value=true, default=3))

We see as we use more and more data points, eventually the full data range is well-represented by interpolating polynomial.
"""

# ╔═╡ 2dfef7e4-23ca-49fa-bb95-f9c7219cdf22
md"""
Let's see how this method performs on our data.
We demonstrate the case for `n_data_monomial = ` $(n_data_monomial),
thus a polynomial degree of $(n_data_monomial - 1)
"""

# ╔═╡ ac94e6f2-ebf0-4ee8-b21f-d9d20c656dd5
begin
	x = temperature[1:n_data_monomial]
	y = rate[1:n_data_monomial]

	V = zeros(n_data_monomial, n_data_monomial)
	for j in 1:n_data_monomial
		for i in 1:n_data_monomial
			V[i, j] = x[i] ^ (j-1)
		end
	end
	V
end

# ╔═╡ 5a1e03f2-07fc-4a31-9715-3b44d12824a2
c = V \ y  # Solve V * c = y

# ╔═╡ f7e3bfc6-1d12-4390-add7-12d0b0f8f4e1
let
	# Use a smart function from Julia to evaluate the polynomial.
	# Essentially computes (3)
	f(x) = evalpoly(x, c)

	p = scatter(temperature, rate; ylims=(1.65, 1.95),
	            label="all data", xlabel="temperature", ylabel="rate")
	scatter!(p, x, y; label="data used in fit")
	plot!(p, f; label=L"Polynomial interpolant $p$", lw=2)
end

# ╔═╡ d6c4e000-c213-40e0-8a6a-e1dfccadb006
md"""
**A word of warning:** Vandermonde matrices are an example of what is known as a *badly conditioned matrix*, that is a matrix where small numerical errors (such as rounding errors due to the finite-precision number format used by the computer) can amplify and lead to very inaccurate solutions. Therefore this method starts to be come very unreliable for $n$ larger than a few tens.

We will explore a better method in the next section.
"""

# ╔═╡ f6fdfe7b-fa92-4168-bc8b-1f3a5a30f8db
md"""
### Lagrange basis

While the monomial basis seemed natural, the fact that the $n$-th degree
interpolating polynomial is *uniquely* determined by the $n+1$ data points
allows us to use *any polynomial basis* to find it.

Both a more practical as well as a numerically more stable way to find
**the** interpolating polynomial is the approach using Lagrange ploynomials.

!!! info "Definition: Lagrange basis"
    The Lagrange polynomials associated to the **nodes**
    $x_1, x_2, \ldots x_{n+1}$ are the $n+1$ polynomials
    ```math
    \tag{5}
    \begin{aligned}
    L_{\textcolor{red}{i}}(x) &= \prod_{\stackrel{j=1}{\textcolor{red}{j\neq i}}}^{n+1} \frac{x-x_j}{\textcolor{red}{x_i} - x_j} \\
           &= \frac{(x-x_1)(x - x_2) \cdots (x - x_{i-1}) (x - x_{i+1}) \cdots (x-x_{n+1})}{
	(\textcolor{red}{x_i} - x_1)
	(\textcolor{red}{x_i} - x_2)
	\cdots
	(\textcolor{red}{x_i} - x_{i-1})
	(\textcolor{red}{x_i} - x_{i+1})
	\cdots
	(\textcolor{red}{x_i} - x_{n})
		   }
    \end{aligned}
    ```
    for $i = 1, \ldots, n+1$.

Each of these polynomials is of degree $n$, moreover they satisfy the **cardinality condition**:
```math
\tag{6}
L_j(x_i) = δ_{ji} = \left\{\begin{array}{ll} 1 &\text{if $i=j$}\\ 0 &\text{if $i\neq j$} \end{array}\right.
```
i.e. they are $1$ only on the nodal point with the same index, but $0$ on all other nodal points.

Visualisation of Lagrange polynomials
on an equally spaced grid between $0$ and $1$:

- Number of nodal points: $(@bind n_lagrange_vis Slider(2:15, default=4, show_value=true))
"""

# ╔═╡ 463fcdfa-89b0-4f9d-aa22-29d522585627
let
	nodes = range(0, 1; length=n_lagrange_vis)
	function L(i, x)
		r = one(x)
		for j in 1:length(nodes)
			i == j && continue
			r = r * (x - nodes[j]) / (nodes[i] - nodes[j])
		end
		r
	end

	if length(nodes) < 6
		examples = 1:length(nodes)
	elseif length(nodes) < 8
		examples = 2:2:length(nodes)
	else
		examples = 2:3:length(nodes)
	end

	p = plot(; ylims=(-2, 2), legend=:bottomright)
 	vline!(p, nodes; c=:grey, label="nodes", ls=:dash)
	
	x = range(0, 1; length=300)
	for (ic, i) in enumerate(examples)
		plot!(p, x, L.(i, x); label=LaTeXString("\$L_{$i}\$"), c=ic+1, lw=2)
	end
	for (ic, i) in enumerate(examples)
		plot!(p, [nodes[i], nodes[i]], [0.0, 1.0];
		      c=ic+1, label="", lw=3, ls=:dash, mark=:4)
	end
	
	hline!(p, [0], c=:grey, label="", ls=:dash)
	p
end

# ╔═╡ 8bfdb2e5-b3bf-4a17-a1cb-161eec581df8
md"""
In particular the nodal property (6) makes it extremely convenient to use a Lagrange polynomial basis to find the interpolating polynomial $p_n$:

!!! note "Proposition 1"
    The $n$-th degree
    polynomial $p_n$ interpolating the data $(x_i, y_i)$ for $i = 1, \ldots, n+1$
    is given by
    ```math
    \tag{7}
    p_n(x) = \sum_{i=1}^{n+1} y_i L_i(x)
    ```

> **Proof:**
> This can be easily verified:
> - Since every linear combination of an $n$-th degree polynomial is itself
>   an $n$-th degree polynomial, $p_n$ (as a linear combination of the $n$-th
>   degree Lagrange polynomials) is a $n$-th degree polynomial.
> - Moreover
>   ```math
>   p_n(x_i) = y_1 \underbrace{L_1(x_i)}_{=0} + y_2 \underbrace{L_2(x_i)}_{=0} + \cdots +  y_i \underbrace{L_i(x_i)}_{=1} + \cdots + y_{n+1} \underbrace{L_{n+1}(x_i)}_{=0} = y_i
>   ```
>   which confirms that $p_n$ is the interpolating polynomial.
"""

# ╔═╡ 2b0638bd-2ab1-49ff-921b-a86c4e49c002
md"""

!!! warning "Example: Using Lagrange polynomials"
    Using Lagrange basis functions, find the interpolating polynomial through the
    points $(x_1 = -1, y_1 = -6)$, $(x_2 = 1, y_2 = 0)$, $(x_3 = 2, y_3 = 6)$:
    - Using equation (7) we obtain
      ```math
      p_2(x) = y_1 L_1(x) + y_2 L_2(x) + y_3 L_3(x) = -6 L_1(x) + 6 L_3(x)
      ```
      where
      ```math
      \begin{aligned}
      L_1(x) &= \frac{(x-x_2)(x-x_3)}{(x_1-x_2)(x_1-x_3)} = \frac{(x-1)(x-2)}{(-1-1)(-1-2)} = \frac{(x-1)(x-2)}{6} \\
      L_3(x) &= \frac{(x-x_1)(x-x_2)}{(x_3-x_1)(x_3-x_2)} = \frac{(x+1)(x-1)}{(2+1)(2-1)} = \frac{(x+1)(x-1)}{3}
      \end{aligned}
      ```
      such that
      ```math
      p_2(x) = -6\frac{(x-1)(x-2)}{6} + 6\frac{(x+1)(x-1)}{3} = -4 + 3x + x^2.
      ```
"""

# ╔═╡ 5ecbe1d5-2250-4bab-8c1e-c319d23a0f4b
md"
For reference: In Julia a convenient way to interpolate a polynomial to a given set of points is provided by the `Polynomials` package, e.g."

# ╔═╡ f40c7c2e-339c-4b6e-91d3-ef5d031091f3
let
	poly = Polynomials.fit(temperature, rate)
	p = scatter(temperature, rate; ylims=(1.65, 1.95), xlims=(245, 355),
	            label="data", xlabel="temperature", ylabel="rate")
	plot!(p, poly; label=L"Polynomial interpolant $p$", lw=2, xlims=(245, 355))
end

# ╔═╡ dada1011-4a09-440b-9b92-ac167c8028fe
md"""
### Error analysis

With the Lagrange polynomials we have a simple approach to find an interpolating polynomial. From a numerical analysis perspective this raises the question **how good a polynomial approximation is**.
For our use case, where we are given $n$ data points $(x_i, y_i)$ this question cannot be answered as such, since in principle any of the (infinitely many) ways one could find an interpolant is equally "good".

So instead we ask another question: Assume we have a true function $f$ and we are allowed to observe $n$ samples $(x_n, f(x_n))$ from it, **would the polynomial approximation converge** to $f$ as we observe more and more samples of $f$. That is we want to check whether
```math
p_n \to f \quad \text{as} \quad n\to \infty.
```
Provided this convergence is indeed the case, we are usually also interested in the convergence rate
--- similar to our discussion about root finding algorithms.
Clearly, the faster $p_n$ converges $f$ the better, since we need less samples
to obtain an accurate reconstruction of the original $f$.

For illustration we contrast two cases, namely the construction of a polynomial interpolation of the functions
"""

# ╔═╡ b1773d1c-4d3c-4f4b-805d-f188623c9571
begin
	fsin(x)   = sin(5x)
	fratio(x) = 1 / (1 + 20x^2)
end

# ╔═╡ e567287b-9244-4b55-9f5b-2e3bf583e46a
md"""
each defined on [-1, 1].

Change the number of samples using the slider and observe the convergence:

- `n_samples_comparison = ` $(@bind n_samples_comparison Slider(2:20; default=10, show_value=true))
"""

# ╔═╡ 31126101-685d-464d-be4a-6233d3803780
let
	samples = range(-1, 1; length=n_samples_comparison)
	fine = range(-1.3, 1.3; length=500)
	
	values_fsin = fsin.(samples)
	poly_fsin   = Polynomials.fit(samples, values_fsin)
	
	values_fratio = fratio.(samples)
	poly_fratio   = Polynomials.fit(samples, values_fratio)
	
	p = plot(fine, poly_fsin.(fine);
	         label="Polynomial fit", xlims=(-1.3, 1.3), lw=2, ylims=(-1.2, 1.2),
	         legend=:topright, title=L"f_{\textrm{sin}}(x) = \sin(5x)",
	         titlefontsize=12)
	plot!(p, fine, fsin.(fine), label="fsin", lw=2, ls=:dash)
	scatter!(p, samples, values_fsin, label="samples")

	q = plot(fine, poly_fratio.(fine);
	         label="Polynomial fit", xlims=(-1.3, 1.3), lw=2, ylims=(-0.5, 2),
	         legend=:topright, title=L"f_{\textrm{ratio}}(x) = 1 / (1 + 20x^2)",
	         titlefontsize=12)
	plot!(q, fratio, label="fratio", lw=2,ls=:dash)
	scatter!(q, samples, values_fratio, label="samples")
	
	plot(p, q)
end

# ╔═╡ af534a92-31d4-4125-bb7c-3eb2b3430852
md"""
We see that the polynomial interpolation of `fsin` (left) converges very well to the original function. Moreover this convergence is very quickly, since already for `n_samples_comparsion = 10` hardly any difference is visible. On the other hand `fratio`, i.e the rather innocent looking function
```math
f_\text{ratio}(x) = \frac{1}{1 + 20x^2}
```
converges very poorly. In particular from about `n_samples_comparsion = 8`
spurious oscillations start to appear towards the end of the domain $[-1, 1]$,
which become more and more pronounced as we *increase* the number of samples
and the polynomial degree.
"""

# ╔═╡ b92ee89f-9241-4b4d-a54f-dafe4900e159
md"""
To understand this behaviour the following error estimate is useful:

!!! info "Theorem 2"
    For a $C^{n+1}$ function $f : [a, b] \to \mathbb{R}$
    and $a = x_1 < x_2 < \cdots < x_{n+1} = b$ *equally distributed*
    nodes in $[a, b]$ the $n$-th degree polynomial interpolant $p_n$
    of the data $(x_i, f(x_i))$ with $i = 1, 2, \ldots, n+1$
    satisfies the estimate
    ```math
    \tag{8}
    \left\|f - p_n\right\|_\infty
    \leq \frac{1}{4(n+1)}\left(\frac{b-a}{n}\right)^{n+1}
    \left\| f^{(n+1)} \right\|_\infty
    ```
    where the infinity norm $\| \phi \|_\infty$ for a function $\phi : D \to \mathbb{R}$
    mapping from a domain $D$ is the expression
    ```math
    \| \phi \|_\infty = \max_{x \in D} |\phi(x)|,
    ```
    i.e. the maximal absolute value the function takes.

In other words the error $\|f - p_n\|_\infty$ only decreases as $n$ increases
when the right-hand side (RHS) of (8) is bounded by a constant. So let's check this for our functions.
- For $f_\text{sin}(x) = \sin(5x)$ we can easily verify
  $|f_\text{sin}^{(n+1)}(x)| = 5 |f_\text{sin}^{(n)}(x)|$ as well as $\max_{[-1,1]} |f_\text{sin}(x)| = \max_{x\in[-1,1]} |\sin(5x)| = 1$, such that
  ```math
  \|f_\text{sin}^{(n+1)}\|_\infty =\max_{x\in[-1,1]} |f_\text{sin}^{(n+1)}(x)| = 5^{n+1}
  ```
  and (8) becomes (using $b = 1$ and $a = -1$):
  ```math
  \left\|f_\text{sin} - p_n \right\|_\infty \leq \frac{10^{n+1}}{4(n+1)n^{n+1}} \longrightarrow 0 \qquad \text{as $n \to \infty$}
  ```
- In contrast for $f_\text{ratio}(x)$ one can show
  ```math
  \begin{aligned}
  \|f_\text{ratio}^{(n+1)}\|_\infty
  &\sim 20^n n! \qquad \text{as $n\to\infty$}
  \end{aligned}
  ```
  and as a result convergence is not guaranteed.

Looking more closely at the error of the problematic case,
that is plotting $|f_\text{ratio}(x) - p_n(x)|$ as a function of $x$,
we observe an interesting pattern:
"""

# ╔═╡ de28f64c-f243-4cb5-a1e6-ae8562cce5bc
let
	fine = range(-1.3, 1.3; length=3000)

	p = plot(yaxis=:log, title=L"Error $|f_\textrm{ratio}(x) - p_n(x)|$ for small degree")
	for n_degree in (4, 8, 12)
		samples = range(-1, 1; length=n_degree+1)
		values_fratio = fratio.(samples)
		poly_fratio   = Polynomials.fit(samples, values_fratio)
		error(x) = abs(poly_fratio(x) - fratio(BigFloat(x)))
		plot!(p, fine, error.(fine), label="error n = $n_degree", lw=2)
	end
	p
end

# ╔═╡ 08d5cb8a-3bf5-4d40-a8d6-d22318e37712
md"""
For small polynomial degrees, the error decreases homogeneously as we increase $n$. Additionally we observe some almost vertical "drops", where the error goes down to machine precision (about $10^{-16}$). This is because at the nodal points --- by construction --- the polynomial is exact and only the floating-point, error remains.
"""

# ╔═╡ c5b33503-973b-468f-95de-acd6a7605a99
let
	fine = range(-1.1, 1.1; length=3000)

	p = plot(yaxis=:log, title=L"Error $|f_\textrm{ratio}(x) - p_n(x)|$ for large degree")
	for n_degree in (24, 34, 50)
		samples = range(-1, 1; length=n_degree+1)
		values_fratio = fratio.(samples)
		poly_fratio   = Polynomials.fit(samples, values_fratio)
		error(x) = abs(poly_fratio(x) - fratio(BigFloat(x)))
		plot!(p, fine, error.(fine), label="error n = $n_degree", lw=2)
	end
	p
end

# ╔═╡ 0a74569f-fd4b-4705-8350-3b47db65f502
md"""
Switching to higher degrees we observe the error to again *increase* at near the boundaries of the interval $[-1, 1]$. However in the central part of the interpolation domain $[-0.5, 0.5]$ the error remains constantly small.
This visual result is also confirmed by a more detailed analysis,
which reveals, that the origin is our choice of a *regular* spacing between the sampling points, an effect known as Runge's phaenomenon.

!!! info "Observation: Runge's phaenomenon"
    Polynomial interpolation employing equally spaced nodal points
    to construct the interpolant $p_n$ may lead to a non-convergence
    as $n \to \infty$. Moreover while for small $n$ the error in the
    infinity norm $\|\,\cdot\,\|_\infty$ still
    decreases, for large $n$ this behaviour can change, such that
    for large $n$ the error *keeps increasing* as $n$ increases.

The solution to this dilemma is to employ **irregularly spaced**
nodal points, in particular the nodal points need to become
*more densely spaced* towards the boundary of the domain.
One especially important node family
are the **Chebyshev extreme points** defined by
```math
x_k = - \cos\left(\frac{k \pi}{n}\right) \qquad \text{for $k = 0, 1, \ldots n$}.
```
Using these to interpolate a degree $n$ polynomial gives a uniform convergence behaviour as $n\to \infty$:
"""

# ╔═╡ c38b9e48-98bb-4b9c-acc4-7375bbd39ade
let
	fine = range(-1.1, 1.1; length=3000)
	p = plot(yaxis=:log, title="Error for Chebyshev interpolants")
	for n_degree in (5, 10, 25, 35, 50, 70)
		t = BigFloat[-cos(π*k/n_degree) for k = 0:n_degree]
		values_fratio = fratio.(t)
		poly_fratio   = Polynomials.fit(t, values_fratio)
		error(x) = abs(poly_fratio(x) - fratio(BigFloat(x)))
		plot!(p, fine, error.(fine), label="degree n = $n_degree", lw=2)
	end
	p
end

# ╔═╡ 479a234e-1ce6-456d-903a-048bbb3de65a
md"""

Notably chebyshev notes enjoy the following convergence result:

!!! info "Theorem 3"
    Let  $f : [-1, 1] \to \mathbb{R}$ be a function, which is
    analytic in an open interval containing $[-1, 1]$,
    that is the Talyor series of $f$ converges to
    $f(x)$ for any $x$ from this open interval.
    Then we can find constants $C > 0$ and $K > 0$
    such that
    ```math
	\| f - p_n \|_\infty = \max_{x\in[-1,1]} \left|f(x) - p_n(x)\right| \leq C K^{-n}
    ```
    where $p_n$ is the unique polynomial of degree $n$ defined by
    interpolation on $n+1$ Chebyshev points.

This is an example of **exponential convergence**: The error of the approximation scheme reduces by a *constant factor* whenever the polynomial degree $n$ is increased by a constant increment. 

The **graphical characterisation** is similar to the iterative schemes we discussed in the previous chapter: We employ a **semilog plot** (using a linear scale for $n$ and a logarithmic scale for the error), where exponential convergence is characterised by a straight line:
"""

# ╔═╡ d4cf71ef-576d-4900-9608-475dbd4d933a
let
	fine = range(-1.0, 1.0; length=3000)
	degrees = 2:2:50
	linf_errors_chebyshev = map(degrees) do n_degree
		t = BigFloat[-cos(π*k/n_degree) for k = 0:n_degree]
		poly_fratio = Polynomials.fit(t, fratio.(t))
		linf_error  = maximum(x -> abs(poly_fratio(x) - fratio(BigFloat(x))), fine)
	end
	linf_errors_equispaced = map(degrees) do n_degree
		t = range(-1.0, 1.0; length=n_degree+1)
		poly_fratio = Polynomials.fit(t, fratio.(t))
		linf_error  = maximum(x -> abs(poly_fratio(x) - fratio(BigFloat(x))), fine)
	end
	
	plot(degrees, linf_errors_chebyshev; yaxis=:log, lw=2, mark=:x, xlabel=L"Polynomial degree $n$", ylabel=L"\Vert f - p_n \Vert_\infty", label="Chebyshev nodes", title="Convergence of polynomial interpolation", ylims=(-Inf, 10))
	plot!(degrees, linf_errors_equispaced; lw=2, mark=:x, label="equally spaced nodes", title="Convergence of polynomial interpolation")
end

# ╔═╡ 56685887-7866-446c-acdb-2c20bd11d4cd
md"""
When designing approximation schemes, **obtaining exponential convergence**
is one of the **desired properties**.

!!! info "Observations"
    - For **exponential convergence**, the error reduces by a constant factor
    - A straight line is obtained when looking at the error norm on a `log`-scale
	  versus an appropriate accuracy parameter (such as the polynomial degree $n$, the spacing of the interpolating notes etc.)

!!! danger "Potential confusion: Convergence terminology"
    When discussing convergences rates of iterative numerical algorithms and
    the accuracy of numerical approximation schemes (interpolation, differentiation, integration, discretisation) unfortunately a different terminology is employed. In the following let $C > 0$ and $K > 0$ denote appropriate constants.
	- **Iterative schemes:** Linear convergence
	  * If the error scales as $C K^{-n}$ where $n$ is the iteration number, we say the scheme has **linear convergence**. (Compare to the last chapter.)
    - **Approximation schemes:** Exponential convergence
	  * If the error scales as $C K^{-n}$ where $n$ is some accuracy parameter (with larger $n$ giving more accurate results), then we say the scheme has **exponential convergence**.
      * For approximation schemes it is often more convenient to instead formulate the method using an approximation parameter $h = O(1/n)$, i.e. which scales inversely to $n$. In this case *smaller* $h$ give more accurate results. In this case exponential convergence is characterised by an error scaling as $C K^h$.
"""

# ╔═╡ 1d61d925-1ceb-439b-aff6-c5e8483fdc31
md"""
### Stability of polynomial interpolation

In the previous discussion we identified the Chebyshev nodal points to provide exponential convergence by avoiding Runge's phaenomenon. Another common question to ask in numerical analysis is referred to as **numerical stability**. The goal of stability analysis is to quantify by how much small perturbations in the input data translate to the obtained result of a numerical algorithm.

We consider the case of interpolating a polynomial $p_n$
to $n+1$ data points $(x_i, y_i)$ for $i = 1, 2, \ldots, n+1$,
where $a = x_1 < x_2 < \cdots < x_{n+1} = b$
and the $y_i$ are drawn from a function $f$, i.e. $y_i = f(x_i)$.
Our goal is effectively to recover $f$ as close as possible.
However, in many practical settings we don't have access to the
true data  $(x_i, y_i)$, since we are only able to obtain data
through a noisy measurement.
This means that the data we are *actually* able to collect
is much rather $\tilde{y}_i = f(x_i) + \varepsilon_i$
where $\varepsilon_i$ represents the measurement error.
We suppose further that $|\varepsilon_i| \leq \varepsilon$ for all $i= 1, 2, \ldots, n+1$ where $\varepsilon > 0$ is a small number,
i.e. that overall the errors are small.
Instead of producing an interpolating polynomial $p_n$
using the exact data $(x_i, y_i)$,
our procedure only has access to the noisy data  $(x_i, \tilde{y}_i)$,
thus producing the polynomial $\tilde{p}_n$.

In **stabilty analysis** we now ask the question:
How different are $p_n$ and $\tilde{p}_n$ given a measurement noise
of order $\varepsilon$.

Let us investigate this using the Lagrange basis, where
```math
p_n(x) = \sum_{i=1}^{n+1} y_i L_i(x)
\qquad \qquad
\tilde{p}_n(x) = \sum_{i=1}^{n+1} (y_i + \varepsilon_i) L_i(x)
```
therefore for all $x \in [a, b]$:
```math
\left| \tilde{p}_n(x) - p_n(x) \right|
= \left|\sum_{i=1}^{n+1} \varepsilon_i L_i(x)\right|
\leq \sum_{i=1}^{n+1}\left|\varepsilon_i L_i(x)\right|
\leq \varepsilon \sum_{i=1}^{n+1}\left| L_i(x)\right|
```
Introducing

!!! info "Definition: Lebesgue's constant"
    Given $n$ nodes $a = x_1 < x_2 < \cdots < x_{n+1} = b$
    from the interval $[a, b]$ we denote Lebesgue's constant
    as the quantity
    ```math
    \Lambda_n = \sup_{x\in[a,b]} \sum_{i=1}^{n+1} |L_i(x)|.
    ```

we can rewrite this as
```math
\tag{9}
\left| \tilde{p}_n(x) - p_n(x) \right| \leq \Lambda_n ε.
```
"""

# ╔═╡ 9766ff2e-2060-4508-8018-8aa6d51219bc
md"""
Notably, if $\Lambda_n$ is small, then small measurement errors $ε$ can only lead to small perturbations in the interpolating polynomial. In that case our polynomial interpolation procedure would be called **stable** or **well-conditioned**. By contrast, if $\Lambda_n$ is very high, then already a small measurement error $ε$ allows for notably deviations in the resulting interpolant
--- we are faced with a **badly conditioned** problem.

In general for numerical problems, we call the factor relating the error in the output quantity (here $\left| \tilde{p}_n(x) - p_n(x) \right|$)
to the error in the input quantity
(here $ε = \max_i |\tilde{x}_i - x_i|$) the **condition number**
of the problem. For polynomial interpolation the condition number is exactly Lebesgue's constant $\Lambda_n$.

Considering the two polynomial interpolation schemes we discussed, one can show
- Equally distributed nodes: $\Lambda_n \sim \frac{2^{n+1}}{e\, n \log n}$ as $n\to \infty$
- Chebyshev nodes: $\frac{2}{\pi} \log(n+1)+a < \Lambda_n < \frac{2}{\pi} \log(n+1) + 1, \qquad a = 0.9625\ldots$

Therefore for Chebyshev nodes the condition number grows only logarithmically with $n$, while for equally spaced nodes it grows *exponentially* !

Thus Chebyshev nodes do not only lead to faster-converging polynomial interpolations,
but also to notably more stable answers. As a result they are one of the standard ingredients in many numerical algorithms.

Let us verify this result visually. We consider again our function $f_\text{sin}(x) = \sin(5x)$ in the interval $[-1, 1]$, which we evaluate at the
distinct nodes $\{x_i\}_{i=1}^{\texttt{n\_nodes\_poly}}$
--- either equally spaced (left plot) or using Chebyshev nodes (right plot). Additional for both cases we consider an exact evaluation, i.e. points $(x_i, y_i) = (x_i, f(x_i))$ as well as noisy evaluations $(x_i, \tilde{y}_i)$ with
```math
\tilde{y}_i = f(x_i) + η_i
```
where $η_i$ is a random number of magnitude $|\eta_i| ≤ 10^\texttt{η\_log\_poly}$.
"""

# ╔═╡ d5de3b11-7781-4100-8a63-d26426685bbc
md"""
- Number of nodes `n_nodes_poly = ` $(@bind n_nodes_poly Slider(2:20; show_value=true, default=10))
- Logarithm of noise amplitude `η_log_poly = ` $(@bind η_log_poly Slider(-3:0.25:0; show_value=true, default=-2))
"""

# ╔═╡ 4a4540a3-b4cf-47fd-a615-a5280505333f
let
	fine = range(-1.1, 1.1; length=500)
	
	n_degree = n_nodes_poly - 1
	values_equispaced = range(-1, 1; length=n_nodes_poly)
	values_chebyshev  = [-cos(π*k/n_degree) for k in 0:n_degree]

	p = plot(; layout=(1, 2))
	for (i, values) in enumerate((values_equispaced, values_chebyshev))
		noise = [10^η_log_poly * (-1 + 2rand()) for _ in values]
		
		fsin_accurate = fsin.(values) #  = [fsin(x) for x in values]
		fsin_noise    = fsin.(values) + noise
		
		poly_accurate = Polynomials.fit(values, fsin_accurate)
		poly_noise    = Polynomials.fit(values, fsin_noise)

		plot!(p, fine, fsin.(fine);
		      label=L"function $f$", subplot=i, lw=2, ls=:dash, ylim=(-2, 2), legend=:topright)
		plot!(p, fine, poly_accurate.(fine); label=L"accurate $p_n$", subplot=i, lw=2)
		plot!(p, fine, poly_noise.(fine); label=L"noisy $\tilde{p}_n$", subplot=i, lw=2)
		scatter!(p, values, fsin_noise; subplot=i, label="noisy data")
	end
	title!(p, "equispaced"; subplot=1)
	title!(p, "Chebyshev"; subplot=2)

	p
end

# ╔═╡ 50a3face-5831-4ef5-b6a5-225b8b7c46a0
md"""
As we increase polynomial order and noise, we see larger discrepancies for the interpolation based on equispaced points then for the Chebyshev points.
"""

# ╔═╡ 2cefad7d-d734-4342-8039-aefdc33c2edd
md"""
## Piecewise linear interpolation

In the previous section we looked at $n$-th degree polynomial interpolation,
which we found to be poorly conditioned, e.g. when equispaced nodes and high 
polynomial degree are employed.
An alternative construction is to employ piecewise polynomials,
i.e. to interpolate a separate polynomial betwen each data point.
The most simple approach in this regard is *linear interpolation*,
i.e. just connecting the dots of each data point by a straight line.

!!! note "Definition: Piecewise linear interpolation"
    Given nodes $x_1 < x_2 < \cdots < x_{n+1}$
    and associated data $(x_i, y_i)$ for $i = 1, \ldots, n+1$
    the **piecewise linear interpolant**
    $p_{1,h}$ is given by
    ```math
    p_{1,h}(x) = y_i + \frac{y_{i+1} - y_i}{x_{i+1} - x_i} (x - x_i)
    \qquad \text{for $x\in [x_i, x_{i+1}]$}
    ```

Instead of using this definition to implement piecewise polynomial interpolation,
a more practical approach is to follow the idea of equation (1)
and construct an appropriate set of basis functions
for piecewise polynomial interpolation.
These are the

!!! note "Definition: Hat function"
    Given nodes $x_1 < x_2 < \cdots < x_{n+1}$
    their associated hat functions are the functions
    ```math
    \tag{10}
    H_i(x) = \left\{ \begin{array}{ll}
    \frac{x - x_{i-1}}{x_i - x_{i-1}} & \text{if $i>1$ and $x\in [x_{i-1}, x_i]$}\\
    \frac{x_{i+1} - x}{x_{i+1} - x_{i}} & \text{if $i\leq n$ and $x\in [x_{i}, x_{i+1}]$}\\
    0 & \text{otherwise}
    \end{array}\right.
    ```
    for $i = 1, \ldots, n+1$.

In code these may be obtained as:
"""

# ╔═╡ 611541de-927a-44fd-bdf6-bedec2befa01
function hatfun(nodes, i)
	# Function to generate the i-th hat function corresponding to nodal
	# points given in the vector nodes. It is assumed that nodes is sorted
	n = length(nodes) - 1

	# Define an inner function, which is returned
	return function (x)
		if i > 1 && nodes[i-1] ≤ x ≤ nodes[i]
			return (x - nodes[i-1]) / (nodes[i] - nodes[i-1])
		elseif i ≤ n && nodes[i] ≤ x ≤ nodes[i+1]
			return (nodes[i+1] - x) / (nodes[i+1] - nodes[i])
		else
			return 0
		end
	end
end

# ╔═╡ 1af1365e-de88-4140-861d-ca527ae67db7
md"""We plot the four hat functions for the example nodal points
$x_1 = 0$, $x_2 = 0.45$, $x_3 = 0.7$, $x_4 = 1.0$:"""

# ╔═╡ 0b21dd78-bdd1-453b-92c5-4b5918cfed85
let
	nodes = [0.0, 0.45, 0.7, 1.0]

	p = plot(; layout=(4, 1), xlabel=L"x", ylims=[-0.1, 1.1], ytick=[1.0])
	for i in 1:length(nodes)
		Hᵢ = hatfun(nodes, i)
		plot!(p, Hᵢ, 0, 1; subplot=i, label=LaTeXString("\$H_$i\$"), lw=2, c=i)
		for (ii, node) in enumerate(nodes)
			scatter!(p, [node], [Hᵢ(node)]; subplot=i, c=ii, label="")
		end
	end
	p
end

# ╔═╡ 5664b4ab-e9b3-4b31-9ef4-ae57c95691f6
md"""
Since each function $H_i$ is globally continuous
and moreover linear inside every interval $[x_i, x_{i+1}]$,
any linear combination $\sum_{i=1}^{n+1} c_i H_i(x)$
of hat functions will have the same property.
Conversely, since the piecewise linear functions form a vector space,
every piecewise linear functions is expressible
in the hat function basis, i.e. in particular
for our piecewise linear interpolant we have
```math
p_{1,h}(x) = \sum_{i=1}^{n+1} c_i H_i(x)
```
for some coefficients $c_i$, $i = 1, \ldots n+1$.

Similar to the Lagrange polynomials the hat functions
satisfy a **cardinality condition**
```math
\tag{11}
H_i(x_j) = δ_{ij} \qquad \text{for $i,j = 1, \ldots n+1$}.
```
With this in mind interpolating a piecewise linear polynomial
becomes analogous to (7) the simple expression
```math
\tag{12}
p_{1,h}(x) = \sum_{i=1}^{n+1} y_i H_i(x)
```
(i.e. the coefficients in above expressions $c_i = y_i$).

As a result piecewise linear interpolation becomes easy to implement:
"""

# ╔═╡ 863e625f-4c9e-4108-af4b-009dadad7809
function pwlinear(x, y)
	# Construct a piecewise linear interpolation for the data values (xᵢ, yᵢ)
	# given in the vectors x = [x₁, …, xₙ₊₁] and y = [y₁, …, yₙ₊₁]
	H = [hatfun(x, i) for i in 1:length(x)]
	function interpolant(xx)
		sum(y[i] * H[i](xx) for i in 1:length(x))  # (12)
	end
end

# ╔═╡ 8009f766-57d5-4bfe-bc90-759a50e863bd
md"""
We try this interpolation algorithm on our previous challenging example
```math
f_\text{ratio}(x) = \frac{1}{1 + 20x^2}
```
using `n_pwlinear` *equispaced* nodes in $[-1, 1]$:
- `n_pwlinear = ` $(@bind n_pwlinear Slider(2:30; default=10, show_value=true))
"""

# ╔═╡ f26fb942-a891-4555-8c1a-9ecc1ffeca85
let
	nodes  = range(-1, 1; length=n_pwlinear)
	p₁ₕ = pwlinear(nodes, fratio.(nodes))
	
	p = plot(fratio; xlims=(-1, 1), ylims=(0.0, 1.1), lw=2, label=L"f_\textrm{ratio}", ls=:dash,
	         title="Piecewise linear interpolation")
	scatter!(p, nodes, fratio.(nodes); label="nodes")
	plot!(p, p₁ₕ, lw=2; label=L"p_{1,h}")
end

# ╔═╡ f6055e59-1ddd-4148-8a85-95f4501e3f9f
let
	p = plot(; xlims=(-1, 1), title="Error for piecewise interpolation")
	for n_nodes in (5, 10, 30, 40)
		nodes = range(-1, 1; length=n_nodes)
		p₁ₕ = pwlinear(nodes, fratio.(nodes))
		
		error(x) = p₁ₕ(x) - fratio(x)
		plot!(p, error, lw=2; label="$(n_nodes) nodes")
	end
	p
end

# ╔═╡ 193131ad-e016-4f2a-b1cb-a544bc497c95
md"""
Unlike the case of polynomial interpolation, we observe a nice convergence
as we increase `n_pwlinear`. Let's analyse this in more detail in the next section.
"""

# ╔═╡ 7c0d831d-3e30-420b-b2c3-7c2658f7c181
md"""### Error analysis

For the error analysis of the piecewise linear interpolation
we restrict ourself to the case of *equidistant* nodes,
which we considered in the most recent exercise.
That is we assume a partitioning of the interval $[a, b]$
with $a = x_1 < x_2 < \cdots < x_{n+x} = b$
and equal nodal distance $h = x_{i+1} - x_i$.

In this setting the error analysis is a
consequence of Theorem 2, equation (8).
Indeed, for every interval $[x_i, x_{i+1}]$
we are constructing a linear interpolation
between the points $(x_i, y_i)$ and $(x_{i+1}, y_{i+1})$.
In Theorem 2 we can thus set $n=1$ and $b-a = x_{i+1} - x_i = h$
and obtain
```math
\max_{x \in [x_{i}, x_{i+1}]}
\left| f(x) - p_{1,h}(x) \right| \leq \frac{h^2}{8} \max_{x \in [x_{i}, x_{i+1}} |f''(x)|
```
therefore
```math
\begin{aligned}
\|f - p_{1,h}\|_\infty
&= \max_{x\in [a,b]} |f(x) - p_{1,h}(x)| \\
&= \max_{i=1,\ldots,n} \max_{x\in[x_{i},x_{i+1}]} |f(x) - p_{1,h}(x)| \\
&\leq \max_{i=1,\ldots,n} \frac{h^2}{8} \max_{x \in [x_{i}, x_{i+1}} |f''(x)| \\
&= \frac{h^2}{8} \| f'' \|_\infty
\end{aligned}.
```
We summarise in a Theorem:

!!! info "Theorem  4"
    Let $f : [a, b] \to \mathbb{R}$ be a $C^2$
    function and $a = x_1 < x_2 < \cdots < x_{n+x} = b$
    with equal nodal distance $h = (b - a) / n$.
    The piecewise linear polynomial interpolating the data
    $(x_i, f(x_i))$ satisfies the error estimate
    ```math
    \tag{13}
    \|f - p_{1,h}\|_\infty \leq C h^2 \| f'' \|_\infty
    ```
    with $C = 1/8$.

Note, that this theorem is only true if the second derivative of $f$
is continuous. Usually the second derivative $f''$
and thus its maximum norm $\|f\|_\infty$ are not known.
However, we obtain that the interpolation error
goes as $O(h^2)$ as $n\to \infty$.
This is an example of **quadratic convergence**.
More generally we define

!!! info "Definition: Algebraic convergence"
    If an approximation has an error with asymptotic
    behaviour $O(h^m)$ as $h\to0$ with $m$ integer
    and $h$ being a discretisation parameter
    (e.g. grid spacing, nodal distance, etc.),
	we say the approximation has **algebraic convergence**.
    If $m$ is the largest such integer
    (i.e. the error is *not* $O(h^{m+1})$)
    then $m$ is the order of accuracy.

    Moreover we often refer to the case $m=1$ as **linear convergence**,
    $m=2$ as **quadratic convergence** and so on.

!!! danger "Potential confusion: Convergence terminology"
    Note again the difference in convergence terminology between **iterative methods** and **approximation schemes**. What is called algebraic convergence for approximation schemes would be sub-linear convergence for iterative schemes.
"""

# ╔═╡ 7e317807-d3ee-4197-91fb-9fc03f7297e8
md"""
Let us finally illustrate the quadratic convergence of piecewise polynomial
interpolation graphically. Since the error $\|f - p_{1,h}\|_\infty \sim C h^2$ as $h\to0$, taking logarithms on both sides we obtain
```math
\log \left(\|f - p_{1,h}\|_\infty\right) \sim \log(C) + 2 \log(h).
```
Therefore the logarithm of the error is a *linear function*
in the logarithm of the discretisation parameter.
Moreover the **slope** gives us the convergence order,
here $m = 2$.
Note, that to obtain this behaviour any logarithm would suffice.
We choose a $log_{10}$-$\log_{10}$ plot as this is most easily realised in `Plots.jl`
and indeed exhibits a slope of $2$, meaning quadratic convergence (note that the $x$-axis is reversed in the plot.)
"""

# ╔═╡ eb96f944-74ac-4cc0-9322-825df83f34f2
let
	fine = range(-1, 1, length=1000)  # Very fine grid
	
	n = 5:5:50   # Number of nodes
	maxerror = Float64[]  # Empty array for Float64 numbers
	for (i, n_nodes) in enumerate(n)
		nodes = range(-1, 1; length=n_nodes)
		p₁ₕ = pwlinear(nodes, fratio.(nodes))

		push!(maxerror, maximum(p₁ₕ.(fine) - fratio.(fine)))
	end

	h = 2 ./ n   # Discretisation parameter  == [2/element for element in n]
	p = plot(h, maxerror; label="error", title="Convergence piecewise linear",
	         xlabel=L"h", xflip=true, xscale=:log10, ylabel=L"|| f-p_{1,h} ||_\infty",
	         yscale=:log10, mark=:o, legend=:topright)

	# Generate guiding slope
	order2 = (h ./ h[1]) .^ 2
	plot!(p, h, order2, label=L"O(h^2)", ls=:dash)

	p
end

# ╔═╡ 9d175a21-71e0-4df8-bbd5-f4eedbeb1384
md"""
To complete the definition of convergence classes,
let us return to the case of Chebyshev polynomial interpolation.
We already mentioned below Theorem 3,
that this method exhibits exponential convergence.
While for this setting the nodal points are not equally
spread, the definition of the nodal points
($x_k = - \cos\left(\frac{k \pi}{n}\right)$
for $k = 0, 1, \ldots n$)
allows to identify a discretisation parameter $h \sim 1/n$,
which scales the distance between the nodes.
We obtain the definition:

!!! info "Definition: Exponential convergence"
	If an approximation has an error scaling as $O(c^{-α h})$
    with $c>0$ and $α > 0$ two constants and $h$ a discretisation parameter,
    we say the approximation has **exponential convergence**.

To **determine the convergence behaviour graphically** we thus need to look at
- a log-log plot ($\log(y)$ versus $\log(x)$) if we suspect *algebraic convergence*.
  In this case the convergence will be a straight line with slope $m$.
- a log-linear plot ($\log(y)$ versus $x$) if we suspect *exponential convergence*.
  Again a straight line with slope $-α$ will result.
"""

# ╔═╡ 3eaddcc1-bddb-45a5-9e0d-0961dba5583c
md"""### Stability analysis

We are again interested in the effect of measurement noise
on the quality of the polynomial interpolation.
The goal is to compare the interpolation of a piecewise linear polynomial $p_{1,h}$
employing the noise-free data $(x_i, y_i)$ with $y_i = f(x_i)$
with the polynomial $\tilde{p}_{1,h}$ obtained from the
noisy data $(x_i, \tilde{y}_i)$ with $\tilde{y}_i = f(x_i) + \varepsilon_i$
and $|\varepsilon_i| \leq ε\ \forall i = 1, \ldots, n+1$.
Due to (12) we can directly write
```math
p_{1,h}(x) = \sum_{i=1}^{n+1} y_i H_i(x)
```
and
```math
\tilde{p}_{1,h}(x) = \sum_{i=1}^{n+1} \tilde{y}_i H_i(x) = \sum_{i=1}^{n+1} (y_i + ε_i) H_i(x).
```
By linearity we note
```math
\tag{14}
|p_{1,h}(x) - \tilde{p}_{1,h}(x) |
= \left| \sum_{i=1}^{n+1} ε_i H_i(x) \right| 
\leq ε \sum_{i=1}^{n+1} |H_i(x)|
= ε \sum_{i=1}^{n+1} H_i(x)
\stackrel{(\ast)}{=} ε 
```
where in the last step $(\ast)$ we used the partitition of unity
property, which you are asked to prove as an exercise.

!!! exercise
    Use (12) to prove the **partition of unity**,
    namely $\sum_{i=1}^{n+1} H_i(x) = 1$.

Based on (14) we can conclude:
- The condition number of piecewise linear interpolation is $1$ independent
  of the data.
- As a result small errors in the input values $x_i$ (e.g. from measurements)
  only introduce small perturbations in the polynomial:
  Piecewise linear interpolation is numerically stable.

We conclude by illustrating this result graphically.
Again we consider the function
$f_\text{sin}(x) = \sin(5x)$ on the interval $[-1, 1]$, 
once evaluated without noise and once taking $\tilde{y}_i = f(x_i) + η_i$
where $|\eta_i| ≤ 10^\texttt{η\_log\_pl}$.
For construct the piecewise linear interpolation
we take $\texttt{n\_nodes\_pl}$ equally spaced nodes.
"""

# ╔═╡ 2096258a-efd9-4de1-8f56-c02295407c0b
md"""
- Number of nodes `n_nodes_pl = ` $(@bind n_nodes_pl Slider(2:2:40; show_value=true, default=20))
- Logarithm of noise amplitude `η_log_pl = ` $(@bind η_log_pl Slider(-3:0.25:0; show_value=true, default=-2))
"""

# ╔═╡ e2b32c49-9c8e-4d32-8bf6-a0517b7caa8a
let
	fine  = range(-1.0, 1.0; length=500)

	nodes = range(-1, 1; length=n_nodes_pl)
	noise = [10^η_log_pl * (-1 + 2rand()) for _ in nodes]

	p₁ₕ_accurate = pwlinear(nodes, fsin.(nodes))
	p₁ₕ_noisy    = pwlinear(nodes, fsin.(nodes) + noise)

	p = plot()
	plot!(p, fine, fsin.(fine); title="Piecewise linear",
	      label=L"function $f$", lw=2, ls=:dash, ylim=(-2, 2), legend=:topright)
	plot!(p, fine, p₁ₕ_accurate.(fine); label=L"accurate $p_{1,h}$", lw=2)
	plot!(p, fine, p₁ₕ_noisy.(fine); label=L"noisy $\tilde{p}_{1,h}$", lw=2)
	scatter!(p, nodes, fsin.(nodes) + noise; label="")

	q = plot(fine, abs.(fsin.(fine) .- p₁ₕ_noisy.(fine)), label=L"noisy fit $\tilde{p}_{1,h}$", yaxis=:log, xlim=(-1, 1), ylims=(10e-5, 1), legend=:topright, c=3, lw=1.5, ylabel="error")
	hline!(q, [10^η_log_pl], label=L"noise $\eta$", c=3, lw=1, ls=:dash)
	plot!(q, fine, abs.(fsin.(fine) .- p₁ₕ_accurate.(fine)), label=L"fit $p_{1,h}$", c=2, lw=1.5)

	plot(p, q; layout=grid(2, 1; heights=[0.75, 0.25]))
end

# ╔═╡ cf296fbf-c900-485c-8cfe-672f02ad3b69
md"""
## Spline interpolation

Over polynomial interpolation
our previously discussed piecewise linear interpolation has the advantage
that it is *always* convergent as the number of nodal points increases,
even for an equidistant scheme. Moreover it enjoys a condition number of $1$,
thus is numerically very stable. However, it suffers from a few flaws,
which make it unsuitable for many applications:
- The interpolant $p_{1,h}$ is globally only continuous since at all nodal points it is not neccessarily differentiable.
- Piecewise linear interpolation only converges quadratically.
  To understand why this is a problem, consider the following convergence curves
  on our example $f_\text{sin} = \sin(5x)$.
  We simply need a lot more nodal points (and thus data / measurements)
  to get the same accuracy !

Show splines: $(@bind show_splines CheckBox())
"""

# ╔═╡ a0a8a8ba-3f25-4245-b25b-84333046af0d
md"""
An alternative approach is *cubic spline interpolation*:

!!! info "Definition: Spline interpolation"
    Let $(x_i, y_i), i = 1, \ldots, n + 1$ denote the available data
    with ordering $a = x_1 < x_2 < \cdots < x_{n+1} = b$ of the nodal points.
    An **interpolating cubic spline** is a function $s_{3,h}$ such that

    1. Between the nodal points, i.e. within the intervals $[x_i, x_{i+1}]$
       the spline $s_{3,h}$ is a polynomial of degree $3$.
    2.   $s_{3,h} \in C^2([a,b])$, that is $s_{3,h}$ is globally twice differentiable.
    3.    $s_{3,h}(x_i) = y_i$, that is $s_{3,h}$ interpolates the data.

Due to condition 1. the restriction $s_{3,h}|_{[x_i, x_{i+1}]}$ to an interval
can be written as
```math
s_{3,h}|_{[x_i, x_{i+1}]}(x) = a_i + b_i x + c_i x^2 + d_i x^3
\qquad \text{for $i = 1, \ldots n$},
```
which makes in total $4n$ unknowns to be determined.

In order to satisfy the continuity condition 2. we need $s_{3,h}$ and its first and second derivative to be continuous at the nodal points, that is
```math
\begin{aligned}
s_{3,h}(x_i^-)   &= s_{3,h}(x_i^+), \qquad i = 2, \ldots, n \\
s'_{3,h}(x_i^-)  &= s'_{3,h}(x_i^+), \qquad i = 2, \ldots, n \\ 
s''_{3,h}(x_i^-) &= s''_{3,h}(x_i^+), \qquad i = 2, \ldots, n \\ 
\end{aligned}
```
where we used the notation $f(x^+) = \lim_{\xi \searrow x} f(\xi)$
and $f(x^-) = \lim_{\xi \nearrow x} f(\xi)$.
This gives $3 (n-1)$ equations.

Additionally the interpolating condition 3. gives us $n+1$ equations,
one for each of the nodes.

Summarising our findings we are thus left with
```math
4n - 3 (n-1) - (n+1) = 2
```
degrees of freedom, which we need to set to determine a unique spline.
For this there are two mjor alternatives:

- **Natural spline:**    $s_1'''(x_1) = s_1'''(x_{n+1}) = 0$
- **Not-a-knot spline:** $s_1'''(x_1) = s_2'''(x_2)$ and $s_1'''(x_n) = s_2'''(x_{n+1})$

Natural splines have properties that make them theoretically more interesting,
but not-a-knot splines give better pointwise accuracy.

Let's check how splines behave for our function
```math
f_\text{ratio}(x) = \frac{1}{1 + 20x^2}
```
using `n_spline` equispaced nodes in $[-1, 1]$.
- `n_spline = ` $(@bind n_spline Slider(2:30; default=10, show_value=true))
"""

# ╔═╡ fd469b84-9093-4b92-ac3e-c1847db3627f
function spinter(t, y)
	# Copied from https://github.com/fncbook/FundamentalsNumericalComputation.jl
    # Available under the MIT licence
	
    n = length(t)-1
    h = [ t[k+1]-t[k] for k in 1:n ]

    # Preliminary definitions.
    Z = zeros(n,n);
    In = I(n);  E = In[1:n-1,:];
    J = diagm(0=>ones(n),1=>-ones(n-1))
    H = diagm(0=>h)

    # Left endpoint interpolation:
    AL = [ In Z Z Z ]
    vL = y[1:n]

    # Right endpoint interpolation:
    AR = [ In H H^2 H^3 ];
    vR = y[2:n+1]

    # Continuity of first derivative:
    A1 = E*[ Z J 2*H 3*H^2 ]
    v1 = zeros(n-1)

    # Continuity of second derivative:
    A2 = E*[ Z Z J 3*H ]
    v2 = zeros(n-1)

    # Not-a-knot conditions:
    nakL = [ zeros(1,3*n) [1 -1 zeros(1,n-2)] ]
    nakR = [ zeros(1,3*n) [zeros(1,n-2) 1 -1] ]

    # Assemble and solve the full system.
    A = [ AL; AR; A1; A2; nakL; nakR ]
    v = [ vL; vR; v1; v2; 0; 0 ]
    z = A\v

    # Break the coefficients into separate vectors.
    rows = 1:n
    a = z[rows]
    b = z[n.+rows];  c = z[2*n.+rows];  d = z[3*n.+rows]
    S = [ Polynomial([a[k],b[k],c[k],d[k]]) for k in 1:n ]

    # This function evaluates the spline when called with a value
    # for x.
    return function (x)
        if x < t[1] || x > t[n+1]    # outside the interval
            return NaN
        elseif x==t[1]
            return y[1]
        else
            k = findlast(x .> t)    # last node to the left of x
            return S[k](x-t[k])
        end
    end
end

# ╔═╡ 22c33635-dd92-4a8d-9f99-74b222a44f2e
let
	fine = range(-1, 1, length=1000)  # Very fine grid
	p = plot(title="Approximation error", xlabel=L"n", ylabel=L"||p - f||_\infty",
		     yscale=:log10, legend=:topright)
	
	n = 8:4:80   # Number of nodes
	maxerror_linear = Float64[]
	for (i, n_nodes) in enumerate(n)
		nodes = range(-1, 1; length=n_nodes)
		p₁ₕ = pwlinear(nodes, fratio.(nodes))
		push!(maxerror_linear, maximum(p₁ₕ.(fine) - fratio.(fine)))
	end
	p = plot!(p, n, maxerror_linear; mark=:o, lw=2, label="Piecewise linear")

	maxerror_cheb = Float64[]
	for (i, n_nodes) in enumerate(n)
		n_degree = n_nodes - 1
		t = BigFloat[-cos(π*k/n_degree) for k = 0:n_degree]
		poly = Polynomials.fit(t, fratio.(t))
		push!(maxerror_cheb, maximum(poly.(fine) - fratio.(fine)))
	end
	plot!(p, n, maxerror_cheb; label="Chebyshev polynomial", mark=:o, lw=2)
	
	maxerror_spline = Float64[]
	for (i, n_nodes) in enumerate(n)
		nodes = range(-1, 1; length=n_nodes)
		spline = spinter(nodes, fratio.(nodes))
		push!(maxerror_spline, maximum(spline.(fine) - fratio.(fine)))
	end

	if show_splines
		plot!(p, n, maxerror_spline; label="Cubic splines", mark=:o, lw=2)
	end
		
	p
end

# ╔═╡ 5987bb0f-c784-42ab-9c34-d6ff76cef2ea
let
	nodes = range(-1, 1; length=n_spline)
	s₃ₕ   = spinter(nodes, fratio.(nodes))
	
	p = plot(fratio; xlims=(-1, 1), ylims=(0.0, 1.1), lw=2,
	         label=L"f_\textrm{ratio}", ls=:dash,
	         title="Cubic spline interpolation")
	scatter!(p, nodes, fratio.(nodes); label="nodes")
	plot!(p, s₃ₕ, lw=2; label=L"s_{3,h}")
end

# ╔═╡ 33dd12c8-f3d5-4238-986f-df58639a6f3a
md"""
Already around 18 nodal points there is hardly any difference between the function and the interpolant visible.

For and example as well as more details on the computation of splines see
[chapter 5.3](https://tobydriscoll.net/fnc-julia/localapprox/splines.html)
of Driscoll, Brown: *Fundamentals of Numerical Computation*.
"""

# ╔═╡ d4226809-9f2a-4651-9778-8f8f66d03720
md"""
### Error analysis

For cubic splines the following convergence result is known:

!!! info "Theorem 5"
    Let $f \in C^4([a, b])$ be $4$ times differentiable and let
    $a < x_1 < x_2 < \cdots < x_{n+1} = b$ be equispaced nodal points
    in $[a,b]$. The cubic spline $s_{3,h} interpolating the data $(x_i, f(x_i))$
    satisfies the error bounds
    ```math
    \begin{aligned}
    \|f - s_{3,h}\|_\infty &≤ C_0 \, h^4 \, \|f^{(4)}\|_\infty\\
    \|f' - s'_{3,h}\|_\infty &≤ C_1 \, h^3 \, \|f^{(4)}\|_\infty\\
    \|f'' - s''_{3,h}\|_\infty &≤ C_2 \, h^2 \, \|f^{(4)}\|_\infty
    \end{aligned}
    ```
    where $h = \frac{b-a}{n}$ is the length of each interval
    and $C_0$, $C_1$, $C_2$ are constants not depending on $h$.


We would thus expect a 4th order convergence.
Let us investigate this visually using the function
$f_\text{ratio}(x) = \frac{1}{1 + 20x^2}$.
Since we expect algebraic convergence, we employ a log-log scale:
"""

# ╔═╡ a806cafa-8ee9-4dd6-9b90-100a831c009d
let
	fine = range(-1, 1, length=1000)  # Very fine grid
	
	n = 2 .^ (4:9)
	maxerror = Float64[]  # Empty array for Float64 numbers
	for (i, n_nodes) in enumerate(n)
		nodes = range(-1, 1; length=n_nodes)
		s₃ₕ = spinter(nodes, fratio.(nodes))
		push!(maxerror, maximum(s₃ₕ.(fine) - fratio.(fine)))
	end

	h = 2 ./ n
	p = plot(h, maxerror; label="error", title="Convergence piecewise linear",
	         xlabel=L"h", xflip=true, xscale=:log10, ylabel=L"|| f-s_{3,h} ||_\infty",
	         yscale=:log10, mark=:o, legend=:topright, lw=2)

	# Generate guiding slope
	order4 = (h ./ h[1]) .^ 4
	plot!(p, h, order4, label=L"O(h^4)", ls=:dash, lw=2)

	p
end

# ╔═╡ ef07aa33-4127-4c30-891a-69e5154bc1a8
md"""
## Regression and curve fitting

Consider again the case where we have acquired $n$ data points
$(x_i, \tilde{y}_i)$, $i=1, 2, \ldots, n$ where 
```math
\tilde{y}_i = f(x_i) + ε_i
```
and the measurement noice $ε_i$ is substantial.
If we were to use an *interpolation technique*
--- as discussed in the previous sections ---
our goal would be to obtain a (potentially piecewise)
polynomial $p$ such that $p(x_i) = \tilde{y}_i$.
However, since $\tilde{y}_i$ is only a *noisy* observation
it makes very little sense to force the fitted polynomial
to go through the data *exactly*.

In this section we will now consider a generalisation to interpolation,
where we give up on the interpolation condition $q(x_i) = \tilde{y}_i$.
Instead we will now seek the representative $q \in \mathbb{P}$
taken from a model space $\mathbb{P}$,
which best represents the data.

For example, in **least-squares** **linear regression** one seeks
the straight line $q(x) = ax + b$, which yields the lowest squared error
```math
\tag{15}
e_\text{LS} = \sum_{i=1}^{n} |p_1^\text{LS}(x_i) - y_i|^2.
```
If we denote by $\mathbb{P}_m$ the space of all polynomials of degree
at most $m$, then linear regression seeks the $q \in \mathbb{P}_1$,
which gives smallest $e_\text{LS}$, or more mathematically
```math
p_1^\text{LS} = \textrm{argmin}_{q\in \mathbb{P}_1} |q(x_i) - y_i|^2.
```

A generalisation to higher-order polynomials can be formally defined as

!!! info "Definition: Least-squares m-th degree polynomial regression"
    Given $n$ data points $(x_i, \tilde{y}_i)$, $i = 1, 2, \ldots, n$
    the $m$-th degree polynomial $p_m^\text{LS} \in \mathbb{P}_m$
    satisfying
    ```math
    \tag{16}
    p_m^\text{LS} = \textrm{argmin}_{q\in \mathbb{P}_m} |q(x_i) - y_i|^2.
    ```
    is called $m$-th degree polynomial least squares approximation
    of the data. In other words $p_m^\text{LS}$ is such that
    ```math
    \sum_{i=1}^{n} \left( y_i - p_m^\text{LS}(x_i) \right)^2
    \leq \sum_{i=1}^{n}\left( y_i - q(x_i) \right)^2
    \qquad \forall q \in \mathbb{P}_m.
    ```

Least-squares regression is usually much better suited to capture trends
in data. To see this consider the following setting, where we fit
polynomials of various degree.
Our example data is the $5$-year averages of the worldwide temperature anomaly as compared to the 1951--1980 average:
"""

# ╔═╡ 5c5b212d-0451-4334-b5f7-273640414e86
begin
	year = 1955:5:2000
	temp = [ -0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
		      0.1180, 0.2100, 0.3320, 0.3340, 0.4560 ]
end;

# ╔═╡ 7bd4720f-3ec8-459a-91be-17d5f2acaa5c
md"In the plot we show interpolation (i.e. $m = n-1 = 9$ in this case)
as well as polynomial regression with $m = 1$, $m=3$ and $m = 5$."

# ╔═╡ 66cfe319-c19e-478f-8f36-645fa1ff8b3b
let
	fine = range(1950, 2005, length=500)
	
	p = plot()
	for m in (1, 3)
		poly = Polynomials.fit(BigFloat.(year), BigFloat.(temp), m)
		plot!(p, fine, poly.(fine); label="Regression m = $m", lw=2, c=m)
	end
	poly = Polynomials.fit(BigFloat.(year), BigFloat.(temp))
	plot!(p, fine, poly.(fine); label="Interpolation", lw=2, ls=:dash, c=2)
	scatter!(p, year, temp; label="data", c=1)

	xlims!(p, 1950, 2005)
	ylims!(p, -0.2, 0.5)
	p
end

# ╔═╡ f6916b49-4fce-429c-9f55-d9a51dc46937
md"""
As expected for such noisy data: The lower-degree polynomials seem to do a much better job.

Note that in polynomial regression
the choice of the error metric (15) essentially determines which of the members of the model space $\mathbb{P}_m$ is considered to be the "best fit"
--- in the sense of being the minimiser of the minimisation problem (16).
In this lecture we will only consider **least squares regression**,
for which the error metric is the **sum of squares error**
```math
e_\text{LS} = \sum_{i=1}^{n} |q(x_i) - y_i|^2.
```
Other choices are for example to employ
-   $\max_i |q(x_i) - y_i|$, the maximal elementwise deviation
-    $\sum_{i=1}^{n+1} |q(x_i) - y_i|$, absolute deviations
However, we will not consider these further.
"""

# ╔═╡ 2122f3c9-a901-4488-bd1b-4bd7e49b7c8a
md"""
### Least squares problems
We will now consider how to solve polynomial regression problems (16).
For this we introduce some linear algebra.

Recall from equation (1) at the start of discussion of polynomial interpolation,
that each polynomial $q \in \mathbb{P}_m$ can be written as a linear combination
```math
q(x) = \sum_{j=0}^m c_j φ_j(x)
```
where $φ_j$ are suitable basis functions, such as the monomials $φ_j(x) = x^j$.
Inserting this, the least-squares error expression can be rewritten as
```math
e_\text{LS} = \sum_{i=1}^{n} \left| y_i - \sum_{j=0}^m c_j φ_j(x_i) \right|^2
```
and further by introducing the matrix / vector notation
```math
\textbf{V} = \left(\begin{array}{cccc}
φ_0(x_1) & φ_1(x_1) & \cdots & φ_m(x_1) \\
φ_0(x_2) & φ_1(x_2) & \cdots & φ_m(x_2) \\
& \vdots\\
φ_0(x_{n}) & φ_1(x_{n}) & \cdots & φ_m(x_{n})
\end{array}\right) \qquad
\textbf{y} = \left(\begin{array}{c} y_1 \\ y_2 \\ \vdots \\ y_n \end{array}\right)
\qquad
\textbf{c} = \left(\begin{array}{c} c_1 \\ c_2 \\ \vdots \\ c_n \end{array}\right)
```
as
```math
\tag{16}
e_\text{LS} = \| \mathbf{y} - \mathbf{V} \mathbf{c} \|_2^2
```
where $\| v \|_2 = \sqrt{ \sum_{i=1}^n v_i^2 }$ is the Euklidean norm.
If we are using the monomial basis
the matrix $\mathbf{V}$ is again equal to a Vandermode matrix,
same as in polynomial interpolation.
However, since for regression problems we usually have that $n > m + 1$
it is rectangular in this case:
```math
\mathbf{V} = \left(\begin{array}{ccc}
1 & x_1 & \ldots & x_1^m \\
1 & x_2 & \ldots & x_2^m \\
\vdots \\
1 &x_n & \ldots & x_n^m \\
\end{array}\right) \in \mathbb{R}^{n\times (m+1)}
```

In polynomial regression our job is now to minimise expression (16),
which means that we want to find the coefficient vector $\mathbf{c}$,
which minimises $\|\mathbf{y} - \mathbf{V} \mathbf{c} \|_2^2$.
A procedure to obtain this coefficient vector from
of $\mathbf{y}$ and $\mathbf{V}$ is obtained by noting
that polynomial regression is just a special case of a
general least-squares problem, defined as:

!!! note "Definition: Least-squares problem"
    Given $\mathbf{A} \in \mathbb{R}^{k\times l}$ and $\mathbf{b} \in \mathbb{R}^k$
    with $k > l$ find
    ```math
    \tag{17}
    \text{argmin}_{\mathbf{x} \in \mathbb{R}^l} \| \mathbf{b} - \mathbf{A} \mathbf{x} \|_2^2.
    ```
    The vector $\textbf{r} = \mathbf{b} - \mathbf{A} \mathbf{x} \in \mathbb{R}^k$
    is referred to as the **residual** of the least-squares problem.
"""

# ╔═╡ 30899dd5-b307-415d-bb94-efd09e7d0864
md"""
Least squares problems arise in many applications (statistics, machine learning, ...)
as well as polynomial regression. For the latter case the Vandermonde matrix $\mathbf{V}$ plays the role of $\mathbf{A}$,
the vector of $y_i$-values $\mathbf{y}$ the role of $\textbf{b}$
and the vector of coefficients $\textbf{c}$ the role of the unknown $\textbf{x}$.

Solving (17) can in fact be achieved by a concise explicit expression:

!!! note "Theorem 6"
    Let $\mathbf{A} \in \mathbb{R}^{k\times l}$ and $\mathbf{b} \in \mathbb{R}^k$
    with $k > l$.
    If $\mathbf{x} \in \mathbb{R}^l$ satisfies $\textbf{A}^T (\textbf{A}\textbf{x} - \textbf{b}) = \textbf{0}$ then $\textbf{x}$ solves the least-squares problem, i.e.
    $\textbf{x}$ is the minimiser of $\|\textbf{A} \textbf{x} - \textbf{b}\|$.

> **Proof:** We first consider the elementary identity
> for the sum of two vectors $\textbf{u}$ and $\textbf{v}$:
> ```math
> \tag{18}
> \begin{aligned}
> \|\textbf{u} + \textbf{v}\|^2 &=
> (\textbf{u} + \textbf{v})^T(\textbf{u} + \textbf{v})
> \\&= \textbf{u}^T\textbf{u} + \textbf{u}^T\textbf{v} + \textbf{v}^T\textbf{u} +  \textbf{v}^T\textbf{v} \\
> &= \|\textbf{u}\|^2 + 2 \textbf{u}^T\textbf{v} + \|\textbf{v}\|^2
> \end{aligned}
> ```
>
> Now let $\textbf{y} \in \mathbb{R}^l$ be an arbitrary vector
> and set $\textbf{u} = \textbf{A}\textbf{x} - \textbf{b}$
> and $\textbf{v} = \textbf{A} \textbf{y}$ in the above development.
> This results in
> ```math
> \begin{aligned}
> \|\textbf{A}(\textbf{x} + \textbf{y}) - \textbf{b}\|^2
> &= \|(\textbf{A}\textbf{x} - \textbf{b}) + \textbf{A}\textbf{y}\|^2 \\
> &\stackrel{(18)}{=} \|\textbf{A}\textbf{x} - \textbf{b}\|^2 + 2 (\textbf{A}\textbf{y})^T (\textbf{A}\textbf{x} - \textbf{b})
> + \|\textbf{A}\textbf{y}\|^2\\
> &\geq \|\textbf{A}\textbf{x} - \textbf{b}\|^2 + 2 \textbf{y}^T \textbf{A}^T (\textbf{A}\textbf{x} - \textbf{b})
> \end{aligned}
> ```
> Therefore, if $\textbf{A}^T (\textbf{A}\textbf{x} - \textbf{b}) = 0$, then
> ```math
> \|\textbf{A}(\textbf{x} + \textbf{y}) - \textbf{b}\|^2 \geq \|\textbf{A}\textbf{x} - \textbf{b}\|^2 \qquad \text{for all $\textbf{y} \in \mathbb{R}^l$},
> ```
> which implies that $\textbf{x}$ is the minimiser of $\|\textbf{A}\textbf{x} - \textbf{b}\|$.
"""

# ╔═╡ d846d3f2-8cfc-4fc5-a22b-d55762f88f45
md"""
Due to the overall importance of least-squares problems
the solution equation $\textbf{A}^T (\textbf{A}\textbf{x} - \textbf{b}) = \textbf{0}$
is usually referred to as the normal equation.

!!! note "Definition: Normal equations"
    Given $\mathbf{A} \in \mathbb{R}^{k\times l}$ and $\mathbf{b} \in \mathbb{R}^k$,
    which define a least-squares problem $\text{argmin}_{\textbf{x}\in \mathbb{R}^l} \|\mathbf{b} - \mathbf{A} \mathbf{x}\|$, the solution equation
    ```math
    \textbf{A}^T (\textbf{A}\textbf{x} - \textbf{b}) = \textbf{0}
    ```
    or equivalently
    ```math
    \tag{19}
    \textbf{A}^T \textbf{A} \textbf{x} = \textbf{A}^T \textbf{b}.
    ```
    are called the **normal equations**.

For a geometric interpretation of the normal equations, see [chapter 3.2](https://tobydriscoll.net/fnc-julia/leastsq/normaleqns.html) of Driscoll, Brown: *Fundamentals of Numerical Computation*.

Based on this strategy we can now perform polynomial regression.
We again employ the word temperature data between 1950 and 2000.
You can change the polynomial degree using the slider:

- `m_poly = ` $(@bind m_poly Slider(1:9; default=2, show_value=true))
"""

# ╔═╡ 0b8f2bcc-0948-47ed-bd62-4a0ceeee9156
let
	x = year
	y = temp
	
	# Build Vandermonde matrix (using BigFloat to avoid numerical issues)
	V = ones(BigFloat, length(x), m_poly + 1)
	for k in 2:m_poly+1
		V[:, k] = V[:, k-1] .* x
	end

	# Solve normal equations
	c = (V'V) \ (V' * y)

	# Construct polynomial
	polynomial(x) = evalpoly(x, c)

	# Plot result
	p = scatter(year, temp; label="data")
	plot!(p, polynomial; label="Regression m = $m_poly", lw=2)

	p
end

# ╔═╡ 1eccef57-f63b-4fba-a8f6-458652686fb1
TODO("First show the interactive plot and use it to illustrate the key points as a motivation, then explain the theory.")

# ╔═╡ bb5988bf-21eb-4078-b8b7-5a8dedf2ff6e
md"""
### Error analysis

In the general setting the analysis of least-squares problems is rather involved. Here, we will only limit ourselves to stating some key results, which illustrate the key properties of the method.
We recall the setting: We are given $n$ data points 
$(x_i, \tilde{y}_i)$, $i=1, 2, \ldots, n$ where 
```math
\tilde{y}_i = f(x_i) + ε_i
```
are the noisy observations of a function $f$.
Further we assume that $\mathbb{E}(ε_i) = 0$, i.e. the noise follows a distribution
with mean zero and $σ = \sqrt{\mathbb{V}\text{ar}(ε_i)}$ is the square root of the variance
of the noise.

- In the **absence of measurement noise**, i.e. $σ=0$, least-squares can recover polynomials exactly. Imagine $f$ to be a polynomial of degree $m \ll n$. In this case the polynomial $f$ is the *only* $m$-th degree polynomial to give a zero least-squares error, i.e.
  ```math
  0 = \sum_{i=1}^{n} |p_m^\text{LS}(x_i) - f(x_i)|^2 \text{ with $p_m^\text{LS}\in\mathbb{P}_m$}
  \quad \Leftrightarrow \quad p_m^\text{LS} = f.
  ```
  and our result is exact.
- With **non-zero measurement errors** the approximation error $\|f - p_m^\text{LS}\|_\infty$ is proportional to $\sqrt{\frac{m+1}{n}}σ$. Thus while every measurement has been perturbed by an error $ε_i$ of order around $σ$, the approximation error is much smaller of order $\sqrt{(m+1)/n}\,σ$ and moreover becomes smaller as the number of samples increases.

  However, if the degree of the fitting polynomial $m$ becomes too large --- say comparable to $n$ ---, then the least-squares polynomial becomes similar to the
  interpolating polynomial. This we saw to become unstable as $m$ gets large.
  Therefore the least-squares approximation in general does similarly deteriorate
  when $m \to n$.
"""

# ╔═╡ 4562b2f4-3fcb-4c60-83b6-cafd6e4b3144
md"""
- Number of samples `n = ` $(@bind n Slider(10:5:50; show_value=true, default=20))
- Polynomial degree `m = ` $(@bind m Slider(1:20; show_value=true, default=2))
- Logarithm of noise amplitude `η_log = ` $(@bind η_log Slider(-3:0.25:0; show_value=true, default=-2))
"""

# ╔═╡ 7829afd0-2693-40dc-b2d5-c8da72e96454
let
	fine  = range(-1.0, 1.0; length=500)

	nodes = range(BigFloat(-1), BigFloat(1); length=n)
	noise = [10^η_log * (-1 + 2rand()) for _ in 1:n]
	data  = fsin.(nodes) + noise

	pₘᴸˢ = Polynomials.fit(nodes, data, m+1)

	p = plot()
	plot!(p, fine, fsin.(fine); title="Polynomial regression",
	      label=L"function $f$", lw=2, ls=:dash, ylim=(-2, 2), legend=:topright)
	scatter!(p, nodes, data; label="data")
	plot!(p, fine, pₘᴸˢ.(fine); label=L"regression $p_m^\textrm{LS}$", lw=2.5)

	q = plot(fine, abs.(fsin.(fine) .-  pₘᴸˢ.(fine)), label="error", yaxis=:log, xlim=(-1, 1), ylims=(10e-5, 1), legend=:bottomright, c=3, lw=1.5)
	hline!(q, [10^η_log], label=L"noise $\eta$", c=3, lw=1, ls=:dash)

	plot(p, q; layout=grid(2, 1; heights=[0.75, 0.25]))
end

# ╔═╡ d6e89a4c-e59b-462d-8db6-0bc1a98ac381
TODO("Summarise Common features and differences between polynomial interpolation and least-squares regression.")

# ╔═╡ 2240f8bc-5c0b-450a-b56f-2b53ca66bb03
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 530)
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
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Plots = "~1.40.0"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.55"
Polynomials = "~1.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "49189c2d7c19149f56d539d92c1dff2ba497fc23"

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
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

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
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "d476eaeddfcdf0de15a67a948331c69a585495fa"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.47.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "8cfa272e8bdedfa88b6aefbbca7c19f1befac519"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.0"

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
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

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
git-tree-sha1 = "8e2d86e06ceb4580110d9e716be26658effc5bfd"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "da121cbdc95b065da07fbb93638367737969693f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.8+0"

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
git-tree-sha1 = "abbbb9ec3afd783a7cbd82ef01dcd088ea051398"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.1"

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

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "ac0aaa807ed5eaf13f67afe188ebc07e828ff640"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.10.0"

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
git-tree-sha1 = "04663b9e1eb0d0eabf76a6d0752e0dac83d53b36"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.28"

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
git-tree-sha1 = "8c57307b5d9bb3be1ff2da469063628631d4d51e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.21"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    DiffEqBiologicalExt = "DiffEqBiological"
    ParameterizedFunctionsExt = "DiffEqBase"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    DiffEqBiological = "eb300fae-53e8-50a0-950c-e21f52c2b7e0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

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
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "4cc0c5a83933648b615c36c2b956d94fda70641e"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.7"

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

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

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
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

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
git-tree-sha1 = "38a748946dca52a622e79eea6ed35c6737499109"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.0"

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
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[deps.Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

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

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

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

[[deps.SpecialFunctions]]
deps = ["OpenSpecFun_jll"]
git-tree-sha1 = "d8d8b8a9f4119829410ecd706da4cc8594a1e020"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "0.10.3"

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
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "d39314cdbaf5b90a047db33858626f8d1cc973e1"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.0.0+2023c"

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

[[deps.TimeZones]]
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "89e64d61ef3cd9e80f7fc12b7d13db2d75a23c03"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.13.0"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

[[deps.TranscodingStreams]]
git-tree-sha1 = "1fbeaaca45801b4ba17c251dd8603ef24801dd84"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.2"
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

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

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

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

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
# ╠═679517fb-6090-4259-9b93-571ded17ab88
# ╟─61e5ef66-a213-4b23-9406-9cc63a58104c
# ╟─8c00d7cc-c050-11ee-23f8-e90afbb1a5e9
# ╟─eb43ef87-343e-4a9b-9a0f-3a18364d6b4b
# ╠═9e933ebe-1805-46fa-9e0a-a64c03a71c5c
# ╠═f9f4ddec-fd3d-4624-ace4-f316b9cfa1e4
# ╠═8be732e2-927b-4c56-8cab-9079c7598e86
# ╟─b73ea6d2-5751-4926-a28d-418deb60da2a
# ╟─2dfef7e4-23ca-49fa-bb95-f9c7219cdf22
# ╟─2a088353-0b78-45a8-b354-6ced35b1c7f9
# ╠═ac94e6f2-ebf0-4ee8-b21f-d9d20c656dd5
# ╟─4c1412f2-c1ac-4f95-afec-2c3ef1b09dd2
# ╠═5a1e03f2-07fc-4a31-9715-3b44d12824a2
# ╟─317058ec-7690-4eea-8fed-b9cd8190cc7b
# ╟─f7e3bfc6-1d12-4390-add7-12d0b0f8f4e1
# ╟─2580774f-b4be-4665-9123-bc7c22210010
# ╟─d6c4e000-c213-40e0-8a6a-e1dfccadb006
# ╟─f6fdfe7b-fa92-4168-bc8b-1f3a5a30f8db
# ╟─463fcdfa-89b0-4f9d-aa22-29d522585627
# ╟─8bfdb2e5-b3bf-4a17-a1cb-161eec581df8
# ╟─2b0638bd-2ab1-49ff-921b-a86c4e49c002
# ╟─5ecbe1d5-2250-4bab-8c1e-c319d23a0f4b
# ╠═f40c7c2e-339c-4b6e-91d3-ef5d031091f3
# ╟─dada1011-4a09-440b-9b92-ac167c8028fe
# ╠═b1773d1c-4d3c-4f4b-805d-f188623c9571
# ╟─e567287b-9244-4b55-9f5b-2e3bf583e46a
# ╟─31126101-685d-464d-be4a-6233d3803780
# ╟─af534a92-31d4-4125-bb7c-3eb2b3430852
# ╟─b92ee89f-9241-4b4d-a54f-dafe4900e159
# ╟─de28f64c-f243-4cb5-a1e6-ae8562cce5bc
# ╟─08d5cb8a-3bf5-4d40-a8d6-d22318e37712
# ╟─c5b33503-973b-468f-95de-acd6a7605a99
# ╟─0a74569f-fd4b-4705-8350-3b47db65f502
# ╟─c38b9e48-98bb-4b9c-acc4-7375bbd39ade
# ╟─479a234e-1ce6-456d-903a-048bbb3de65a
# ╟─d4cf71ef-576d-4900-9608-475dbd4d933a
# ╟─56685887-7866-446c-acdb-2c20bd11d4cd
# ╟─1d61d925-1ceb-439b-aff6-c5e8483fdc31
# ╟─9766ff2e-2060-4508-8018-8aa6d51219bc
# ╟─d5de3b11-7781-4100-8a63-d26426685bbc
# ╠═4a4540a3-b4cf-47fd-a615-a5280505333f
# ╟─50a3face-5831-4ef5-b6a5-225b8b7c46a0
# ╟─2cefad7d-d734-4342-8039-aefdc33c2edd
# ╠═611541de-927a-44fd-bdf6-bedec2befa01
# ╟─1af1365e-de88-4140-861d-ca527ae67db7
# ╠═0b21dd78-bdd1-453b-92c5-4b5918cfed85
# ╟─5664b4ab-e9b3-4b31-9ef4-ae57c95691f6
# ╠═863e625f-4c9e-4108-af4b-009dadad7809
# ╟─8009f766-57d5-4bfe-bc90-759a50e863bd
# ╠═f26fb942-a891-4555-8c1a-9ecc1ffeca85
# ╠═f6055e59-1ddd-4148-8a85-95f4501e3f9f
# ╟─193131ad-e016-4f2a-b1cb-a544bc497c95
# ╟─7c0d831d-3e30-420b-b2c3-7c2658f7c181
# ╟─7e317807-d3ee-4197-91fb-9fc03f7297e8
# ╠═eb96f944-74ac-4cc0-9322-825df83f34f2
# ╟─9d175a21-71e0-4df8-bbd5-f4eedbeb1384
# ╟─3eaddcc1-bddb-45a5-9e0d-0961dba5583c
# ╟─2096258a-efd9-4de1-8f56-c02295407c0b
# ╠═e2b32c49-9c8e-4d32-8bf6-a0517b7caa8a
# ╟─cf296fbf-c900-485c-8cfe-672f02ad3b69
# ╟─22c33635-dd92-4a8d-9f99-74b222a44f2e
# ╟─a0a8a8ba-3f25-4245-b25b-84333046af0d
# ╠═5987bb0f-c784-42ab-9c34-d6ff76cef2ea
# ╟─fd469b84-9093-4b92-ac3e-c1847db3627f
# ╟─33dd12c8-f3d5-4238-986f-df58639a6f3a
# ╟─d4226809-9f2a-4651-9778-8f8f66d03720
# ╠═a806cafa-8ee9-4dd6-9b90-100a831c009d
# ╟─ef07aa33-4127-4c30-891a-69e5154bc1a8
# ╠═5c5b212d-0451-4334-b5f7-273640414e86
# ╟─7bd4720f-3ec8-459a-91be-17d5f2acaa5c
# ╟─66cfe319-c19e-478f-8f36-645fa1ff8b3b
# ╟─f6916b49-4fce-429c-9f55-d9a51dc46937
# ╟─2122f3c9-a901-4488-bd1b-4bd7e49b7c8a
# ╟─30899dd5-b307-415d-bb94-efd09e7d0864
# ╟─d846d3f2-8cfc-4fc5-a22b-d55762f88f45
# ╠═0b8f2bcc-0948-47ed-bd62-4a0ceeee9156
# ╠═1eccef57-f63b-4fba-a8f6-458652686fb1
# ╟─bb5988bf-21eb-4078-b8b7-5a8dedf2ff6e
# ╟─4562b2f4-3fcb-4c60-83b6-cafd6e4b3144
# ╠═7829afd0-2693-40dc-b2d5-c8da72e96454
# ╠═d6e89a4c-e59b-462d-8db6-0bc1a98ac381
# ╟─2240f8bc-5c0b-450a-b56f-2b53ca66bb03
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
