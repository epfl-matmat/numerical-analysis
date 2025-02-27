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

# ╔═╡ 46b46b8e-b388-44e1-b2d8-8d7cfdc3b475
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/05_Interpolation.pdf)
"""

# ╔═╡ 61e5ef66-a213-4b23-9406-9cc63a58104c
TableOfContents()

# ╔═╡ 345043c8-2e69-4a38-a3e9-c7a95f755f4c
md"""
# Interpolation

In this chapter we will return to one of the problems we already briefly discussed in the introductory lecture, namely:

!!! info "Definition: Interpolation problem"
    Suppose we are given data $(x_i, y_i)$ with $i = 1, 2, 3, \ldots n$,
    where the $x_i$ are all distinct. Find a function $p$
    such that $p(x_i) = y_i$.

The function $p$ is usually called the **interpolant** or **interpolating function**. 
"""

# ╔═╡ 263fa0f1-71c4-4794-9cc1-330ce50c282f
md"""
One can also consider the task of finding such an interpolating function $p$ as one simple example for a **data-driven method**:
- Let us assume the data $(x_i, y_i)$ originates from a complex statistical process (e.g. a lab experiment or an involved simulation)
- Having observed $n$ observed data points $(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)$ we thus want to obtain a cheaper **model** to predict a measurement $(x_{n+1}, y_{n+1})$ for a so far unseed $x_{n+1}$.
- By the means of interpolation we effectively **train** such a model, namely the interpolant $p$. Based on the interpolated $p$ we can **make a prediction** of the data point $(x_{n+1}, y_{n+1})$ as $(x_{n+1}, p(x_{n+1}))$.
"""

# ╔═╡ 0a08db0d-41ec-471c-a4fe-ca4205635c40
Foldable("Data-driven methods",
md"""
**Data-driven methods** refer to mathematical approaches, where one strives to come up with powerful predictive models from oberving only a few data points. Examples include methods from statistical learning such as Bayesian regression or neural networks.
""")

# ╔═╡ 0852f479-cf4d-4293-9a7a-0b7aa6e10d81
md"""
Of course we have not really discussed precisely *how* to obtain such an $p$. Indeed the interpolation condition $p(x_i) = y_i$ is rather weak and in practice and leaves considerable freedom:

- Show polynomial interpolant: $(@bind show_poly CheckBox(; default=true))
- Show piecewise linear interpolant: $(@bind show_pw CheckBox(; default=true))
- Show periodic interpolant: $(@bind show_periodic CheckBox(; default=true))
"""

# ╔═╡ 25258c1f-2685-4df8-a292-8708e013284a
md"""
Which of the above interpolants is the *best* is strongly dependent on the scientific context, for example: 
- If we know our data to have some periodicity, than the periodic interpolant (purple) may be better.
- If we need differentiability, the piecewise linear interpolant (green) may not be a good choice.

Therefore, in practical applications the quality of the fit usually depends on how well our method to construct $p$ agrees with the behaviour of the data itself.
"""

# ╔═╡ 6e46dbc2-7cd0-4fd7-a6bb-55355c53d861
md"""
The **typical ansatz** for an interpolation problem is 
 to take $p$ from a family of suitable functions
(e.g. polynomials, trigonometric functions etc),
which form a vector space.
If $φ_1, φ_2, \ldots φ_n$ is a basis for for this vector space,
then we can write $p$ as the linear combination
```math
\tag{1}
p(x) = \sum_{j=1}^n c_j φ_j(x)
```
of $n$ basis functions. Employing further the condition
that $p$ needs to pass through the data points, i.e. $p(x_i) = y_i$ for all $i$,
then we get $n$ equations
```math
p(x_i) = \sum_{j=1}^n c_j φ_j(x_i) = y_i \qquad i = 1, \ldots n
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

We start with one of the simplest forms of interpolation: Namely fitting polynomials through the given data.

Let us come back to our earlier problem of finding an interpolating function $p$
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

# ╔═╡ 944fa679-77d2-4fe5-9a79-ef63248105c7
md"""
### Error analysis

With the Lagrange polynomials we have a simple approach to find an interpolating polynomial. From a numerical analysis perspective this raises the question **how good a polynomial approximation is**.
Or to put it into the language of **data-driven methods**: Having obtained an interpolation model, how well does it generalise to unseen $x_n$ ?

To study this mathematically we consider the following setup:
 Assume we have a true function $f$ and we are allowed to observe $n$ samples $(x_n, f(x_n))$ from it. We want to understand the following questions:
- **Would the polynomial approximation converge** to $f$ as we observe more and more samples of $f$? 
- **How fast** is this convergence ?
- **What is our error** if we only observe $n$ samples $(x_n, f(x_n))$ and are not given any more information ?
"""

# ╔═╡ 055cb883-a140-4078-bedb-666f7e6259d0
md"""
In a mathematical language, we want to check whether
```math
p_n \to f \quad \text{as} \quad n\to \infty.
```
Provided this convergence is indeed the case, we are usually also interested in the **convergence rate**
--- similar to our discussion about fixed-point algorithms.
Clearly, the faster $p_n$ converges $f$ the better, since we need **less samples**
to obtain an **accurate reconstruction** of the original $f$.

"""

# ╔═╡ 5f1abde6-f4fc-4c5e-9a19-abf28ca3584d
md"""
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

# ╔═╡ 08fe0c39-f73d-46d0-b435-17abb3f60583
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

Note that in this theorem the error $\|f - p_n\|_\infty$ effectively measures how badly our polynomial interpolation model $p_n$ generalises for unseen datapoints $(x_{n+1}, f(x_{n+1}))$ with $x \in D$: If this error $\|f - p_n\|_\infty$ is small, $p_n$ is a very good statistical model for $f$. If this error is large, it is a rather inaccurate model.
"""

# ╔═╡ 5f63abce-f843-4c33-9db6-0770323b55ac
md"""
The **key conclusion of the previous theorem** is that
the error $\|f - p_n\|_\infty$ only decreases as $n$ increases
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

# ╔═╡ 890b203c-d380-45ae-bb08-739dd0f4a1da
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
"""

# ╔═╡ 25b82572-b27d-4f0b-9be9-323cd4e3ce7a
md"""
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
	yticks!(p, 10.0 .^ (-10:2:4))
	p
end

# ╔═╡ 479a234e-1ce6-456d-903a-048bbb3de65a
md"""

Notably Chebyshev nodes enjoy the following convergence result:

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

# ╔═╡ eaaf2227-1a19-4fbc-a5b4-45503e832280
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

As a graphical visualisation:
"""

# ╔═╡ d9e318e7-5ffc-4433-9b0c-f393948c10ef
begin
	ns = 3:13
	L_equi(n) = 2^(n+1) / (n * ℯ * log(n))
	L_Chebyshev(n) = 2π * log(n+1)
	plot(ns, L_equi.(ns), lw=2, mark=:o, label=L"$Λ_n$ for equally spaced nodes")
	plot!(ns, L_Chebyshev.(ns), lw=2, mark=:o, label=L"$Λ_n$ for Chebyshev nodes")
end

# ╔═╡ 64ef7e43-c07b-42de-9f71-89302433c178
md"""
Therefore for **Chebyshev nodes** the **condition number** grows only **logarithmically** with $n$, while for **equally spaced nodes** it grows **exponentially** !

Thus **Chebyshev nodes** do **not only** lead to **faster-converging polynomial interpolations**,
but also to notably **more stable answers**. As a result they are one of the standard ingredients in many numerical algorithms.

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
- Logarithm of noise amplitude `η_log_poly = ` $(@bind η_log_poly Slider(-3:0.25:-0.5; show_value=true, default=-2))
"""

# ╔═╡ 4a4540a3-b4cf-47fd-a615-a5280505333f
let
	fine = range(-1.1, 1.1; length=500)
	
	n_degree = n_nodes_poly - 1
	values_equispaced = range(-1, 1; length=n_nodes_poly)
	values_chebyshev  = [-cos(π*k/n_degree) for k in 0:n_degree]

	p = plot(; layout=grid(2, 2; heights=[0.75, 0.25]))
	for (i, values) in enumerate((values_equispaced, values_chebyshev))
		noise = [10^η_log_poly * (-1 + 2rand()) for _ in values]
		
		fsin_accurate = fsin.(values) #  = [fsin(x) for x in values]
		fsin_noise    = fsin.(values) + noise
		
		poly_accurate = Polynomials.fit(values, fsin_accurate)
		poly_noise    = Polynomials.fit(values, fsin_noise)

		plot!(p, fine, fsin.(fine);
		      label=L"function $f$", subplot=i, lw=2, ls=:dash,
		      ylim=(-2.75, 2), legend=:bottomright)
		plot!(p, fine, poly_accurate.(fine); label=L"accurate $p_n$", subplot=i, lw=2)
		plot!(p, fine, poly_noise.(fine); label=L"noisy $\tilde{p}_n$", subplot=i, lw=2)
		scatter!(p, values, fsin_noise; subplot=i, label="noisy data")

		sub_noise = i + 2
		plot!(p, fine, abs.(fsin.(fine) .- poly_noise.(fine)), label="noisy fit", yaxis=:log, xlim=(-1, 1), ylims=(10e-5, 1), legend=:topright, c=3, lw=1.5, ylabel="error", subplot=sub_noise)
		hline!(p, [10^η_log_poly], label=L"noise $\eta$", c=3, lw=1, ls=:dash, subplot=sub_noise)
		# plot!(p, fine, abs.(fsin.(fine) .- poly_accurate.(fine)), label="fit", c=2, lw=1.5, subplot=sub_noise)
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
*Note:* We will only discuss the high-level results of this part in the lecture.

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

# ╔═╡ 321c3303-824b-459f-b30c-f56bcb6d8a56
let
	x = [0.05648884843647836, 1.0015220577537893, 2.1104608983667728, 2.885496598146311, 4.127461807466144, 5.074047771027794, 6.116880303070477, 6.932437329389779, 7.951621443647719, 9.146205522650643, 10.027023988855284]
	y = [0.9495736137923013, -0.4274209000654566, -0.4317414334706035, 0.8845385764965122, -0.3616427771927499, -0.7418367372339307, 0.9012312093079629, 0.3034312370811605, -1.0285796376668908, 0.8420835485339082, 0.41544701770481224]

	f_poly = Polynomials.fit(x, y)
	f_pl   = pwlinear(x, y)
	
	V_cos  = reduce(hcat, [[cos((i-1) * x[j]) for i in 1:length(y)] for j in 1:length(x)])'
	c_cos  = V_cos \ y
	f_cos(x) = sum(c_cos[i] * cos((i-1) * x) for i in 1:length(c_cos))


	xfine = range(0, 20, length=400)

	p = plot([20, 10, 10, 20, 20], [-2, -2, 2, 2, -2]; fill=(0, :grey, 0.3), c=:grey, label="Extrapolative regime")

	
	p = scatter!(p, x, y, label="Data", ylims=(-2, 2), legend=:bottomright, xlims=(0, 20), c=1, markersize=6)
	if show_poly
		p = plot!(p, xfine, f_poly.(xfine), label="polynomial", lw=3, c=2)
	end
	if show_pw
		p = plot!(p, xfine, f_pl.(xfine), label="piecewise", lw=3, ls=:dot, c=3)
	end
	if show_periodic
		p = plot!(p, xfine, f_cos.(xfine), label="periodic", lw=3, ls=:dash, c=4)
	end
	p
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

# ╔═╡ d0a8b733-af40-41f7-a900-bf0ec6b17ed3
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
\left| f(x) - p_{1,h}(x) \right| \leq \frac{h^2}{8} \max_{x \in [x_{i}, x_{i+1}]} |f''(x)|
```
therefore
```math
\begin{aligned}
\|f - p_{1,h}\|_\infty
&= \max_{x\in [a,b]} |f(x) - p_{1,h}(x)| \\
&= \max_{i=1,\ldots,n} \max_{x\in[x_{i},x_{i+1}]} |f(x) - p_{1,h}(x)| \\
&\leq \max_{i=1,\ldots,n} \frac{h^2}{8} \max_{x \in [x_{i}, x_{i+1}]} |f''(x)| \\
&= \frac{h^2}{8} \| f'' \|_\infty
\end{aligned}.
```
"""

# ╔═╡ 2c4256d0-bcd7-42b4-ad85-aedea7063aae
md"""
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
- The **condition number** of piecewise linear interpolation is $1$ **independent
  of the data**.
- As a result **small errors in the input values** $x_i$ (e.g. from measurements)
  **only introduce small perturbations** in the polynomial:
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
## Optional: Spline interpolation

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

# ╔═╡ 923ed582-e187-4f58-b945-93973bdd358c
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
"""

# ╔═╡ fe1f46df-1b42-4380-80bd-e38d5eb37e00
md"""
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

For examples as well as more details on the computation of splines see
[chapter 5.3](https://tobydriscoll.net/fnc-julia/localapprox/splines.html)
of Driscoll, Brown: *Fundamentals of Numerical Computation*.
"""

# ╔═╡ d4226809-9f2a-4651-9778-8f8f66d03720
md"""
### Optional: Error analysis

For cubic splines the following convergence result is known:

!!! info "Theorem 5"
    Let $f \in C^4([a, b])$ be $4$ times differentiable and let
    $a < x_1 < x_2 < \cdots < x_{n+1} = b$ be equispaced nodal points
    in $[a,b]$. The cubic spline $s_{3,h}$ interpolating the data $(x_i, f(x_i))$
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
	p = plot(h, maxerror; label="error", title="Convergence spline interpolation",
	         xlabel=L"h", xflip=true, xscale=:log10, ylabel=L"|| f-s_{3,h} ||_\infty",
	         yscale=:log10, mark=:o, legend=:topright, lw=2)

	# Generate guiding slope
	order4 = (h ./ h[1]) .^ 4
	plot!(p, h, order4, label=L"O(h^4)", ls=:dash, lw=2)

	p
end

# ╔═╡ b80bdd34-fc55-4e68-a10a-f3e4cff6451a
md"""
## Regression and curve fitting

Consider again the case where we have acquired $n$ data points
$(x_i, \tilde{y}_i)$, $i=1, 2, \ldots, n$ where 
```math
\tilde{y}_i = f(x_i) + ε_i
```
and the **measurement noice $ε_i$ is substantial**.
If we were to use an *interpolation technique*
--- as discussed in the previous sections ---
our goal would be to obtain a (potentially piecewise)
polynomial $p$ such that $p(x_i) = \tilde{y}_i$.
However, since $\tilde{y}_i$ is only a *noisy* observation
it makes **little sense to force the fitted polynomial
to go through the data exactly**.

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
"""

# ╔═╡ 36c8221e-8a40-4d2f-abdb-f7ed96b89d47
md"""
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

**Least-squares regression** is usually much **better suited to capture trends**
in data than polynomial interpolation.
To see this consider the following setting, where we fit
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
    Let $\mathbf{A} \in \mathbb{R}^{k\times \ell}$ and $\mathbf{b} \in \mathbb{R}^k$
    with $k > \ell$.
    If $\mathbf{x} \in \mathbb{R}^\ell$ satisfies $\textbf{A}^T (\textbf{A}\textbf{x} - \textbf{b}) = \textbf{0}$ then $\textbf{x}$ solves the least-squares problem, i.e.
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

# ╔═╡ d3fdad97-275d-49d8-a4aa-b1f6438a01b8
md"""
### Error analysis

In the general setting the analysis of least-squares problems is rather involved. 
Here, we will restrict ourselves to exploring the parameter space a little using an interactive visualisation.
"""

# ╔═╡ 14bb1d8e-dd73-4795-9bb9-0213e83020d9
Foldable("Some interesting experiments with the visualisation below.",
md"""
**Polynomial interpolation regime: Sanity check with our results from before:**
  - `n = 20` and `η_log` small: Increase degree `m` slowly up to `20`: you should observe Runge's phaenomenon for large `m` as before.
  - Keep `n = 20` and `m = 20` and increase the noise `η_log`: The error drastically increases (well beyond the original noise level `η_log`); the fit is very unstable as we are essentially doing polynomial interpolation with equispaced points.
  - In general for `m = n` (polynomial interpolation) the situation is similar, try for example `n=15` and `m=15` with varying noise.

**The regime $n \gg m$ of least-squares regression:**
  - Set `m=15`, `m=15` and `η_log = -1`, so polynomial interpolation with large noise. The errors are fairly large. Now increase `n` slowly to `50`: All of a sudden the errors are in control and basically never exceed `η_log`. Play with `η_log` by making it even larger: You should see the error remains well below `η_log` in all cases.
  - Set `m=40` and `η_log = -0.5` and slide the degree `m` between `20` and `10`. At `m=10` the error is notacibly smaller than at `m=20`. This suggests that somehow the ratio $\frac{m}{n}$ has some influence to what extend measurement noise $η$ translates to an error in the outcome of polynomial regression.
  - Set `n=20` and `η_log = -1.5` and slide the polynomial degree. We realise that there is sweet spot aronud `m=10` when the error is overall smallest. At too low polynomial degrees $m$ the model is not rich enough to approximate the sine, while at too large degrees $m$ the ratio $\frac{m}{n}$ gets large we get into the regime of a polynomial interpolation. Therefore the numerical noise begins to amplify, resulting in larger and larger errors as we keep increasing $m$.
""")

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

# ╔═╡ c2431601-afc7-4b37-9113-5ae85d4e5549
md"""
To close off we state  some key results,
which illustrate the key properties of the least-squares problem.
We recall the setting: We are given $n$ data points 
$(x_i, \tilde{y}_i)$, $i=1, 2, \ldots, n$ where 
```math
\tilde{y}_i = f(x_i) + ε_i
```
are the noisy observations of a function $f$.
Further we assume that $\mathbb{E}(ε_i) = 0$, i.e. the noise follows a distribution
with mean zero and $σ = \sqrt{\mathbb{V}\text{ar}(ε_i)}$ is the square root of the variance
of the noise.

- In the **absence of measurement noise**, i.e. $σ=0$, **least-squares can recover polynomials exactly**. Imagine $f$ to be a polynomial of degree $m \ll n$. In this case the polynomial $f$ is the *only* $m$-th degree polynomial to give a zero least-squares error, i.e.
  ```math
  0 = \sum_{i=1}^{n} |p_m^\text{LS}(x_i) - f(x_i)|^2 \text{ with $p_m^\text{LS}\in\mathbb{P}_m$}
  \quad \Leftrightarrow \quad p_m^\text{LS} = f.
  ```
  and our result is exact.
- With **non-zero measurement errors** the approximation error $\|f - p_m^\text{LS}\|_\infty$ is proportional to $\sqrt{\frac{m+1}{n}}σ$. Thus while every measurement has been perturbed by an error $ε_i$ of order around $σ$, the **approximation error** is much smaller of order $\sqrt{\frac{m+1}{n}}σ$ and moreover **becomes smaller** as the **number of samples increases**.

  However, if the **degree of the fitting polynomial $m$ becomes too large** --- say comparable to $n$ ---, then the least-squares polynomial becomes similar to the
  interpolating polynomial. This we saw to become unstable as $m$ gets large.
  Therefore the **quality of the least-squares approximation** does **deteriorate
  when $m \to n$**.
"""

# ╔═╡ ebf6b5f4-b963-4296-9255-47e426d2d12d
md"""
## Summary
- **Interpolation** and **least-squares regression** techniques is one common approach to **extract a model** $p$ **from** $n$ **observed data** points $(x_i, y_i)$ (with $i=1,2,\ldots,n$). Based on this model one can make predictions about unseen $x_{n+1}$ namely as the points $(x_{n+1}, p(x_{n+1}))$.

- Polynomial interpolation on **equally spaced data points** leads to Runge's phaenomenon as the polynomial degree $m$ is growing.
  * Moreover this problem is **ill-conditioned**, i.e. extremely susceptible to numerical or experimental noise in the training data $(x_i, y_i)$.

- **One solution**: Keep the polynomial degree $m$ low, for example:
  * Use piecewise polynomial (e.g. **piecewise linear**) interpolation techniques.
  * Use **more observations $n$ than $m$**, i.e. perform **least-squares regression** with $n \gg m$
  * This generally leads to methods with **algebraic convergence** when approximating a smooth function $f$.
  * Moreover such problems are generally **well-conditioned**. E.g. for piecewise linear interpolation the condition number is $1$ independent of the employed data --- the best-possible value.

- The **other solution** is to use **non-equally spaced points**:
  * The typical approach are **Chebyshev nodes**
  * These lead to **exponential convergence**
"""

# ╔═╡ 2240f8bc-5c0b-450a-b56f-2b53ca66bb03
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 500)
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
LaTeXStrings = "~1.4.0"
Plots = "~1.40.0"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.55"
Polynomials = "~4.0.17"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "09ea56855cccbf151b4ad3576b58d9161b143e43"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

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
git-tree-sha1 = "4bf4b400a8234cff0f177da4a160a90296159ce9"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.41"

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
git-tree-sha1 = "cd714447457c660382fe634710fb56eb255ee42e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.6"

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

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "0973615c3239b1b0d173b77befdada6deb5aa9d8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.17"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

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

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

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
git-tree-sha1 = "ee6f41aac16f6c9a8cab34e2f7a200418b1cc1e3"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+0"

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
# ╟─46b46b8e-b388-44e1-b2d8-8d7cfdc3b475
# ╠═679517fb-6090-4259-9b93-571ded17ab88
# ╟─61e5ef66-a213-4b23-9406-9cc63a58104c
# ╟─345043c8-2e69-4a38-a3e9-c7a95f755f4c
# ╟─263fa0f1-71c4-4794-9cc1-330ce50c282f
# ╟─0a08db0d-41ec-471c-a4fe-ca4205635c40
# ╟─0852f479-cf4d-4293-9a7a-0b7aa6e10d81
# ╟─321c3303-824b-459f-b30c-f56bcb6d8a56
# ╟─25258c1f-2685-4df8-a292-8708e013284a
# ╟─6e46dbc2-7cd0-4fd7-a6bb-55355c53d861
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
# ╟─944fa679-77d2-4fe5-9a79-ef63248105c7
# ╟─055cb883-a140-4078-bedb-666f7e6259d0
# ╟─5f1abde6-f4fc-4c5e-9a19-abf28ca3584d
# ╠═b1773d1c-4d3c-4f4b-805d-f188623c9571
# ╟─e567287b-9244-4b55-9f5b-2e3bf583e46a
# ╟─31126101-685d-464d-be4a-6233d3803780
# ╟─af534a92-31d4-4125-bb7c-3eb2b3430852
# ╟─08fe0c39-f73d-46d0-b435-17abb3f60583
# ╟─5f63abce-f843-4c33-9db6-0770323b55ac
# ╟─de28f64c-f243-4cb5-a1e6-ae8562cce5bc
# ╟─08d5cb8a-3bf5-4d40-a8d6-d22318e37712
# ╟─c5b33503-973b-468f-95de-acd6a7605a99
# ╟─890b203c-d380-45ae-bb08-739dd0f4a1da
# ╟─25b82572-b27d-4f0b-9be9-323cd4e3ce7a
# ╟─c38b9e48-98bb-4b9c-acc4-7375bbd39ade
# ╟─479a234e-1ce6-456d-903a-048bbb3de65a
# ╟─d4cf71ef-576d-4900-9608-475dbd4d933a
# ╟─56685887-7866-446c-acdb-2c20bd11d4cd
# ╟─1d61d925-1ceb-439b-aff6-c5e8483fdc31
# ╟─eaaf2227-1a19-4fbc-a5b4-45503e832280
# ╟─d9e318e7-5ffc-4433-9b0c-f393948c10ef
# ╟─64ef7e43-c07b-42de-9f71-89302433c178
# ╟─d5de3b11-7781-4100-8a63-d26426685bbc
# ╟─4a4540a3-b4cf-47fd-a615-a5280505333f
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
# ╟─d0a8b733-af40-41f7-a900-bf0ec6b17ed3
# ╟─2c4256d0-bcd7-42b4-ad85-aedea7063aae
# ╟─7e317807-d3ee-4197-91fb-9fc03f7297e8
# ╠═eb96f944-74ac-4cc0-9322-825df83f34f2
# ╟─9d175a21-71e0-4df8-bbd5-f4eedbeb1384
# ╟─3eaddcc1-bddb-45a5-9e0d-0961dba5583c
# ╟─2096258a-efd9-4de1-8f56-c02295407c0b
# ╟─e2b32c49-9c8e-4d32-8bf6-a0517b7caa8a
# ╟─cf296fbf-c900-485c-8cfe-672f02ad3b69
# ╟─22c33635-dd92-4a8d-9f99-74b222a44f2e
# ╟─923ed582-e187-4f58-b945-93973bdd358c
# ╟─fe1f46df-1b42-4380-80bd-e38d5eb37e00
# ╠═5987bb0f-c784-42ab-9c34-d6ff76cef2ea
# ╟─fd469b84-9093-4b92-ac3e-c1847db3627f
# ╟─33dd12c8-f3d5-4238-986f-df58639a6f3a
# ╟─d4226809-9f2a-4651-9778-8f8f66d03720
# ╠═a806cafa-8ee9-4dd6-9b90-100a831c009d
# ╟─b80bdd34-fc55-4e68-a10a-f3e4cff6451a
# ╟─36c8221e-8a40-4d2f-abdb-f7ed96b89d47
# ╠═5c5b212d-0451-4334-b5f7-273640414e86
# ╟─7bd4720f-3ec8-459a-91be-17d5f2acaa5c
# ╟─66cfe319-c19e-478f-8f36-645fa1ff8b3b
# ╟─f6916b49-4fce-429c-9f55-d9a51dc46937
# ╟─2122f3c9-a901-4488-bd1b-4bd7e49b7c8a
# ╟─30899dd5-b307-415d-bb94-efd09e7d0864
# ╟─d846d3f2-8cfc-4fc5-a22b-d55762f88f45
# ╠═0b8f2bcc-0948-47ed-bd62-4a0ceeee9156
# ╟─d3fdad97-275d-49d8-a4aa-b1f6438a01b8
# ╟─14bb1d8e-dd73-4795-9bb9-0213e83020d9
# ╟─4562b2f4-3fcb-4c60-83b6-cafd6e4b3144
# ╠═7829afd0-2693-40dc-b2d5-c8da72e96454
# ╟─c2431601-afc7-4b37-9113-5ae85d4e5549
# ╟─ebf6b5f4-b963-4296-9255-47e426d2d12d
# ╟─2240f8bc-5c0b-450a-b56f-2b53ca66bb03
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
