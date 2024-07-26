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

# ╔═╡ 4949225a-ccf2-11ee-299b-9b834eb6bd42
begin
	using LinearAlgebra
	using PlutoUI
	using PlutoTeachingTools
	using Plots
	using LaTeXStrings
	using HypertextLiteral
end

# ╔═╡ 13298dc4-9800-476d-9474-182359a7671b
TableOfContents()

# ╔═╡ 377d9d74-2826-4c53-bc9a-8bfe5a087a84
md"""
# Eigenvalue problems

In the previous notebooks we already noted a connection between the eigenspectrum of the iteration matrix and the convergence properties of fixed-point problems or the convergence properties of solving linear systems. 

In this notebook we will discuss some simple iterative methods for actually computing  eigenpairs. However, the topic is vast and we will only scratch the surface. 
Readers interested in a more in-depth treatment
of eigenvalue problems are encouraged to attend the master class
[MATH-500: Error control in scientific modelling](https://staging-edu.epfl.ch/coursebook/en/error-control-in-scientific-modelling-MATH-500).

A more comprehensive treatment of the topic can also be found in the book [Numerical Methods for Large Eigenvalue Problems](https://epubs.siam.org/doi/book/10.1137/1.9781611970739) by Youssef Saad as well as the [Lecture notes on Large Scale Eigenvalue Problems](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf) by Peter Arbenz.
"""

# ╔═╡ 71d74b6a-96b2-4855-b20f-59b00c8f560b
md"""
## Power iteration

We start with a simple question: What happens if we apply a matrix multiple times ?
Here, we choose the matrix
"""

# ╔═╡ 6db511ac-a217-49e2-a508-2590b64ea636
A = [ 0.5  1/3;
      0.5  2/3]

# ╔═╡ 814e602d-1afa-44ba-bb25-ed32dd66e2b9
md"""
and a random 2-element vector:
"""

# ╔═╡ 58a88c32-cdf6-4cb9-af08-0ef9a186fe1c
md"Applying the matrix once seems to be very innocent:"

# ╔═╡ 7da94235-399e-4c13-8a02-25d8fa795e17
md"""
But if we appy many times we start to see something ...
"""

# ╔═╡ 5a4038ac-1905-42e5-88ee-ced49c5375ad
@bind dummy Button("Regenerate random matrix and rerun experiment")

# ╔═╡ 0dabc016-9cb8-4d1b-bd44-fb9493b8cf28
begin
	dummy  # For rerun to work
	x = randn(2)
end

# ╔═╡ aa47c74d-2202-413a-9cad-17cad20832e1
x

# ╔═╡ 5d990087-901a-4457-a9ec-ba3b2d60f470
A * x

# ╔═╡ 272e7bce-e05e-412f-b2f5-2742d33d0b6f
begin
	for j in 1:6
		x = A * x
	end
	maxerror = maximum(abs, (x - A * x))
	@show maxerror
	[x A * x]
end

# ╔═╡ 076bbe2b-130a-44e5-9218-6fdcc1998cf5
TODO("Graphical visualisation")

# ╔═╡ db0c875b-9955-4bc1-bdb1-ec2e1a44942d
md"""Note how the iterations stabilise, i.e. that $\textbf{x}$ and $\textbf{A} \textbf{x}$ start to be alike. In other words we seem to achieve $\mathbf{A} \mathbf{x} = \mathbf{x}$, which is nothing else than saying that $\mathbf{x}$ is an eigenvector of $\mathbf{A}$ with eigenvalue $1$.
"""

# ╔═╡ 4d1eb1b1-7634-4f38-a564-98f64923e377
TODO(md"Warning box about eigenvalue ordering: Note explicitly the order of Julia and the orderings we use here. Note that we use different eigenvalue orderings depending on the algorithm which is discussed and the context that one focuses on. Note in partiuclar that $λ_1$ and $λ_2$ throughout these notes may mean different things !")

# ╔═╡ 73c62d23-ac2a-48be-b3c7-0d51ffce773c
md"""
Note, that in this developent we have cheated a little, namely we allowed ourself to 
normalise the columns (the `M = M ./ sum(M; dims=1)` above).

Let us understand what happened in this example in detail now. We consider the case $\mathbf{A} \in \mathbb{R}^{n \times n}$ diagonalisable and let further
```math
\tag{1}
|λ_1| > |λ_2| ≥ |λ_3| ≥ \cdots ≥ |λ_n|
```
be its eigenvalues with corresponding eigenvectors $\mathbf v_1$, $\mathbf v_2$, $\ldots, \mathbf v_n$,
which we collect column-wise in a unitary matrix $\mathbf{V}$, i.e.
```math
\mathbf{V} = \begin{pmatrix} \mathbf v_1 & \mathbf v_2 & \cdots & \mathbf v_n \\
				\downarrow & \downarrow & \cdots & \downarrow
             \end{pmatrix}
```
Note that $|λ_1| > |λ_2|$, such that $λ_1$ is non-degenerate and the absolutely largest eigenvalue. We call it the **dominant eigenvalue** of $\mathbf{A}$.
If $k > 0$ is a positive integer, then we have
```math
\mathbf{A}^k = \left(\mathbf{V} \, \mathbf{D}\, \mathbf{V}^{-1}\right)^k = \mathbf{V} \, \mathbf{D}^k \,\mathbf{V}^{-1} 
```
and we can expand any random starting vector $\mathbf{x}^{(1)}$ in terms of the eigenbasis
of $\mathbf A$, i.e.
```math
\mathbf{x}^{(1)} = \sum_{i=1}^n z_i \mathbf v_i = \mathbf{V} \mathbf{z},
```
where $\mathbf{z}$ is the vector collecting the coefficients $z_i$.
Notably this implies $\mathbf z = \mathbf{V}^{-1} \mathbf{x}^{(1)}$.
"""

# ╔═╡ 16d6559c-977b-4fd9-aecb-2b6f61290f04
md"""
We now consider applying $\mathbf{A}$ a $k$ number of times to $\mathbf{x}^{(1)}$,
which yields
```math
\tag{2}
\begin{aligned}
\mathbf{A}^k \mathbf{x}^{(1)} &= \mathbf{V} \, \mathbf{D}^k \,\underbrace{\mathbf{V}^{-1} \mathbf{x}^{(1)}}_{= \mathbf{z}} = \mathbf{V} \begin{pmatrix}
λ_1^k z_1 \\
λ_2^k z_2 \\
\vdots\\
λ_n^k z_n
\end{pmatrix} \\ 
&= λ_1^k \left( z_1 \mathbf v_1 + \left( \frac{λ_2}{λ_1} \right)^k z_2 \mathbf v_2 + \cdots + \left( \frac{λ_n}{λ_1} \right)^k z_n \, \mathbf v_n   \right).
\end{aligned}
```
Therefore if $z_1 \neq 0$ then
```math
\tag{3}
\left\| \frac{\mathbf{A}^k \mathbf{x}^{(1)}}{λ_1^k} - z_1 \textbf v_1 \right\| ≤ 
\left|\frac{λ_2}{λ_1}\right|^k \, |z_2| \, \|\textbf v_2\| +
\cdots
+ \left|\frac{λ_n}{λ_1}\right|^k \, |z_n| \, \|\textbf v_n\|
```
Now, each of the terms $\left|\frac{λ_j}{λ_1}\right|^k$
for $j = 2, \ldots, n$ goes to zero for $k \to \infty$ due to the
descending eigenvalue ordering of Equation (1).
Therefore overall
```math
\left\| \frac{\mathbf{A}^k \mathbf{x}^{(1)}}{λ_1^k} - z_1 \textbf v_1 \right\|
\to 0 \qquad \text{as $k \to \infty$}.
```
In other words $\frac{1}{λ_1^k} \mathbf{A}^k \mathbf{x}^{(1)}$ eventually becomes a multiple of the eigenvector associated with the dominant eigenvalue.
"""

# ╔═╡ 139ad60c-e04f-4dd6-b18c-a6de140a0154
TODO(md"
First show that generalising to a 2x2 matrix without the rescaling trick fails because the entries in the iterated vector $x$ get larger and larger.

Then suggests to rescale the entries somehow, one natural idea is to simply use the largest by magnitude entry. Show that this fixes the issue.

Then motivate why β is a good estimate for the eigenvalue, show the algorithm, show its convergence and then show its convergence proof.
")

# ╔═╡ 1f512d18-1f9c-4730-9676-ca333a7531de
md"""
### Power iteration algorithm

In order to turn this idea into an algorithm there are two aspects missing:
-   $λ_1^k$ becomes either extremely small (for $λ_1 < 1$) or extremely large (for $λ_1 > 1$) if $k \to \infty$. As a result representing the matrix entries of $\mathbf A^k \mathbf{x}$ and performing the associated matrix-vector products numerically accurately in finite-precision floating-point arithmetic is difficult. Some form of *normalisation* at each step of our iteration techinque is thus required.
- For the normalisation we cannot simply divide by $λ_1$ itself, since this eigenvalue is unknown to us.

But other choices of normalisation exist. In this case we will take the infinity norm:
```math
\|x\|_\infty = \max_{i=1,\ldots n} |x_i|,
```
which simply selects the largest absolute entry of the vector.
Based on this idea we obtain the algorithm

!!! info "Algorithm 1: Power iterations"
    Given a diagonalisable matrix $\mathbf A \in \mathbb R^{n\times n}$
    and an initial guess $\mathbf x^{(1)} \in \mathbb R^n$ we iterate
    for $k = 1, 2, \ldots$:
    1. Set $\mathbf y^{(k)} = \mathbf A \, \mathbf x^{(k)}$
    2. Find the index $m$ such that $\left|y_m^{(k)}\right|$ is the largest (by magnitude) element of $\textbf{y}^{(k)}$, i.e.
       ```math
       m = \textrm{argmax}_i \left|y_i^{(k)}\right|
       ```
    3. Set $α^{(k)} = \frac{1}{y^{(k)}_m}$ and $β^{(k)} = \frac{y^{(k)}_m}{x^{(k)}_m}$.
    4. Set $\mathbf x^{(k+1)} = α^{(k)} \mathbf y^{(k)}$

Note for this algorithm we can write
```math
\tag{4}
\textbf{x}^{(k)} = α^{(1)} \cdot α^{(2)} \cdot \, \cdots \, \cdot \, α^{(k)} \cdot \, \textbf{A}^k \textbf x^{(1)}
= \prod_{i=1}^k α^{(i)} \cdot \mathbf A^k \mathbf x^{(1)}
```
and by construction we always have $\left\| \textbf{x}^{(k)} \right\|_\infty = 1$. Algorithm 1 thus only presents a minor modification over the scheme we discussed in the previous section.
"""

# ╔═╡ 2073afcb-7e3c-4d3f-924e-cc768681a8fa
md"""
Finally we note that if $\textbf{x}^{(k)}$ is almost an eigenvector of the dominant eigenvalue, then $\mathbf A\, \mathbf{x}^{(k)}$ is almost $\lambda_1 \mathbf{x}^{(k)}$
and thus one would assume $β^{(k)} = \frac{y^{(k)}_m}{x^{(k)}_m}$ to be a reasonable eigenvalue estimate. Indeed note that
```math
\begin{aligned}
β^{(k)} &= \frac{y^{(k)}_m}{x^{(k)}_m}
= \frac{\left(\mathbf{A}\, \mathbf{x}^{(k)}\right)_m}{x^{(k)}_m}\\
&\stackrel{(4)}{=} \frac{
\left( \prod_{i=1}^k α^{(i)} \cdot \mathbf A^{k+1} \mathbf x^{(1)} \right)_\textcolor{red}{m}
}{
\left( \prod_{i=1}^k α^{(i)} \cdot \mathbf A^k \mathbf x^{(1)} \right)_\textcolor{red}{m}
}\\
&\stackrel{(2)}{=}
\frac{
λ_1^{k+1} \left( z_1 \mathbf v_1 + \left( \frac{λ_2}{λ_1} \right)^{k+1} z_2 \mathbf v_2 + \cdots + \left( \frac{λ_n}{λ_1} \right)^{k+1} z_n \, \mathbf v_n   \right)_\textcolor{red}{m}
}{
λ_1^{k} \left( z_1 \mathbf v_1 + \left( \frac{λ_2}{λ_1} \right)^{k} z_2 \mathbf v_2 + \cdots + \left( \frac{λ_n}{λ_1} \right)^k z_n \, \mathbf v_n   \right)_\textcolor{red}{m}
}\\
&= 
\frac{
λ_1 \left( b_1 + \left( \frac{λ_2}{λ_1} \right)^{k+1} b_2 + \cdots + \left( \frac{λ_n}{λ_1} \right)^{k+1} b_n   \right)
}{
b_1 + \left( \frac{λ_2}{λ_1} \right)^{k} b_2 + \cdots + \left( \frac{λ_n}{λ_1} \right)^k b_n
}
\end{aligned}
```
where $_m$ denotes that we take the $m$-th element of the vector in the brackets
and we further defined $b_i = z_i \left(\mathbf{v}_i\right)_m$, that is $z_i$ multiplied by the $m$-th element of the eigenvector $\mathbf{v}_i$.
Again we assume $z_1 \neq 0$, which implies $b_1 \neq 0$.
Under this assumption
```math
\tag{5}
β^{(k)} = \frac{y^{(k)}_m}{x^{(k)}_m}=
λ_1 \, \frac{
1 + \left( \frac{λ_2}{λ_1} \right)^{k+1} \frac{b_2}{b_1} + \cdots + \left( \frac{λ_n}{λ_1} \right)^{k+1} \frac{b_n}{b_1}
}{
1 + \left( \frac{λ_2}{λ_1} \right)^{k} \frac{b_2}{b_1} + \cdots + \left( \frac{λ_n}{λ_1} \right)^k \frac{b_n}{b_1}
}
\to λ_1 \qquad \text{as $k\to\infty$}.
```
Keeping in mind that $\frac{λ_i}{λ_1} < 1$ for all $2 ≤ i ≤ n$ we indeed
observe $β^{(k)}$ to converge to $λ_1$.

An implementation of this power method algorithm in:
"""

# ╔═╡ 01186225-602f-4637-99f2-0a6dd569a703
function power_method(A, x; maxiter=100)
    n = size(A, 1)
    x = normalize(x, Inf)  # Normalise initial guess

	history = Float64[]  # Record a history of all βs (estimates of eigenvalue)
    for k in 1:maxiter
        y = A * x
        m = argmax(abs.(y))
		α = 1 / y[m]
		β = y[m] / x[m]
		push!(history, β)
        x = α * y
    end

	(; x, λ=last(history), history)
end

# ╔═╡ 25f2b565-cc6a-4919-ac82-37d6dd62ba16
md"""
Let's try this on a lower-triangular matrix as a test. Note that the eigenvalues of a lower-triangular matrix are exactly the diagonal entries. We set
"""

# ╔═╡ 22b9c8b5-02cd-481e-aca7-bb2fe4e85c3c
begin
	λref = [1, -0.75, 0.6, -0.4, 0]
	# Make a triangular matrix with eigenvalues on the diagonal.
	M = triu(ones(5,5),1) + diagm(λref)
end

# ╔═╡ 81a6721f-f7e5-42c9-b389-518a4c458ae9
md"""
By construction the largest eigenvalue of this matrix is $1$. So let's track convergence:
"""

# ╔═╡ 033b6a9b-2187-4053-bc69-bbc1f03494fa
begin
	λlargest = 1.0  # Largest eigenvalue of M (by construction)

	x1 = ones(size(M, 2))
	results = power_method(M, x1; maxiter=70)
	error = abs.(results.history .- λlargest)
	plot(error; mark=:o, label="", yaxis=:log,
	     title="Convergence of power iteration",
	     ylabel=L"|β^{(k)} - λ_{\textrm{ref}}|", xlabel=L"k")
end

# ╔═╡ cc3495af-e53c-4ee8-8ec0-f4b597e17f53
md"""
!!! note "Remark on z₁ ≠ 0"
    In the theoretical discussions so far we always assumed $z₁ ≠ 0$,
    i.e. that the initial guess $\mathbf{x}^{(1)}$ has a non-zero component along
    the eigenvector $\mathbf v_1$ of the dominant eigenvalue. This may not
    always be the case. However even if $z₁ ≠ 0$,
    computing $\mathbf{A} \mathbf{x}^{(1)}$
    in finite-precision floating-point arithmetic is associated with
    a small error, which means that 
    which implies that $\mathbf{A} \mathbf{x}^{(1)}$ and thus
    $\mathbf{x}^{(2)}$ does usually
    *not* have a zero component along $\mathbf v_1$,
    such that in practical calculations assuming $z₁ ≠ 0$
    is not a restriction.
"""

# ╔═╡ c1b01e10-235c-4424-8524-0836a07106ff
md"""
### Convergence of power method

The above plot already suggests a **linear convergence** towards
the exact eigenvalue. Indeed a detailed analysis shows that the convergence rate
can be computed as
```math
\tag{6}
r_{\textrm{power}} = \lim_{k\to\infty} \frac{|β^{(k+1)} - λ_1|}{|β^{(k)} - λ_1|}
= \left|\frac{λ_2}{λ_1}\right|.
```
"""

# ╔═╡ 38015a4f-17af-47c5-a5b6-5e19c429c841
details("Detailed derivation",
md"""
To show the rate $r_{\textrm{power}}$ we revisit equation (5)
```math
β^{(k)} = \frac{y^{(k)}_m}{x^{(k)}_m}=
λ_1 \frac{
1 + \left( \frac{λ_2}{λ_1} \right)^{k+1} \frac{b_2}{b_1} + \cdots + \left( \frac{λ_n}{λ_1} \right)^{k+1} \frac{b_n}{b_1}
}{
1 + \left( \frac{λ_2}{λ_1} \right)^{k} \frac{b_2}{b_1} + \cdots + \left( \frac{λ_n}{λ_1} \right)^k \frac{b_n}{b_1}
}
```
and introduce the shorthand notation
```math
B_i = b_i/b_1 \qquad\text{and}\qquad
r_i = λ_i/λ_1
```
such that
```math
β^{(k)} =
λ_1 \frac{
1 +  r_2^{k+1} B_2 + \cdots + r_n^{k+1} B_n
}{
1 +  r_2^{k} B_2 + \cdots + r_n^k B_n.
}
```
We will make the *additional assumption* $|λ_2| < |λ_3|$,
such that $|λ_3 / λ_2| < 1$. This is not strictly speaking neccessary,
but it allows to simplify our developments,
since the common expression
```math
1 + r_2^{k} B_2 + \cdots + r_n^{k} B_n
= 1 + r_2^k \left[
B_2 + \left( \frac{λ_3}{λ_2} \right)^k B_3
+ \cdots + \left( \frac{λ_n}{λ_2} \right)^k B_n
\right]
```
clearly approaches $1 + r_2^k B_2$ as $k \to \infty$.
Now noting for $k \to \infty$ that
```math
\frac{1}{1 + r_2^k B_2} = 1 - r_2^k B_2 + O(r_2^{2k})
```
we can conclude in this limit
```math
β^{(k)} - λ_1 = λ_1 \left(1 +  r_2^{k+1} B_2\right)
\left( 1 - r_2^k B_2 + O(r_2^{2k}) \right) - λ_1
= (r_2 - 1) r_2^k B_2 + O(r_2^{2k}).
```
This is indeed linear convergence with convergence rate
```math
\tag{6}
\frac{|β^{(k+1)} - λ_1|}{|β^{(k)} - λ_1|}
\to \frac{|(r_2 - 1) r_2^{k+1} B_2|}{|(r_2 - 1) r_2^{k} B_2|} = |r_2| = \left|\frac{λ_2}{λ_1}\right|
\quad \text{as $k \to \infty$}
```
""")

# ╔═╡ d6a56c08-374a-4665-9968-d538fe2c4c59
md"""
So the smaller the ratio between $|λ_2|$ and $|λ_1|$,
the faster the convergence.
"""

# ╔═╡ a8dfd1bc-5146-4218-9fad-0676c7406f16
md"and we will put these on the diagonal of a triangular matrix, such that the eigenvalues of the resulting matrix $M$ are exactly the $λ_δ$:"

# ╔═╡ 385b07d6-15ef-40ee-8d2e-b3b6505a4171
md"""
We check this: We introduce a slider to tune the $λ_2$ of equation (6):
- `λ₂ = `$(@bind λ₂ Slider(-0.6:-0.05:-0.95; show_value=true, default=-0.9) )
"""

# ╔═╡ 79afdc06-a39e-4187-b91d-0af65a62613b
md"Let us take a case where $λ_1 = 1.0$. Therefore the size of $λ_2$ determines the rate of convergence. Let's take $λ_2$ = $(round(λ₂; sigdigits=2))."

# ╔═╡ 59b21a54-5cac-4e7f-b800-89c1f5f2d01f
begin
	λ₁ = 1
	λ_δ = [λ₁, λ₂, 0.6, -0.4, 0]   # The reference eigenvalues we will use
end

# ╔═╡ e5bf05e5-b064-4e43-b64b-22e5e45ddcd0
M_δ = triu(ones(5, 5), 1) + diagm(λ_δ)

# ╔═╡ 5d4697b4-2f3e-452b-afb8-78443efa1f25
let
	λlargest = maximum(abs.(λref))  # Largest eigenvalue of M_δ (by construction)

	x1 = ones(size(M_δ, 2))
	results = power_method(M_δ, x1; maxiter=40)
	error = abs.(results.history .- λlargest)
	p = plot(error; mark=:o, label="", yaxis=:log, ylims=(1e-10, 3),
	     title="Convergence of power iteration",
	     ylabel=L"|β^{(k)} - λ_{\textrm{ref}}|",
	     xlabel=L"k", legend=:bottomleft, lw=2)

	r_power = abs(λ₂ / λ₁)	
	plot!(p, k -> r_power^k; ls=:dash, label=L"expected rate $| λ_2 / λ_1|$", lw=2)
end

# ╔═╡ 0e4a443b-e99d-436f-8151-d6cfed9b773b
md"""
## Spectral transformations

The above procedure provides us with a simple algorithm to compute the *largest* eigenvalue of a matrix and its associated eigenvector.
A natural question is now how one could compute the smallest value or some eigenvalue in between.
In this section we discuss an extension to power iteration, which makes this feasible.
We only need a small ingredient from linear algebra that makes this feasible,
the so-called spectral transformations.

We explore based on a few examples. Consider
"""

# ╔═╡ 4fbed8eb-a929-4afb-a971-d37a75377d8f
Ashift = [ 0.4 -0.6 0.2;
	      -0.3 0.7 -0.4;
	      -0.1 -0.4 0.5]

# ╔═╡ 8290d531-bad1-40b8-a1e9-49cb3a9dae4c
md"Its eigenvalues and eigenvectors are:"

# ╔═╡ 2f97e0aa-e9b1-4bdf-bc81-968a3217b2e6
eigen(Ashift)

# ╔═╡ ad005fc4-04fd-400a-9113-757e52488beb
md"Now we can add a multiple of the identity matrix, e.g.:"

# ╔═╡ 2d5eabf9-ba09-43f8-8ee6-9b35bd826640
σ = 2  # Shift

# ╔═╡ 2dbd51f3-6629-456d-9f45-bf6aeea07be4
eigen(Ashift + σ * I)  # Add 2 * identity matrix

# ╔═╡ b0a4cf39-749d-470f-8a25-f59cda1a740c
md"Notice, how the eigenvectors are the same and only the eigenvalues have been shifted by $σ$.

Similarly we notice:
"

# ╔═╡ 850f1909-b5ab-402a-b4f4-d11decc4d457
let
	A⁻¹ = inv(Ashift + σ * I)
	eigen(A⁻¹)
end

# ╔═╡ a2f76acd-0f52-4d39-8c4a-c106c8db25ca
md"Notice how this matrix still has the same eigenpairs (albeit in a different order)  and simply the *inverted* eigenvalues of $\mathbf A + σ \mathbf I$."

# ╔═╡ 0775cca6-48e4-49eb-b63b-8e183ae2e2ed
1 ./ eigvals(Ashift + σ * I)

# ╔═╡ 677bd623-d1cd-497f-9155-ffd8dba28054
md"""
We now aim to formalise our observation more rigourously,
recall that if $λ_i$ is an eigenvalue of $\mathbf A$ with (non-zero) eigenvector $\mathbf x_i$ if and only if
```math
\tag{7}
\mathbf{A} \mathbf x_i = λ_i \mathbf x_i.
```
We usually refer to the pair $(λ_i, \mathbf x_i)$ as an **eigenpair** of $\mathbf A$.

The statement of the **spectral transformations** is:

!!! info "Theorem 1: Spectral transformations"
    Assume $\mathbf A \in \mathbb R^{n \times n}$ is diagonalisable. 
	Let $(λ_i, \mathbf x_i)$ be an eigenpair of $\mathbf A$, then
	  1. If $A$ is invertible, then $(\frac{1}{λ_i}, x_i)$ is an eigenpair of $A^{-1}$.
	  2. For every $σ \in \mathbb{R}$ we have that
         $(λ_i - σ, x_i)$ is an eigenpair of $\mathbf A - σ \mathbf I$.
	  3. If $\mathbf A - σ \mathbf I$ is invertible, then $(\frac{1}{λ_i - σ}, x_i)$ is an eigenpair of $(\mathbf A - σ \mathbf I)^{-1}$.

> **Proof:** 
> The result follows from a few calculations on top of (7). We proceed 
> in order of the statements above.
> 
> 1. Let $(λ_i, \mathbf x_i)$ be an arbitrary eigenpair of $\mathbf A$.
>    Since $\mathbf A$ is invertible by assumption, $λ_i \neq 0$.
>    Therefore for all eigenvalues $λ_i$ of $\mathbf A$ the fraction
>    $\frac{1}{λ_i}$ is meaningful and we can show:
>    ```math
>    \frac{1}{λ_i} \mathbf x_i = \frac{1}{λ_i} \mathbf I \mathbf x_i
>    = \frac{1}{λ_i} \left(\mathbf A^{-1} \mathbf{A} \right) \mathbf x_i
>    = \frac{1}{λ_i} \mathbf A^{-1} \left(\mathbf{A} \mathbf x_i \right)
>    \stackrel{(6)}{=} \mathbf A^{-1} \mathbf x_i,
>    ```
>    which indeed is the statement that $(\frac{1}{λ_i}, x_i)$ is an eigenpair of $A^{-1}$.
> 
> 2. The argument is similar to 1 and we leave it as an exercise.
> 
> 3. This statement follows by combining statements 1. and 2.

!!! tip "Exercise 1"
    Assume $\mathbf A \in \mathbb R^{n \times n}$ is diagonalisable.
    Prove statement 2. of Theorem 1, that is for all $σ \in \mathbb{R}$:
    If $(λ_i, \mathbf x_i)$ is an eigenpair of $\mathbf A$,
    then $\mathbf x_i$ is also an eigenvector
    of $\mathbf A - σ \mathbf I$ with eigenvalue $λ_i - σ$.
"""

# ╔═╡ 2a2192b5-0c49-4bdd-9470-23c81e3912f7
md"""
Consider point 1. of Theorem 1 and assume that $\mathbf A$ has a *smallest* eigenvalue
```math
|λ_1| ≥ |λ_2| ≥ \cdots |λ_{n-1}| > |λ_n| > 0,
```
then $\mathbf A^{-1}$ has the eigenvalues
```math
|λ_1^{-1}| ≤ |λ_2^{-1}| ≤ \cdots |λ_{n-1}^{-1}| < |λ_n^{-1}|,
```
thus a dominating (i.e. largest) eigenvalue $|λ_n^{-1}|$.
Applying the `power_method` function (Algorithm 1) to $\mathbf A^{-1}$
we thus converge to $|λ_n^{-1}|$ from which we can deduce $|λ_n|$,
the eigenvalue of $\mathbf A$ closest to zero.

Now consider point 3. and assume $σ$ has been chosen such that
```math
\tag{8}
|λ_n - σ| ≥ |λ_{n-1} - σ| ≥ \cdots |λ_{2} - σ| > |λ_1 - σ| > 0,
```
i.e. such that $σ$ is closest to the eigenvalue $λ_1$.
It follows
```math
|λ_n-σ|^{-1} ≤ |λ_{n-1}-σ|^{-1} ≤ \cdots |λ_2-σ|^{-1} < |λ_1-σ|^{-1},
```
such that the `power_method` function converges to $|λ_1-σ|^{-1}$.
From this value we can deduce $λ_1$, since we know $σ$.
In other words by applying Algorithm 1 to the
shift-and-invert matrix $(\mathbf A - σ \mathbf I)^{-1}$
enables us to find the eigenvalue of $\mathbf A$ closest to $σ$.

Note, that a naive application of Algorithm 1 would to compute
```math
\mathbf y^{(k)} = (\mathbf A - σ \mathbf I)^{-1} \mathbf x^{(k)},
```
which is numerically unstable due to the computation of the matrix inverse.
To avoid this computation we obtain $\mathbf y^{(k)}$ by solving
a linear system. We arrive at the following algorithm,
where the changes compared to Algorithm 1 are marked in red.

!!! info "Algorithm 2: Inverse iterations"
     Given
     - a diagonalisable matrix $\mathbf A \in \mathbb R^{n\times n}$,
     - a shift $σ\in \mathbb R$, such that $\mathbf A - σ \mathbf I$ is
       invertible
     - an initial guess $\mathbf x^{(1)} \in \mathbb R^n$
     we iterate for $k = 1, 2, \ldots$:
     1.    $\textcolor{brown}{\text{Solve }(\mathbf A - σ \mathbf I) \mathbf y^{(k)} = \mathbf x^{(k)}
        \text{ for }\mathbf y^{(k)}}$.
     2. Find the index $m$ such that $y^{(k)}_m = \|\mathbf y^{(k)}\|$
     3. Set $α^{(k)} = \frac{1}{y^{(k)}_m}$ and $\textcolor{brown}{β^{(k)} = σ + \frac{x^{(k)}_m}{y^{(k)}_m}}$.
     4. Set $\mathbf x^{(k+1)} = α^{(k)} \mathbf y^{(k)}$

Note the additional change in step 3: In Algorithm 1 we obtained the estimate of the dominant eigenvalue as $y^{(k)}_m / x^{(k)}_m$.
Here this estimate approximates $\frac 1 {λ_1 - σ}$, the dominant eigenvalue of $(\mathbf A - σ \mathbf I)^{-1}$.
Therefore an estimate of $λ_1$ itself is obtained by solving
```math
\frac {y^{(k)}_m} {x^{(k)}_m} = \frac 1 {λ_1 - σ} \quad\Longrightarrow\quad
λ_1 = σ + \frac{x^{(k)}_m}{y^{(k)}_m}
```
for $λ_1$, which yields exactly the expression shown in step 3 of Algorithm 2.

Before we take a look at the implementation, one final point. Each iteration of Algorithm 2 requires the solution of a linear system with the same system matrix $\mathbf B = \mathbf A - σ \mathbf I$. Here, we will attempt this using direct methods, that is PLU factorisation. Since $\mathbf B$ does not change between iterations, the factorisation itself only needs to be computed once, which we
do right before entering the loop. Since for dense matrices
computing the factorisation scales $O(n^3)$,
but solving linear systems based on the factorisation only scales $O(n^2)$
(recall chapter 6), this reduces the cost per iteration.
"""

# ╔═╡ 4873fe3b-d411-4e06-92b3-8cfbed493a12
function inverse_iterations(A, σ, x; maxiter=100)
	# A:  Matrix
	# σ:  shift
	# x:  initial guess
	
    n = size(A, 1)
    x = normalize(x, Inf)

	fact = lu(A - σ*I)   # Compute LU factorisation of A-σI
	
	history = Float64[]
    for k in 1:maxiter
        y = fact \ x
        m = argmax(abs.(y))
		α = 1 / y[m]
		β = σ + x[m] / y[m]
		push!(history, β)
        x = α * y
    end

	(; x, λ=last(history), history)
end

# ╔═╡ 59493fbc-1b97-4ade-bbba-30eea2b5c91a
md"""
Inserting $\mathbf A - σ \mathbf I$ in the place of $\mathbf A$ in (6) and recalling the eigenvalue ordering (8), i.e.
$|λ_n - σ| ≥ |λ_{n-1} - σ| ≥ \cdots |λ_{2} - σ| > |λ_1 - σ| > 0$
one can show similarly that inverse iterations converge again linearly with rate
```math
\tag{9}
r_\text{inviter} =
\lim_{k\to\infty} \frac{|β^{(k+1)} - λ_1|}{|β^{(k)} - λ_1|}
= \left|\frac{λ_1 - σ}{λ_2 - σ}\right|
\qquad \text{as $k\to \infty$}.
```
This implies that convergence is fastest if $σ$ is closest to the targeted
eigenvalue $λ_1$ than to all other eigenvalues of $\mathbf A$.
"""

# ╔═╡ 93173364-1a49-4f87-b814-f95cd00d657f
md"""To demonstrate this rate we again consider the case of diagonalising a triangular matrix
"""

# ╔═╡ c734b379-ab45-4270-8bea-461332887640
begin
	λT = [1, 0.75, 0.6, -0.4, 0]   # The reference eigenvalues we will use
	T  = triu(ones(5, 5), 1) + diagm(λT)
end

# ╔═╡ 6ffd6918-3b4d-4e1f-ba7e-5c825f96681d
md"""
This time the slider allows to tune the value of $σ$ (which we store in the variable `σT`):
- `σT = `$(@bind σT Slider(-0.5:0.01:1.1; default=0.4, show_value=true))
"""

# ╔═╡ 43d637d5-0b62-4c60-a363-e1616f517eb7
md"This value of $σ$ results in a rate:"

# ╔═╡ 908ae6f5-23fe-4ca8-bceb-e484589441a2
begin
	# compute |λ_5 - σ| > |λ_4 - σ| > |λ_3 - σ| > |λ_2 - σ| > |λ_1 - σ|
	# We need sort to get the data in that order
	diff_sorted = sort(abs.(λT .- σT); rev=true)

	# take last element (|λ_1 - σ|) and last-but-one (|λ_2 - σ|)
	r_inviter = diff_sorted[end] / diff_sorted[end-1]
end

# ╔═╡ c0bcba34-3435-4d0a-a52d-a2fbbfd7519c
TODO(md"Give an example using a 3x3 tridiagonal matrix using power iterations, inverse power iterationsn ($σ=0$), inverse power iterations with non-zero shift. One case should not be converging due to two dominant eigenvalues. In each case indicate the eigenvalue to which one converges and the rate of convergence.")

# ╔═╡ 03e47af0-840c-425e-aa82-d265e862da7a
md"""
## Optional: Dynamic shifting

In the discussion in the previous section we noted that the convergence of inverse iterations is best if $σ$ is chosen close to the eigenvalue of the matrix $\mathbf A$. Since inverse iterations (Algorithm 2) actually produce an estimate for the eigenvalue in each iteration (the $β^{(k)}$), a natural idea is to *update* the shift. 

In other words instead of using the same $σ$ in each iteration of Algorithm 2,
we dynamically update $σ = β^{(k)}$ once the new eigenvalue estimate $β^{(k)}$ has been computed.
This also implies that for solving the linear system in step 1, i.e.
```math
(\mathbf A - σ \mathbf I) \, \mathbf y^{(k)} = \mathbf x^{(k)}
```
the system matrix is different for each $k$, such that pre-computing the PLU factorisation once before entering the iteration loop is no longer possible.

We arrive at the following implementation for an **inverse iteration algorithm with dynamic shifting**:
"""

# ╔═╡ f0ecec01-dde6-4d51-9039-cdfc648d16e9
function dynamic_shifting(A, σ, x; maxiter=100, tol=1e-8)
	# A:  Matrix
	# σ:  shift
	# x:  initial guess
	
    n = size(A, 1)
    x = normalize(x, Inf)
		
	history = Float64[]
    for k in 1:maxiter
        y = (A - σ*I) \ x   # PLU-factorise (A - σ*I) and solve system
        m = argmax(abs.(y))
		α = 1 / y[m]
		β = σ + x[m] / y[m]
		push!(history, β)
        x = α * y

		if abs(σ - β) < tol  # We are converged, so exit the iterations.
			break
		end
		σ = β
    end

	(; x, λ=last(history), history)
end

# ╔═╡ 1af4bb50-277f-4368-97b8-2cfcdfea69bc
md"""
A consequence of the fact that we need to compute the PLU factorisation in each iteration is that the overall cost of `dynamic_shifting` is larger than the cost of `inverse_iterations`, but the convergence is also improved:
We go from linear to **quadratic convergence**
(see [chapter 8.3](https://tobydriscoll.net/fnc-julia/krylov/inviter.html#dynamic-shifting) of Driscoll, Brown: Fundamentals of Numerical Computation).

This is easily checked numerically. We use the same matrix
"""

# ╔═╡ b4599ca9-3515-4d91-8e8f-0da11870937e
T

# ╔═╡ f5bf3801-9daa-4abc-a552-561fcce3e48e
md"""
and introduce another slider to tune the value of the *initial* shift $σ$,
which we store in the variable `σD`:
- `σD = `$(@bind σD Slider(-0.5:0.01:1.1; default=0.4, show_value=true))
"""

# ╔═╡ 63b946a7-fa76-4367-98bf-0c52230ef158
md"""The much faster quadratic convergence is clearly visible."""

# ╔═╡ 234362c4-7607-45ec-b4b9-a7c5a08bf350
xstart = randn(size(T, 2));

# ╔═╡ 57709e05-2ca5-43d4-8d7f-19c04393dcb3
let
	# The eigenvalue the iteration targets
	i = argmin(abs.(λT .- σT))
	λtarget = λT[i]

	p = plot(title="Convergence",
		ylabel=L"|β^{(k)} - λ_{\textrm{ref}}|", xlabel=L"k", legend=:bottomleft,
		yaxis=:log, ylims=(1e-10, 3))

	q = plot(title="History", xlabel=L"k", ylabel=L"β^{(k)}",
	         legend=:bottomleft, ylims=(-0.6, 1.2))
	hline!(q, λT, ls=:dash, label=L"eigenvalues of $T$", c=2, lw=1.5)
	hline!(q, [σT], ls=:dash, label=L"shift $σ$", c=3, lw=1.5)

	if abs(σT - λtarget) > 1e-10  # Guard to prevent numerical issues
		results = inverse_iterations(T, σT, xstart; maxiter=20)
		error = abs.(results.history .- λtarget)

		plot!(p, error; mark=:o, label="", lw=1.5)
		plot!(p, k -> r_inviter^k; ls=:dash, label="expected rate", lw=2)

		plot!(q, results.history, mark=:o, c=1, label="", lw=1.5)
	end
	plot(p, q; layout=(1,2))
end

# ╔═╡ 9a72aec8-81fc-4b76-a163-d1fce03bc3d8
let
	# The eigenvalue the iteration targets
	i = argmin(abs.(λT .- σD))
	λtarget = λT[i]

	p = plot(title="Convergence",
		ylabel=L"|β^{(k)} - λ_{\textrm{ref}}|", xlabel=L"k", legend=:bottomleft,
		yaxis=:log, ylims=(1e-10, 3))

	q = plot(title="History", xlabel=L"k", ylabel=L"β^{(k)}",
	         legend=:bottomleft, ylims=(-0.6, 1.2))
	hline!(q, λT, ls=:dash, label=L"eigenvalues of $T$", c=2, lw=1.5)
	hline!(q, [σD], ls=:dash, label=L"initial shift $σ$", c=3, lw=1.5)

	if abs(σD - λtarget) > 1e-10  # Guard to prevent numerical issues
		results = dynamic_shifting(T, σD, xstart; maxiter=7)
		error = abs.(results.history .- λtarget)

		plot!(p, error; mark=:o, label="", lw=1.5)
		plot!(q, results.history, mark=:o, c=1, label="", lw=1.5)
	end
	plot(p, q; layout=(1, 2))
end

# ╔═╡ 1601c0a2-877c-43ed-b692-0754363c7ec4
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
LaTeXStrings = "~1.3.1"
Plots = "~1.40.1"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.56"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "1318db7425f7ff74ae1d5667e50d47e4061fdaf7"

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
git-tree-sha1 = "75bd5b6fc5089df449b5d35fa501c846c9b6549b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.12.0"
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
git-tree-sha1 = "211cdf570992b0d977fda3745f72772e0d5423f2"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.56"

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
# ╠═4949225a-ccf2-11ee-299b-9b834eb6bd42
# ╟─13298dc4-9800-476d-9474-182359a7671b
# ╟─377d9d74-2826-4c53-bc9a-8bfe5a087a84
# ╟─71d74b6a-96b2-4855-b20f-59b00c8f560b
# ╠═6db511ac-a217-49e2-a508-2590b64ea636
# ╟─814e602d-1afa-44ba-bb25-ed32dd66e2b9
# ╠═0dabc016-9cb8-4d1b-bd44-fb9493b8cf28
# ╠═aa47c74d-2202-413a-9cad-17cad20832e1
# ╟─58a88c32-cdf6-4cb9-af08-0ef9a186fe1c
# ╠═5d990087-901a-4457-a9ec-ba3b2d60f470
# ╟─7da94235-399e-4c13-8a02-25d8fa795e17
# ╠═272e7bce-e05e-412f-b2f5-2742d33d0b6f
# ╟─5a4038ac-1905-42e5-88ee-ced49c5375ad
# ╠═076bbe2b-130a-44e5-9218-6fdcc1998cf5
# ╟─db0c875b-9955-4bc1-bdb1-ec2e1a44942d
# ╟─4d1eb1b1-7634-4f38-a564-98f64923e377
# ╟─73c62d23-ac2a-48be-b3c7-0d51ffce773c
# ╟─16d6559c-977b-4fd9-aecb-2b6f61290f04
# ╟─139ad60c-e04f-4dd6-b18c-a6de140a0154
# ╟─1f512d18-1f9c-4730-9676-ca333a7531de
# ╟─2073afcb-7e3c-4d3f-924e-cc768681a8fa
# ╠═01186225-602f-4637-99f2-0a6dd569a703
# ╟─25f2b565-cc6a-4919-ac82-37d6dd62ba16
# ╠═22b9c8b5-02cd-481e-aca7-bb2fe4e85c3c
# ╟─81a6721f-f7e5-42c9-b389-518a4c458ae9
# ╠═033b6a9b-2187-4053-bc69-bbc1f03494fa
# ╟─cc3495af-e53c-4ee8-8ec0-f4b597e17f53
# ╟─c1b01e10-235c-4424-8524-0836a07106ff
# ╟─38015a4f-17af-47c5-a5b6-5e19c429c841
# ╟─d6a56c08-374a-4665-9968-d538fe2c4c59
# ╟─79afdc06-a39e-4187-b91d-0af65a62613b
# ╠═59b21a54-5cac-4e7f-b800-89c1f5f2d01f
# ╟─a8dfd1bc-5146-4218-9fad-0676c7406f16
# ╠═e5bf05e5-b064-4e43-b64b-22e5e45ddcd0
# ╟─385b07d6-15ef-40ee-8d2e-b3b6505a4171
# ╠═5d4697b4-2f3e-452b-afb8-78443efa1f25
# ╟─0e4a443b-e99d-436f-8151-d6cfed9b773b
# ╠═4fbed8eb-a929-4afb-a971-d37a75377d8f
# ╟─8290d531-bad1-40b8-a1e9-49cb3a9dae4c
# ╠═2f97e0aa-e9b1-4bdf-bc81-968a3217b2e6
# ╟─ad005fc4-04fd-400a-9113-757e52488beb
# ╠═2d5eabf9-ba09-43f8-8ee6-9b35bd826640
# ╠═2dbd51f3-6629-456d-9f45-bf6aeea07be4
# ╟─b0a4cf39-749d-470f-8a25-f59cda1a740c
# ╠═850f1909-b5ab-402a-b4f4-d11decc4d457
# ╟─a2f76acd-0f52-4d39-8c4a-c106c8db25ca
# ╠═0775cca6-48e4-49eb-b63b-8e183ae2e2ed
# ╟─677bd623-d1cd-497f-9155-ffd8dba28054
# ╟─2a2192b5-0c49-4bdd-9470-23c81e3912f7
# ╠═4873fe3b-d411-4e06-92b3-8cfbed493a12
# ╟─59493fbc-1b97-4ade-bbba-30eea2b5c91a
# ╟─93173364-1a49-4f87-b814-f95cd00d657f
# ╠═c734b379-ab45-4270-8bea-461332887640
# ╟─6ffd6918-3b4d-4e1f-ba7e-5c825f96681d
# ╟─43d637d5-0b62-4c60-a363-e1616f517eb7
# ╠═908ae6f5-23fe-4ca8-bceb-e484589441a2
# ╟─57709e05-2ca5-43d4-8d7f-19c04393dcb3
# ╠═c0bcba34-3435-4d0a-a52d-a2fbbfd7519c
# ╟─03e47af0-840c-425e-aa82-d265e862da7a
# ╠═f0ecec01-dde6-4d51-9039-cdfc648d16e9
# ╟─1af4bb50-277f-4368-97b8-2cfcdfea69bc
# ╠═b4599ca9-3515-4d91-8e8f-0da11870937e
# ╟─f5bf3801-9daa-4abc-a552-561fcce3e48e
# ╟─9a72aec8-81fc-4b76-a163-d1fce03bc3d8
# ╟─63b946a7-fa76-4367-98bf-0c52230ef158
# ╠═234362c4-7607-45ec-b4b9-a7c5a08bf350
# ╟─1601c0a2-877c-43ed-b692-0754363c7ec4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
