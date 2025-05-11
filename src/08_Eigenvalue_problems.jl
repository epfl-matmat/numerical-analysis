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

# ╔═╡ 4949225a-ccf2-11ee-299b-9b834eb6bd42
begin
	using LinearAlgebra
	using PlutoUI
	using PlutoTeachingTools
	using Plots
	using LaTeXStrings
	using HypertextLiteral
	using Printf
end

# ╔═╡ 34beda8f-7e5f-42eb-b32c-73cfc724062e
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/08_Eigenvalue_problems.pdf)
"""

# ╔═╡ 13298dc4-9800-476d-9474-182359a7671b
TableOfContents()

# ╔═╡ a138fb39-aae0-41b4-bd7b-d2f7eaad7a53
md"""
# Eigenvalue problems

Recall that the eigenpairs of a matrix $\mathbf A$ are the pairs $(λ_i, \mathbf{v}_i)$
of eigenvalues $λ_i$ and eigenvectors $\mathbf{v}_i$ such that
```math
\mathbf{A} \mathbf{v}_i = λ_i \mathbf{v}_i.
```
Geometrically speaking the eigenvectors provide special directions in space along which forming the matrix-vector-product $\mathbf{A}\mathbf{v}_i$ is particularly simple, namely it just *scales* the vector $\mathbf{v}_i$ by a number.

But more generally if $\mathbf x$ is an arbitrary vector
and if for simplicity we assume $\mathbf{A}$ to be symmetric
and positive definite,
then we find
```math
\| \mathbf A \mathbf x \| \leq \| \mathbf A \| \, \| \mathbf x \| \leq 
\sqrt{λ_\text{max}(\mathbf A^T \mathbf A)} \, \| \mathbf x \|
= λ_\text{max}(\mathbf A) \, \| \mathbf x \|
```
where we used the inequalities introduced at the end of  [Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.html).
We note that the **largest eigenvalue of $\mathbf A$**
provides a **bound to the action of $\mathbf{A}$**.
"""

# ╔═╡ a702d4b6-e70c-417d-ac98-92c534d52770
md"""
This may sound technical, but as **matrices are common
in physics and engineering**
and since their **eigenpairs characterise the action of these matrices**,
the computation of
eigenpairs often carries a **physical interpretation**.

For example, in the classical mechanics of rotating objects, the eigenvectors of the **[Moment of inertia](https://en.wikipedia.org/wiki/Moment_of_inertia)** tensor are the **principle axes** along which an object spins without coupling to other rotational degrees of freedom.

In engineering the **eigenvalues of the Hessian matrix (the Laplacian)**
of an object's total energy describe the **resonance frequences**.
That is the frequencies at which the object best absorbs energy from an outside excitation. When these frequencies coincide
with the motion of humans or cars this can lead to a [Resonance disaster](https://en.wikipedia.org/wiki/Mechanical_resonance#Resonance_disaster), e.g. the famous [Tacoma Narrows Bridge Collapse](https://en.wikipedia.org/wiki/Tacoma_Narrows_Bridge_(1940)#Collapse) or the [Millenium Bridge Lateral Excitation Problem](https://en.wikipedia.org/wiki/Millennium_Bridge,_London#Resonance)).
Analysing the eigenfrequencies of bridges is nowadays required as part of the procedure to obtain the permissons for construction.
"""

# ╔═╡ 73889935-2e4c-4a15-a281-b155cf1ca1c9
md"""
In this notebook we will discuss some simple iterative methods for actually computing  eigenpairs. However, the topic is vast and we will only scratch the surface. 
Readers interested in a more in-depth treatment
of eigenvalue problems are encouraged to attend the master class
[MATH-500: Error control in scientific modelling](https://teaching.matmat.org/error-control/).
Some recommended further reading can also be found in the book [Numerical Methods for Large Eigenvalue Problems](https://epubs.siam.org/doi/book/10.1137/1.9781611970739) by Youssef Saad as well as the [Lecture notes on Large Scale Eigenvalue Problems](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf) by Peter Arbenz.
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
@bind dummy Button("Regenerate random vector and rerun experiment")

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
	history = [x]
	for j in 1:6
		x = A * x
		push!(history, x)
		@printf "Iteration %i:  x = [%12.4f, %12.4f]\n" j x[1] x[2]
	end
	maxerror = maximum(abs, (x - A * x))
	@show maxerror
end;

# ╔═╡ 2dbe6cbf-6b27-49c1-b28c-e03a9b11e09d
let
	l = norm(history[1])	
	p = plot(xlims=(-l, l), ylims=(-l, l), legend=:topleft)
	for (i, x) in enumerate(history)
		plot!([0, x[1]], [0, x[2]], arrow=true, label="Iteration $(i-1)", lw=2,
			  colour=cgrad(:speed)[i/length(history)])
	end

	ev = [0.5547, -0.83205]
	plot!([-3ev[1], 3ev[1]], [3ev[2], -3ev[2]], ls=:dash, lw=1.5, c=:grey, label="final direction")
	p
end

# ╔═╡ db0c875b-9955-4bc1-bdb1-ec2e1a44942d
md"""Note how the iterations stabilise, i.e. that $\textbf{x}$ and $\textbf{A} \textbf{x}$ start to be alike. In other words we seem to achieve $\mathbf{A} \mathbf{x} = \mathbf{x}$, which is nothing else than saying that $\mathbf{x}$ is an eigenvector of $\mathbf{A}$ with eigenvalue $1$.
"""

# ╔═╡ 73c62d23-ac2a-48be-b3c7-0d51ffce773c
md"""
Let us understand what happened in this example in detail. We consider the case $\mathbf{A} \in \mathbb{R}^{n \times n}$ diagonalisable and let further
```math
\tag{1}
|λ_1| ≤ |λ_2| ≤ |λ_3| ≤ \cdots ≤ |λ_{n-1}| \textcolor{red}{<} |λ_n|
```
be its eigenvalues with corresponding eigenvectors $\mathbf v_1$, $\mathbf v_2$, $\ldots, \mathbf v_n$,
which we collect column-wise in a unitary matrix $\mathbf{V}$, i.e.
```math
\mathbf{V} = \begin{pmatrix} \mathbf v_1 & \mathbf v_2 & \cdots & \mathbf v_n \\
				\downarrow & \downarrow & \cdots & \downarrow
             \end{pmatrix}
```
Note that $|λ_{n-1}| < |λ_n|$, such that $λ_n$ is non-degenerate and the absolutely largest eigenvalue. We call it the **dominant eigenvalue** of $\mathbf{A}$.
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
Consider applying $\mathbf{A}$ a $k$ number of times to $\mathbf{x}^{(1)}$.
This yields
```math
\tag{2}
\begin{aligned}
\mathbf{A}^k \mathbf{x}^{(1)} &= \mathbf{V} \, \mathbf{D}^k \,\underbrace{\mathbf{V}^{-1} \mathbf{x}^{(1)}}_{= \mathbf{z}} = \mathbf{V} \begin{pmatrix}
λ_1^k z_1 \\
λ_2^k z_2 \\
\vdots\\
λ_n^k z_n
\end{pmatrix} \\ 
&= λ_n^k \left( \left( \frac{λ_1}{λ_n} \right)^k z_1 \mathbf v_1 + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^k z_{n-1} \, \mathbf v_{n-1} + z_n \mathbf v_{n}   \right).
\end{aligned}
```
If $z_n \neq 0$ then
```math
\tag{3}
\left\| \frac{\mathbf{A}^k \mathbf{x}^{(1)}}{λ_n^k} - z_n \textbf v_n \right\| ≤ 
\left|\frac{λ_1}{λ_n}\right|^k \, |z_1| \, \|\textbf v_1\| +
\cdots
+ \left|\frac{λ_{n-1}}{λ_n}\right|^k \, |z_{n-1}| \, \|\textbf v_{n-1}\|
```
Now, each of the terms $\left|\frac{λ_j}{λ_n}\right|^k$
for $j = 1, \ldots, {n-1}$ goes to zero for $k \to \infty$ due to the
ascending eigenvalue ordering of Equation (1).
Therefore overall
```math
\left\| \frac{\mathbf{A}^k \mathbf{x}^{(1)}}{λ_n^k} - z_n \textbf v_n \right\|
\to 0 \qquad \text{as $k \to \infty$}.
```
In other words $\frac{1}{λ_n^k} \mathbf{A}^k \mathbf{x}^{(1)}$ **eventually becomes a multiple of the eigenvector** associated with the dominant eigenvalue.
"""

# ╔═╡ 7fc2ff92-bd21-498e-8e1e-4e8a18355501
md"""
!!! danger "Conventions of eigenpair ordering"
	In the literature as well as across linear algebra codes
	there are multiple conventions regarding the ordering
	of eigenvalues. E.g. whether eigenvalues are ordered from
	smallest to largest, from largest to smallest, whether
	one considers the absolute value of the eigenvalues or
	the signed values etc.
	Wherever this is possible and makes sense we will employ the
	ordering
	```math
	|λ_1| ≤ |λ_2| ≤ |λ_3| ≤ \cdots ≤ |λ_n|
	```
	i.e. that eigenvalues increase in magnitude,
	but their signs may differ.
"""

# ╔═╡ cd72b974-37ab-40ad-b9a6-e930efa3e767
md"""
### Power iteration algorithm
Let's try applying our idea from the earlier section to the matrix
"""

# ╔═╡ a3377f9f-7fba-4e6b-b014-e85de15ca3f0
md"""
- `λₙ = ` $(@bind λₙ Slider([0.2, 1.0, 5.0, 10.0]; default=5.0, show_value=true))
"""

# ╔═╡ 612948f2-abd5-4350-9327-bbfa51cebc57
B = [0.1 5.0;
	 0.0 λₙ]

# ╔═╡ bfefdb2f-f42a-4867-9a5e-94d22226645b
md"""which has eigenvalues $0.1$ and `λₙ = ` $λₙ. The latter is the dominant one, which can further be changed by this slider:""" 

# ╔═╡ 1e461d1e-bcc7-4c9d-8398-f87a51311fdd
md"We again run 6 subsequent applications of $\mathbf B$ and hope the iterations to stabilise:"

# ╔═╡ d4b22b38-04f7-44a7-834b-1bdbad64182d
let
	x = randn(2)
	for j in 1:6
		x = B * x
		@printf "Iteration %i:  x = [%12.4f, %12.4f]\n" j x[1] x[2]
	end
end;

# ╔═╡ 56352bf0-5837-4c38-89c1-62a940e80b3c
md"""
But this time this does not work ... unless `λₙ` happens to be 1.0. 

- This can be understood looking at equation (2): if $λ_n > 1.0$,
  then $λ_n^k$ becomes extremely large, such that $\mathbf{A}^k \mathbf{x}^{(1)}
  ≈ λ_n^k z_n \mathbf v^k$, which grows significantly from one iteration to the next,
  such that no stabilisation is achieved.
  Similarly if $λ_n < 1.0$ then as $k\to \infty$ the $λ_n^k$ becomes smaller
  and smaller.
- Apart from not converging, this also poses difficulties from a numerical point of view, since accurately representing the vector entries of $\mathbf A^k \mathbf{x}$ and performing the associated matrix-vector products becomes more and more difficult the larger the range of numbers that have to be represented.


Naively one could take a look at (3) and just **normalise in each iteration** by dividing by $λ_n$. And indeed this works:
"""

# ╔═╡ 60bbeaee-dc3a-4c44-b016-02d8b1b5623b
let
	x = randn(2)
	for j in 1:6
		x = x / λₙ   # Normalise
		x = B * x
		@printf "Iteration %i:  x = [%12.4f, %12.4f]\n" j x[1] x[2]
	end
end;

# ╔═╡ b728388d-e6e4-475a-adce-149b2b53a1d4
md"""
The problem here is that for general problems **we don't know $λ_n$**,
so we cannot use it for normalisation.
Fortunately it turns out, however, that pretty much any normalisation of $\mathbf x$ works. For example we can use the infinity norm:
```math
\|x\|_\infty = \max_{i=1,\ldots n} |x_i|,
```
which simply selects the largest absolute entry of the vector.
Again this works as expected:
"""

# ╔═╡ 2425ae53-d861-47f7-9f86-750cf4f363c1
let
	x = randn(2)
	for j in 1:6
		x = x / maximum(abs.(x))  # Normalise 
		x = B * x
		@printf "Iteration %i:  x = [%12.4f, %12.4f]\n" j x[1] x[2]
	end
	@printf "Estimate for eigenvalue: %12.6f" maximum(abs.(x))
end;

# ╔═╡ a11825e7-9ee5-416c-b170-edd1c5eb746c
md"We also note in passing that $\|x\|_\infty$ seems to converge to the dominant eigenvalue."

# ╔═╡ b69a8d6c-364a-4951-afd9-24588ac10b64
md"""
Based on this idea we formulate the algorithm

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

"""

# ╔═╡ 4ebfc860-e179-4c76-8fc5-8c1089301078
md"""
Note that in this algorithm
$y_m^{(k)} = \| y^{(k)} \|_\infty$
and $α^{(k)} = 1 / \| y^{(k)} \|_\infty$.
Step 4 is thus performing the normalisation we developed above.
Furthermore instead of employing $\| x \|_\infty$
as the eigenvalue estimate it employs
$β^{(k)}$, which is just a scaled version of $\| x \|_\infty$.
To see why this should be expected to be an estimate of the
dominant eigenvalue $λ_n$ assume that
$\textbf{x}^{(k)}$ is already close to the associated eigenvector $\mathbf v_n$.
Then $\mathbf A\, \mathbf{x}^{(k)}$ is almost $\lambda_n \mathbf{x}^{(k)}$,
such that the ratio $β^{(k)} = \frac{y^{(k)}_m}{x^{(k)}_m}$
becomes close to $λ_n$ itself.

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

# ╔═╡ 25836ba5-5850-44d5-981b-3cc4dd5c7fdc
md"""
We want to demonstrate why $β^{(k)} \to λ_n$ as $k\to\infty$ more rigorously.

Let us first note that for our algorithm
```math
\tag{4}
\textbf{x}^{(k)} = α^{(1)} \cdot α^{(2)} \cdot \, \cdots \, \cdot \, α^{(k)} \cdot \, \textbf{A}^k \textbf x^{(1)}
= \prod_{i=1}^k α^{(i)} \cdot \mathbf A^k \mathbf x^{(1)}
```
and by construction we always have $\left\| \textbf{x}^{(k)} \right\|_\infty = 1$. 
We further note
```math
\begin{aligned}
β^{(k)} &= \frac{y^{(k)}_m}{x^{(k)}_m}
= \frac{\left(\mathbf{A}\, \mathbf{x}^{(k)}\right)_\textcolor{red}{m}}{x^{(k)}_m}\\
&\stackrel{(4)}{=} \frac{
\left( \prod_{i=1}^k α^{(i)} \cdot \mathbf A^{k+1} \mathbf x^{(1)} \right)_\textcolor{red}{m}
}{
\left( \prod_{i=1}^k α^{(i)} \cdot \mathbf A^k \mathbf x^{(1)} \right)_\textcolor{red}{m}
}\\
&\stackrel{(2)}{=}
\frac{
λ_n^{k+1} \left( 
\left( \frac{λ_1}{λ_n} \right)^{k+1} z_1 \mathbf v_1 + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k+1} z_{n-1} \, \mathbf v_{n-1} + z_n \mathbf v_{n} 
\right)_\textcolor{red}{m}
}{
λ_n^{k} \left( 
\left(\frac{λ_1}{λ_n} \right)^{k} z_1 \mathbf v_1 + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k} z_{n-1} \, \mathbf v_{n-1} + z_n \mathbf v_{n}  
\right)_\textcolor{red}{m}
}\\
&= 
\frac{
λ_n \left( 
\left( \frac{λ_1}{λ_n} \right)^{k+1} b_1 + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k+1} b_{n-1} + b_{n} 
\right)
}{
\left( 
\left(\frac{λ_1}{λ_n} \right)^{k} b_1 + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k} b_{n-1} + b_{n}  
\right)
}
\end{aligned}
```
where $_\textcolor{red}{m}$ denotes that we take the $m$-th element of the vector in the brackets
and we further defined $b_i = z_i \left(\mathbf{v}_i\right)_m$, that is $z_i$ multiplied by the $m$-th element of the eigenvector $\mathbf{v}_i$.
Again we assume $z_n \neq 0$, which implies $b_n \neq 0$.
Under this assumption
```math
\tag{5}
β^{(k)} = \frac{y^{(k)}_m}{x^{(k)}_m}=
λ_n \, \frac{
\left( \frac{λ_1}{λ_n} \right)^{k+1} \frac{b_1}{b_n} + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k+1} \frac{b_{n-1}}{b_n} + 1
}{
\left( \frac{λ_1}{λ_n} \right)^{k} \frac{b_1}{b_n} + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k} \frac{b_{n-1}}{b_n} + 1
}
\to λ_n \qquad \text{as $k\to\infty$}.
```
Keeping in mind that $\frac{λ_j}{λ_n} < 1$ for all $1 ≤ i ≤ n-1$ we indeed
observe $β^{(k)}$ to converge to $λ_n$.
"""

# ╔═╡ cc3495af-e53c-4ee8-8ec0-f4b597e17f53
Foldable("A remark on zₙ ≠ 0", 
md"""
In the theoretical discussions so far we always assumed $z_n ≠ 0$,
i.e. that the initial guess $\mathbf{x}^{(1)}$ has a non-zero component along
the eigenvector $\mathbf v_n$ of the dominant eigenvalue. This may not
always be the case. However even if $z_n ≠ 0$,
computing $\mathbf{A} \mathbf{x}^{(1)}$
in finite-precision floating-point arithmetic is associated with
a small error, which means that 
$\mathbf{A} \mathbf{x}^{(1)}$ and thus
$\mathbf{x}^{(2)}$ does usually
*not* have a zero component along $\mathbf v_n$:
Assuming $z_n ≠ 0$ is not a restriction for practical calculations.
""")

# ╔═╡ c1b01e10-235c-4424-8524-0836a07106ff
md"""
### Convergence of power method

The above plot already suggests a **linear convergence** towards
the exact eigenvalue. Indeed a detailed analysis shows that the convergence rate
can be computed as
```math
\tag{6}
r_{\textrm{power}} = \lim_{k\to\infty} \frac{|β^{(k+1)} - λ_n|}{|β^{(k)} - λ_n|}
= \left|\frac{λ_{n-1}}{λ_n}\right|.
```
"""

# ╔═╡ 38015a4f-17af-47c5-a5b6-5e19c429c841
details("Optional: Detailed derivation",
md"""
To show the rate $r_{\textrm{power}}$ we revisit equation (5)
```math
β^{(k)} =
λ_n \, \frac{
\left( \frac{λ_1}{λ_n} \right)^{k+1} \frac{b_1}{b_n} + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k+1} \frac{b_{n-1}}{b_n} + 1
}{
\left( \frac{λ_1}{λ_n} \right)^{k} \frac{b_1}{b_n} + \cdots + \left( \frac{λ_{n-1}}{λ_n} \right)^{k} \frac{b_{n-1}}{b_n} + 1
}
```
and introduce the shorthand notation
```math
B_i = b_i/b_n \qquad\text{and}\qquad
r_i = λ_i/λ_n
```
such that
```math
β^{(k)} =
λ_n \frac{
r_1^{k+1} B_1 + r_2^{k+1} B_2 + \cdots + r_{n-1}^{k+1} B_{n-1} + 1
}{
r_1^{k}   B_1 + r_2^{k}   B_2 + \cdots + r_{n-1}^{k}   B_{n-1} + 1
}
```
We will make the *additional assumption* $|λ_{n-2}| < |λ_{n-1}|$,
such that $|λ_{n-2} / λ_{n-1}| < 1$.
We could also complete this development without this assumption,
but the resulting expressions are involved.
Now
```math
\begin{aligned}
r_1^{k}   B_1 &+ r_2^{k}   B_2 + \cdots + r_{n-1}^k   B_{n-1} + 1\\
&= 1 + r_{n-1}^k \left[
  \left(\frac{λ_1}{λ_{n-1}}\right)^k B_1
+ \cdots 
+ \left(\frac{λ_{n-2}}{λ_{n-1}}\right)^k B_{n-2}
+ B_{n-1}
\right]
\end{aligned}
```
clearly approaches $1 + r_{n-1}^k B_{n-1}$ as $k \to \infty$.
Noting for $k \to \infty$ that
```math
\frac{1}{1 + r_{n-1}^k B_{n-1}} = 1 - r_{n-1}^k B_{n-1} + O(r_{n-1}^{2k})
```
we can conclude in this limit
```math
\begin{aligned}
β^{(k)} - λ_n &= λ_n \left(1 +  r_{n-1}^{k+1} B_{n-1}\right)
\left[ 1 - r_{n-1}^k B_{n-1} + O(r_{n-1}^{2k}) \right] - λ_n\\
&= (r_{n-1} - 1) r_{n-1}^k B_{n-1} + O(r_{n-1}^{2k}).
\end{aligned}
```
This is indeed linear convergence with convergence rate
```math
\frac{|β^{(k\textcolor{red}{+1})} - λ_n|}{|β^{(k)} - λ_n|}
\to \frac{|(r_{n-1} - 1) r_{n-1}^{k\textcolor{red}{+1}} B_{n-1}|}{|(r_{n-1} - 1) r_{n-1}^{k} B_{n-1}|} = |r_{n-1}| = \left|\frac{λ_{n-1}}{λ_n}\right|
\quad \text{as $k \to \infty$}
```
""")

# ╔═╡ d6a56c08-374a-4665-9968-d538fe2c4c59
md"""
So the smaller the ratio between $|λ_{n-1}|$ and $|λ_n|$,
the faster the convergence.
"""

# ╔═╡ a8dfd1bc-5146-4218-9fad-0676c7406f16
md"and we will put these on the diagonal of a triangular matrix, such that the eigenvalues of the resulting matrix $M$ are exactly the $λ_δ$:"

# ╔═╡ 385b07d6-15ef-40ee-8d2e-b3b6505a4171
md"""
We introduce a slider to tune the $λ_{n-1}$ of equation (6):
- `λₙ₋₁ = `$(@bind λₙ₋₁ Slider(-0.6:-0.05:-0.95; show_value=true, default=-0.9) )
"""

# ╔═╡ 79afdc06-a39e-4187-b91d-0af65a62613b
md"Let us take a case where $λ_n = 1.0$. Therefore the size of $λ_{n-1}$ determines the rate of convergence. Let's take $λ_{n-1}$ = $(round(λₙ₋₁; sigdigits=2))."

# ╔═╡ 59b21a54-5cac-4e7f-b800-89c1f5f2d01f
begin
	λ_δ = [0.0, -0.4, 0.6, λₙ₋₁, 1.0]   # The reference eigenvalues we will use
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

	λₙ = 1.0
	r_power = abs(λₙ₋₁ / λₙ)	
	plot!(p, k -> r_power^k; ls=:dash, label=L"expected rate $| λ_{n-1} / λ_n|$", lw=2)
end

# ╔═╡ 0e4a443b-e99d-436f-8151-d6cfed9b773b
md"""
## Spectral transformations

The power method provides us with a simple algorithm to compute the *largest* eigenvalue of a matrix and its associated eigenvector.
But what if one actually wanted to compute the smallest eigenvalue
or an eigenvalue somewhere in the middle ?

In this section we discuss an extension to power iteration,
which makes this feasible.
We only need a small ingredient from linear algebra
the **spectral transformations**.

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
md"Now we add a multiple of the identity matrix, e.g.:"

# ╔═╡ 2d5eabf9-ba09-43f8-8ee6-9b35bd826640
σ = 2  # Shift

# ╔═╡ 2dbd51f3-6629-456d-9f45-bf6aeea07be4
eigen(Ashift + σ * I)  # Add 2 * identity matrix

# ╔═╡ b0a4cf39-749d-470f-8a25-f59cda1a740c
md"Notice, how the eigenvectors are the same and only the eigenvalues have been shifted by $σ$.

Similarly:
"

# ╔═╡ 850f1909-b5ab-402a-b4f4-d11decc4d457
let
	A⁻¹ = inv(Ashift + σ * I)
	eigen(A⁻¹)
end

# ╔═╡ a2f76acd-0f52-4d39-8c4a-c106c8db25ca
md"Notice how this matrix still has the same eigenvectors (albeit in a different order) and the *inverted* eigenvalues of $\mathbf A + σ \mathbf I$."

# ╔═╡ 0775cca6-48e4-49eb-b63b-8e183ae2e2ed
1 ./ eigvals(Ashift + σ * I)

# ╔═╡ b25fb4e5-127b-4939-a44b-001cc68ea811
md"""
We want to formalise our observation more rigourously.
Recall that $λ_i$ is an eigenvalue of $\mathbf A$ with (non-zero) eigenvector $\mathbf x_i$ if and only if
```math
\tag{7}
\mathbf{A} \mathbf x_i = λ_i \mathbf x_i.
```
We usually refer to the pair $(λ_i, \mathbf x_i)$ as an **eigenpair** of $\mathbf A$.

The statement of the **spectral transformations** is:
"""

# ╔═╡ c19c7118-37e5-4bfd-adc9-434a96358abd
md"""
!!! info "Theorem 1: Spectral transformations"
    Assume $\mathbf A \in \mathbb R^{n \times n}$ is diagonalisable. 
	Let $(λ_i, \mathbf x_i)$ be an eigenpair of $\mathbf A$, then
	  1. If $\mathbf A$ is invertible, then $(\frac{1}{λ_i}, x_i)$ is an eigenpair of $A^{-1}$.
	  2. For every $σ \in \mathbb{R}$ we have that
         $(λ_i - σ, x_i)$ is an eigenpair of $\mathbf A - σ \mathbf I$.
	  3. If $\mathbf A - σ \mathbf I$ is invertible, then $(\frac{1}{λ_i - σ}, x_i)$ is an eigenpair of $(\mathbf A - σ \mathbf I)^{-1}$.
"""

# ╔═╡ f5014b54-79ea-4506-8ddf-8262bd09d9f9
md"""
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
"""

# ╔═╡ 5ea67b6d-a297-403e-bb3e-e69f647f90a8
md"""
!!! tip "Exercise 1"
    Assume $\mathbf A \in \mathbb R^{n \times n}$ is diagonalisable.
    Prove statement 2. of Theorem 1, that is for all $σ \in \mathbb{R}$:
    If $(λ_i, \mathbf x_i)$ is an eigenpair of $\mathbf A$,
    then $\mathbf x_i$ is also an eigenvector
    of $\mathbf A - σ \mathbf I$ with eigenvalue $λ_i - σ$.
"""

# ╔═╡ 9b32cc52-b74d-49fe-ba3f-192f5efc4ef2
md"""
Consider point 1. of Theorem 1 and assume that $\mathbf A$
has a *smallest* eigenvalue
```math
0 < |λ_1| < |λ_2| ≤ |λ_3| ≤ \cdots ≤ |λ_n|
```
then $\mathbf A^{-1}$ has the eigenvalues
```math
|λ_1^{-1}| > |λ_2^{-1}| ≥ |λ_3^{-1}| ≥ \cdots ≥ |λ_n^{-1}|,
```
thus a dominating (i.e. largest) eigenvalue $|λ_1^{-1}|$.
**Applying** the `power_method` function (Algorithm 1) **to $\mathbf A^{-1}$**
we thus converge to $λ_1^{-1}$ from which we can deduce $λ_1$,
the **eigenvalue of $\mathbf A$ closest to zero**.
"""

# ╔═╡ 36dc90e6-fd8a-4c88-a816-140b663532ef
md"""
Now consider point 3. and assume $σ$ has been chosen such that
```math
\tag{8}
|λ_n - σ| ≥ |λ_{n-1} - σ| ≥ \cdots |λ_{2} - σ| > |λ_1 - σ| > 0,
```
i.e. such that $σ$ is closest to the eigenvalue $λ_1$.
It follows
```math
|(λ_n-σ)^{-1}| ≤ |(λ_{n-1}-σ)^{-1}| ≤ \cdots |(λ_2-σ)^{-1}| < |(λ_1-σ)^{-1}|,
```
such that the `power_method` function converges to $(λ_1-σ)^{-1}$.
From this value we can deduce $λ_1$, since we know $σ$.
In other words by **applying Algorithm 1** to the
**shift-and-invert matrix** $(\mathbf A - σ \mathbf I)^{-1}$
enables us to find the **eigenvalue of $\mathbf A$ closest to $σ$**.
"""

# ╔═╡ 29a4d3a6-29e0-4067-b703-43da5171f7bf
md"""
A naive application of Algorithm 1 would first compute $\mathbf{P} = (\mathbf A - σ \mathbf I)^{-1}$ and then apply 
```math
\mathbf y^{(k)} = (\mathbf A - σ \mathbf I)^{-1} \mathbf x^{(k)} \mathbf{P} \mathbf x^{(k)}
```
in each step of the power iteration.
However, for many problems the explicit computation of the inverse
$\mathbf{P}$ is numerically unstable.
Instead of computing $\mathbf{P}$ explicitly
one instead obtains $\mathbf y^{(k)}$ by **solving a linear system**
```math
(\mathbf A - σ \mathbf I) \, \mathbf y^{(k)} = \mathbf x^{(k)}
```
for $\mathbf y^{(k)}$,
which is done **using LU factorisation**.
We arrive at the following algorithm,
where the changes compared to Algorithm 1 are marked in red.
"""

# ╔═╡ 49838e6f-63c2-4494-a5d8-9e3dd25554d3
md"""
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
"""

# ╔═╡ 3e0e12ad-4d09-4128-9175-14e71c7de5a2
md"""
Note the additional change in step 3: In Algorithm 1 we obtained the estimate of the dominant eigenvalue as $y^{(k)}_m / x^{(k)}_m$.
Here this estimate approximates $\frac 1 {λ_1 - σ}$, the dominant eigenvalue of $(\mathbf A - σ \mathbf I)^{-1}$.
Therefore an estimate of $λ_1$ itself is obtained by solving
```math
\frac {y^{(k)}_m} {x^{(k)}_m} = \frac 1 {λ_1 - σ} \quad\Longrightarrow\quad
λ_1 = σ + \frac{x^{(k)}_m}{y^{(k)}_m}
```
for $λ_1$, which yields exactly the expression shown in step 3 of Algorithm 2.

An implementation of this algorithm is given in:
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

# ╔═╡ 67d258d1-9959-46ad-9f13-2bc9abb4da37
md"""
Notice that in this implementation we make use of the fact that
once an **LU factorisation is computed** it can be **re-used
for solving multiple linear systems** with changing right-hand sides:
in our case we can expect to solve many linear systems
involving the matrix $\mathbf A - σ \mathbf I$.
Therefore we **compute the LU factorisation** only once,
namely at the beginning of the algorithm and before entering
the iterative loop.

Since for dense matrices
computing the factorisation scales $O(n^3)$,
but solving linear systems based on the factorisation only scales $O(n^2)$
([recall chapter 6](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.html)), this reduces the cost per iteration.
"""

# ╔═╡ 8e01a98d-c49f-43b3-9681-07d8e4b7f12a
md"""
To investigate the possibilities enabled by inverse iterations,
we consider a few examples using the following triangular matrix
"""

# ╔═╡ 19e64038-afb3-4df6-ad47-16bad5a98a7b
C = [-4.0  3.0  4.0;
	  0.0  4.0  5.0;
	  0.0  0.0  2.0]

# ╔═╡ 75ab9929-200d-4a69-bb84-187d962662eb
md"""which has eigenvalues $λ_1 = -4$, $λ_2 = 4$ and $λ_3 = 2$.

Since it has no unique dominant eigenvalue ($|λ_1| = |λ_2| = 4$)
**plain power iterations** do not converge, but rather oscillate:
"""

# ╔═╡ fe15a206-e650-4d59-962c-a00bc226bf48
let
	xinit = randn(size(C, 2))
	res = power_method(C, xinit; maxiter=30)

	plot(res.history; mark=:o, c=1, label="", lw=1.5, ylabel=L"β^{(k)}", xlabel=L"k",
		 title="Power iterations on C")
end

# ╔═╡ a152d4e2-afde-4b57-8a22-33f5f26f6880
md"""
In contrast for **inverse iterations** (with $σ=0$) what matters are the eigenvalues of $\mathbf{C}^{-1}$, which are $1 / λ_i$ or
```math
\frac{1}{λ_1} = - \frac 14
\qquad
\frac{1}{λ_2} = \frac 14
\qquad
\frac{1}{λ_3} = \frac 12
```
Therefore there is a single dominant eigenvalue ($\frac 12$) and inverse iterations will converge to the eigenvalue $2$:
"""

# ╔═╡ d467f287-0e29-48fe-adb5-71c52aa66cfc
let
	σ = 0.0  # Just perform plain inverse iterations
	
	xinit = randn(size(C, 2))
	res = inverse_iterations(C, σ, xinit; maxiter=30)

	plot(res.history; mark=:o, c=1, label="", lw=1.5, ylabel=L"β^{(k)}", xlabel=L"k",
		 title="Inverse iterations on C with σ = $σ")
end

# ╔═╡ 9345e24b-f19b-46ce-8e35-fb8acbaabb53
md"""
Now suppose we knew that $\mathbf{C}$ has one eigenvalue near $-6$, such that we take $σ = -6$. Then what matters for convergence in the inverse iterations are the eigenvalues of $(\mathbf C - (-6) \mathbf I)^{-1}$, which are
$\mu_i = \frac1{λ_i-σ}$ or
```math
\mu_1 = \frac{1}{-4 + 6} = \frac12 \qquad \mu_2 = \frac{1}{4+6} = \frac1{10} \qquad \mu_3 = \frac1{2+6} = \frac{1}{8},
```
such that $\mu_1 = \frac12$ is the dominating eigenvalue. Therefore
with $σ = -6$ inverse iterations converge to $-6 + 1 / \mu_1 = -6 + 1 / (\frac{1}2) = -4$:
"""

# ╔═╡ 10558bcb-43d4-41cd-b39f-63ec58e86b1b
let
	σ = -6.0
	
	xinit = randn(size(C, 2))
	res = inverse_iterations(C, σ, xinit; maxiter=30)

	plot(res.history; mark=:o, c=1, label="", lw=1.5, ylabel=L"β^{(k)}", xlabel=L"k",
		 title="Inverse iterations on C with σ = $σ")
end

# ╔═╡ 59493fbc-1b97-4ade-bbba-30eea2b5c91a
md"""
### Convergence of inverse iterations

Inserting $\mathbf A - σ \mathbf I$ in the place of $\mathbf A$ in (6) and recalling the eigenvalue ordering (8), i.e.
$|λ_n - σ| ≥ |λ_{n-1} - σ| ≥ \cdots |λ_{2} - σ| > |λ_1 - σ| > 0$
one can show similarly that **inverse iterations converge linearly** with rate
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

# ╔═╡ d998c91f-0d3c-4846-a6df-9d2c407f5836
md"""
Finally let us consider a case where we use a slider to be able to
change the applied shift:
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

# ╔═╡ 03e47af0-840c-425e-aa82-d265e862da7a
md"""
## Optional: Dynamic shifting

In the discussion in the previous section we noted that the convergence of inverse iterations is best if $σ$ is chosen close to the eigenvalue of the matrix $\mathbf A$. Since **inverse iterations** (Algorithm 2) actually **produce an estimate for the eigenvalue** in each iteration (the $β^{(k)}$),
a natural idea is to **update the shift**:
instead of using the same $σ$ in each iteration,
we select a different $σ = β^{(k)}$ *dynamically*
based on the best eigenvalue estimate currently available to us.

This also implies that for solving the linear system in step 1, i.e.
$(\mathbf A - σ \mathbf I) \, \mathbf y^{(k)} = \mathbf x^{(k)}$
the system matrix is different for each $k$.
Therefore pre-computing the LU factorisation
before entering the iteration loop is no longer possible.

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
        y = (A - σ*I) \ x   # LU-factorise (A - σ*I) and solve system
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
Since we now need to compute a fresh LU factorisation
in each iteration, the overall cost of `dynamic_shifting`
is larger than the cost of `inverse_iterations`.
On the upside, convergence is improved:
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
	Sidebar(Markdown.parse(read("sidebar.md", String)), 320)
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
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.1"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.56"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "ed4b990b71c74481d1b3d103b825239737d31586"

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
# ╟─34beda8f-7e5f-42eb-b32c-73cfc724062e
# ╠═4949225a-ccf2-11ee-299b-9b834eb6bd42
# ╟─13298dc4-9800-476d-9474-182359a7671b
# ╟─a138fb39-aae0-41b4-bd7b-d2f7eaad7a53
# ╟─a702d4b6-e70c-417d-ac98-92c534d52770
# ╟─73889935-2e4c-4a15-a281-b155cf1ca1c9
# ╟─71d74b6a-96b2-4855-b20f-59b00c8f560b
# ╠═6db511ac-a217-49e2-a508-2590b64ea636
# ╟─814e602d-1afa-44ba-bb25-ed32dd66e2b9
# ╠═0dabc016-9cb8-4d1b-bd44-fb9493b8cf28
# ╠═aa47c74d-2202-413a-9cad-17cad20832e1
# ╟─58a88c32-cdf6-4cb9-af08-0ef9a186fe1c
# ╠═5d990087-901a-4457-a9ec-ba3b2d60f470
# ╟─7da94235-399e-4c13-8a02-25d8fa795e17
# ╠═272e7bce-e05e-412f-b2f5-2742d33d0b6f
# ╟─2dbe6cbf-6b27-49c1-b28c-e03a9b11e09d
# ╟─5a4038ac-1905-42e5-88ee-ced49c5375ad
# ╟─db0c875b-9955-4bc1-bdb1-ec2e1a44942d
# ╟─73c62d23-ac2a-48be-b3c7-0d51ffce773c
# ╟─16d6559c-977b-4fd9-aecb-2b6f61290f04
# ╟─7fc2ff92-bd21-498e-8e1e-4e8a18355501
# ╟─cd72b974-37ab-40ad-b9a6-e930efa3e767
# ╠═612948f2-abd5-4350-9327-bbfa51cebc57
# ╟─bfefdb2f-f42a-4867-9a5e-94d22226645b
# ╟─a3377f9f-7fba-4e6b-b014-e85de15ca3f0
# ╟─1e461d1e-bcc7-4c9d-8398-f87a51311fdd
# ╠═d4b22b38-04f7-44a7-834b-1bdbad64182d
# ╟─56352bf0-5837-4c38-89c1-62a940e80b3c
# ╠═60bbeaee-dc3a-4c44-b016-02d8b1b5623b
# ╟─b728388d-e6e4-475a-adce-149b2b53a1d4
# ╠═2425ae53-d861-47f7-9f86-750cf4f363c1
# ╟─a11825e7-9ee5-416c-b170-edd1c5eb746c
# ╟─b69a8d6c-364a-4951-afd9-24588ac10b64
# ╟─4ebfc860-e179-4c76-8fc5-8c1089301078
# ╠═01186225-602f-4637-99f2-0a6dd569a703
# ╟─25f2b565-cc6a-4919-ac82-37d6dd62ba16
# ╠═22b9c8b5-02cd-481e-aca7-bb2fe4e85c3c
# ╟─81a6721f-f7e5-42c9-b389-518a4c458ae9
# ╠═033b6a9b-2187-4053-bc69-bbc1f03494fa
# ╟─25836ba5-5850-44d5-981b-3cc4dd5c7fdc
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
# ╟─b25fb4e5-127b-4939-a44b-001cc68ea811
# ╟─c19c7118-37e5-4bfd-adc9-434a96358abd
# ╟─f5014b54-79ea-4506-8ddf-8262bd09d9f9
# ╟─5ea67b6d-a297-403e-bb3e-e69f647f90a8
# ╟─9b32cc52-b74d-49fe-ba3f-192f5efc4ef2
# ╟─36dc90e6-fd8a-4c88-a816-140b663532ef
# ╟─29a4d3a6-29e0-4067-b703-43da5171f7bf
# ╟─49838e6f-63c2-4494-a5d8-9e3dd25554d3
# ╟─3e0e12ad-4d09-4128-9175-14e71c7de5a2
# ╠═4873fe3b-d411-4e06-92b3-8cfbed493a12
# ╟─67d258d1-9959-46ad-9f13-2bc9abb4da37
# ╟─8e01a98d-c49f-43b3-9681-07d8e4b7f12a
# ╠═19e64038-afb3-4df6-ad47-16bad5a98a7b
# ╟─75ab9929-200d-4a69-bb84-187d962662eb
# ╠═fe15a206-e650-4d59-962c-a00bc226bf48
# ╟─a152d4e2-afde-4b57-8a22-33f5f26f6880
# ╠═d467f287-0e29-48fe-adb5-71c52aa66cfc
# ╟─9345e24b-f19b-46ce-8e35-fb8acbaabb53
# ╠═10558bcb-43d4-41cd-b39f-63ec58e86b1b
# ╟─59493fbc-1b97-4ade-bbba-30eea2b5c91a
# ╟─d998c91f-0d3c-4846-a6df-9d2c407f5836
# ╠═c734b379-ab45-4270-8bea-461332887640
# ╟─6ffd6918-3b4d-4e1f-ba7e-5c825f96681d
# ╟─43d637d5-0b62-4c60-a363-e1616f517eb7
# ╠═908ae6f5-23fe-4ca8-bceb-e484589441a2
# ╠═57709e05-2ca5-43d4-8d7f-19c04393dcb3
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
