### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0501f270-c39c-11ee-2077-7f7f409b6c66
begin
	using LinearAlgebra
	using PlutoUI
	using PlutoTeachingTools
	using Plots
	using LaTeXStrings
	using HypertextLiteral
end

# ╔═╡ 63bb7fe9-750f-4d2f-9d18-8374b113373e
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/07_Iterative_methods.pdf)
"""

# ╔═╡ 7d9c9392-3aec-4efd-a9ba-d8965687b163
TableOfContents()

# ╔═╡ 4e0e5311-3428-44f5-a616-68cf7634214f
md"""
# Iterative methods for linear systems

In the previous notebook we looked at direct methods for solving linear systems
$\mathbf{A} \mathbf{x} = \mathbf{b}$ based on LU factorisation of the system matrix $\mathbf{A}$. We saw that even for sparse matrices we may need $O(n^3)$ computational time to perform the factorisation and $O(n^2)$ memory and to store the $\mathbf{L}$ and $\mathbf{U}$ factors due to fill in.

Both can make the solution of linear systems prohibitively expensive for large matrices. In this notebook we will develop *iterative methods*.
These build a sequence of solution vectors $\mathbf{x}^{(k)}$,
which is designed to converge to the solution $\mathbf{x}$,
that is $\lim_{k \to \infty} \mathbf{x}^{(k)} = \mathbf{x}$.
Typically these methods are based on performing
matrix-vector products $\textbf{A} \textbf{v}$
in each iteration step.
Exactly this difference to the direct methods
--- that is *not* employing the matrix $\textbf{A}$, but only the matrix-vector product ---
is what often leads to an overall reduction of the computational cost.

"""

# ╔═╡ b5bf49df-f6ff-49f9-b016-7f0159b4af7d
md"""
## Richardson iterations

A general framework for building iterative methods for solving linear systems
is Richardson's method. The idea is to introduce a *matrix splitting* $\mathbf{A} = \mathbf{P} - \mathbf{R}$ where $\mathbf{P}$ is an invertible matrix
and $\mathbf{R} = \mathbf{P} - \mathbf{A}$.
With this idea we rewrite $\mathbf{A} \mathbf{x} = \mathbf{b}$ to
```math
(\mathbf{P} - \mathbf{R})\mathbf{x} = \mathbf{b}
\qquad \Leftrightarrow \qquad
\mathbf{P} \mathbf{x} = (\mathbf{P} - \mathbf{A}) \mathbf{x} + \mathbf{b}
\qquad \Leftrightarrow \qquad
\mathbf{x} = \mathbf{P}^{-1} (\mathbf{P} - \mathbf{A}) \mathbf{x} + \mathbf{b}
```
If we define $g(\mathbf{x}) = \mathbf{P}^{-1} (\mathbf{P} - \mathbf{A}) \mathbf{x} + \mathbf{b}$ we observe this to be exactly the setting of
the multi-dimensional fixed-point methods we developed
in chapter 4,
i.e. our goal is to obtain a solution $\mathbf{x}_\ast$
with $\mathbf{x}_\ast = g(\mathbf{x}_\ast) = \mathbf{P}^{-1} (\mathbf{P} - \mathbf{A}) \mathbf{x}_\ast + \mathbf{b}$.

Starting from an initial vector $\mathbf{x}^{(0)}$ we thus iterate
```math
\tag{2}
\mathbf{P} \mathbf{x}^{(k+1)} = (\mathbf{P} - \mathbf{A}) \mathbf{x}^{(k)} + \mathbf{b},
\qquad k = 0, 1, \ldots,
```
that is in each step we solve the linear system $\mathbf{P} \mathbf{x}^{(k+1)} = \widetilde{\mathbf{b}}^{(k)}$ where $\widetilde{\mathbf{b}}^{(k)} = (\mathbf{P} - \mathbf{A}) \mathbf{x}^{(k)} + \mathbf{b}$. Assuming this sequence converges,
i.e. that $\lim_{k\to\infty} \textbf{x}^{(k)} = \textbf{x}_\ast$,
then
```math
\mathbf{P} \mathbf{x}_\ast = (\mathbf{P} - \mathbf{A}) \mathbf{x}_\ast + \mathbf{b}
\qquad \Leftrightarrow \qquad \textbf{0} = \mathbf{b} - \mathbf{A}\mathbf{x}_\ast,
```
i.e. the limit $\mathbf{x}_\ast$ is indeed the solution to the original equation.
"""

# ╔═╡ d3f1bc8e-f14a-48dd-8f24-c19cd4f31844
md"""
The matrix $\textbf{P}$ is usually referred
to as the **preconditioner**. Since we need to solve a linear
system $\textbf{P} \textbf{u}^{(k)} = \textbf{r}^{(k)}$ in *each
iteration* $k$, an important criterion for a good preconditioner
is that solving such linear systems is cheap.

Let us rewrite equation (2) as
```math
\mathbf{P} \underbrace{\left(\textbf{x}^{(k+1)} - \textbf{x}^{(k)} \right)}_{=\textbf{u}^{(k)}} = \underbrace{\textbf{b} - \textbf{A} \textbf{x}^{(k)}}_{=\textbf{r}^{(k)}}.
```
The quantity $\textbf{r}^{(k)} = \textbf{b} - \textbf{A} \textbf{x}^{(k)}$
appearing in this expression is called the **residual** at iteration $k$.
It measures how far the iterate $\textbf{x}^{(k)}$
is from being a solution to $\mathbf{A} \mathbf{x} = \mathbf{b}$.
This can be understood by observing that 
$\textbf{r}^{(k)} = \textbf{0}$ indeed implies
$\mathbf{A} \mathbf{x}^{(k)} = \mathbf{b}$,
i.e. that $\textbf{x}^{(k)}$ is the solution.

We formulate Richardson's iterations:

!!! info "Algorithm 1: Richardson iteration"
    Let $\mathbf{A} \in \mathbb{R}^{n\times n}$ be a system matrix,
    right-hand side $\mathbf{b} \in \mathbb{R}^n$, an invertible
    preconditioner $\mathbf{P} \in \mathbb{R}^{n\times n}$
    as well as the desired convergence tolerance $ε$.
    For $k = 0, 1, \ldots$ iterate:
    1. Compute residual $\mathbf{r}^{(k)} = \mathbf{b} - \mathbf{A} \mathbf{x}^{(k)}$
    2. Solve linear system $\mathbf{P} \mathbf{u}^{(k)} = \mathbf{r}^{(k)}$
       for the update $\mathbf{u}^{(k)}$.
    3. Compute new $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \mathbf{u}^{(k)}$.
    The iteration is stopped as soon as
    $\frac{\| \mathbf{r}^{(k)} \|}{\| \mathbf{b} \|} < ε$.

We will discuss the rationale behind the stopping condition
$\| \mathbf{r}^{(k)} \| / \| \mathbf{b} \| < ε$ 
in the subsection on error control below.

Note, that Algorithm 1 can be directly brought into the setting of the
fixed point methods we studied in chapter 4.
Essentially Algorithm 1 provides an iteration
$\mathbf{x}^{(k+1)} = g(\mathbf{x}^{(k)})$ with the map $g$ given by
```math
\tag{3}
g(\mathbf{x}) = \mathbf{x} + \mathbf{P}^{-1} \underbrace{\left( \textbf{b} - \mathbf{A} \mathbf{x} \right)}_{=\mathbf{r}}.
```
Indeed, its fixed point $\mathbf{x}_\ast = g(\mathbf{x}_\ast)$ necessarily satisfies $\mathbf{0} = \mathbf{r}_\ast = \textbf{b} - \mathbf{A} \mathbf{x}_\ast$ and therefore is a solution to the linear system.

An implementation of Algorithm 1 in Julia is given by
"""

# ╔═╡ 94a56970-0c97-4bc4-bbd5-d297cd54823e
function richardson(A, b, P; x=zero(b), tol=1e-6, maxiter=100)
	history = Float64[]  # Track relative residual norm
	for k in 1:maxiter
		r = b - A * x

		relnorm = norm(r) / norm(b)
		push!(history, relnorm)
		if relnorm < tol
			break
		end
		
		u = P \ r
		x = x + u
	end
	(; x, history)  # Return current iterate and history
end

# ╔═╡ 76dcdc97-b0e3-4f3f-8ebd-615df0cea57f
md"""Let us consider the test problem $\textbf{A} \textbf{x} = \textbf{b}$ with"""

# ╔═╡ 1c4ba2a0-553f-4b35-8b45-3de1c9e3a4e8
# Generate a random matrix, which has large entries on the diagonal
A = randn(100, 100) + Diagonal(10 .+ 50rand(100))

# ╔═╡ 49134b41-d3ca-4286-ae9a-5be0d1569241
b = rand(100);

# ╔═╡ e4a1e62c-646b-4378-b2a5-86cfe36bada9
md"""In this case the problem is sufficiently small that it's reference can still be computed by a direct method, e.g. the LU factorisation performed by default when employing Julia's `\` operator:"""

# ╔═╡ 9a4c02fb-6cbf-4b25-9a3a-74066ba14666
x_reference = A \ b

# ╔═╡ 65bdbf4a-b304-43d7-a3b7-a5c932b23964
md"""
To apply Richeardson iterations, we need a preconditioner $\mathbf P$.
Due to the construction of $\mathbf A$ we chose
it has its largest entries in each row on the diagonal.
A particularly simple preconditioner, which typically works well for such matrices
is to use the diagonal of $\mathbf A$ as $\mathbf P$, i.e.
"""

# ╔═╡ c9419d03-d32b-4811-aef7-644803e7d971
P = Diagonal(A)

# ╔═╡ 9efa5334-f2e9-43d6-899a-8cb90fff8c6a
md"""
Applying this preconditioner within the Richardson iterations, yields the correct result to the chosen tolerance:
"""

# ╔═╡ 33ebd4d0-c369-4994-9aab-455a4a4e9c4b
begin
	richardson_result = richardson(A, b, P; tol=1e-6)
	maximum(abs, richardson_result.x - x_reference)
end

# ╔═╡ 0b374ba9-24af-43a2-86cb-c3a4cf5d999c
let
	plot(richardson_result.history;
	     yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
	     title="Richardson convergence", label=L"$P = \textrm{Diagonal(}A\textrm{)}$")
end

# ╔═╡ 2d925aaf-5ffc-4df4-bd7e-1a45cdfcc442
md"""
## Computational cost

Our main motivation for considering iterative methods was to overcome the $O(n^3)$ computational cost of Gaussian elimination / LU factorisation for solving linear systems.

Let us revist Richardson iteration (Algorithm 1) and discuss the computational cost
for the case where $\mathbf{A} \in \mathbb{R}^{n\times n}$ is a full matrix.
Recall that the cost of LU factorisation scaled as $O(n^3)$ for this setting.

1. In the first step the most costly operation is the computation of the matrix-vector product $\mathbf{A} \mathbf{x}^{(k)}$. For a full matrix this costs $O(n^2)$.
2. In the second step we solve $\mathbf{P} \mathbf{u}^{(k)} = \mathbf{r}^{(k)}$. If $\mathbf{P}$ is a diagonal matrix (like we just considered), this costs $O(n)$,
   whereas for a triangular matrix (where one would perform backward or forward substitution) the cost is $O(n^2)$.
3. The computation of the new $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \mathbf{u}^{(k)}$ again costs only $O(n)$.

Overall for *full matrices* the cost of Richardson iterations is thus only $O(n^2)$.

For *sparse matrices* $\mathbf{A} \mathbf{x}^{(k)}$ can even be done in $O(n)$ computational cost. Similarly using clever preconditioners step 2 can be done in $O(n)$ time. Examples would be a sparsity-adapted forward substitution algorithm or again a diagonal preconditioner. Overall Richardson iterations thus only cost $O(n)$ --- in stark contrast to the $O(n^3)$ cost also required for LU factorisation in this setting.

!!! note "Observation: Computational cost of iterative methods"
    When solving linear systems using iterative methods,
    the computational cost (per iteration) scales
    - for **full matrices** $\mathbf{A} \in \mathbb{R}^{n\times n}$ as $O(n^2)$
    - for **sparse matrices** as $O(n)$.
    A recurring theme is that the cost of the matrix-vector product
    $\mathbf{A} \mathbf{x}^{(k)}$ essentially determis the
    cost of the iterative scheme.
"""

# ╔═╡ 4c88dc9e-3629-4188-bac7-55e3e8276091
md"""
## Convergence analysis

Having established above that Richardson iteration indeed leads to numerically the correct answer, we proceed to analyse its convergence.

Following our previous remark that Algorithm 1 can be brought into the setting of a fixed point method $\mathbf{x}^{(k+1)} = g(\mathbf{x}^{(k)})$ with
```math
g(\mathbf{x}) = \mathbf{x} + \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x} \right),
```
we study the evolution of the error $\mathbf{e}^{(k)} = \mathbf{x}^{(k)} - \mathbf{x}_\ast$, where $\mathbf{x}_\ast = g(\mathbf{x}_\ast)$ is the fixed point.
Subtracting $\mathbf{x}_\ast$ from both sides of $\mathbf{x}^{(k+1)} = g(\mathbf{x}^{(k)})$  we find
```math
\tag{4}
\begin{aligned}
\mathbf{e}^{(k+1)} &= g(\mathbf{x}^{(k)}) - \mathbf{x}_\ast \\
&= g(\mathbf{x}^{(k)}) - g(\mathbf{x}_\ast) \\
&= \mathbf{x}^{(k)} + \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x}^{(k)} \right) - \mathbf{x}_\ast - \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x}_\ast \right) \\
&= \mathbf{e}^{(k)} +  \mathbf{P}^{-1} \left(
\textbf{b} - \mathbf{A} \mathbf{x}^{(k)} 
-  \textbf{b} + \mathbf{A} \mathbf{x}_\ast \right) \\
&= \mathbf{e}^{(k)} -  \mathbf{P}^{-1} \mathbf{A} \left(\mathbf{x}^{(k)} - \mathbf{x}_\ast \right) \\
&= \underbrace{ \left(\mathbf{I} - \mathbf{P}^{-1} \mathbf{A} \right) }_{= \mathbf{B}}\,\mathbf{e}^{(k)},
\end{aligned}
```
where $\mathbf{I}$ is the identity matrix and the matrix $\mathbf{B}$ is usually called the **iteration matrix**.

Recall the definition of the matrix norm for a matrix $M \in \mathbb{R}^{n\times n}$,
namely
```math
\|\mathbf{M}\| = \max_{\mathbf{0}\neq\mathbf{v}\in\mathbb{R}^n} \frac{\| \mathbf{M}\mathbf{v}\|}{\|\mathbf{v}\|}.
```
Using this and (4) we thus have
```math
\begin{aligned}
\|\mathbf{e}^{(k+1)}\| &= \left\|\mathbf{B}\mathbf{e}^{(k)}\right\| \\
	&≤ \left\|\mathbf{B}\right\| \left\|\mathbf{e}^{(k)}\right\| \\
	&≤ \left\|\mathbf{B}\right\|^2 \left\|\mathbf{e}^{(k-1)}\right\| \\
	&\ \vdots \\
	&≤ \left\|\mathbf{B}\right\|^{k+1} \left\|\mathbf{e}^{(0)}\right\| \\
\end{aligned}
```
"""

# ╔═╡ 5b517ed9-cbcf-4b52-bdcf-c8269690736b
md"""
Provided $\|\mathbf{B}\| < 1$ the iteration thus converges. We summarise in a theorem:

!!! info "Theorem 1"
    Given a system matrix $\mathbf{A} \in \mathbb{R}^{n \times n}$,
    RHS $\mathbf{b} \in \mathbb{R}^n$ and preconditioner
    $\mathbf{P} \in \mathbb{R}^{n \times n}$ invertible,
    the Richardson method with iteration matrix
    $\mathbf{B} = \mathbf{I} - \mathbf{P}^{-1} \mathbf{A}$
    converges for *any* initial guess $\mathbf{x}^{(0)}$
    if $\left\|\mathbf{B}\right\| < 1$.
    Moreover
    ```math
    \|\mathbf{x}_\ast - \mathbf{x}^{(k)}\| ≤ \left\|\mathbf{B}\right\|^k \|\mathbf{x}_\ast - \mathbf{x}^{(0)}\|.
    ```

Notice that $\mathbf{B}$ is actually just the **Jacobian**
of the fixed-point map $g$.
In fact we can thus view the above theorem in the same light
as the condition $|g'(x)| < 1$ we previously found for the convergence
of the fixed-point methods for non-linear equations.
However, the notable difference is that Richardson iterations
converge for *any initial guess*, which for fiexd-point methods
is not true in general.

The norm of iteration matrix
  $\|\mathbf{B}\|$ is the crucial quantity to determine not only *if*
Richardson iterations converge, but also at which speed.
To understand it a little better and in particular to understand
the role of the preconditioner $\mathbf{P}$ we introduce the spectral radius:

!!! info "Definition: Spectral radius"
    Let $\mathbf{B} \in \mathbb{R}^{n \times n}$ be a diagonalisable matrix
    and denote its eigenvalues as $λ_i(\mathbf{B})$, $i = 1, \ldots n$.
    The **spectral radius** of $\mathbf{B}$ is the quantity
    ```math
    \tag{5}
    ρ(\mathbf{B}) = \max_{i=1,\ldots,n}\left|λ_i(\mathbf{B})\right|
    ```

Note that one can show that
```math
   \|\mathbf{B}\|^2 = ρ(\mathbf{B}),
```
which implies

!!! info "Theorem 2"
    Richardson iterations (Algorithm 1) converge if and only if
    $ρ(\mathbf{B}) < 1$,
    where $\mathbf{B} = \mathbf{I} - \mathbf{P}^{-1} \mathbf{A}$.

With this we can understand why one would even precondition at all.
If we were not to perform any preconditioning, i.e. $\mathbf{P} = \mathbf{I}$,
then $ρ(\mathbf{B}) = ρ\left(\mathbf{I} - \mathbf{A}\right)$
becomes
"""

# ╔═╡ 54cb29d7-4933-4ab4-a3b5-26431d8432f9
let
	B = I - A    # Iteration matrix for P = I, i.e. no preconditioner
	ρ = maximum(abs.(eigvals(B)))
end

# ╔═╡ af630960-a734-4124-81d3-a2657456f6c2
md"""
whereas when using a diagonal preconditioner `P = Diagonal(A)`,
i.e. $ρ(\mathbf{B}) = ρ\left(\mathbf{I} - \mathbf{P}^{-1} \mathbf{A}\right)$,
then
"""

# ╔═╡ 51030114-04d7-49eb-9ac1-042941681446
let
	P = Diagonal(A)
	B = I - inv(P) * A  # Iteration matrix for P = Diagonal(A)
	ρ = maximum(abs.(eigvals(B)))
end

# ╔═╡ 114bc28c-7888-4c6c-85bd-8c6232624491
md"""
the result is that only the preconditioned iterations converge:
"""

# ╔═╡ d8e0e584-985d-4ff0-8a5c-b21c83c226d9
let
	P = Diagonal(A)
	richardson_diagonal = richardson(A, b, P; tol=1e-6, maxiter=30)
	richardson_no_preconditioner = richardson(A, b, I; tol=1e-6, maxiter=30)
	
	plot(richardson_diagonal.history;
	     yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
	     label=L"$P = \textrm{Diagonal(}A\textrm{)}$", legend=:bottomright)
	plot!(richardson_no_preconditioner.history;
	      label=L"$P = I$", lw=2, mark=:o)
end

# ╔═╡ eaee3d23-bee2-4f22-bfcf-4b1e890033b2
md"""
## Choosing a good preconditioner

From the above experiments we notice:

!!! info "Observation: Preconditioners"
    A good preconditioner $P$ in the Richardson iterations, satisfies the following properties:
    1. It is *cheap to invert*, that is linear systems $\mathbf{P} \mathbf{u} = \mathbf{r}$ are cheap to solve.
    2. The spectral radius $ρ(\mathbf I - \mathbf P^{-1} \mathbf A)$ is as small as possible,
      definitely smaller than one (to ensure convergence).

Clearly condition 2. suggests that the *ideal preconditioner* is $P = A$, such that $\mathbf{P}^{-1} = \mathbf{A}^{-1}$. In this setting $ρ(\mathbf I - \mathbf{P}^{-1} \mathbf{A}) = ρ(\mathbf I - \mathbf I) = 0$. The trouble of this choice is that step 2 of Algorithm 1 (Richardson iteration) now effectively requires to solve the system $\mathbf A \mathbf u^{(k)} = \mathbf r^{(k)}$ for $\mathbf u^{(k)}$, i.e. we are back to solving the original problem. 

On the other hand condition 1. suggests to use $\mathbf P^{-1} = \mathbf I^{-1}$ (since the identity is cheapest to invert --- nothing needs to be done). However, this does not really do anything and does thus not improve the spectral radius.

In practice we thus need to find a compromise between the two. As mentioned above standard choices for the preconditioner are:
- the diagonal of $\mathbf A$
- the lower triangle of $\mathbf A$
- another common choice is the *incomplete LU factorisation*, i.e. where one seeks a factorisation $\tilde{\mathbf L} \tilde{\mathbf U} \approx \mathbf A$ with a lower triangular matrix $\tilde{\mathbf L}$ and an upper triangular $\tilde{\mathbf U}$, which only *approximately* factorises $\mathbf A$. This gives additional freedom in designing the factorisations, in particular to avoid fill-in for sparse matrices.
- in many physics or engineering applications $\mathbf A$ results from a physical model. A preconditioner $\mathbf P$ can thus be resulting from an approximation to the employed physical model (e.g. by dropping the modelling of some physical effects).

We will not discuss the art of designing good preconditioners any further at this stage. Interested readers are referred to the excellent book [Youssef Saad *Iterative methods for Sparse Linear Systems*, SIAM (2003)](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf).
"""

# ╔═╡ 1dca7d0a-53f7-4903-b390-b5cc34a1f319
md"""
### Error control and stopping criterion

Let us return to clarifying the choice of the stopping criterion in the Richardson iterations (Algorithm 1).

As usual our goal is to control the error $\|\mathbf{x}_\ast - \mathbf{x}^{(k)}\|$
based on our stopping criterion,
but without having access to the exact solution $x_\ast$.
However, we know that $\mathbf{A} \mathbf{x}_\ast = \mathbf{b}$,
since $\mathbf{x}_\ast$ is just the solution to the linear system
we want to solve. 

Similarly $\mathbf{r}^{(k)} =  \mathbf{b} - \mathbf{A} \mathbf{x}^{(k)}$
directly implies
```math
\mathbf{A} \mathbf{x}^{(k)} = \underbrace{\mathbf{b} - \mathbf{r}^{(k)}}_{=\widehat{\mathbf{b}}}
```
where we defined the "perturbed" RHS $\widehat{\mathbf{b}}$.
Notably $\mathbf{x}^{(k)}$ is the *exact* solution of a linear system
involving $\mathbf{A}$ and this perturbed RHS.
Using further that
```math
\tag{6}
\|\mathbf{b}\| = \|\mathbf{A} \mathbf{x}_\ast\| ≤ \|\mathbf{A}\| \, \|\mathbf{x}_\ast\| 
\quad
\Rightarrow
\quad \|\mathbf{x}_\ast\| ≥ \frac{\|\mathbf{b}\|}{\|\mathbf{A}\|}
```
we develop for the relative error
```math
\begin{aligned}
\frac{\|\mathbf{x}_\ast - \mathbf{x}^{(k)} \|}{\| \mathbf{x}_\ast \|}
\stackrel{(6)}{≤} \frac{ \| \mathbf{A} \| \, \left\|\mathbf{A}^{-1} \left(\mathbf{b} - \widehat{\mathbf{b}}\right) \right\|}{\| \mathbf{b} \|}
≤ \| \mathbf{A} \| \, \| \mathbf{A}^{-1} \| \, \frac{ \left\| \mathbf{b} - \widehat{\mathbf{b}} \right\| }{\| \mathbf{b} \|}
=  \| \mathbf{A} \| \, \| \mathbf{A}^{-1} \| \, \frac{ \left\| \mathbf{r}^{(k)} \right\| }{\| \mathbf{b} \|}.
\end{aligned}
```
Recall the previously discussed condition number
```math
κ(\mathbf{A}) = \| \mathbf{A} \| \, \| \mathbf{A}^{-1} \|
= \sqrt{ \frac{λ_\text{max}(\mathbf{A}^T\mathbf{A} ) }{ λ_\text{min}(\mathbf{A}^T \mathbf{A})  }  }
```
as well as our stopping criterion from Algorithm 1, that is
```math
\frac{\left\|\mathbf{r}^{(k)} \right\|}{\| \mathbf{b} \|} < ε,
```
With these quantities we rewrite
```math
\frac{\|\mathbf{x}_\ast - \mathbf{x}^{(k)} \|}{\| \mathbf{x}_\ast \|}
≤ κ(\mathbf{A}) \frac{\left\|\mathbf{r}^{(k)} \right\|}{\| \mathbf{b} \|} = κ(\mathbf{A}) \, ε.
```
In other words our stopping criterion ensures that the relative error of the returned solution is no larger than $κ(\mathbf{A})$ times the chosen tolerance.

If the matrix is well-conditioned, i.e. $κ(\mathbf{A})$ is close to $1$,
then the relative residual $\frac{\left\|\mathbf{r}^{(k)} \right\|}{\| \mathbf{b} \|}$
provides a good estimate of the true error and our stopping criterion is appropriate.
However, if $κ(\mathbf{A})$ is large, even a small residual may imply a large error in the returned solution.
"""

# ╔═╡ 858b7c87-80de-4c1a-8598-03657608e323
md"""
## Jacobi and Gauss-Seidel method

We will now briefly discuss the Jacobi and Gauss-Seidel methods,
which can be seen as particular cases of Richardson iterations.
"""

# ╔═╡ d3f4ef3a-1a17-4078-b7a8-9e9b767ee541
md"""
Recall equation (2)
```math
\mathbf{P} \mathbf{x}^{(k+1)} = (\mathbf{P} - \mathbf{A}) \mathbf{x}^{(k)} + \mathbf{b},
\qquad k = 0, 1, \ldots.
```

In the **Jacobi method** the key idea is to use the *diagonal* of the matrix $\mathbf{A}$ as the preconditioner
```math
\mathbf{P} = \begin{pmatrix}
A_{11} &        &        & \\
       & A_{22} &        & \\
       &        & \ddots &  \\
       &        &        & A_{nn}
\end{pmatrix}.
```
and to explicitly insert this expression into equation (2).
Rewriting we obtain in the $(k+1)$-st iteration of Jacobi's method
```math
\begin{pmatrix}
A_{11} &        &        & \\
       & A_{22} &        & \\
       &        & \ddots &  \\
       &        &        & A_{nn}
\end{pmatrix}
\begin{pmatrix}
x_1^{(k+1)} \\ x_2^{(k+1)} \\ \vdots \\ x_n^{(k+1)}
\end{pmatrix}
= -\begin{pmatrix}
0 & A_{12} & \cdots & A_{1n} \\
A_{21}       & 0      &        & \vdots \\
\vdots       &        & \ddots &  \\
A_{n1}       & \cdots & A_{n,n-1} & 0
\end{pmatrix}
\begin{pmatrix}
x_1^{(k)} \\ x_2^{(k)} \\ \vdots \\ x_n^{(k)}
\end{pmatrix}
+
\begin{pmatrix}
b_1 \\ b_2 \\ \vdots \\ b_n
\end{pmatrix}
```
The diagonal structure of $\textbf{P}$ allows the explicit
computation of $\textbf{x}^{(k+1)}$. Its $i$-th component can be obtained as
```math
\tag{7}
x^{(k+1)}_i = \frac{1}{A_{ii}} \left( b_i - \sum_{\stackrel{j=1}{j\neq i}}^n A_{ij} \, x^{(k)}_j \right)
```
for all $i = 1, \ldots, n$ and as long as $A_{ii} \neq 0$.
"""

# ╔═╡ 77f1c8d5-8058-47cb-a0b5-436b4259aec9
function jacobi(A, b; x=zero(b), tol=1e-6, maxiter=100)
	history = Float64[]  # Track relative residual norm
	
	n  = length(x)
	xᵏ = x
	for k in 1:maxiter
		xᵏ⁺¹ = zeros(n)
		for i in 1:n
			Aᵢⱼxⱼ = 0.0
			for j in 1:n
				if i ≠ j
					Aᵢⱼxⱼ += A[i, j] * xᵏ[j]
				end
			end  # Loop j
			xᵏ⁺¹[i] = (b[i] - Aᵢⱼxⱼ) / A[i, i]
		end  # Loop i
		
		relnorm_rᵏ = norm(b - A * xᵏ⁺¹) / norm(b)  # Relative residual norm
		push!(history, relnorm_rᵏ)  # Push to history
		if relnorm_rᵏ < tol         # Check convergence
			break
		end

		xᵏ = xᵏ⁺¹
	end  # Loop k
	(; x=xᵏ, history)
end

# ╔═╡ 96352c57-7f86-4af8-ae6c-a3c20f190461
md"""
The **Gauss-Seidel method** employs the lower triangular part of $A$ as the preconditioner $\mathbf{P}$ (in Julia `P = LowerTriangular(A)`):
```math
\mathbf{P} = \begin{pmatrix}
A_{11} &        &        & \\
A_{21} & A_{22} &        & \\
\vdots &        & \ddots &  \\
A_{n1} & A_{n2} & \cdots & A_{nn}
\end{pmatrix}.
```
By inserting the lower-triangular $\textbf{P}$ into (2) we obtain in the $(k+1)$-st iteration of Gauss-Seidel:
```math
\begin{pmatrix}
A_{11} &        &        & \\
A_{21} & A_{22} &        & \\
\vdots &        & \ddots &  \\
A_{n1} & A_{n2} & \cdots & A_{nn}
\end{pmatrix} 
\begin{pmatrix}
x_1^{(k+1)} \\ x_2^{(k+1)} \\ \vdots \\ x_n^{(k+1)}
\end{pmatrix}
= -\begin{pmatrix}
0 & A_{12} & \cdots & A_{1n} \\
           & 0 &  \ddots      & \vdots \\
           &   & \ddots & A_{n-1,n} \\
           &   &        & 0
\end{pmatrix}
\begin{pmatrix}
x_1^{(k)} \\ x_2^{(k)} \\ \vdots \\ x_n^{(k)}
\end{pmatrix}
+
\begin{pmatrix}
b_1 \\ b_2 \\ \vdots \\ b_n
\end{pmatrix}.
```
Using forward substitution to solve this linear system
leads to the following form for the $i$-th component
of the vector $\mathbf{x}^{(k+1)}$.
```math
\tag{8}
x^{(k+1)}_i = \frac{1}{A_{ii}} \left( b_i - \sum_{j=1}^{i-1} A_{ij} \, x^{(k+1)}_j - \sum_{j=i+1}^n A_{ij} \, x^{(k)}_j \right)
```
for all $i = 1, \ldots, n$ and as long as $A_{ii} \neq 0$.

Notice that Gauss-Seidel is very similar to Jacobi's method, just with the small difference that for computing the new component $x^{(k+1)}_i$ we use the components $x^{(k+1)}_j$ for $j<i$, which have already been updated to their new values.
"""

# ╔═╡ d376fb74-6ebd-4b07-9c5c-c34e7e16c25a
function gauss_seidel(A, b; x=zero(b), tol=1e-6, maxiter=100)
	history = Float64[]  # Track relative residual norm
	
	n  = length(x)
	xᵏ = x
	for k in 1:maxiter
		xᵏ⁺¹ = zeros(n)
		for i in 1:n
			Aᵢⱼxⱼ = 0.0
			for j in 1:i-1
				Aᵢⱼxⱼ += A[i, j] * xᵏ⁺¹[j]
			end  # Loop j
			for j in i+1:n
				Aᵢⱼxⱼ += A[i, j] * xᵏ[j]
			end  # Loop j
			xᵏ⁺¹[i] = (b[i] - Aᵢⱼxⱼ) / A[i, i]
		end  # Loop i
		
		relnorm_rᵏ = norm(b - A * xᵏ⁺¹) / norm(b)  # Relative residual norm
		push!(history, relnorm_rᵏ)  # Push to history
		if relnorm_rᵏ < tol         # Check convergence
			break
		end

		xᵏ = xᵏ⁺¹
	end  # Loop k
	(; x=xᵏ, history)
end

# ╔═╡ de24e92e-faab-41fe-afc2-287835b1e015
md"""
These two cases just provide two examples of the many flavours
of Richardson iterations, which are used in practice.
A good overview provides chapter 4 of
[Youssef Saad *Iterative methods for Sparse Linear Systems*, SIAM (2003)](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf).
"""

# ╔═╡ b36bda99-fddb-422c-9ac4-cab7d3495334
md"We return to our example problem $\mathbf{A} \mathbf{x} = \mathbf{b}$ with"

# ╔═╡ bc76f57d-22ba-488b-9c57-99e0b89c09f0
A

# ╔═╡ 08ca7c51-34c9-4cff-8e11-d875f75c899e
md"and"

# ╔═╡ f4484dda-074e-4323-844a-319c69a478b0
b

# ╔═╡ 5e217a36-ceb1-4dc9-9e3c-ad89a4c83503
md"On this problem the convergence is as follows:"

# ╔═╡ d8b7f223-20fc-4052-8614-28e5f226014d
let
	result_jacobi       = jacobi(A, b; tol=1e-10, maxiter=30)
	result_gauss_seidel = gauss_seidel(A, b; tol=1e-10, maxiter=30)

	
	P = Diagonal(A)
	richardson_diagonal = richardson(A, b, P; tol=1e-6, maxiter=30)
	richardson_no_preconditioner = richardson(A, b, I; tol=1e-6, maxiter=30)
	
	plot(result_jacobi.history;
	     yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
	     label="Jacobi", legend=:topright)
	plot!(result_gauss_seidel.history;
	      label="Gauss Seidel", lw=2, mark=:o)
end

# ╔═╡ 08dc1f90-2fc6-4328-a39d-8856f215db33
md"""
## Linear systems involving symmetric positive-definite matrices

In many applications in engineering and the sciences
the arising system matrices $\mathbf{A}$ have additional properties,
which can be exploited to obtain better-suited algorithms.
An important case are the symmetric positive definite matrices.

This part of the notes only provides a condensed introduction.
A more detailed, but very accessible introduction to the steepest descent and conjugate gradient methods discussed in this section provides
[Jonathan Shewchuk *An Introduction to the Conjugate Gradient Method Without the Agonizing Pain* (1994)](https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf).

Let us first recall the definition of symmetric positive definite matrices.

!!! info "Definition: Positive definite"
    A matrix $\textbf{A} \in \mathbb{R}^{n\times n}$ is called **positive definite**
    if
    ```math
    \forall\, \mathbf{0} \neq \mathbf{v} \in \mathbb{R}^n
    \qquad \mathbf{v}^T \mathbf{A} \mathbf{v} > 0
    ```
"""

# ╔═╡ 33cb7e52-abe1-40f0-81d0-f8b00468a9ff
md"""
For **symmetric matrices** we additionally have that $\mathbf{A}^T = \mathbf{A}$.
Symmetric positive definite matrices (often abbreviated as s.p.d. matrices)
arise naturally in physical systems when looking for configurations
minimising their energy. A result underlining this construction is

!!! info "Theorem 2"
    An matrix $\textbf{A} \in \mathbb{R}^{n\times n}$ is
    s.p.d. if and only if all its eigenvalues are real and positive.
"""

# ╔═╡ b46abcc8-f739-4da0-9d63-f6673ed6884f
md"""
To every linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$
involving an s.p.d. system matrix $\mathbf{A}$ we can associate an
**energy function** $\phi$
```math
ϕ\, :\, \mathbb{R}^n \to \mathbb{R}\, :\, ϕ(\mathbf{v}) = \frac12 \mathbf{v}^T\mathbf{A}\mathbf{v} - \mathbf{v}^T \mathbf{b}
```
Using basic vector calculus one shows that
```math
\tag{9}
\begin{aligned}
\text{Gradient} && \nabla ϕ(\mathbf{x}) &= \mathbf{A} \mathbf{x} - \mathbf{b} = -\mathbf{r} && \text{\textit{(negative residual)}} \\
\text{Hessian}  && H_ϕ(\mathbf{x}) &= \mathbf{A} && \text{\textit{(system matrix)}}
\end{aligned}
```
Setting the gradient to zero, we obtain the stationary points.
The corresponding equation
```math
\mathbf{0} = \nabla ϕ(\mathbf{x}) = \mathbf{A} \mathbf{x} - \mathbf{b}
```
only has a single solution, provided that $\mathbf{A}$ is non-singular.
Since the Hessian $H_ϕ(\mathbf{x}) = \mathbf{A}$ is positive definite,
this stationary point is a minimum. As a result we obtain

!!! info "Proposition 3"
    Given an s.p.d. matrix $\mathbf{A} \in \mathbb{R}^{n \times n}$,
    the solution of the system $\mathbf{A} \mathbf{x}_\ast = \mathbf{b}$
    is the unique minimum of the function $ϕ$, i.e.
    ```math
    \tag{10}
    \mathbf{x}_\ast = \text{argmin}_{\textbf{v} \in \mathbb{R}^{n}} \, ϕ(\mathbf{v}).
    ```
"""

# ╔═╡ 590e8877-0600-4bbb-8e88-7033548815b5
md"""
## Steepest descent method
"""

# ╔═╡ fa7c4b56-c322-4662-a190-07a700947dd8
TODO("More details")

# ╔═╡ c330ec39-bcce-4b97-b5c9-11c13967de56
md"""
A Julia implementation of steepest descent is given by
"""

# ╔═╡ 0159a195-627c-4e07-834e-222248665da8
function steepest_descent_simple(A, b; x=zero(b), tol=1e-6, maxiter=100)
	history = Float64[]  # Track relative residual norm
	
	r = b - A * x
	for k in 1:maxiter
		relnorm = norm(r) / norm(b)
		push!(history, relnorm)
		if relnorm < tol
			break
		end
		
		Ar = A * r
		α = dot(r, r) / dot(r, Ar)
		x = x + α * r
		r = r - α * Ar
	end
	(; x, history)
end

# ╔═╡ 9e34bd04-a4b7-4b60-b93b-dcd37b25287b
begin
	# Generate a random s.p.d. matrix
	U, Σ, V = svd(randn(100, 100))
	Aₛ = U * Diagonal(abs.(Σ)) * U'
end

# ╔═╡ 45cfdf30-8738-4cb6-a49a-534392128475
bₛ = rand(100);

# ╔═╡ 646ebf55-d081-4884-8524-15e5b0ecebcc
md"The system matrix has a condition number larger well larger than $1$:"

# ╔═╡ d7eedfae-90ae-4917-b629-fc56f00bc3da
cond(Aₛ)

# ╔═╡ 7481b32b-d8c4-405b-9a97-dfb17ce87f53
md"""Still, as the theory predicts, the method converges. However, convergence gets slow after a while and approaches the theoretical rate"""

# ╔═╡ 1a675fdb-885d-462e-bc77-9ea074d49614
rate_steep_desc = (cond(Aₛ) - 1) / (cond(Aₛ) + 1)

# ╔═╡ f0c0071f-5499-4bfa-b838-4f3d72ae3b3d
let
	descent = steepest_descent_simple(Aₛ, bₛ; tol=1e-6, maxiter=40)
	
	plot(descent.history;
		yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||", ylims=(5e-2, 1.3),
		title="Steepest descent convergence", label="Observed convergence", c=2)
	
	plot!(1:40, x -> 0.15rate_steep_desc^x, lw=2, ls=:dash, label="Theoretical rate", c=2)
end

# ╔═╡ 8090fe02-249a-492b-b04a-b3d1409253c4
md"""
## Preconditioned steepest descent
"""

# ╔═╡ 7207e407-a42e-42c3-919c-69f8d8114f74
TODO("More details")

# ╔═╡ db4ab984-439d-4c43-8d02-1c0f288abdae
md"""
An implementation is
"""

# ╔═╡ fdbe3ba7-d1a3-43b0-bd56-ad0bc6f2200d
function steepest_descent(A, b, P; x=zero(b), tol=1e-6, maxiter=100)
	history = Float64[]
	
	r = b - A * x
	for k in 1:maxiter
		relnorm = norm(r) / norm(b)
		push!(history, relnorm)
		if relnorm < tol
			break
		end

		z = P \ r  # Solve Pz = r  for z, the main addition compared to
		           # steepest_descent_simple in this function
		Az = A * z
		α = dot(r, r) / dot(r, Az)
		x = x + α * r
		r = r - α * r
	end
	(; x, history)
end

# ╔═╡ 15f921ac-227c-49e9-8777-86ae502ebc6c
md"""
We  consider the convergence of both the preconditioned and the non-preconditioned steepest descent method on the model problem $\mathbf{A}_\text{s} \mathbf{x} = \mathbf{b}_\text{s}$.
"""

# ╔═╡ bf028cd7-c911-4329-86b0-442061a8e8c8
let
    P = Diagonal(Aₛ)
    descent_diagonal = steepest_descent(Aₛ, bₛ, P; tol=1e-6, maxiter=30)
    descent_no_preconditioner = steepest_descent_simple(Aₛ, bₛ; tol=1e-6, maxiter=30)
       
    plot(descent_diagonal.history;
         yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
         label=L"$P = \textrm{Diagonal(}A\textrm{)}$", legend=:bottomright,
         title="Steepest descent convergence", ylims=(1e-4, 1.3))
    plot!(descent_no_preconditioner.history;
          label=L"$P = I$", lw=2, mark=:o)
end

# ╔═╡ ca833c65-eba0-4344-ac9f-a1e360acc7dc
md"""
## Conjugate gradient method
"""

# ╔═╡ d698fe83-68a6-4988-b1e0-9e2dace2c6fd
md"""
We only show the conjugate gradient method *without preconditioning*.
"""

# ╔═╡ 7050ffeb-91db-456c-8f07-3ee748f5c395
TODO("Show a classic case where steepest descent fails and CG works, e.g. a Rosenbrock or a cartesian product of a steep and a shallow parabola.")

# ╔═╡ 3f5576db-8ee6-451c-b967-6a3a97626422
TODO("More details")

# ╔═╡ 9e3340b7-1ecb-4f6d-96a0-089932ec306d
md"""
For this example the theoretical rate is:
"""

# ╔═╡ 7cbc6819-875a-498f-a81a-f63652c3dc8c
rate_cg = (sqrt(cond(Aₛ)) - 1) / (sqrt(cond(Aₛ)) + 1)

# ╔═╡ cf8ae086-da8a-4fd7-83cd-167b8b8d7dd4
md"""
A (non-optimised) CG implementation is:
"""

# ╔═╡ 5291c5f3-7d79-42b4-b704-24697fcc9782
function conjugate_gradient_simple(A, b; x=zero(b), tol=1e-6, maxiter=100)
	history = Float64[]  # Track relative residual norm
	
	r = b - A * x
	p = r
	for k in 1:maxiter
		relnorm = norm(r) / norm(b)
		push!(history, relnorm)
		if relnorm < tol
			break
		end

		# Descent along conjugate direction p (instead of along r)
		Ap = A * p
		α  = dot(r, r) / dot(p, Ap)
		x  = x + α * p
		r_new = r - α * Ap

		# Update conjugate direction p
		β = dot(r_new, r_new) / dot(r, r)
		p = r_new + β * p
		r = r_new
	end
	(; x, history)
end

# ╔═╡ d0dfcbe6-6115-4f87-a310-5db427d90346
let
	P = Diagonal(A)
	descent_no_preconditioner = steepest_descent_simple(Aₛ, bₛ; tol=1e-6, maxiter=40)
	cg_no_preconditioner = conjugate_gradient_simple(Aₛ, bₛ; tol=1e-6, maxiter=40)
	
	plot(descent_no_preconditioner.history;
	     yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
	     label="Steepest descent", legend=:bottomleft, c=2,  ylims=(1e-4, 1.3))
	plot!(cg_no_preconditioner.history; label="CG", lw=2, mark=:o, c=1)


	plot!(1:40, x -> 0.15rate_steep_desc^x, lw=2, ls=:dash, label="Rate steep. desc.", c=2)
	plot!(1:40, x -> 0.5rate_cg^x, lw=2, ls=:dash, label="Rate CG", c=1)
end

# ╔═╡ e566f53e-e841-40d4-914d-de7f660a1e67
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 550)
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
Plots = "~1.40.0"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "52bb0330f5adf3a5307f53920d9dad77ae913952"

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
git-tree-sha1 = "522b8414d40c4cbbab8dee346ac3a09f9768f25d"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.5+0"

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
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

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
# ╟─63bb7fe9-750f-4d2f-9d18-8374b113373e
# ╠═0501f270-c39c-11ee-2077-7f7f409b6c66
# ╟─7d9c9392-3aec-4efd-a9ba-d8965687b163
# ╟─4e0e5311-3428-44f5-a616-68cf7634214f
# ╟─b5bf49df-f6ff-49f9-b016-7f0159b4af7d
# ╟─d3f1bc8e-f14a-48dd-8f24-c19cd4f31844
# ╠═94a56970-0c97-4bc4-bbd5-d297cd54823e
# ╟─76dcdc97-b0e3-4f3f-8ebd-615df0cea57f
# ╠═1c4ba2a0-553f-4b35-8b45-3de1c9e3a4e8
# ╠═49134b41-d3ca-4286-ae9a-5be0d1569241
# ╟─e4a1e62c-646b-4378-b2a5-86cfe36bada9
# ╠═9a4c02fb-6cbf-4b25-9a3a-74066ba14666
# ╟─65bdbf4a-b304-43d7-a3b7-a5c932b23964
# ╠═c9419d03-d32b-4811-aef7-644803e7d971
# ╟─9efa5334-f2e9-43d6-899a-8cb90fff8c6a
# ╠═33ebd4d0-c369-4994-9aab-455a4a4e9c4b
# ╠═0b374ba9-24af-43a2-86cb-c3a4cf5d999c
# ╟─2d925aaf-5ffc-4df4-bd7e-1a45cdfcc442
# ╟─4c88dc9e-3629-4188-bac7-55e3e8276091
# ╟─5b517ed9-cbcf-4b52-bdcf-c8269690736b
# ╠═54cb29d7-4933-4ab4-a3b5-26431d8432f9
# ╟─af630960-a734-4124-81d3-a2657456f6c2
# ╠═51030114-04d7-49eb-9ac1-042941681446
# ╟─114bc28c-7888-4c6c-85bd-8c6232624491
# ╠═d8e0e584-985d-4ff0-8a5c-b21c83c226d9
# ╟─eaee3d23-bee2-4f22-bfcf-4b1e890033b2
# ╟─1dca7d0a-53f7-4903-b390-b5cc34a1f319
# ╟─858b7c87-80de-4c1a-8598-03657608e323
# ╟─d3f4ef3a-1a17-4078-b7a8-9e9b767ee541
# ╠═77f1c8d5-8058-47cb-a0b5-436b4259aec9
# ╟─96352c57-7f86-4af8-ae6c-a3c20f190461
# ╠═d376fb74-6ebd-4b07-9c5c-c34e7e16c25a
# ╟─de24e92e-faab-41fe-afc2-287835b1e015
# ╟─b36bda99-fddb-422c-9ac4-cab7d3495334
# ╠═bc76f57d-22ba-488b-9c57-99e0b89c09f0
# ╟─08ca7c51-34c9-4cff-8e11-d875f75c899e
# ╠═f4484dda-074e-4323-844a-319c69a478b0
# ╟─5e217a36-ceb1-4dc9-9e3c-ad89a4c83503
# ╠═d8b7f223-20fc-4052-8614-28e5f226014d
# ╟─08dc1f90-2fc6-4328-a39d-8856f215db33
# ╟─33cb7e52-abe1-40f0-81d0-f8b00468a9ff
# ╟─b46abcc8-f739-4da0-9d63-f6673ed6884f
# ╟─590e8877-0600-4bbb-8e88-7033548815b5
# ╠═fa7c4b56-c322-4662-a190-07a700947dd8
# ╟─c330ec39-bcce-4b97-b5c9-11c13967de56
# ╠═0159a195-627c-4e07-834e-222248665da8
# ╠═9e34bd04-a4b7-4b60-b93b-dcd37b25287b
# ╠═45cfdf30-8738-4cb6-a49a-534392128475
# ╟─646ebf55-d081-4884-8524-15e5b0ecebcc
# ╠═d7eedfae-90ae-4917-b629-fc56f00bc3da
# ╟─7481b32b-d8c4-405b-9a97-dfb17ce87f53
# ╠═1a675fdb-885d-462e-bc77-9ea074d49614
# ╠═f0c0071f-5499-4bfa-b838-4f3d72ae3b3d
# ╟─8090fe02-249a-492b-b04a-b3d1409253c4
# ╠═7207e407-a42e-42c3-919c-69f8d8114f74
# ╟─db4ab984-439d-4c43-8d02-1c0f288abdae
# ╠═fdbe3ba7-d1a3-43b0-bd56-ad0bc6f2200d
# ╟─15f921ac-227c-49e9-8777-86ae502ebc6c
# ╠═bf028cd7-c911-4329-86b0-442061a8e8c8
# ╟─ca833c65-eba0-4344-ac9f-a1e360acc7dc
# ╠═d698fe83-68a6-4988-b1e0-9e2dace2c6fd
# ╠═7050ffeb-91db-456c-8f07-3ee748f5c395
# ╠═3f5576db-8ee6-451c-b967-6a3a97626422
# ╟─9e3340b7-1ecb-4f6d-96a0-089932ec306d
# ╠═7cbc6819-875a-498f-a81a-f63652c3dc8c
# ╟─cf8ae086-da8a-4fd7-83cd-167b8b8d7dd4
# ╠═5291c5f3-7d79-42b4-b704-24697fcc9782
# ╠═d0dfcbe6-6115-4f87-a310-5db427d90346
# ╟─e566f53e-e841-40d4-914d-de7f660a1e67
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
