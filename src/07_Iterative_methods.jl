### A Pluto.jl notebook ###
# v0.20.5

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

Both can make the **solution of linear systems prohibitively expensive for large matrices**. In this notebook we will develop *iterative methods*.
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
is Richardson's method. The idea is to introduce a **matrix splitting** $\mathbf{A} = \mathbf{P} - \mathbf{R}$ where $\mathbf{P}$ is an invertible matrix
and $\mathbf{R} = \mathbf{P} - \mathbf{A}$.
With this idea we rewrite $\mathbf{A} \mathbf{x} = \mathbf{b}$ to
```math
(\mathbf{P} - \mathbf{R})\,\mathbf{x} = \mathbf{b}
\qquad \Leftrightarrow \qquad
\mathbf{P} \mathbf{x} = (\mathbf{P} - \mathbf{A})\, \mathbf{x} + \mathbf{b}
\qquad \Leftrightarrow \qquad
\mathbf{x} = \mathbf{P}^{-1} \left[(\mathbf{P} - \mathbf{A}) \mathbf{x} + \mathbf{b}\right]
```
If we define $g(\mathbf{x}) = \mathbf{P}^{-1} \left[(\mathbf{P} - \mathbf{A}) \mathbf{x} + \mathbf{b}\right]$ we observe this to be exactly the setting of
the multi-dimensional fixed-point methods we developed
in chapter 4,
i.e. our goal is to obtain a solution $\mathbf{x}_\ast$
with $\mathbf{x}_\ast = g(\mathbf{x}_\ast) = \mathbf{P}^{-1} \left[(\mathbf{P} - \mathbf{A}) \mathbf{x}_\ast + \mathbf{b}\right]$.

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

# ╔═╡ 1292abe9-5d88-45b0-8904-b89a03009968
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
appearing in this expression is the **residual** at iteration $k$.
As usual it measures how far the iterate $\textbf{x}^{(k)}$
is from being a solution to $\mathbf{A} \mathbf{x} = \mathbf{b}$:
if $\textbf{r}^{(k)} = \textbf{0}$, then indeed
$\mathbf{A} \mathbf{x}^{(k)} = \mathbf{b}$,
i.e. that $\textbf{x}^{(k)}$ is the solution.
"""

# ╔═╡ c6617ee2-be59-49b2-b0bc-7423b49c2fce
md"""
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

Comparing with our discussion on the 
[fixed point methods we studied in chapter 4](https://teaching.matmat.org/numerical-analysis/04_Nonlinear_equations.html)
we notice that Algorithm 1 is essentially fixed-point iteration
$\mathbf{x}^{(k+1)} = g(\mathbf{x}^{(k)})$ with the map $g$ given by
```math
\tag{3}
g(\mathbf{x}) = \mathbf{x} + \mathbf{P}^{-1} \underbrace{\left( \textbf{b} - \mathbf{A} \mathbf{x} \right)}_{=\mathbf{r}}
```
and as discussed above the fixed point  $\mathbf{x}_\ast = g(\mathbf{x}_\ast)$
necessarily satisfies $\mathbf{0} = \mathbf{r}_\ast = \textbf{b} - \mathbf{A} \mathbf{x}_\ast$ and therefore is a solution to the linear system.

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
A = randn(100, 100) + Diagonal(15 .+ 50rand(100))

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

# ╔═╡ 03984b5d-3fe9-4557-9a0b-ddcc7d0a6c34
md"""
## Computational cost

Our main motivation for considering iterative methods was to overcome the $O(n^3)$ computational cost of Gaussian elimination / LU factorisation for solving linear systems.

Let us revist Richardson iteration (Algorithm 1) and discuss the computational cost
for the case where $\mathbf{A} \in \mathbb{R}^{n\times n}$ is a full matrix.
Recall that the cost of LU factorisation scaled as $O(n^3)$ for this setting.

1. In the **first step** the most costly operation is the **computation of the matrix-vector product** $\mathbf{A} \mathbf{x}^{(k)}$. For a full matrix this costs $O(n^2)$.
2. In the **second step** we solve $\mathbf{P} \mathbf{u}^{(k)} = \mathbf{r}^{(k)}$. If $\mathbf{P}$ is a **diagonal matrix** (like we just considered), this costs $O(n)$,
   whereas for a **triangular matrix** (where one would perform backward or forward substitution) the cost is $O(n^2)$.
3. The computation of the new $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \mathbf{u}^{(k)}$ again costs only $O(n)$.

Overall **for full matrices** the **cost of Richardson iterations** is thus only $O(n^2)$.
"""

# ╔═╡ 5913dace-5e27-4d4d-8bc3-0b4118b1b093
md"""
For **sparse matrices** $\mathbf{A} \mathbf{x}^{(k)}$ can computed in $O(n)$ computational cost. Similarly using clever preconditioners step 2 can be done in $O(n)$ time. Examples would be a sparsity-adapted forward substitution algorithm or again a diagonal preconditioner. Overall Richardson iterations thus **only cost $O(n)$** --- in stark contrast to the $O(n^3)$ cost also required for LU factorisation in this setting.
"""

# ╔═╡ 7b8bb857-35d3-4014-851b-9ab263f48967
md"""
!!! note "Observation: Computational cost of iterative methods"
    When solving linear systems using iterative methods,
    the computational cost (per iteration) scales
    - for **full matrices** $\mathbf{A} \in \mathbb{R}^{n\times n}$ as $O(n^2)$
    - for **sparse matrices** as $O(n)$.
    A recurring theme is that the cost of the matrix-vector product
    $\mathbf{A} \mathbf{x}^{(k)}$ essentially determis the
    cost of the iterative scheme.
"""

# ╔═╡ 1f7765a2-4986-4ef2-849e-ca1efb59abcf
md"""
## Convergence analysis

Having established above that Richardson iteration indeed leads to numerically the correct answer, we proceed to analyse its convergence.

As an iterative method Algorithm 1 can be brought into the setting of a fixed point method $\mathbf{x}^{(k+1)} = g(\mathbf{x}^{(k)})$ by choosing
```math
g(\mathbf{x}) = \mathbf{x} + \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x} \right).
```
In line with our discussion in [Root finding and fixed-point problems](https://teaching.matmat.org/numerical-analysis/04_Nonlinear_equations.html) the convergence depend on the properties of the **Jacobian** (i.e the derivative) of this map at the fixed point. The term $\mathbf{P}^{-1} \mathbf{b}$ does not depend on $\mathbf{x}$ and thus drops, for the rest we compute:
```math
\begin{align*}
g'(\mathbf{x}_\ast) &= \underbrace{\mathbf{I} - \mathbf{P}^{-1} \mathbf{A}}_{=\mathbf{B}}
\end{align*}
```
where the matrix $\mathbf{B}$ is usually referred to as the **iteration matrix**.
"""

# ╔═╡ ba77a115-f9fb-4f80-a0ea-eb50ddd9ae65
md"""
Previously we considered scalar-valued functions where the condition for convergence was that $|g'(x_\ast)| < 1$. A careful analysis in fact reveals that for vector-valued fixed-point maps like in the case of Richardson iterations this condition becomes
```math
\| g'(\mathbf{x}_\ast) \| = \| \mathbf{B} \| < 1,
```
i.e. the modulus is simply replaced by the matrix norm.

As a reminder, recall that the definition of the matrix norm for a matrix $M \in \mathbb{R}^{n\times n}$ is
```math
\|\mathbf{M}\| = \max_{\mathbf{0}\neq\mathbf{v}\in\mathbb{R}^n} \frac{\| \mathbf{M}\mathbf{v}\|}{\|\mathbf{v}\|}.
```
"""

# ╔═╡ 6f3d76f5-aeec-4a05-8aea-a4f97d625a4a
md"""
We summarise in a theorem:

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
	This is **linear convergence with rate**  $\| \mathbf{B}\|$.


Notice that the theorem mentions that Richardson iterations
converges for *any initial guess*, which for general fixed-point methods
is not true.
"""

# ╔═╡ 28d6cb69-71ab-493f-a9a1-60f2fe37c1ed
md"""
This is in fact a consequence of the fact that $g$ is a linear function,
i.e. that its Taylor expansion
```math
g(\mathbf x^{(k)}) = g(\mathbf{x}_\ast) + g'(\mathbf{x}_\ast) (\mathbf{x}^{(k)} - \mathbf{x}_\ast)
```
terminates after the linear term as we easily verify by inserting our obtained expressions into the right-hand side:
```math
\begin{aligned}
g(\mathbf{x}_\ast) + g'(\mathbf{x}_\ast) (\mathbf{x}^{(k)} - \mathbf{x}_\ast)
&= \mathbf{x}_\ast + \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x}_\ast \right) + \left(\mathbf{I} - \mathbf{P}^{-1} \mathbf{A} \right) (\mathbf{x}^{(k)} - \mathbf{x}_\ast) \\
&=\mathbf{x}_\ast + \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x}_\ast \right) + \left(\mathbf{x}^{(k)} - \mathbf{x}_\ast - \mathbf{P}^{-1} \mathbf{A}\mathbf{x}^{(k)} + \mathbf{P}^{-1} \mathbf{A}\mathbf{x}_\ast\right) \\
&= \mathbf{x}^{(k)} +  \mathbf{P}^{-1} \mathbf{b} - \mathbf{P}^{-1} \mathbf{A}\mathbf{x}^{(k)} \\
&= g(\mathbf{x}^{(k)})
\end{aligned}
```
Denoting the error in the $k$-th iteration as
$\mathbf{e}^{(k)} = \mathbf{x}^{(k)} - \mathbf{x}_\ast$
one can thus show that
```math
\tag{4}
\mathbf{e}^{(k+1)} = \underbrace{ \left(\mathbf{I} - \mathbf{P}^{-1} \mathbf{A} \right) }_{= \mathbf{B}  \,=\, g'(\mathbf{x}_\ast)} \,\mathbf{e}^{(k)},
```
such that if $\left\|\mathbf{B}\right\| < 1$ the error is *guaranteed to shrink* in an iteration, independent on the current iterate $\mathbf{x}^{(k)}$.
"""

# ╔═╡ 6a5346df-6c21-47fb-9f99-d0bc75532135
md"""
This is in contrast to the Convergence analysis discussion in [Root finding and fixed-point problems](https://teaching.matmat.org/numerical-analysis/04_Nonlinear_equations.html)
where for the case of *non-linear* fixed point functions we found that
```math
e^{(k+1)} = g'(x_\ast) e^{(k)} + O(|e^{(k)}|^2) 
```
where thus this was *not* an *exact equality* and one can thus only ensure
the error to shrink for $k$ to $k+1$ if the higher-order terms are small enough,
i.e. if $|e^{(k)}|$ is small enough, i.e. if one starts sufficiently close to the fixed point $x_\ast$.
"""

# ╔═╡ f2d4fbe9-e19a-41e0-9f3b-2f77ef0b12a7
md"""
!!! exercise
	Show that for the fixed-point map of Richardson iteration
	$g(\mathbf{x}) = \mathbf{x} + \mathbf{P}^{-1} \left( \textbf{b} - \mathbf{A} \mathbf{x} \right)$ we indeed have that
	$\mathbf{e}^{(k+1)} = \mathbf{B} \,\mathbf{e}^{(k)}$
	where $\mathbf{B} = \mathbf{I} - \mathbf{P}^{-1} \mathbf{A}$.
	Making use of the inequality
	```math
	\| \mathbf{M}\mathbf{v}\| \leq \|\mathbf{M}\| \, \|\mathbf{v}\|
	```
	and by adapting our arguments in the
	convergence analysis discussion of [Root finding and fixed-point problems](https://teaching.matmat.org/numerical-analysis/04_Nonlinear_equations.html)
	show that for Richardson iteration the error always tends to zero
	as the iteration progresses, i.e. that $\lim_{k\to\infty} \|\mathbf{e}^{(k)}\| = 0$.
"""

# ╔═╡ ef1ecb29-3ca0-4505-b96f-a56abd961979
md"""
From Theorem 1 we take away that the norm of iteration matrix
$\|\mathbf{B}\|$ is the crucial quantity to determine not only *if*
Richardson iterations converge, but also *at which rate*.
Recall in Lemma 4 of [Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.html)
we had the result that for any matrix  $\mathbf{B} \in \mathbb{R}^{m \times n}$
```math
\tag{5}
\|\mathbf{B}\| = \sqrt{ \lambda_\text{max}(\mathbf{B}^T\mathbf{B}) }.
```
With this and the condition $\|\mathbf{B}\| < 1$ for convergence
in Theorem we can understand the **role of preconditioning**.
If we were not to perform any preconditioning, i.e. $\mathbf{P} = \mathbf{I}$,
then $\|\mathbf{B}\| = \|\mathbf{I} - \mathbf{A}\|$
becomes
"""

# ╔═╡ e277ff4c-87a8-440d-bba4-f30a20b47788
A

# ╔═╡ 54cb29d7-4933-4ab4-a3b5-26431d8432f9
let
	B = I - A    # Iteration matrix for P = I, i.e. no preconditioner
	norm_b = sqrt(maximum(eigvals( B' * B )))
end

# ╔═╡ af630960-a734-4124-81d3-a2657456f6c2
md"""
However, when using a diagonal preconditioner `P = Diagonal(A)`,
i.e. $\|\mathbf{B}\| = \|\mathbf{I} - \mathbf{P}^{-1} \mathbf{A}\|$,
then
"""

# ╔═╡ 51030114-04d7-49eb-9ac1-042941681446
let
	P = Diagonal(A)
	B = I - inv(P) * A  # Iteration matrix for P = Diagonal(A)
	norm_b = sqrt(maximum(eigvals( B' * B )))
end

# ╔═╡ 114bc28c-7888-4c6c-85bd-8c6232624491
md"""
Which is much smaller, in particular less than $1$.
Therefore only the preconditioned iterations converge:
"""

# ╔═╡ d8e0e584-985d-4ff0-8a5c-b21c83c226d9
let
	P = Diagonal(A)
	richardson_diagonal = richardson(A, b, P; tol=1e-10, maxiter=30)
	richardson_no_preconditioner = richardson(A, b, I; tol=1e-10, maxiter=30)
	
	plot(richardson_diagonal.history;
	     yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
	     label=L"$P = \textrm{Diagonal(}A\textrm{)}$", legend=:bottomright,
		 ylims=(1e-10, 1e5))
	plot!(richardson_no_preconditioner.history;
	      label=L"$P = I$", lw=2, mark=:o)
end

# ╔═╡ 96f2e900-d8e0-45d5-a79f-8df6d654fa3f
md"""
## Choosing a good preconditioner

From the above experiments we notice:

!!! info "Observation: Preconditioners"
    A good preconditioner $P$ in the Richardson iterations, satisfies the following properties:
    1. It is *cheap to invert*, that is linear systems $\mathbf{P} \mathbf{u} = \mathbf{r}$ are cheap to solve.
    2. The iteration matrix norm $\|\mathbf I - \mathbf P^{-1} \mathbf A\|$ is as small as possible,
      definitely smaller than one (to ensure convergence).

Clearly condition 2. suggests that the *perfect preconditioner* is $\mathbf{P} = \textbf{A}$, such that $\mathbf{P}^{-1} = \mathbf{A}^{-1}$. In this setting $\|\mathbf I - \mathbf{P}^{-1} \mathbf{A}\| = \|\mathbf I - \mathbf I\| = 0$,
i.e. we converge in a single Richardson step !
The trouble of this choice is that step 2 of Algorithm 1 (Richardson iteration) now effectively requires to solve the system $\mathbf A \mathbf u^{(k)} = \mathbf r^{(k)}$ for $\mathbf u^{(k)}$, i.e. we are back to solving the original problem. 

On the other hand condition 1. suggests to use $\mathbf P^{-1} = \mathbf I^{-1}$ (since the identity is cheapest to invert --- nothing needs to be done). However, this does not really do anything and does thus not reduce the value of $\|\mathbf{B}\|$. It will just be $\|\mathbf I - \mathbf A\|$.
"""

# ╔═╡ a7ebacfb-b163-40aa-99c8-f5873478249b
md"""
In practice we thus need to find a compromise between the two. As mentioned above standard choices for the preconditioner are:
- the diagonal of $\mathbf A$
- the lower triangle of $\mathbf A$
- another common choice is the *incomplete LU factorisation*, i.e. where one seeks a factorisation $\tilde{\mathbf L} \tilde{\mathbf U} \approx \mathbf A$ with a lower triangular matrix $\tilde{\mathbf L}$ and an upper triangular $\tilde{\mathbf U}$, which only *approximately* factorises $\mathbf A$. This gives additional freedom in designing the factorisations, in particular to avoid fill-in for sparse matrices.
- in many physics or engineering applications $\mathbf A$ results from a physical model. A preconditioner $\mathbf P$ can thus be resulting from an approximation to the employed physical model (e.g. by dropping the modelling of some physical effects).

We will not discuss the art of designing good preconditioners any further at this stage. Interested readers are referred to the excellent book [Youssef Saad *Iterative methods for Sparse Linear Systems*, SIAM (2003)](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf).
"""

# ╔═╡ d9bb93ce-5fac-4f3a-98b4-3024a23bcfb1
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
\mathbf{A} \mathbf{x}^{(k)} = \underbrace{\mathbf{b} - \mathbf{r}^{(k)}}_{=\widetilde{\mathbf{b}}}
```
where we defined the "perturbed" right-hand side $\widetilde{\mathbf{b}}$.
Notably $\mathbf{x}^{(k)}$ is thus *exact* solution of a linear system
involving $\mathbf{A}$ and this perturbed RHS,
while we actually care to find the true solution $\mathbf{x}_\ast$
obtained by solving the system $\mathbf{A} \mathbf{x}_\ast = \mathbf{b}$ employing an "unperturbed" right-hand side $\mathbf{b}$.

We are thus in exactly the same setting as our
final section on *Numerical stability* in our discussion
on [Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.html)
where instead of solving $\mathbf{A} \mathbf{x}_\ast = \mathbf{b}$
we are only able to solve the perturbed system
$\mathbf{A} \widetilde{\textbf{x}} = \widetilde{\mathbf{b}}$.
"""

# ╔═╡ 55a69e52-002f-40dc-8830-7fa16b7af081
md"""
We can thus directly apply Theorem 2
from [Direct methods for linear systems](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.html), which states that
```math
\frac{\|\mathbf{x}_\ast - \widetilde{\mathbf{x}} \|}{\| \mathbf{x}_\ast \|}
≤ κ(\mathbf{A}) 
\frac{ \left\| \mathbf{b} - \widetilde{\mathbf{b}} \right\| }{\| \mathbf{b} \|}
```
Keeping in mind that here $\widetilde{\mathbf{x}} = \mathbf{x}^{(k)}$
and $\mathbf{b} - \widetilde{\mathbf{b}} = \mathbf{r}^{(k)}$
we thus obtain
```math
\frac{\|\mathbf{x}_\ast - \mathbf{x}^{(k)} \|}{\| \mathbf{x}_\ast \|}
≤ κ(\mathbf{A}) 
\frac{ \left\| \mathbf{r}^{(k)} \right\| }{\| \mathbf{b} \|}
```
Combining this with our stopping criterion from **Algorithm 1**, that is
```math
\frac{\left\|\mathbf{r}^{(k)} \right\|}{\| \mathbf{b} \|} < ε,
```
we finally obtain
```math
\frac{\|\mathbf{x}_\ast - \mathbf{x}^{(k)} \|}{\| \mathbf{x}_\ast \|}
< κ(\mathbf{A}) \, ε.
```
In other words **our stopping criterion ensures** that the **relative error of the returned solution** is smaller than $κ(\mathbf{A})$ times the chosen tolerance.

If the matrix is well-conditioned, i.e. $κ(\mathbf{A})$ is close to $1$,
then the relative residual $\frac{\left\|\mathbf{r}^{(k)} \right\|}{\| \mathbf{b} \|}$
provides a good estimate of the true error and our stopping criterion is appropriate.
However, if $κ(\mathbf{A})$ is large, even a small residual may imply a large error in the returned solution.
"""

# ╔═╡ a5bedab1-b6a3-4539-82af-c18de2da4e92
Foldable("Alternative derivation without applying the previous Theorem 2",
md"""
We start by noting that $\mathbf{A} \mathbf{x}_\ast = \mathbf{b}$
and as discussed above that
$\mathbf{A} \mathbf{x}^{(k)} = \mathbf{b} - \mathbf{r}^{(k)} = \widetilde{\mathbf{b}}$. Now using the inequality
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
\stackrel{(6)}{≤} \frac{ \| \mathbf{A} \| \, \left\|\mathbf{A}^{-1} \left(\mathbf{b} - \widetilde{\mathbf{b}}\right) \right\|}{\| \mathbf{b} \|}
≤ \| \mathbf{A} \| \, \| \mathbf{A}^{-1} \| \, \frac{ \left\| \mathbf{b} - \widetilde{\mathbf{b}} \right\| }{\| \mathbf{b} \|}
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
With these quantities we rewrite to again obtain
```math
\frac{\|\mathbf{x}_\ast - \mathbf{x}^{(k)} \|}{\| \mathbf{x}_\ast \|}
≤ κ(\mathbf{A}) \frac{\left\|\mathbf{r}^{(k)} \right\|}{\| \mathbf{b} \|} = κ(\mathbf{A}) \, ε.
```
""")

# ╔═╡ 858b7c87-80de-4c1a-8598-03657608e323
md"""
## Jacobi and Gauss-Seidel method

*Note:* We will only discuss the high-level ideas of this part in the lecture. You can expect that there will not be any detailed exam questions on Jacobi and Gauss-Seidel  without providing you with the formulas and algorithms.


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

# ╔═╡ 41a269c8-fe64-4660-aa66-054b8af8af20
md"""
Importantly there is thus a **relation between optimisation problems** and **solving linear systems** if the system matrix $\textbf{A}$ is s.p.d.
"""

# ╔═╡ bf9a171a-8aa4-4f21-bde3-56ccef40de24
md"""
SPD matrices are not unusual. For example, recall that in polynomial regression problems (see least-squares problems in [Interpolation](https://teaching.matmat.org/numerical-analysis/05_Interpolation.html)),
where we wanted to find the best polynomial through the points
$(x_i, y_i)$ for $i=1, \ldots n$ by minimising the least-squares error,
we had to solve the *normal equations*
```math
\mathbf{V}^T \mathbf{V} \mathbf{c} = \mathbf{V}^T \mathbf{y}
```
where $\mathbf{c} \in \mathbf{R}^m$ are the unknown coefficients of the polynomial $\sum_{j=0}^m c_j x^j$, $\mathbf{y}$ is the vector
collecting all $y_i$ values and the $\mathbf{V}$ is the Vandermonde matrix
```math
\mathbf{V} = \left(\begin{array}{ccc}
1 & x_1 & \ldots & x_1^m \\
1 & x_2 & \ldots & x_2^m \\
\vdots \\
1 &x_n & \ldots & x_n^m \\
\end{array}\right) \in \mathbb{R}^{n\times (m+1)}.
```
In this case the system matrix $\mathbf{A} = \mathbf{V}^T \mathbf{V}$ is *always* s.p.d.
"""

# ╔═╡ 590e8877-0600-4bbb-8e88-7033548815b5
md"""
## Steepest descent method
"""

# ╔═╡ fff80ab4-d1df-4bef-a1d2-2aa08c78ce54
md"""
If we view the problem of solving linear systems as an optimisation problem,
a relatively simple idea is to construct an iterative method, where at every step $k$ we try to decrease the energy function $\phi$.

To guide our thoughts we consider a 2D problem which is easy to visualise. We take as an example
```math
\mathbf A^\text{2D} = \begin{pmatrix} 2& 0 \\ 0& 60 \end{pmatrix}
\qquad \mathbf b^\text{2D} = \begin{pmatrix}0\\0\end{pmatrix}
```
such that
```math
\phi(\mathbf x) = \frac12 \mathbf x^T \mathbf A^\text{2D} \mathbf x - \mathbf x^T \mathbf b^\text{2D} = x_1^2 + 30 x_2^2
```
"""

# ╔═╡ bb45f71b-b4cc-4b9a-abcb-7b0647f32a9b
A2d = [1.0  0.0;
	   0.0 20.0]

# ╔═╡ 53cc323c-d742-4cc6-8b96-e966e2b1ccb4
b2d = [0.0;
       0.0]

# ╔═╡ e22c0abe-caa8-4e65-a3c1-8b05e820b7e7
ϕ(x₁, x₂) = x₁^2 + 20 * x₂^2

# ╔═╡ 46cf4acb-5883-4cf0-bc15-3150f401a401
md"""
This problem is visualised as a contour plot below:
"""

# ╔═╡ 92035249-948e-4f2f-a00e-8b4ae0ec2742
begin
	function plot_contour(; show_grad=false)
		x = -6:0.1:6
		y = -2:0.1:2
		p = contour(x, y, ϕ; levels=5:5:100, clim=(0, 100), xlabel=L"x_1", ylabel=L"x_2", title="Contour plot of ϕ", aspect_ratio=1, legend=:topleft)

		scatter!(p, [5], [1]; mark=:o, label=L"x^{(k)}", markersize=5, c=1)

		if show_grad
			α = 1/10
			grad = A2d * [5; 1]
			plot!(p, [5, 5 - α*grad[1]], [1, 1 - α*grad[2]],
					 arrow=true, lw=3, c=3, label=L"$-\frac{1}{10} \, \nabla \phi(x^{(k)})$ "
					)
		end
		
		scatter!(p, [0], [0]; mark=:utriangle, label=L"Minimum $x_\ast$", markersize=6, c=1)
	end

	plot_contour(; show_grad=true)
end

# ╔═╡ 32b2a6e4-31ca-477f-9839-71cd940949b8
md"""
Now suppose we have an estimate $\mathbf{x}^{(k)}$ of the solutin of $\mathbf{A}\mathbf{x} = \mathbf{b}$ shown above by the blue dot.
This $\mathbf{x}^{(k)}$ is also an estimate of the minimum of $\phi$.
Our goal is thus to find a $\textbf{x}^{(k+1)}$ satisfying
$\phi(\textbf{x}^{(k+1)}) < \phi(\textbf{x}^{(k)})$,
i.e. to get closer to the minimum (blue triangle).

The gradient $\nabla \phi(\textbf{x}^{(k)})$ provides
the slope of the function $\phi$ at $\textbf{x}^{(k)}$
in the *upwards* direction.
Therefore, taking a step along the direction of $-\nabla \phi$
takes us downhill (see green arrow). We thus propose an algorithm
```math
\tag{11}
\begin{aligned}
\textbf{x}^{(k+1)} &= \textbf{x}^{(k)} - \alpha_k \nabla \phi(\textbf{x}^{(k)}) \\&= \textbf{x}^{(k)} + \alpha_k \textbf r^{(k)},
\end{aligned}
```
where we used that $\nabla \phi(\textbf{x}^{(k)}) = - \textbf r^{(k)}$,
and where $\alpha_k > 0$ is a parameter determining
how far to follow along the direction of the negative gradient.

Note, that this method is quite related to Richardson's iterations:
for $α_k = 1$ we actually recover **Algorithm 1** with $\mathbf{P} = \mathbf{I}$.
"""

# ╔═╡ f8f62fa0-6a68-4083-909f-4a3367157f12
Foldable("Rationale for introducing αₖ. Or: why not just take αₖ = 1 ?",
md"""
The rationale of having this parameter $\alpha_k$
is that taking too large steps from $\mathbf x^{(k)}$
to $\mathbf x^{(k+1)}$ is in general
not useful. For example in the example visualised above the green arrow
represents only $1/40$ times the actual gradient.
Therefore simply
updating $\mathbf x^{(k+1)} = \mathbf x^{(k)} + \mathbf r^{(k)}$
with $\alpha_k = 1$ will substantially overshoot,
such that $\mathbf x^{(k+1)}$ will be even further from the minimum
than $\mathbf x^{(k)}$.
Therefore we employ $0 < \alpha_k < 1$ to down-scale this step.
""")

# ╔═╡ cf1af660-1326-4c57-900f-ea12e9836f85
md"""
Since overall our goal is to find the minimum of $\phi$,
a natural idea to determine $\alpha_k$ exactly such that
we make the value of $\phi(\textbf{x}^{(k+1)})$ as small as possible.
We compute

```math
\begin{aligned}
ϕ(\mathbf{x}^{(k+1)}) &= ϕ(\mathbf{x}^{(k)} + α_k \mathbf{r}^{(k)}) \\
&= \frac12 \left(\mathbf{x}^{(k)} + α_k \mathbf{r}^{(k)}\right)^T \mathbf{A}
\left(\mathbf{x}^{(k)} + α_k \mathbf{r}^{(k)}\right)
- \left(\mathbf{x}^{(k)} + α_k \mathbf{r}^{(k)}\right)^T\mathbf{b}\\
&= \underbrace{\frac12 (\mathbf{r}^{(k)})^T\mathbf A\mathbf{r}^{(k)}}_{=a} α_k^2
+ \underbrace{\left[
	\frac12 (\mathbf{x}^{(k)})^T\mathbf A\mathbf{r}^{(k)}
	+ \frac12 (\mathbf{r}^{(k)})^T\mathbf A\mathbf{x}^{(k)}
	- (\mathbf{r}^{(k)})^T \mathbf{b}
\right]}_{=b} α_k^{(k)} \\
	&\hspace{8.8em} + \underbrace{\left[
	\frac12 (\mathbf{x}^{(k)})^T\mathbf A\mathbf{x}^{(k)}
	- (\mathbf{x}^{(k)})^T \mathbf{b}
\right]}_{=c}
\end{aligned}
```
where
```math
\begin{aligned}
a &= \frac12 (\mathbf{r}^{(k)})^T\mathbf A\mathbf{r}^{(k)}&
c &= 
	\frac12 (\mathbf{x}^{(k)})^T\mathbf A\mathbf{x}^{(k)}
	- (\mathbf{x}^{(k)})^T \mathbf{b}
\end{aligned}
```
and
```math
\begin{aligned}
b &= 
	\frac12 (\mathbf{x}^{(k)})^T\mathbf A\mathbf{r}^{(k)}
	+ \frac12 (\mathbf{r}^{(k)})^T\mathbf A\mathbf{x}^{(k)}
	- (\mathbf{r}^{(k)})^T \mathbf{b}\\
	&= (\mathbf{x}^{(k)})^T\mathbf A\mathbf{r}^{(k)}- (\mathbf{r}^{(k)})^T \mathbf{b}\\
	&= -(\mathbf{r}^{(k)})^T \, \mathbf{r}^{(k)}
\end{aligned}
```
"""

# ╔═╡ 599d5e8c-8e36-401e-a576-1d47ef045c18
md"""
We notice that $ϕ(\mathbf{x}^{(k+1)}) = a α_k^2 + bα_k + c$ is a second-degree polynomial in $α_k$. Setting its gradient to zero gives us the $α_k$
which minimises the $ϕ$ as much as possible in this step. From
```math
0 = \frac{dϕ}{dα_k}(\mathbf{x}^{(k+1)}) = 2aα_k + b
```
we obtain the optimal step size as
```math
\tag{12}
α_k = \frac{-b}{2a} = \frac{(\mathbf{r}^{(k)})^T \, \mathbf{r}^{(k)}}{(\mathbf{r}^{(k)})^T\mathbf A\,\mathbf{r}^{(k)}}.
```
In summary we obtain the algorithm:
"""

# ╔═╡ 6d4b968d-9cdc-4ca0-971a-6f2e5a8d06d8
md"""
!!! info "Algorithm 2: Steepest descent method"
    Let $\mathbf{A} \in \mathbb{R}^{n\times n}$ be an s.p.d. system matrix,
    right-hand side $\mathbf{b} \in \mathbb{R}^n$, 
	an initial guess $\mathbf{x}^{(0)} \in \mathbb{R}^n$
	and convergence threshold $ε$.

	Compute the inital residual $\mathbf{r}^{(0)} = \mathbf{b} - \mathbf{A} \mathbf{x}^{(0)}$.
	Then for $k = 0, 1, \ldots$ iterate:
	1. Compute $\mathbf{w}^{(k)} = \mathbf{A} \mathbf{r}^{(k)}$
	2. Step size: $α_k = \frac{(\mathbf{r}^{(k)})^T \, \mathbf{r}^{(k)}}{(\mathbf{r}^{(k)})^T\textcolor{red}{\mathbf{w}^{(k)}}}$ $\qquad$ *(Optimal because of (12) )*
	3. Take step: $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + α_k \mathbf{r}^{(k)}$.
	4. Update residual $\mathbf{r}^{(k+1)} = \mathbf{r}^{(k)} - α_k \mathbf{w}^{(k)}$
    The iteration is stopped as soon as
    $\frac{\| \mathbf{r}^{(k)} \|}{\| \mathbf{b} \|} < ε$.
"""

# ╔═╡ 3e4785cc-4120-4f69-b6ad-27da5bb3a17e
md"""
Notice that in this algorithm we used the trick
```math
\mathbf{r}^{(k+1)} = \mathbf{b} - \mathbf{A}\mathbf{x}^{(k+1)}
= \mathbf{b} - \mathbf{A}\left(\mathbf{x}^{(k)} + α_k \mathbf{r}^{(k)} \right)
= \mathbf{r}^{(k)} - α_k \mathbf{w}^{(k)}
```
to compute the residual $\mathbf{r}^{(k+1)}$ from $\mathbf{r}^{(k)}$ and $\mathbf{w}^{(k)}$, that is without computing another matrix-vector product. As matrix-vector products scale as $O(n^2)$, but vector operations only as $O(n)$ this saves computational cost.

An implementation of this algorithm is given by the following function:
"""

# ╔═╡ 0159a195-627c-4e07-834e-222248665da8
function steepest_descent_simple(A, b; x=zero(b), tol=1e-6, maxiter=100)
	resnorms = Float64[]  # Track relative residual norm
	history = [x]
	
	r = b - A * x
	for k in 1:maxiter
		relnorm = norm(r) / norm(b)
		push!(resnorms, relnorm)
		if relnorm < tol
			break
		end
		
		Ar = A * r
		α = dot(r, r) / dot(r, Ar)
		x = x + α * r
		r = r - α * Ar

		push!(history, x)
	end
	(; x, history, resnorms)
end

# ╔═╡ c330ec39-bcce-4b97-b5c9-11c13967de56
md"""
Running it on our 2x2 example problem plotted above it produces the following steps
(first 6 steps shown). 

We realise that convergence is steadily towards the minimum, but seems to oscillate
around the "best possible path". We will see in a second why this is.
"""

# ╔═╡ 0b7aad76-2e9d-41de-851d-5075d9f840b3
let
	p = plot_contour()
	res = steepest_descent_simple(A2d, b2d; x=[5.0, 1.0], maxiter=10)
	for i in 1:10
		label = i == 1 ? "Steepest Decent iterations" : ""
		plot!(p, first.(res.history)[i:i+1], last.(res.history)[i:i+1]; c=:red, arrow=i < 6, label, lw=2)
	end
	p
end

# ╔═╡ 3c09a9c8-015b-4d17-9643-1f270b8a1d9b
md"""
In line with our previous discussion we can again view **Algorithm 2**
as a fixed-point method $\mathbf{x}^{(k+1)} = g(\mathbf{x}^{(k)})$,
this time with fixed-point function
```math
g(\mathbf{x}) = \mathbf{x} + α(\mathbf{x}) \, \left( \mathbf{b} - \mathbf{A} \mathbf{x} \right)
```
where $α(\mathbf{x})$ is meant to indicate the step size determined according to equation (12). Clearly a fixed-point of this function satisfies
$\mathbf{x} =  \mathbf{x} + α(\mathbf{x}) \, \left( \mathbf{b} - \mathbf{A} \mathbf{x} \right)$
and thus $\mathbf{b} = \mathbf{A} \mathbf{x}$.
"""

# ╔═╡ 7207b521-e642-4ffd-a590-74b0f474f7a8
md"""
While the details are more involved and beyond our scope,
applying the usual convergence theory on fixed-point methods
yields the following observation:

!!! note "Theorem 4: Convergence of steepest descent"
    Given an s.p.d. matrix $\mathbf{A} \in \mathbb{R}^{n \times n}$
	and a right-hand side $\mathbf{b}$ the solution of the linear system
	$\mathbf{A} \mathbf{x} = \mathbf{b}$ can be obtained by applying
	steepest descent starting from *any* initial position $\mathbf{x}^{(0)}$.
	The following error estimate holds:
	```math
	\tag{13}
	\|\mathbf{x} - \mathbf{x}^{(k)}\| \leq C \left(\frac{κ(A) - 1}{κ(A) + 1} \right)^k \|\mathbf{x} - \mathbf{x}^{(0)}\|
	```
	where $C$ is a positive constant.
	This is **linear convergence** with rate $\frac{κ(A) - 1}{κ(A) + 1}$.

The main result is thus that steepest descent **always converges**.
However, if $κ(A)$ is large the convergence can be slow.
Here is a plot of the rate with increasing $κ$.
Recall, that smaller rate means faster convergence
and a rate of $1$ is essentially stagnation.
"""

# ╔═╡ 7cc7488b-fe17-496c-8c59-c074819d81f0
plot(κ -> (κ-1)/(κ+1), xaxis=:log, xlims=(1, 10^3), label="", xlabel=L"κ(A)", ylabel="Rate", lw=3)

# ╔═╡ 95b48015-48d4-4742-94dd-fac1dbe5f041
md"We try steepest descent on the following random SPD matrix."

# ╔═╡ 9e34bd04-a4b7-4b60-b93b-dcd37b25287b
begin
	# Generate a random s.p.d. matrix
	U, Σ, V = svd(randn(100, 100))
	Aₛ = U * Diagonal(abs.(Σ)) * U'
end

# ╔═╡ 45cfdf30-8738-4cb6-a49a-534392128475
bₛ = ones(100);

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
	
	plot(descent.resnorms;
		yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||", ylims=(5e-2, 1.3),
		title="Steepest descent convergence", label="Observed convergence", c=2)
	
	plot!(1:40, x -> 0.15rate_steep_desc^x, lw=2, ls=:dash, label="Theoretical rate", c=2)
end

# ╔═╡ cc4b5c31-cd6a-4df8-b1cf-068e209aaa33
md"""
## Optional: Preconditioned steepest descent

We noticed above that a condition number $κ(\mathbf{A}) \approx 1$ is providing the highest convergence rate for steepest descent.
For cases where the condition number is large, a potential cure is to
apply *preconditioning*. The idea is similar to Richardson iteration.
Assume we have a matrix $\mathbf{P} \approx \mathbf{A}$.
Then instead of applying steepest descent to $\mathbf{A} \mathbf{x} = \mathbf{b}$
we instead apply it to
```math
\mathbf{P}^{-1} \mathbf{A} \mathbf{x} = \mathbf{P}^{-1} \mathbf{b},
```
such that this new method will now converge as
```math
\|\mathbf{x} - \mathbf{x}^{(k)}\| \leq C \left(\frac{κ(\mathbf{P}^{-1}\mathbf A) - 1}{κ(\mathbf{P}^{-1} \mathbf A) + 1} \right)^k \|\mathbf{x} - \mathbf{x}^{(0)}\|
```
If we use a good preconditioner, then $κ(\mathbf{P}^{-1} \mathbf A) \ll κ(\mathbf A)$ and convergence accelerates.
"""

# ╔═╡ 1c819b4a-d6d4-408f-a30f-aa4f6bb2871e
md"""
Without going into many details, we remark that there are a few subtleties that apply:
- Even if $\mathbf{P}$ *and* $\mathbf{A}$ are s.p.d. matrices, the matrix $\mathbf{P}^{-1} \mathbf A$ *may not be s.p.d.*. Therefore a practical implementation of preconditioned gradient descent is not realised by blindly applying the plain algorithm to $\mathbf{P}^{-1} \mathbf{A}$, but instead by "inlining" the application of $\mathbf{P}^{-1}$ into the actual algorithm. At the $k$-th iteration we thus want to update the solution $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + α_k \mathbf{z}^{(k)}$ using the *preconditioned residual*
  ```math
  \mathbf{z}^{(k)} = \mathbf{P}^{-1}\mathbf{b} - \mathbf{P}^{-1} \mathbf{A} \mathbf{x}^{(k)} = \mathbf{P}^{-1} \mathbf{r}^{(k)}.
  ```
  Following through with this change also in the computation of the optimal
  step size leads to $α_k = \frac{(\mathbf{z}^{(k)})^T \, \mathbf{r}^{(k)}}{(\mathbf{z}^{(k)})^T\mathbf{A} \mathbf{z}^{(k)}}$, thus the following algorithm:
  !!! info "Algorithm 3: Preconditioned steepest descent method"
      Let $\mathbf{A} \in \mathbb{R}^{n\times n}$ be an s.p.d. system matrix,
      right-hand side $\mathbf{b} \in \mathbb{R}^n$, 
      preconditioner  $\mathbf{P} \in \mathbb{R}^{n\times n}$,
      an initial guess $\mathbf{x}^{(0)} \in \mathbb{R}^n$
      and convergence threshold $ε$.

      Compute the inital residual $\mathbf{r}^{(0)} = \mathbf{b} - \mathbf{A} \mathbf{x}^{(0)}$.
      Then for $k = 0, 1, \ldots$ iterate:
      1. Solve $\mathbf{P} \mathbf{z}^{(k)} = \mathbf{r}^{(k)}$ for $\mathbf{z}^{(k)}$ $\qquad$ *(This is the new step)*
      2. Compute $\mathbf{w}^{(k)} = \mathbf{A} \mathbf{z}^{(k)}$  $\qquad$ *(This has been changed)*
      3. Step size: $α_k = \frac{(\mathbf{z}^{(k)})^T \, \mathbf{r}^{(k)}}{(\mathbf{z}^{(k)})^T\mathbf{w}^{(k)}}$
      4. Take step: $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + α_k \mathbf{z}^{(k)}$.
      5. Update residual $\mathbf{r}^{(k+1)} = \mathbf{r}^{(k)} - α_k \mathbf{w}^{(k)}$
      The iteration is stopped as soon as
      $\frac{\| \mathbf{r}^{(k)} \|}{\| \mathbf{b} \|} < ε$.
- Similar to our discussion in the context of Richardson iterations, the choice of the preconditioner is delicate. The "perfect" preconditioner $\mathbf{P} = \mathbf{A}$ brings our condition number down to $κ(\mathbf A^{-1} \mathbf A) = κ(\mathbf I) = 1$, but requires us to solve the full problem in the first step of **Algorithm 3**. On the other hand no preconditioning, i.e. $\mathbf{P} = \mathbf{I}$, gives a cheap computation of $\mathbf{z}^{(k)}$ from $\mathbf{P} \mathbf{z}^{(k)} = \mathbf{r}^{(k)}$, but does nothing to reduce the conditioning. Often a simple preconditioner, such as the diagonal of the matrix already provides noteworthy improvements. For example:
"""

# ╔═╡ fdbe3ba7-d1a3-43b0-bd56-ad0bc6f2200d
function steepest_descent(A, b, P; x=zero(b), tol=1e-6, maxiter=100)
	resnorms = Float64[]
	
	r = b - A * x
	for k in 1:maxiter
		relnorm = norm(r) / norm(b)
		push!(resnorms, relnorm)
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
	(; x, resnorms)
end

# ╔═╡ bf028cd7-c911-4329-86b0-442061a8e8c8
let
    P = Diagonal(Aₛ)
    descent_diagonal = steepest_descent(Aₛ, bₛ, P; tol=1e-6, maxiter=30)
    descent_no_preconditioner = steepest_descent_simple(Aₛ, bₛ; tol=1e-6, maxiter=30)
       
    plot(descent_diagonal.resnorms;
         yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
         label=L"$P = \textrm{Diagonal(}A\textrm{)}$", legend=:bottomright,
         title="Steepest descent convergence", ylims=(1e-4, 1.3))
    plot!(descent_no_preconditioner.resnorms;
          label=L"$P = I$", lw=2, mark=:o)
end

# ╔═╡ d698fe83-68a6-4988-b1e0-9e2dace2c6fd
md"""
## Optional: Conjugate gradient method

Recall that the key idea of steepest descent was to employ the update
```math
\textbf{x}^{(k+1)} = \textbf{x}^{(k)} - α_k \nabla \phi(\textbf{x}^{(k)})
=  \textbf{x}^{(k)} + α_k \mathbf{r}^{(k)}
```
that is to follow in the *direction of maximal descent*
--- i.e. along the direction of the negative gradient of $\phi$.

While this is certainly a natural choice, we also saw that this can lead to a rather unsteady convergence behaviour:
"""

# ╔═╡ 653dd5fd-2771-492e-9508-2b473d448be4
let
	p = plot_contour()
	res = steepest_descent_simple(A2d, b2d; x=[5.0, 1.0], maxiter=10)
	for i in 1:10
		label = i == 1 ? "Steepest Decent iterations" : ""
		plot!(p, first.(res.history)[i:i+1], last.(res.history)[i:i+1]; c=:red, arrow=i < 6, label, lw=2)
	end
	p
end

# ╔═╡ 6671fe27-7b0f-458a-a012-747c5fcb76eb
md"""
One way to cure this behaviour is in fact make a rather different choice for the update direction. While for the first update we keep $\mathbf{p}^{(1)} = \mathbf{r}^{(k)}$ we subsequently choose directions $\mathbf{p}^{(k)}$ with the property
```math
0 = \mathbf{p}^{(k)} \mathbf{A} \mathbf{p}^{(j)} \qquad \text{with $j = 0, 1, \ldots, k-1$}.
```
This property of the vectors $\mathbf{p}^{(k)}$ is usually called $\mathbf{A}$-orthogonality. Based on this update direction we iterate as
```math
\textbf{x}^{(k+1)} = \textbf{x}^{(k)} + α_k \mathbf{p}^{(k)}
```
with the optimal step size now being given by
```math
α_k = \frac{(\mathbf{p}^{(k)})^T \, \mathbf{r}^{(k)}}{(\mathbf{p}^{(k)})^T\mathbf{A} \mathbf{p}^{(k)}}.
```
The next $\mathbf{p}^{(k+1)}$ is found by $A$-orthogonalising $\mathbf{r}^{(k+1)}$ against all previous $\mathbf{p}^{j}$ with $j = 1, \ldots, k$. Perhaps surprisingly this can be achieved by the following recurrent algorithm
```math
\mathbf{p}^{(k+1)} = \mathbf{r}^{(k+1)} - β_k\, \mathbf{p}^{k}
\qquad \text{with} \quad β_k = \frac{(\mathbf{A} \mathbf{p}^{(k)})^T \, \mathbf{r}^{(k)}}{(\mathbf{p}^{(k)})^T\mathbf{A} \mathbf{p}^{(k)}}.
```

This is called the [Conjugate gradient method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) and despite its invention in 1952 is still the state-of-the-art method for solving linear systems or optimisation problems involving s.p.d. matrices.
"""

# ╔═╡ 1d3ea21f-0ebd-491b-a165-2e5cf7cc583f
md"""
Studying its convergence leads to the following strong result:

!!! note "Theorem 4: Convergence of Conjugate Gradient (CG)"
    Given an s.p.d. matrix $\mathbf{A} \in \mathbb{R}^{n \times n}$
	and a right-hand side $\mathbf{b}$ the solution of the linear system
	$\mathbf{A} \mathbf{x} = \mathbf{b}$ can be obtained by applying
	the conjugate gradient algorithm
	starting from *any* initial position $\mathbf{x}^{(0)}$
	**in at most $n$ steps** in exact arithmetic.
	Moreover the following error estimate holds:
	```math
	\tag{14}
	\|\mathbf{x} - \mathbf{x}^{(k)}\| \leq C \left(\frac{\sqrt{κ(A)} - 1}{\sqrt{κ(A)} + 1} \right)^k \|\mathbf{x} - \mathbf{x}^{(0)}\|
	```
	where $C$ is a positive constant.
	This is **linear convergence** with rate $\frac{\sqrt{κ(A)} - 1}{\sqrt{κ(A)} + 1}$.

Notice (1) that this method is **guaranteed to converge after $n$ steps** and (2)  that the rate deteriorates much slower as the condition number increases (recall that smaller rates are better):
"""

# ╔═╡ 9061df24-9a36-4612-9b6f-926d8641231a
let
	p = plot(κ -> (κ-1)/(κ+1), xaxis=:log, xlims=(1, 10^3), label="Steepest descent", xlabel=L"κ(A)", ylabel="Rate", lw=3)
	plot!(p, κ -> (sqrt(κ)-1)/(sqrt(κ)+1), lw=3, label="Conjugate gradient", ls=:dash)
end

# ╔═╡ 9e3340b7-1ecb-4f6d-96a0-089932ec306d
md"""
Indeed for our 2x2 matrix example two steps of CG are enough:
"""

# ╔═╡ 96d9ceb9-2917-44bd-a7af-724d7dd50414
md"""
Furthermore we note our random matrix example to converge noticably faster:
"""

# ╔═╡ 7cbc6819-875a-498f-a81a-f63652c3dc8c
rate_cg = (sqrt(cond(Aₛ)) - 1) / (sqrt(cond(Aₛ)) + 1)

# ╔═╡ cf8ae086-da8a-4fd7-83cd-167b8b8d7dd4
md"""
A (non-optimised) CG implementation is:
"""

# ╔═╡ 5291c5f3-7d79-42b4-b704-24697fcc9782
function conjugate_gradient_simple(A, b; x=zero(b), tol=1e-6, maxiter=100)
	resnorms = Float64[]  # Track relative residual norm
	history  = [x]
	
	r = b - A * x
	p = r
	for k in 1:maxiter
		relnorm = norm(r) / norm(b)
		push!(resnorms, relnorm)
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

		push!(history, x)
	end
	(; x, resnorms, history)
end

# ╔═╡ 2ad6d547-88d0-4ad2-b9d5-09719212b0cc
let
	p = plot_contour()

	res = conjugate_gradient_simple(A2d, zeros(2); x=[5.0, 1.0], maxiter=4)
	for i in 1:2
		label = i == 1 ? "Conjugate Gradient" : ""
		plot!(p, first.(res.history)[i:i+1], last.(res.history)[i:i+1]; c=:blue, arrow=true, label, lw=4)
	end
	
	res = steepest_descent_simple(A2d, zeros(2); x=[5.0, 1.0], maxiter=10)
	for i in 1:10
		label = i == 1 ? "Steepest Decent" : ""
		plot!(p, first.(res.history)[i:i+1], last.(res.history)[i:i+1]; c=:red, arrow=i < 6, label, lw=2)
	end
	
	p
end

# ╔═╡ d0dfcbe6-6115-4f87-a310-5db427d90346
let
	P = Diagonal(A)
	descent_no_preconditioner = steepest_descent_simple(Aₛ, bₛ; tol=1e-6, maxiter=40)
	cg_no_preconditioner = conjugate_gradient_simple(Aₛ, bₛ; tol=1e-6, maxiter=40)
	
	plot(descent_no_preconditioner.resnorms;
	     yaxis=:log, mark=:o, lw=2, ylabel=L"||r|| / ||b||",
	     label="Steepest descent", legend=:bottomleft, c=2,  ylims=(1e-4, 1.3))
	plot!(cg_no_preconditioner.resnorms; label="CG", lw=2, mark=:o, c=1)


	plot!(1:40, x -> 0.15rate_steep_desc^x, lw=2, ls=:dash, label="Rate steep. desc.", c=2)
	plot!(1:40, x -> 0.5rate_cg^x, lw=2, ls=:dash, label="Rate CG", c=1)
end

# ╔═╡ e566f53e-e841-40d4-914d-de7f660a1e67
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 420)
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
LaTeXStrings = "~1.4.0"
Plots = "~1.40.0"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "70b8b0e1f666a1a5a6c73ee25fe9849578771333"

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
git-tree-sha1 = "a434e811d10e7cbf4f0674285542e697dca605d0"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.42"

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
git-tree-sha1 = "24be21541580495368c35a6ccef1454e7b5015be"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.11"

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
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

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
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

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
# ╟─63bb7fe9-750f-4d2f-9d18-8374b113373e
# ╠═0501f270-c39c-11ee-2077-7f7f409b6c66
# ╟─7d9c9392-3aec-4efd-a9ba-d8965687b163
# ╟─4e0e5311-3428-44f5-a616-68cf7634214f
# ╟─b5bf49df-f6ff-49f9-b016-7f0159b4af7d
# ╟─1292abe9-5d88-45b0-8904-b89a03009968
# ╟─c6617ee2-be59-49b2-b0bc-7423b49c2fce
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
# ╟─03984b5d-3fe9-4557-9a0b-ddcc7d0a6c34
# ╟─5913dace-5e27-4d4d-8bc3-0b4118b1b093
# ╟─7b8bb857-35d3-4014-851b-9ab263f48967
# ╟─1f7765a2-4986-4ef2-849e-ca1efb59abcf
# ╟─ba77a115-f9fb-4f80-a0ea-eb50ddd9ae65
# ╟─6f3d76f5-aeec-4a05-8aea-a4f97d625a4a
# ╟─28d6cb69-71ab-493f-a9a1-60f2fe37c1ed
# ╟─6a5346df-6c21-47fb-9f99-d0bc75532135
# ╟─f2d4fbe9-e19a-41e0-9f3b-2f77ef0b12a7
# ╟─ef1ecb29-3ca0-4505-b96f-a56abd961979
# ╠═e277ff4c-87a8-440d-bba4-f30a20b47788
# ╠═54cb29d7-4933-4ab4-a3b5-26431d8432f9
# ╟─af630960-a734-4124-81d3-a2657456f6c2
# ╠═51030114-04d7-49eb-9ac1-042941681446
# ╟─114bc28c-7888-4c6c-85bd-8c6232624491
# ╠═d8e0e584-985d-4ff0-8a5c-b21c83c226d9
# ╟─96f2e900-d8e0-45d5-a79f-8df6d654fa3f
# ╟─a7ebacfb-b163-40aa-99c8-f5873478249b
# ╟─d9bb93ce-5fac-4f3a-98b4-3024a23bcfb1
# ╟─55a69e52-002f-40dc-8830-7fa16b7af081
# ╟─a5bedab1-b6a3-4539-82af-c18de2da4e92
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
# ╟─41a269c8-fe64-4660-aa66-054b8af8af20
# ╟─bf9a171a-8aa4-4f21-bde3-56ccef40de24
# ╟─590e8877-0600-4bbb-8e88-7033548815b5
# ╟─fff80ab4-d1df-4bef-a1d2-2aa08c78ce54
# ╠═bb45f71b-b4cc-4b9a-abcb-7b0647f32a9b
# ╠═53cc323c-d742-4cc6-8b96-e966e2b1ccb4
# ╠═e22c0abe-caa8-4e65-a3c1-8b05e820b7e7
# ╟─46cf4acb-5883-4cf0-bc15-3150f401a401
# ╟─92035249-948e-4f2f-a00e-8b4ae0ec2742
# ╟─32b2a6e4-31ca-477f-9839-71cd940949b8
# ╟─f8f62fa0-6a68-4083-909f-4a3367157f12
# ╟─cf1af660-1326-4c57-900f-ea12e9836f85
# ╟─599d5e8c-8e36-401e-a576-1d47ef045c18
# ╟─6d4b968d-9cdc-4ca0-971a-6f2e5a8d06d8
# ╟─3e4785cc-4120-4f69-b6ad-27da5bb3a17e
# ╠═0159a195-627c-4e07-834e-222248665da8
# ╟─c330ec39-bcce-4b97-b5c9-11c13967de56
# ╟─0b7aad76-2e9d-41de-851d-5075d9f840b3
# ╟─3c09a9c8-015b-4d17-9643-1f270b8a1d9b
# ╟─7207b521-e642-4ffd-a590-74b0f474f7a8
# ╟─7cc7488b-fe17-496c-8c59-c074819d81f0
# ╟─95b48015-48d4-4742-94dd-fac1dbe5f041
# ╠═9e34bd04-a4b7-4b60-b93b-dcd37b25287b
# ╠═45cfdf30-8738-4cb6-a49a-534392128475
# ╟─646ebf55-d081-4884-8524-15e5b0ecebcc
# ╠═d7eedfae-90ae-4917-b629-fc56f00bc3da
# ╟─7481b32b-d8c4-405b-9a97-dfb17ce87f53
# ╠═1a675fdb-885d-462e-bc77-9ea074d49614
# ╠═f0c0071f-5499-4bfa-b838-4f3d72ae3b3d
# ╟─cc4b5c31-cd6a-4df8-b1cf-068e209aaa33
# ╟─1c819b4a-d6d4-408f-a30f-aa4f6bb2871e
# ╟─fdbe3ba7-d1a3-43b0-bd56-ad0bc6f2200d
# ╟─bf028cd7-c911-4329-86b0-442061a8e8c8
# ╟─d698fe83-68a6-4988-b1e0-9e2dace2c6fd
# ╟─653dd5fd-2771-492e-9508-2b473d448be4
# ╟─6671fe27-7b0f-458a-a012-747c5fcb76eb
# ╟─1d3ea21f-0ebd-491b-a165-2e5cf7cc583f
# ╠═9061df24-9a36-4612-9b6f-926d8641231a
# ╟─9e3340b7-1ecb-4f6d-96a0-089932ec306d
# ╟─2ad6d547-88d0-4ad2-b9d5-09719212b0cc
# ╟─96d9ceb9-2917-44bd-a7af-724d7dd50414
# ╠═7cbc6819-875a-498f-a81a-f63652c3dc8c
# ╟─cf8ae086-da8a-4fd7-83cd-167b8b8d7dd4
# ╟─d0dfcbe6-6115-4f87-a310-5db427d90346
# ╟─5291c5f3-7d79-42b4-b704-24697fcc9782
# ╟─e566f53e-e841-40d4-914d-de7f660a1e67
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
