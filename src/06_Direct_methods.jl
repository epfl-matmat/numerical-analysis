### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3295f30c-c1f4-11ee-3901-4fb291e0e4cb
begin
	using LinearAlgebra
	using SparseArrays
	using PlutoUI
	using PlutoTeachingTools
	using HypertextLiteral
end

# ╔═╡ ca2c949f-a6a0-485f-bd52-5dae3b050612
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/06_Direct_methods.pdf)
"""

# ╔═╡ 21c9a859-f976-4a93-bae4-616122712a24
TableOfContents()

# ╔═╡ 8e4312b5-f4ce-440b-a487-689959851098
md"""
# Direct methods for linear systems

In the previous chapter on polynomial interpolation we were already
confronted with the need to solve linear systems, that is a system
of equations of the form
```math
\tag{1}
\mathbf{A} \mathbf{x} = \mathbf{b},
```
where we are given a matrix $\mathbf{A} \in \mathbb{R}^{n\times n}$
as well as a right-hand side $\mathbf{b} \in \mathbb{R}^n$.
As the solution we seek the unknown $\mathbf{x} \in \mathbb{R}^n$.

Solving such equations is a standard exercise in linear algbra
and you probably have already done it in previous courses.
However, you might also know that solving such linear systems is
not always equally easy.

This one, for example, is particularly difficult:
"""

# ╔═╡ fafa019e-8a47-4a80-bfd8-a017ae31ab1f
M_difficult = [π          ℯ;
			   355/113    23225/8544]   # Rational approximations to π and e

# ╔═╡ e45952b2-714b-427b-8979-98d11c330294
b_difficult = [1.0, 1.0]

# ╔═╡ 52d8634e-116b-48e6-9673-2ee2a97423c0
md"""
Solving it in 64-bit `Float64` gives this result:
"""

# ╔═╡ 5fb995cc-b338-4608-a01c-fc0e84d9dfe9
M_difficult \ b_difficult

# ╔═╡ 4d9d19da-7a4f-49ba-9c6f-890dad00672c
md"""
while using 32-bit `Float32` gives this:
"""

# ╔═╡ 08ec21d9-5474-4aae-8fb1-62fde1c7e864
let
	M_32 = convert(Matrix{Float32}, M_difficult)
	b_32 = convert(Vector{Float32}, b_difficult)
	M_32 \ b_32
end

# ╔═╡ a8559c5e-744a-4c25-a75c-5394558bc672
md"""
These clearly differ and we would probably trust the `Float64` result more.
But is this justified ? 
And what are conditions to quantify that a linear system is more
challenging than another ?

This is what we want to explore in this part of the course.
"""

# ╔═╡ 0b8bccd7-28ad-4ebc-8bcc-71e6adda71e3
md"""
## What Julia's `\` does under the hood

For solving linear systems, we so far employed Julia's backslash `\` operator.
So let's find out what it actually does under the hood. We take the problem
"""

# ╔═╡ 9230561b-26f3-4aa0-9c54-2f15aa231d80
begin
	A = Float64[-4  3 -1;
	             2  1  0;         
	             4 -3  4]
	b = [2, 4, -2]
end;

# ╔═╡ c74de7e2-b1b0-4f7f-9ac8-31d216caf987
md"as an example. Normally we would now just do `A \ b`, so let's check what this calls"

# ╔═╡ 128af99c-4c34-4845-9fe8-876f5a003412
methods(\, (Matrix, Vector))

# ╔═╡ 0d040c4e-f306-4a84-93ac-9c35bc924130
md"""
Ok, so this calls into Julia's linear algebra library. The code is [located here](https://github.com/JuliaLang/julia/blob/f2e168a0a8bbd6be82b4a0dfa205bc8080a48824/stdlib/LinearAlgebra/src/generic.jl#L1132). Essentially it performs
an **LU factorisation**:
"""

# ╔═╡ fd7a0025-1b77-45ab-bcbb-df39e38b3ffd
LU = lu(A)

# ╔═╡ 5471b1c1-224c-4464-a4a3-7719347048fc
md"""
which as we can see *factorises* the matrix $\mathbf A$ into a **lower triangular** matrix $\mathbf L$ and an **upper triangular** matrix $\mathbf U$. More formally:

!!! note "Definition: LU factorisation"
    Given a matrix $\mathbf A \in \mathbb{R}^{n\times n}$ find a lower triangular matrix
    $\mathbf L \in \mathbb{R}^{n\times n}$ and an upper triangular matrix
    $\mathbf U \in \mathbb{R}^{n\times n}$,
    such that
	```math
    \tag{2}
	\mathbf{A} = \mathbf{L} \mathbf{U}
	```

We can easily check that by multiplying $\mathbf L$ and  $\mathbf U$
we indeed revover $\mathbf A$:
"""

# ╔═╡ 2a1a8730-795b-4659-b688-6e320ffa81f8
LU.L * LU.U

# ╔═╡ a067a699-5cfe-405e-bae2-4e1033f9f9c9
A - LU.L * LU.U

# ╔═╡ ececcd18-56a3-41f7-8b15-bdff974080ac
md"""
Now we are a step closer to what happens, but why is this useful ?
"""

# ╔═╡ 06716c35-b551-40a6-bb14-906bc7569434
md"""
## Solving triangular systems

The factorisation into two triangular matrices is useful,
because it is especially easy to solve a system where
the matrix is triangular. For example consider the *lower* triangular system

```math
  \begin{pmatrix}
    4 & 0 & 0 & 0 \\
    3 & -1 & 0 & 0 \\
    -1 & 0 & 3 & 0 \\
    1 & -1 & -1 & 2
  \end{pmatrix} \, \mathbf{x} =
  \begin{pmatrix}
    8 \\ 5 \\ 0 \\ 1
  \end{pmatrix}.
```
The first row of this system simply states $4 x_1 = 8$,
which is very easy to solve, namely $x_1 = 8 / 4 = 2$.
The second row states $3x_1-x_2=5$.
However, $x_1$ is already known and can be inserted
to find $x_2 = -(5-3\cdot 2)=1$.
Following the same idea the third row gives
$x_3=(0+1\cdot 2)/3 = 2/3$ and the last row
$x_4=(1-1\cdot 2 + 1\cdot 1 + 1\cdot 2/3)/2 = 1/3$.
In total we found
```math
  \mathbf{x} =
  \begin{pmatrix} 2 \\ 1 \\ 2/3 \\ 1/3
  \end{pmatrix}.
```
Generalising to an arbitrary 4x4 lower-rectangular matrix
```math
  \begin{pmatrix}
    L_{11} & 0 & 0 & 0 \\
    L_{21} & L_{22} & 0 & 0 \\
    L_{31} & L_{32} & L_{33} & 0 \\
    L_{41} & L_{42} & L_{43} & L_{44}
  \end{pmatrix} \, \mathbf{x} =
  \begin{pmatrix}
    b_1 \\ b_2 \\ b_3 \\ b_4
  \end{pmatrix}
```
this solution algorithm can be expressed as
```math
\tag{3}
\begin{aligned}
  x_1 &= \frac{b_1}{L_{11}}, \\
  x_2 &= \frac{b_2 - L_{21}x_1}{L_{22}}, \\
  x_3 &= \frac{b_3 - L_{31}x_1 - L_{32}x_2}{L_{33}}, \\
  x_4 &= \frac{b_4 - L_{41}x_1 - L_{42}x_2 - L_{43}x_3}{L_{44}}.
\end{aligned}
```
We obtain

!!! info "Algorithm 1: Forward substitution"
    Given a lower-triangular matrix $\mathbf{L} \in \mathbb{R}^{n\times n}$
    a linear system $\mathbf{L} \mathbf x = \mathbf b$ can be solved
    by looping from $i = 1, \ldots n$ and computing
    ```math
    x_i = \frac{1}{L_{ii}} \left(b_i - \sum_{j=1}^{i-1} L_{ij} x_j\right)
    ```

or in form of an implementation
"""

# ╔═╡ ef0babcf-d7a8-4fe3-9d64-a54504d8da77
function forward_substitution(L, b)
	n = size(L, 1)
    x = zeros(n)
    x[1] = b[1] / L[1, 1]
    for i in 2:n
        row_sum = sum(L[i, j] * x[j] for j in 1:i-1)
        x[i] = (b[i] - row_sum) / L[i, i]
    end
    x
end

# ╔═╡ 6e6fce09-d14c-459d-9dd6-5f8fd1ee8248
md"""
For upper triangular matrices we solve linear systems proceeds
following the same idea
--- just in this case we start from the last row and not the first.
For example a system
```math
  \begin{pmatrix}
    U_{11} & U_{12} & U_{13} & U_{14} \\
    0 & U_{22} & U_{23} & U_{24} \\
    0 & 0 & U_{33} & U_{34} \\
    0 & 0 & 0 & U_{44}
  \end{pmatrix} \, \mathbf{x} =
  \begin{pmatrix}
    b_1 \\ b_2 \\ b_3 \\ b_4
  \end{pmatrix}
```
we solve by starting with $x_4$ and working our way forward to $x_1$:
```math
\tag{4}
\begin{aligned}
  x_4 &= \frac{b_4}{U_{44}}, \\
  x_3 &= \frac{b_3 - U_{34}x_4}{U_{33}}, \\
  x_2 &= \frac{b_2 - U_{23}x_3 - U_{24}x_4}{U_{22}}, \\
  x_1 &= \frac{b_1 - U_{12}x_2 - U_{13}x_3 - U_{14}x_4}{U_{11}}.
\end{aligned}
```
More formally this algorithm is called
!!! info "Algorithm 2: Backward substitution"
    Given an upper-triangular matrix $\mathbf{U} \in \mathbb{R}^{n\times n}$
    a linear system $\mathbf{U} \mathbf x = \mathbf b$ can be solved
    by looping in reverse order $i = n, \ldots, 1$ and computing
    ```math
    x_i = \frac{1}{U_{ii}} \left(b_i - \sum_{j=i+1}^{n} U_{ij} x_j\right)
    ```
or in code
"""

# ╔═╡ 8e31aaec-2512-4ad9-88a4-3130734d6e15
md"""
Based on our discussion we now understand why it is advantageous to
perform an LU factorisation when solving a linear system:
For both the resulting triangular matrix $\mathbf{L}$ as well as the matrix
$\mathbf{U}$, simple solution algorithms are available.
In summary we thus obtain

!!! note "Algorithm 3: Solving linear systems by LU factorisation / Gaussian elimination"
    We are given a matrix $\mathbf{A} \in \mathbb{R}^{n\times n}$,
    a right-hand side $\mathbf b \in \mathbb{R}^{n}$. This algorithm
    computes the solution $\mathbf x \in \mathbb{R}^{n}$
    to $\mathbf{A} \mathbf x = \mathbf b$:
    1. Factorise $\mathbf{A} = \mathbf{L} \mathbf{U}$.
    2. Solve $\mathbf{L} \mathbf{z} = \mathbf b$ for $\mathbf z$
       using *forward* substitution.
    3. Solve $\mathbf U \mathbf x = \mathbf z$ for $\mathbf x$
       using *backward* substitution.

A few remarks are in order:
- We have so far not discussed how to even compute the LU factorisation.
  As we will see in the next section this step will actually be accomplished
  by the Gaussian elimination algorithm, which probably already know
  from your linear algebra lectures.
- A severe drawback of this naive algorithm is immediately apparent from
  our 4x4 example problem, where (3) and (4) provide explicit solution
  expressions. If by any chance an element $L_{ii}$ or $U_{ii}$ happens
  to be zero, the algorithm cannot work. This already points to a problem
  we will discuss in the section on *pivoting*.
- Note, that steps 2. and 3. of the algorithm do *not depend* on the
  original matrix $\mathbf A$, but only on the right hand side $\mathbf b$.
  In other words as soon as the factorisation $\mathbf{A} = \mathbf{L} \mathbf{U}$
  has been computed, solving $\mathbf A \mathbf x = \mathbf b$ for many
  right hand sides $\mathbf b$ only requires the execution of steps 2. and 3.
  As we will see later, step 1. is the most expensive step. Therefore
  once $\mathbf{A}$ has been factorised, the cost for solving an
  *arbitrary* linear systems involving $\mathbf{A}$ is reduced as its
  factorised form $\mathbf{L} \mathbf{U}$ can be used in its place.
"""

# ╔═╡ d1c0d1ea-43c9-442e-899c-e83155f14841
md"""
## Gaussian elimination and LU factorisation

Final question and the missing puzzle piece is how to compute LU factorization.
As it turns out the Gaussian elimination algorithm you already learned in previous
linear algebra classes provides exactly this missing ingredient.

Recall that **Gaussian elimination** allows us to
transform a linear system into an upper
triangular system using elementary operations,
which combine the different equations.

!!! warning "Example: Manual Gaussian elimination"
	We will discuss one example, namely the linear system
    ```math
	\underbrace{\begin{pmatrix}
	2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4
	\end{pmatrix}}_{=\textbf A}
    \underbrace{\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix}}_{=\textbf x}
    = \underbrace{\begin{pmatrix}
	4 \\ 2 \\ -2
	\end{pmatrix}}_{= \textbf b}.
    ```
    Let us call $r_1^{(1)}$, $r_2^{(1)}$ and $r_3^{(1)}$
    the three equations of the system:
	```math
	\begin{aligned}
		& r_1^{(1)}: &&2x_1 &+&x_2 && &=& 4& \\
		& r_2^{(1)}: &-& 4x_1 &+&3x_2 &-&x_3 &=&2& \\
		& r_3^{(1)}: &&4x_1 &-&3x_2 &+&4x_3 &=&-2&
	\end{aligned}
	```
	on which we perform the following operations.

	**1st step:** Zero out all but the first row in the first column.
	```math
	\begin{aligned}
	    r_1^{(2)} &\leftarrow r_1^{(1)}  && \Longrightarrow 2x_1 &+&x_2 && &=& 4\\
	    r_2^{(2)} &\leftarrow r_2^{(1)} - \underbrace{\left(\boldsymbol{\frac{-4}{2}}\right)}_{L_{21}}r_1^{(1)} &&\Longrightarrow &&5x_2&-&x_3 &=& 10 \\
	    r_3^{(2)} &\leftarrow r_3^{(1)} - \underbrace{\left(\boldsymbol{\frac{4}{2}}\right)}_{L_{31}}r_1^{(1)} && \Longrightarrow &-&5x_2 &+&4x_3 &=& -10.
	\end{aligned}
	```

	**2nd step:** Zero out all but the first and second row in the second column.
	```math
	\begin{aligned}
    r_1^{(3)} &\leftarrow r_1^{(2)}  && \Longrightarrow 2x_1 &+&x_2 && &=& 4\\
    r_2^{(3)} &\leftarrow r_2^{(2)}  &&\Longrightarrow &&5x_2&-&x_3 &=& 10 \\
    r_3^{(3)} &\leftarrow r_3^{(2)} - \underbrace{\left(\boldsymbol{\frac{-5}{5}}\right)}_{L_{32}}r_2^{(2)} && \Longrightarrow && &&3x_3 &=& 0.
	\end{aligned}
	```

	After this step we have thus obtained an upper-triangular system,
    which we would solve using backward substitution.
    Note that during the iterations the system matrix
    got slowly reduced to an upper triangular matrix,
    which we can identify with $\mathbf{U}$:
	```math
	\underbrace{\textbf A^{(1)}}_{=\textbf A}=\begin{pmatrix} 2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4\end{pmatrix}  \; \Rightarrow \;
	\textbf A^{(2)}=\begin{pmatrix} 2 & 1 & 0 \\ \textcolor{grey}{0} & 5 & -1 \\ \textcolor{grey}{0} & -5 & 4\end{pmatrix}  \; \Rightarrow \;
	\textbf A^{(3)}=\begin{pmatrix} 2 & 1 & 0 \\ \textcolor{grey}{0} & 5 & -1 \\ \textcolor{grey}{0} & \textcolor{grey}{0} & 3\end{pmatrix}=\textbf U
    ```
	In turn the matrix $\textbf L$ is obtained by collecting the factors
    we employed to multiply the rows. We add diagonal entries of $1$
    to indicate that the $i$-th row has been left unchanged in the $i$-th step.
    Accumulating step by step we obtain
    ```math
	\mathbf L^{(1)} =   \begin{pmatrix} 1 &\phantom{-1}&\phantom{1} \\ -\frac{4}{2} &  &  \\ \frac{4}{2} &  &  \end{pmatrix} \; \Rightarrow \;
	\mathbf L^{(2)} =   \begin{pmatrix} 1 &\phantom{-1}& \\ -2 & 1 &  \\ 2 & -\frac{5}{5} &  \end{pmatrix} \; \Rightarrow \;
	\mathbf L^{(3)} =   \begin{pmatrix} 1 &\phantom{-1}&\phantom{1} \\ -2 &  1&  \\ 2 & -1 & 1 \end{pmatrix} = \mathbf L
    ```
    Finally by multiplying out the matrices we easily verify
    ```math
    \mathbf{L} \mathbf U = \begin{pmatrix} 1 &\textcolor{grey}{0}&\textcolor{grey}{0} \\ -2 &  1& \textcolor{grey}{0} \\ 2 & -1 & 1 \end{pmatrix}\begin{pmatrix}
	2 & 1 & 0 \\ \textcolor{grey}{0}& 5 & -1 \\ \textcolor{grey}{0} & \textcolor{grey}{0} & 3
	\end{pmatrix} = \begin{pmatrix}
	2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4
	\end{pmatrix} = \mathbf A
    ```

So effectively Gaussian elimination provides us with an approach to compute the LU factorisation. As a result Algorithm 3 is nothing but a formalised form of the procedure you used in your linear algebra lectures to solve linear systems using Gaussian elimination. We summarise:
"""

# ╔═╡ 1bcd8f38-5ae9-459e-8f91-8c9ad7cde407
md"""
!!! info "Algorithm 4: LU factorisation"
    **Input:** $\textbf A \in \mathbb{R}^{n\times n}$ 

    **Output:** $\mathbf U \in \mathbb{R}^{n\times n}$, $\mathbf L \in \mathbb{R}^{n\times n}$
    -  $\textbf A^{(1)} = \textbf A$
    - for $k = 1, \ldots, n-1$ *(algorithm steps)*
      *  $L_{kk} = 1$
      * for $i = k+1, \ldots, n$  *(Loop over rows)*
        -   $L_{ik} = \frac{A_{ik}^{(k)}}{A^{(k)}_{kk}}$
        - for $j = k+1, \ldots n$  *(Loop over columns)*
          *   $A_{ij}^{(k+1)} = A_{ij}^{(k)} - L_{ik} A_{kj}^{(k)}$
    -  $\textbf U = \textbf A^{(n)}$ 


An implementation of LU factorisation could be realised as such:
"""

# ╔═╡ 9692c144-531a-4995-9057-60af2b91ecfa
function factorise_lu(A)
    n = size(A, 1)
    L = diagm(ones(n))   # Initialise L by ones on diagonal, zeros elsewhere
	U = zeros(n, n)
    Aᵏ = float(copy(A))  # Make a copy of A and ensure that all entries
	                     # are converted to floating-point numbers

    for k in 1:n-1          # Algorithm steps
		U[k, :] = Aᵏ[k, :]  # Copy k-th row to U, since this row is in its final form.
		                    # From now on in the algorithm use U[k, :] instead of Aᵏ,
		                    # since Aᵏ is updated in-place to become Aᵏ⁺¹
		for i in k+1:n      # Loop over rows
			L[i, k] = Aᵏ[i, k] / U[k, k]
			for j = k+1:n   # Loop over columns
				Aᵏ[i, j] = Aᵏ[i, j] - L[i, k] * U[k, j]	 # Update Aᵏ in-place
			end
		end
    end
	U[n, n] = Aᵏ[n, n]

	# Return L and U using Julia datastructures to indicate
	# their special lower-triangular and upper-triangular form.
    return LowerTriangular(L), UpperTriangular(U)
end

# ╔═╡ 4da6a534-6a8f-464d-bbd4-9a6eec942688
md"""
Let us test our algorithms `factorise_lu`, `forward_substitution` and `backward_substitution` on our earlier example problem $\mathbf A \mathbf x = \mathbf b$ with
"""

# ╔═╡ 995a73db-6289-423b-8b5c-0fca4240472e
A

# ╔═╡ fa46c270-9eff-4fe4-9716-cf4cbc86d9c0
b

# ╔═╡ 5b6d6431-3b77-4b2e-ae19-69000602935b
md"""
We follow Algorithm 3. First we need to find the LU factorisation:
"""

# ╔═╡ 8f92f95a-9596-4226-8167-7db772b01ff4
L, U = factorise_lu(A)

# ╔═╡ d6a1c5e7-0921-4c4e-bc4a-fb2ebc94e45c
function backward_substitution(U, b)
	n = size(L, 1)
    x = zeros(n)
    x[n] = b[n] / U[n, n]
    for i in n-1:-1:1
        row_sum = sum(U[i, j] * x[j] for j in i+1:n)
        x[i] = (b[i] - row_sum) / U[i, i]
    end
    x
end

# ╔═╡ e66acbf7-a93b-491f-8b2a-9554abd33a6c
md"""
Next we forward substitute:
"""

# ╔═╡ 2e1043ec-f822-4c68-b173-fdaac26502db
z = forward_substitution(L, b)

# ╔═╡ 097763de-9f4b-4fa1-a64d-3938c46607c7
md"""
Finally we backward substitute:
"""

# ╔═╡ 890f4d32-04ff-4931-bc3c-e69005221183
x = backward_substitution(U, z)

# ╔═╡ 57279ca8-5824-445a-a837-6a7312a366a9
md""" and check the result:"""

# ╔═╡ f84cedc8-07b2-4744-81c1-a7b98a86a664
A * x - b

# ╔═╡ 50d0ccfb-78e4-432d-8f0a-abc18ccddacc
md"""
Hooray ! This suceeded. But does it always work ?
The answer to this is unfortunately now as we will see in the next section.
"""

# ╔═╡ 2baa3355-f51b-464e-a55a-3471fb7e0100
md"""
## Gaussian eliminitation with pivoting

Consider the matrix
"""

# ╔═╡ 036177c4-08dd-4c22-b3fb-ca0ee807b307
AA = [1 2 3;
      2 4 5;
	  7 8 9]

# ╔═╡ e35fe3ba-ae7d-407c-b034-5430449d414b
md"""the LU factorisation gives strange results"""

# ╔═╡ a5e8f3d1-246e-4b33-bceb-c682f3bb4fbc
let
	L, U = factorise_lu(AA)
	L
end

# ╔═╡ c9a3dffd-c418-48a0-ade1-4c6ac3df9a40
md"""The source of the problem is easily found. After the first step of the algorithm
 we get the matrices
```math
  \textbf A^{(2)} = \begin{pmatrix} 1 & 2 & 3\\0 & \textcolor{red}{0} & -1\\ 0 &-6 &-12\end{pmatrix}, \qquad
  \textbf L^{(1)} = \begin{pmatrix} 1 & & \\ 2 &\phantom{1} & \\ 7 & &\phantom{1} \end{pmatrix}.
```
As a result the element `Aᵏ[k, k]` is zero and if we blindly continue
in step `k = 2` a division
by zero is encountered in the computation of `L[k, k]`
--- which reflects above in the introduction of a `-Inf`
in the matrix. From this point our computations are numerical garbage.
Due their central role in the Gaussian elimination algorithm
the elements $\left(A^{(k)}\right)_{kk}$
--- respectively `Aᵏ[k, k]` in the implementation ---
are usually referred to as **pivots**.

To overcome the encountered problem we will consider
a strategy called **row pivoting**, i.e. where we allow ourselves
some flexibility by selecting the pivot amongst the $n-k$ lowermost
columns in the $k$-th step of Gaussian elimination.
This is usually performed by *permuting* (swapping the order)
of the $n-k$ lowermost columns.

Looking at $\textbf A^{(2)}$ it seems very reasonable to just swap the second and the
third column in $\textbf A^{(2)}$ and thus move the $-6$ to become the new pivot.
For consistency we not only need to swap $A^{(2)}$, but also $\textbf L^{(2)}$.
This gives us
```math
  \widetilde{\textbf{A}}^{(2)} = \begin{pmatrix} 1 & 2 & 3\\0 & \textcolor{red}{-6} &-12 \\
0 & 0 & -1\end{pmatrix}, \qquad
  \widetilde{\textbf{L}}^{(1)} = \begin{pmatrix} 1 & & \\ 7 &\phantom{1} & \\ 2 & &\phantom{1} \end{pmatrix}.
```
If we continue the algorithm now, the $-6$ sits in the position of the pivot
$\left(A^{(k)}\right)_{kk}$ and the division by zero is avoided.

In fact in this particular case there is nothing left to do and the matrix
$\widetilde{\textbf{A}}^{(2)}$ already is in upper triangular form.
We obtain from our algorithm
```math
\textbf U = \textbf A^{(3)} = \widetilde{\textbf{A}}^{(2)} = \begin{pmatrix} 1 & 2 & 3\\\textcolor{grey}{0} & -6 &-12 \\
\textcolor{grey}{0} & \textcolor{grey}{0} & -1\end{pmatrix}
\qquad
\textbf L = \textbf L^{(3)} = \begin{pmatrix}1& \textcolor{grey}{0} &\textcolor{grey}{0} \\ 7 & 1 & \textcolor{grey}{0} \\ 2 & 0& 1 \end{pmatrix}.
```
Since the additional row permutation we performed
not taken into account in the $\textbf L$ and $\textbf U$ factors,
multiplying them out will not yield $\textbf A$,
but a matrix where the second and third rows of $\textbf A$
are swapped:
```math
\mathbf L \mathbf U = 
\begin{pmatrix}1& \textcolor{grey}{0} &\textcolor{grey}{0} \\ 7 & 1 & \textcolor{grey}{0} \\ 2 & 0& 1 \end{pmatrix}
\begin{pmatrix} 1 & 2 & 3\\\textcolor{grey}{0} & -6 &-12 \\
\textcolor{grey}{0} & \textcolor{grey}{0} & -1\end{pmatrix} 
= \begin{pmatrix} 
 1 & 2 & 3 \\
 7 & 8 & 9 \\
 2 & 4 & 5 \\
\end{pmatrix}
```
To resolve this, we introduce the **permutation matrix**
```math
\textbf P = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0\end{pmatrix},
```
which is obtained by performing exactly the same row permutations
we did during the algorithm to the identity matrix
--- i.e. we swap the second and third row.
With this matrix we can easily check that
```math
\textbf L \textbf U = \textbf P \textbf A.
```

If we thus extend our `factorise_lu` function to additionally
perform such pivoting permutations,
the Gaussian elimination algorithm would terminate successfully.
But pivoting brings additional opportunities.
As it turns out numerical stability of LU factorisation
improves if one permutes not only if a pivot $\left(A^{(k)}\right)_{kk}$ is zero,
but in general if one *always* swaps the order of the rows
such that the pivot is taken as large as possible.
In other words in the $k$-th step of LU factorisation
we always exchange row $k$ with the row $l$ where $l$ satisfies
```math
\left|\left(A^{(k)}\right)_{lk}\right| \leq \left|\left(A^{(k)}\right)_{ik}\right|
\qquad
\text{for all $i = k, \ldots n$}.
```
The appropriate swaps are tracked and returned as well.
Note that instead of returning a permutation matrix
$\textbf P$ (which requries to store $n^2$ elements)
it is usually more convenient to store a permutation vector $\textbf p$,
which has only $n$ elements.
This vector tracks the indices of the rows of $\mathbf A$
in the order they are used as pivots. In other words if
"""

# ╔═╡ 6754911a-81dd-41ce-8b59-2a0743b324c0
idmx = diagm(ones(4))

# ╔═╡ c11a03f1-4aa2-4b6d-abf7-3e8c94d50498
md"""is the identity matrix and """

# ╔═╡ a67a89e9-858d-4501-834b-18f1ffd2ad0e
perm = [1, 3, 2, 4]

# ╔═╡ 5e1903be-3fd7-40dd-be3f-c8327c3c6033
md"""the permutation vector, then"""

# ╔═╡ 9dcfe511-f4a9-49ae-8571-33b402ca8d64
idmx[perm, :]

# ╔═╡ a25724a6-097a-48b5-a500-f4c5f2fe5388
md"""
is the permutation matrix.

The code below presents an implementation of row-pivoted LU factorisation.
"""

# ╔═╡ d09898da-dc12-4e0d-b548-72c8f32458e0
function factorise_lu_pivot(A)
    n = size(A, 1)
    L = zeros(n, n)
	U = zeros(n, n)
	p = fill(0, n)
    Aᵏ = float(copy(A))  # Make a copy of A and ensure that all entries
	                     # are converted to floating-point numbers

    for k in 1:n-1             # Algorithm steps
		p[k] = argmax(abs.(Aᵏ[:, k]))  # Find row with maximal pivot
		
		U[k, :] = Aᵏ[p[k], :]  # Copy pivot row to U, use U now instead of Aᵏ,
		                       # which is again updated in-place
		for i in 1:n           # Row loop: Note the full range as any row may
		                       #           be non-zero
			L[i, k] = Aᵏ[i, k] / U[k, k]
			for j = 1:n        # Column loop: Again full range
				Aᵏ[i, j] = Aᵏ[i, j] - L[i, k] * U[k, j]
			end
		end
    end
	p[n] = argmax(abs.(Aᵏ[:, n]))
	U[n, n] = Aᵏ[p[n], n]
	L[:, n] = Aᵏ[:, n] / U[n, n]

	# To simplify assembling L we so far kept the rows in the same order
	# as in A. To make the matrix upper triangular we also apply the column
	# permutation p before returning the results.
	LowerTriangular(L[p, :]), UpperTriangular(U), p
end

# ╔═╡ 24df3f6a-8ac3-492b-b4e7-80029af6918d
md"""We conclude with a theoretical result:

!!! info "Theorem 1"
    Every non-singular matrix $\textbf A \in \mathbb{R}^{n\times n}$ admits a factorisation
    ```math
    \textbf P \textbf A = \textbf L \textbf U
    ```
    where $\textbf L$ is lower-triangular, $\textbf U$ upper-triangular
    and $\textbf P$ is a permutation matrix.

Note, that this implies that `factorise_lu_pivot` can be employed to construct a robust method for solving linear systems, namely Gaussian elimination with pivoting.
The main idea is that using pivoted LU factorisation we can write:
```math
\textbf A \textbf x = \textbf b
\quad \Longleftrightarrow \quad
\textbf P \textbf A \textbf x = \textbf P  \textbf b
\quad \Longleftrightarrow \quad
\textbf L \textbf U \textbf x = \textbf P  \textbf b
```

!!! note "Algorithm 5: Solving linear systems with pivoted LU factorisation"
    Given a matrix $\mathbf{A} \in \mathbb{R}^{n\times n}$,
    a right-hand side $\mathbf b \in \mathbb{R}^{n}$
    the linear system $\mathbf{A} \mathbf x = \mathbf b$
    can be solved for  $\mathbf x \in \mathbb{R}^{n}$ as follows:
    1. Factorise $\textbf P \mathbf{A} = \mathbf{L} \mathbf{U}$.
    2. Solve $\mathbf{L} \mathbf{z} = \textbf P \mathbf b$ for $\mathbf z$
       using *forward* substitution.
    3. Solve $\mathbf U \mathbf x = \mathbf z$ for $\mathbf x$
       using *backward* substitution.


To conclude we return to our starting question:
Understanding the `\` operator in Julia.
- We already saw that julia essentially solves `A \ b` by replacing it by `lu(A) \ b`.
- Julia's `lu` by default also performs a row-pivoted LU factorisation.
- The returned object (an `LinearAlgebra.LU` object) internally stores
  both the `L`, the `U` and the permutation vector `p`
- When we use `\` with an `LU` object it, Julia immediately exploits the factorised form and performs first (permuted) forward substitution (step 2.) then backward substitution (step 3.).
- Therefore overall `A \ b` exactly performs Algorithm 5 above.

We can check against our implementation:
"""

# ╔═╡ 29d974a4-a748-4818-9c94-ca0a7b71e423
fac = lu(AA)

# ╔═╡ ae82d253-4937-4e72-b516-6f3da8f19df1
fac.p

# ╔═╡ b8045050-09e2-42ea-a7aa-c702ec323fd3
fac.L

# ╔═╡ f79f92bc-e50b-4e5c-a82e-555269fa80b5
factorise_lu_pivot(AA)[3]  # p

# ╔═╡ 8f06e886-6ebb-4a16-a1c5-8123e09ef723
factorise_lu_pivot(AA)[1]  # L

# ╔═╡ dd9023ad-cf79-4540-9259-61a28efa1174
TODO("Explain how to disable pivoting in julia and how to get the permuation matrix and so on")

# ╔═╡ 0f8c8b5c-33e0-4952-8d89-f26606a8503e
MM = [2 3  -4;
	  2 2  -2;
	  6 4 -11]

# ╔═╡ 2abb0086-0375-4633-8c1f-467b995ddd33
lu(MM).P

# ╔═╡ cee62fc9-1d60-4356-b170-d9928e7de807
([0 0 1; 0 1 0; 1 0 0] * [1 0 0; 0 0 1; 0 1 0])'

# ╔═╡ 8bd7c41a-112c-4ba7-b69e-af8e8d621c26
lu(MM, NoPivot())

# ╔═╡ bd574872-321d-4315-9f0f-60c95933226f
lu(MM)

# ╔═╡ 3e110a7e-1272-4c50-80f6-290f61844952
md"""
## Memory usage and fill-in

### Sparse matrices

The matrices and linear systems one encounters in physics and engineering applications not infrequently reach huge sizes. For example in the numerical solution of partial differential equations using finite difference or finite element methods the vector of unknowns frequently has a dimension in the order of $n \simeq 10^6$. This implies that matrices, i.e. entities from $\mathbb{R}^{n\times n}$ could hold around $10^{12}$ elements.

In double precision, i.e. $64$-bit numbers, each entry requires $8$ bytes. As a result storing all elements of such a matrix explicitly in memory reqires
"""

# ╔═╡ 1f34c503-6815-4639-a932-9fe2d7e1b4c2
10^12 * 8 / 1024^4  # TiB

# ╔═╡ ef60d36e-02c7-44e1-8910-241f99a89a38
md"""
tibibyes of storage. This is a challenge even for standard compute clusters, where typically a node has around $512$ GiB. Laptops nowadays feature around $32$ GiB, so would be completely out of the question to solve such problems.

However, in many applications the arising matrices are sparse,
meaning that they contain a large number of zero elements,
which do not need to be stored. Let us discuss a few examples.
"""

# ╔═╡ cea0c83c-6c6d-45a9-815c-6a112c3b76a9
md"""
!!! info "Definition: Full matrix"
    A matrix $A \in \mathbb{R}^{n\times n}$ is called **full** if the number
    of non-zero elements is at the order of $n^2$. For such matrices almost all
    elements are non-zero and need to be stored.

Full matrices are the "standard case" and for them memory constraints usually set the limit of the size of problems, which can be tackled. For example an $32$ GiB memory laptop can store around $4 \cdot 10^9$ double-precision numbers, which means that
linear problems with more than around $n \simeq 60000$ unknows cannot be treated.

!!! info "Definition: Sparse matrix"
    A matrix $A \in \mathbb{R}^{n\times n}$ is called **sparse** if the number
    of non-zero elements is at the order of $n$.

If we know which elements are zero, we can thus save memory by only storing those elements, which are non-zero. For this purpose the [`SparseArrays` Julia package](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) implements many primitives
for working with sparse arrays. 

This includes, generating random sparse arrays:
"""

# ╔═╡ 4a0dc42d-31a8-49bc-8db0-5f77474d8785
Asp = sprand(100, 100, 0.15)

# ╔═╡ d0a852cb-0d68-4807-861f-5302a5aeecd4
md"or simply sparsifying a dense array using the `sparse` function, which converts a full matrix to a sparse matrix by dropping the explicit storage of all zero entries."

# ╔═╡ 8353d1fd-75cc-41e0-866c-5234634219d5
M = [ 1 0 0  5;
     -2 0 1  0;
      0 0 6  0;
      0 1 0 -1]

# ╔═╡ d6b61d9a-817e-4c09-9a83-fde7ca591d23
sparse(M)

# ╔═╡ fb87aa13-1fd6-4c13-8ba4-68caef5f775b
md"""
Using the `SparseArray` data structure from `SparseArrays` consistently allows to fully exploit sparsity when solving a problem. As a result the storage costs scale only as $O(n)$. With our laptop of 32 GiB memory we can thus tackle problems with around $n\simeq 4 \cdot 10^9$ unknowns --- much better than the $60000$ we found when using full matrices.

Finally we introduce a special kind of sparse matrix:

!!! info "Definition: Band matrix"
    A matrix $A \in \mathbb{R}^{n\times n}$ is called a **band matrix**
    with bandwidth $d$ if $A_{ij} = 0$ when $|j - i| > d$.
    Every line of the matrix contains at most $2d + 1$ non-zero elements
    and the number of non-zeros thus scales as $O(nd)$.

An example for a banded matrix with bandwidth $5$ is:
"""

# ╔═╡ 5dc919c8-01e1-4b23-b406-492562c9e338
band = spdiagm(-5 =>  -ones(100),
	           -3 =>  3ones(102),
	           -2 =>   ones(103),
	            0 => -10ones(105),
	            1 =>  -ones(104),
	            3 =>  3ones(102),
	            4 =>   ones(101),
	            5 =>  -ones(100))

# ╔═╡ 74709cb2-0be8-4b24-aa9d-926ef059fe2d
md"""
### LU factorisation of sparse matrices

Since the amount of available memory can put hard constraints on the size of linear systems which can be solved, we now want to investigate the memory requirement of LU factorisation $\textbf A = \textbf L \textbf U$.

If $\textbf A$ is a full matrix, then we know that $\textbf L$ and $\textbf U$ are triangular, thus contain less non-zero elements than  $\textbf A$ itself.

Let's contrast this with the sparse matrices we have considered above.
First the random sparse matrix:
"""

# ╔═╡ 68396ec1-444b-4c56-966b-5c70a2208d34
Asp

# ╔═╡ 8b51fcc1-6287-43cd-a44c-26c7423e8d7a
lu(Asp).L  # U looks similar

# ╔═╡ 2429c458-7108-4f13-aad2-dd04a3f1c3cf
md"""Storing both L and U thus requires a total number of"""

# ╔═╡ ba2a91f4-7f7a-4393-90da-3d5848d379e1
nnz(lu(Asp).L) + nnz(lu(Asp).U)

# ╔═╡ d7063c34-6418-4e0b-80df-ddd07566e578
md"""non-zero elements, which is about **4 times as much** as the original matrix !"""

# ╔═╡ 0033e298-a85f-4718-acfa-fe9fd5a67746
md"Now the banded matrix:"

# ╔═╡ 84a76d06-bf4e-461b-904c-40614c5b9626
band

# ╔═╡ 39f1d6d1-5d53-4c92-b4f0-8abd55a4fc57
lu(band).L  # U looks similar

# ╔═╡ 3c384fb9-fa39-441c-81d6-7306df13657d
md"At least the banded structure is kept, but still we require"

# ╔═╡ 1cb81418-f8cb-4b4d-9343-586b5950f0da
nnz(lu(band).L) + nnz(lu(band).U)

# ╔═╡ 60ca2879-ee47-4163-b014-4a0af90b449f
md"nonzeros, which is about 1.5 times than the original in this case."

# ╔═╡ 191fb08e-fade-43cd-8d8f-2273b594da68
md"In both cases storing the $\mathbf L$ and the $\mathbf U$ factors require more non-zero elements than the original matrix. This phenomenon is usually referred to as **fill in**.

In particular for the sparse matrix `Asp` LU factorisation did not preserve the structure at all. In contrast it almost resulted in a full matrix in the lower right corner of the $\textbf L$. For $\mathbf U$ this looks exactly the same. In fact one can show that in general even for sparse matrices the **memory usage of LU factorisation** is $O(n^2)$. Therefore while one may be able to store a sparse matrix $\mathbf A$ with a huge size in memory (due to the $O(n)$ memory cost), one may actually run out of memory while computing the factorisation.
"

# ╔═╡ f3fee115-f239-4578-81e0-a3153027056f
md"""
Let us add that for banded matrices the situation is slightly better as one can show that in this case the fill-in takes place at most inside the band. As a result the memory requirement stays at $O(n d)$.
"""

# ╔═╡ ca3eca34-3fc6-4488-ae74-69e9453c05e1
md"""
## Computational cost of LU factorisation

In this section we want to investigate how long it will take a computer to perform LU factorisation on a large matrix.

Modern laptop computers nowadays have a clock frequence fo a few Gigahertz (GHz). This means they are able to perform about $10^9$ operations per second, where for simplicity we assume that "one operation" is an elementary addition, multiplication, division and so on. If we can determine how many such operations are needed to perform LU factorisation we can estimate the computational cost.

Before we look at the LU algorithm (Algorithm 4) we first understand a few simpler cases from linear algebra:
"""

# ╔═╡ c0b02aa0-642c-4a6f-ac21-cdb8e27d9279
md"""
### Scalar product
Given two vectors $x, y \in \mathbb{R}^n$ consider computing
the scalar product
```math
x \cdot y = x^T y = \sum_{i=1}^n x_i \, y_i
```
which in code is achieved as
"""

# ╔═╡ 48a19e05-ba2e-421a-b059-07d994bed73c
function scalar_product(x, y)
	result = 0.0
	for i in 1:length(x)
		result = result + x[i] * y[i]
	end
	result
end

# ╔═╡ e91b91ec-77b4-4c77-b2ce-0e022f54edb7
md"""
In the loop for each iteration we require $1$ multiplication and $1$ addition. In total the `scalar_product` function therefore requires $n$ multiplications and $n$ additions. The number of elementary operations is thus $2n$. We say the computational cost is $O(n)$ (i.e. on the order of $n$ or linear in $n$).

If we take the dimensionality $n=1000$ the number of operations is $O(1000)$. On a 1 GHz computer (with $10^9$ operations per second) this therefore takes about $10^{-6}$ seconds ... hardly noticable.
"""

# ╔═╡ 44b90b25-bdcd-4f38-804d-53e466ec6374
md"""
### Matrix-vector product

Given a matrix $A \in \mathbb{R}^{n\times n}$ and a vector $x \in \mathbb{R}^n$
the matrix-vector product $y = Ax$ is computed as
```math
y_i = \sum_{j=1}^{n} A_{ij} x_j \qquad \text{for $i = 1, \ldots, n$}
```
or in code
"""

# ╔═╡ 0f8b23ab-c797-4513-b70c-5827376f6093
function matrix_vector_product(A, x)
	y = zeros(size(A, 1))  # Create a vector of zeros (#)
	
	for i in 1:size(A, 1)
		for j in 1:size(A, 2)
			y[i] = y[i] + A[i, j] * y[j]  # (*)
		end
	end
	result
end

# ╔═╡ 756a5787-6f64-4308-b879-da7676d39a8c
md"""
In the innermost loop we observe again that we require
$1$ addition and $1$ multiplication per iteration.
This instruction is performed once for each combination of $i$ and $j$,
so we look at the limits of each of the nested loops.
The loop over $i$ runs over $n$ values (the number of rows of $A$)
and the loop over $j$ over $n$ values as well (the number of columns of $A$).
In total the inner most instruction `(*)` is thus run $n^2$ times,
each times costing $2$ operations.
The total cost is thus $O(n^2)$.

For our case with $n = 1000$ and a 1 GHZ computer we thus now need
$(10^3)^2 / 10^{9} s = 10^{-3} s = 1ms$, which is again a rather short time.

"""

# ╔═╡ 01a66e18-bc56-4204-b4b8-9bf46155aec1
Foldable("Expert information: What about the allocation in `(#)` ?",
md"""
You may wonder about the line marked `(#)` in the above code example,
which allocates the output vector $y$ with space for storing $n$ elements.

In fact allocating memory also costs computational time, which usually scales linearly with the size of the vector. So step `(#)` also costs $O(n)$,
such that the precise scaling of the cost of `matrix_vector_product` is actually
$O(n + n^2)$. Since the $O$-notation actually makes a statement about the *asymptotic complexity* it turns out that all lower terms can be neglected.
In fact $O(n^2 + n) = O(n^2)$, see [Big-O notation](https://en.wikipedia.org/wiki/Big_O_notation) for details.

For this course if you are asked to determine the asymptotic complexity,
you will not asked to deal with these details. For our problems it will be sufficient to only look at the most deeply nested loop and ignore the rest of the code as part of a complexity analysis.
""")

# ╔═╡ 4fac9508-bc29-447c-9188-a58da07abeb3
md"""
### LU factorisation

We revisit Algorithm 4:

!!! info ""
    **Input:** $\textbf A \in \mathbb{R}^{n\times n}$ 

    **Output:** $\mathbf U \in \mathbb{R}^{n\times n}$, $\mathbf L \in \mathbb{R}^{n\times n}$
    -  $\textbf A^{(1)} = \textbf A$
    - for $k = 1, \ldots, n-1$ *(algorithm steps)*
      *  $L_{kk} = 1$
      * for $i = k+1, \ldots, n$  *(Loop over rows)*
        -   $L_{ik} = \frac{A_{ik}^{(k)}}{A^{(k)}_{kk}}$
        - for $j = k+1, \ldots n$  *(Loop over columns)*
          *   $A_{ij}^{(k+1)} = A_{ij}^{(k)} - L_{ik} A_{kj}^{(k)}$
    -  $\textbf U = \textbf A^{(n)}$ 

We notice that the instruction at the innermost loop level is
```math
  A_{ij}^{(k+1)} = A_{ij}^{(k)} - L_{ik} A_{kj}^{(k)},
```
which costs again $1$ addition and $1$ multiplication
and which is executed once for each tuple $(i, j, k)$.
In the worst case $k$ is takes $n - 1$ elements,
$i$ also $n-1$ elements (assume $k=1$)
and $j$ again in its worst case takes $n-1$ elements.
In total the cost is thus $(n - 1)^3 = O(n^3)$.

For our case with $n = 1000$ and a 1 GHZ computer the computation
thus needs approximately $(10^3)^3 / 10^{9} s = 1 s$,
which starts to be a noticable amount of time.
If we even consider $n = 10^5$ all of a sudden it takes $10^3 s$,
which is about $15$ minutes.
For large matrices the cost of computing
LU factorisation can thus become important !

For **banded matrices** the computational scaling is a little improved.
E.g. consider LU factorisation on a banded matrix $A$ with band width $d$.
Then bothe the loop over $i$ (row loop) as well as the loop over $j$ (column loop) in the LU factorisation algorithm can be truncated and at most run over $d$ elements. As a result the computational cost at worst becomes $O(n d^2)$ which for a small band width (i.e. $d < n$) is substantially less than $O(n^3)$.

!!! info "Overview of LU factorisation computational cost and memory cost"
    We summarise the cost of LU factorisation in a table:

	|  type of $n\times n$ matrix   | computational cost | memory usage |
	| ------------------------- | ------------------ | ------------ |
	|  full matrix  | $O(n^3)$           | $O(n^2)$      |
	|  banded matrix, band width $d$  | $O(n\,d^2)$    | $O(n\,d)$      |

"""

# ╔═╡ e8de3cd9-05bb-439e-b6fe-39d9c82f0ba6
md"""
## Numerical stability

To close off our discussion of Gaussian elimination we consider the question of its numerical stability, i.e. how it is effected by the round off errors introduced in the algorithm by finite-precision floating-point arithmetic.

In this chapter our goal has been to solve a linear system of the form
```math
\tag{5}
\mathbf{A} \mathbf{x} = \mathbf{b},
```
where we are given a matrix $\mathbf{A} \in \mathbb{R}^{n\times n}$
as well as a right-hand side $\mathbf{b} \in \mathbb{R}^n$
and we wish to solve for $\mathbf{x}$
We will later refer to this system as the **exact system**.

Since the computer in general is unable to represent the matrix and right-hand side exactly in finite-precision floating point arithmetic,
even providing this problem to a computer is already associated with making a small error: in practice the computer is only able to obtain a solution
 $\widetilde{\mathbf{x}}$ to a slighly modified linear system
```math
\tag{6}
\widetilde{\mathbf{A}} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}},
```
where $\widetilde{\mathbf{A}}$ is an approximation to $\mathbf{A}$,
$\widetilde{\mathbf{b}}$ is an approximation to $\mathbf{b}$.
This linear system we will refer to as the **perturbed system**.

In standard floating-point arithmetic the *relative error* in representing
any number is around $10^{-16}$. Mathematically we can find the following
relationships between the elements of $\widetilde{\mathbf{b}}$
and $\mathbf{b}$ as well as $\widetilde{\mathbf{A}}$
and $\mathbf{A}$, respectively:
```math
\begin{aligned}
\widetilde{b}_i    &= b_i   \, (1 + \varepsilon_i) && \text{where $\varepsilon_i$ is on the order of $10^{-16}$} \\
\widetilde{A}_{ij} &= A_{ij} \, (1 + \eta_{ij}) && \text{where $\eta_{ij}$ is on the order of $10^{-16}$}
\end{aligned}
```
To understand the numerical stability when solving linear systems
our goal is now to understand how far the solution
$\widetilde{\mathbf{x}}$ --- obtained on a computer ---
differs from $\mathbf{x}$ --- the solution to (5), the true linear system
we want to solve.

To simplify matters in this course we will not attempt to answer this question in full generality, but concentrate on the case $\eta_{ij} = 0$, i.e. a setting where we only consider the rounding-error in the right-hand side $\widetilde{\mathbf{b}}$
and $\mathbf{b}$, but assume that the system matrix is accessible to the computer without error ($\widetilde{\mathbf{A}} = \mathbf{A}$).
For an alternative and more detailed discussion see also
[chapter 2.8](https://tobydriscoll.net/fnc-julia/linsys/condition-number.html) of Driscoll, Brown: *Fundamentals of Numerical Computation*.
"""

# ╔═╡ 97e1379a-b0bf-4c6f-8375-299e1f42899f
md"""
To provide a motivation why this type of analysis matters in practice,
we first consider the following example:

!!! warning "Example: Rounding error in a 2x2 system"
    Consider the *exact* linear system
    ```math
	\mathbf{A} \mathbf{x} = \mathbf{b} \qquad \Leftrightarrow \qquad
	\begin{pmatrix} 1 & 10^{-16} \\ 1 &0\end{pmatrix}
	\begin{pmatrix} x_1 \\ x_2 \end{pmatrix}
	= \begin{pmatrix} 1\\1 \end{pmatrix}
    ```
	with solution $(x_1, x_2) = (1, 0)$ and the *perturbed* linear system
	```math
	\mathbf{A} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}} \qquad \Leftrightarrow \qquad
	\begin{pmatrix} 1 & 10^{-16} \\ 1 &0\end{pmatrix}
	\begin{pmatrix} \widetilde{x}_1 \\ \widetilde{x}_2 \end{pmatrix}
	= \begin{pmatrix} 1 + 10^{-16}\\1 \end{pmatrix},
	```
	where only the first component of the right-hand side has been perturbed
	by $10^{-16}$. This second system has solution
	$(\widetilde{x}_1, \widetilde{x}_2) = (1, 1)$.
	The second component of this solution is thus *completely wrong*!

	The small perturbation of $10^{-16}$,
	which may readily occurr just by inputing the data to a computer,
	can already have an impactful effect on the solution and in fact
	render the solution we obtain from our our computation
	(i.e. $\widetilde{\mathbf{x}}$) rather far from the answer we are actually after
	(i.e. $\mathbf{x}$).
"""

# ╔═╡ f11b6f84-d0c4-4855-a5c6-31d837880407
md"""
### Matrix norm and condition numbers

To analyse the error $\mathbf{x} - \widetilde{\mathbf{x}}$ mathematically, we first have to introduce some notation.

Recall that for a vector $\mathbf{v} \in \mathbb{R}^n$ its Euclidean norm
is $\|\mathbf{v}\| = \sqrt{\sum_{i=1}^n v_i^2}$.

We define the following:

!!! info "Definition: Relative error of linear systems"
    The **relative error** in $\widetilde{\mathbf{x}}$,
    the solution to the perturbed system (6),
    is the quantity
	```math
	\frac{\|\mathbf{x} - \widetilde{\mathbf{x}}\|}{\|\mathbf{x}\|},
	```
	where $\mathbf{x}$ is the exact solution,
    i.e. the solution to the exact linear system (5) one actually wanted to solve.

We further need a generalisation of norms to matrices:

!!! info "Definition: Matrix norm"
    Given $\mathbf{M} \in \mathbb{R}^{m \times n}$ a real matrix (not necessarily square), we define the matrix norm of $\mathbf{M}$ as
	```math
		\|\mathbf{M}\| = \max_{\stackrel{\mathbf{v} \in \mathbb{R}^n}{\mathbf{v}\neq\mathbf{0}}} \frac{\|\mathbf{M}\mathbf{v}\|}{\|\mathbf{v}\|}.
	```

The previous definition implies in particular that
```math
\|\mathbf{M}\mathbf{v}\| \leq \|\mathbf{M}\| \, \|\mathbf{v}\|
\quad \forall \mathbf{x} \in \mathbb{R}^n.
\tag{7}
```
"""

# ╔═╡ 0e4e510c-9662-4962-aeea-2a31f6c027fe
md"""
With this definitions in place we begin our analysis for (5) and (6) where $\mathbf{A} \in \mathbb{R}^{n\times n}$.
From the definition of the exact and perturbed linear systems
```math
\mathbf{A} \mathbf{x} = \mathbf{b} \qquad\text{and}\qquad \mathbf{A} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}}
```
we obtain by subtraction and since $\mathbf{A}$ is invertible
```math
\mathbf{A} \left( \mathbf{x} -  \widetilde{\mathbf{x}} \right)
= \mathbf{b} -  \widetilde{\mathbf{b}} 
\qquad \Rightarrow \qquad
\mathbf{x} -  \widetilde{\mathbf{x}} = \mathbf{A}^{-1} \left( \mathbf{b} -  \widetilde{\mathbf{b}}  \right)
```
Using (7) we therefore have
```math
\|\mathbf{x} -  \widetilde{\mathbf{x}}\| \leq \left\| \mathbf{A}^{-1}\right\| \left\| \mathbf{b} -  \widetilde{\mathbf{b}}  \right\|.
```
Furthermore we have from applying (7) to $\mathbf{A} \mathbf{x} = \mathbf{b}$:
```math
\|\mathbf{b}\| \leq \|\mathbf{A}\| \, \|\mathbf{x}\|
\qquad \Rightarrow \qquad
\frac{1}{\|\mathbf{x}\|} \leq \frac{\|\mathbf{A}\|}{\|\mathbf{b}\| }
```
such that with the above result we can bound the relative error as
```math
\tag{8}
\frac{\|\mathbf{x} -  \widetilde{\mathbf{x}}\|}{\|\mathbf{x}\|}
\leq \left\| \mathbf{A}^{-1}\right\| \left\| \mathbf{b} -  \widetilde{\mathbf{b}}  \right\| \cdot \frac{\|\mathbf{A}\|}{\|\mathbf{b}\|}
= \underbrace{\left\| \mathbf{A}^{-1}\right\| \, \|\mathbf{A}\|}_{=\kappa(\mathbf{A})}
\, \frac{\left\| \mathbf{b} -  \widetilde{\mathbf{b}} \right\|}{\|\mathbf{b}\|}.
```

From as stability analysis point of view the quantity $\kappa(\mathbf{A}) = \left\| \mathbf{A}^{-1}\right\| \, \|\mathbf{A}\|$ tus relates the relative error in the right-hand side (input quantity)
to the relative error in the solution $\mathbf{x}$ (output quantity).
It is thus the condition number for a linear system:

!!! info "Definition: Condition number of a linear system"
	Given a linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$
	with an invertible system matrix $\mathbf{A} \in \mathbb{R}^{n \times n}$
	its **condition number** is defined as
	```math
	\kappa(\mathbf{A}) = \left\| \mathbf{A}^{-1}\right\| \, \left\| \mathbf{A}\right\|.
	```
	Due to the importance of solving linear systems in linear algebra
	one usually calls $\kappa(\mathbf{A})$ also the **condition number
	of the matrix $\mathbf{A}$**.

But notably (8) can be interpreted even more generally as the way
it the solution of linear systems are changed as we change the right-hand side
as is summarised below:

!!! info "Theorem 2: Change of solution when changing the right-hand side"
	Given a linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$,
	which we solve for $\mathbf{x}$ and a related
	linear system $\mathbf{A} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}}$
	which we solve for $\widetilde{\mathbf{x}}$,
	then the two solutions are related as
	```math
	\frac{\|\mathbf{x} -  \widetilde{\mathbf{x}}\|}{\|\mathbf{x}\|}
	\leq \kappa(\mathbf{A})
	\, \frac{\left\| \mathbf{b} -  \widetilde{\mathbf{b}} \right\|}{\|\mathbf{b}\|}.
	```
"""

# ╔═╡ d326776b-9a89-4464-b720-92d66e93e83c
md"""
Before we continue to develop expression (8) we pause for a second to discuss how condition numbers can be computed in practice.

First let us introduce the following notation:

!!! info "Definition: Minimal and maximal eigenvalues"
	Given a diagonalisable matrix $\mathbf{M} \in \mathbb{R}^{n\times n}$
	we denote by $\lambda_\text{max}(\mathbf{M})$ the *largest* eigenvalue
	of $\mathbf{M}$ and by $\lambda_\text{min}(\mathbf{M})$ the *smallest* eigenvalue
	of $\mathbf{M}$, i.e.
	```math
	\lambda_\text{max}(\mathbf{M}) = \max_{i=1,\ldots,n} \lambda_i(\mathbf{M})
	\qquad\text{and}\qquad
	\lambda_\text{min}(\mathbf{M}) = \min_{i=1,\ldots,n} \lambda_i(\mathbf{M})
	```
	where $\lambda_i(\mathbf{M})$ for $i=1,\ldots,n$ are the eigenvalues of $\mathbf{M}$.

It turns out that the minimal and maximal eigenvalues provide
a convenient way to compute matrix norms:

!!! info "Lemma 3: Computing matrix norms"
	For any matrix  $\mathbf{M} \in \mathbb{R}^{m \times n}$
	```math
	\|\mathbf{M}\| = \sqrt{ \lambda_\text{max}(\mathbf{M}^T\mathbf{M}) }.
	```
	If $\mathbf{M}$ is moreover square and
	invertible, then
	```math
	\|\mathbf{M}^{-1}\| = \sqrt{ \lambda_\text{max}(\mathbf{M}^{-T}\mathbf{M}^{-1}) } = \frac 1 {\sqrt{ \lambda_\text{min}(\mathbf{M}^T\mathbf{M}) }}.
	```
"""

# ╔═╡ ab532d5b-4396-4792-b5e5-232216d6e4ca
md"""
Based on this we can deduce a few useful formulas for computing condition numbers:

- For any invertible $\mathbf{A} \in \mathbb{R}^{n\times n}$ the condition number
  can be computed as
  ```math
  \kappa(\mathbf{A}) = \left\| \mathbf{A}^{-1}\right\| \, \left\| \mathbf{A}\right\|
  = \frac{ \sqrt{ \lambda_\text{max}(\mathbf{A}^T\mathbf{A}) }}
  { \sqrt{ \lambda_\text{min}(\mathbf{A}^T\mathbf{A}) }}.
  ```
- If $\mathbf{A} \in \mathbb{R}^{n\times n}$ is invertible and **symmetric**
  then $\mathbf{A}^T = \mathbf{A}$. As a result we observe that
  ```math
  \lambda_{i}(\mathbf{A}^T\mathbf{A}) = \lambda_{i}(\mathbf{A}^2) = \lambda_{i}(\mathbf{A})^2 \qquad \text{for all $i=1,\ldots,n$}.
  ```
  As a result
  ```math
  \kappa(\mathbf{A})
  = \frac{ \sqrt{ \lambda_\text{max}(\mathbf{A}^T\mathbf{A}) }}
  { \sqrt{ \lambda_\text{min}(\mathbf{A}^T\mathbf{A}) }}
  = \frac{ \sqrt{ \max_{i=1,\ldots,n} \lambda_{i}(\mathbf{A})^2}}
     { \sqrt{\min_{i=1,\ldots,n} \lambda_{i}(\mathbf{A})^2}}
  = \frac{ \max_{i=1,\ldots,n} |\lambda_{i}(\mathbf{A})|}
     { \min_{i=1,\ldots,n} |\lambda_{i}(\mathbf{A})|}
  ```
- If $\mathbf{A} \in \mathbb{R}^{n\times n}$ is **both symmetric
  and positive definite** (i.e. all eigenvalues are strictly positive),
  then we can drop the moduli in the above expression and obtain:

!!! info "Lemma 4: Condition number of symmetric, positive-definite matrices"
	If $\mathbf{A} \in \mathbb{R}^{n\times n}$ is symmetric and positive definite,
    then
    ```math
	\kappa(\mathbf{A}) = \frac{\lambda_\text{max}(\mathbf{A})}
     { \lambda_\text{min}(\mathbf{A})}
	```
"""

# ╔═╡ f2e5abd8-a544-4c21-88fe-af0edec938ab
md"""
!!! warning "Rounding error in a 2x2 system (continued)"
	We continue our example from earlier and evaluate the condition number
	of the matrix
	```math
	\mathbf{A} = \begin{pmatrix} 1 & 10^{-16} \\ 1 &0\end{pmatrix}.
	```
	This matrix is not symmetric and we therefore use the expression
	```math
	\kappa(\mathbf{A}) = \frac { \sqrt{ \lambda_\text{max}(\mathbf{A}^T\mathbf{A}) }}  { \sqrt{ \lambda_\text{min}(\mathbf{A}^T\mathbf{A}) }},
	```
	which requires us to compute the eigenvalues of $\mathbf{A}^T\mathbf{A}$.
	We note
	```math
	\mathbf{A}^T\mathbf{A})
	= \begin{pmatrix} 1 & 1 \\ 10^{-16} &0\end{pmatrix}
	\begin{pmatrix} 1 & 10^{-16} \\ 1 &0\end{pmatrix}
	= \begin{pmatrix} 2 & 10^{-16}  \\ 10^{-16} & 10^{-32} \end{pmatrix}
	```
	and
	```math
	\det(\mathbf{A}^T\mathbf{A}) - \lambda I)
	= \lambda^2 - (2 +  10^{-32}) \lambda +  10^{-32},
	```
	such that the two eigenvalues of $\mathbf{A}^T\mathbf{A})$ are
	```math
	\lambda_{\pm}(\mathbf{A}^T\mathbf{A})) = \frac12\left( 2 +  10^{-32} \pm \sqrt{4 +  10^{-64}} \right).
	```
	We conclude that
	```math
	\lambda_\text{max}(\mathbf{A}^T\mathbf{A}) \approx 2
	\qquad\text{and}\qquad
	\lambda_\text{min}(\mathbf{A}^T\mathbf{A}) \approx \frac{10^{-32}}{2}
	```
	therefore
	```math
	\kappa(\mathbf{A}) \approx 2 \cdot 10^{16},
	```
	which is huge !
"""

# ╔═╡ 3eb70fdb-457a-4b5b-94d8-13b062895cd3
md"""
In Julia code the condition number of a matrix is computed programmatically using the `cond` function. It automatically makes a guess to use the most appropriate of teh above formulas. For the above example we can compute:
"""

# ╔═╡ 6cb42cac-faa0-4bc4-bb2c-54195f07fd9a
let
	A = [1 1e-16;
	     1  0   ]
	κ = cond(A)
end

# ╔═╡ 9571044c-7d41-4dd0-8601-fd5f37f05f0c
md"""
### Stability result
With these computational tools let us return to the previously obtained
relationship (8) between relative error in the right-hand side
and relative error in the solution of the exact and perturbed systems:

```math
\frac{\|\mathbf{x} -  \widetilde{\mathbf{x}}\|}{\|\mathbf{x}\|}
\leq \kappa(\mathbf{A})
\, \frac{\left\| \mathbf{b} -  \widetilde{\mathbf{b}} \right\|}{\|\mathbf{b}\|}.
```

In our setup we said that storing $\mathbf{b}$ on the computer
(which leads to $\widetilde{\mathbf{b}}$)
introduces a small relative error at each of the elements of the right-hand side,
leading to a relationship
```math
\widetilde{b}_i = b_i \, (1 + \varepsilon_i).
```
Let us further define
```math
\epsilon = \max_{i=1,\ldots,n} |\varepsilon_i|.
```
For standard double-precision floating-point number
as discusset at the beginning of this section,
we have that $\epsilon \approx 10^{-16}$.

With this setup we can simplify the
expression of the absolute error of right-hand side:
```math
\left\| \mathbf{b} -  \widetilde{\mathbf{b}}  \right\|
= \sqrt{\sum_{i=1}^n b_i^2 \varepsilon_i^2}
\leq \sqrt{\sum_{i=1}^n b_i^2 \underbrace{\max_{i=1,\ldots,n} |\varepsilon_i|^2}_{=\epsilon^2}}
= \epsilon \, \| \mathbf{b} \|.
```

This leads to the following stability result:

!!! info "Theorem 5: Stability of solving linear systems"
	Given a linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$
	and $\widetilde{\mathbf{x}}$ the solution to a perturbed
	linear system $\mathbf{A} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}}$,
	where $\widetilde{b}_i = b_i \, (1 + \varepsilon_i)$ for $i=1,\ldots,n$
	and $\epsilon = \max_{i=1,\ldots,n} |\varepsilon_i|$, then
	```math
	\frac{\left\| \mathbf{x} -  \widetilde{\mathbf{x}}  \right\|}
	{\| \mathbf{x}\|} \leq \kappa(\mathbf{A}) \, \epsilon
	```

This inequality shows that the condition number of the matrix $\mathbf{A}$ plays the role of *amplifying* the round-off error introduced by the floating-point representation of the right-hand side.

!!! warning "Rounding error in a 2x2 system (continued)"
	Using Theorem 5 we can finally understand the behaviour we observe
	in our example.
	Previously we found that a small perturbation
	$1 \to 1 + 10^{-16}$ in one of the elements of the right hand side,
	introduces a change in the solution from
	$(x_1, x_2) = (1, 0)$ to $(\widetilde{x}_1, \widetilde{x}_2) = (1, 1)$.
	This is a relative error of
	```math
	\frac{\left\| \mathbf{x} -  \widetilde{\mathbf{x}}  \right\|}
	{\| \mathbf{x}\|} = \frac{1}{1} = 1.
	```
	Note that in this example $\epsilon \approx 10^{-16}$ just like
	for the case of a rounding error in the right-hand side.

	As we found previously
	```math
	\kappa(\mathbf{A}) \approx 2 \cdot 10^{16},
	```
	such that this large relative error is explained by the large
	condition number of this matrix.
"""

# ╔═╡ 45df2cc2-44d7-4b39-872d-b3a49cd59dfe
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 525)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
HypertextLiteral = "~0.9.5"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "8e611390a0aadd23aca741d761c4cd296c83671d"

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

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "5d3a5a206297af3868151bb4a2cf27ebce46f16d"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.33"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "5b0d630f3020b82c0775a51d05895852f8506f50"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.4"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "0b898aba6cb0b01fb96245fa5375accb651a241a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.0.0"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "677b65e17aeb6b4a0be1982e281ec03b0f55155c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.16"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─ca2c949f-a6a0-485f-bd52-5dae3b050612
# ╠═3295f30c-c1f4-11ee-3901-4fb291e0e4cb
# ╟─21c9a859-f976-4a93-bae4-616122712a24
# ╟─8e4312b5-f4ce-440b-a487-689959851098
# ╠═fafa019e-8a47-4a80-bfd8-a017ae31ab1f
# ╠═e45952b2-714b-427b-8979-98d11c330294
# ╟─52d8634e-116b-48e6-9673-2ee2a97423c0
# ╠═5fb995cc-b338-4608-a01c-fc0e84d9dfe9
# ╟─4d9d19da-7a4f-49ba-9c6f-890dad00672c
# ╠═08ec21d9-5474-4aae-8fb1-62fde1c7e864
# ╟─a8559c5e-744a-4c25-a75c-5394558bc672
# ╟─0b8bccd7-28ad-4ebc-8bcc-71e6adda71e3
# ╠═9230561b-26f3-4aa0-9c54-2f15aa231d80
# ╟─c74de7e2-b1b0-4f7f-9ac8-31d216caf987
# ╠═128af99c-4c34-4845-9fe8-876f5a003412
# ╟─0d040c4e-f306-4a84-93ac-9c35bc924130
# ╠═fd7a0025-1b77-45ab-bcbb-df39e38b3ffd
# ╟─5471b1c1-224c-4464-a4a3-7719347048fc
# ╠═2a1a8730-795b-4659-b688-6e320ffa81f8
# ╠═a067a699-5cfe-405e-bae2-4e1033f9f9c9
# ╟─ececcd18-56a3-41f7-8b15-bdff974080ac
# ╟─06716c35-b551-40a6-bb14-906bc7569434
# ╠═ef0babcf-d7a8-4fe3-9d64-a54504d8da77
# ╟─6e6fce09-d14c-459d-9dd6-5f8fd1ee8248
# ╠═d6a1c5e7-0921-4c4e-bc4a-fb2ebc94e45c
# ╟─8e31aaec-2512-4ad9-88a4-3130734d6e15
# ╟─d1c0d1ea-43c9-442e-899c-e83155f14841
# ╟─1bcd8f38-5ae9-459e-8f91-8c9ad7cde407
# ╠═9692c144-531a-4995-9057-60af2b91ecfa
# ╟─4da6a534-6a8f-464d-bbd4-9a6eec942688
# ╠═995a73db-6289-423b-8b5c-0fca4240472e
# ╠═fa46c270-9eff-4fe4-9716-cf4cbc86d9c0
# ╟─5b6d6431-3b77-4b2e-ae19-69000602935b
# ╠═8f92f95a-9596-4226-8167-7db772b01ff4
# ╟─e66acbf7-a93b-491f-8b2a-9554abd33a6c
# ╠═2e1043ec-f822-4c68-b173-fdaac26502db
# ╟─097763de-9f4b-4fa1-a64d-3938c46607c7
# ╠═890f4d32-04ff-4931-bc3c-e69005221183
# ╟─57279ca8-5824-445a-a837-6a7312a366a9
# ╠═f84cedc8-07b2-4744-81c1-a7b98a86a664
# ╟─50d0ccfb-78e4-432d-8f0a-abc18ccddacc
# ╟─2baa3355-f51b-464e-a55a-3471fb7e0100
# ╠═036177c4-08dd-4c22-b3fb-ca0ee807b307
# ╟─e35fe3ba-ae7d-407c-b034-5430449d414b
# ╠═a5e8f3d1-246e-4b33-bceb-c682f3bb4fbc
# ╟─c9a3dffd-c418-48a0-ade1-4c6ac3df9a40
# ╠═6754911a-81dd-41ce-8b59-2a0743b324c0
# ╟─c11a03f1-4aa2-4b6d-abf7-3e8c94d50498
# ╠═a67a89e9-858d-4501-834b-18f1ffd2ad0e
# ╟─5e1903be-3fd7-40dd-be3f-c8327c3c6033
# ╠═9dcfe511-f4a9-49ae-8571-33b402ca8d64
# ╟─a25724a6-097a-48b5-a500-f4c5f2fe5388
# ╠═d09898da-dc12-4e0d-b548-72c8f32458e0
# ╟─24df3f6a-8ac3-492b-b4e7-80029af6918d
# ╠═29d974a4-a748-4818-9c94-ca0a7b71e423
# ╠═ae82d253-4937-4e72-b516-6f3da8f19df1
# ╠═b8045050-09e2-42ea-a7aa-c702ec323fd3
# ╠═f79f92bc-e50b-4e5c-a82e-555269fa80b5
# ╠═8f06e886-6ebb-4a16-a1c5-8123e09ef723
# ╠═dd9023ad-cf79-4540-9259-61a28efa1174
# ╠═0f8c8b5c-33e0-4952-8d89-f26606a8503e
# ╠═2abb0086-0375-4633-8c1f-467b995ddd33
# ╠═cee62fc9-1d60-4356-b170-d9928e7de807
# ╠═8bd7c41a-112c-4ba7-b69e-af8e8d621c26
# ╠═bd574872-321d-4315-9f0f-60c95933226f
# ╟─3e110a7e-1272-4c50-80f6-290f61844952
# ╠═1f34c503-6815-4639-a932-9fe2d7e1b4c2
# ╟─ef60d36e-02c7-44e1-8910-241f99a89a38
# ╟─cea0c83c-6c6d-45a9-815c-6a112c3b76a9
# ╠═4a0dc42d-31a8-49bc-8db0-5f77474d8785
# ╟─d0a852cb-0d68-4807-861f-5302a5aeecd4
# ╠═8353d1fd-75cc-41e0-866c-5234634219d5
# ╠═d6b61d9a-817e-4c09-9a83-fde7ca591d23
# ╟─fb87aa13-1fd6-4c13-8ba4-68caef5f775b
# ╟─5dc919c8-01e1-4b23-b406-492562c9e338
# ╟─74709cb2-0be8-4b24-aa9d-926ef059fe2d
# ╠═68396ec1-444b-4c56-966b-5c70a2208d34
# ╠═8b51fcc1-6287-43cd-a44c-26c7423e8d7a
# ╟─2429c458-7108-4f13-aad2-dd04a3f1c3cf
# ╠═ba2a91f4-7f7a-4393-90da-3d5848d379e1
# ╟─d7063c34-6418-4e0b-80df-ddd07566e578
# ╟─0033e298-a85f-4718-acfa-fe9fd5a67746
# ╠═84a76d06-bf4e-461b-904c-40614c5b9626
# ╠═39f1d6d1-5d53-4c92-b4f0-8abd55a4fc57
# ╟─3c384fb9-fa39-441c-81d6-7306df13657d
# ╠═1cb81418-f8cb-4b4d-9343-586b5950f0da
# ╟─60ca2879-ee47-4163-b014-4a0af90b449f
# ╟─191fb08e-fade-43cd-8d8f-2273b594da68
# ╟─f3fee115-f239-4578-81e0-a3153027056f
# ╟─ca3eca34-3fc6-4488-ae74-69e9453c05e1
# ╟─c0b02aa0-642c-4a6f-ac21-cdb8e27d9279
# ╠═48a19e05-ba2e-421a-b059-07d994bed73c
# ╟─e91b91ec-77b4-4c77-b2ce-0e022f54edb7
# ╟─44b90b25-bdcd-4f38-804d-53e466ec6374
# ╠═0f8b23ab-c797-4513-b70c-5827376f6093
# ╟─756a5787-6f64-4308-b879-da7676d39a8c
# ╟─01a66e18-bc56-4204-b4b8-9bf46155aec1
# ╟─4fac9508-bc29-447c-9188-a58da07abeb3
# ╟─e8de3cd9-05bb-439e-b6fe-39d9c82f0ba6
# ╟─97e1379a-b0bf-4c6f-8375-299e1f42899f
# ╟─f11b6f84-d0c4-4855-a5c6-31d837880407
# ╟─0e4e510c-9662-4962-aeea-2a31f6c027fe
# ╟─d326776b-9a89-4464-b720-92d66e93e83c
# ╟─ab532d5b-4396-4792-b5e5-232216d6e4ca
# ╟─f2e5abd8-a544-4c21-88fe-af0edec938ab
# ╟─3eb70fdb-457a-4b5b-94d8-13b062895cd3
# ╠═6cb42cac-faa0-4bc4-bb2c-54195f07fd9a
# ╟─9571044c-7d41-4dd0-8601-fd5f37f05f0c
# ╟─45df2cc2-44d7-4b39-872d-b3a49cd59dfe
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
