### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ b3cb31aa-c982-4454-8882-5b840c68df9b
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
"""

# ╔═╡ 419d11bf-2561-49ca-a6e7-40c8d8b88b24
md"""
- `nmax = ` $(@bind nmax Slider([5, 10, 12, 15]; default=10, show_value=true))
"""

# ╔═╡ be5d3f98-4c96-4e69-af91-fa2ae5f74af5
md"""
## A difficult example

Solving linear equations is a standard exercise in linear algebra
and you probably have already done it in previous courses.
However, you might also know that solving such linear systems is
not always equally easy.

Let us consider a family of innocent-looking matrices, which are famously ill-conditioned, the **Hilbert matrices**. Here we show the $(nmax) by $(nmax) case.
Feel free to increase the `nmax` to make the problem even more challenging:
"""

# ╔═╡ 011c25d5-0d60-4729-b200-cdaf3dc89faf
M_difficult = [ 1/(i+j) for i=1:nmax, j=1:nmax ]

# ╔═╡ 8f6bffb7-e07a-45cf-a616-de6b3f0739b3
md"We take the simple right-hand side of all ones:"

# ╔═╡ e45952b2-714b-427b-8979-98d11c330294
b_difficult = ones(nmax)

# ╔═╡ 52d8634e-116b-48e6-9673-2ee2a97423c0
md"""
And solve the system using `\`:
"""

# ╔═╡ 5fb995cc-b338-4608-a01c-fc0e84d9dfe9
x_difficult = @time M_difficult \ b_difficult

# ╔═╡ 4d9d19da-7a4f-49ba-9c6f-890dad00672c
md"""
Looks like a reasonable answer, **but is it ?**

Let's check against a computation using Julia's `BigFloat` number. This is usually between 10 and 100 times more expensive, so not useful for practical computations. But it will give us a reference answer to much higher precision.
"""

# ╔═╡ fec623d0-08d3-434c-85af-266abde46da1
begin
	M_big = [ 1/(big(i)+big(j)) for i=1:nmax, j=1:nmax ]
	x_big = @time M_big \ b_difficult
end

# ╔═╡ 6e2c0370-a2a9-450c-b973-14925754edf8
md"""Looking at the second entry we already see some **significant deviations in the standard `Float64` answer**, which actually get way worse as we increase `nmax`:
"""

# ╔═╡ fff3d0d9-f018-42ac-ac36-d12e11ce9362
x_big - x_difficult

# ╔═╡ a8559c5e-744a-4c25-a75c-5394558bc672
md"""
While for `nmax = 5` the answer still kind of ok, **the result of the standard `\`-operator become numerical garbage** from `nmax = 12` onwards.

The related questions we want to ask here are:
- How does Julia's `\`-operator work ?
- How can we quantify when a linear system is more
  challenging to solve than another ?
- Based on this: When can we trust the results we get ?

This is what we want to explore in this part of the course.
"""

# ╔═╡ 0b8bccd7-28ad-4ebc-8bcc-71e6adda71e3
md"""
## Motivation: The \ (backslash) operator

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
Ok, so this calls into Julia's linear algebra library. The code is [located here](https://github.com/JuliaLang/julia/blob/8561cc3d68d3551a5728b40f782c244834fd3348/stdlib/LinearAlgebra/src/generic.jl#L1118). Essentially it performs
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

# ╔═╡ a4aab2c0-b254-4e24-a833-95a47abd0b25
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
- The first row of this system simply states $4 x_1 = 8$,
  which is very easy to solve, namely $x_1 = 8 / 4 = 2$.
- The second row states $3x_1-x_2=5$.
  However, $x_1$ is already known and can be inserted
  to find $x_2 = -(5-3\cdot 2)=1$.
- Following the same idea the third row gives
  $x_3=(0+1\cdot 2)/3 = 2/3$ and the last row
  $x_4=(1-1\cdot 2 + 1\cdot 1 + 1\cdot 2/3)/2 = 1/3$.
- In total we found
  ```math
  \mathbf{x} =
  \begin{pmatrix} 2 \\ 1 \\ 2/3 \\ 1/3
  \end{pmatrix}.
  ```
"""

# ╔═╡ 4396fbdd-12d7-4d73-bb98-948d20303292
md"""
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
        row_sum = 0.0
		for j in 1:i-1
			row_sum += L[i, j] * x[j]
		end
		
        x[i] = 1 / L[i, i] * (b[i] - row_sum)
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

# ╔═╡ d6a1c5e7-0921-4c4e-bc4a-fb2ebc94e45c
function backward_substitution(U, b)
	n = size(U, 1)
    x = zeros(n)
    x[n] = b[n] / U[n, n]
    for i in n-1:-1:1  # Note that this loop goes *backwards* n-1, n-2, ..., 1
        row_sum = 0.0
		for j in i+1:n
			row_sum += U[i, j] * x[j]
		end
        x[i] = (b[i] - row_sum) / U[i, i]
    end
    x
end

# ╔═╡ 7c66ecc0-f194-489a-bdbf-7b537f2d8567
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
       using **forward substitution** (Algorithm 1).
    3. Solve $\mathbf U \mathbf x = \mathbf z$ for $\mathbf x$
       using **backward substitution** (Algorithm 2).
"""

# ╔═╡ 4f017115-fd86-4087-8ce4-8ef896fa4959
md"""
A few remarks are in order:
- We have so far not discussed how to even compute the LU factorisation.
  As we will see in the next section this step will actually be accomplished
  by the **Gaussian elimination algorithm**, which probably already know
  from your linear algebra lectures.
- A severe **drawback of this naive algorithm** is immediately apparent from
  our 4x4 example problem, where (3) and (4) provide explicit solution
  expressions. If by any chance an **element $L_{ii}$ or $U_{ii}$ happens
  to be zero**, the **algorithm cannot work**. So one needs to come up with
  ways around this. We will mention a few ideas in the section on **pivoting**.
- Note, that **steps 2. and 3.** of the algorithm **do not depend on the
  original matrix $\mathbf A$**, but only on the right hand side $\mathbf b$.
  In other words once the factorisation $\mathbf{A} = \mathbf{L} \mathbf{U}$
  has been computed, solving $\mathbf A \mathbf x = \mathbf b$ for different
  right hand sides $\mathbf b$ only requires the execution of steps 2. and 3.
  As we will see later, step 1. is the most expensive step. Therefore
  **once $\mathbf{A}$ has been factorised**, the **cost for solving an
  *arbitrary* linear system** involving $\mathbf{A}$ **is reduced** as its
  factorised form $\mathbf{L} \mathbf{U}$ can be used in its place.
"""

# ╔═╡ 95448773-88f3-4e52-a2d3-5cd7e4a3e28f
md"""
## LU factorisation

The missing piece in our discussion is how the LU factorisation can be computed. As it turns out this is the Gaussian elimination algorithm, which you already learned in previous linear algebra classes. Indeed, this algorithm reduces a matrix to upper triangular form --- now we just need to be careful with the book-keeping to also extract the factor $\mathbf L$.
"""

# ╔═╡ 1bcd8f38-5ae9-459e-8f91-8c9ad7cde407
begin
	alg_LU = md"""
!!! info "Algorithm 4: LU factorisation"
    **Input:** $\textbf A \in \mathbb{R}^{n\times n}$ 

    **Output:** $\mathbf U \in \mathbb{R}^{n\times n}$, $\mathbf L \in \mathbb{R}^{n\times n}$
	-  $\textbf U = \textbf A$ 
    - for $k = 1, \ldots, n-1$ $\quad$  *(algorithm steps)*
      *  $L_{kk} = 1$
      * for $i = k+1, \ldots, n$ $\quad$   *(Loop over rows)*
        -   $L_{ik} = \frac{U_{ik}}{U_{kk}}$
        - for $j = k, \ldots n$ $\quad$   *(Loop over columns)*
          *   $U_{ij} \leftarrow U_{ij} - L_{ik} U_{kj}$  $\qquad$ *(modify $U_{ij}$ to $U_{ij} - L_{ik} U_{kj}$)*
	-  $L_{nn} = 1$  $\quad$ *(the loop above only runs until $n-1$)*
	"""

	alg_LU
end

# ╔═╡ c092eb61-b7f6-4868-8dc4-ca053c697a92
md"""
!!! warning "Example: Manual LU factorisation"
	We will run this algorithm manually for solving the linear system
    ```math
	\underbrace{\begin{pmatrix}
	2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4
	\end{pmatrix}}_{=\textbf A}
    \underbrace{\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix}}_{=\textbf x}
    = \underbrace{\begin{pmatrix}
	4 \\ 2 \\ -2
	\end{pmatrix}}_{= \textbf b}.
    ```
	That is for factorising the matrix $\mathbf A$. Before we start the loop over `k`, the matrix $\mathbf L$ is empty and $\mathbf U$ is just equal to $\mathbf A$:
	```math
	\mathbf L =   \begin{pmatrix} \phantom{-1} &\phantom{-1}&\phantom{1} \\ \phantom{-1} &  &  \\ \phantom{-1} &  &  \end{pmatrix} \qquad
	\mathbf U = \begin{pmatrix} \textcolor{red}{2} & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4\end{pmatrix}
	```
	- **`k=1` (Step 1):** In Gaussian elimination we would first zero out the **2nd and 3rd row** of the **1st column** of matrix $\mathbf A$. Here we do the same in a loop over rows starting at `k+1 = 2`.
	  * **`i = 2` (Row 2):** To zero out the first entry of the second row by subtraction, we need to multiply the first row with this factor:
	    -  $L_{21} = \frac{U_{21}}{\textcolor{red}{U_{11}}} = \frac{-4}{\textcolor{red}{2}} = \mathbf{-2}$
	    - The *loop over columns* now just uses this factor to update the second
	      row by subtracting $L_{21} = -2$ times the first --- or equivalently adding 2 times the first.
	    - After this step we have:
	      ```math
	      \mathbf L =   \begin{pmatrix} 1 &\phantom{-1}&\phantom{1} \\ \mathbf{-2} &  &  \\ \phantom{-1} &  &  \end{pmatrix} \qquad
	      \mathbf U = \begin{pmatrix} \textcolor{red}{2} & 1 & 0 \\ \textbf{0} & \textbf{5} & \textbf{-1} \\ 4 & -3 & 4\end{pmatrix}
	      ```
          where **bold** highlights the elements, that have changed.
	  * **`i = 3` (Row 3):** Here we zero out $U_{31}$. We determine the factor
	    -  $L_{31} = \frac{U_{31}}{\textcolor{red}{U_{11}}} = \frac{4}{\textcolor{red}{2}} = \mathbf{2}$
	    - The *loop over columns* updates the third row by subtracting $L_{31} = 2$ times the first row.
	    - After this step:
	      ```math
	      \mathbf L =   \begin{pmatrix} 1 &\phantom{-1}&\phantom{1} \\ -2 &  &  \\ \textbf{2} &  &  \end{pmatrix} \qquad
	      \mathbf U = \begin{pmatrix} \textcolor{red}{2} & 1 & 0 \\ 0 & 5 & -1 \\ \textbf{0} & \textbf{-5} & \textbf{4}\end{pmatrix}
	      ```
	    - In *loop over rows* `i` only runs until `n = 3`, so we are done with it.
	- **`k=2` (Step 2):** Our goal is now to zero out the 2nd column in all rows below the diagonal. We thus start another loop over rows, this time starting at `k+1 = 3`:
	  * **`i = 3` (Row 3):** Our factor is now
	    -  $L_{32} = \frac{U_{31}}{\textcolor{green}{U_{22}}} = \frac{-5}{\textcolor{green}{5}} = \mathbf{-1}$
	    - After subtracting $L_{32} = -1$ times the second row from the 3rd in the *loop over columns*, i.e. we add 2nd and 3rd row, we get
	      ```math
	      \mathbf L =   \begin{pmatrix} 1 &\phantom{-1}&\phantom{1} \\ -2 & 1 &  \\ 2 & \textbf{-1} &  \end{pmatrix} \qquad
	      \mathbf U = \begin{pmatrix} 2 & 1 & 0 \\ 0 & \textcolor{green}{5} & -1 \\ 0 & \textbf{0} & \textbf{3}\end{pmatrix}
	      ```
	  * Again we have reached the end of the *loop over rows* as `i = n = 3`.
	- Since `k = n-1 = 2` we also reached the end of the loop over `k`
    - Finally we set the missing $L_{33} = 1$ to obtain the final result:
	  ```math
	  \mathbf L =   \begin{pmatrix} 1 &\phantom{-1}&\phantom{1} \\ -2 & 1 &  \\ 2 & -1 & \textbf{1} \end{pmatrix} \qquad
	  \mathbf U = \begin{pmatrix} 2 & 1 & 0 \\ 0 & 5 & -1 \\ 0 & 0 & 3\end{pmatrix}
	  ```
    This is the LU factorisation of $\mathbf A$, which is easily verified by multiplying out the matrices:
    ```math
    \mathbf{L} \mathbf U = \begin{pmatrix} 1 &\textcolor{grey}{0}&\textcolor{grey}{0} \\ -2 &  1& \textcolor{grey}{0} \\ 2 & -1 & 1 \end{pmatrix}\begin{pmatrix}
	  2 & 1 & 0 \\ \textcolor{grey}{0}& 5 & -1 \\ \textcolor{grey}{0} &  \textcolor{grey}{0} & 3
	\end{pmatrix} = \begin{pmatrix}
	2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4
	\end{pmatrix} = \mathbf A
    ```
"""

# ╔═╡ 29bd39a3-9c0e-4038-a5ac-783fc9ac1629
Foldable("Expand for an explict Gaussian elimination procedure, which works directly working on the equations in x", 
md"""
!!! warning "Example: Gaussian elimination by working on the equations."
	We contrast the above algorithmic picture by a more naive approach working directly on the equations. In this case, too, the factors $L_{21}$, $L_{31}$ and $L_{32}$ appear.

	Let us run standard Gaussian eliminition manually for solving the linear system
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

    Notice that as before system matrix $\mathbf{A}$
    got slowly reduced to an upper triangular matrix,
    which we can identify with $\mathbf{U}$:
	```math
	\underbrace{\textbf A^{(1)}}_{=\textbf A}=\begin{pmatrix} 2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4\end{pmatrix}  \; \Rightarrow \;
	\textbf A^{(2)}=\begin{pmatrix} 2 & 1 & 0 \\ \textcolor{grey}{0} & 5 & -1 \\ \textcolor{grey}{0} & -5 & 4\end{pmatrix}  \; \Rightarrow \;
	\textbf A^{(3)}=\begin{pmatrix} 2 & 1 & 0 \\ \textcolor{grey}{0} & 5 & -1 \\ \textcolor{grey}{0} & \textcolor{grey}{0} & 3\end{pmatrix}=\textbf U
    ```
""")

# ╔═╡ af2e5735-b60d-4407-94fc-7ac3b14ac6a5
md"""
We see that indeed the LU factorisation algorithm can be seen as a formalisation of the Gaussian elimination procedure, reducing the matrix to triangular form $\mathbf{U}$.

Finally, for completeness, we provide a Julia implementation of LU factorisation (Algorithm 4):
"""

# ╔═╡ 9692c144-531a-4995-9057-60af2b91ecfa
function factorise_lu(A)
    n = size(A, 1)
    L = zeros(n, n)     # Initialise L and U by zeros
    U = float(copy(A))  # Make a copy of A and ensure that all entries
	                    # are converted to floating-point numbers

    for k in 1:n-1          # Algorithm steps
		L[k, k] = 1.0
		for i in k+1:n      # Loop over rows
			L[i, k] = U[i, k] / U[k, k]
			for j in k:n    # Loop over columns
				U[i, j] = U[i, j] - L[i, k] * U[k, j]  # Update U in-place
			end
		end
    end
	L[n, n] = 1.0           # Since the loop only runs until n-1
		
	# Return L and U using Julia datastructures to indicate
	# their special lower-triangular and upper-triangular form.
    return LowerTriangular(L), UpperTriangular(U)
end

# ╔═╡ ee45b963-fd91-4a7f-94b9-062bf45d5d7e
md"""We stay with the example where we performed manual Gaussian elimination:"""

# ╔═╡ cd1bbd16-ce88-49a4-bf1f-a9f42bc8920b
A_manual = Float64[ 2  1  0;
	               -4  3 -1;
                    4 -3  4]

# ╔═╡ 039b93ac-0a0e-45c3-acf9-70ec49b077c3
md"""With this slider you can stop `factorise_lu` after 1, 2 or 3 steps, checking that agrees with the steps we computed manually. The matrices displayed below show the situation after `k = nstep_lu_A` in **Algorithm 4** has finished.

- `nstep_lu_A = ` $(@bind nstep_lu_A Slider(0:3; default=0, show_value=true))
"""

# ╔═╡ c918c14f-1532-4f40-9857-127554aaaf42
begin
	function factorise_lu_steps(A; nstep=size(A, 1))
	    n = size(A, 1)
	    L = zeros(n, n)     # Initialise L and U by zeros
	    U = float(copy(A))  # Make a copy of A and ensure that all entries
	        	            # are converted to floating-point numbers
	
	    for k in 1:n-1          # Algorithm steps
			if k > nstep
				break
			end
			
			L[k, k] = 1.0
			for i in k+1:n      # Loop over rows
				L[i, k] = U[i, k] / U[k, k]
				for j in k:n   # Loop over columns
					U[i, j] = U[i, j] - L[i, k] * U[k, j]
				end
			end
	    end
		if nstep ≥ n
			L[n, n] = 1.0      # Since the loop only runs until n-1
		end
			
	    return (; U, L)
	end

	factorise_lu_steps(A_manual; nstep=nstep_lu_A)
end

# ╔═╡ 4da6a534-6a8f-464d-bbd4-9a6eec942688
md"""
### Running Algorithm 3

With this we have the missing ingredient to run Algorithm 3 and
numerically solve a linear system.

We stay with the problem
```math
\underbrace{\begin{pmatrix}
2 & 1 & 0 \\ -4 & 3 & -1 \\ 4 & -3 & 4
\end{pmatrix}}_{=\textbf A_\text{manuel}}
\underbrace{\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix}}_{=\textbf x}
= \underbrace{\begin{pmatrix}
4 \\ 2 \\ -2
\end{pmatrix}}_{= \textbf b_\text{manuel}}
```
which we solved manually beforehand:
"""

# ╔═╡ 995a73db-6289-423b-8b5c-0fca4240472e
A_manual

# ╔═╡ fa46c270-9eff-4fe4-9716-cf4cbc86d9c0
b_manual = [4, 2, -2]

# ╔═╡ 5b6d6431-3b77-4b2e-ae19-69000602935b
md"""
We follow **Algorithm 3**. First we need to find the LU factorisation:
"""

# ╔═╡ 8f92f95a-9596-4226-8167-7db772b01ff4
L, U = factorise_lu(A_manual)

# ╔═╡ e66acbf7-a93b-491f-8b2a-9554abd33a6c
md"""
which agrees what we obtained manually. Next we forward substitute:
"""

# ╔═╡ 2e1043ec-f822-4c68-b173-fdaac26502db
z = forward_substitution(L, b_manual)

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

# ╔═╡ fa504b66-4b10-4ba8-a535-895aa0b0c658
md"""
Hooray ! This suceeded!

Unfortunately **this simple algorithm does not always work**.

Consider the matrix:
"""

# ╔═╡ 87473cf7-68f6-4e34-b3b3-6780dc355506
D = [1 2 3;
     2 4 5;
	 7 8 9]

# ╔═╡ 75d4b513-4260-48dd-aa10-089ada83c093
md"""
If we apply our `factorise_lu` to this matrix we obtain:
"""

# ╔═╡ c1bb18ef-b58d-40a3-9d84-f16d6ccf6930
let
	L, U = factorise_lu(D)
	L
end

# ╔═╡ 8c0d6843-6bde-4c14-8a98-6cf7cb9f244a
md"""This `-Inf` is suspicious and points to a problem in the algorithm,
which we will investigate in the next section."""

# ╔═╡ 2baa3355-f51b-464e-a55a-3471fb7e0100
md"""
## LU factorisation with pivoting

We stay with the problem we identified at the end of the previous section. That is factorising the matrix

```math
\mathbf{D} = \begin{pmatrix} 1 & 2 & 3\\\textcolor{blue}{2} & \textcolor{blue}{4} & \textcolor{blue}{5}\\ \textcolor{orange}{7} &\textcolor{orange}{8} &\textcolor{orange}{9}\end{pmatrix}
```
"""

# ╔═╡ fbbff089-03ff-45ca-9ee7-2644d9fa8489
md"""
Applying **Algorithm 4** the first step ($k=1$) will zero out the first column in the second and third row by subtracting the 2 times (7 times) the first row from the second (third) row. When entering the loop for $k=2$ this results in:
```math
\text{at beginning of $k=2$:}\qquad
  \textbf U = \begin{pmatrix} 1 & 2 & 3\\0 & \textcolor{red}{0} & -1\\ 0 &-6 &-12\end{pmatrix}, \qquad
  \textbf L = \begin{pmatrix} 1 & & \\ 2 &\phantom{1} & \\ 7 & &\phantom{1} \end{pmatrix}.
```
The first step of the $k=2$ iteration in **Algorithm 4** will be to compute $L_{ik} = \frac{U_{ik}}{U_{kk}}$. However, the element $U_{kk}$ is zero as marked $\textcolor{red}{\text{in red}}$. We thus divide by zero, which is exactly what leads to the introduction of a `-Inf`
in the matrix `L`
"""

# ╔═╡ 8c4bc8a8-06bf-48ce-9d76-bf0c03f5a79d
md"""
To see this we use the slider to advance the algorithm step by step,
the matrices below show the situtation after the `k = nstep_lu_D` iteration has finished.

- `nstep_lu_D = ` $(@bind nstep_lu_D Slider(0:3; default=0, show_value=true))
"""

# ╔═╡ 716a0a8d-2bf4-4e3d-8044-09d542152cdc
factorise_lu_steps(D; nstep=nstep_lu_D)

# ╔═╡ 9d5aab1f-3156-4c73-9e2f-b91f3ebe3984
md"""
Due their central role in the Gaussian elimination algorithm
the denominators $U_{kk}$ in the computation of the $L_{ik}$
are usually referred to as **pivots**.

In summary we observe that from step $k \geq 2$ and onwards our computations are numerical garbage because we have used a $0$ pivot.
"""

# ╔═╡ 58c13b75-006d-48e4-8ddf-290df272488b
md"""
However, if instead of $\textbf D$ we factorise the **permuted matrix**
```math
\mathbf{P}\mathbf{D} = \begin{pmatrix} 1 & 2 & 3\\
\textcolor{orange}{7} &\textcolor{orange}{8} &\textcolor{orange}{9}\\
\textcolor{blue}{2} & \textcolor{blue}{4} & \textcolor{blue}{5}
\end{pmatrix}
```
which is the matrix $\textbf D$ in which the last two rows are swapped,
the algorithm goes through as expected:
"""

# ╔═╡ dca2b3b4-2b15-4428-91d0-2b949442b6bf
md"""
- `nstep_lu_PD = ` $(@bind nstep_lu_PD Slider(0:3; default=0, show_value=true))
"""

# ╔═╡ d194baba-556c-4669-a345-255c03081965
md"""
Let us also check that the LU factorisation indeed yields P * U:
"""

# ╔═╡ 49f66db0-3757-431d-85d7-4fe75096f9cc
md"We remark that the permutation matrix $\mathbf P$ is given by ..."

# ╔═╡ 6f59c78d-7966-45ab-8b79-8443b82ba1df
P = [1 0 0;
	 0 0 1;
	 0 1 0]

# ╔═╡ a85a3ad2-a89d-4e7c-8275-a2f3cd2921b0
P * D

# ╔═╡ f6cc916b-78d0-47d4-8b89-f4ab18926c1b
factorise_lu_steps(P * D; nstep=nstep_lu_PD)

# ╔═╡ 0ea63086-f090-4b5b-a7d5-d6a1b21b1fd8
let
	L, U = factorise_lu(P * D)
	L * U - P * D   # show that L * U = P * D
end

# ╔═╡ 2d461a32-c40f-4670-9252-09baa5f3a6d5
md"... and exactly achieves the task of swapping the last two rows, but leaving the rest of $\mathbf D$ intact as we saw above."

# ╔═╡ c42c3c63-b96c-4af0-a276-72ebd96fe23c
md"""
We notice that even though $\mathbf D$ cannot be permuted it is possible to obtain an  LU factorisation $\mathbf{L} \mathbf{U} = \mathbf{P} \mathbf{D}$
if we additionally allow the freedom to cleverly permute the rows of $\mathbf D$.

This is in fact a general result:

!!! info "Theorem 1"
    Every non-singular matrix $\textbf A \in \mathbb{R}^{n\times n}$ admits a factorisation
    ```math
    \textbf P \textbf A = \textbf L \textbf U
    ```
    where $\textbf L$ is lower-triangular, $\textbf U$ upper-triangular
    and $\textbf P$ is a permutation matrix.
"""

# ╔═╡ d85ab4d1-730b-4a6e-bd4d-46e450261f64
md"""
That is to say, that while in general **factorising a matrix $\mathbf A$ can fail**,
**it always suceeds** for non-singular matrices if we **give ourselves the freedom to permute the rows** of $\mathbf A$.

Finding a suitable $\mathbf P$ can be achieved by a small
modification of **Algorithm 4**.
Essentially this modification boils down to **selecting a permutation
of rows** of the factorised matrix on the fly,
**such that the pivots** of the permuted matrix $\mathbf{P} \mathbf{A}$ **are always nonzero**.
The precise way how this is achieved is called **pivoting strategy**
and the details are beyond the scope of this course.
The interested reader can find some discussion in the optional
subsection below.

Note, that Julia's implementation of LU factorisation
(the `lu` Julia function) does indeed implement one such
pivoting strategies, such that it works flawlessly on the problematic matrix `D`:
"""

# ╔═╡ 0de87398-22d4-4ba6-8831-bcf30e67a6db
facD = lu(D);

# ╔═╡ d46ee923-0912-4e85-bc33-0ddc6abe2fa5
md"Notably beyond the factors $\mathbf L$ and $\mathbf U$"

# ╔═╡ 12a63b8a-360b-421b-bf4c-03f923f46cc4
facD.L

# ╔═╡ edc42947-af9f-48f9-817c-037d9244f07a
facD.U

# ╔═╡ 7681fb06-703a-4b30-8973-20d45ca50732
md"this factorisation object also contains the employed permutation matrix $\mathbf P$:"

# ╔═╡ 9cf9c710-c76a-48b2-8799-efaaf7f1a82b
facD.P

# ╔═╡ dbf95d60-0333-407c-bc89-97914fb08437
md"Such that as expected $\mathbf L * \mathbf U = \mathbf P * \mathbf D$:"

# ╔═╡ ff388011-a9ed-4d1c-9f13-0b70c72187ae
facD.L * facD.U - facD.P * D

# ╔═╡ adae52fe-4862-47f6-a8e4-70ae7cc563b1
md"However, we notice that Julia's pivoting strategy did end up with a different permutation than our example."

# ╔═╡ 169ee211-c734-44d8-a034-b342aa5393d3
md"""
## Solving linear systems based on LU factorisation

The result of Theorem 1 is clearly that a factorisation
$\textbf P \textbf A = \textbf L \textbf U$
can always be achieved if the linear system $\mathbf A \mathbf x = \mathbf b$
has a unique solution, that is that $\mathbf A$ is non-singular.
To employ pivoted LU factorisation to solve this linear system we note that
we can always multiply both left and right hand sides by $\mathbf P$,
therefore:

```math
\textbf A \textbf x = \textbf b
\quad \Longleftrightarrow \quad
\textbf P \textbf A \textbf x = \textbf P  \textbf b
\quad \Longleftrightarrow \quad
\textbf L \textbf U \textbf x = \textbf P  \textbf b
```

which leads to the following algorithm:

!!! note "Algorithm 5: Solving linear systems with pivoted LU factorisation"
    Given a matrix $\mathbf{A} \in \mathbb{R}^{n\times n}$,
    a right-hand side $\mathbf b \in \mathbb{R}^{n}$
    the linear system $\mathbf{A} \mathbf x = \mathbf b$
    can be solved for  $\mathbf x \in \mathbb{R}^{n}$ as follows:
    1. Factorise $\textbf P \mathbf{A} = \mathbf{L} \mathbf{U}$, that is obtain $\mathbf P$, $\mathbf L$ and $\mathbf U$ $\qquad$ *(The `lu` Julia function)*.
    2. Solve $\mathbf{L} \mathbf{z} = \textbf P \mathbf b$ for $\mathbf z$
       using *forward* substitution.
    3. Solve $\mathbf U \mathbf x = \mathbf z$ for $\mathbf x$
       using *backward* substitution.

When we use Julia's backslash `\`-operator, effectively this **Algorithm 5** is executed under the hood.
"""

# ╔═╡ 40f6f4c5-dc86-4834-a6c2-51852d87a3bb
md"""Let us understand this algorithm by executing the steps manually for solving the problem

```math
\mathbf D \, \mathbf x_D = \mathbf b
\qquad \text{with} \quad \mathbf{D} = \begin{pmatrix} 1 & 2 & 3\\\textcolor{blue}{2} & \textcolor{blue}{4} & \textcolor{blue}{5}\\ \textcolor{orange}{7} &\textcolor{orange}{8} &\textcolor{orange}{9}\end{pmatrix}
 \quad \text{and}\quad \mathbf b = \begin{pmatrix} 2 \\ 4 \\ -2 \end{pmatrix}
```
We already defined defined $\mathbf{D}$ and $\mathbf{b}$ before:
"""

# ╔═╡ 24b7c334-28d4-41d1-b8f2-23fbc80ec2f3
D

# ╔═╡ e1d08709-d4bc-452f-a766-e43d25406fc2
b

# ╔═╡ 61f16e01-f7a3-420b-940f-ce0c6be93a99
md"**Step 1:** Perform LU factorisation of $\mathbf D$:"

# ╔═╡ ce383027-483d-424a-84c3-e8102475c43e
begin
	fac = lu(D)
	(; fac.P, fac.L, fac.U)
end

# ╔═╡ dc00020e-2fa0-4816-8b13-d457aee00823
md"**Step 2:** Compute $\mathbf z$ by forward substitution of $\mathbf P \mathbf b$ wrt. $\mathbf L$"

# ╔═╡ 43f22f5c-1fab-459e-a67e-fe67f6310d86
zD = forward_substitution(fac.L, fac.P * b)

# ╔═╡ 22d2cc3c-b50a-4915-8ba4-d7f6221ebd61
md"**Step 3:** Backward-substitute $\mathbf z$ wrt. $\mathbf U$:"

# ╔═╡ e13c4c6e-b55c-4360-9ec4-501042490d59
xD = backward_substitution(fac.U, zD)

# ╔═╡ 38deeaeb-0596-48d5-b671-2d1fc2e73a27
md"""**Verify result:** This solves the problem $\mathbf{D} \mathbf{x}_D = \mathbf{b}$ as desired:"""

# ╔═╡ 25f18e3c-a00b-4d3a-9b5e-60d40eeb5596
D * xD - b

# ╔═╡ 48979215-458c-40f5-82e4-a4485390bca4
md"""
## Optional: More details on pivoting

In this section we provide some details of one common form of
pivoting from LU factorisation, namely **row pivoting**.

In this approach we allow ourselves 
some flexibility in **Algorithm 4** by 
allowing ourselves to change the order in the last few rows
and thus *choose* the pivot amongst the entries
$U_{ik}$ with $i \geq k$ in each step $k$.
The resulting change of row order will then define the permutation
matrix $\mathbf{P}$ we employed in our above discussion.
"""

# ╔═╡ 827a75c2-2e3c-4da4-8deb-6c8a999596f2
md"""
As a guiding example we again consider the
problem of factorising
```math
\mathbf{D} = \begin{pmatrix} 1 & 2 & 3\\\textcolor{blue}{2} & \textcolor{blue}{4} & \textcolor{blue}{5}\\ \textcolor{orange}{7} &\textcolor{orange}{8} &\textcolor{orange}{9}\end{pmatrix}
```
where we saw previously non-pivoted LU factorisation to fail.
Without any pivoting / row swapping after the first LU factorisation step
(i.e when $k=2$ will start) the situation is

```math
\text{at beginning of $k=2$:}\qquad
  \textbf U = \begin{pmatrix} 1 & 2 & 3\\0 & \textcolor{red}{0} & -1\\ \textcolor{orange}{0} &\textcolor{orange}{-6} &\textcolor{orange}{-12}\end{pmatrix}, \qquad
  \textbf L = \begin{pmatrix} 1 & & \\ \textcolor{blue}{2} &\phantom{1} & \\ \textcolor{orange}{7} & &\phantom{1} \end{pmatrix}.
```

Looking at this $\textbf U$ it seems very reasonable to just swap the second and the
third column in and thus move the $-6$ to become the new pivot.
For consistency we not only need to swap $\textbf U$, but also $\textbf L$.
This gives us
```math
	\text{after row swap:}\qquad
  \textbf{U} = \begin{pmatrix} 1 & 2 & 3\\\textcolor{orange}{0} & \textcolor{red}{-6} &\textcolor{orange}{-12} \\
0 & 0 & -1\end{pmatrix}, \qquad
  \textbf{L} = \begin{pmatrix} 1 & & \\ \textcolor{orange}{7} &\phantom{1} & \\ \textcolor{blue}{2} & &\phantom{1} \end{pmatrix}.
```
If we continue the algorithm now, the $-6$ sits in the position of the pivot
$\left(A^{(k)}\right)_{kk}$ and the division by zero is avoided.
"""

# ╔═╡ f4a02392-ca42-44ac-bd1c-13cb3d6fafa2
md"""
In fact in this particular case the matrix $\mathbf U$ is already in upper triangular form after step $k=2$, such that in the step $k=3$ nothing will change, in fact we will just get $L_{32} = 0$. After $k=3$ we thus obtain from our algorithm:
```math
\textbf U = \begin{pmatrix} 1 & 2 & 3\\\textcolor{grey}{0} & -6 &-12 \\
\textcolor{grey}{0} & \textcolor{grey}{0} & -1\end{pmatrix}
\qquad
\textbf L = \begin{pmatrix}1& \textcolor{grey}{0} &\textcolor{grey}{0} \\ 7 & 1 & \textcolor{grey}{0} \\ 2 & 0& 1 \end{pmatrix}.
```
Due to the additional row permutation we performed
multiplying out $\textbf L \textbf U$ will not yield $\textbf A$,
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
"""

# ╔═╡ 2708ff44-5f30-43b3-b89a-19b9d8954100
md"""
This is not surprising as with pivoting we expect to get
$\textbf L \textbf U = \textbf P \textbf A$, so our missing piece is to
find the permutation matrix $\mathbf P$.

Here the correct matrix is 
```math
\textbf P = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0\end{pmatrix},
```
as we saw before. We can now understand how this matrix has been constructed:
Namely if we take the identity matrix and perform exactly the same row permutations as during the pivoted LU factorisation.
That is we swap the second and third row.
With this matrix we can easily check that
```math
\textbf L \textbf U = \textbf P \textbf A.
```

If we thus extend our `factorise_lu` function to additionally
perform such pivoting permutations,
the Gaussian elimination algorithm would always terminate successfully.
"""

# ╔═╡ 59083552-006d-4376-9c8a-2f5b85c6bb44
md"""
But pivoting brings additional opportunities.
As it turns out numerical stability of LU factorisation
can improve if one permuts the rows of $\textbf U$
not only if a pivot $U_{kk}$ is zero,
but in fact during each iteration $k$, ensuring that the 
**pivot** is **as large as possible**.

In other words in the $k$-th step of LU factorisation
we always exchange row $k$ with the row $l$ where $l$ satisfies
```math
\left|\,U_{lk}\,\right| \leq \left|\,U_{ik}\,\right|
\qquad
\text{for all $i = k, \ldots n$}.
```
The appropriate swaps are tracked and returned as well.

In practical algorithms instead of returning a permutation matrix
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
    Ak = float(copy(A))  # Make a copy of A and ensure that all entries
	                     # are converted to floating-point numbers

    for k in 1:n-1             # Algorithm steps
		p[k] = argmax(abs.(Ak[:, k]))  # Find row with maximal pivot
		
		U[k, :] = Ak[p[k], :]  # Copy pivot row to U, use U now instead of Ak,
		                       # which is again updated in-place
		for i in 1:n           # Row loop: Note the full range as any row may
		                       #           be non-zero
			L[i, k] = Ak[i, k] / U[k, k]
			for j = 1:n        # Column loop: Again full range
				Ak[i, j] = Ak[i, j] - L[i, k] * U[k, j]
			end
		end
    end
	p[n] = argmax(abs.(Ak[:, n]))
	U[n, n] = Ak[p[n], n]
	L[:, n] = Ak[:, n] / U[n, n]

	# To simplify assembling L we so far kept the rows in the same order
	# as in A. To make the matrix upper triangular we also apply the column
	# permutation p before returning the results.
	(; L=LowerTriangular(L[p, :]), U=UpperTriangular(U), p=p)
end

# ╔═╡ 24df3f6a-8ac3-492b-b4e7-80029af6918d
md"""
We check our implementation of pivoted LU factorisation against Julia's `RowMaximum` pivoting strategy, which implements the same algorithm.
"""

# ╔═╡ f55927a7-93de-4795-a1a9-734eca0ee60f
reference = lu(D, RowMaximum())

# ╔═╡ b659f816-a443-4be1-8359-776105d6e01e
md"The resulting L, U and pivot vector values are:"

# ╔═╡ 6d0a9fcf-f8e9-48fb-8167-eadd17cd0082
fac.L

# ╔═╡ ae82d253-4937-4e72-b516-6f3da8f19df1
fac.U

# ╔═╡ b8045050-09e2-42ea-a7aa-c702ec323fd3
fac.p

# ╔═╡ e5bf55d6-2698-4124-bd1c-2ea1d7e98ebd
md"In contrast we obtain:"

# ╔═╡ 8f06e886-6ebb-4a16-a1c5-8123e09ef723
factorise_lu_pivot(D).L

# ╔═╡ 88000701-f480-4efd-93b8-50b60c9d5d95
factorise_lu_pivot(D).U

# ╔═╡ f79f92bc-e50b-4e5c-a82e-555269fa80b5
factorise_lu_pivot(D).p

# ╔═╡ 3e110a7e-1272-4c50-80f6-290f61844952
md"""
## Memory usage and fill-in

### Sparse matrices

The matrices and linear systems one encounters in physics and engineering applications not infrequently reach **huge sizes**. For example in the numerical solution of partial differential equations using finite difference or finite element methods the vector of unknowns frequently has a **dimension on the order of $n \simeq 10^6$**.

Taking $n = 10^6$ thus as an example this implies that matrices,
i.e. entities from $\mathbb{R}^{n\times n}$,
hold around **$10^{12}$ elements**.
In double precision, i.e. $64$-bit numbers, each entry requires $8$ bytes. As a result **storing all elements** of such a matrix explicitly in memory reqires
"""

# ╔═╡ d0771eb3-6b89-4282-adfa-31c9a525f13c
10^12 * 8

# ╔═╡ 4647f5ea-9c17-4fd0-8f50-a734d8a1f840
md"bytes or"

# ╔═╡ 1f34c503-6815-4639-a932-9fe2d7e1b4c2
10^12 * 8 / 1024^4  # TiB

# ╔═╡ ef60d36e-02c7-44e1-8910-241f99a89a38
md"""
tibibyes of storage. **This is a challenge** even for modern compute clusters, where typically a node has around $512$ GiB. Laptops nowadays feature around $32$ GiB, so would be completely out of the question to solve such problems.

However, **in many applications** the arising **matrices are sparse**,
meaning that they contain a large number of zero elements,
which do not need to be stored. Let us discuss a few examples.
"""

# ╔═╡ 3ab0ccfe-df7c-4aed-b83f-f0113001610b
md"""
!!! info "Definition: Full matrix"
    A matrix $A \in \mathbb{R}^{n\times n}$ is called **full** if the number
    of non-zero elements is at the order of $n^2$. For such matrices almost all
    elements are non-zero and need to be stored.

Full matrices are the "standard case" and for them memory constraints usually set the limit of the size of problems, which can be tackled. For example an $32$ GiB memory laptop can store around $4 \cdot 10^9$ double-precision numbers, which means that
linear problems with more than around $n \simeq 60000$ unknows cannot be treated.
"""

# ╔═╡ 946fb6e6-e7e5-4aee-b566-b7c33ead7789
md"""
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
Using the `SparseArray` data structure from `SparseArrays` consistently allows to fully exploit sparsity when solving a problem. As a result the **storage costs scale only as $O(n)$**. With our laptop of 32 GiB memory we can thus tackle problems with around $n\simeq 4 \cdot 10^9$ unknowns --- much better than the $60000$ we found when using full matrices.

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
### Memory: LU factorisation of full matrices

Since the amount of available memory can put hard constraints on the size of linear systems which can be solved, we now want to investigate the memory requirement of LU factorisation $\textbf A = \textbf L \textbf U$.

If $\textbf A$ is a **full matrix**, then we know that $\textbf L$ and $\textbf U$ are triangular, thus **L and U** contain less non-zero elements or **more zero elements** than  $\textbf A$ itself and thus **require together as much memory** to be stored in memory as $\mathbf A$ itself.

### Memory: LU factorisation of sparse matrices

Let's **contrast** this **with the sparse matrices** we have considered above.
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
md"""non-zero elements, which is about **$(
round( (nnz(lu(Asp).L) + nnz(lu(Asp).U) ) / nnz(Asp), sigdigits=2)
) times as much** as the original matrix !"""

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
md"nonzeros, which is about $(
round( (nnz(lu(band).L) + nnz(lu(band).U) ) / nnz(band), sigdigits=2)
) times than the original in this case."

# ╔═╡ 191fb08e-fade-43cd-8d8f-2273b594da68
md"In both cases storing the $\mathbf L$ and the $\mathbf U$ factors require more non-zero elements than the original matrix. This phenomenon is usually referred to as **fill in**.

In particular for the sparse matrix `Asp` LU factorisation did not preserve the structure at all. In contrast it almost resulted in a full matrix in the lower right corner of the $\textbf L$. For $\mathbf U$ this looks exactly the same. In fact one can show that in general even for sparse matrices the **memory usage of LU factorisation** is $O(n^2)$. Therefore while one may be able to store a sparse matrix $\mathbf A$ with a huge size in memory (due to the $O(n)$ memory cost), one may actually run out of memory while computing the factorisation.
"

# ╔═╡ f3fee115-f239-4578-81e0-a3153027056f
md"""
Let us add that for banded matrices the situation is slightly better as one can show that in this case the fill-in takes place at most inside the band. As a result the memory requirement stays at $O(n d)$.
"""

# ╔═╡ e4aba0e4-adfc-4633-800e-3e2915300839
md"""
!!! info "Overview of LU factorisation memory cost"
    We summarise the cost of LU factorisation in a table:

	|  type of $n\times n$ matrix   | memory usage | comment |
	| ------------------------- |  ------------ | ------------ |
	|  full matrix  | $O(n^2)$             |
    |  general sparse matrix | $O(n)$ | fill-in
	|  banded matrix, band width $d$     | $O(n\,d)$     |  stays block-diagonal

"""

# ╔═╡ ca3eca34-3fc6-4488-ae74-69e9453c05e1
md"""
## Computational cost of LU factorisation

In this section we want to investigate how long it will take a computer to perform LU factorisation on a large matrix.

**Modern laptop computers** nowadays have a clock frequency of a few Gigahertz (GHz). This means they are **able to perform about $10^9$ operations per second**, where for simplicity we assume that "one operation" is an elementary addition, multiplication, division and so on. If we can determine how many such operations are needed to perform LU factorisation we can estimate the computational cost.

Before we look at the LU algorithm (Algorithm 4) we first understand a few simpler cases from linear algebra:
"""

# ╔═╡ c0b02aa0-642c-4a6f-ac21-cdb8e27d9279
md"""
### Scalar product
Given two vectors $\textbf x, \textbf y \in \mathbb{R}^n$ consider computing
the scalar product
```math
\textbf x \cdot \textbf y = \textbf x^T \textbf y = \sum_{i=1}^n x_i \, y_i
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
In the loop for each iteration we require $1$ multiplication and $1$ addition. In total the `scalar_product` function therefore requires $n$ multiplications and $n$ additions. The number of elementary operations is thus $2n$. We say the **computational cost is $O(n)$** (i.e. on the order of $n$ or linear in $n$).

If we take the dimensionality $n=1000$ the number of operations is $O(1000)$. On a 1 GHz computer (with $10^9$ operations per second) this therefore takes about $10^{-6}$ seconds ... hardly noticable.
"""

# ╔═╡ 44b90b25-bdcd-4f38-804d-53e466ec6374
md"""
### Matrix-vector product

Given a matrix $\textbf A \in \mathbb{R}^{n\times n}$ and a vector $\textbf x \in \mathbb{R}^n$
the matrix-vector product $\textbf y = \textbf A \textbf x$ is computed as
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
			y[i] = y[i] + A[i, j] * x[j]  # (*)
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
The loop over $i$ runs over $n$ values (the number of rows of $\textbf A$)
and the loop over $j$ over $n$ values as well (the number of columns of $\textbf A$).
In total the inner most instruction `(*)` is thus run $n^2$ times,
each times costing $2$ operations.
The total **cost is thus $O(n^2)$**.

For our case with $n = 1000$ and a 1 GHZ computer we thus now need
$(10^3)^2 / 10^{9} s = 10^{-3} s = 1\text{ms}$, which is again a rather short time.

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

# ╔═╡ 334d83a2-4809-4e85-808f-1213a56be9b7
md"""
!!! info "Overview of computational cost"
	For vectors $\textbf v, \textbf w \in \mathbb{R}^n$ and matrices $\mathbf A \in \mathbb R^{n\times n}$ the computational cost is

	| operation  | cost  |
	| ------------------------- |  ------------ |
	|  dot product $\mathbf v ^T \mathbf w$ | $O(n)$             |
    |  matriv-vector product $\mathbf A \mathbf v$ | $O(n^2)$
    |  matrix-matrix multiplication $\mathbf A \mathbf B$ | $O(n^3)$

"""

# ╔═╡ 3de7544d-a5e2-43fc-9af0-d2e37386b72a
md"""
Finally a few rough guidelines:

!!! info "General guideline to estimate computational cost"
	1. Determine the instructions at the innermost (most deeply nested) loop level
	2. Find the most expensive of these instructions and determine its cost in terms of *number of operations*
	3. Determine the ranges of all loops in terms of the dimensionality of your problem. Typically for each loop level this is $n$
	4. Multiply the result of 2. and all index ranges of 3 to get the total scaling. Typically for a single loop nesting the cost is $O(n)$ for a doubly nested loop $O(n^2)$ and so on.
"""

# ╔═╡ 7afabc20-a2a8-4b6d-87f8-09b84206a6dd
md"""
!!! warning "Example Matrix-matrix multiplication"
	Let us code up an algorithm how to compute the product of two matrices
	$\mathbf A, \mathbf B \in \mathbb{R}^{n\times n}$ and analyse its complexity.
"""

# ╔═╡ 1863e20c-ef35-4cea-8d3c-a33a0a42fc2e
function matmul(A, B)
	C = zeros(size(A, 1), size(B, 2))

	# loops ...
	
	C
end

# ╔═╡ 2d98da64-481b-4e86-9daa-f199294ed417
md"""
### LU factorisation

We revisit Algorithm 4:

$(alg_LU)
"""

# ╔═╡ 41ebce0a-1fde-4704-b889-b1a8391ccac7
md"""
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

For our case with $n = 1000$ and a 1 GHz computer the computation
thus needs approximately $(10^3)^3 / 10^{9} s = 1 s$,
which starts to be a noticable amount of time.
If we even consider $n = 10^5$ all of a sudden it takes $10^3 s$,
which is about $15$ minutes.
For **large matrices the cost of computing
LU factorisation can thus become important** !
"""

# ╔═╡ c34c30e8-b777-43b7-aae6-cd88d7179015
md"""
For **banded matrices** the computational scaling is a little improved.
E.g. consider LU factorisation on a banded matrix $A$ with band width $d$.
Then bothe the loop over $i$ (row loop) as well as the loop over $j$ (column loop) in the LU factorisation algorithm can be truncated and at most run over $d$ elements. As a result the computational cost at worst becomes $O(n d^2)$ which for a small band width (i.e. $d < n$) is substantially less than $O(n^3)$.
"""

# ╔═╡ 491d4fd8-1b22-4dac-9f2a-dcef9f3fb034
md"""
!!! info "Overview of LU factorisation computational cost and memory cost"
    We summarise the cost of LU factorisation in a table:

	|  type of $n\times n$ matrix   | computational cost | memory usage |
	| ------------------------- | ------------------ | ------------ |
	|  full matrix  | $O(n^3)$           | $O(n^2)$      |
    |  general sparse matrix
	|  banded matrix, band width $d$  | $O(n\,d^2)$    | $O(n\,d)$      |

"""

# ╔═╡ 51ab9105-5f8b-4370-867e-9042962fc777
md"""
## Numerical stability

To close off our discussion of Gaussian elimination we consider the question of its numerical stability, i.e. how is the **result of LU factorisation effected by small round-off errors** introduced when performing a calculation on a computer with its finite precision to represent numbers.
"""

# ╔═╡ fe34bed1-702b-49f3-b5f2-0c9961ca4e72
md"""
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
even **just providing this problem to a computer**
will usually be associated with **making a small error**:
in practice the computer only ever sees the **perturbed system**
```math
\tag{6}
\widetilde{\mathbf{A}} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}},
```
where $\widetilde{\mathbf{A}}$ is an approximation to $\mathbf{A}$,
$\widetilde{\mathbf{b}}$ is an approximation to $\mathbf{b}$.
It will thus yield the solution $\widetilde{\mathbf{x}}$
by solving the perturbed system instead of providing $\mathbf{x}$.

To understand the numerical stability when solving linear systems
our goal is thus to understand how far the solution
$\widetilde{\mathbf{x}}$ --- obtained on a computer ---
differs from $\mathbf{x}$ --- the solution to (5), the true linear system
we want to solve.
"""

# ╔═╡ 97e1379a-b0bf-4c6f-8375-299e1f42899f
md"""
To provide a motivation **why this type of analysis matters in practice**,
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
	can already have a very noticable effect on the solution and in fact
	render the solution we obtain from our our computation
	(i.e. $\widetilde{\mathbf{x}}$) rather far from the answer we are actually after
	(i.e. $\mathbf{x}$).
"""

# ╔═╡ fdad2eae-1aa5-42da-a276-3be34c598c19
md"""
In standard floating-point arithmetic the **relative error** in representing
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
However, to simplify the treatment in this course we will not attempt to discuss the effect of such relative errors on the numerical stability in full generality. Much rather we will employ the **error model**:
"""

# ╔═╡ 596b167c-6704-40bb-ac8a-99409466e831
md"""
!!! info "Error model in this notebook"
	Between the exact system and the perturbed system we will assume the
	foollowing relationship
	```math
	\begin{aligned}
	\widetilde{b}_i    &= b_i   \, (1 + \varepsilon) && \text{where $\varepsilon \approx 10^{-16}$} \\
	\widetilde{A}_{ij} &= A_{ij}
	\end{aligned}
	```

In other words we assume that the system matrix $\mathbf{A}$ is exactly represented
on the computer and that round-off errors only affect the right-hand side $\widetilde{\mathbf{b}}$. Moreover we assume that the relative error for all elements of $\mathbf b$ is identical to $\varepsilon$.

For standard double-precision floating-point numbers we have
$\epsilon \approx 10^{-16}$.

For an alternative and more detailed discussion see also
[chapter 2.8](https://tobydriscoll.net/fnc-julia/linsys/condition-number.html) of Driscoll, Brown: *Fundamentals of Numerical Computation*.
"""

# ╔═╡ f11b6f84-d0c4-4855-a5c6-31d837880407
md"""
### Matrix condition numbers and stability result

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

# ╔═╡ 25e0ebf6-ca45-4491-b5bc-c29f2725aa7b
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
"""

# ╔═╡ 679bfd40-ec6a-4fda-b2ab-b15d6acdbd32
md"""
From as stability analysis point of view the quantity $\kappa(\mathbf{A}) = \left\| \mathbf{A}^{-1}\right\| \, \|\mathbf{A}\|$ thus relates the relative error in the right-hand side (input quantity)
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

# ╔═╡ c4b82a47-1e66-4dfe-9160-ec3422420ce0
md"""

In this notebook we employed an error model
where storing $\mathbf{b}$ on the computer
(which leads to $\widetilde{\mathbf{b}}$)
introduces a small relative error:
```math
\widetilde{b}_i = b_i \, (1 + \varepsilon).
```
This enables to simplify the expression of the absolute error
of right-hand side:
```math
\left\| \mathbf{b} -  \widetilde{\mathbf{b}}  \right\|
= \sqrt{\sum_{i=1}^n b_i^2 \varepsilon^2}
= ε \, \sqrt{\sum_{i=1}^n b_i^2}
= \epsilon \, \| \mathbf{b} \|,
```
which leads to the following result:
"""

# ╔═╡ 630368a5-6934-4304-bbac-6995b819af8d
md"""
!!! info "Theorem 3: Stability of solving linear systems"
	Given a linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$
	and $\widetilde{\mathbf{x}}$ the solution to a perturbed
	linear system $\mathbf{A} \widetilde{\mathbf{x}} = \widetilde{\mathbf{b}}$,
	where $\widetilde{b}_i = b_i \, (1 + \varepsilon)$ for $i=1,\ldots,n$, then
	```math
	\frac{\left\| \mathbf{x} -  \widetilde{\mathbf{x}}  \right\|}
	{\| \mathbf{x}\|} \leq \kappa(\mathbf{A}) \, \varepsilon
	```

This inequality shows that the condition number of the matrix $\mathbf{A}$ plays the role of *amplifying* the round-off error introduced by the floating-point representation of the right-hand side.
"""

# ╔═╡ 73e6a875-69bd-4bef-84d6-3ade74e30015
md"""
Since for standard double-precision floating-point numbers we have
$\epsilon \approx 10^{-16}$, this means that condition numbers
of $10^{16}$ lead to a relative error of $1$, i.e. $100\%$ error
--- all precision is lost.
In general condition numbers above $10^8$ start to become problematic
as in this case the relative error in the solution is at most $10^{-8}$,
i.e. roughly speaking no more than $8$ digits.
"""

# ╔═╡ 3eb70fdb-457a-4b5b-94d8-13b062895cd3
md"""
**Numerically computing the condition number.**
In Julia code the condition number of a matrix is computed using the `cond` function. It effectively uses the above formulas for its computation.

For example for our running example we obtain:
"""

# ╔═╡ 6cb42cac-faa0-4bc4-bb2c-54195f07fd9a
let
	A = [1 1e-16;
	     1  0   ]
	κ = cond(A)
end

# ╔═╡ e9454cae-fe9a-4f35-ae53-2e723ea5a6b0
md"""Notice, that this is $2 \cdot 10^{16}$, which is a huge number !"""

# ╔═╡ 52280937-94b1-4259-93f8-d0d9214c27fc
md"""
!!! warning "Rounding error in a 2x2 system (continued)"
	Using Theorem 3 we can finally understand the behaviour we observe
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

	As we computed above
	```math
	\kappa(\mathbf{A}) \approx 2 \cdot 10^{16},
	```
	such that this large relative error is explained by the large
	condition number of this matrix.
"""

# ╔═╡ 86954ae7-a870-4695-8162-d9436cd19a33
md"""
### Computing condition numbers

Using the `cond` function of Julia enables us to easily compute condition numbers in practice.

However, the expression $\kappa(\mathbf{A}) = \left\| \mathbf{A}^{-1}\right\| \, \left\| \mathbf{A}\right\|$ is not extremly handy to provide a good intuition for what the condition number actually means or how it varies as the matrix $\mathbf A$ is changed.

In this section we thus discuss a few standard techniques how to compute condition numbers of matrices. These techniques are in fact exactly the algorithms that `cond` uses under the hood to do its computation.
"""

# ╔═╡ 8c327ba8-ea35-4acf-8774-a1f45126f040
md"""
First we need some notation:

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
"""

# ╔═╡ fab51c98-55c8-4bdc-a206-3456de4a324e
md"""
It turns out that the minimal and maximal eigenvalues provide
a convenient way to compute matrix norms:

!!! info "Lemma 4: Computing matrix norms"
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

# ╔═╡ 9be40396-9ace-4088-ac7b-8e0a0117e4f4
md"""
Based on this we can deduce a few useful formulas for computing condition numbers.
We start with the general expression:

!!! info "Corollary 5: Condition number for invertible matrices"
	If $\mathbf{A} \in \mathbb{R}^{n\times n}$ is invertible
	the condition number can be computed as
	```math
	\kappa(\mathbf{A}) = \left\| \mathbf{A}^{-1}\right\| \, \left\| \mathbf{A}\right\|
	= \frac{ \sqrt{ \lambda_\text{max}(\mathbf{A}^T\mathbf{A}) }}
	{ \sqrt{ \lambda_\text{min}(\mathbf{A}^T\mathbf{A}) }}.
	```
"""

# ╔═╡ 6f077373-cc53-4e9d-8849-c59c45b40178
md"""
Furthemore if $\mathbf{A} \in \mathbb{R}^{n\times n}$ is invertible and **symmetric**
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

If $\mathbf{A} \in \mathbb{R}^{n\times n}$ is **both symmetric
and positive definite** (i.e. all eigenvalues are strictly positive),
then we can drop the moduli in the above expression and obtain:

!!! info "Lemma 6: Condition number of symmetric, positive-definite matrices"
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
	As an example we now compute the condition number
	of the matrix analytically
	```math
	\mathbf{A} = \begin{pmatrix} 1 & 10^{-16} \\ 1 &0\end{pmatrix}.
	```
	Recall that the value computed by Julia's `cond` function was roughly $2 \cdot 10^{16}$.

	This matrix is not symmetric and we therefore use the expression
	```math
	\kappa(\mathbf{A}) = \frac { \sqrt{ \lambda_\text{max}(\mathbf{A}^T\mathbf{A}) }}  { \sqrt{ \lambda_\text{min}(\mathbf{A}^T\mathbf{A}) }},
	```
	which requires us to compute the eigenvalues of $\mathbf{A}^T\mathbf{A}$.
	We note
	```math
	\mathbf{A}^T\mathbf{A}
	= \begin{pmatrix} 1 & 1 \\ 10^{-16} &0\end{pmatrix}
	\begin{pmatrix} 1 & 10^{-16} \\ 1 &0\end{pmatrix}
	= \begin{pmatrix} 2 & 10^{-16}  \\ 10^{-16} & 10^{-32} \end{pmatrix}
	```
	and
	```math
	\det(\mathbf{A}^T\mathbf{A} - \lambda I)
	= \lambda^2 - (2 +  10^{-32}) \lambda +  10^{-32},
	```
	such that the two eigenvalues of $\mathbf{A}^T\mathbf{A}$ are
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
	--- i.e. the same huge number we obtained before.
"""

# ╔═╡ 45df2cc2-44d7-4b39-872d-b3a49cd59dfe
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 615)
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
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "a05f52688e4369db5594db7e90cf9dd4733214fc"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "a434e811d10e7cbf4f0674285542e697dca605d0"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.42"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

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
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "8252b5de1f81dc103eb0293523ddf917695adea1"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.1"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─ca2c949f-a6a0-485f-bd52-5dae3b050612
# ╠═3295f30c-c1f4-11ee-3901-4fb291e0e4cb
# ╟─21c9a859-f976-4a93-bae4-616122712a24
# ╟─b3cb31aa-c982-4454-8882-5b840c68df9b
# ╟─be5d3f98-4c96-4e69-af91-fa2ae5f74af5
# ╟─419d11bf-2561-49ca-a6e7-40c8d8b88b24
# ╠═011c25d5-0d60-4729-b200-cdaf3dc89faf
# ╟─8f6bffb7-e07a-45cf-a616-de6b3f0739b3
# ╠═e45952b2-714b-427b-8979-98d11c330294
# ╟─52d8634e-116b-48e6-9673-2ee2a97423c0
# ╠═5fb995cc-b338-4608-a01c-fc0e84d9dfe9
# ╟─4d9d19da-7a4f-49ba-9c6f-890dad00672c
# ╠═fec623d0-08d3-434c-85af-266abde46da1
# ╟─6e2c0370-a2a9-450c-b973-14925754edf8
# ╠═fff3d0d9-f018-42ac-ac36-d12e11ce9362
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
# ╟─a4aab2c0-b254-4e24-a833-95a47abd0b25
# ╟─4396fbdd-12d7-4d73-bb98-948d20303292
# ╠═ef0babcf-d7a8-4fe3-9d64-a54504d8da77
# ╟─6e6fce09-d14c-459d-9dd6-5f8fd1ee8248
# ╠═d6a1c5e7-0921-4c4e-bc4a-fb2ebc94e45c
# ╟─7c66ecc0-f194-489a-bdbf-7b537f2d8567
# ╟─4f017115-fd86-4087-8ce4-8ef896fa4959
# ╟─95448773-88f3-4e52-a2d3-5cd7e4a3e28f
# ╟─1bcd8f38-5ae9-459e-8f91-8c9ad7cde407
# ╟─c092eb61-b7f6-4868-8dc4-ca053c697a92
# ╟─29bd39a3-9c0e-4038-a5ac-783fc9ac1629
# ╟─af2e5735-b60d-4407-94fc-7ac3b14ac6a5
# ╠═9692c144-531a-4995-9057-60af2b91ecfa
# ╟─ee45b963-fd91-4a7f-94b9-062bf45d5d7e
# ╠═cd1bbd16-ce88-49a4-bf1f-a9f42bc8920b
# ╟─039b93ac-0a0e-45c3-acf9-70ec49b077c3
# ╟─c918c14f-1532-4f40-9857-127554aaaf42
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
# ╟─fa504b66-4b10-4ba8-a535-895aa0b0c658
# ╟─87473cf7-68f6-4e34-b3b3-6780dc355506
# ╟─75d4b513-4260-48dd-aa10-089ada83c093
# ╠═c1bb18ef-b58d-40a3-9d84-f16d6ccf6930
# ╟─8c0d6843-6bde-4c14-8a98-6cf7cb9f244a
# ╟─2baa3355-f51b-464e-a55a-3471fb7e0100
# ╟─fbbff089-03ff-45ca-9ee7-2644d9fa8489
# ╟─8c4bc8a8-06bf-48ce-9d76-bf0c03f5a79d
# ╟─716a0a8d-2bf4-4e3d-8044-09d542152cdc
# ╟─9d5aab1f-3156-4c73-9e2f-b91f3ebe3984
# ╟─58c13b75-006d-48e4-8ddf-290df272488b
# ╠═a85a3ad2-a89d-4e7c-8275-a2f3cd2921b0
# ╟─dca2b3b4-2b15-4428-91d0-2b949442b6bf
# ╟─f6cc916b-78d0-47d4-8b89-f4ab18926c1b
# ╟─d194baba-556c-4669-a345-255c03081965
# ╠═0ea63086-f090-4b5b-a7d5-d6a1b21b1fd8
# ╟─49f66db0-3757-431d-85d7-4fe75096f9cc
# ╠═6f59c78d-7966-45ab-8b79-8443b82ba1df
# ╟─2d461a32-c40f-4670-9252-09baa5f3a6d5
# ╟─c42c3c63-b96c-4af0-a276-72ebd96fe23c
# ╟─d85ab4d1-730b-4a6e-bd4d-46e450261f64
# ╠═0de87398-22d4-4ba6-8831-bcf30e67a6db
# ╟─d46ee923-0912-4e85-bc33-0ddc6abe2fa5
# ╠═12a63b8a-360b-421b-bf4c-03f923f46cc4
# ╠═edc42947-af9f-48f9-817c-037d9244f07a
# ╟─7681fb06-703a-4b30-8973-20d45ca50732
# ╠═9cf9c710-c76a-48b2-8799-efaaf7f1a82b
# ╟─dbf95d60-0333-407c-bc89-97914fb08437
# ╠═ff388011-a9ed-4d1c-9f13-0b70c72187ae
# ╟─adae52fe-4862-47f6-a8e4-70ae7cc563b1
# ╟─169ee211-c734-44d8-a034-b342aa5393d3
# ╟─40f6f4c5-dc86-4834-a6c2-51852d87a3bb
# ╠═24b7c334-28d4-41d1-b8f2-23fbc80ec2f3
# ╠═e1d08709-d4bc-452f-a766-e43d25406fc2
# ╟─61f16e01-f7a3-420b-940f-ce0c6be93a99
# ╠═ce383027-483d-424a-84c3-e8102475c43e
# ╟─dc00020e-2fa0-4816-8b13-d457aee00823
# ╠═43f22f5c-1fab-459e-a67e-fe67f6310d86
# ╟─22d2cc3c-b50a-4915-8ba4-d7f6221ebd61
# ╠═e13c4c6e-b55c-4360-9ec4-501042490d59
# ╟─38deeaeb-0596-48d5-b671-2d1fc2e73a27
# ╠═25f18e3c-a00b-4d3a-9b5e-60d40eeb5596
# ╟─48979215-458c-40f5-82e4-a4485390bca4
# ╟─827a75c2-2e3c-4da4-8deb-6c8a999596f2
# ╟─f4a02392-ca42-44ac-bd1c-13cb3d6fafa2
# ╟─2708ff44-5f30-43b3-b89a-19b9d8954100
# ╟─59083552-006d-4376-9c8a-2f5b85c6bb44
# ╠═6754911a-81dd-41ce-8b59-2a0743b324c0
# ╟─c11a03f1-4aa2-4b6d-abf7-3e8c94d50498
# ╠═a67a89e9-858d-4501-834b-18f1ffd2ad0e
# ╟─5e1903be-3fd7-40dd-be3f-c8327c3c6033
# ╠═9dcfe511-f4a9-49ae-8571-33b402ca8d64
# ╟─a25724a6-097a-48b5-a500-f4c5f2fe5388
# ╠═d09898da-dc12-4e0d-b548-72c8f32458e0
# ╟─24df3f6a-8ac3-492b-b4e7-80029af6918d
# ╠═f55927a7-93de-4795-a1a9-734eca0ee60f
# ╟─b659f816-a443-4be1-8359-776105d6e01e
# ╠═6d0a9fcf-f8e9-48fb-8167-eadd17cd0082
# ╠═ae82d253-4937-4e72-b516-6f3da8f19df1
# ╠═b8045050-09e2-42ea-a7aa-c702ec323fd3
# ╟─e5bf55d6-2698-4124-bd1c-2ea1d7e98ebd
# ╠═8f06e886-6ebb-4a16-a1c5-8123e09ef723
# ╠═88000701-f480-4efd-93b8-50b60c9d5d95
# ╠═f79f92bc-e50b-4e5c-a82e-555269fa80b5
# ╟─3e110a7e-1272-4c50-80f6-290f61844952
# ╠═d0771eb3-6b89-4282-adfa-31c9a525f13c
# ╟─4647f5ea-9c17-4fd0-8f50-a734d8a1f840
# ╠═1f34c503-6815-4639-a932-9fe2d7e1b4c2
# ╟─ef60d36e-02c7-44e1-8910-241f99a89a38
# ╟─3ab0ccfe-df7c-4aed-b83f-f0113001610b
# ╟─946fb6e6-e7e5-4aee-b566-b7c33ead7789
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
# ╟─e4aba0e4-adfc-4633-800e-3e2915300839
# ╟─ca3eca34-3fc6-4488-ae74-69e9453c05e1
# ╟─c0b02aa0-642c-4a6f-ac21-cdb8e27d9279
# ╠═48a19e05-ba2e-421a-b059-07d994bed73c
# ╟─e91b91ec-77b4-4c77-b2ce-0e022f54edb7
# ╟─44b90b25-bdcd-4f38-804d-53e466ec6374
# ╠═0f8b23ab-c797-4513-b70c-5827376f6093
# ╟─756a5787-6f64-4308-b879-da7676d39a8c
# ╟─01a66e18-bc56-4204-b4b8-9bf46155aec1
# ╟─334d83a2-4809-4e85-808f-1213a56be9b7
# ╟─3de7544d-a5e2-43fc-9af0-d2e37386b72a
# ╟─7afabc20-a2a8-4b6d-87f8-09b84206a6dd
# ╠═1863e20c-ef35-4cea-8d3c-a33a0a42fc2e
# ╟─2d98da64-481b-4e86-9daa-f199294ed417
# ╟─41ebce0a-1fde-4704-b889-b1a8391ccac7
# ╟─c34c30e8-b777-43b7-aae6-cd88d7179015
# ╟─491d4fd8-1b22-4dac-9f2a-dcef9f3fb034
# ╟─51ab9105-5f8b-4370-867e-9042962fc777
# ╟─fe34bed1-702b-49f3-b5f2-0c9961ca4e72
# ╟─97e1379a-b0bf-4c6f-8375-299e1f42899f
# ╟─fdad2eae-1aa5-42da-a276-3be34c598c19
# ╟─596b167c-6704-40bb-ac8a-99409466e831
# ╟─f11b6f84-d0c4-4855-a5c6-31d837880407
# ╟─25e0ebf6-ca45-4491-b5bc-c29f2725aa7b
# ╟─679bfd40-ec6a-4fda-b2ab-b15d6acdbd32
# ╟─c4b82a47-1e66-4dfe-9160-ec3422420ce0
# ╟─630368a5-6934-4304-bbac-6995b819af8d
# ╟─73e6a875-69bd-4bef-84d6-3ade74e30015
# ╟─3eb70fdb-457a-4b5b-94d8-13b062895cd3
# ╠═6cb42cac-faa0-4bc4-bb2c-54195f07fd9a
# ╟─e9454cae-fe9a-4f35-ae53-2e723ea5a6b0
# ╟─52280937-94b1-4259-93f8-d0d9214c27fc
# ╟─86954ae7-a870-4695-8162-d9436cd19a33
# ╟─8c327ba8-ea35-4acf-8774-a1f45126f040
# ╟─fab51c98-55c8-4bdc-a206-3456de4a324e
# ╟─9be40396-9ace-4088-ac7b-8e0a0117e4f4
# ╟─6f077373-cc53-4e9d-8849-c59c45b40178
# ╟─f2e5abd8-a544-4c21-88fe-af0edec938ab
# ╟─45df2cc2-44d7-4b39-872d-b3a49cd59dfe
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
