### A Pluto.jl notebook ###
# v0.20.21

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
	using Plots
end

# ╔═╡ ca2c949f-a6a0-485f-bd52-5dae3b050612
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/05_Direct_methods.pdf)
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
Ok, so this calls into Julia's linear algebra library. The code is [located here](https://github.com/JuliaLang/LinearAlgebra.jl/blob/master/src/generic.jl#L1226). Essentially it performs an **LU factorisation** for square matrices:
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

# ╔═╡ 85175a00-d0cf-46ca-98b0-14da532d9f52
Foldable("Teacher hint: Live coding",
md"""
```julia
function forward_substitution(L, b)
	n = size(L, 1)
    x = zeros(n)
    x[1] = # ...
    for i in 2:n
        row_sum = 0.0
		for j in 1:i-1
			row_sum += # ....
		end
        x[i] = # ...
    end
    x
end
```
""")

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

# ╔═╡ 458e6259-af57-4053-ac92-fe502b04f370
md"""
To make this more quantitative, we remind ourselves:

!!! info "Definition: Absolute and relative error"
	Given two vectors $\mathbf{x} \in \mathbb{R}^n$ and $\tilde{\mathbf{x}} \in \mathbb{R}^n$, where $\tilde{\mathbf{x}}$ is thought to be an approximation of $\mathbf{x}$. Then the absolute error of $\tilde{\mathbf{x}}$ is
	$\Vert\mathbf{x} - \tilde{\mathbf{x}}\Vert$, while the relative error is
	$\frac{\Vert\mathbf{x} - \tilde{\mathbf{x}}\Vert}{\Vert\mathbf{x}\Vert}$.

!!! warning "Rounding error in a 2x2 system (continued)"
	We saw above that a small perturbation
	$1 \to 1 + 10^{-16}$ in one of the elements of the right hand side,
	introduces a change in the solution from
	$(x_1, x_2) = (1, 0)$ to $(\widetilde{x}_1, \widetilde{x}_2) = (1, 1)$.
	In the solution this is an absolute error of
	```math
	\| \mathbf{x} -  \widetilde{\mathbf{x}} \| = \left\| \begin{pmatrix} 0 \\ -1  \end{pmatrix} \right\| = 1
	```
	and thus a relative error of 
	```math
	\frac{\left\| \mathbf{x} -  \widetilde{\mathbf{x}}  \right\|}
	{\| \mathbf{x}\|} = \frac{1}{1} = 1.
	```
	while in the right-hand side this only was an absolute error of
	```math
	\| \mathbf{b} -  \widetilde{\mathbf{b}} \| = 10^{-16}
	```
	and a relative error of
	```math
	\frac{\| \mathbf{b} -  \widetilde{\mathbf{b}} \|}{\| \mathbf{b}\|} = \frac{10^{-16}}{\sqrt{2}}
	```
"""

# ╔═╡ fdad2eae-1aa5-42da-a276-3be34c598c19
md"""
As outlined [in the introduction](https://teaching.matmat.org/numerical-analysis/01_Introduction.html), in standard floating-point arithmetic
the **relative error** in representing
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
on the computer and that round-off errors only affect the right-hand side $\widetilde{\mathbf{b}}$. Moreover we assume that the relative error for all elements of $\mathbf b$ is identical to $\varepsilon$, i.e. the overall relative error in the right-hand side is
```math
\frac{\|\mathbf b - \widetilde{\mathbf{b}}\|}{\|\mathbf b\|} =
\frac{|1 - (1+ε)| \, \|\mathbf b\|}{\|\mathbf b\|} = ε
```
as well.

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
We call this the condition number for a linear system:

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

This inequality shows that the condition number of the matrix $\mathbf{A}$ plays the role of *a worst-case amplification* of the round-off error introduced by the floating-point representation of the right-hand side.
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
	We also found this to be a relative error of
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

!!! info "Definition: Absolutely minimal and maximal eigenvalues"
	Given a diagonalisable matrix $\mathbf{M} \in \mathbb{R}^{n\times n}$
	we denote by $\lambda_\text{max}^\text{abs}(\mathbf{M})$ the *largest absolute* eigenvalue
	of $\mathbf{M}$ and by $\lambda_\text{min}^\text{abs}(\mathbf{M})$ the *smallest absolute* eigenvalue
	of $\mathbf{M}$, i.e.
	```math
	\lambda_\text{max}^\text{abs}(\mathbf{M}) = \max_{i=1,\ldots,n} |\lambda_i(\mathbf{M})|
	\qquad\text{and}\qquad
	\lambda_\text{min}^\text{abs}(\mathbf{M}) = \min_{i=1,\ldots,n} |\lambda_i(\mathbf{M})|
	```
	where $\lambda_i(\mathbf{M})$ for $i=1,\ldots,n$ are the eigenvalues of $\mathbf{M}$.
"""

# ╔═╡ d6311d6e-f774-4cdd-8ce9-3504d15ddc9f
md"""
!!! warning "Examples"
	To develop some intuition about $\lambda^\text{abs}_\text{max}$
	we consider the diagonal matrices:
	```math
	\begin{aligned}
	\mathbf{A} &= \begin{pmatrix} 5 & 0 \\ 0 & 10 \end{pmatrix}
	& \lambda_\text{max}^\text{abs}(\mathbf{A}) &= 10
	& \lambda_\text{min}^\text{abs}(\mathbf{A}) &= 5 \\
	%
	\mathbf{B} &= \begin{pmatrix} -5 & 0 \\ 0 & 10 \end{pmatrix}
	& \lambda_\text{max}^\text{abs}(\mathbf{B}) &= 10
	& \lambda_\text{min}^\text{abs}(\mathbf{B}) &= 5 \\
	%
	\mathbf{C} &= \begin{pmatrix} -500 & 0 \\ 0 & 10 \end{pmatrix}
	& \lambda_\text{max}^\text{abs}(\mathbf{C}) &= 500
	& \lambda_\text{min}^\text{abs}(\mathbf{C}) &= 10 \\
	\end{aligned}
	```	
"""

# ╔═╡ fab51c98-55c8-4bdc-a206-3456de4a324e
md"""
It turns out that the minimal and maximal absolute eigenvalues provide
a convenient way to compute matrix norms:

!!! info "Lemma 4: Computing matrix norms"
	For any matrix  $\mathbf{M} \in \mathbb{R}^{m \times n}$
	```math
	\|\mathbf{M}\| = \sqrt{ \lambda_\text{max}^\text{abs}(\mathbf{M}^T\mathbf{M}) }.
	```
	If $\mathbf{M}$ is moreover square and
	invertible, then
	```math
	\|\mathbf{M}^{-1}\| = \sqrt{ \lambda_\text{max}^\text{abs}(\mathbf{M}^{-T}\mathbf{M}^{-1}) } = \frac 1 {\sqrt{ \lambda_\text{min}^\text{abs}(\mathbf{M}^T\mathbf{M}) }}.
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
	= \frac{ \sqrt{ \lambda_\text{max}^\text{abs}(\mathbf{A}^T\mathbf{A}) }}
	{ \sqrt{ \lambda_\text{min}^\text{abs}(\mathbf{A}^T\mathbf{A}) }}.
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
= \frac{ \sqrt{ \lambda_\text{max}^\text{abs}(\mathbf{A}^T\mathbf{A}) }}
{ \sqrt{ \lambda_\text{min}^\text{abs}(\mathbf{A}^T\mathbf{A}) }}
= \frac{ \sqrt{ \max_{i=1,\ldots,n} |\lambda_{i}(\mathbf{A})|^2}}
   { \sqrt{\min_{i=1,\ldots,n} |\lambda_{i}(\mathbf{A})|^2}}
= \frac{ \max_{i=1,\ldots,n} |\lambda_{i}(\mathbf{A})|}
   { \min_{i=1,\ldots,n} |\lambda_{i}(\mathbf{A})|}
```

We obtain:

!!! info "Lemma 6: Condition number of symmetric matrices"
	If $\mathbf{A} \in \mathbb{R}^{n\times n}$ is symmetric, then
    ```math
	\kappa(\mathbf{A}) = \frac{\lambda_\text{max}^\text{abs}(\mathbf{A})}
     { \lambda_\text{min}^\text{abs}(\mathbf{A})}
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
	\kappa(\mathbf{A}) = \frac { \sqrt{ \lambda_\text{max}^\text{abs}(\mathbf{A}^T\mathbf{A}) }}  { \sqrt{ \lambda_\text{min}^\text{abs}(\mathbf{A}^T\mathbf{A}) }},
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
	\lambda_\text{max}^\text{abs}(\mathbf{A}^T\mathbf{A}) \approx 2
	\qquad\text{and}\qquad
	\lambda_\text{min}^\text{abs}(\mathbf{A}^T\mathbf{A}) \approx \frac{10^{-32}}{2}
	```
	therefore
	```math
	\kappa(\mathbf{A}) \approx 2 \cdot 10^{16},
	```
	--- i.e. the same huge number we obtained before.
"""

# ╔═╡ 86bc2d9e-903e-4069-be88-1b8cf1f737ff
md"""
## Lessons to learn about numerical stability

In the our previous discussion we noted that for a linear problem $\mathbf{A} \mathbf{x} = \mathbf{b}$ the condition number of the matrix $κ(\mathbf{A})$ provides the factor by which noise in $\mathbf{b}$ is at most amplified in the solution $\mathbf{x}$.

This concept of "condition number" or "conditioning" of a numerical problem
is in fact more general as defined below:

!!! info "Definition: Condition number of an algorithm"
	Consider an algorithm $f$, which relates an input quantity $\mathbf b$
	to an output quantity $\mathbf x$, i.e. $\mathbf{x} = f(\mathbf{b})$.
	The (relative) **condition number** $k$ of this algorithm is
	defined by the largest constant $K$ satisfying the relation
	```math
	\frac{\|f(\widetilde{\mathbf b}) - f(\mathbf b)\|}{\|f(\mathbf b)\|}
	≤ K \frac{\|\widetilde{\mathbf{b}} - \mathbf{b}\|}{\|\mathbf{b}\|}
	```
	for all valid inputs $\mathbf{b}$ and $\tilde{\mathbf{b}}$ to $f$.

In our case of the linear system we had that $f$ is the LU factorisation algorithm 3 to obtain the solution $\mathbf{x}$ of the linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$,
i.e. $\widetilde{\mathbf{x}} = f(\widetilde{\mathbf b})$ and ${\mathbf{x}} = f(\mathbf b)$. The condition number $K$ of this algorithm (see Theorem 2) is than just $\kappa(\mathbf{A})$, the condition number of the matrix $\mathbf{A}$.

!!! danger "Potential confusion: Two related condition number concepts"
	Note that the term *condition number* both refers to the *condition number of a matrix $\mathbf{A}$* and more generally to the *condition number of an algorithm*. While **for solving linear systems the two happen to coincide**, this is **not the case in general**.

Finally let us note:

!!! info "Definiton: Well-conditioned"
	We call an algorithm *well-conditioned* if its condition number is small (i.e. close to $1$). Else we call a problem badly conditioned.
"""

# ╔═╡ 66f6453b-9ad9-48ff-b5f5-60787c5cb381
md"""
## Overview of matrix factorisations

Algorithms to factorise matrices is a **fundamental building block** of numerical methods. We so far discussed the pivoted LU factorisation
$\mathbf{L} \mathbf{U} = \mathbf{P} \mathbf{A}$,
which applies to all square matrices $\mathbf{A}$.

You know at least one more matrix factorisation, namely the eigendecomposition for diagonalisable matrices
```math
\mathbf{U} \mathbf{\Lambda} \mathbf{U}^\dagger = \mathbf{A}
```
where $\mathbf{U}$ is the matrix of the eigenvectors of $\mathbf{A}$ as columns
and $\mathbf{\Lambda}$ is a diagonal matrix of eigenvalues.

In fact beyond these two factorisations, there are at least 8-10 more,
which you will encounter if you keep studying numerical methods for
scientific or engineering problems, see the documentation of Julia functions such as `qr`, `cholesky`, `svd`, `ldlt`, `bunchkaufman` if you are curious.

Each of these factorisations have there respective advantages and disadvantages.
Some are more computationally expensive to compute, some require additional properties of the matrix (e.g. positive-definite), explaining why there is such a zoo of methods to choose from.

Beyond LU and eigendecomposition we will only consider one additional
matrix factorisation method in this class,
namely **QR factorisation**.
This factorisation applies to **rectangular matrices** $\mathbf{A} \in \mathbb{R}^{n\times m}$ with $n > m$ and factorises this matrix as $\mathbf{Q} \mathbf{R} = \mathbf{A}$ with $\mathbf{Q}$ being a matrix of $m$ orthonormal vectors
and $\mathbf{R}$ being an $m \times m$ upper-triangular matrix.
We will discuss this method in detail in the context of solving [Regression and curve fitting problems](https://teaching.matmat.org/numerical-analysis/07_Interpolation.html#QR-factorisation),
where it is most commonly employed.

Comparing QR to LU factorisation the key difference is that LU only applies to square matrices, while QR can also be employed for rectangular matrices.
In fact coming back to our
question what Julia's `\` (backslash) operator does also make use of `qr`
(see the last line of [the implementation](https://github.com/JuliaLang/LinearAlgebra.jl/blob/master/src/generic.jl#L1226))
--- exactly for the case where $\mathbf{A}$ is not a square matrix.
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
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
HypertextLiteral = "~1.0.0"
Plots = "~1.41.4"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "83edc4557a1e92dbcab9548dbaecf1a63209c225"

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
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

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
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

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
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

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
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

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
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

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
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "f2fd0ae89599c11473fbd846ea1b9ed8c24613fb"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.21"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "55965bc822094d9e5d314d2845fb77f1b74fc8f7"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.21+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "6b4d2dc81736fe3980ff0e8879a9fc7c33c44ddf"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.2+0"

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
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

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
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

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
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

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
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39a11854f0cba27aa41efaedf43c77c5daa6be51"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.0+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
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
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "063ef757a1e0e15af77bbe92be92da672793fd4e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.4"

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

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

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
git-tree-sha1 = "34f7e5d2861083ec7596af8b8c092531facf2192"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "8f528b0851b5b7025032818eb5abbeb8a736f853"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

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
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

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
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

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
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

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
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

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
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
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
# ╟─85175a00-d0cf-46ca-98b0-14da532d9f52
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
# ╟─458e6259-af57-4053-ac92-fe502b04f370
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
# ╟─d6311d6e-f774-4cdd-8ce9-3504d15ddc9f
# ╟─fab51c98-55c8-4bdc-a206-3456de4a324e
# ╟─9be40396-9ace-4088-ac7b-8e0a0117e4f4
# ╟─6f077373-cc53-4e9d-8849-c59c45b40178
# ╟─f2e5abd8-a544-4c21-88fe-af0edec938ab
# ╟─86bc2d9e-903e-4069-be88-1b8cf1f737ff
# ╟─66f6453b-9ad9-48ff-b5f5-60787c5cb381
# ╟─45df2cc2-44d7-4b39-872d-b3a49cd59dfe
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
