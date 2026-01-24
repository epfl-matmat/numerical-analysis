### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ a8e3750d-62bf-4c1f-8309-3ef8981b5d26
begin
	using HypertextLiteral
	using PlutoUI
	using PlutoTeachingTools
end

# ╔═╡ c805c721-16ca-4f4f-ac76-c6382c1760b7
md"""
!!! info ""
    [Click here to view the PDF version.](https://teaching.matmat.org/numerical-analysis/03_Preliminaries.pdf)
"""

# ╔═╡ 2ff23584-4b35-11ef-0c69-a98bfb3ca906
md"""
# Revision and preliminaries


This chapter provides some revision on concepts you should already know from previous classes. We will need these throughout the course. It is strongly advised you read in carefully in your own time to remind yourself.
"""

# ╔═╡ 514be580-6ba5-45bc-9a8f-cfcbfaedf9df
TableOfContents()

# ╔═╡ 4c126a61-f2dd-42ca-bddc-6db64eb10dc5
md"""
## Taylor expansion and approximation

Given an infinitely differentiable function $f$ its **Taylor series** at the point **$a$** is the infinite sum
```math
\tag{1}
\begin{align*}
f(x) &= f(a) + f'(a) \cdot (x-a) +  \frac{1}{2!} f''(a) \cdot (x-a)^2
+ \frac{1}{3!} f^{(3)}(a) (x-a)^3 + \cdots \\
&= \sum_{n=0}^\infty \frac{1}{n!} f^{(n)}(a) (x-a)^n
\end{align*}
```
where $f^{(n)}$ is a short-hand for differentiating $f$ $n$ number of times.
"""

# ╔═╡ 316a1688-1315-4542-999d-a1dc4289a79f
md"""
### Big O notation

Truncating this sum early (i.e. not summing all the way to infinity)
is a frequent approximation technique often termed **Taylor approximation**.

Here we choose the example of a **second-degree approximation**,
i.e. truncating the sum at $n=2$ when we take at most two derivatives.
We would thus get an approximation
```math
f(x) \approx f(a) + (x-a) \  f'(a) + \frac12 f''(a) \  (x-a)^2.
```

Comparing with (1) we see
```math
\tag{2}
f(x) = f(a) + (x-a) \  f'(a) + \frac12 f''(a) \  (x-a)^2 + R(x)
```
where $R(x) = \sum_{n=3}^\infty \frac{1}{n!} f^{(n)}(a) (x-a)^n$
is usually called the **remainder term**.

As $x \to a$ the remainder $R(x)$ vanishes.
More precisely it becomes small *at least as fast as $(x-a)^3$*.
The idea is that as $x$ gets closer to $a$ then
$|x-a|$ < 1, such that $|x-a|^n < |x-a|^3$ for all $n > 3$.
While some of the derivatives $f^{(n)}(a)$ may be large,
there still is a point when $x$ is so close to $a$,
such that $|x-a| \ll 1$ and still
$f^{(n)}(a) |x-a|^n < |x-a|^3$ for all $n > 3$.

Mathematically one writes as:
> There exist positive constants $M, δ > 0$ such that
> ```math
> |R(x)| ≤ M |x - a|^3 \qquad \forall 0 < |x-a| < δ,
> ```

or more compactly using **big O notation** as

> ```math
> R(x) = O(|x - a|^3) \qquad \text{as $x \to a$}.
> ```

Employing this within (2) we arrive at

```math
\tag{3}
f(x) = f(a) + (x-a) \  f'(a) + \frac12 f''(a) \  (x-a)^2 + O(|x-a|^3)
```

This **big O notation** is thus a mathematically precise way of saying
that there are more terms in the expansion (3) that we don't show
but their order is at most $|x - a|^3$.
Or more generally

```math
\tag{4}
f(x) = \sum_{k=0}^n \frac{1}{k!} f^{(k)}(a) (x-a)^k + O(|x-a|^{n+1})
```
"""

# ╔═╡ 1405fdc1-614e-4f22-915e-1707e3cfec8b
md"""
### Lagrange form of the reminder

An important result of the study of Taylor approximations is, that the remainder term can in fact be also expressed in terms of the next term of the Taylor approximation itself.
Considering for example the general approximation (4)
we can find a $ξ$ between $x$ and $a$
--- i.e. assuming $x < a$ a number $ξ \in (x, a)$ ---
such that
```math
f(x) = \sum_{k=0}^n \frac{1}{k!} f^{(k)}(a) (x-a)^k + 
\textcolor{red}{
\frac{1}{(n+1)!} f^{(n+1)}(ξ) \ (x-a)^{n+1}
}
```
holds. This is the so-called **Lagrange form** of the remainder.
Note, how the red part is the next term of the Taylor polynomial,
just with the derivative evaluated at a (generally unknown) position $\xi$.

!!! warning "Second-degree approximation (continued)"
    To return to our example of a second-degree approximation
    and assuming $x < a$ we would thus be able to write
	```math
	f(x) = f(a) + (x-a) \cdot f'(a) + \frac12 f''(a) \cdot (x-a)^2 + \frac1{6} f'''(ξ) (x-a)^3
	```
	for some $ξ \in (x, a)$.
"""

# ╔═╡ 9ace76f9-1601-411d-ac70-dc20e3c8c55d
md"""
### Typical Taylor variants

In the above example we considered the expansion of $f(x)$ around a point $x = a$
leading to an expansion of the form
```math
f(x) = \sum_{k=0}^n \frac{1}{k!} f^{(k)}(a) (x-a)^k + O(|x-a|^{n+1}).
```
Another common setting is to expand a function $f(x+h)$ about $x$.
Replacing $x$ by $x+h$ in the above expression directly yields
```math
f(x+h) = \sum_{k=0}^n \frac{1}{k!} f^{(k)}(x) \, h^k + O(|h|^{n+1})
```
!!! warning "Second-degree approximation (continued)"
	Again in our second-degree example we would obtain
	```math
	f(x+h) = f(x) + h \ f'(x) + \frac12 f''(x) \  h^2 + O(|h|^3)
	```
	or using the Lagrange remainder
	```math
	f(x+h) = f(x) + h \ f'(x) + \frac12 f''(x) \  h^2 + \frac16 f'''(ξ) \ h^3
	```
	for some $ξ \in (x, x+h)$.
"""

# ╔═╡ ca35aec5-f9f1-44cd-acac-6b97454cefed
md"""
## Vector spaces

From a physical point of view the concept of a vector space $V$ boils down
to a well-defined set of objects, which are closed under linear combination.
That is to say for two vectors $\textbf{v} \in V$ and $\textbf{w} \in V$
and real scalars $α, β \in \mathbb{R}$ we want to have that any linear combination
is also a vector, i.e.
```math
α \textbf{v} + β \textbf{w} \in V.
```
Provided that ask the scalar multiplication $α \textbf{v}$
and the vector addition $\textbf{v} + \textbf{w}$ to satisfy
a few basic axioms, such as
-   $\textbf{v} + \textbf{w} = \textbf{w} + \textbf{v}$ ($+$ is commutative)
-   $α \textbf{v} + α \textbf{w} = α (\textbf{v} + \textbf{w})$
  and $α \textbf{v} + β \textbf{v} = (α+β) \textbf{v}$
  ($+$ and scalar multiplication are distributive)
- ... (see [Wikipedia](https://en.wikipedia.org/wiki/Vector_space) for the full list)
"""

# ╔═╡ 1dbdf498-01fa-4075-b6dd-907bcdc49352
md"""
we get a range of important properties, such as:

- Concept of linear independence,
  i.e. $n$ vectors $\textbf{v}_1, \textbf{v}_2, \ldots, \textbf{v}_n$
  are linearly independent if the only scalars $α_1, α_2, \ldots, α_n$
  to obtain
  ```math
  \textbf{0} = \sum_{i=1}^n α_i \textbf{v}_i
  ```
  are the zeros $α_1 =  α_2 =  \cdots =  α_n = 0$.

- Existance of a **basis**, i.e. a set of $n$ linearly independent vectors
  $\textbf{v}_1, \textbf{v}_2, \ldots, \textbf{v}_n$
  such that all elements of $V$ can be generated as linear combinations
  of these vectors, i.e. for all $\textbf{w} \in V$ we can find scalars
  $α_1, α_2, \ldots, α_n \in \mathbb{R}$ such that
  ```math
  \tag{5}
  \textbf{w} = \sum_{k=1}^n α_k \textbf{v}_k
  ```
"""

# ╔═╡ 61e6edb3-6ab3-4894-9447-4b2803f2f714
md"""
## Scalar product and norm

Frequently vector spaces feature an additional operation, the **inner product** or **scalar product**. For example for two Euclidean vectors $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$ the scalar product is
```math
\textbf x \cdot \textbf y = \textbf x^T \textbf y = \sum_{i=1}^n x_i \, y_i
```
Notable are the properties for all $\mathbf{x}, \mathbf{y}, \mathbf{z} \in \mathbb{R}^n$ and $α, β \in \mathbb{R}$:

-  $\textbf x \cdot \textbf y = \textbf y \cdot \textbf x$   *(Symmetry)*
-  *Linearity:*
```math
 ( \alpha \textbf x  + \beta \textbf y ) \cdot \textbf z  = α \textbf{x}\cdot \textbf{z} + β \textbf{y} \cdot \textbf{z} 
```
-   $\textbf x \cdot \textbf x \geq 0$ *(Non-negativity)*

Derived from the scaler product is the **vector norm**, measuring the length of a vector:
```math
\|\textbf x\| = \sqrt{\textbf x \cdot \textbf x} = \sqrt{\textbf x^T \textbf x}
= \sqrt{ \sum_{i=1}^n x_i^2 }
```
"""

# ╔═╡ 3579e3fb-e540-4c9c-94c5-d336047258c5
md"""
### Euclidean vector spaces

In this lecture when talking about vectors $\textbf{v}$ without specifying any further details we usually refer to elements of the Euclidean vector space $\mathbb{R}^n$, that is elements
```math
	\textbf{v} = \left( \begin{array}{c} v_1 \\ \vdots \\ v_n \end{array}\right)
```
"""

# ╔═╡ bbc55b67-c0fd-48a5-9c87-3cf1a4d5fd07
md"""
### Polynomials as vector spaces

Any polynomial $p(x)$ of degree less than $n$ can be written as
```math
\tag{6}
p(x) = \sum_{k=1}^n a_k x^k
```
where $a_k \in \mathbb{R}$. Comparing with (5) immediately suggests the
set of all polynomials of degree less than $n$ to be a vector space.
This indeed can be shown to be the case.

- In comparing with (5) we notice the **monomials** $1 = x^0, x = x^1, x^2, x^3, \ldots, x^n$ to be a possible choice for a basis.
- Similar to Euclidean vector spaces this is not the only choice of basis and in fact many families of polynomials are known, which are frequently employed as basis functions (e.g. [Lagrange polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial), [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials), [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials), ...)
- One basis we will discuss in the context of [polynomial interpolation](https://teaching.matmat.org/numerical-analysis/07_Interpolation.html) are Lagrange polynomials, which have the form 
  ```math
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
"""

# ╔═╡ 8e287200-7bb7-4a6a-a20f-2fdfaf450e33
md"""
## Taylor expansions of multi-dimensional functions
### Scalar-valued functions
We consider a function $f : \mathbb{R}^n \to \mathbb{R}$.
It's **second-order Taylor expansion** around $\mathbf{a} \in \mathbb{R}^n$
can be written as
```math
f(\mathbf{x}) = f(\mathbf{a})
+ \left(\nabla f(\mathbf{a})\right)^T \ (\mathbf{x} - \mathbf{a})
+ \frac12 (\mathbf{x} - \mathbf{a})^T \ \mathbf{H}_f(\mathbf{a}) \ (\mathbf{x} - \mathbf{a})
+ O(\|\mathbf{x} - \mathbf{a}\|^3),
```
where we introduced the gradient of $f$ at $\mathbf{x} = \mathbf{a}$
```math
\nabla f(\mathbf{a}) = \left( \begin{array}{c}
\left.\frac{\partial f}{\partial x_1}\right|_{\mathbf{x} = \mathbf{a}} \\
\left.\frac{\partial f}{\partial x_2}\right|_{\mathbf{x} = \mathbf{a}} \\
\vdots\\
\left.\frac{\partial f}{\partial x_n}\right|_{\mathbf{x} = \mathbf{a}}
\end{array}\right) \in \mathbb{R}^n
```
which is a vector of first partial derivatives.
Similarly the Hessian of $f$ at $\mathbf{x} = \mathbf{a}$
is the matrix of all second partial derivatives
```math
\mathbf{H}_f(\mathbf{a}) = \left( \begin{array}{cccc}
\left.\frac{\partial^2 f}{\partial x_1 \partial x_1}\right|_{\mathbf{x} = \mathbf{a}} &
\left.\frac{\partial^2 f}{\partial x_1 \partial x_2}\right|_{\mathbf{x} = \mathbf{a}} &
\ldots &
\left.\frac{\partial^2 f}{\partial x_1 \partial x_n}\right|_{\mathbf{x} = \mathbf{a}} \\
%
\left.\frac{\partial^2 f}{\partial x_2 \partial x_1}\right|_{\mathbf{x} = \mathbf{a}} &
\left.\frac{\partial^2 f}{\partial x_2 \partial x_2}\right|_{\mathbf{x} = \mathbf{a}} &
\ldots &
\left.\frac{\partial^2 f}{\partial x_2 \partial x_n}\right|_{\mathbf{x} = \mathbf{a}} \\
%
\vdots & \vdots && \vdots \\
%
\left.\frac{\partial^2 f}{\partial x_n \partial x_1}\right|_{\mathbf{x} = \mathbf{a}} &
\left.\frac{\partial^2 f}{\partial x_n \partial x_2}\right|_{\mathbf{x} = \mathbf{a}} &
\ldots &
\left.\frac{\partial^2 f}{\partial x_n \partial x_n}\right|_{\mathbf{x} = \mathbf{a}}
\end{array}\right) \in \mathbb{R}^{n\times n}
```
"""

# ╔═╡ e363da18-0001-4870-adb5-ef7d73514a36
md"""
!!! warning "Example"
	We consider the function $f : \mathbb{R}^3 \to \mathbb{R}$ defined as
	```math
	f(\textbf{x}) = f(\, (x_1, x_2, x_3)^T\,) = x_1^3 + x_1 x_3^2 + x_2^2 x_3 + 3 x_3^2
	```
	and compute the second-order Taylor expansion at
	```math
	\textbf{a} = \left(\begin{array}{c} 1 \\ 0 \\ 1 \end{array}\right)
	```
	First we compute the gradient $\nabla f(\textbf{a})$
	```math
	\nabla f(\mathbf{a}) = \left( \begin{array}{c}
	\left.\frac{\partial f}{\partial x_1}\right|_{\mathbf{x} = \mathbf{a}} \\
	\left.\frac{\partial f}{\partial x_2}\right|_{\mathbf{x} = \mathbf{a}} \\
	\left.\frac{\partial f}{\partial x_3}\right|_{\mathbf{x} = \mathbf{a}}
	\end{array}\right)
	= \left( \begin{array}{c}
	3a_1^2 + a_3^2 \\
	2a_2a_3 \\
	2a_1 a_3 + a_2^2 + 6a_3
	\end{array}\right)
	= \left( \begin{array}{c}
	4 \\
	0 \\
	8
	\end{array}\right)
	```
	Next the Hessian matrix $\mathbf{H}_f(\textbf{a})$
	```math
	\begin{align*}
	\mathbf{H}_f(\mathbf{a}) &= \left( \begin{array}{ccc}
	\left.\frac{\partial^2 f}{\partial x_1 \partial x_1}\right|_{\mathbf{x} = \mathbf{a}} &
	\left.\frac{\partial^2 f}{\partial x_1 \partial x_2}\right|_{\mathbf{x} = \mathbf{a}} &
	\left.\frac{\partial^2 f}{\partial x_1 \partial x_3}\right|_{\mathbf{x} = \mathbf{a}} \\
	%
	\left.\frac{\partial^2 f}{\partial x_2 \partial x_1}\right|_{\mathbf{x} = \mathbf{a}} &
	\left.\frac{\partial^2 f}{\partial x_2 \partial x_2}\right|_{\mathbf{x} = \mathbf{a}} &
	\left.\frac{\partial^2 f}{\partial x_2 \partial x_3}\right|_{\mathbf{x} = \mathbf{a}} \\
	\left.\frac{\partial^2 f}{\partial x_3 \partial x_1}\right|_{\mathbf{x} = \mathbf{a}} &
	\left.\frac{\partial^2 f}{\partial x_3 \partial x_2}\right|_{\mathbf{x} = \mathbf{a}} &
	\left.\frac{\partial^2 f}{\partial x_3 \partial x_3}\right|_{\mathbf{x} = \mathbf{a}}
	\end{array}\right) \\[0.3em]
	&= \left( \begin{array}{ccc}
	6a_1 & 0 & 2 a_3 \\
	0    & 2a_3 & 2a_2 \\
	2 a_3 & 2a_2 & 6 \\
	\end{array}\right)
	= \left( \begin{array}{ccc}
	6 & 0 & 2 \\
	0  & 2 & 0 \\
	2  & 0 & 6 \\
	\end{array}\right)
	\end{align*}
	```
	Combining the results we get
	```math
	f(\mathbf{x}) = f(\textbf{a})
	+ \left( \begin{array}{c}
	4 \\
	0 \\
	8
	\end{array}\right)^T (\mathbf{x} - \mathbf{a})
	+ \frac12 (\mathbf{x} - \mathbf{a})^T \ \left( \begin{array}{ccc}
	6 & 0 & 2 \\
	0  & 2 & 0 \\
	2  & 0 & 6 \\
	\end{array}\right) \ (\mathbf{x} - \mathbf{a})
	+ O(\|\mathbf{x} - \mathbf{a}\|^3),
	```
	which simplifies (after some algebra) to
	```math
	f(\mathbf{x}) = 4
	- \left( \begin{array}{c}
	10 \\
	0 \\
	0
	\end{array}\right)^T \mathbf{x}
	%
	+ 
	\mathbf{x}^T \ \left( \begin{array}{ccc}
	3  & 0 & 1 \\
	0  & 1 & 0 \\
	1  & 0 & 3 \\
	\end{array}\right) \ \mathbf{x} 
	+ O(\|\mathbf{x} - \mathbf{a}\|^3),
	```
"""

# ╔═╡ bd53a113-a153-405e-9340-da78a61b775a
md"""
### Vector-valued function
Now we consider a vector-valued function $\mathbf{f} : \mathbb{R}^n \to \mathbb{R}^m$.
It's **first-order Taylor expansion** around $\mathbf{a} \in \mathbb{R}^n$
can be written as
```math
\mathbf{f}(\mathbf{x}) = \mathbf{f}(\mathbf{a})
+ \textbf{J}_\textbf{f}(\mathbf{a}) \ (\mathbf{x} - \mathbf{a})
+ O(\|\mathbf{x} - \mathbf{a}\|^2).
```
In this the **Jacobian matrix**
$\textbf{J}_\textbf{f}(\textbf{x}) \in \mathbb{R}^{m\times n}$
is the collection of all partial derivatives of $\textbf{f}$, i.e.
```math
\textbf{J}_\textbf{f}(\mathbf{a}) = \left(\begin{array}{cccc}
\left.\frac{\partial f_1}{\partial x_1}\right|_{\textbf{x} = \textbf{a}} &
\left.\frac{\partial f_1}{\partial x_2}\right|_{\textbf{x} = \textbf{a}} &
\ldots &
\left.\frac{\partial f_1}{\partial x_n}\right|_{\textbf{x} = \textbf{a}} \\

\left.\frac{\partial f_2}{\partial x_1}\right|_{\textbf{x} = \textbf{a}} &
\left.\frac{\partial f_2}{\partial x_2}\right|_{\textbf{x} = \textbf{a}} &
\ldots &
\left.\frac{\partial f_2}{\partial x_n}\right|_{\textbf{x} = \textbf{a}} \\

\vdots & & \ddots & \vdots\\

\left.\frac{\partial f_m}{\partial x_1}\right|_{\textbf{x} = \textbf{a}} &
\left.\frac{\partial f_m}{\partial x_2}\right|_{\textbf{x} = \textbf{a}} &
\ldots &
\left.\frac{\partial f_m}{\partial x_n}\right|_{\textbf{x} = \textbf{a}}
\end{array}\right) \in \mathbb{R}^{m\times n}.
```
"""

# ╔═╡ f8720ba3-9f6e-476c-adf2-bf66e07c3755
md"""
!!! warning "Example"
	We consider the function $\textbf{f} : \mathbb{R}^3 \to \mathbb{R}^2$
	defined as
	```math
	\left\{ \begin{aligned}
	f_1(x_1, x_2, x_3) &= -x_1 \cos(x_2) - 1\\
	f_2(x_1, x_2, x_3) &= x_1 \, x_2 + x_3\\
	\end{aligned} \right.
	```
	and compute its first-order Taylor expansion at 
	```math
	\textbf{a} = \left(\begin{array}{c} 1 \\ 0 \\ 1 \end{array}\right)
	```
	The Jacobian is
	```math
	\begin{align*}
	\textbf{J}_\textbf{f}(\mathbf{a}) &= \left(\begin{array}{ccc}
	\left.\frac{\partial f_1}{\partial x_1}\right|_{\textbf{x} = \textbf{a}} &
	\left.\frac{\partial f_1}{\partial x_2}\right|_{\textbf{x} = \textbf{a}} &
	\left.\frac{\partial f_1}{\partial x_3}\right|_{\textbf{x} = \textbf{a}} \\
	%
	\left.\frac{\partial f_2}{\partial x_1}\right|_{\textbf{x} = \textbf{a}} &
	\left.\frac{\partial f_2}{\partial x_2}\right|_{\textbf{x} = \textbf{a}} &
	\left.\frac{\partial f_2}{\partial x_3}\right|_{\textbf{x} = \textbf{a}} \\
	\end{array}\right) \\
	&= \left(\begin{array}{ccc}
	-\cos(a_2) & a_1 \sin(a_2) & 0 \\
	a_2 & a_1 & 1
    \end{array}\right)\\
	&= \left(\begin{array}{ccc}
	-1 & 0 & 0 \\
	0 & 1 & 1
    \end{array}\right)
	\end{align*}
	```
	such that we arrive at
	```math
	\mathbf{f}(\mathbf{x}) = \mathbf{f}(\mathbf{a})
	+  \left(\begin{array}{ccc}
	-1 & 0 & 0 \\
	0 & 1 & 1
    \end{array}\right) \ (\mathbf{x} - \mathbf{a})
	+ O(\|\mathbf{x} - \mathbf{a}\|^2)
	```
	Again some algebra this is simplified to
	```math
	\mathbf{f}(\mathbf{x}) = \left(\begin{array}{c} -3 \\ 2 \end{array}\right)
	+  \left(\begin{array}{ccc}
		-1 & 0 & 0 \\
		0 & 1 & 1
	    \end{array}\right)\mathbf{x} + O(\|\mathbf{x} - \mathbf{a}\|^2) \\
	```
"""

# ╔═╡ c9b1c58c-b5a9-4780-8403-bed1ab74f0db
md"""
## Triangle inequality

An important inequality between numbers and vectors is the **triangle inequality**.
For example between two vectors $\textbf x$, $\textbf y$ it reads:
```math
\| \textbf x + \textbf y \| ≤ \| \textbf x \| + \| \textbf y \|.
```
Similarly between real numbers $α_1, α_2 \in \mathbb{R}$:
```math
| α_1 + α_2 | ≤ |α_1| + |α_2|.
```
Notably, this generalises to arbitrarily long sums:
```math
\left|  \sum_{i=1}^n α_i \right| ≤ \sum_{i=1}^n |α_i|
```
"""

# ╔═╡ 831aef2e-e129-4183-95c2-b6f98afb6b44
let
	RobustLocalResource("https://teaching.matmat.org/numerical-analysis/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(Markdown.parse(read("sidebar.md", String)), 435)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~1.0.0"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "ab9eab837fe12cbc0a139e7c843824a816db4625"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

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

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

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

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

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

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

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

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─c805c721-16ca-4f4f-ac76-c6382c1760b7
# ╟─2ff23584-4b35-11ef-0c69-a98bfb3ca906
# ╟─a8e3750d-62bf-4c1f-8309-3ef8981b5d26
# ╟─514be580-6ba5-45bc-9a8f-cfcbfaedf9df
# ╟─4c126a61-f2dd-42ca-bddc-6db64eb10dc5
# ╟─316a1688-1315-4542-999d-a1dc4289a79f
# ╟─1405fdc1-614e-4f22-915e-1707e3cfec8b
# ╟─9ace76f9-1601-411d-ac70-dc20e3c8c55d
# ╟─ca35aec5-f9f1-44cd-acac-6b97454cefed
# ╟─1dbdf498-01fa-4075-b6dd-907bcdc49352
# ╟─61e6edb3-6ab3-4894-9447-4b2803f2f714
# ╟─3579e3fb-e540-4c9c-94c5-d336047258c5
# ╟─bbc55b67-c0fd-48a5-9c87-3cf1a4d5fd07
# ╟─8e287200-7bb7-4a6a-a20f-2fdfaf450e33
# ╟─e363da18-0001-4870-adb5-ef7d73514a36
# ╟─bd53a113-a153-405e-9340-da78a61b775a
# ╟─f8720ba3-9f6e-476c-adf2-bf66e07c3755
# ╟─c9b1c58c-b5a9-4780-8403-bed1ab74f0db
# ╟─831aef2e-e129-4183-95c2-b6f98afb6b44
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
