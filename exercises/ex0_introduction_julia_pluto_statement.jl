### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 33490d1f-ef11-4159-823f-0040c096539c
# Install some packages
begin
	using LinearAlgebra
	using PlutoUI
	using PlutoTeachingTools
	using Plots
end

# ╔═╡ 0d3aec92-edeb-11ea-3adb-cd0dc17cbdab
md"# Exercice Session 0: Getting Started with Julia

This first session is meant for you to get acquainted with the [Julia programming language](https://julialang.org/). 

The following tutorial will give you an interactive tour of the main Julia functionalities we will need for the rest of the semester. 
By reading and running this notebook, you should be able to get a decent overview of Julia.

**References**: 
- When in doubt about syntax, have a look at the [Julia Cheatsheet](https://cheatsheet.juliadocs.org/)
- This introduction is inspired by the MIT lecture [Introduction to Computational Thinking](https://computationalthinking.mit.edu) where you can find more examples of Julia code and Pluto notebooks.
- Perhaps you also find the [Comparative cheatsheet Python <-> Julia <-> Matlab](https://cheatsheets.quantecon.org/) useful.
"


# ╔═╡ 4dbd5dd5-02d9-4243-a0fc-ac7089fc588a
md"""# Pluto 

[Pluto](https://plutojl.org/) is a browser-based notebook framework for Julia. 
All programming for this lecture will be done within Pluto notebooks. 

Pluto allows for running code in an interactive fashion while also presenting computation results in a readable way, also integrating formatted text.


### Cells
All Pluto inputs and outputs appear in **cells**. Cells can contain **text** or **code** and can be **visible** or **invisible**.

A cell can be run by **clicking the play button** below it or by hitting `Shift + Enter`. Its **visibility** can be toggled by clicking the *eye* button on its left (note that this only hides the cell's input, not its output).


#### Text cells
Text can be formatted by using [Markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax). 
Markdown cells begin with `\"md ` and end with `\"`.

!!! exercise "Exercise"
    **Modify** the cell below, **run** it and look at the output. Then **toggle its visibility**.
"""

# ╔═╡ 5e963571-a2bf-4554-9082-1ab5ec48d089
md" # Title
## Subtitle 

Text can be put in **bold** or *italic*.

- List element 1
- List element 2 

`inline quote` 

```julia
# code block 
using Plots
plot(sin)
```

[Test link](https://epfl.ch)
"

# ╔═╡ 462e68ff-5c86-4403-bcc5-907dbca15f2e
md"""
# Exercises
In the notebooks, there will be many interactive exercises, denoted by a green box such as this:

!!! exercise
    This is what an exercise looks like!

Many exercises are interactive:
- The statement is in a green box.
- There is a cell below where you have to complete the code.
- There is an immediate feedback box below.

Here is a first exercise to get you started:

!!! exercise
    Change the following line to `i_am_ready = true` then run the cell.
"""

# ╔═╡ 6424fb26-3b04-48d9-9c6f-0a1db75ef35f
i_am_ready = false

# ╔═╡ a249a972-8911-4c74-884d-5e30fd43dbec
if i_am_ready
	correct(text=md"""Good job! Now that you understand interactive exercises, let's continue.""")
else
	almost(md"""This is an interactive feedback cell.
	
	**You have not completed this exercise yet.** Change the line above to `i_am_ready = true` then run the cell.""")
end

# ╔═╡ 3b038ee0-edeb-11ea-0977-97cc30d1c6ff
md"# Julia Tutorial
**Now on to the Julia language itself.** The rest of this notebook demonstrates all the basic Julia functionalities we will need for the lecture. It is written in an interactive and conversational way, the goal being that by the time you have read and run the provided commands you will have a basic working understanding of how to use Julia.


**Please read and run all the cells below. Feel free to experiment by modifying the content of the cells.**
"

# ╔═╡ 6af9c3a0-2534-4bf1-9eb8-f5e34a091c64
md"""
## Expressions
Julia supports common math operations:
"""

# ╔═╡ 28c2ac78-2462-452e-bdc5-8319576fb33a
2 + 3

# ╔═╡ 0c8d79af-d762-4165-9992-688532185d48
2 - 3

# ╔═╡ 48ebe8f9-7ee1-478a-8ea3-f12f8b8b25e6
2 * 3

# ╔═╡ 9195da81-4465-4fc1-9831-b9aea1697bbe
2 / 3

# ╔═╡ 76f3caa3-4541-4268-aac6-6fb4d7aa6436
2 ^ 3

# ╔═╡ 5e062a24-edeb-11ea-256a-d938f77d7815
md"By default Julia displays the output of the last operation. (You can suppress the output by adding `;` (a semicolon) at the end.)
"

# ╔═╡ 166e8c0c-73d7-4a4d-b27f-13f9b4005dd2
md"""
Some more complex examples:
"""

# ╔═╡ 52f0d7e9-5eb2-42e6-b5ce-be8c6dd0ea4f
2 + 3 * 4

# ╔═╡ 025abeca-9f34-4921-9578-c571234c0616
(2 + 3) * 4

# ╔═╡ 4f9f0cc7-f6cf-443b-94af-bb31f213503e
2 + 3 * 10 ^ 2

# ╔═╡ f608bad0-43a6-4796-bfc8-7f30499797ad
2 + 3 * 10^2

# ╔═╡ b5143875-138b-413b-8073-1dc613555803
md"""
## Variables

We can define a variable using `=` (assignment). Then we can use its value in other expressions:
"""

# ╔═╡ 3e8e0ea0-edeb-11ea-22e0-c58f7c2168ce
x = 3

# ╔═╡ 59b66862-edeb-11ea-2d62-71dcc79dbfab
y = 2 * x

# ╔═╡ e938820b-daad-4cc9-9fd0-c8adc84adf46
md"""Pluto is **reactive**, meaning that when a variable changes, all other variables depending on it are automatically updated. 

!!! exercise "Exercise"
    Try modifying `z` in the below cell and run it. See how the value of `w` in the cell below is automagically updated.
"""

# ╔═╡ b6a3144d-eaba-4dff-802f-d20becff06e6
z = 10

# ╔═╡ 8ad688a5-2831-4017-ad39-2c98e5a0efd9
w = 2 * z

# ╔═╡ 0d0f4edf-a411-4a51-9c51-24256050b397
md"""
## A note on multiple expressions per cell
In Pluto, multiple expressions per cell are not allowed, as show by the error below:
"""

# ╔═╡ 2a4ec691-0944-4b3a-8def-a83eb8bf2400
var1 = 2
var2 = 3

# ╔═╡ 72e5b9d3-df98-4f6e-b283-52f928e7a80c
md"""
One solution is to split the code across multiple cells:
```julia
var1 = 2
```
```julia
var2 = 3
```

Another solution is to put the code inside a `begin`...`end` block:
```julia
begin
    var1 = 2
    var2 = 3
end
# var1 and var2 are available in other cells
```

Either way, the variables `a` and `b` will be available in future cells.

You might also encounter `let`...`end` blocks, which do not make the variables available to future cells:
```julia
let
    var1 = 2
    var2 = 3
end
# var1 and var2 are NOT available outside of let...end
```
"""

# ╔═╡ 44986e44-3371-477f-967e-9af708d103b3
md"""
!!! exercise
    In the following cell, assign `2` to the `var1` variable and `3` to the `var2` variable
"""

# ╔═╡ 35887694-5ad2-463a-b4a9-9f115542c9e2
# Set var1 to 2 and var2 to 3

# ╔═╡ 5e381128-b6f6-4594-8baa-f585d1e44061
#=╠═╡
if !isdefined(@__MODULE__, :var1)
	var_not_defined("var1")
elseif !isdefined(@__MODULE__, :var2)
	var_not_defined("var2")
elseif var1 == 2 && var2 == 3
	correct()
else
	almost(md"The variables exist but their values are not correct.")
end
  ╠═╡ =#

# ╔═╡ 1b2086c1-9b29-4b7d-8e88-7e8fbd86bde3
md"""
## Types
"""

# ╔═╡ 7e46f0e8-edeb-11ea-1092-4b5e8acd9ee0
md"""
In Julia, every value has a type that determines what can be done with it.
We can find the type of values using `typeof`:
"""

# ╔═╡ 5111713e-8c05-4161-b55a-814d72d32bc9
typeof(10)

# ╔═╡ 799938f0-8a4d-4505-bdcf-d5de7f70dcbb
typeof(10.0)

# ╔═╡ 57165a23-204d-4e68-8943-10bee7d68fa5
typeof(2.5)

# ╔═╡ 89336d7d-ebb8-46cc-92b2-dea09e3282c0
md"""
We can also ask for the type of a variable:
"""

# ╔═╡ 8a695b86-edeb-11ea-08cc-17263bec09df
typeof(y)

# ╔═╡ 769ed521-7faa-459a-b77e-15d5fbeef046
md"""
`Int64` means that this variable contains a signed 64-bit [integer](https://en.wikipedia.org/wiki/Integer_(computer_science)). In practice, this means that it is number without any decimal part.

`Float64` means that this variable contains a 64-bit [floating point number](https://en.wikipedia.org/wiki/Floating-point_arithmetic). In practice, this means a number that can have a decimal part. 
"""

# ╔═╡ 685174e3-21ca-46f7-bd6d-467718a6850e
md"""
As you learn Julia, you will discover more and more types. Sometimes, it can get quite complicated:
"""

# ╔═╡ 71bcb101-8715-4754-82ce-1d24fea8a512
typeof((; x=[[[[[1]]]]]))

# ╔═╡ d3287858-b154-42f6-9046-7030e6b8d722
md"""
## Strings
"""

# ╔═╡ 9072a34d-862f-4d6a-bbc4-7fce7a55c575
md"""
Another common type is `String`, which is used for text:
"""

# ╔═╡ 19c956cb-776b-4628-911d-a119f51c62d8
typeof("Some text")

# ╔═╡ cefeb394-9322-432d-8f25-6f70f1d15432
md"Variables can be inserted in a string using `$`:"

# ╔═╡ b9e1c787-3594-4107-b1f0-0293a666d565
let
	variable = 10
	"variable is $variable"
end

# ╔═╡ d3b94264-ff93-455e-a44c-3798cdacc660
md"To add a `$` to a string, escape it with `\` as follows:"

# ╔═╡ ae409c89-ec72-4d7f-8d39-79aa8946e533
println("0.02 \$")

# ╔═╡ 4a4dcb69-89b9-4622-b43f-d806bb27e861
md"To combine strings, use the `*` operator:"

# ╔═╡ 0af19cc8-d47c-4e01-afe9-d4f5d9225396
"part 1 " * "part 2 " * "part 3"

# ╔═╡ b7684e96-5a32-4e12-bbca-92ce7f57ce70
md"""
## Symbols
You may occasionally encounter symbols. They are typically used to represent identifiers and start with `:`.
"""

# ╔═╡ f7b4f385-78ed-4191-b2a9-fa4423b58f72
typeof(:name)

# ╔═╡ 273e451c-17fc-4816-831c-82f0c1f2f2ad
md"Keep in mind they are different from strings:"

# ╔═╡ 8ca29458-10dc-401a-b681-e63564552cf4
:name == "name"

# ╔═╡ 8e2dd3be-edeb-11ea-0703-354fb31c12f5
md"## Functions"

# ╔═╡ 82a0aae9-3937-4f60-9f75-532feb59e3f7
md"""
Julia comes with many built-in functions, for example `exp(x)` to compute $e^x$, `cos(x)` to compute $\cos(x)$, and many others...

Typing the function's name gives some basic information about the function:
"""

# ╔═╡ 80074898-5414-445b-8831-4c94fd8c9d52
exp

# ╔═╡ 1d945be1-3664-462c-ba4a-29e5e72724ea
cos

# ╔═╡ 9f28c628-dc42-4d4d-bf5d-d85a39cc1271
md"""
To call a function we must use parentheses:
"""

# ╔═╡ b8acc676-3e5b-4048-8120-ad6a182d2315
exp(0)

# ╔═╡ 2a3f7d44-0fb6-4929-9248-5265d0680619
exp(1)

# ╔═╡ dc29fcf1-ed6c-44a4-9e89-57d3f635a542
cos(0)

# ╔═╡ 7db7d2fd-4004-41bd-bfe9-b1ead83bdfea
cos(π)

# ╔═╡ 13a0268c-d815-4844-8b90-b2a475e9a4fa
md"""
### Creating new functions
"""

# ╔═╡ 96b5a28c-edeb-11ea-11c0-597615962f54
md"Of course, we can also write our own functions. For simple functions, we can use a short-form, one-line function definition:"

# ╔═╡ a7453572-edeb-11ea-1e27-9f710fd856a6
f(x) = 2 + x

# ╔═╡ b341db4e-edeb-11ea-078b-b71ac00089d7
md"As before, typing the function's name gives some basic information about the function."

# ╔═╡ 23f9afd4-eded-11ea-202a-9f0f1f91e5ad
f

# ╔═╡ d50fdb40-62eb-454d-baf3-42d1bb34f082
md"To call it we must use parentheses:"

# ╔═╡ cc1f6872-edeb-11ea-33e9-6976fd9b107a
f(10)

# ╔═╡ ce9667c2-edeb-11ea-2665-d789032abd11
md"For longer functions we use the following syntax with the `function` keyword and `end`:"

# ╔═╡ d73d3400-edeb-11ea-2dea-95e8c4a6563b
function g(x, y)
	z = x + y
	return z^2
end

# ╔═╡ e04ccf10-edeb-11ea-36d1-d11969e4b2f2
g(1, 2)

# ╔═╡ 2c446801-df67-4912-942c-aaa664302b28
md"""
**Note that the final `return` is not necessary.** In functions, the last expression is automatically returned:
"""

# ╔═╡ fc0bc6bd-b01c-461d-8cf8-ad3912ab2cbd
function h(x, y)
    z = x + y
	z ^ 2
end

# ╔═╡ b5e5cee9-dcf2-4e23-afce-72e136cbd988
h(1, 2)

# ╔═╡ d15f35df-23af-4793-bb4a-5878a4a52b23
md"""
!!! exercise
    Complete the function `line(a, b, x)` below, which should return $ax + b$.
"""

# ╔═╡ 971be2b7-d1e6-4e16-be24-a6a8d3092772
#=╠═╡
function line(a, b, x)
	nothing
end
  ╠═╡ =#

# ╔═╡ 7cfb56d4-b978-4e02-a64e-96431db95b61
#=╠═╡
let res = line(0.2885, 0.1026, 0.3452)
	if isnothing(res)
		still_nothing()
	elseif res == 0.2021902
		correct()
	else
		keep_working()
	end
end
  ╠═╡ =#

# ╔═╡ 93a231f4-edec-11ea-3b39-299b3be2da78
md"## Conditionals: `if`, `elseif`, `else`"

# ╔═╡ 82e63a24-eded-11ea-3887-15d6bfabea4b
md"""
We can evaluate whether a condition is true or not by using comparison operators:
- `<`: smaller than
- `<=`: smaller or equal
- `>`: greater than
- `>=`: greater or equal
- `==`: equal
- `!=`: not equal
"""

# ╔═╡ 9b339b2a-eded-11ea-10d7-8fc9a907c892
a = 3

# ╔═╡ 9535eb40-eded-11ea-1651-e33c9c23dbfb
a < 5

# ╔═╡ 73f81dc2-1138-4f44-be8c-fc0ee55aa2e2
a != 5

# ╔═╡ 58053efc-8656-4b9c-abfc-a15c29cb10b3
a >= 10

# ╔═╡ 530f923d-c820-4a45-8215-a9a57e375a18
typeof(a >= 10)

# ╔═╡ a16299a2-eded-11ea-2b56-93eb7a1010a7
md"We see that conditions have a boolean (`true` or `false`) value. The corresponding Julia type is `Bool`.

We can then use `if` to control what we do based on that value:"

# ╔═╡ bc6b124e-eded-11ea-0290-b3760cb81024
if a < 5
	"small"
else
	"big"
end

# ╔═╡ cfb21014-eded-11ea-1261-3bc30952a88e
md"""Note that the `if` also returns the last value that was evaluated, in this case the string `"small"` or `"big"`. Since Pluto is reactive, changing the definition of `a` above will automatically cause this to be reevaluated!"""

# ╔═╡ ca7089ae-0394-4d52-82cc-b893e2300894
md"""
Intermediate checks can be added using `elseif`:
"""

# ╔═╡ ee0510a9-95b3-4e43-addb-8ffdb724d3ad
if a < 0
	"negative"
elseif a == 0
	"zero"
elseif a < 10
	"1 to 9"
else
	"10 or larger"
end

# ╔═╡ b96761f5-b235-47b0-b3b4-f350feea46ac
md"""
## Conditionals: logical operators
Comparisons can be combined using logical operators:
- `a && b`: checks that `a` **and** `b` are `true`.
- `a || b`: checks that at least `a` **or** `b` is `true`.
- `!a`: checks that `a` is **not** `true`.

For example:
"""

# ╔═╡ 22bd1bc0-c2fa-489a-a2c9-b1c75b36f7ba
1 < 2 && 2 < 3

# ╔═╡ f2426c25-32c6-4408-80c3-d234a0fd6a65
1 < 2 || 2 < 1

# ╔═╡ b8201bdd-585a-4ac2-ae44-bb2c1b378dca
1 < 2 && 2 < 1

# ╔═╡ 18adab01-e74e-4537-8549-fadad03f30b7
!(1 == 2)

# ╔═╡ 77d9277c-54bc-4d80-acc7-d8d866753e27
md"""
!!! exercise
    Complete the function `is_leap_year(year)` below. The function should return `true` if `year` is a leap year, and `false` if it is not. Here is a reminder of the algorithm:
    - Every year divisible by 4 is a leap year, except that:
    - Every year divisible by 100 **is not** a leap year, except that:
    - Every year divisible by 400 **is** a leap year.

    To check if the year is divisible by some number `n`, check that `year % n == 0`. (`%` is called the [modulo operator](https://en.wikipedia.org/wiki/Modulo)).
"""

# ╔═╡ 54a7f75f-9835-497e-a1da-ecc36ac61641
#=╠═╡
function is_leap_year(year)
	nothing
end
  ╠═╡ =#

# ╔═╡ 4053a6a7-54a0-4b05-aa44-f6ffaf38cfe7
#=╠═╡
if isnothing(is_leap_year(2000))
	still_nothing()
elseif is_leap_year(133) != false || is_leap_year(2025) != false
	keep_working(md"Non-multiples of 4 are not correct yet.")
elseif is_leap_year(2004) != true || is_leap_year(2008) != true || is_leap_year(64) != true
	keep_working(md"Multiples of 4 are not correct yet.")
elseif is_leap_year(2100) != false || is_leap_year(1500) != false
	keep_working(md"Multiples of 100 are not correct yet.")
elseif is_leap_year(1600) != true || is_leap_year(2000) != true
	keep_working(md"Multiples of 400 are not correct yet.")
else
	correct()
end
  ╠═╡ =#

# ╔═╡ e297c5cc-edeb-11ea-3bdd-090f415685ab
md"## For loops"

# ╔═╡ ec751446-edeb-11ea-31ba-2372e7c71b42
md"Use `for` to loop through a pre-determined set of values:"

# ╔═╡ fe3fa290-edeb-11ea-121e-7114e5c573c1
let
    s = 0
	
	for i in 1:10
		s += i    # Equivalent to s = s + i
	end
	
	s
end

# ╔═╡ 394b0ec8-eded-11ea-31fb-27392068ef8f
md"Here, `1:10` is a **range** representing the numbers from 1 to 10:"

# ╔═╡ 4dc00908-eded-11ea-25c5-0f7b2b7e18f9
typeof(1:10)

# ╔═╡ 6c44abb4-edec-11ea-16bd-557800b5f9d2
md"Above we used a `let` block to define a new local variable `s`. 
But blocks of code like this are usually better inside functions, so that they can be reused. For example, we could rewrite the above as follows:
"

# ╔═╡ 683af3e2-eded-11ea-25a5-0d90bf099d98
function mysum(n)
	s = 0
	
	for i in 1:n
		s += i    
	end
	
	return s
end

# ╔═╡ 76764ea2-eded-11ea-1aa6-296f3421de1c
mysum(100)

# ╔═╡ bd39bd6d-152a-48d2-bc7b-43f426eb1165
md"""
`for` loops work over many other things than ranges. Here is an example of looping over a vector (see below). `$element` includes the value of the `element` variable in a string.
"""

# ╔═╡ ddcac2dc-f5b4-47fe-81f2-262b1b4718b2
for element in [1, 10, 100, 1000]
	println("Looping over $element")
end

# ╔═╡ 1b49f71c-a8b4-440b-a46a-17dfb4220fa5
md"""
## Exercise: Fibonacci
Before we move on, let us implement a classic exercise.
!!! exercise
    Complete the below cell by writing a function `fibonacci` that accepts an integer `` n `` and returns the `` n ``-th term of the [Fibonacci sequence](https://en.wikipedia.org/wiki/Fibonacci_sequence) as an integer,
    that is the sequence $F_n$ with
    ```math
    \begin{aligned}
    F_0 &= 0 \\
    F_1 &= 1 \\
    F_2 &= 1 \\
    F_3 &= 2 \\
    \ldots\\
    F_n &= F_{n-1} + F_{n-2}
    \end{aligned}
    ```
"""

# ╔═╡ 0cfc3b0d-c662-4399-9c60-3bc851890338
#=╠═╡
function fibonacci(n)
	nothing
end
  ╠═╡ =#

# ╔═╡ e3c394dc-b8f7-408a-a129-df42b5e0914d
#=╠═╡
if isnothing(fibonacci(0))
	still_nothing()
elseif (fibonacci(0) == 0 && fibonacci(1) == 1 && fibonacci(2) == 1 && fibonacci(7) == 13)
	correct()
elseif (fibonacci(0) == 0 && fibonacci(1) == 1)
	almost(md"Good start `fibonacci(0)` and `fibonacci(1)` are correct. Keep going!")
elseif !type_eq(fibonacci(1), Int)
	almost(md"Return type is not integer. Make sure all numbers in your implementation are integers. E.g. use integers like `1` and not floats like `1.0`.")
else
	keep_working()
end
  ╠═╡ =#

# ╔═╡ ffee7d80-eded-11ea-26b1-1331df204c67
md"## Arrays"

# ╔═╡ 903b52a1-47dd-4924-bd10-4d4573f7e7a8
md"""
An array is a collection of elements. They can be one-dimensional, corresponding to a list or vector. They can be two-dimensional, corresponding to a grid of numbers or matrix. They can have more dimensions.

Arrays have the type `Array{T, N}` where:
- `T` is the type of element inside the array.
- `N` is the number of dimensions in the array: `1` for a vector, `2` for a matrix, etc...

`Vector{T}` is an alias for `Array{T, 1}` and `Matrix{T}` is an alias for `Array{T, 1}`.

!!! tip
    If you are already familiar with Python or Matlab, you should check out this comparative cheatsheet already linked in the introduction, which covers a lot of Julia's array syntax: [Comparative cheatsheet Python <-> Julia <-> Matlab](https://cheatsheets.quantecon.org/).
"""

# ╔═╡ cae4137e-edee-11ea-14af-59a32227de1b
md"### 1D arrays (`Vector`s)"

# ╔═╡ 714f4fca-edee-11ea-3410-c9ab8825d836
md"We can make a `Vector` (1-dimensional, or 1D array) using square brackets:"

# ╔═╡ 82cc2a0e-edee-11ea-11b7-fbaa5ad7b556
v = [1, 2, 3]

# ╔═╡ 85916c18-edee-11ea-0738-5f5d78875b86
typeof(v)

# ╔═╡ 881b7d0c-edee-11ea-0b4a-4bd7d5be2c77
md"""
The type tells us that this is a 1D array of integers.
**Don't forget the commas (`,`)** between the elements, or you will create a matrix instead:
"""

# ╔═╡ 3384d534-9871-484d-aa7f-c07ecb8b0c49
[1 2 3]

# ╔═╡ 839463f9-0fbd-456e-bca5-7d844048fd78
md"""
We access elements using square brackets:
"""

# ╔═╡ 547120f6-0b1d-4a17-8eea-ddbb5e2e1999
v[1]

# ╔═╡ a298e8ae-edee-11ea-3613-0dd4bae70c26
v[2]

# ╔═╡ 8e28a294-df42-4f85-a262-8333272a9510
md"""
The special syntax **`end`** can be used to refer to the end of the array. For example, `v[end]` will access the last element, and `v[end-1]` the one before.
"""

# ╔═╡ ae0c73b7-3148-4137-ad4f-163ad66d36a6
v[end-1]

# ╔═╡ 50d57c04-57a3-40bd-bd54-4bfefdf26e8b
md"""
**In Julia, arrays start at `1`. Accessing an array at index `0` will cause an error:**
"""

# ╔═╡ ffa187da-5e70-4317-b88f-922ea31cc8c7
v[0]

# ╔═╡ 17998b62-91df-4122-b4fa-a50e4bf6aaef
md"""
Arrays can be modified:
"""

# ╔═╡ a5ebddd6-edee-11ea-2234-55453ea59c5a
v[2] = 10

# ╔═╡ 1cec47b2-d37b-4f04-bc80-c5a06a925083
md"""
However types matter! We cannot store a decimal number in an array of `Int64`:
"""

# ╔═╡ 9ef18d2c-8f3d-4e47-97bf-50ff096d7f74
v[2] = 2.5

# ╔═╡ a9b48e54-edee-11ea-1333-a96181de0185
md"Note that Pluto does not automatically update cells when you modify elements of an array, but the value does change."

# ╔═╡ 68c4ead2-edef-11ea-124a-03c2d7dd6a1b
md"A nice way to create `Vector`s following a certain pattern is to use an **array comprehension**:"

# ╔═╡ 84129294-edef-11ea-0c77-ffa2b9592a26
v2 = [i^2 for i in 1:10]

# ╔═╡ 872cbae6-7dd6-4567-b582-2ce55109f46d
md"""
## Element-wise array operations
Arrays can be added together, or subtracted, using the regular `+` and `-` operators.
"""

# ╔═╡ 15e4b26d-6d7e-45bf-939a-a498c11bf741
vec1 = [1, 2, 3]

# ╔═╡ cfa2a8d5-467e-4084-aea3-40275a2e9470
vec2 = [1, 4, 9]

# ╔═╡ a8401aa0-34d7-449c-8877-a96cc7875393
vec1 + vec2

# ╔═╡ aeb60671-3d68-47fb-b6d6-4aa846024558
vec1 - vec2

# ╔═╡ 114e304e-8642-4b00-a566-7c014b7cf326
md"""
This adds and substracts corresponding elements of each array. We call this an **element-wise** operation.

In general, other mathematical operator or functions will not work on arrays. To apply operators or functions to an array **element-wise**, we use special syntax:
- Adding `.` before an operator will apply the operator element-wise. For example, `.*` and `./` will perform element-wise multiplication.
- Adding `.` after a function will apply the function element-wise. For example, `cos.` will compute the cosine element-wise.
"""

# ╔═╡ ec85d0a2-7f82-4a36-b2a8-429ec7015140
vec1 .* vec2

# ╔═╡ e84ad8a0-853c-4886-8133-87c9c1fe6fa3
vec1 ./ vec2

# ╔═╡ 5ccc054f-7f31-4ca3-8769-6d4aa28ce8d0
vec1 .^ 2

# ╔═╡ 612aee5c-648d-469f-9ec9-c37300005b11
vec1.^2 == vec2

# ╔═╡ 0cbf039b-bb06-44af-a349-55b614c94e5b
exp.(vec1)

# ╔═╡ 9632c36a-8fe2-4e11-b949-8213a6ac4bdb
cos.([0, π, 2π])

# ╔═╡ 39bc95f9-74e8-468d-b7e4-700438675bf0
md"""
Let's put this into practice with the Fibonacci function you wrote above:

!!! exercise
    Complete the following function that will compute multiple Fibonacci numbers at once. It will receive a vector of integers, and should return a vector with the corresponding Fibonacci numbers.
    For example, `many_fibonacci([1, 2, 4])` should return `[1, 1, 3]`.

    Your answer should be very short and use the `fibonacci` function that you wrote above.
"""

# ╔═╡ f538c5f7-f5c3-4c47-9db3-b785d8b7b9ca
#=╠═╡
many_fibonacci(ns) = fibonacci(0)
  ╠═╡ =#

# ╔═╡ 8f51ebf4-7dc7-4441-a0b6-9225acbbbfb0
#=╠═╡
if isnothing(many_fibonacci([0]))
	still_nothing(text=md"First, implement the `fibonacci` function in the exercise a few sections above.")
elseif many_fibonacci([10]) == 0
	still_nothing(md"Replace `fibonacci(0)` in the cell above.")
elseif many_fibonacci([1, 2, 4]) != [1, 1, 3]
	keep_working(md"The example `many_fibonacci([1, 2, 4])` is not correct yet.")
elseif many_fibonacci([1, 2, 4, 7, 10]) != [1, 1, 3, 13, 55]
	keep_working()
else
	correct()
end
  ╠═╡ =#

# ╔═╡ d364fa16-edee-11ea-2050-0f6cb70e1bcf
md"## 2D arrays (matrices)"

# ╔═╡ db99ae9a-edee-11ea-393e-9de420a545a1
md"We can make small matrices (2D arrays) with square brackets too:"

# ╔═╡ 04f175f2-edef-11ea-0882-712548ebb7a3
M = [1 2;
	 3 4]

# ╔═╡ 0a8ac112-edef-11ea-1e99-cf7c7808c4f5
typeof(M)

# ╔═╡ 1295f48a-edef-11ea-22a5-61e8a2e1d005
md"The `2` in the type confirms that this is a 2D array."

# ╔═╡ 3e1fdaa8-edef-11ea-2f03-eb41b2b9ea0f
md"This won't work so easily for larger matrices, though. For that we can use e.g."

# ╔═╡ 48f3deca-edef-11ea-2c18-e7419c9030a0
zeros(5, 5)

# ╔═╡ a8f26af8-edef-11ea-2fc7-2b776f515aea
md"Note that `zeros` gives `Float64`s by default. We can also specify a type for the elements:"

# ╔═╡ b595373e-edef-11ea-03e2-6599ef14af20
zeros(Int, 4, 5)

# ╔═╡ c6676b63-aa1d-4698-ad33-e8a354534164
md"Same as the arrays themselves. E.g. contrast"

# ╔═╡ b1d27139-6a41-414a-8bc2-6b54373f5404
Float64[1, 2, 3]

# ╔═╡ f85a45d7-2467-470b-b688-5d9f0f69b133
md"which creates a `Float64` array versus"

# ╔═╡ 9449e4ad-bef7-4342-9623-5583308f4d28
Int[1, 2, 3]  # or just [1, 2, 3]

# ╔═╡ 4cb33c04-edef-11ea-2b35-1139c246c331
md"which creates an integer array. We can then fill in the values we want by manipulating the elements, e.g. with a `for` loop."

# ╔═╡ 54e47e9e-edef-11ea-2d75-b5f550902528
md"A nice alternative syntax to create matrices following a certain pattern is an array comprehension with a *double* `for` loop:"

# ╔═╡ 6348edce-edef-11ea-1ab4-019514eb414f
[i + j for i in 1:5, j in 1:6]

# ╔═╡ 2f4d9a36-d93d-4cf3-aa73-23442fad0bac
md"""
To access matrix elements directly, we need to use two indices:
"""

# ╔═╡ 237bc17e-ae95-4d7f-b6e3-fa54ad49234f
M[1, 2]

# ╔═╡ 75341b42-172a-4a18-8bd6-1d68859a0380
M[2, 2]

# ╔═╡ 953a86dc-0013-4e80-9810-ffcae96444e9
md"""
## Element-wise vs matrix operations
Element-wise operations also work on matrices, using the `.` syntax explained above:
"""

# ╔═╡ 75825042-cc74-40f9-8b36-9493b6dbeb45
exp.(M)

# ╔═╡ 3e9b1407-f4ca-47b8-be67-0b1b402cf6a3
md"""
**Many functions can be applied to the whole square matrix, giving a different result than applying the function element-wise.**
"""

# ╔═╡ 3acae1a9-562f-4d09-b4e0-dd973478d095
exp(M)

# ╔═╡ 09855245-b441-4af9-a1fe-f399d9e24b1c
md"""
Oops, we computed the [matrix exponential](https://en.wikipedia.org/wiki/Matrix_exponential) instead of the element-wise exponential.
"""

# ╔═╡ 768b7bdd-0348-4cdf-bd0b-159638449f7f
md"""
As another example, **matrix multiplication** is performed with `*` whereas element-wise multiplication is performed with `.*`:
"""

# ╔═╡ 91d34d35-0b1d-4331-9622-8b3619b0b20c
M * M

# ╔═╡ 2e0aecbc-8368-42ff-abc7-62e96bbb44a2
M .* M

# ╔═╡ f3cd579b-dd2a-4677-b410-bd9824f815f6
md"""
As a final example, let's revisit two variations of the famous $\cos(x)^2 + \sin(x)^2 = 1$ equation:
"""

# ╔═╡ 6f952243-5b13-45d5-9d78-b0e9fca84c25
cos(M)^2 + sin(M)^2

# ╔═╡ 2b9177e4-49cf-4687-9677-fb3939d39b83
cos.(M).^2 + sin.(M).^2

# ╔═╡ 8b14d642-c0a6-430e-9050-3d31544e8703
md"""
Element-wise, we indeed get `1` for every entry. For the whole matrix, we get an identity matrix which is indeed the `1` of 2x2 matrices.

Note: `-1.66533e-16` means $-1.66533 \times 10^{-16}$ which is a very small number. It happens due to the imprecision of floating-point arithmetic. It would be `0` if computer math were exact.
"""

# ╔═╡ e7958456-9f46-4273-96ef-014ee631efe3
md"""
# The End
This concludes the tutorial section of this first notebook.

Feel free to keep reading for a look at a more advanced example, or skip the following section and come back to it later.
In any case, a second notebook awaits you on Moodle with an introduction to plotting in Julia.
"""

# ╔═╡ ff2060d7-017f-43f4-af6f-ca3c921fe975
md"""
# Optional: Step-by-step replication of the Arrhenius fit from the lecture
_(This part is intented as a deeper dive into Julia's capabilities.)_


We here reproduce in full detail the Arrhenius fit example from the lecture.

- We were given the data:
"""

# ╔═╡ d0501dcf-8cc0-4e80-bab9-7a8421e082cb
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


# ╔═╡ ff49ef7e-ec35-42b9-b267-a4dfc4380a2e
md"""
This data is in plain text form, so we need to preprocess it a little in order to be able to plot it graphically.

First we split the data into lines:
"""

# ╔═╡ 4543612a-7a80-456c-810e-0a1b344941c9
lines = split(data, "\n")

# ╔═╡ 821a8995-2a94-41fb-89e1-f305b37074f4
md"""
The head line
"""

# ╔═╡ 791e8815-45ba-469c-8849-66a63f455687
lines[1]

# ╔═╡ c669ca4a-06c6-4f0f-8324-929ebe95c228
md"""
is not needed, and similarly the last line
"""

# ╔═╡ f0afcf68-9796-47d5-8616-d6e4b90310b9
lines[end]

# ╔═╡ 1c18d499-5cb2-4afa-9fb8-056cfe538c40
md"""
is empty, so we strip them using Julia's array masks:
"""

# ╔═╡ 09b9eb1c-eb18-4c05-8dd6-ae4686da32f0
lines_relevant = lines[2:end-1]

# ╔═╡ cb18cacb-3d77-4822-b211-cdcba564ede7
md"""Repeating the spliting on the lines, we can obtain the temperature and rate data in separate arrays:"""

# ╔═╡ 021021b6-87db-4ea6-9d0c-dcff04d8a9f8
temperature_string = [split(line)[1] for line in lines_relevant]

# ╔═╡ ba4ed27b-d910-4994-b885-6e1d1aadcbb2
rate_string = [split(line)[2] for line in lines_relevant]

# ╔═╡ 42dc5569-00db-490b-91de-e6ad2cf892ec
md"""To get floating-point numbers out of these strings, we `parse` them:"""

# ╔═╡ ba1b06f9-2bac-4270-885a-3ab7b617f6fa
[parse(Float64, string) for string in temperature_string]

# ╔═╡ f952e9e3-3bd6-4ea3-a717-e330ccbc9234
md"""More compactly we could have written this as"""

# ╔═╡ d6296330-0f76-4a74-9fc9-a5ffa0b370f1
begin
	temperature = [parse(Float64, split(line)[1]) for line in lines[2:end-1]]
	rate        = [parse(Float64, split(line)[2]) for line in lines[2:end-1]]
end;

# ╔═╡ bf1b8bc6-651e-4aaa-9449-8ecbe9dcac86
md"""where `begin ... end` allows to combine multiple statements in one cell.

Finally we plot:"""

# ╔═╡ 2eece0f2-ffec-4f55-be39-86178519f8cf
scatter(temperature, rate)

# ╔═╡ 6b82c743-7be1-4980-9c1b-ef2f86eb474f
md"""
The best fit Arrhenius equation is given by the function
"""

# ╔═╡ 67c63a8e-88c2-4d17-9225-30819f0886a3
function arrhenius(T)
	5exp(-300 / T)
end

# ╔═╡ cfb8897a-4312-450b-920a-e5ec90e88565
md"""
Let's first investigate which values this function would take at these points
and plot it in the same graph:
"""

# ╔═╡ 69432575-0503-4331-8ac1-2b1f90f1f202
begin
	predicted_rates = [arrhenius(t) for t in temperature]

	plt1 = scatter(temperature, rate, label="Measured data")
	scatter!(plt1, temperature, predicted_rates, label="Predicted rate (Arrhenius)")
	plt1
end

# ╔═╡ d009e504-dca9-43e0-8668-273a303edd79
md"Note that here the first `scatter` call returns an object, which represents the plotting canvas. Using a second `scatter!` call we can add other entities for plotting to it.

If we want to plot a continuous graph instead of data points, we can use `plot` or `plot!`, e.g.:"

# ╔═╡ 71d64b5a-4e0c-4dea-bc87-479e8e6aa7e0
begin
	# Note: predicted_rates already defined above
	
	plt2 = scatter(temperature, rate, label="Measured data")
	plot!(plt2, temperature, predicted_rates, label="Predicted rate (Arrhenius)",
	      linewidth=2)
	plt2
end

# ╔═╡ ca2ea37d-21cd-404e-bcb0-11b12c6af9b7
md"For convenience, functions like `arrhenius` can also be plotted directly, without evaluating them explicitly:"

# ╔═╡ f23a8a39-0cd1-4c18-b598-41fc3c0581ad
begin
	plt3 = scatter(temperature, rate, label="Measured data")
	plot!(plt3, arrhenius, label="Predicted rate (Arrhenius)",
	      linewidth=2)
	plt3
end

# ╔═╡ 8ff0d83c-9a5a-419d-9a0d-efe9f0034fc2
md"""
Note, that also a similar `let ... end` block exists to combine multiple statements.
The point of `let` is to *hide* variable names and content from the outside the cell.
This is needed because in Pluto notebooks each variable name may only be used a single time in the public context, e.g. defining `a` twice leads to a warning
and the disabling of one cell.
"""

# ╔═╡ ab70be10-5821-44ac-8a5b-aeb32f24d1cf
md"""
This is a safety feature to avoid overwriting computational results.

Sometimes (especially for setting up plots), we have no interest in using the data outside of the cell anyway. In this case using `let` ... `end` is usually better:
"""

# ╔═╡ 562622cc-d910-4cc5-a7f5-f5b919923a47
let
	p = scatter(temperature, rate, label="Measured data")
	plot!(p, arrhenius, label="Predicted rate (Arrhenius)",
	      linewidth=2)
	p
end

# ╔═╡ 34c12be4-e036-4a1b-a3f8-72530073c6d4
# ╠═╡ disabled = true
#=╠═╡
variable = 4
  ╠═╡ =#

# ╔═╡ 7c3a5972-a0b8-44ba-8525-ade8c0bc4135
#=╠═╡
variable = 5
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.41.4"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "6479b20828288c84692da3cd72d5bbf8e4c6ebab"

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
git-tree-sha1 = "6949039c4db076a5967228bf9d7b6f32f2324328"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.20"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "04a6eb1e57dc3a6992b5af5aa753a2f8d3a1028e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.20+0"

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
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

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
# ╠═33490d1f-ef11-4159-823f-0040c096539c
# ╟─0d3aec92-edeb-11ea-3adb-cd0dc17cbdab
# ╟─4dbd5dd5-02d9-4243-a0fc-ac7089fc588a
# ╠═5e963571-a2bf-4554-9082-1ab5ec48d089
# ╟─462e68ff-5c86-4403-bcc5-907dbca15f2e
# ╠═6424fb26-3b04-48d9-9c6f-0a1db75ef35f
# ╟─a249a972-8911-4c74-884d-5e30fd43dbec
# ╟─3b038ee0-edeb-11ea-0977-97cc30d1c6ff
# ╟─6af9c3a0-2534-4bf1-9eb8-f5e34a091c64
# ╠═28c2ac78-2462-452e-bdc5-8319576fb33a
# ╠═0c8d79af-d762-4165-9992-688532185d48
# ╠═48ebe8f9-7ee1-478a-8ea3-f12f8b8b25e6
# ╠═9195da81-4465-4fc1-9831-b9aea1697bbe
# ╠═76f3caa3-4541-4268-aac6-6fb4d7aa6436
# ╟─5e062a24-edeb-11ea-256a-d938f77d7815
# ╟─166e8c0c-73d7-4a4d-b27f-13f9b4005dd2
# ╠═52f0d7e9-5eb2-42e6-b5ce-be8c6dd0ea4f
# ╠═025abeca-9f34-4921-9578-c571234c0616
# ╠═4f9f0cc7-f6cf-443b-94af-bb31f213503e
# ╠═f608bad0-43a6-4796-bfc8-7f30499797ad
# ╟─b5143875-138b-413b-8073-1dc613555803
# ╠═3e8e0ea0-edeb-11ea-22e0-c58f7c2168ce
# ╠═59b66862-edeb-11ea-2d62-71dcc79dbfab
# ╟─e938820b-daad-4cc9-9fd0-c8adc84adf46
# ╠═b6a3144d-eaba-4dff-802f-d20becff06e6
# ╠═8ad688a5-2831-4017-ad39-2c98e5a0efd9
# ╟─0d0f4edf-a411-4a51-9c51-24256050b397
# ╠═2a4ec691-0944-4b3a-8def-a83eb8bf2400
# ╟─72e5b9d3-df98-4f6e-b283-52f928e7a80c
# ╟─44986e44-3371-477f-967e-9af708d103b3
# ╠═35887694-5ad2-463a-b4a9-9f115542c9e2
# ╟─5e381128-b6f6-4594-8baa-f585d1e44061
# ╟─1b2086c1-9b29-4b7d-8e88-7e8fbd86bde3
# ╟─7e46f0e8-edeb-11ea-1092-4b5e8acd9ee0
# ╠═5111713e-8c05-4161-b55a-814d72d32bc9
# ╠═799938f0-8a4d-4505-bdcf-d5de7f70dcbb
# ╠═57165a23-204d-4e68-8943-10bee7d68fa5
# ╟─89336d7d-ebb8-46cc-92b2-dea09e3282c0
# ╠═8a695b86-edeb-11ea-08cc-17263bec09df
# ╟─769ed521-7faa-459a-b77e-15d5fbeef046
# ╟─685174e3-21ca-46f7-bd6d-467718a6850e
# ╠═71bcb101-8715-4754-82ce-1d24fea8a512
# ╟─d3287858-b154-42f6-9046-7030e6b8d722
# ╟─9072a34d-862f-4d6a-bbc4-7fce7a55c575
# ╠═19c956cb-776b-4628-911d-a119f51c62d8
# ╟─cefeb394-9322-432d-8f25-6f70f1d15432
# ╠═b9e1c787-3594-4107-b1f0-0293a666d565
# ╟─d3b94264-ff93-455e-a44c-3798cdacc660
# ╠═ae409c89-ec72-4d7f-8d39-79aa8946e533
# ╟─4a4dcb69-89b9-4622-b43f-d806bb27e861
# ╠═0af19cc8-d47c-4e01-afe9-d4f5d9225396
# ╟─b7684e96-5a32-4e12-bbca-92ce7f57ce70
# ╠═f7b4f385-78ed-4191-b2a9-fa4423b58f72
# ╟─273e451c-17fc-4816-831c-82f0c1f2f2ad
# ╠═8ca29458-10dc-401a-b681-e63564552cf4
# ╟─8e2dd3be-edeb-11ea-0703-354fb31c12f5
# ╟─82a0aae9-3937-4f60-9f75-532feb59e3f7
# ╠═80074898-5414-445b-8831-4c94fd8c9d52
# ╠═1d945be1-3664-462c-ba4a-29e5e72724ea
# ╟─9f28c628-dc42-4d4d-bf5d-d85a39cc1271
# ╠═b8acc676-3e5b-4048-8120-ad6a182d2315
# ╠═2a3f7d44-0fb6-4929-9248-5265d0680619
# ╠═dc29fcf1-ed6c-44a4-9e89-57d3f635a542
# ╠═7db7d2fd-4004-41bd-bfe9-b1ead83bdfea
# ╟─13a0268c-d815-4844-8b90-b2a475e9a4fa
# ╟─96b5a28c-edeb-11ea-11c0-597615962f54
# ╠═a7453572-edeb-11ea-1e27-9f710fd856a6
# ╟─b341db4e-edeb-11ea-078b-b71ac00089d7
# ╠═23f9afd4-eded-11ea-202a-9f0f1f91e5ad
# ╟─d50fdb40-62eb-454d-baf3-42d1bb34f082
# ╠═cc1f6872-edeb-11ea-33e9-6976fd9b107a
# ╟─ce9667c2-edeb-11ea-2665-d789032abd11
# ╠═d73d3400-edeb-11ea-2dea-95e8c4a6563b
# ╠═e04ccf10-edeb-11ea-36d1-d11969e4b2f2
# ╟─2c446801-df67-4912-942c-aaa664302b28
# ╠═fc0bc6bd-b01c-461d-8cf8-ad3912ab2cbd
# ╠═b5e5cee9-dcf2-4e23-afce-72e136cbd988
# ╟─d15f35df-23af-4793-bb4a-5878a4a52b23
# ╠═971be2b7-d1e6-4e16-be24-a6a8d3092772
# ╟─7cfb56d4-b978-4e02-a64e-96431db95b61
# ╟─93a231f4-edec-11ea-3b39-299b3be2da78
# ╟─82e63a24-eded-11ea-3887-15d6bfabea4b
# ╠═9b339b2a-eded-11ea-10d7-8fc9a907c892
# ╠═9535eb40-eded-11ea-1651-e33c9c23dbfb
# ╠═73f81dc2-1138-4f44-be8c-fc0ee55aa2e2
# ╠═58053efc-8656-4b9c-abfc-a15c29cb10b3
# ╠═530f923d-c820-4a45-8215-a9a57e375a18
# ╟─a16299a2-eded-11ea-2b56-93eb7a1010a7
# ╠═bc6b124e-eded-11ea-0290-b3760cb81024
# ╟─cfb21014-eded-11ea-1261-3bc30952a88e
# ╟─ca7089ae-0394-4d52-82cc-b893e2300894
# ╠═ee0510a9-95b3-4e43-addb-8ffdb724d3ad
# ╟─b96761f5-b235-47b0-b3b4-f350feea46ac
# ╠═22bd1bc0-c2fa-489a-a2c9-b1c75b36f7ba
# ╠═f2426c25-32c6-4408-80c3-d234a0fd6a65
# ╠═b8201bdd-585a-4ac2-ae44-bb2c1b378dca
# ╠═18adab01-e74e-4537-8549-fadad03f30b7
# ╟─77d9277c-54bc-4d80-acc7-d8d866753e27
# ╠═54a7f75f-9835-497e-a1da-ecc36ac61641
# ╟─4053a6a7-54a0-4b05-aa44-f6ffaf38cfe7
# ╟─e297c5cc-edeb-11ea-3bdd-090f415685ab
# ╟─ec751446-edeb-11ea-31ba-2372e7c71b42
# ╠═fe3fa290-edeb-11ea-121e-7114e5c573c1
# ╟─394b0ec8-eded-11ea-31fb-27392068ef8f
# ╠═4dc00908-eded-11ea-25c5-0f7b2b7e18f9
# ╟─6c44abb4-edec-11ea-16bd-557800b5f9d2
# ╠═683af3e2-eded-11ea-25a5-0d90bf099d98
# ╠═76764ea2-eded-11ea-1aa6-296f3421de1c
# ╟─bd39bd6d-152a-48d2-bc7b-43f426eb1165
# ╠═ddcac2dc-f5b4-47fe-81f2-262b1b4718b2
# ╟─1b49f71c-a8b4-440b-a46a-17dfb4220fa5
# ╠═0cfc3b0d-c662-4399-9c60-3bc851890338
# ╟─e3c394dc-b8f7-408a-a129-df42b5e0914d
# ╟─ffee7d80-eded-11ea-26b1-1331df204c67
# ╟─903b52a1-47dd-4924-bd10-4d4573f7e7a8
# ╟─cae4137e-edee-11ea-14af-59a32227de1b
# ╟─714f4fca-edee-11ea-3410-c9ab8825d836
# ╠═82cc2a0e-edee-11ea-11b7-fbaa5ad7b556
# ╠═85916c18-edee-11ea-0738-5f5d78875b86
# ╟─881b7d0c-edee-11ea-0b4a-4bd7d5be2c77
# ╠═3384d534-9871-484d-aa7f-c07ecb8b0c49
# ╟─839463f9-0fbd-456e-bca5-7d844048fd78
# ╠═547120f6-0b1d-4a17-8eea-ddbb5e2e1999
# ╠═a298e8ae-edee-11ea-3613-0dd4bae70c26
# ╟─8e28a294-df42-4f85-a262-8333272a9510
# ╠═ae0c73b7-3148-4137-ad4f-163ad66d36a6
# ╟─50d57c04-57a3-40bd-bd54-4bfefdf26e8b
# ╠═ffa187da-5e70-4317-b88f-922ea31cc8c7
# ╟─17998b62-91df-4122-b4fa-a50e4bf6aaef
# ╠═a5ebddd6-edee-11ea-2234-55453ea59c5a
# ╟─1cec47b2-d37b-4f04-bc80-c5a06a925083
# ╠═9ef18d2c-8f3d-4e47-97bf-50ff096d7f74
# ╟─a9b48e54-edee-11ea-1333-a96181de0185
# ╟─68c4ead2-edef-11ea-124a-03c2d7dd6a1b
# ╠═84129294-edef-11ea-0c77-ffa2b9592a26
# ╟─872cbae6-7dd6-4567-b582-2ce55109f46d
# ╠═15e4b26d-6d7e-45bf-939a-a498c11bf741
# ╠═cfa2a8d5-467e-4084-aea3-40275a2e9470
# ╠═a8401aa0-34d7-449c-8877-a96cc7875393
# ╠═aeb60671-3d68-47fb-b6d6-4aa846024558
# ╟─114e304e-8642-4b00-a566-7c014b7cf326
# ╠═ec85d0a2-7f82-4a36-b2a8-429ec7015140
# ╠═e84ad8a0-853c-4886-8133-87c9c1fe6fa3
# ╠═5ccc054f-7f31-4ca3-8769-6d4aa28ce8d0
# ╠═612aee5c-648d-469f-9ec9-c37300005b11
# ╠═0cbf039b-bb06-44af-a349-55b614c94e5b
# ╠═9632c36a-8fe2-4e11-b949-8213a6ac4bdb
# ╟─39bc95f9-74e8-468d-b7e4-700438675bf0
# ╠═f538c5f7-f5c3-4c47-9db3-b785d8b7b9ca
# ╟─8f51ebf4-7dc7-4441-a0b6-9225acbbbfb0
# ╟─d364fa16-edee-11ea-2050-0f6cb70e1bcf
# ╟─db99ae9a-edee-11ea-393e-9de420a545a1
# ╠═04f175f2-edef-11ea-0882-712548ebb7a3
# ╠═0a8ac112-edef-11ea-1e99-cf7c7808c4f5
# ╟─1295f48a-edef-11ea-22a5-61e8a2e1d005
# ╟─3e1fdaa8-edef-11ea-2f03-eb41b2b9ea0f
# ╠═48f3deca-edef-11ea-2c18-e7419c9030a0
# ╟─a8f26af8-edef-11ea-2fc7-2b776f515aea
# ╠═b595373e-edef-11ea-03e2-6599ef14af20
# ╟─c6676b63-aa1d-4698-ad33-e8a354534164
# ╠═b1d27139-6a41-414a-8bc2-6b54373f5404
# ╟─f85a45d7-2467-470b-b688-5d9f0f69b133
# ╠═9449e4ad-bef7-4342-9623-5583308f4d28
# ╟─4cb33c04-edef-11ea-2b35-1139c246c331
# ╟─54e47e9e-edef-11ea-2d75-b5f550902528
# ╠═6348edce-edef-11ea-1ab4-019514eb414f
# ╟─2f4d9a36-d93d-4cf3-aa73-23442fad0bac
# ╠═237bc17e-ae95-4d7f-b6e3-fa54ad49234f
# ╠═75341b42-172a-4a18-8bd6-1d68859a0380
# ╟─953a86dc-0013-4e80-9810-ffcae96444e9
# ╠═75825042-cc74-40f9-8b36-9493b6dbeb45
# ╟─3e9b1407-f4ca-47b8-be67-0b1b402cf6a3
# ╠═3acae1a9-562f-4d09-b4e0-dd973478d095
# ╟─09855245-b441-4af9-a1fe-f399d9e24b1c
# ╟─768b7bdd-0348-4cdf-bd0b-159638449f7f
# ╠═91d34d35-0b1d-4331-9622-8b3619b0b20c
# ╠═2e0aecbc-8368-42ff-abc7-62e96bbb44a2
# ╟─f3cd579b-dd2a-4677-b410-bd9824f815f6
# ╠═6f952243-5b13-45d5-9d78-b0e9fca84c25
# ╠═2b9177e4-49cf-4687-9677-fb3939d39b83
# ╟─8b14d642-c0a6-430e-9050-3d31544e8703
# ╟─e7958456-9f46-4273-96ef-014ee631efe3
# ╟─ff2060d7-017f-43f4-af6f-ca3c921fe975
# ╠═d0501dcf-8cc0-4e80-bab9-7a8421e082cb
# ╟─ff49ef7e-ec35-42b9-b267-a4dfc4380a2e
# ╠═4543612a-7a80-456c-810e-0a1b344941c9
# ╟─821a8995-2a94-41fb-89e1-f305b37074f4
# ╠═791e8815-45ba-469c-8849-66a63f455687
# ╟─c669ca4a-06c6-4f0f-8324-929ebe95c228
# ╠═f0afcf68-9796-47d5-8616-d6e4b90310b9
# ╟─1c18d499-5cb2-4afa-9fb8-056cfe538c40
# ╠═09b9eb1c-eb18-4c05-8dd6-ae4686da32f0
# ╟─cb18cacb-3d77-4822-b211-cdcba564ede7
# ╠═021021b6-87db-4ea6-9d0c-dcff04d8a9f8
# ╠═ba4ed27b-d910-4994-b885-6e1d1aadcbb2
# ╟─42dc5569-00db-490b-91de-e6ad2cf892ec
# ╠═ba1b06f9-2bac-4270-885a-3ab7b617f6fa
# ╟─f952e9e3-3bd6-4ea3-a717-e330ccbc9234
# ╠═d6296330-0f76-4a74-9fc9-a5ffa0b370f1
# ╟─bf1b8bc6-651e-4aaa-9449-8ecbe9dcac86
# ╠═2eece0f2-ffec-4f55-be39-86178519f8cf
# ╟─6b82c743-7be1-4980-9c1b-ef2f86eb474f
# ╠═67c63a8e-88c2-4d17-9225-30819f0886a3
# ╟─cfb8897a-4312-450b-920a-e5ec90e88565
# ╠═69432575-0503-4331-8ac1-2b1f90f1f202
# ╟─d009e504-dca9-43e0-8668-273a303edd79
# ╠═71d64b5a-4e0c-4dea-bc87-479e8e6aa7e0
# ╟─ca2ea37d-21cd-404e-bcb0-11b12c6af9b7
# ╠═f23a8a39-0cd1-4c18-b598-41fc3c0581ad
# ╟─8ff0d83c-9a5a-419d-9a0d-efe9f0034fc2
# ╠═34c12be4-e036-4a1b-a3f8-72530073c6d4
# ╠═7c3a5972-a0b8-44ba-8525-ade8c0bc4135
# ╟─ab70be10-5821-44ac-8a5b-aeb32f24d1cf
# ╠═562622cc-d910-4cc5-a7f5-f5b919923a47
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
