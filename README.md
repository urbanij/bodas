## bodas
[![Downloads](https://pepy.tech/badge/bodas)](https://pepy.tech/project/bodas)

Asymptotic Bode plots in Python.

![](https://github.com/urbanij/bodas/blob/main/docs/example1.png?raw=true)

### Installation
`pip install bodas`


### Simple usage example
```python
In [1]: import bodas 

In [2]: import sympy                # import [SymPy](https://www.sympy.org) the Python 
                                    # library for symbolic mathematics

In [3]: s = sympy.Symbol('s')       # define `s` as symbol

In [4]: H = (1+s/23)/(1+s/123)**2   # assign to `H` the function you want to plot

In [5]: sympy.pretty_print(H)
  s
  ── + 1
  23
──────────
         2
⎛ s     ⎞
⎜─── + 1⎟
⎝123    ⎠


In [6]: bodas.plot(H)               # call the `plot` function defined in the bodas library
```

---

### Todo
See [TODO.md](https://github.com/urbanij/bodas/blob/main/TODO.md)
