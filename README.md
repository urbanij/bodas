## bodas

Asymptotic Bode plots in Python.

![](https://github.com/urbanij/bodas/blob/main/docs/example1.png?raw=true)

### Installation
`pip install bodas`


### Basic usage
```python
import bodas 
import sympy            # import [SymPy](https://www.sympy.org) the Python 
                        # library for symbolic mathematics

s = sympy.Symbol('s')   # define `s` as symbol

H = 1/(1 + s/120)       # assign to `H` the function you want to plot

bodas.plot( str(H) )    # call the `plot` function defined in the bodas library
```

---

### Todo
See [TODO.md](https://github.com/urbanij/bodas/blob/main/TODO.md)
