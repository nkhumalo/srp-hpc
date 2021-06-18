# Optimizing the parameters

The code we are using currently uses a Monte Carlo method to optimize the
wave function. While this is slow it does converge and produce reliable
results. The natural orbital functional parameters do not converge
very well with Monte Carlo methods. However, maybe a quasi-Newton-Raphson
method works for these parameters.

First of all note that we are trying to minimize the error in the potential
energy surface produced by the functional. This error is defined as the 
Infinity norm of the difference between the Full-CI potential energy
surfaces and the WFN1 potential energy surfaces. 

Assume that the function we want to minimize can be aproximated as a
quadratic function:
```
        f(x) = a*x^2 + b*x + c
```
If we choose three points and compute the corresponding `f(x)`
```
        x1 = -d   f(x1) = f(-d)
        x2 =  0   f(x2) = f( 0)
        x3 = +d   f(x3) = f(+d)
```
then the parameters `a`, `b`, and `c` can be obtained as
```
        c = f(x2)
        b = (f(x3)-f(x1))/(2*d)
        a = (f(x3)+f(x1)-2*f(x2))/(2*d^2)
```
and the optimal step in `x` (let us call that `dx`) is
```
        dx = -b/a
```
