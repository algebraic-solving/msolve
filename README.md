# Multivariate polynomial system solver 

[https://msolve.lip6.fr](https://msolve.lip6.fr)

# Install Instructions

See INSTALL.

# Input File Format

More informations are given in the documentation (see [https://msolve.lip6.fr](https://msolve.lip6.fr))

`msolve` input files need to be of the following format:

**1st line**: variables as commata separated list, e.g. `x1,x2,x3,x4,y1,y2`.<br/>
**2nd line**: field characteristic, e.g. `0`.<br/>
**following lines**: generating polynomials, all but the last one need
to terminate by a `,`, e.g.
```
x1,x2,x3,x4,y1,y2
101
x1+x2+x3+x4,
2*y1-145*y2
```
Polynomials may be multiline, thus `,` as a separator.

Coefficients can be rational, using `/`, e.g. `-2/3*x2*y1^2+...`.
