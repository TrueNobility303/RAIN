# RAIN (Recursive Anchoing IterRation)

Offical codes for "Near-optimal algorithms for making the gradient small in stochastic minimax optimization."

We conduct numerical experiments on the following three hard instances:
*  a convex-concave bilinear function (in the folder `src/bilinear_func`).
*  a convex-concave $(\delta,\nu)$ function with additive Guassian noise (in the folder `src/delta_func`). 
*  a nonconvex-nonconcave function under the comonotonicity condition (in the folder `src/rho_func`). 
 
To reimplement the expeiments in our paper, please run 
```
exp_gnorm.m
```
in the corresponding folders.
