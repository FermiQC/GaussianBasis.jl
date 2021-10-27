# GaussianBasis.jl
[![CI](https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml)

GaussianBasis offers a high-level wrap around the integral library [libcint](https://github.com/sunqm/libcint). Install it by simply runnig

```
pkg> add GaussianBasis
```

Everything should work automatically, but sometimes there might be some difficulties building `libcint`. The file deps/build.jl contains simple commands to clone and build this library, you might need to modify it to better suit your system. If you do, rerun the build step using `pkg> build GaussianBasis`. Please reach out if you encounter any problem.

Current features include:

- Standard basis set files (`gbs` format)
- Basis set parsing 
- One-electron integral (1e)
- Two-electron two-center integral (2e2c)
- Two-electrons three-center integral (2e3c)
- Two-electrons four-center integral (2e4c)

In all cases, a dense array is returned containing all integral values. For the 2e4c integral, a sparse version exists where only unique non-zero elents are returned, along with its indices. 

## TODO

1) Integral gradients are being constructed, currently we have some finite differences options (mainly used for debugging) and functional, but not optimal, gradients for 1e and 2e4c. Thus the current goals are

- Optimize gradients for 1e
- Optimize gradients for 2e4c
- Implement gradients for 2e2c and 2e3c
- Hessians

2) While we have light wraps around basic `libcint` functions that should allow you to do anything, it would be nice to have higher-abstracted functions for batched integral computations.
