# VibronicToolkit

Tools for vibronic Hamiltonians.
This package aims for correctness and readability rather than performance, so these tools are mainly useful for verifying other implementations.

Tested with Julia 1.0.


## Installation

```
pkg> add https://github.com/0/VibronicToolkit.jl.git
```

In order to run the driver scripts in `bin/`, you will also need to
```
pkg> add ArgParse
```


## Units

* The input energy parameters define the unit of energy `E`.
* The input frequency parameters are in units of `E` (i.e. `hbar = 1`).
* Positions are dimensionless, so the other input parameters are also in units of `E`.
* Reciprocal temperatures are in units of reciprocal energy (i.e. `k_B = 1`).
* Heat capacities are dimensionless (in units of `k_B`).


## Examples

To run the following examples, you should set the project (e.g. using `--project` or `JULIA_PROJECT`) to a Julia project that has the prerequisites installed.

* `julia bin/analytical.jl --conf examples/s2m2_uncoupled.json --beta 43.21`
* `julia bin/sos.jl --conf examples/s2m2_uncoupled.json --beta 43.21 --basis-size 30`
* `julia bin/trotter.jl --conf examples/s2m2_uncoupled.json --beta 43.21 --basis-size 30 --num-links 20`
* `julia bin/sampling.jl --conf examples/s2m2_uncoupled.json --beta 43.21 --num-links 20 --num-samples 1000000`
* `julia bin/sampling.jl --conf examples/s2m2_uncoupled.json --beta 43.21 --dbeta 1e-4 --num-links 20 --num-samples 1000000`


## Publications

The following publications contain data created using this package:

* Neil Raymond, Dmitri Iouchtchenko, Pierre-Nicholas Roy, and Marcel Nooijen. **A path integral methodology for obtaining thermodynamic properties of nonadiabatic systems using Gaussian mixture distributions.** The Journal of Chemical Physics 148, 194110 (2018). [doi:10.1063/1.5025058](https://aip.scitation.org/doi/abs/10.1063/1.5025058), [arXiv:1805.05971](https://arxiv.org/abs/1805.05971).


## Acknowledgements

Thanks to Neil Raymond for designing the parameter file format and helping to verify this implementation!


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
