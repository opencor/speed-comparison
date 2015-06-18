The idea behind this speed comparison is to show the difference in computational time between OpenCOR and MATLAB when used out-of-the-box, i.e. with whatever solver that comes with either of those packages.

Testing environment:
 - MacBook Pro (Retina, Mid 2012).
 - Processor: 2.6 GHz Intel Core i7.
 - Memory: 16 GB 1600 MHz DDR3.
 - Operating system: OS X Yosemite 10.10.3.

[OpenCOR](http://www.opencor.ws/):
 - Version: 0.4.1.
 - Solver: CVODE with its default settings, except for its ```Maximum step``` parameter, which is to be set to a model's stimulation duration, if needed.

[MATLAB](http://www.mathworks.com/products/matlab/):
 - Version: R2013a.
 - Solver: ode15s (i.e. a solver suitable for stiff problems and which has low to medium order of accuracy) with both its ```RelTol``` and ```AbsTol``` parameters set to ```1e-7``` and its ```MaxStep``` parameter set to the stimulation duration, if needed.

Tested models (with its simulation duration in seconds and, if needed, its stimulation duration in milliseconds):
 - [Bondarenko et al. 2004](http://models.cellml.org/e/41) (```10 s``` / ```0.5 ms```);
 - [Courtemanche et al. 1998](http://models.cellml.org/exposure/0e03bbe01606be5811691f9d5de10b65) (```100 s``` / ```2 ms```; note: the value of ```membrane.stim_end``` was increased so as to get action potentials for the duration of the simulation);
 - [Faber & Rudy 2000](http://models.cellml.org/exposure/55643f2114a2a463ada007deb9fc3913) (```50 s``` / ```2 ms```);
 - [Garny et al. 2003](http://models.cellml.org/exposure/d71105df45dd7030b3c99b2b1e95b8c0) (```100 s```);
 - [Luo & Rudy 1991](http://models.cellml.org/exposure/2d2ce7737b42a4f72d6bf8b67f6eb5a2) (```200 s``` / ```2 ms```; note: the value of ```membrane.stim_end``` was increased so as to get action potentials for the duration of the simulation);
 - [Noble 1962](http://models.cellml.org/exposure/812eeafbc8ebe97bef435340c80cfcce) (```1,000 s```);
 - [Noble et al. 1998](http://models.cellml.org/exposure/a40c4434423c0436e2789a2d457b7ab2) (```100 s``` / ```3 ms```);
 - [Nygren et al. 1998](http://models.cellml.org/exposure/ad761ce160f3b4077bbae7a004c229e3) (```100 s``` / ```6 ms```); and
 - [ten Tusscher & Panfilov 2006](http://models.cellml.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d) (```100 s``` / ```1 ms```).
