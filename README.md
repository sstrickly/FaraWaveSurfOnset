# FaraWaveSurfOnset
ON THE ONSET OF SURFACTANT-COVERED FARADAY WAVES

The code in this repository is linked to a research paper published in the IMA Journal of Applied Mathematics:
"The role of higher-order viscous and interfacial effects on the onset of surfactant-covered faraday waves"
by Stephen L. Strickland*a; Karen E. Danielsb; Michael Shearerc

The code in this repository can be run in either Matlab or Octave.

In this repository is a MasterExampleScript.m which shows how one would use the code.
The .m files named  Aon...  are functions which will calculate the onset acceleration for different situations.
  0Surf - no surfactant
  GSurf - with surfactant
  HInf - Infinite depth
  HFin - Finite depth
  PInf - Infinite Peclet number (no surface diffusion)
  OrdX - the highest order of the series expansion which the function will calcualte
    Ord1 will only calculate the first order term of the series expansion.
    Ord2 will calculate the sum of the first AND second order terms of the series expansion.
    etc.

The other files in this respository assist in calculating the onset acceleration.
  KelvinDispersionRelationSolver.m - numerically solves the finite-depth Kelvin dispersion relation, calculating the wave number given the frequency, depth, density, surface tension, etc.  This code uses a Levenberg–Marquardt algorithm.
  KelvinDispersionRelationTester.m - an auxiliarly function which assists the KelvinDispersionRelationSolver.m
  LMBuildYKelvin.m - an auxiliarly function which assists the KelvinDispersionRelationSolver.m


