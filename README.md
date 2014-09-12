WIM1d
=====

Code for 1d version of WIM

1. OcMod2013-paper has the code for Williams et al (2013a,b).
   NB please don't modify this branch!!

2. Please make new branches for new developments/testing.
   (These can be merged with the master branch if they work out.)

3. test-youngs to allow attenuation to depend on Youn's modulus

Things possibly on the agenda:
- make attenuation depend on Young's modulus as well
- set a definite time scale for the breaking to happen at
  *currently it happens at the advection/attenuation time step;
   possibly the behaviour with Courant number could change if we do this
  *this may also impact on the way we (effectively) remove energy from
   the waves if breaking occurs - ie we may need to do it explicitly.
- include generation of waves by wind
- numerical description of floe size distribution (ie have floe size bins)  &/or evolution equation for it
- refreezing
