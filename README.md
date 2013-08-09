WIM1d
=====

Code for 1d version of WIM

OcMod2013_paper_code has the code for Williams et al (2013a,b).

Please make new directories for new developments.

Things possibly on the agenda:
- make attenuation depend on Young's modulus as well
- set a definite time scale for the breaking to happen at
  *currently it happens at the advection/attenuation time step;
   possibly the behaviour with Courant number could change if we do this
  *this may also impact on the way we (effectively) remove energy from
   the waves if breaking occurs - ie we may need to do it explicitly.
- include generation of waves by wind
- numerical desciption of floe size distribution &/or evolution equation for it
- refreezing
