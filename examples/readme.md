---
layout: page
title: Examples   
---

The `examples` directory contains all input and output files related to an `MP_tools` runs on brief sequences of five snapshot dumps from more extended and MD simulations. Obviously, such short MD sequences do not provide any meaningful access to functionalities concerning inelastic scattering - these call for sequences of 10^2 to 10^4 MD snapshots (do not hesitate to contact the author for assistance).

The first example is part of a [`DL_POLY`](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx) run on a baryum titanate box consisting of 50 x 50 x 50 elementary cells (for details cf. Pasciak et al., PRL 120, 167601 (2018)), in the second case a [`LAMMPS`](https://www.lammps.org) run on a copper box of 48 x 48 x 48 cells. In the `LAMMPS` case the file numbering together with the elementary time step (microstep) of 0.0002 ps and the dumps after each 50 microsteps result in a time sequence with a 0.01 ps step.  

In the case of `LAMMPS` data the whole run can be reproduced, once the `MP_tools` are properly installed. As a guide may serve the `mp_tools_lammps_cu_48_log.txt` file containing a dump of the related command line dialogue. Because of data volume limitations the `DL_POLY` example contains only a stub of the trajectory ASCII dump file, yet the binary data again can be used to reproduce the other `MP_TOOLS` operations. Again the command line log `mp_tools_dl_poly_BTO500c_log.txt` is available for reference.

<b>In both cases a particular importance have the parameter (`.PAR`) files `BTO500c.par` and `cu_c48.par`, which can be used - after a careful examination and update of the options - as templates for any further work.</b>