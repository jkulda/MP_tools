
*****************************************************************************************

Use of MPTOOLS to treat single crystal MD data on a cubic BaZrO3 perovskite

*****************************************************************************************

- this directory contains the MP_TOOLS parameter files, input data, command line logs and test results on single crystal MD simulation data kindly provided by M. Pasciak (FZU Prague), in addition to the parameter file for the normal crystal lattice (BZO200d.par, CELL type data) used in the demos we give also the alternative BULK type files (BZO200b32. par, BZO200b1.par that would be of interest for systems with additional disorder)

- the MP_DBIN tool is first employed to convert the original ASCII data (only a stub is provided because of data volume limits) to the proprietary binary format of MP_TOOLS

- the MP_INSP tool permits to inspect the converted data

- the MP_PDF tool first accumulates directly the g(r) pair-distribution function by a projective MC sampling, avoiding the r^2 divergence of the radial distribution function, and procedes then by the conventional transformations towards G(r), S(Q) etc.; mind the difference of PDF_PIX between the cases of normal lattice (0.5) here and of "atoms in the box" (.02) in the amorphous example - in both case the period of the overlay grid is expressed as a fraction of the A_CELL lattice parameter

- the MP_DOS tool provides plots of the phonon density of states and its partials; a few hundred snapshots covering at least 10 ps of continuous trajectory in time are needed for reasonable quality; a sequence of 20 frames already provides a crude insight into the overall dynamics of the system with low-energy Ba vibrations quite neatly separated from the higher energy dynamics of the network of the ZrO3 octahedra

- the MP_SQL tool calculates the full 3D NUFFT (non-uniform fast Fourier transform, cf . references https://mptools.fr) and obtains then a spherically averaged total scattering intensity function I(Q) that can be compared to data from a scattering experiment; this procedure can prove efficient to compare high-quality model and experimental data over a restricted Q-range, which would not allow for a correct Q->R Fourier transform (for the subsequent analysis of the model data in terms of g(r) etc. the  MP_PDF can be used)

- the MP_SQOM tool calculates 2D NUFFT slices of the reciprocal lattice; a few of them can provide a crude insight in the total scattering distribution,; hundreds up to thousand are required for high quality total - S(Q) - and/or energy-resolved - S(Q,w) - scattering function maps; the output based on a sequence of 20 frames already gives a hint of what may be the difference between total (bzo200d_sq_0001_0020_07.png) and energy-resolved scattering function (bzo200d_sq_0001_0020_11.png)

IMPORTANT NOTE: 
    the PGPLOT library is incapable to produce BITMAP output (eg. using the PNG driver) meeting standards of 21st century and the POST_SCRIPT format is not suited to this purpose (trouble with filling the big pixels); two work-arounds are possible:
    1/ choose one of the postscript drivers and set the screen plot as large as possible via P_SIZE=13 (inches), make a window/screen shot and resize it afterwards using some viewer application
    2/ use the tabulated output in the corresponding .TXT file as input to some modern-age plotting tool (Matlab, Origin, ...)

- the ASCII data files (.txt) accompanying the graphics output contain tabulated MP_TOOLS results that can be used as input to more advanced plotting software

- the included command line logs may prove useful as a guide to first test sessions, which should be run directly from this directory (the presensence of the BZO200**.par files and the ./data subdirectory (with its *.dat contents) are crucial

- the description on the https://mptools.fr web site applies only partly (update foreseen by January 2024)

- the displayed execution times correspond to a model consisting of a 32^3 cell box on a dated 3.1 GHz Dual-Core Intel Core i7 (2015) processor
    

