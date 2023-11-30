******   MP_TOOLS code version 1.56 **** (C) Jiri Kulda, Grenoble/Prague 2017-2023  ******


Computing tools for visualisation and analysis of atomic scale models of crystalline solids. Fast generation of elastic and inelastic neutron and Xray diffuse scattering intensity maps in the momentum-energy space (Q,w) based the non-uniform fast Fourier transform of input data from molecular dynamics (MD), phasefield (Ferrodo) and lattice defect (Discus) simulations. Complementary pair distribution function (PDF) and vibrational density of states (DOS) computation. For more information and updated documentation consult https://mptools.fr


Disclaimer:
The MP_TOOLS code is (and probably remains forever) at a public_beta level of maturity and most probably contains bugs - please contribute to its development by (possibly) fixing them yourself and (anyways) by reporting them to the author CONTACT(AT)MPTOOLS.FR.

Code specifications:
FORTRAN code in free form (F90) with extensive use of F95 evolutions (pointers, dynamical variable allocation, vector operations) as implemented in the GNU GCC package (GFORTRAN).
  
Contributed elements: MP_TOOLS depend on the FINUFFT (non-uniform fast Fourier transform) 
  and on the FFTW (standard fast Fourier transform) libraries and use also the code of  
  Singleton's classical FFT implementation and the PGPLOT graphical output library.


1/ Contents of the MP_TOOLS v.1.56 package

  README.txt            this file
  LICENCE.txt           licence declarations
  compile_mp_156.txt    installation script
  COPYING.txt				    GNU GPL licence

  code [dir]            package code & related info
  contrib [dir]         licensing info on work of others
  doc  [dir]            brief descriptions of individual tools
  ref[dir]              reference data tables
  src  [dir]            the Fortran source files and the compile script
  examples [dir]        data examples (stubs), results & logs of related runs
  PGplot_install [dir]  installation comments and a debugged, high-res PNG driver


2/ Compiling and installing MP_TOOLS (linux, Mac OS X Darwin)

Hardware requirements:
  - for static models with < 10^6 atoms (cf. examples and tasks of corresponding size)
    a middle-class notebook with an intel i5_core2 processor (or equivalent) and 8Gb
    memory is sufficient (the time tags in the Examples logfiles correspond to such a
    configuration)
    
  - internal memory of 16Gb is a rock bottom minimum to explore large static supercells
    (10^7 atoms and more) and/or lattice dynamics with sequences of 10^2-10^3 frames; a
    workstation with 10-20 cores and 64-128Gb memory will bring comfort and boost overall
    efficiency    

Software requirements (all within the reach of standard link and execution pathes):
  - the GNU GCC suite including GFORTRAN (cf. the COMPILE_MP_154.TXT script)
  
  - the FFTW (Fastest Fourier Transform in the West) library (http://fftw.org/) 
  
  - the OpenMP code paralellisation library (https://www.openmp.org) 
    
  - the FINUFFT non-uniform FFT library (https://github.com/flatironinstitute/finufft)

  - the PGPLOT libraries (http://www.astro.caltech.edu/~tjp/pgplot/) installed in their
    standard location /usr/local/pgplot; you may wish to implement the corrected, high-res
    .PNG driver attached to the present distribution (useful for bitmap output)
  
Installing MP_TOOLS (unix, linux, Mac OS_X Darwin, may need the 'sudo' privileges):
  - copy the distribution archive into a suitable location (outside /usr/local) and
    uncompress it
  - create the '/usr/local/mp_tools' directory (and get privileges to access it)
  - check that the '/usr/local/bin' directory exists, otherwise create it and include it 
    into your execution path
  - set '<whatever_path>/MP_tools/code' as your current directory
  - execute the installation script 'source compile_mp_156.txt'; feel free to edit the
    script according to whether you prefer to overwrite the binaries of the previous
    version or not
  - in case of need (segmentation faults) add the '-fbounds-check' option to the compiler 
    calls to help localise the fault origins; similarly, in case of stable operation feel
    free to attempt more aggressive code optimisation by replacing the '-O' option by the
    '-O3' and by adding '-march=native' (or a specific architecture indentifier, cf. the 
    GCC documentation)
  
For the moment there is no makefile as the MP_TOOLS are a suite of short standalone programs; in case of modifications they can be recompiled individually copying just a single line from the compiling script to the command line. The executables will be automatically placed into /usr/local/bin, which is normally part of the execution path.

The EXAMPLES section contains data together with terminal logs and output files giving an idea of how to use MP_TOOLS, check DOCS/8_MP_TOOLS_USE_v154.txt for more details.