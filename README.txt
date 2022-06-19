******   MP_TOOLS code suite **** (C) Jiri Kulda, Grenoble/Prague 2017-2022  ******


Computing tools for visualisation and analysis of atomic scale models of crystalline solids.
Fast generation of elastic and inelastic neutron and Xray diffuse scattering intensity maps 
in the momentum-energy space (Q,w) based the non-uniform fast Fourier transform of input data 
from molecular dynamics (MD), phasefield (Ferrodo) and lattice defect (Discus) simulations.
Complementary pair distribution function (PDF) and vibrational density of states (DOS) computation.


Disclaimer:
  The MP_TOOLS code is (and probably remains forever) at a public_beta level of maturity,
  most probably containing bugs - please contribute to its development by (possibly) fixing
  them yourself and (anyways) reporting them to the author at mp_tools(at)free.fr.

Code specifications:
  Fortran code in fixed form (F77) with extensive use of F95 evolutions (pointers, dynamical 
  variable allocation, vector operations) as implemented in the GNU GCC package (GFORTRAN).
  
Contributed elements: MP_TOOLS depend on the FINUFFT (non-uniform fast Fourier transform) 
  and on the FFTW (standard fast Fourier transform) libraries and use also the code of  
  Singleton's classical FFT implementation and the PGPLOT graphical output library.


1/ Contents of the MP_TOOLS v.1.5 package

  README.txt        this file
  LICENCE.txt       licence declarations

  code [dir]
    contrib [dir]     licensing info on work of others
    doc  [dir]        brief descriptions of individual tools
    ref[dir]          reference data tables
    src  [dir]        the Fortran source files and the compile script
    compile_mp_15.txt installation script
    COPYING.txt				GNU GPL licence

  examples [dir]    data & results of a short BaZrO3 MD run


2/ Compiling and installing MP_TOOLS (linux, Mac OS X Darwin)

Hardware requirements:
  - for static models with < 10^6 atom (attached examples and tasks of corresponding size)
    a middle class notebook with an intel i5_core2 processor (or equivalent) and 8Gb memory 
    is sufficient (the time tags in the Examples logfiles correspond to such a configuration)
    
  - internal memory of 16Gb opens the way to explore large static supercells and/or lattice 
    dynamics with sequences of 10^2-10^3 frames  is a rock bottom minimum and a workstation 
    with 10-20 Xeon cores and 64-128Gb memory will bring comfort and overall efficiency    

Software requirements (all within the reach of standard link and execution pathes):
  - the GNU GCC suite including GFORTRAN; place the following alias into your alias folder
    (logout and login back again to bring it into effect)
    alias ff='gfortran -funderscoring -O -ffixed-form -ffixed-line-length-132 -fimplicit-none –fPIC'
  
  - the FFTW (Fastest Fourier transform in the West) library (http://fftw.org/) 
  
  - the OpenMP code paralellisation library (https://www.openmp.org) 
    
  - the FINUFFT non-uniform FFT library (https://github.com/flatironinstitute/finufft)

  - the PGPLOT libraries (http://www.astro.caltech.edu/~tjp/pgplot/) installed in their standard 
    location /usr/local/pgplot
  
Installing MP_TOOLS (unix, linux, Mac OS_X Darwin) may need the 'sudo' privileages:
  - copy the distribution archive into a suitable location and uncompress it
  - create the '/usr/local/mp_tools' directory (and get privileges to access it)
  - check that the '/usr/local/bin' directory exists, otherwise create it and include it 
    into your execution path
  - set '<path>/MP_tools/code' as your current directory
  - execute the installation script 'source compile_mp_15.txt'
  - in case of need (segmentation faults) add the '-fbounds-check' option to the compiler 
    calls to help localise the fault origins; similarly, in case of stable operation feel free 
    to attempt more aggressive code optimisation by replacing the '-O' option by the '-O3' 
    and by adding '-march=native' (or a specific architecture indentifier, cf. the GCC 
    documentation)
  
For the moment there is no makefile as the MP_TOOLS are a suite of short standalone programs; 
in case of modifications they can be recompiled individually copying just a single line from
the compiling script to the command line. The executables will be automatically placed into 
/usr/local/bin, which is normally part of the execution path.


