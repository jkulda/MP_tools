
MP_TOOLS use the PGPLOT code to create graphics output. 

On Mac OS X (Darwin) systems the PGPLOT libraries have to be built locally from source code downloadable from various sites (cf. B).

For LINUX systems precompiled PGPLOT libraries are obtainable as software package, eg. via the APT system of UBUNTU:
      >apt install pgplot5

Nevertheless, in order to be able to record high-resolution bitmap graphics, it is necessary, for the moment, to perform the manual installation as well.


A/ The PNG (portable network graphics) driver

Historically, there have been some issues concerning its driver handling output of bitmap graphics in the .PNG file format, which is of central importance for recording scattering intensity maps in our case. Patches applicable to the PNDRIV.C source code can be found easily (in fact they concern a single line in the middle of the code) and the binary PGPLOT5 distributions already contain this correction.
      
Nevertheless, these corrected PNG drivers still use a low resolution of 85 ppi (points-per-inch) and produce bitmaps of only ≈ 680x680 pixels. The PNDRIV_CORRECTED_HIRES.C driver, contained in this distribution, resolves this issue by defining a default 300 ppi pitch yielding cartoons of 2048x2048 pixel size and introducing a new environment variable PGPLOT_PNG_PPI, which permits to adjust this default value according to actual needs by 
  >export PGPLOT_PNG_PPI=300   (in case of BASH shell)

Known issues:
1/ on smaller computer systems the default 300 ppi value may result in an only partially rendered image field - in such cases the resolution can be adapted by setting a lower value via the above mentioned instruction

2/ sometimes it may happen that the PNG file is generated, but not saved to the disk:
a new message on STDOUT warns when the PNG was not saved, usually repeating the operation once again resolves the problem


B/ Manual installation of PGPLOT libraries

The local (manual) installation is similar on any unix-like system and is described step-by-step in what follows:

0/ (prerequisite) have the GNU Compiler Collection (GCC) installed - eg. via APT on the the linux UBUNTU/DEBIAM systems or via HOMEBREW on Mac and use it to compile PGPLOT (for consistency reasons)

1/ get a recent PGPLOT distribution (5.2.2 or 5.3.1)

2/ unpack the tarball in any suitable place (<dir>; eg. /usr/local/src)

3/ replace its MAKEMAKE file by MAKEMAKE_MP.SAVE from this distribution (it contains the patches advertised over internet)
     >cp ./makemake_mp.save <dir>/pgplotsrc/makemake

4/ replace the PNG driver PNDRIV.C in the distribution by PNDRIV_CORRECTED_HIRES.C
     >cp ./pndriv_corrected_hires.c <dir>/pgplotsrc/drivers/pndriv.c

5/ inspect the list of system subdirectories (SYS_UBUNTU etc.) and their *.CONF files for a suitable one, using GFORTRAN and CC as compilers 
  
  a) for most LINUX environments the right choice should be SYS_UBUNTU (with GCC.CONF
  inside) from the present MP_tools distribution, provided you have installed the GNU GCC
  compiler suite on your system, copy this directory into the PGPLOT distribution
     >cp -r ./sys_ubuntu <dir>/pgplotsrc

  b) for DARWIN on Mac OS X systems the SYS_MACOSX contains configuration files adapted to
  Intel MACs (GFORTRAN_GCC_64.CONF) and to more resent Apple Silicon (M1,M2) systems .....  
  
6/ create the PGPLOT directory and make it your working directory
     >sudo mkdir /usr/local/pgplot
     >cd /usr/local/pgplot
     
7/ run from there the MAKEMAKE script with suitable arguments specifying your system
  configuration
      >sudo <dir>/pgplotsrc/makemake <dir>/pgplotsrc/sys_ubuntu gcc
  
8/ execute MAKE from your /usr/local/pgplot directory 
      >sudo make
     
  ... don't be scared by the possible WARNINGS! only ERRORS are FATAL ...
  
9/ now you are able to build MP_TOOLS using the compile script in the distribution, 
    -L/usr/local/pgplot   includes /usr/local/pgplot into the link path
    -lpgplot              indicates that libpglot* has to be looked for