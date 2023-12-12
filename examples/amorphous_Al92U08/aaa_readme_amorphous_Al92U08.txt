
Use of MPTOOLS to treat RMC model data on an amorphous Al92U08 alloy

- this directory contains the MP_TOOLS parameter file, input data, command line logs and test results on amorphous AlU input data kindly provided by Jozef Bednarcik (Kosice/DESY), which served basis to the publication 
S Michalik, J Bednarcik et al., J. Phys.: Condens. Matter 22 (2010) 404209 (6pp)

- the MP_LBIN tool is first employed to convert the original ASCII data to the proprietary binary format of MP_TOOLS

- the MP_INSP tool permits to inspect the converted data

- the MP_PDF tool first accumulates directly the g(r) pair-distribution function by a projective MC sampling, avoiding the r^2 divergence of the radial distribution function, and procedes then by the conventional transformations towards G(r), S(Q) etc. 

- the MP_SQL tool calculates first the full 3D NUFFT (non-uniform fast Fourier transform, cf . references https://mptools.fr) and obtains then a spherically averaged total scattering intensity function I(Q) and transforms it to the other conventional scattering functions S(Q), F(Q) etc.; this procedure can prove efficient to compare high-quality model and experimental data over a restricted Q-range, which would not allow for a correct FT on the opposite direction (for the subsequent analysis of the model data in terms of g(r) etc. can be used MP_PDF)

- the ASCII data files (.txt) contain the MP_TOOLS output corresponding to color curves in the PGplot graphs; in principle the MP_PDF and MP_SQL tools have the option to import external data (absent in the present examples) to compare with the MP_TOOLS output, these data would appear as broad grey curves underlying the coloured MP_TOOLS output

- the included command line logs may prove useful as a guide to first test sessions, which should be run directly from this directory (the presensence of the Al_U_large.par file and the ./data subdirectory (with its *.dat contents) are crucial

- the description on the https://mptools.fr web site applies only partly (update foreseen by January 2024)

- the execution times correspond to a model consisting of 30000 atoms in box: 27600 Al + 2400 U; those shown in the logs are quite long, corresponding to a dated 3.1 GHz Dual-Core Intel Core i7 (2015) processor; up-to-date performance would be as follows:   

                    Apple Mini M1 (16Gb) & OMP 4 threads    

     ## NOTE: first run ever may take longer to settle the memory footprint
******************************************************************************************

                                   FINUFFT  AVERAGE
    MMP_SQL5 Q_range = 17    6Gb     3.46    4.38                                     
                             6Gb     3.50    4.43

    MMP_SQL5 Q_range = 20    8Gb     8.57    10.03    
                             8Gb     9.51    11.51
                             8Gb     9.44    11.58                             
                             8Gb     9.47    11.61
                             8Gb     9.65    11.65
                             8Gb     9.54    11.30
                             
    MMP_SQL5 Q_range = 25   15Gb    64.00    74.40    

     ## NOTE: memory requirements of MP_SQL explode with Q_range**3 
     ## NOTE:for short Q_ranges NUFFT overheads dominate

******************************************************************************************

    MMP_PDF5 n_h = 450                               
                          EFF1[%]    EFF2[%]  t[s]    EFF1 - hits to empty overlay cells
              J_ACC =2     67.4       4.1     57.9    EFF1 - hits outside PDF range (both EFF determined by data set & J_RAND)
                           67.4       4.1     53.0
                           67.4       4.1     49.0
                           67.4       4.1     48.7
              J_ACC =1     70.6      12.7     56.5 
                           70.2      12.7     52.7
                           70.2      12.7     56.6
                               
       
       n_h = 100 J_ACC =2  67.4       4.1     11.8    
                           67.4       4.1     10.8 
                           67.4       4.1     12.5 
                           67.4       4.1     10.8 
                              
       ## NOTE:  the choice of the overlay grid (PDF_PIX in the .PAR file) period needs
                 some expertise and testing: too small values reduce the EFF1 efficiency,
                 the too big ones leed to double occupancies of the cells (a few ones can
                 be tolerated(=skipped) by typing return in the dialogue)
                  
