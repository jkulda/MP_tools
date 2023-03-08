  
      program mp_pdf55

C *************************************************************************************
C *****
C ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   calculates the pair distribution functions (PDF) for simulated supercell data
C *****
C**********					Copyright (C) 2022  Jiri Kulda, Grenoble/Prague          **********
C**	
C** This file is part MP_TOOLS developed and maintained by Jiri Kulda <jkulda@free.fr>
C**
C**	MP_TOOLS are free software: you can use it, redistribute it and/or modify it 
C**	under the terms of the GNU General Public License as published by the Free Software 
C**	Foundation, either version 3 of the License, or (at your option) any later version, 
C**	for details see <https://www.gnu.org/licenses/>
C**
C**	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
C**	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
C**	See the GNU General Public License for more details.
C**
C ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C ***** Ver. 1.1 The correlation is asymmetric: N_HSUM random cell against the complete supercell
C ***** Ver. 1.2 the whole supercell is read into memory once for ever
C ***** Ver. 1.21 the correlation is symmetric N_HSUM random cells against N_HSUM random cells
C *****
C ***** Ver. 1.31 between two delayed snapshot frame uses the 1.21 symmetric algorithm
C ***** Ver. 1.31 now the delay may be 0, so 1.31 replaces also ver. 1.21 and will override it
C *****
C ***** Ver. 1.32 the IDOM dialog moved after reading first input file
C ***** Ver. 1.32 the TIME_LOOP replaces the MASTER_LOOP
C ***** Ver. 1.32 normalization of correlation functions by true number of active cells n_norm = 45**3
C ***** 
C ***** Ver. 1.33 the reference snapshot and the first one in the time-series defined independently
C ***** Ver. 1.33 hence instantaneous (lag 0) correlations can be calculated explicitly
C ***** 
C ***** Ver. 1.4 	- flexible cell structure
C ***** 				 	- use of auxiliary <file_title.par> file (correlation pairs etc.)
C ***** Ver. 1.41	- corrected inconsistency in output <rdf>.txt file name
C *****						- correct histogram identification by j_hist_type
C ***** 
C ***** Ver. 1.42	- lighter & parallelied version (no correlation histograms, no IDOM)
C *****						- integrates the PGPLOT graphics of mp_rplot.f
C ***** 					- only PDF
C ***** 					- for correlation histograms (refer to mp_corr)
C ***** 					- modified .PAR file format, consequent .LOG
C ***** 
C ***** Ver. 1.50	- start of new series (identical to mp_pdf 1.42) as MP_PDF
C ***** 
C ***** Ver. 1.51	- all data arrays allocatable, no predefined array size limits
C ***** 					- supercell format in .PAR
C *****						- number of correlation pairs read from .PAR (don't use ver. 1.50!)
C *****
C ***** Ver 1.53  - record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****           - CELL and QUICK data type (BULK not eligible for PDF)
C *****						- secondary sampling volume clipped for highly anisotropic or 2D supercells
C *****
C ***** 
C ***** Ver 1.54  - NAMELIST I/O for .PAR files and .DAT headers implemented 
C *****           - PNG driver for PGPLOT included
C *****           
C ***** Ver 1.55  - completely rewritten MC accumulation part  
C *****           - consistent MC results between OMP and non-OMP runs
C *****           - estimate of MC statistical accuracy
C *****           - automatic accumulation of the complete correlation landscape
C *****           - pseudoatoms to represent group correlations
C *****           - overlay with adjustable partial PDF scales
C ***** 
C ***** file name convention:
C *****     <master_file>_n<snapshot number>.dat
C *****     example: H05_10K_n0001.dat from H05_10K.txt
C *****
C ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
C ***** record numbers are directly related to cell positions and atom types
C *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
C *****						1 going for the header record
C *****
C ***** atom positions are converted from reduced lattice coordinates (x) to real space distances
C ***** 
  
  		use omp_lib
  		
CC  		use mp_nopgplot					! uncomment when not able to use PGPLOT, compile and link together with the module mp_nopgplot.f

      integer,parameter :: l_rec  =  1024		    !record length in real(4)
      real,parameter    :: pi = 3.14159

      character(4),allocatable  ::	at_name_par(:),at_label(:),at_name_ext(:),at_name_pseudo(:),at_name_plot(:)
      character(10),allocatable ::	corr_name(:)

      integer,allocatable ::  at_ind_file(:,:,:),at_mask(:),j_corr(:,:),mult_corr(:),ind_pseudo(:,:),numbers(:) 

			real,allocatable ::  at_pos_file(:,:,:),rdf_tot(:),rdf_tot_plot(:),r(:),rdf_p2(:,:,:),rdf_p2_n(:,:,:,:),rdf_p2_plot(:,:,:),rdf_err(:)
      real,allocatable ::  b_coh(:),at_weight(:),at_weight_matrix(:,:),at_base(:,:),at_pos(:,:),at_pos2(:,:)
	  
			integer ::       i_rdf,n_pdf,j_acc,n_int,n_int2,n_pseudo,i_ref,j_pdf,j_rand
      real    ::       rdf_dist,rdf_range,rdf_sum,rdf_tot_sum,rdf_norm,at_weight_sum
      real 		:: 			 arg,c_min,c_max,c_max2,c1,c2,r_start,r_end,pdf_step,c_smooth,f_smooth(51)
      
      character(4)   :: c_int(2),c_fil(2),version,head,atom,ps_out(2),size_out(2)
      character(10)  :: at_weight_scheme(3),pg_out,string,section,c_date,c_time,c_zone,c_nfile_min,c_nfile,c_jfile
      character(16)  :: sim_type_par,data_type,string16,filter_name
      character(40)  :: subst_name,file_master,file_inp,file_out,time_stamp,int_mode,x_file_name
      character(60)  :: file_dat,file_dat_t0,file_res,file_ps,file_log,line,masks
      character(128) :: cwd_path,plot_header,plot_title
      character(l_rec):: header_record
      
			logical ::  nml_in,found,found_txt,found_ps,t_single

      integer(8) ::  n_in,n_out,n_inbox,n_outbox,n_pdf_range(3),n_mc,n_mc_max,n_norm
      integer ::  j_proc,proc_num,proc_num_in,thread_num,thread_num_max,m,j_name,jm,j1m,j2m,j_mode
      integer ::  i_start,i_end,ind_part(2,4),n_step,nfile_step,n_hsum,n_h,j_head_in,hist_ind(3),j_verb,n_tot,n_head
      integer ::  i,j,k,ii,jj,ind,ind2,jat,i_rec,ind_rec,nrec,ncell,nrow,nlayer,nsuper,nfile,nfile_t0,nfile_min,nfile_max,jfile
      integer ::  j1,j2,j_plane,j_grid,j_logsc,j_ps,j_out,j_txt,i_seed,i_r1,i_r11,i_r12,i_time(8),j_sr,j_pb,j_sr2,j_pb2
      integer ::  at_no,at_ind(3),at_ind2(3),d_ind(3),j_dom,j_base(3),n_plot,n_smooth,n_smooth_fwhm,n_part,n_part_max,n_atom_tot,n_pseudo_max
      integer ::  j_atom,j_weight,j_xray,j_edit,j_mask,ifile,ifile0,jfil_step,jint,jint2,n_atom,sc_c2,sc_c1,ier,ios
      
      real :: t_dump,filter_fwhm,tt0,tt,b_sum,rand1,rand2,sc_r,at_displ,p_size
      integer :: rand1_seed(8),rand2_seed(8),seed_size

      real :: diff_pos(3),x_plot,y_plot,part_scale(4)

C **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
      character(4),allocatable :: at_name_out(:)
      integer(4),allocatable   :: at_ind_in(:),nsuper_r(:)	

      real(4),allocatable ::	at_pos_in(:),at_occup_r(:)

      character(16)  :: sim_type,dat_type,input_method,file_par,dat_source
      integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,n_tot_in,idum,iran0
      real(4)        :: rec_zero(l_rec),t_ms,t_step,t0,t1,a_par(3),angle(3),a_par_pdf(3),temp

      namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,
     1    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
      namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
     
      namelist /mp_gen/ j_verb,j_proc       
      namelist /mp_out/ j_weight,j_logsc,j_ps,j_txt,p_size,j_grid,pg_out       
      										!general rule: namelists of tools should only contain their local parameters
                          !what is of global interest they should pass into data_header
			namelist /mp_pdf/ n_pdf,pdf_step,a_par_pdf,j_rand,n_mc,j_acc,j_pdf,n_part_max,n_pseudo_max,j_out
			
C **** PGPLOT stuff
      INTEGER :: PGOPEN,j_xserv


CCC ********************* Test of RANDOM_NUMBER generator (start) *******************************      
CCC        
CCC        - remove the CC to check possible local and processor-dependent effects
CCC        - default seed size 8, can be changed by
CCC        
CCC     CALL RANDOM_SEED(SIZE=K)               !Puts size of seed into K
CCC     CALL RANDOM_SEED(PUT = SEED (1 : K))   ! Define a new seed
CCC     CALL RANDOM_SEED(GET = old_seed)       ! Read the current seed
CC
CC      integer ::  old_seed(8),new_seed(8)		
CC      	
CC			old_seed = 0
CC			new_seed = (/1,2,3,4,5,6,7,8/)
CC			CALL RANDOM_SEED         !Processor initialization
CC      do i=1,4
CC      CALL RANDOM_SEED(GET = old_seed)  ! Read current seed
CC      call random_number(rand1)
CC      write(*,*) 'default',rand1,old_seed
CC      enddo
CC
CC      CALL RANDOM_SEED(PUT = new_seed) ! Define seed
CC      do i=1,4
CC      CALL RANDOM_SEED(GET = old_seed)  ! Read current seed
CC      call random_number(rand1)
CC      write(*,*) 'new',rand1,old_seed
CC      enddo
CC
CC      CALL RANDOM_SEED(PUT = new_seed) ! check the reproducibility
CC      do i=1,4
CC      CALL RANDOM_SEED(GET = old_seed)  ! Read current seed
CC      call random_number(rand1)
CC      write(*,*) 'new_again',rand1,old_seed
CC      enddo
CC 
CC C *** Output of test:
CC 
CC default  0.774232209      -328864355 -1812850565  1795335799  -423129627   258875061   746672867 -1013455889  1109822438
CC default  0.499739051     -1772520199  1041085153  -345388489  -275392182 -1581757171   246595309  -441046729 -1125993601
CC default  0.551364720       512484964  1716921112  -363648824 -1071663748 -1457211347  -619233579   394076467  -840249998
CC default  0.671572745      -592193939 -1415474521  -203200136 -1628592313  1380558838  1345367254  1710769802   823960254
CC new  0.471070886               1           2           3           4           5           6           7           8
CC new  0.117344737       -16064360  1080567055  2124207785  -289480891   105059873 -2032997881  -987088397  1102541325
CC new  0.357547939     -1494515791  2126487000  -379374510 -1783270911 -1145639077 -1358320318  -112593222  -960835836
CC new  0.318134785       957948723  -227147975 -1341437354  -415669684 -1377879149  2060778058   557735840  1455132063
CC new_again  0.471070886               1           2           3           4           5           6           7           8
CC new_again  0.117344737       -16064360  1080567055  2124207785  -289480891   105059873 -2032997881  -987088397  1102541325
CC new_again  0.357547939     -1494515791  2126487000  -379374510 -1783270911 -1145639077 -1358320318  -112593222  -960835836
CC new_again  0.318134785       957948723  -227147975 -1341437354  -415669684 -1377879149  2060778058   557735840  1455132063
CC     
CCC ********************* Test of RANDOM_NUMBER generator (end) *******************************      




C ********************* Initialization *******************************      

			write(*,*) '*** Program MP_PDF 1.55 ** Copyright (C) Jiri Kulda (2019,2021,2022) ***'
			write(*,*)
			
C ********************* Get a time stamp and open a .LOG file *******************************
      call getcwd(cwd_path)
	    call date_and_time(c_date,c_time,c_zone,i_time)
	    write(time_stamp,110) c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
110		format(a,'/',a,'/',a,' ',a,':',a,':',a)	    
	    write(file_log,111) trim(c_date)
111   format('mp_tools_',a,'.log')
			inquire(file=file_log,exist=found_txt)
			if(found_txt)then
		    open (9,file=file_log,position='append')
		  else
		    open (9,file=file_log)
		  endif

			write(9,*) 
			write(9,*) 
			write(9,*) trim(time_stamp),'  MP_PDF 1.55  ',trim(cwd_path)
			write(9,*) 
      
C *********************  OpenMP initialization start  *******************************      
C
C    cf. line ≈410 (using input from .PAR)
C
C ********************* OpenMP initialization end *******************************      

C *** other initialisations
			c_int = (/'GS','MC'/)
			c_fil = (/'g','m'/)
			ps_out = (/'OFF ','ON  '/)			!PGPLOT
			size_out = (/'S   ','XXL '/)			!PGPLOT
			
C *** Generate data file access
			write(*,*) 'Input data file_master: '
			read(*,*) file_master 
						
      write(*,*) 'Read data files number min, max (0 0 no numbers, single file): '
      read(*,*)   nfile_min,nfile_max
      nfile_step = 1
			t_single = nfile_min==0.and.nfile_max==0

			if(t_single)then
				nfile_min = 1		!for conformity with loop control conventions
				nfile_max = 1
			endif
			
      nfile = ((nfile_max-nfile_min)/nfile_step)+1
      
			if(t_single)then
				write(file_dat_t0,'("./data/",a,".dat")') trim(file_master)
			else
        if(nfile_min<=9999) then
          write(file_dat_t0,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min
        elseif(nfile_min>=10000) then
          write(string,'(i8)') nfile_min
          file_dat_t0 = './data/'//trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
        endif
			endif

	    open (1,file=file_dat_t0,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_dat_t0),' not found! Stop execution.'
				stop
			endif
			
C *** Read the header record 

      i_rec = 1   			
      read(1,rec=i_rec) header_record
      head = header_record(1:4)
      call up_case(head)
      if(head == 'MP_T') then      !new structure with namelist
        nml_in = .true.
				read(1,rec=i_rec) dat_source,version,string16
				read(string16,*) n_head
        write(*,*)		'Input data:  ',dat_source,version,n_head
        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_1)	
        t0 = t_dump
CC        write(*,nml=data_header_1)
      elseif(head.eq.'TIME'.or.head.eq.'STAT') then                                  !old w/o structure
        write(*,*)		'Input data:  ','old header format'
        nml_in = .false.
       read(1,rec=1) 
     1		  sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,n_cond,n_rec					
         n_head = 1
         call up_case(sim_type)
         input_method = 'CELL'
				 if(index(sim_type,'BULK')/=0) input_method = 'BULK'
      else
      	write(*,*) 'header record wrong or missing'
      	write(*,*) trim(header_record)
      	stop
      endif 
      
      if(input_method == 'BULK') then
        write(*,*) 'BULK data are not suitable for this PDF algorithm, works with CELL data only.'
        stop
      endif

      allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom))
			allocate(at_label(n_atom),at_name_par(n_atom),at_name_ext(n_atom))
			allocate(at_base(n_atom,3),at_weight(n_atom),at_weight_matrix(n_atom,n_atom),at_mask(n_atom),rdf_err(n_atom))
			allocate(b_coh(n_atom),SOURCE=.0)												!we need this to read .par
      at_mask = 1

      if(nml_in) then      !new structure with namelist
        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_2)			
      else                                  !old w/o structure
			  read(1,rec=1) 
     1		sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell_out,n_cond,n_rec,n_tot,			
     2    at_name_out,at_occup_r,nsuper_r																		
     		if(head.eq.'TIME') sim_type = 'TIMESTEP'
     		if (head.eq.'STAT') sim_type = 'STATIC'
      endif 
      close(1)
      
      
C **** Read the auxiliary file <file_par.par> with structure parameters, atom names and further info
      write(*,*) 'Parameter file name (.par to be added) (confirm or type other name): ', file_par
      read(*,*) file_par
      file_inp = trim(file_par)//'.par'

      open(4,file=file_inp,action='read',status ='old',iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_inp),' not found! Stop execution.'
				stop
			endif

      write(9,*) 'Read the parameter file:  ',trim(file_inp)

      j_verb = 0
      j_proc = 0
      read(4,nml=mp_gen)

      j_ps = 0               !defaults for &mp_out
      j_txt = 0
      j_weight = 1
      j_logsc = 0
      j_grid = 0
      pg_out = 'ps'
      p_size = 7.
      rewind(4)
      read(4,nml=mp_out)

      at_weight_scheme(1) = 'Uniform'
			at_weight_scheme(2) = 'Neutron'
			at_weight_scheme(3) = 'Xray'

      n_pdf = 1024
      pdf_step = 0.06
      a_par_pdf = 1.
      n_part_max = 4         !defaults for &mp_pdf
      n_pseudo_max = 4
      j_out = 0
      j_rand = 1
			n_mc = 1000000					!unit for counting MC cycles (n_h)
      rewind(4)
      read(4,nml=mp_pdf)

      allocate(ind_pseudo(n_atom,n_pseudo_max),at_name_pseudo(n_pseudo_max),ind_part(2,n_part_max))
      ind_pseudo = 0
      ind_part = 0
      at_name_pseudo = ''

C *** Read the atom positions       
      rewind(4)
      section = 'atoms'
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
          stop
        endif
        if(string(1:5).eq.section) exit	!find the mp_simple part of the .par file
      enddo
			do j=1,n_atom
        read(4,*) at_label(j),at_name_par(j),arg,arg,arg,arg,at_name_ext(j)	!at_base & conc come from the header, for BULK the at_base and conc are not significant
        if(at_name_ext(j)=='n'.or.at_name_ext(j)=='N') at_name_ext(j)=''
			enddo

C *** Read the partial PDF definitions       
      if(n_part_max>0) then
        rewind(4)
        section = 'partial'
        do
          read(4,'(a)',iostat=ios) string
          if(ios/=0) then
            write(*,*) 'Section title: PARTIAL_PDF  not found (can be added in dialogue)'    !n_part,n_pseudo
            n_part = 0
            found = .false.
            exit
          endif
          found = (string(1:7).eq.section) 
          if(found) exit	!found the 'partial_pdf' part of the .par file
        enddo
        
        if(found) then
          do j=1,n_part_max
            read(4,*,iostat=ios) string,ind_part(:,j)	!indices of partial PDFs to be displayed reading until end of the list
            if(ios/=0) then
              n_part = j-1
              exit
            endif
            n_part = j
          enddo
        endif
      endif
            
C *** Read the pseudo atom definitions       
      if(n_pseudo_max>0) then
        rewind(4)
        section = 'pseudo'
        do
          read(4,'(a)',iostat=ios) string
          if(ios/=0) then
            write(*,*) 'Section title: PSEUDO_ATOMS not found (can be added in dialogue)'
            n_pseudo = 0
            found = .false.
            exit
          endif
          found = (string(1:6)==section)
          if(found) exit	                  !found the 'pseudo_atom' part of the .par file
        enddo
        
        if(found) then
          do j=1,n_pseudo_max
            read(4,*,iostat=ios) at_name_pseudo(j),ind_pseudo(1:n_atom,j)	              ! pseudo_atom name and indeces
            if(ios/=0) then
              n_pseudo = j-1
              exit
            endif
            n_pseudo = j
          enddo
        endif
      endif

			close(4)
			
C *** Check the .PAR input       
			do j=1,n_atom
				if(at_name_par(j)/=at_name_out(j)) then
					write(*,*) 'Not-matching  atom names in .PAR and .DAT: ',j,at_name_par(j),at_name_out(j)
					write(*,*) 'Prefer .DAT? (1/0)'
					read(*,*) ii
					if(ii==1) at_name_par = at_name_out
					exit 
				endif
			enddo			

C *** read neutron scattering lengths (always)
      open(4,file='neutron_xs.txt',action='read',status ='old',iostat=ios)
      if(ios.ne.0) then
        open(4,file='/usr/local/mp_tools/ref/neutron_xs.txt',action='read',status ='old',iostat=ios)
        if(ios.ne.0) then
          do
            write(*,*) 'File neutron_xs.txt not found, type in valid access path/filename'
            read(*,'(a)') x_file_name
            open(4,file=trim(x_file_name),action='read',status ='old',iostat=ios)
            if(ios==0) exit
            write(*,*) 'File',trim(x_file_name),' not found, try again ...'
          enddo
        endif
      endif
      atom_loop: do j=1,n_atom
        rewind(4)
        do i=1,210
          read(4,*) atom
          if(atom==trim(at_label(j))) then
            backspace(4)
            read(4,*) atom,b_coh(j)
            cycle atom_loop
          endif
        enddo
        write(*,*) 'b_coh for ',trim(at_label(j)),' not found,'
        write(*,*) 'check your spelling and the neutron_xs.txt table; use unit weights'
      enddo atom_loop
      close(4)
      b_coh = .1*b_coh  !convert b_coh from FM to 10^12 cm
CC      write(*,*)
CC      write(*,*) 'Neutron scattering lengths:'
CC      do j=1,n_atom
CC        write(*,*) at_label(j),b_coh(j)
CC      enddo			

C *** read Xray formfactor parameters 
C     open(4,file='xray_ff.txt',action='read',status ='old',iostat=ios)
C     if(ios.ne.0) then
C       open(4,file='/usr/local/mp_tools/ref/xray_ff.txt',action='read',status ='old',iostat=ios)
C       if(ios.ne.0) then
C         do
C           write(*,*) 'File xray_ff.txt not found, type valid access path/filename'
C           read(*,'(a)') x_file_name
C           open(4,file=trim(x_file_name),action='read',status ='old',iostat=ios)
C           if(ios==0) exit
C           write(*,*) 'File',trim(x_file_name),' not found, try again ...'
C         enddo
C       endif
C     endif
C     atom_loop: do j=1,n_atom
C       rewind(4)
C       do i=1,210
C         read(4,*) atom
C         if(atom==trim(at_label(j))//trim(at_name_ext(j))) then
C           backspace(4)
C           read(4,*) atom,x_ffpar(j,1:9)
C           cycle atom_loop
C         endif
C       enddo
C       write(*,*) 'Xray formfactor for ',trim(at_label(j))//trim(at_name_ext(j)),' not found,'
C       write(*,*) 'check your spelling and the neutron_xs.txt table; use unit weights'
C     enddo atom_loop
C     close(4)
C

C *** write overview of atom data
			write(*,*)
      write(*,*) 'Substance name: ',subst_name	  
			write(*,*) 'Atoms from ',trim(file_inp)
      do j=1,n_atom
        write(*,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j)
        write(9,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j)
      enddo
      b_sum = sum(at_occup_r(1:n_atom)*b_coh(1:n_atom))

			write(*,*) 

C **** use the .PAR value of a_par in case of input in LATTICE units      
      if(a_par(1)*a_par(2)*a_par(3)==1.and.a_par_pdf(1)*a_par_pdf(2)*a_par_pdf(3)/=0.) a_par = a_par_pdf
      write(*,*) 'Setting a_par to',a_par
      
C *** define references for lattice search	
			nrow = n_row(1)
      nlayer = n_row(1)*n_row(2) 						! 48**2 gives   2304
      nsuper = n_row(1)*n_row(2)*n_row(3) 	! 48**3 gives 110592  
      n_tot = n_atom*nsuper   
CC    	write(*,*) 'n_row',n_row
      rdf_range = n_pdf*pdf_step
      do j=1,3
      	n_pdf_range(j) = min(nint(rdf_range/a_par(j)),n_row(j))		!clip the integration range for strongly anisotropic supercells
      enddo
      n_norm = 8*n_pdf_range(1)*n_pdf_range(2)*n_pdf_range(3)		!number of unit cells within the PDF reach from R1
      n_mc_max = nsuper*(n_norm-1) !number of independent site correlation pairs
CC      write(*,*) 'n_pdf_range,n_norm,n_mc_max',n_pdf_range,n_norm,n_mc_max


C *********************  OpenMP initialization start  *******************************      
C
C     write(*,*) 'proc_num_in?'
C     read(*,*) proc_num_in
			proc_num_in = j_proc
C		write (*,*) 'OMP proc_num_in          = ',proc_num_in
			thread_num_max = omp_get_max_threads( )								!this gives maximum number of threads available (limited by concurrent tasks??)
			proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
			if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
			call omp_set_num_threads(proc_num_in)
			thread_num = omp_get_num_threads( )				!this gives threads available at this moment: here always = 1 as we are outside of a PARALLEL range


 			if(proc_num_in.gt.1) then			
				write(*,*) 
CC				write (*,*) 'OMP processors available = ', proc_num
				write(*,*) 'OMP threads max          = ', thread_num_max
				write(*,*) 'OMP processes requested  = ', proc_num_in										
C			write(*,*) 'OMP processes allocated 		= ', proc_num
C			write(*,*) 'OMP threads allocated       = ', thread_num
				write(9,*) 
				write(9,*) 'OMP threads max          = ', thread_num_max
				write(9,*) 'OMP processes requested  = ', proc_num_in
				write(9,*) 'OMP processes allocated  = ', proc_num
				write(9,*) 'OMP threads allocated    = ', thread_num
			else
				write(9,*) 'OpenMP not in use'
				write(*,*) 'OpenMP not in use'
			endif
			write(*,*) 
			write(9,*) 
C
C ********************* OpenMP initialization end *******************************      


C *********************      Initial info output   ************************
      
CC *** for the moment the MC integration strategy is imposed!

			j_acc = 2
			
C *** initialise the random_number initialisation
      call random_seed(size=seed_size)               !Puts size of seed into K
      allocate(numbers(seed_size))
      if(j_verb==1) then
        call random_seed(get=numbers)               !Gets actual seeds 
        write(*,*) 'random_number seed size:',seed_size
        write(*,*) 'seeds:',numbers
        write(*,*) 
      endif

C *** initialise the MC integration
			n_h = n_mc_max/1000000
			write(string,"(i8)") n_h
			write(*,*) 'MC sampling pairs per frame ([x 10^6], ',trim(adjustl(string)),' max)'
			read(*,*)   n_h
			n_hsum = n_h*n_mc !number of MC cycles per snapshot
			if(n_hsum.gt.n_mc_max) then
				write(*,*) 'WARNING: n_hsum exceeds max number of cell pairs ',n_mc_max
				write(*,*) 'input reduced n_h (<',n_mc_max/n_mc,'):'
				read(*,*)   n_h
				n_hsum = n_h*n_mc !number of MC cycles per snapshot
			endif
			n_int = n_hsum	! for MC cycle over n_hsum randomly chosen cell pairs
			n_int2 = 1      ! each MC step deals with a single R1-R2 site pair
			call random_seed	!initiate the random number sequence for MC sampling
			int_mode = 'Monte Carlo'
			write(*,*) 'Monte Carlo integration over',n_int*n_int2,'cell pairs'
			
C *** Allocate and clear the histogram arrays for accumulation across several snapshots	       
			allocate (at_pos_file(4,nsuper,n_atom),at_ind_file(4,nsuper,n_atom),at_pos_in(4*n_tot),at_ind_in(4*n_tot))
			allocate(at_pos(n_atom,3),at_pos2(n_atom,3))
			allocate(r(n_pdf),rdf_tot(n_pdf),rdf_p2(n_atom,n_atom,n_pdf),rdf_p2_n(n_atom,n_atom,n_pdf,n_h))

      rdf_p2 = .0
      rdf_p2_n = .0

      call cpu_time(t0)
      
      if(t_single)then
      	nfile_min = 1
      	nfile_max = 1
      	nfile_step = 1
      endif
      
C *************  the file loop: cycle over the tt files to augment statistics
      CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)

      file_loop: do ifile=nfile_min,nfile_max,nfile_step									

C ***  open the t0 file (binary MD snapshot file)
			if(t_single)then
				write(file_dat,'("./data/",a,".dat")') trim(file_master)
			else
        if(nfile_min<=9999) then
          write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min+(ifile-nfile_min)*nfile_step
        elseif(nfile_min>=10000) then
          write(string,'(i8)') nfile_min+(ifile-nfile_min)*nfile_step
          file_dat = './data/'//trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
        endif
			endif

      CALL SYSTEM_CLOCK (COUNT = sc_c2)
 			write(*,*)
			write(*,*)'input: ',trim(file_dat)
			write(9,*)'input: ',trim(file_dat)

		  open(1,file=file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_dat),' not opened! IOS =',ios
				write(*,*) 'Skip(1), stop execution(0)?'
				read(*,*) jj
				if(jj==1) exit file_loop
				if(jj==0) stop
			endif

			i_rec = n_head	
			do j=1,n_rec-1
				i_rec = i_rec+1
				read(1,rec=i_rec) at_ind_in((j-1)*l_rec+1:j*l_rec)			
			enddo	
			i_rec = i_rec+1
			read(1,rec=i_rec) at_ind_in((n_rec-1)*l_rec+1:4*n_tot)			
			
			do j=1,n_rec-1
				i_rec = i_rec+1
				read(1,rec=i_rec) at_pos_in((j-1)*l_rec+1:j*l_rec)			
			enddo	
			i_rec = i_rec+1
			read(1,rec=i_rec) at_pos_in((n_rec-1)*l_rec+1:4*n_tot)	
			
			at_ind_file(:,:,:) = reshape(at_ind_in,(/4,nsuper,n_atom/))
			
			at_pos_file(:,:,:) = reshape(at_pos_in,(/4,nsuper,n_atom/))
			
			close(1)
			
			write(*,*) 'Accumulating the PDFs ...'     

C *************  the integration loop: cycle over the site pairs to accumulate the PDFs ******************

			n_inbox = 0
			n_outbox = 0
			n_in = 0
			n_out = 0

                                            !if j_rand=0 do nothing, the system will supply unique, machine dependent seeds each time this code runs
      if(j_rand>1) then                     !if j_rand>1 generate a seed for later reference & numerical reproducibility checks
        idum = j_rand
        numbers(1) = iran0(idum)            !a dry call to initialise ran0
        do i=1,seed_size
          numbers(i) = iran0(idum)          !use a trivial random number generator to produce the seeds (they could even be all the same small ones, but ...)
        enddo
        call random_seed(put=numbers)       !Produce a seed to start, this permits to reproduce exactly the same results on the same system 
        do i=1,5
          call random_number(r(i))
        enddo
        if(j_verb==1) then
          write(*,*) 'Random_seed:',numbers
          write(*,*) 'Reference 1st 5 random numbers:',(r(i),i=1,5)
        endif
      endif


C *** MonteCarlo integration

!$omp parallel shared(at_pos_file,rdf_p2_n,a_par,pdf_step,nsuper,nlayer,n_row,n_atom,n_pdf_range,n_h,j_rand) private(at_pos,at_pos2,diff_pos,at_ind,at_ind2,d_ind,j_base,ind2,i_r1,i_rdf,rand1,rand2,numbers,r)
!$omp do

				do k=1,n_h
				  if(j_rand==1) then                 ! if j_rand=1 produce k-dependent standard seeds for each of the threads and cycles to test consistence of OMP_on and OMP_off results
            idum = j_rand+k
            numbers(1) = iran0(idum)            !a dry call to initialise ran0
            do i=1,seed_size
              numbers(i) = iran0(idum)
            enddo
            call random_seed(put=numbers)
            do i=1,5
              call random_number(r(i))
            enddo
            if(j_verb==1) then
              write(*,*) 'k =',k
              write(*,*) 'Random_seed:',numbers
              write(*,*) 'Reference 1st 5 random numbers:',(r(i),i=1,5)
            endif
				  endif
          
					integration_loop_2: do jint = 1,n_mc/1000
						call random_number(rand1)			!the first cell
						i_r1 = nsuper*rand1+1	!i_r1 may become 1 to nsuper

						do j = 1,n_atom
							at_pos(j,:) = at_pos_file(1:3,i_r1,j)
						enddo

						at_ind(3) = (i_r1-1)/nlayer+1
						at_ind(2) = mod((i_r1-1),nlayer)/n_row(1)+1
						at_ind(1) = mod(i_r1-1,n_row(1))+1
					
					do m=1,1000
            do i=1,3
              call random_number(rand2)		!find the second cell within a cube of 2*n_pdf_range around the first one
              d_ind(i) = anint(2.*n_pdf_range(i)*(rand2-.5))
              if(n_pdf_range(i)==1) d_ind(i) = 0
            enddo

C *** skip selfcorrelations
            if(d_ind(1)==0.and.d_ind(2)==0.and.d_ind(3)==0) then
C					n_outbox = n_outbox+25
C					n_out = n_out+13
              cycle
            endif

            at_ind2 = at_ind+d_ind
          
C *** treat atoms out of the supercell by applying cyclic boundary conditions
            j_base = 0
            do ii =1,3
              if(at_ind2(ii).lt.1) then
                at_ind2(ii) = at_ind2(ii)+n_row(ii)
                j_base(ii) = -n_row(ii)
              endif

              if(at_ind2(ii).gt.n_row(ii)) then
                at_ind2(ii) = at_ind2(ii)-n_row(ii)
                j_base(ii) = n_row(ii)
              endif
            enddo

C *** retrieve the atom positions
            ind2 = nlayer*(at_ind2(3)-1)+n_row(1)*(at_ind2(2)-1)+at_ind2(1)
            if(ind2<1.or.ind2>nsuper) then
              write(*,*) 'out of range: at_ind,at_ind2,ind2,nsuper',at_ind,at_ind2,ind2,nsuper
            endif
          
            do j = 1,n_atom
              at_pos2(j,:) = at_pos_file(1:3,ind2,j)
              if(at_pos2(j,1)/=.0.and.at_pos2(j,2)/=.0.and.at_pos2(j,3)/=0.) at_pos2(j,:) = at_pos2(j,:)+j_base	
            enddo

C *** accumulate all the PDF at once
            do ii=1,n_atom
              if(sum(abs(at_pos(ii,:)))==.0) cycle     !handle unoccupied positions
              do jj=1,n_atom
                if(sum(abs(at_pos2(jj,:)))==.0) cycle
                diff_pos = at_pos(ii,:)-at_pos2(jj,:) 
                i_rdf = anint(norm2(diff_pos*a_par)/pdf_step)
              
                if(i_rdf.ge.1.and.i_rdf.le.n_pdf) then
                  rdf_p2_n(ii,jj,i_rdf,k) = rdf_p2_n(ii,jj,i_rdf,k)+1								
C 							n_inbox = n_inbox+1
                else
C							n_outbox = n_outbox+1
                endif
              enddo	!jj
            enddo		!ii
					enddo !m

					enddo integration_loop_2
				enddo		!k=1,n_h
!$omp end do
!$omp end parallel

C *** sum up contributions from individual cycles and estimate statistical error at the a_par(1) position
        do k=1,n_h                     
          rdf_p2 = rdf_p2 + rdf_p2_n(:,:,:,k)
        enddo
        rdf_p2_n = .0        
      enddo file_loop

 			CALL SYSTEM_CLOCK (COUNT = sc_c2)
 			write(*,*) 
 			write(*,*) 'Accumulation finished         ',(sc_c2-sc_c1)/sc_r,' sec'
 			write(*,*) 

C *** estimate statistical error (Poisson statistics, sqrt(n)**1) at the a_par(1) position
      i_ref = nint(a_par(1)/pdf_step)
      do j=1,n_atom
        rdf_err(j) = (1./sqrt(rdf_p2(j,j,i_ref)))
      enddo
C *** scale by 1/i**2 and normalise the accumulated PDF

      do i=1,n_pdf
        rdf_p2(:,:,i) = rdf_p2(:,:,i)/i**2     ! radial scaling
        rdf_tot(i) = sum(rdf_p2(:,:,i))
        rdf_p2(:,:,i) = rdf_p2(:,:,i)+transpose(rdf_p2(:,:,i))    ! transform to LD/UD form w/o putting zeros beneath/above diagonal, the SUM must
        do ii=1,n_atom
          rdf_p2(ii,ii,i) = rdf_p2(ii,ii,i)*.5
        enddo
      enddo

      rdf_norm = 4.*pi*pdf_step**3*n_int*nfile*(sum(at_occup_r))**2 /(n_norm*a_par(1)*a_par(2)*a_par(3)) ! PDF norm from "first principles" 
                                                                          ! n_int*nfile/n_norm = density of MC events per Å^3
                                                                          ! (sum(at_occup_r))**2/a_par(1)*a_par(2)*a_par(3) = ro_0 atomic/scatt density   
      rdf_p2 = rdf_p2/rdf_norm
      rdf_tot = rdf_tot/rdf_norm
      
C     write(*,*) 'rdf_tot(195:205)',rdf_tot(195:205)
      
C *** Accumulated! now plot PDFs with different atom at_weights

			at_mask = 1
      if(j_weight==1) then
        at_weight(1:n_atom) = 1.
      elseif(j_weight==2) then
        at_weight(1:n_atom) = b_coh(1:n_atom)
      endif

C *** set initial plot range and partial PDFs

      do i=1,n_pdf
				r(i) = i*pdf_step
			enddo

      r_start = pdf_step
      r_end = 10.
      i_start = (r_start-r(1))/pdf_step + 1
      i_end = (r_end-r(1))/pdf_step+1
      n_plot = i_end-i_start+1
        
      if(n_part==0) ind_part = 0
      part_scale = 1.
      
C     write(*,*) 'n_part,ind_part',n_part,ind_part
      
			at_weights_loop: do		
			
C *** extend the RDF matrix by pseudo_atom rows&columns 
       n_atom_tot = n_atom+n_pseudo             
       allocate(rdf_p2_plot(n_atom_tot,n_atom_tot,n_pdf),rdf_tot_plot(n_pdf),at_name_plot(n_atom_tot)) 

       rdf_p2_plot = .0
       rdf_tot_plot = .0

       at_name_plot(1:n_atom) = at_name_par       
       if(n_pseudo>0) at_name_plot(n_atom+1:n_atom+n_pseudo) = at_name_pseudo(1:n_pseudo)
				
								
				write(9,*)
				write(*,'(1x,"Atoms:         ",50(1x,a8))')  (at_name_plot(i),i=1,n_atom_tot)
				write(9,'(1x,"Atoms:         ",50(1x,a8))')  (at_name_plot(i),i=1,n_atom_tot)
				write(*,'(1x,"Atoms no.:  ",50(1x,i8))') (i,i=1,n_atom_tot)
				write(9,'(1x,"Atoms no.:  ",50(1x,i8))') (i,i=1,n_atom_tot)
				write(*,113) trim(at_weight_scheme(j_weight)),(at_weight(i),i=1,n_atom),(sum(ind_pseudo(:,j)*at_weight),j=1,n_pseudo)
				write(9,113) trim(at_weight_scheme(j_weight)),(at_weight(i),i=1,n_atom),(sum(ind_pseudo(:,j)*at_weight),j=1,n_pseudo)
113     format(1x,a,' weights: ',50(1x,f8.4))	
        write(*,'(1x,"Rel_std_error:   ",50(1x,f8.4))') (rdf_err(i),i=1,n_atom)							
				write(*,*) 
				write(*,*) 'Actual masks:'
				write(*,'((50i3))') (at_mask(i),i=1,n_atom)        
				write(*,*) 
        
        if(minval(at_mask(1:n_atom))==1) then  !if all masks =1 don't put them into plot title
          masks = ''
        else
          write(masks,'("  Masks: ",50i1.1)') (at_mask(i),i=1,n_atom)			!to be used in plot titles
        endif
         
C			at_weight = at_weight*at_mask
C			at_weight_sum = sum(at_weight*at_occup_r)/sum(at_occup_r*at_mask)
				at_weight_sum = sum(at_weight*at_mask*at_occup_r)/sum(at_occup_r)
				
        do jj=1,n_atom
          do ii=1,n_atom
						at_weight_matrix(ii,jj) = at_weight(ii)*at_mask(ii)*at_weight(jj)*at_mask(jj)
				  enddo
				enddo
				at_weight_matrix = at_weight_matrix/at_weight_sum**2
				
C			write(*,*) 'at_weight_matrix'
C			do ii=1,n_atom
C					write(*,*) ii,at_weight_matrix(ii,:) 
C			enddo
				       
C *** fill the normal atom part of the RDF matrix  
        do i=1,n_pdf                                                                
          rdf_p2_plot(1:n_atom,1:n_atom,i) = rdf_p2(:,:,i)*at_weight_matrix         !also puts zeros below diagonal
        enddo

        do j=1,n_atom
          do i=1,j  
             rdf_tot_plot(:) = rdf_tot_plot(:)+rdf_p2_plot(i,j,:)
          enddo
        enddo

C *** fill the pseudo_atom part of the RDF matrix, rdf_tot is not concerned          
        if(n_pseudo>0) then
          do j=1,n_pseudo                               ! atom - pseudo_atom corr
            do i=1,n_atom
              do k=1,n_atom
                rdf_p2_plot(i,n_atom+j,:) = rdf_p2_plot(i,n_atom+j,:)+ind_pseudo(k,j)*rdf_p2_plot(i,k,:)
              enddo
            enddo
          enddo
          
          do j=1,n_pseudo                              ! pseudo_atom - pseudo_atom corr
            do i=1,j
              do k=1,n_atom
                rdf_p2_plot(n_atom+i,n_atom+j,:) = rdf_p2_plot(n_atom+i,n_atom+j,:)+ind_pseudo(k,i)*rdf_p2_plot(k,n_atom+i,:)
              enddo
            enddo
          enddo
        endif 

C *** open the PGPLOT graphics window (X11)
      j_xserv = PGOPEN('/xserv')   
      if (j_xserv.LE.0) then    
      	write(*,*) 'Could not open PGPLOT /xserv'
      	STOP
      endif
      
      CALL PGPAP(10.0,.6)     ! define the plot area as landscape
	  	call PGSUBP(1,1)				! PGPLOT window of 1x1 no of panes
      CALL PGASK(.FALSE.)     ! would not ask for <RET>
	  	CALL PGSCRN(0, 'white', IER)	!plot on white background
	  	CALL PGSCRN(1, 'black', IER)
	  	
C *** Set my colors for line plots								  
          CALL PGSHLS (20,.0,.3,.0)     !dark grey
          CALL PGSHLS (21,.0,.4,.7)     !my blue
          CALL PGSHLS (22,120.,.5,1.)   !my red
          CALL PGSHLS (23,240.,.35,.8)  !my green
          CALL PGSHLS (24,60.,.4,.9)    !my violet
          CALL PGSHLS (25,170.,.5,.9)   !my yellow
          CALL PGSHLS (26,320.,.4,.9)   !my turquoise
          CALL PGSHLS (27,.0,.7,.0)     !light grey


        c_min = .0
        c_max = .0

        plot_loop: do
 
C *** Prepare the data to be plotted and smooth the PDF
        if(c_max==.0) then          !only do this at the real start of plot_loop
          if(j_weight==1) then
            c_min = .0
          else
            c_min = minval(rdf_tot_plot(i_start:i_end))
            c_min = min(c_min,-.1)
          endif
          c_max = maxval(rdf_tot_plot(i_start:i_end))
          c_max = .1*(int(10*c_max)+1)
        endif
      
        write(*,*) 'c_min,c_max',c_min,c_max
        
        write(plot_header,'("PDF plot ",a,"    ",a," weights  ")') trim(file_dat_t0),trim(at_weight_scheme(j_weight))
        plot_header = trim(plot_header)//trim(masks)

				scale_loop: do
          write(plot_title,114) trim(subst_name)
114       format(a)   
					CALL PGSLCT(j_xserv)
          CALL PGSCI (1)  !white
          CALL PGSCH(1.)
          CALL PGSLW(2)
          CALL PGSLS (1)  !full
          CALL PGENV(r_start,r_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
          CALL PGLAB('r(A)', 'g(r)',trim(plot_header))  !put the axis labels
          CALL PGSCI (1)  !white
          CALL PGSLW(5)			!operates in steps of 5
          x_plot = r_start+.75*(r_end-r_start)
          y_plot = .9*c_max
          CALL PGLINE(n_plot,r(i_start:i_end),rdf_tot_plot(i_start:i_end))  !plots the curve
          CALL PGSTBG(0)																				 !erase graphics under text
          CALL PGTEXT (x_plot,y_plot,plot_title)
          CALL PGSLW(2)

          do j=1,n_part
            if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
           	write(plot_title,115) at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
          	y_plot = (.9-.06*j)*c_max
115       	format(a,a,' x',f4.1)   
            CALL PGSCI (j+20)  !my red-green-blue
            CALL PGLINE(n_plot,r(i_start:i_end),part_scale(j)*rdf_p2_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title)
          	CALL PGSLW(2)			
          enddo
          
          if(j_weight==1) then              !unit weights
            write(*,*) 'Adjust vertical scale (max) (0 EXIT, -1 adjust horizontal scale)'
            c_min = .0
            c2 = c_max
            read(*,*) c2
            if(c2==.0) then
              exit scale_loop
            elseif(c2<.0) then
              write(*,*) 'Confirm/adjust plot range r_1, r_2, max =',n_pdf*pdf_step
              read(*,*) r_start,r_end 
              if(r_start.lt.pdf_step) r_start = pdf_step
              i_start = (r_start-r(1))/pdf_step + 1
              i_end = (r_end-r(1))/pdf_step+1
              n_plot = i_end-i_start+1
            else
              c_max = c2
            endif
          else
            write(*,*) 'Adjust vertical scale (min, max) (0 0 EXIT, -1 -1 adjust horizontal scale)'
            c1 = c_min
            c2 = c_max
            read(*,*) c1,c2
            if(c1==.0.and.c2==.0) then
              exit
            elseif(c1==-1.and.c2==-1) then
              write(*,*) 'Confirm/adjust plot range r_1, r_2, max =',n_pdf*pdf_step
              read(*,*) r_start,r_end 
              if(r_start.lt.pdf_step) r_start = pdf_step
              i_start = (r_start-r(1))/pdf_step + 1
              i_end = (r_end-r(1))/pdf_step+1
              n_plot = i_end-i_start+1
            else
              c_min = c1
              c_max = c2
            endif
          endif
        enddo scale_loop
         
C **** Prepare and plot the same on .PS, look for existing output files in order not overwrite them				
				if(j_ps==1) then

					jfile = 1
					do						!look for existing .ps files to continue numbering
            if(j_name==0.and.t_single)then
              write(file_ps,1041) trim(file_master),trim(c_fil(j_acc)),jfile,trim(pg_out)
1041   		format(a,'_rdf',a,'_',i2.2,'.',a)      
              write(file_res,1042) trim(file_master),trim(c_fil(j_acc)),jfile
1042   		format(a,'_rdf',a,'_',i2.2,'.txt')      
            else
              write(c_jfile,'("_",i2.2)') jfile
              if(nfile_min<=9999) then
                write(c_nfile_min,'(i4.4)') nfile_min
              elseif(nfile_min>=10000) then
                write(c_nfile_min,'(i8)') nfile_min
              endif
              c_nfile_min = '_'//adjustl(c_nfile_min)
    
              if(nfile<=9999) then
                write(c_nfile,'(i4.4)') nfile
              elseif(nfile>=10000) then
                write(c_nfile,'(i8)') nfile
              endif
              c_nfile = '_'//adjustl(c_nfile)
    
              file_res = trim(file_master)//'_pdf'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//'.txt'							
              file_ps  = trim(file_master)//'_pdf'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//'.'//trim(pg_out)
            endif

						inquire(file=file_ps,exist=found_ps)
						inquire(file=file_res,exist=found_txt)
						if(.not.found_txt.and.(.not.found_ps)) exit				
						jfile = jfile+1
						if(jfile==100) then
							write(*,*)'Tidy up .txt/.ps files to restart count from 01 and type [RET]'
							read(*,*)
							jfile = 1
						endif	
					enddo						
							
					ier = PGOPEN(file_ps//'/CPS')
					IF (ier.LE.0) STOP
					CALL PGASK(.FALSE.)     ! would not ask for <RET>
					CALL PGPAP(11.0,.7)     ! define the plot area as landscape
					CALL PGSUBP(1,1)				! PGPLOT window of 1x1 no of panes
					CALL PGSCRN(0, 'white', IER)	!plot on white background
					CALL PGSCRN(1, 'black', IER)

					plot_title = file_ps//'  '//trim(at_weight_scheme(j_weight))//' weights  '//trim(masks)	
 
          CALL PGSCI (1)  !white
          CALL PGSCH(1.)
          CALL PGSLW(2)
          CALL PGSLS (1)  !full
          CALL PGENV(r_start,r_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes w/o grid
          CALL PGLAB('r(A)', 'g(r)',trim(plot_title))  !put the axis labels
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW(5)			!operates in steps of 5
				  CALL PGLINE(n_plot,r(i_start:i_end),rdf_tot_plot(i_start:i_end))  !plots the curve
					y_plot = .9*c_max
					CALL PGSTBG(0)																				 !erase graphics under text
					CALL PGTEXT (x_plot,y_plot,'PDF_tot '//trim(subst_name)) !the x_plot,y_plot are in world (axis units) coordinates
          CALL PGSLW(2)

          do j=1,n_part
            if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
            write(plot_title,115) at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
          	y_plot = (.9-.06*j)*c_max
            CALL PGSCI (j+1)  !red-green-blue
            CALL PGLINE(n_plot,r(i_start:i_end),part_scale(j)*rdf_p2_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title)				!the x_plot,y_plot are in world (axis units) coordinates
          	CALL PGSLW(2)
          enddo    
					CALL PGCLOS
					write(*,*) ' Postscript output written to: ',file_ps	
					write(9,*)
					write(9,*) '  ',trim(int_mode),' integration',n_int*n_int2,' pairs'
					write(9,*) '  ','Masks:',(at_mask(i),i=1,n_atom)				
					if(n_smooth>1) write(9,*) '  ','Smoothing FWHM:',n_smooth_fwhm	
					write(9,*) ' Postscript output written to: ',file_ps

C *** save the PDF results into an ASCII file (each line corresponds to a distance point)
C     look for existing output files in order not overwrite them				
          if(j_txt==1) then
            open (4,file=file_res)

C *** write PDF file header
            write(4,*) 'Substance:',trim(subst_name),'       ',time_stamp
            write(4,*) 'Input files: ',trim(file_dat_t0),' to ',trim(file_dat),' step',nfile_step									   
            write(4,*) 'Supercell size:',nrow																				
				    write(4,*) 'OMP processes  = ', proc_num_in									
            write(4,*) 'Integration ',trim(int_mode),'  ',n_int*n_int2,'cell pairs, j_rand',j_rand
            write(4,*) 'Unit cell parameter:',a_par																				
            write(4,*) 'Atoms/unit_cell:',n_atom
            write(4,*) 	 'Atoms  :',(('    '//at_name_plot(i)),i=1,n_atom)
            write(4,'(1x,"Atom no. ",50i8)') (i,i=1,n_atom)
            write(4,125) (at_occup_r(i),i=1,n_atom)
            write(4,126) (at_weight(i),i=1,n_atom)
            write(4,127) (at_mask(i),i=1,n_atom)
125     format(1x,"Occup's: ",50f8.4)								
126     format(1x,'Weights: ',50f8.4)								
127     format(1x,'Masks	: ',50i8)								
            write(4,*)	

C *** write the PDFs
            if(j_out==1) then
              write(4,107) ((ii,jj,ii=1,jj),jj=1,n_atom)
              write(4,106) ((trim(at_name_plot(ii)),trim(at_name_plot(jj)),ii=1,jj),jj=1,n_atom)
              do i=1,n_pdf
                write(4,105) r(i),rdf_tot_plot(i),((rdf_p2_plot(ii,jj,i),ii=1,jj),jj=1,n_atom)
              enddo
            else
              write(4,107) (ind_part(1,j),ind_part(2,j),j=1,n_part)
              write(4,106) (trim(at_name_plot(ind_part(1,j))),trim(at_name_plot(ind_part(2,j))),j=1,n_part)
              do i=i_start,i_end
                write(4,105) r(i),rdf_tot_plot(i),(rdf_p2_plot(ind_part(1,j),ind_part(2,j),i),j=1,n_part)
              enddo
            endif
105   	format(1x,f7.3,1x,f8.3,32(2x,5f8.3))	
107   	format(1x,"Partials index  ",32(2x,5(2x,2i3)))	
106   	format('    r[Å]  PDF_tot   ',32(2x,5(a,'_',a,3x)))	
            close(4)
            write(*,*) ' Text output written to: ',file_res	  
            write(9,*) ' Text output written to: ',file_res	  

				  endif 																							!j_txt
				endif 																							!j_ps

C ***   all done, now decide what comes next
			
			 way_point: do
					write(*,*) 'Choose a PDF plot option (FILE output is ',trim(ps_out(j_ps+1)),' size ',trim(size_out(j_out+1)),'):'
	        write(*,*) '       1   REPLOT the PDFs '
	        write(*,'("        2   select max ",i2," partial PDFs & replot")') n_part_max
	        write(*,*) '       3   adjust partial PDF scales & replot '
          write(*,*) '       4   modify atom WEIGHTS (',trim(at_weight_scheme(j_weight)),')'
          write(*,*) '       5   edit atom MASKS ',trim(masks(8:))
          write(*,'("        6   create/modify max ",i2," PSEUDO_ATOMS ")') n_pseudo_max
					write(*,*) '       7   toggle FILE output ',trim(ps_out(mod(j_ps+1,2)+1)),' (mind the J_TXT switch in .PAR)'
					write(*,*) '       8   toggle .TXT output SIZE to ',trim(size_out(mod(j_out+1,2)+1))
          write(*,*) '       9   RESTART with updated weights, masks & pseudo_atoms '
          
          write(*,*) '       0   EXIT'  

          read(*,*) jj
          
          select case(jj)
						case(1) 
							cycle plot_loop

						case(2) 
              write(*,*) '("Confirm/modify up to ", i2," pairs of partial/pseudo PDF indices (-1 -1 skip the rest):")',n_part_max
              do j=1,n_part_max
                j1 = ind_part(1,j)
                j2 = ind_part(2,j) 
                write(*,'(i2,": [",i2,",",i2,"]  ")',advance='no') j,j1,j2
                read(*,*) j1,j2
                if(j1==-1.and.j2==-1) exit
                ind_part(1,j) = j1
                ind_part(2,j) = j2
                if (ind_part(1,j)>ind_part(2,j)) ind_part(:,j) = cshift(ind_part(:,j),1)      !the 2nd one should be bigger  
              enddo
            
              j1 = 0
              do j=1,n_part_max         !compact the partials list
                if(ind_part(1,j)/=0.and.ind_part(2,j)/=0) then
                  j1 = j1+1
                  ind_part(:,j1) = ind_part(:,j)
                endif              
              enddo
              n_part = j1
                     
              cycle plot_loop

						case(3) 
              write(*,"('Confirm/modify partial scale factors ',4f4.1)") part_scale(1:n_part)
              read(*,*) part_scale(1:n_part)
              cycle plot_loop

						case(4) 
              do
                write(*,*) 'Atom weights (1= uniform, 2= neutron b_c^2)'	
                read(*,*) j_weight
                if(j_weight==1.or.j_weight==2) exit
              enddo
              if(j_weight==1) then
                at_weight(1:n_atom) = 1.
              elseif(j_weight==2) then
                at_weight(1:n_atom) = b_coh(1:n_atom)
              endif
              cycle way_point
              
						case(5) 
              write(*,103) (at_mask(i),i=1,n_atom)
103           format(' Actual masks:'(50i3))
              write(*,*)'Type in new ones (0/1):'
              do
                read(*,*)(at_mask(i),i=1,n_atom)
                if(any(at_mask(1:n_atom).ge.0).and.any(at_mask(1:n_atom).le.1)) exit
                write(*,*) 'Input out of range, repeat ...'
              enddo
              cycle way_point

						case(6) 
						  pseudo_loop: do
                if(n_pseudo==0) then
                  write(*,*) 'No actual pseudo_atoms '
                else
                  write(*,*) 'Actual pseudo_atoms:'
                  do i=1,n_pseudo
                    write(*,*) trim(at_name_pseudo(i)),i,'  masks:',ind_pseudo(1:n_atom,i)
                  enddo
                endif
                
                write(*,'("Modify (number), create (max+1) or exit (0): ")',advance = 'no')
                read(*,*) j
                
                if(j>n_pseudo_max) then
                  write(*,*) 'Maximum number of pseudo-atoms reached, consider modifying existing ones'
                  write(*,'("Modify (number): ")',advance = 'no')
                  read(*,*) j
                endif
                  
                if(j<=0) then
                  exit pseudo_loop
                elseif(j<=n_pseudo) then
                  write(*,'("New indices: ")',advance = 'no')
                  read(*,*)ind_pseudo(1:n_atom,j)
                  cycle pseudo_loop
                else
                  write(*,'("Type pseudo_atom name (<=4 char): ")',advance = 'no')
                  read(*,*) at_name_pseudo(n_pseudo+1)
                  write(*,'("Type pseudo_atom indices (n_atom): ")',advance = 'no')
                  read(*,*)ind_pseudo(1:n_atom,n_pseudo+1)
                  n_pseudo = n_pseudo+1
                  cycle pseudo_loop
                endif 
              enddo pseudo_loop                
              cycle way_point

						case(7) 
              j_ps = j_ps+1
              j_ps = mod(j_ps,2)
              cycle way_point

						case(8) 
              j_out = j_out+1
              j_out = mod(j_out,2)
              cycle way_point

						case(9) 
              exit plot_loop

						case(0) 
              exit at_weights_loop
					end select

				enddo way_point
			
        enddo plot_loop

        deallocate(rdf_p2_plot,rdf_tot_plot,at_name_plot)  			
      enddo at_weights_loop

			flush(9)
			close(9)
			            
      end program mp_pdf55
      
C **** minimal random number generator (Numerical Recipees))
C           
C     Minimal” random number generator of Park and Miller. 
C     Returns a uniform random deviate between 0.0 and 1.0. 
C     Set or reset idum to any integer value (except the unlikely value MASK) to initialize the sequence; 
C     idum must not be altered between calls for successive deviates in a sequence.
C
      FUNCTION iran0(idum)
      INTEGER(4) idum,IA,IM,IQ,IR,MASK
      PARAMETER (IA=16807,IM=2147483647,IQ=127773,IR=2836,MASK=123459876)
      INTEGER(4) k,iran0
C       write(*,*) 'R: idum0',idum
      idum=ieor(idum,MASK)
C      write(*,*) 'R: idum1',idum
      k=idum/IQ
C     write(*,*) 'R: k',k
      idum=IA*(idum-k*IQ)-IR*k
C     write(*,*) 'R: idum2',idum
      iran0 = idum-IM/2
C     write(*,*) 'R: iran0',iran0
      idum=ieor(idum,MASK)
C     write(*,*) 'R: idum2',idum
C     write(*,*) 
      return
      END      
      
      
  
C **** string conversion to all upper case
C     
 			subroutine up_case (string)

			character(*), intent(inout)	:: string
			integer											:: j, nc
			character(len=26), parameter	:: lower = 'abcdefghijklmnopqrstuvwxyz'
			character(len=26), parameter	:: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

			do j = 1, len(string)
				nc = index(lower, string(j:j))
				if (nc > 0) string(j:j) = upper(nc:nc)
			end do

			end subroutine up_case	     


 
      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set the RAINBOW palette of colors to be used by PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL JKL(9), JKR(9), JKG(9), JKB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA JKL / 0,   0.1,  0.2,  0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA JKR / 1.0, 0.5,  0.1,  0.1,  0.3,  1.0,  1.0, 1.0, 1.0/
      DATA JKG / 1.0, 0.7,  0.4,   .8,  1.,   1.0,  0.6, 0.0, 1.0/
      DATA JKB / 1.0, 1.0,  0.9,   .6,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.6) THEN
C        -- JK rainbow
         CALL PGCTAB(JKL, JKR, JKG, JKB, 9, CONTRA, BRIGHT)
      END IF
      END
      
   
  
  
     