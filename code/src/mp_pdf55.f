  
      program mp_pdf55

C *************************************************************************************
C *****
C ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   calculates the pair distribution functions (PDF) for simulated supercell data
C *****
C**********					Copyright (C) 2023  Jiri Kulda, Grenoble/Prague          **********
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
C ***** Ver. 1.1 The correlation is asymmetric: n_int random cell against the complete supercell
C ***** Ver. 1.2 the whole supercell is read into memory once for ever
C ***** Ver. 1.21 the correlation is symmetric n_int random cells against n_int random cells
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
C ***** Ver. 1.42	- lighter & parallelised version (no correlation histograms, no IDOM)
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
			use singleton

C *** notes on FT & SINGLETON:
C                 - check whether automatic scaling by 1/sqrt(N) inside SINGLETON is blocked or not - now it is!
C                 - if YES  you need to do it explicitly within the code	
C                 - the scaled FT transform restitutes properly f(x) after going forth and back
C                 - NEVERTHELESS the FT value will depend on the sampling frequency as SQRT(N)
C                 - if the FT value itself is important and should be invariant w.r.t. the sampling frequency the sample sequence has to be scaled by 1/SQRT(N) as well
  		
CC  		use mp_nopgplot					! uncomment when not able to use PGPLOT, compile and link together with the module mp_nopgplot.f

      integer,parameter :: l_rec  =  1024		    !record length in real(4)
      integer,parameter :: n_mc  =  1e6		    !record length in real(4)
      real,parameter    :: pi = 3.14159
      character(4),allocatable  ::	at_name_par(:),at_label(:),at_name_ext(:),at_name_pseudo(:),at_name_plot(:),pdf_out(:)
      character(10),allocatable ::	corr_name(:),curve_label(:),x_label(:),y_label(:)

      integer,allocatable ::  at_mask(:),ind_part(:,:),ind_pseudo(:,:),numbers(:),ind_hist(:,:,:),ind_pdf(:,:,:),ind_at1(:) 

			real,allocatable ::  w(:),rdf_p2(:,:,:),rdf_p2_n(:,:,:,:),rdf_p2_ws(:,:,:),rdf_p2_plot(:,:,:),rdf_err(:),rdf_fft(:)
      real,allocatable ::  at_base(:,:),b_coh(:),at_weight(:),x_ffpar(:,:),x_ffq(:,:),at_weight_matrix(:,:),f_smooth(:)
      real,allocatable,target ::  r(:),q(:)
      real,pointer ::      x(:)
	  
			integer ::       i_rdf,n_pdf,j_acc,n_int,n_pseudo,i_ref,j_pdf,j_rand,n_ind(3),n_skip1,n_skip2
      real    ::       rdf_dist,pdf_range,rdf_sum,rdf_tot_sum,rdf_norm,at_weight_av,at_pos(3),at_pos2(3)
      real 		:: 			 arg,c_min,c_max,c_max2,c1,c2,x_start,x_end,pdf_step,q_step,ro_0,c_smooth,rnd(5)
      
      character(4)   :: version,head,atom,ps_out(2),size_out(2)
      character(10)  :: at_weight_scheme(3),pg_out,pg_ext,string,section,c_date,c_time,c_zone,c_nfile_min,c_nfile,c_jfile
      character(16)  :: sim_type_par,data_type,string16,filter_name
      character(40)  :: subst_name,file_master,file_inp,file_out,time_stamp,int_mode,x_file_name
      character(60)  :: file_dat,file_dat_t0,file_res,file_ps,file_log,line,masks,smooth
      character(128) :: cwd_path,plot_header,plot_title
      character(l_rec):: header_record
      
			logical ::  nml_in,found,found_txt,found_ps,t_single

      integer(8) ::  n_pdf_range(3),n_norm,n_mc_max                !this number may be biiiiig!   
      integer ::  j_proc,proc_num,proc_num_in,thread_num,thread_num_max,m,j_name,jm,j1m,j2m,j_mode,n_mode
      integer ::  i_start,i_end,n_step,nfile_step,n_h,j_head_in,hist_ind(3),j_verb,n_tot,n_head
      integer ::  i,j,k,ii,jj,ind,ind2,jat,i_rec,ind_rec,nrec,nfile,nfile_min,nfile_max,jfile,i_at1
      integer ::  j1,j2,j_plane,j_grid,j_logsc,j_ps,j_out,j_txt,i_seed,i_r1,i_r11,i_r12,i_time(8)
      integer ::  at_ind(3),at_ind2(3),d_ind(3),d_ind_shift(3),j_dom,j_base(3),n_plot,n_smooth,n_smooth_fwhm,n_part,n_part_max,n_atom_tot,n_pseudo_max,n_pdf_grid(3)
      integer ::  j_atom,j_weight,j_xray,j_edit,j_mask,ifile,jint,n_atom,sc_c2,sc_c1,ier,ios,n_pix_step
      
      real :: t_dump,filter_fwhm,tt0,tt,b_sum,at_sum,rand1,rand2,sc_r,at_displ,p_size,rdf_tot_err,pdf_pix,pdf_pix_shift(3),a_par_grid(3)
      integer :: rand1_seed(8),rand2_seed(8),seed_size

      real :: pdf_grid_min(3),pdf_grid_max(3),diff_pos(3),diff_pos_norm,x_plot,y_plot,part_scale(4),d_pos(3)

C **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
      character(4),allocatable :: at_name_out(:)
      integer(4),allocatable   :: at_ind_in(:),nsuper_r(:)	

      real(4),allocatable ::	at_pos_in(:),at_occup_r(:)

      character(16)  :: sim_type,dat_type,input_method,file_par,dat_source
      integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,idum,iran0,mult,j_pgc
      real(4)        :: rec_zero(l_rec),t_ms,t_step,t0,t1,a_par(3),angle(3),a_par_pdf(3),temp

      namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,
     1    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
      namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
     
      namelist /mp_gen/ j_verb,j_proc       
      namelist /mp_out/ j_weight,j_logsc,j_ps,j_txt,p_size,j_grid,pg_out,j_out,j_pgc       
      										!general rule: namelists of tools should only contain their local parameters
                          !what is of global interest they should pass into data_header
			namelist /mp_pdf/ pdf_range,pdf_step,x_end,a_par_pdf,pdf_pix,pdf_pix_shift,j_rand,n_h,j_acc,j_pdf,n_part_max,n_pseudo_max,n_cond    !j_acc not used
			
C **** PGPLOT stuff
      INTEGER :: PGOPEN,j_xserv

      j_xserv = 0
      
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

			write(*,*) '*** Program MP_PDF 1.55 ** Copyright (C) Jiri Kulda (2019-2023) ***'
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
      n_mode = 3
      allocate(pdf_out(n_mode),curve_label(n_mode),x_label(n_mode),y_label(n_mode))
      curve_label = (/'g_tot','G_tot','S_tot'/) 
			pdf_out = (/'g(r)','G(r)','S(Q)'/)			
      y_label = pdf_out
      x_label = (/'r[A]  ','r[A]  ','Q[A-1]'/)
			ps_out = (/'OFF ','ON  '/)			!PGPLOT
			size_out = (/'S   ','XXL '/)		!PGPLOT
			
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
      
      allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom))
			allocate(at_label(n_atom),at_name_par(n_atom),at_name_ext(n_atom))
			allocate(at_base(n_atom,3),at_weight(n_atom),at_weight_matrix(n_atom,n_atom),at_mask(n_atom),rdf_err(n_atom))
			allocate(b_coh(n_atom),x_ffpar(n_atom,9),SOURCE=.0)												!we need this to read .par
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

      j_pdf = 0
      j_txt = 0               !defaults for &mp_out
      j_out = 0
      j_weight = 1
      j_logsc = 0
      j_grid = 0
      j_ps = 0 
      pg_out = 'ps'
      p_size = 7.
      rewind(4)
      read(4,nml=mp_out)
      
      call down_case(pg_out) 
      if(index(pg_out,'png')/=0) then
        pg_ext = '.png'
      else
        pg_ext = '.ps'
      endif
      
      if(j_weight==3) then
        write(*,*) 'Xray weights not yet implemented, setting to Uniform'
        j_weight = 1
      endif

      at_weight_scheme(1) = 'Uniform'
			at_weight_scheme(2) = 'Neutron'
			at_weight_scheme(3) = 'Xray'

      n_h = 0
      pdf_step = 0.02
      pdf_pix = 1.
      a_par_pdf = 1.
      n_part_max = 4         !defaults for &mp_pdf
      n_pseudo_max = 0
      j_rand = 1
      x_end = 10.
      rewind(4)
      read(4,nml=mp_pdf)
      
      allocate(ind_pseudo(n_atom,n_atom+n_pseudo_max+1),at_name_pseudo(n_pseudo_max+1),ind_part(2,n_part_max))       !1st pseudo is TOT by default
      ind_pseudo = 0
      ind_part = 0
      at_name_pseudo = ''
      at_name_pseudo(1) = 'TOT'
      ind_pseudo = 0
      do j=1,n_atom
        ind_pseudo(j,j) = 1
      enddo
      ind_pseudo(:,n_atom+1) = 1

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
            n_pseudo = 1
            found = .false.
            exit
          endif
          found = (string(1:6)==section)
          if(found) exit	                  !found the 'pseudo_atom' part of the .par file
        enddo
        
        if(found) then
          do j=2,n_pseudo_max
            read(4,*,iostat=ios) at_name_pseudo(j),ind_pseudo(1:n_atom,n_atom+j)	              ! pseudo_atom name and indices
            if(ios/=0) then
              n_pseudo = j-1
              exit
            endif
            n_pseudo = j
          enddo
        endif
      endif
			close(4)
			
C     do i=1,n_pseudo
C       write(*,*) trim(at_name_pseudo(i)),i,'  masks:',ind_pseudo(:,n_atom+i)
C     enddo

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
      bc_loop: do j=1,n_atom
        rewind(4)
        do i=1,210
          read(4,*) atom
          if(atom==trim(at_label(j))) then
            backspace(4)
            read(4,*) atom,b_coh(j)
            cycle bc_loop
          endif
        enddo
        write(*,*) 'b_coh for ',trim(at_label(j)),' not found,'
        write(*,*) 'check your spelling and the neutron_xs.txt table; use unit weights'
      enddo bc_loop
      close(4)
      b_coh = .1*b_coh  !convert b_coh from FM to 10^12 cm
            
C *** read Xray formfactor parameters 
      open(4,file='xray_ff.txt',action='read',status ='old',iostat=ios)
      if(ios.ne.0) then
        open(4,file='/usr/local/mp_tools/ref/xray_ff.txt',action='read',status ='old',iostat=ios)
        if(ios.ne.0) then
          do
            write(*,*) 'File xray_ff.txt not found, type valid access path/filename'
            read(*,'(a)') x_file_name
            open(4,file=trim(x_file_name),action='read',status ='old',iostat=ios)
            if(ios==0) exit
            write(*,*) 'File',trim(x_file_name),' not found, try again ...'
          enddo
        endif
      endif
      xff_loop: do j=1,n_atom
        rewind(4)
        do i=1,210
          read(4,*) atom
          if(atom==trim(at_label(j))//trim(at_name_ext(j))) then
            backspace(4)
            read(4,*) atom,x_ffpar(j,1:9)
            cycle xff_loop
          endif
        enddo
        write(*,*) 'Xray formfactor for ',trim(at_label(j))//trim(at_name_ext(j)),' not found,'
        write(*,*) 'check your spelling and the neutron_xs.txt table; use unit weights'
      enddo xff_loop
      close(4)
 

C *** write overview of atom data
			write(*,*)
      write(*,*) 'Substance name: ',subst_name	  
			write(*,*) 'Atoms from ',trim(file_inp)
      do j=1,n_atom
        write(*,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j)
        write(9,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j)
      enddo
      b_sum = sum(at_occup_r(1:n_atom)*b_coh(1:n_atom))
      at_sum = sum(at_occup_r(1:n_atom))

			write(*,*) 

C **** use the .PAR value of a_par in case of input in LATTICE units      
      if(a_par(1)*a_par(2)*a_par(3)==1.and.a_par_pdf(1)*a_par_pdf(2)*a_par_pdf(3)/=0.) then
        a_par = a_par_pdf
        write(*,*) 'Setting a_par to',a_par
      endif
 
 
C **** check the PDF range    
      if(n_cond>0) then                         !box with periodic bound_cond
        if(pdf_range>minval(a_par*n_row)) then
          write(*,*) 'PDF range',pdf_range,' exceeds box size',minval(a_par*n_row)
          write(*,*) 'Setting PDF range to',minval(a_par*n_row)
          pdf_range = minval(a_par*n_row)        
        endif
      else                                      !box with non-periodic bound_cond
        write(*,*) 'ATTENTION: non-periodic boundary conditions!'
        if(pdf_range> .4*minval(a_par*n_row)) then
          write(*,*) 'Box with non-periodic boundary conditions:'
          write(*,*) 'PDF range',pdf_range,' exceeds 40% box size',.4*minval(a_par*n_row)
          pdf_range = int(.4*minval(a_par*n_row))        
          write(*,*) 'Setting PDF range to',pdf_range
        endif      
      endif
      
      n_pdf = anint(pdf_range/pdf_step)
      if(n_cond==0) then
        n_mc_max =  real(n_tot)**2/(500.*n_mc)       !more exactly: (.008*n_tot)*n_atom*(4./3.)*pi*product(.4*n_row)=.00215*n_tot**2
      else
        n_mc_max = real(n_tot)**2/(2.*n_mc)    !(1+n_tot/n_mc)*n_atom*product(n_row)/2 product(n_row)/2 is roughly the volume of sphere with radius of pdf_range_max !n_mc=1e6
      endif
 
 
C **** establish the overlay PDF_GRID 
      n_pdf_grid = anint(n_row/pdf_pix)     
C     d_pos = pdf_pix/a_par               !d_pos is indexing mesh cell size in lattice units (default corresponding to 1Å)
      d_pos = pdf_pix                     !d_pos is indexing mesh cell size in lattice units, has to be commensurate with the box size, should be commensurate with the cell size
      a_par_grid = a_par*pdf_pix           !a_par_grid is its lattice parameter in Å
!      pdf_pix_shift                      !pdf_pix_shift specified in .PAR is offset of pdf_grid to lattice to get most atoms in centres of its cells

      if(j_verb==1) write(*,*) 'Fine grid size, step and offset:',n_pdf_grid,d_pos,pdf_pix_shift
      if(j_verb==1) write(*,*) 

      n_pix_step = anint(sqrt(3.)*pdf_pix/pdf_step)       !we shall avoid correlations within the 1st pixel range
     
      
C *********************  OpenMP initialization start  *******************************      
C
			proc_num_in = j_proc
			thread_num_max = omp_get_max_threads( )			!this gives maximum number of threads available (limited by concurrent tasks??)
			proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
			if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
			call omp_set_num_threads(proc_num_in)
			thread_num = omp_get_num_threads( )				  !this gives threads available at this moment: here always = 1 as we are outside of a PARALLEL range


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

C *** initialise the random_number generator
      call random_seed(size=seed_size)               !Puts size of seed into seed_size
      allocate(numbers(seed_size))

      if(j_rand==0) then                     
        write(*,*) 'Random_number: j_rand =',j_rand,'  the system will supply unique, machine dependent seeds each time this code runs'
        call random_seed(get=numbers)               !Gets actual seeds 
        if(j_verb==1) then
          write(*,*) 'Random_number seed size:',seed_size
          write(*,*) 'Random_seed:',numbers
          write(*,*) 'Reference 1st 5 random numbers:',(rnd(i),i=1,5)
        endif
      elseif(j_rand==1) then                     !if j_rand>1 generate a seed for later reference & numerical reproducibility checks
        write(*,*) 'Random_number: j_rand =',j_rand,'  the system will supply k-dependent standard seeds for each of the OMP threads (use only for testing the consistence of OMP_on/OMP_off results)'
      elseif(j_rand>1) then                     !if j_rand>1 generate a seed for later reference & numerical reproducibility checks
        write(*,*) 'Random_number: j_rand =',j_rand,'  this seeding reference can be used to exactly reproduce this MC-run later on'
       idum = j_rand
        numbers(1) = iran0(idum)            !a dry call to initialise ran0
        do i=1,seed_size
          numbers(i) = iran0(idum)          !use a trivial random number generator to produce the seeds (they could even be all the same small ones, but ...)
        enddo
        call random_seed(put=numbers)       !Produce a seed to start, this permits to reproduce exactly the same results on the same system 
        do i=1,5
          call random_number(rnd(i))
        enddo
        if(j_verb==1) then
          write(*,*) 'Random_seed:',numbers
          write(*,*) 'Reference 1st 5 random numbers:',(rnd(i),i=1,5)
        endif
      endif
      write(*,*)

C *** initialise the MC integration
			if(n_h==0) then
        write(*,*) 'MC sampling pairs per frame ([x 10^6], ',trim(adjustl(string)),' max):'
        read(*,*)   n_h
      endif
      if(n_h.gt.n_mc_max) then                   !n_mc_max is in units of n_mc to avoid overflow for large boxes
        write(*,*) 'WARNING: n_h exceeds max number of 10^6 atom pairs ',n_mc_max
        write(*,*) 'type in a reduced n_h (<',n_mc_max,'):'
        read(*,*)   n_h
      endif
			n_int = n_h*n_mc                           !number of MC cycles per snapshot
			int_mode = 'Monte Carlo'
			write(*,*) trim(int_mode),' integration over',n_h,'*10^6 cell pairs'

C *** Allocate and clear the histogram arrays for accumulation across several snapshots	       
			allocate (at_pos_in(4*n_tot),at_ind_in(4*n_tot),ind_at1(n_tot))
			allocate(ind_pdf(n_pdf_grid(1),n_pdf_grid(2),n_pdf_grid(3)))
			allocate(rdf_p2(n_atom,n_atom,n_pdf),rdf_p2_n(n_atom,n_atom,n_pdf,n_h))

      ind_pdf = 0
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
      
      nfile = 0
 
      file_loop: do ifile=nfile_min,nfile_max,nfile_step									

        n_skip1 = 0

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

C     CALL SYSTEM_CLOCK (COUNT = sc_c2)
        write(*,*)
        write(*,*)'Input: ',trim(file_dat)
        write(9,*)'Input: ',trim(file_dat)

        open(1,file=file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
        if(ios.ne.0) then
          write(*,*) 'File ',trim(file_dat),' not opened! IOS =',ios
          write(*,*) 'Skip(1), stop execution(0)?'
          read(*,*) jj
          if(jj==1) exit file_loop
          if(jj==0) stop
        endif
        nfile = nfile+1

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
        close(1)
C         write(*,*) 'Read finished'


C **** searching for atom_1 candidates in a non-periodic box
        if(n_cond==0) then
          ind_at1 = 0
          i_at1 = 0
          do i=1,n_tot
            if(abs(at_pos_in(4*(i-1)+1))<.1*n_row(1).and.abs(at_pos_in(4*(i-1)+2))<.1*n_row(2).and.abs(at_pos_in(4*(i-1)+3))<.1*n_row(3)) then
              i_at1 = i_at1+1
              ind_at1(i_at1) = i
            endif
          enddo
          write(*,*) 'Atom_1 pool:',i_at1
          if(n_h>nint((i_at1*1000)/1.e6)) then
            write(*,*) 'WARNING: n_h exceeds max number of 10^6 atom pairs ',nint((i_at1*1000)/1.e6)
            write(*,*) 'restart with an adequate n_h and use a larger number of frames'
            stop
          endif
        endif

C *** reindexing atoms on a finer atomary grid	
        ii =1		
        do i=1,n_tot
          if(maxval(abs(at_pos_in(4*(i-1)+1:4*(i-1)+3)))>.0) then
            d_ind = 1+anint((at_pos_in(4*(i-1)+1:4*(i-1)+3)-pdf_pix_shift+n_row/2)/d_pos)   !+d_ind_shift
C           if(j_verb==1.and.i==ii)then
C             write(*,*)'i,at_pos_in',i,at_ind_in(4*(i-1)+1:4*(i-1)+3),at_pos_in(4*(i-1)+1:4*(i-1)+3)
C             write(*,*)'d_ind convert',d_ind
C           endif

            do j=1,3                                                  !handle basis atoms displaced out of the nominal box
              if(d_ind(j)<=0) then
                d_ind(j) = d_ind(j)+n_pdf_grid(j)
                at_pos_in(4*(i-1)+j) = at_pos_in(4*(i-1)+j)+n_row(j)
              endif
              if(d_ind(j)>n_pdf_grid(j)) then
                d_ind(j) = d_ind(j)-n_pdf_grid(j)
                at_pos_in(4*(i-1)+j) = at_pos_in(4*(i-1)+j)-n_row(j)
              endif
            enddo
C           if(j_verb==1.and.i==ii)then
C             write(*,*)'d_ind shift',d_ind,'next ii?'
C             read(*,*)ii
C           endif
           
           if(ind_pdf(d_ind(1),d_ind(2),d_ind(3))/= 0) then
              write(*,*) 'Fine grid cell already taken (you can accept a few by RETURN, else modify pdf_pix in .PAR):'
              write(*,*) 'at_ind_in,at_pos_in',at_ind_in(4*(i-1)+1:4*(i-1)+4),'  ',at_pos_in(4*(i-1)+1:4*(i-1)+3)
C             write(*,*) 'i,d_ind,ind_pdf(d_ind(1),d_ind(2),d_ind(3))',i,d_ind,ind_pdf(d_ind(1),d_ind(2),d_ind(3))
C             write(*,*) 'OLD: at_ind_in,at_pos_in',at_ind_in(4*(ind_pdf(d_ind(1),d_ind(2),d_ind(3))-1)+1:4*(ind_pdf(d_ind(1),d_ind(2),d_ind(3))-1)+4),
C    1              '  ',at_pos_in(4*(ind_pdf(d_ind(1),d_ind(2),d_ind(3))-1)+1:4*(ind_pdf(d_ind(1),d_ind(2),d_ind(3))-1)+3)
             read(*,*)
             n_skip1 = n_skip1+1
            endif
            ind_pdf(d_ind(1),d_ind(2),d_ind(3)) = i
            at_ind_in(4*(i-1)+1:4*(i-1)+3) = d_ind            !re-use at_ind_in(:,1:3) to store the new indices		
          endif	
        enddo

        CALL SYSTEM_CLOCK (COUNT = sc_c2)
        write(*,*) 
        if(j_verb==1) write(*,*) 'Input finished         ',(sc_c2-sc_c1)/sc_r,' sec'
        
        if(n_skip1>10) then
          write(*,*) 'Number of atoms skipped due to fine grid double occupancy:', n_skip1
          write(*,*) 'Modify .PAR pdf_pix (0), continue (1)?'
          read(*,*) jj
          if(jj==0) stop
        endif
 			
C *************  the integration loop: cycle over the site pairs to accumulate the PDFs ******************

        write(*,*) 'Accumulating the PDFs ...'     
        n_skip1 = 0
        n_skip2 = 0

C *** MonteCarlo integration

CCCCCCC OMP continuation line must have a valid character in column 6  CCCCCCC
!$omp parallel shared(at_pos_in,at_ind_in,ind_pdf,rdf_p2_n,a_par,pdf_step,n_tot,n_pdf_grid,n_row,n_atom,n_h,j_rand,n_skip1,n_skip2,i_at1)
!$omp& private(at_pos,at_pos2,diff_pos,diff_pos_norm,at_ind,at_ind2,d_ind,j_base,i_r1,i_rdf,idum,rand1,rand2,numbers,r,jint,m,ii,jj)
!$omp do

        do k=1,n_h                           ! runs over 1000*1000 atom pairs for each k
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
        
          jint = 0
          integration_loop_2: do              !runs over 1000 "first atoms"
            call random_number(rand1)			!the first cell
            if(n_cond>0) then
              i_r1 = n_tot*rand1+1	!i_r1 may become 1 to n_tot
              at_pos(:) = at_pos_in(4*(i_r1-1)+1:4*(i_r1-1)+3)       !randomly select first atom
              if(sum(abs(at_pos(:)))==.0) cycle
              at_ind = at_ind_in(4*(i_r1-1)+1:4*(i_r1-1)+3)   !at_ind_in(:,1:3) has been re-used to store the fine-grid indices
              ii = at_ind_in(4*(i_r1-1)+4)
            else
              i_r1 = i_at1*rand1+1	!i_r1 may become 1 to i_at1
              if(ind_at1(i_r1)<=0) write(*,*) i_at1,rand1,i_r1,'?'
              if(ind_at1(i_r1)<=0) read(*,*)
              at_pos(:) = at_pos_in(4*(ind_at1(i_r1)-1)+1:4*(ind_at1(i_r1)-1)+3)       !randomly select first atom among their shortlist
              if(sum(abs(at_pos(:)))==.0) cycle
              at_ind = at_ind_in(4*(ind_at1(i_r1)-1)+1:4*(ind_at1(i_r1)-1)+3)   !at_ind_in(:,1:3) has been re-used to store the fine-grid indices
              ii = at_ind_in(4*(ind_at1(i_r1)-1)+4)
            endif 
          

C           write(*,*) 'jint,at_ind,at_pos',jint,at_ind,at_pos
            m = 0
            m_cycle: do                 !runs over 1000 "second atoms"
              do
                do i=1,3
                  call random_number(rand2)		!generate the position within a cube [-1,1]; only points within the unit sphere will be retained
                  diff_pos(i) = 2.*rand2-1.
                enddo
                diff_pos_norm = norm2(diff_pos)
                if(diff_pos_norm>.0.and.diff_pos_norm<=1.) exit     !only diff_pos within unit sphere is kept !throwing away also the unlikely event of diff_pos=.0, not d_ind=0!
              enddo

              diff_pos = diff_pos/diff_pos_norm             !now we project random points onto spherical surface by d_ind/norm2(d_ind) ....
              call random_number(rand2)
              diff_pos = pdf_range*rand2*diff_pos     ! .... and select a random point on their radius (rand2) to emulate scaling 1/i**2, diff_pos is in Å now

              d_ind = anint(diff_pos/a_par_grid)       ! .... and select a random point on their radius (rand2) to emulate scaling 1/i**2   (we needn't switch to l.u. and back)

              at_ind2 = at_ind+d_ind          !all of them are on the fine grid

C *** treat atoms out of the supercell by applying cyclic boundary conditions
              j_base = 0
              do j =1,3                                !n_pdf_grid is the supercell size on the fine grid
                if(at_ind2(j).lt.1) then
                  at_ind2(j) = at_ind2(j)+n_pdf_grid(j)
                  j_base(j) = -n_row(j)
                endif

                if(at_ind2(j).gt.n_pdf_grid(j)) then
                  at_ind2(j) = at_ind2(j)-n_pdf_grid(j)
                  j_base(j) = n_row(j)
                endif
              enddo

C *** retrieve the 2nd atom positions
              i = ind_pdf(at_ind2(1),at_ind2(2),at_ind2(3))
              if(i==0) then
                n_skip1 = n_skip1+1             
C               if(n_skip1 == 100*(n_skip1/100)) write(*,*) 'n_skip1=',n_skip1,at_ind2
                cycle                                !cycle over void indexing mesh cells
              endif
              at_pos2 = at_pos_in(4*(i-1)+1:4*(i-1)+3)             
              
              if(at_pos2(1)/=.0.and.at_pos2(2)/=.0.and.at_pos2(3)/=0.) at_pos2(:) = at_pos2(:)+j_base	!handle possible PBC effects
              jj = at_ind_in(4*(i-1)+4)

              at_pos2 = at_pos2-at_pos

C *** accumulate all the PDF at once
              i_rdf = anint(norm2(at_pos2*a_par)/pdf_step)+1     !we need 0 in the 1st channel for FT
              
C             write(*,*)'d_at_pos2,i_rdf',at_pos2,i_rdf

              if(i_rdf>n_pix_step.and.i_rdf<=n_pdf) then         !we are avoiding the low-R range of self-correlations
                rdf_p2_n(ii,jj,i_rdf,k) = rdf_p2_n(ii,jj,i_rdf,k)+1.      
                m = m+1
                if(m<=0.or.m>1000) then 
                  write(*,*) 'Trouble: m=',m
                  read(*,*)
                endif
              else
                n_skip2 = n_skip2+1             
C               write(*,*) 'skip2',jint,m,d_ind,at_pos,at_pos2,i_rdf,n_skip2
C               read(*,*)
                cycle
              endif
                if(m==1000) exit
            enddo m_cycle   !m
            jint = jint+1
C           if(jint == 100*(jint/100)) write(*,*) 'jint=',jint
            if(jint==n_mc/1000) exit integration_loop_2
					enddo integration_loop_2
				enddo		!k=1,n_h
!$omp end do
!$omp end parallel

        write(*,*)'MC hits to empty fine grid cells [%]:',(100.*n_skip1)/(real(n_h*n_mc)+real(n_skip1))
        write(*,*)'MC hits out of PDF range [%]:',(100.*n_skip2)/(real(n_h*n_mc))

C *** sum up contributions from individual cycles and estimate statistical error at the a_par(1) position
        do k=1,n_h                     
          rdf_p2 = rdf_p2 + rdf_p2_n(:,:,:,k)
        enddo
        
C       if(j_verb==1) write(*,*) 'n_pdf,n_int,sum(rdf_p2)',n_pdf,n_int,sum(rdf_p2)
        
        rdf_p2_n = .0
        ind_pdf = 0
      enddo file_loop

 			CALL SYSTEM_CLOCK (COUNT = sc_c1)
 			write(*,*) 
 			write(*,*) 'Accumulation finished         ',(sc_c1-sc_c2)/sc_r,' sec'
C       if(j_verb==1) write(*,*) 'n_skip2',n_skip2
 			write(*,*) 

C *** normalise the accumulated PDF
      rdf_norm = nfile*n_int/real(n_pdf-n_pix_step)          !PDF_STEP completes the 4*pi*r**2*dr volume element, PDF_RANGE goes with the MC scaling
      rdf_p2 = rdf_p2/rdf_norm

      if(j_verb==1) then    
        write(*,*) 'PDF_norm',rdf_norm
C       write(*,*) 'rdf_tot mean',sum(rdf_p2(:,:,:)/real(n_pdf-n_pix_step))
CC       write(*,*) 'rdf_tot mean1',sum(rdf_p2(:,:,n_pix_step/2+1:n_pdf/4)/real(n_pdf/4-n_pix_step))
C       write(*,*) 'rdf_tot mean1',sum(rdf_p2(:,:,n_pix_step+1:n_pdf/4)/real(n_pdf/4-n_pix_step))
C       write(*,*) 'rdf_tot mean2',sum(rdf_p2(:,:,n_pdf/4+1:n_pdf/2)/real(n_pdf/4))
C       write(*,*) 'rdf_tot mean3',sum(rdf_p2(:,:,n_pdf/2+1:3*n_pdf/4)/real(n_pdf/4))
C       write(*,*) 'rdf_tot mean4',sum(rdf_p2(:,:,3*n_pdf/4+1:n_pdf)/real(n_pdf/4))
C       do i=1,n_atom
C         write(*,*) 'rdf_tot mean4_diag i =',i,sum(rdf_p2(i,i,3*n_pdf/4+1:n_pdf)/real(n_pdf/4))
C       enddo
        write(*,*) 
      endif

C *** suppress termination effects in last pixels
      do i=1,n_atom
        do j=1,n_atom
          rdf_p2(i,j,n_pdf-n_pix_step:n_pdf) = at_occup_r(i)*at_occup_r(j)/sum(at_occup_r)**2        !the asymptote is 1 for rdf_tot: we are generating directly g(r) with unit weights
        enddo
      enddo

C *** estimate statistical error (Poisson statistics, sqrt(n)^-1) at the a_par(1) position
C     i_ref = nint(a_par(1)/pdf_step)
      rdf_tot_err = 1./sqrt(rdf_norm)                      !rdf_norm is the average count rate in the flat part of g(r)
      rdf_err = rdf_tot_err*n_atom/at_occup_r              ! each single-atom partial error is n_atom-times larger
     
C *** Accumulated! set initial weights, plot range and partial PDFs

			allocate(r(n_pdf),q(n_pdf),w(n_pdf),x_ffq(n_atom,n_pdf))


      do i=1,n_pdf
        w(i) = .5*(1.+cos(pi*i/real(n_pdf)))       ! FT window == Hann profile
CC					t_wind(i) = 0.355768-0.487396*cos(twopi*(i)/real(n_int+1))							!max is at n_int/2+1		Nuttall window
CC					t_wind(i) = t_wind(i)+0.144232*cos(2.*twopi*(i)/real(n_int+1))						
CC					t_wind(i) = t_wind(i)-0.012604*cos(3.*twopi*(i)/real(n_int+1))							
C				w(i) = 0.355768+0.487396*cos(pi*(i)/real(n_pdf))							!max is at i=0		Nuttall window
C				w(i) = w(i)+0.144232*cos(2.*pi*(i)/real(n_pdf))						
C				w(i) = w(i)+0.012604*cos(3.*pi*(i)/real(n_pdf))							
      enddo
C     w = 1.
C     write(*,*) 'w',w(1),w(n_pdf/2),w(n_pdf)
      
      q_step = 2.*pi/(n_pdf*pdf_step)
      do i=1,n_pdf
				r(i) = (i-1)*pdf_step
				q(i) = (i-1)*q_step
			enddo
			q(1) = 1.e-8        !to avoid singularities 
			r(1) = 1.e-8        !to avoid singularities 
			
C		write(*,*) 'q_step,q(1),q(n_pdf)',q_step,q(1),q(n_pdf)

C *** generate the Xray form-factor table
      if(j_weight==3) then
        do i=1,n_pdf
          do j=1,n_atom
            x_ffq(j,i) = x_ffpar(j,9)
            do ii=1,4
              x_ffq(j,i) = x_ffq(j,i)+x_ffpar(j,2*ii-1)*exp(-.25*x_ffpar(j,2*ii)*q(i)**2)
            enddo
          enddo
        enddo
      endif

C *** set atom weights
			at_mask = 1
      if(j_weight==1) then
        at_weight(1:n_atom) = 1.
      elseif(j_weight==2) then
        at_weight(1:n_atom) = b_coh(1:n_atom)
      elseif(j_weight==3) then
        at_weight(1:n_atom) = x_ffq(1:n_atom,1)
      endif


      if(n_part==0) ind_part = 0
      part_scale = 1.
      
 			at_weights_loop: do		
			
C *** extend the RDF matrix by pseudo_atom rows&columns 
        n_atom_tot = n_atom+n_pseudo             
        allocate(rdf_p2_ws(n_atom,n_atom,n_pdf),rdf_p2_plot(n_atom_tot,n_atom_tot,n_pdf),at_name_plot(n_atom_tot),rdf_fft(2*n_pdf)) 

        rdf_p2_ws = .0
        rdf_p2_plot = .0

        at_name_plot(1:n_atom) = at_name_par       
        if(n_pseudo>0) at_name_plot(n_atom+1:n_atom_tot) = at_name_pseudo(1:n_pseudo)
												
				write(9,*)
				write(*,'(1x,"Atoms:         ",51(1x,a8))')  (at_name_plot(i),i=1,n_atom+1)
				write(9,'(1x,"Atoms:         ",51(1x,a8))')  (at_name_plot(i),i=1,n_atom+1)
				write(*,'(1x,"Atoms no.:  ",51(1x,i8))') (i,i=1,n_atom+1)
				write(9,'(1x,"Atoms no.:  ",51(1x,i8))') (i,i=1,n_atom+1)
				write(*,'(1x,"Occupancy:       ",51(1x,f8.4))') (at_occup_r(i),i=1,n_atom),1.                      !1. for TOT
				write(9,'(1x,"Occupancy:       ",51(1x,f8.4))') (at_occup_r(i),i=1,n_atom),1.
				write(*,113) trim(at_weight_scheme(j_weight)),(sum(ind_pseudo(:,j)*at_weight),j=1,n_atom+1)
				write(9,113) trim(at_weight_scheme(j_weight)),(sum(ind_pseudo(:,j)*at_weight),j=1,n_atom+1)
113     format(1x,a,' weights: ',50(1x,f8.4))	
        write(*,'(1x,"Mean rel_error:  ",50(1x,f8.4))') (rdf_err(i),i=1,n_atom),rdf_tot_err						
				write(*,*) 
				write(*,*) 'Actual masks:'
				write(*,'((50i3))') (at_mask(i),i=1,n_atom)        
				write(*,*) 
        
        if(minval(at_mask(1:n_atom))==1) then  !if all masks =1 don't put them into plot title
          masks = ''
        else
          write(masks,'("  Masks: ",50i1.1)') (at_mask(i),i=1,n_atom)			!to be used in plot titles
        endif

 			  at_weight_av = sum(at_weight*at_mask*at_occup_r)/sum(at_occup_r)     !average weight (scattering length) per atom
				
        do jj=1,n_atom
          do ii=1,n_atom
						at_weight_matrix(ii,jj) = at_weight(ii)*at_mask(ii)*at_weight(jj)*at_mask(jj)
				  enddo
				enddo
 			  at_weight_matrix = at_weight_matrix/at_weight_av**2
				
        if(j_verb==1) then
          write(*,*) 'at_weight average',at_weight_av
          write(*,*) 'at_weight matrix'
          do ii=1,n_atom
              write(*,*) ii,at_weight_matrix(ii,:) 
          enddo      
          write(*,*) 'sum at_weight_matrix',sum(at_weight_matrix)
          write(*,*) 
        endif
				       
        write(*,*) 'Gaussian smooth FWHM in pdf_steps (1 no smoothing, 4 - 20 useful)'
        read(*,*) n_smooth_fwhm
			
				if(n_smooth_fwhm==1) then
				  n_smooth=1			               !no smoothing
				else
  				n_smooth = 5*n_smooth_fwhm/2+1
				endif
				
				allocate(f_smooth(n_smooth))

				do j=1,n_smooth
					f_smooth(j) = 2.**(-((j-n_smooth/2-1)/(.5*n_smooth_fwhm))**2) 
				enddo
				f_smooth = f_smooth/sum(f_smooth(1:n_smooth))				!profile normalized to unit integral
C			  write(*,*) 'f_smooth', (f_smooth(j),j=1,n_smooth)		

				if(n_smooth>1) write(*,*) 'Applying Gaussian smoothing with FWHM=',n_smooth_fwhm,' steps of',pdf_step,'[A^]'
        write(smooth,'("Smooth FWHM",f5.2," [A]")') pdf_step*n_smooth_fwhm

        if(n_smooth==1) then
          do i = n_pix_step,n_pdf
            rdf_p2_ws(:,:,i) = rdf_p2(:,:,i)*at_weight_matrix         
          enddo               
        else
          do i = n_pix_step,n_pdf
            do j = 1,n_smooth
              if (i-n_smooth/2+j-1.le.0) cycle  !zeros out of range
              if (i-n_smooth/2+j-1.gt.n_pdf) then
                rdf_p2_ws(:,:,i)=rdf_p2_ws(:,:,i)+f_smooth(j)*rdf_p2(:,:,n_pdf)     !replace out-of-range points by the asymptote, placed there at the normalisation
              else
                rdf_p2_ws(:,:,i)=rdf_p2_ws(:,:,i)+rdf_p2(:,:,i-n_smooth/2+j-1)*f_smooth(j)
              endif
            enddo 
            rdf_p2_ws(:,:,i) = rdf_p2_ws(:,:,i)*at_weight_matrix            !now we allow for b_coh in the g(r), it won't change the overall norm, just redistribute intensity between poartials
          enddo
        endif

				deallocate(f_smooth)

        ro_0 = sum(at_mask*at_occup_r)/product(a_par)                  !ro_0 is numeric density
        write(*,*) 'Atom density ro_0 [Å-3]',ro_0
        write(*,*) 

        if(j_pdf>=2) then
          do j=1,n_atom
            do i=1,n_atom  
              rdf_p2_ws(i,j,:) = 4.*pi*ro_0*at_mask(i)*at_occup_r(i)*at_mask(j)*at_occup_r(j)
     1              *r*(rdf_p2_ws(i,j,:)-at_weight_matrix(i,j)/sum(at_weight_matrix))    !for j_pdf>=2 make it G(r)
            enddo
          enddo
        endif

C *** calculate S(Q) components & total
        if(j_pdf==3) then
          
          do j=1,n_atom
            do i=1,n_atom  
                rdf_p2_ws(i,j,:) = fft((0.,1.)*w*rdf_p2_ws(i,j,:),inv=.false.)*pdf_step           !pdf_step = pdf_range/n_pdf with 1/sqrt(n_pdf) coming once from FT and and once from sampling norm
                rdf_p2_ws(i,j,1:n_pdf/2) = at_weight_matrix(i,j)/sum(at_weight_matrix)+rdf_p2_ws(i,j,2:n_pdf/2+1)/q(2:n_pdf/2+1)      !skip Q=0
                rdf_p2_ws(i,j,n_pdf/2+1:n_pdf) = at_weight_matrix(i,j)/sum(at_weight_matrix)
            enddo
          enddo
        endif

C *** produce rdf_p2_plot and sum up symmetric off-diagonal elements         

        do i=1,n_pdf
          rdf_p2_plot(:,:,i) = matmul(matmul(transpose(ind_pseudo(:,1:n_atom_tot)),rdf_p2_ws(:,:,i)),ind_pseudo(:,1:n_atom_tot))
        enddo

        do j=1,n_atom_tot
          do i=1,j  
            if(i/=j) rdf_p2_plot(i,j,:) = rdf_p2_plot(i,j,:)+rdf_p2_plot(j,i,:)
            if(i/=j) rdf_p2_plot(j,i,:) = 0.
            if(i/=j.and.i>n_atom.and.j>n_atom) rdf_p2_plot(i,j,:) = 0.
           enddo
        enddo

C       if(j_verb==1) then
C         do j=1,4
C           write(*,*) 'rdf_tot_plot mean',j,sum(rdf_p2_plot(n_atom+1,n_atom+1,1+(j-1)*n_pdf/4:j*n_pdf/4)/real(n_pdf/4))
C         enddo
C       endif
        
        if(j_pdf<=2) then
          x => r
        else        
          x => q
        endif
 
        x_start = x(1)
        i_start = 1
        i_end = (x_end-x(1))/(x(2)-x(1))
        n_plot = i_end-i_start+1
        
C *** open the PGPLOT graphics window (X11)
      if (j_xserv.LE.0) then          
        j_xserv = PGOPEN('/xserv')
      else
        call PGSLCT(j_xserv)
      endif   
      if (j_xserv.LE.0) then    
      	write(*,*) 'Could not open PGPLOT /xserv'
      	STOP
      endif
      
      CALL PGPAP(10.0,.6)     ! define the plot area as landscape
	  	call PGSUBP(1,1)				! PGPLOT window of 1x1 no of panes
      CALL PGASK(.FALSE.)     ! would not ask for <RET>
	  	CALL PGSCRN(0, 'white', IER)	!plot on white background
	  	CALL PGSCRN(1, 'black', IER)
	  	
C *** Set JK colors for line plots								  
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
 
        if(c_max==.0) then          !only do this at the real start of plot_loop
          c_min = minval(rdf_p2_plot(n_atom+1,n_atom+1,i_start+10:i_end))        !avoid maximum close to the origin
          c_min = .1*(int(10*c_min)-1)
          c_max = maxval(rdf_p2_plot(n_atom+1,n_atom+1,i_start+10:i_end))
          c_max = .1*(int(10*c_max)+1)
        endif
              
        write(*,*) 'Vertical scale c_min,c_max',c_min,c_max
        
        write(plot_header,'(a,"    ",a," weights  ")') trim(file_dat_t0),trim(at_weight_scheme(j_weight))
        plot_header = trim(subst_name)//'  '//trim(pdf_out(j_pdf))//'  '//trim(plot_header)//trim(masks)//'  '//trim(smooth)

				scale_loop: do
				   
					CALL PGSLCT(j_xserv)
          CALL PGSCI (1)  !white
          CALL PGSCH(1.)
          CALL PGSLW(2)
          CALL PGSLS (1)  !full
          CALL PGENV(x_start,x_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
          CALL PGLAB(x_label(j_pdf),y_label(j_pdf),trim(plot_header))  !put the axis labels
          CALL PGSCI (1)  !white
          CALL PGSLW(5)			!operates in steps of 5
          x_plot = x_start+.75*(x_end-x_start)
          y_plot = .9*c_max
          CALL PGLINE(n_plot,x(i_start:i_end),rdf_p2_plot(n_atom+1,n_atom+1,i_start:i_end))  !plots the curve
         
          CALL PGSTBG(0)																				 !erase graphics under text
C         CALL PGTEXT (x_plot,y_plot,trim(curve_label(j_pdf))//'    x 1.0')
          CALL PGTEXT (x_plot,y_plot,'Total     x 1.0')
          CALL PGSLW(2)

          do j=1,n_part
            if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
           	write(plot_title,115) at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
          	y_plot = (.9-.06*j)*c_max
115       	format(a,a,' x',f4.1)   
            CALL PGSCI (j+20)  !my red-green-blue
            CALL PGLINE(n_plot,x(i_start:i_end),part_scale(j)*rdf_p2_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title)
          	CALL PGSLW(2)			
          enddo
          
          write(*,*) 'Adjust vertical scale (min, max) (0 0 EXIT, -1 -1 to adjust plot range)'
          c1 = c_min
          c2 = c_max
          read(*,*) c1,c2
          if(c1==.0.and.c2==.0) then
            exit
          elseif(c1==-1.and.c2==-1) then              
            write(*,*) 'Confirm/adjust plot range, max =',n_pdf*(x(2)-x(1))
            write(*,'("x_start, x_end ",2f7.1,":  ")',advance='no') x_start,x_end 
            read(*,*) x_start,x_end 
            if(x_start.lt.pdf_step) x_start = x(1)
            i_start = (x_start-x(1))/(x(2)-x(1)) + 1
            i_end = anint((x_end-x(1))/(x(2)-x(1)))
            n_plot = i_end-i_start+1
          else
            c_min = c1
            c_max = c2
          endif
        enddo scale_loop
         
C **** Prepare and plot the same on .PS, look for existing output files in order not overwrite them				
				if(j_ps==1) then
					jfile = 1
					do						!look for existing .ps files to continue numbering
            if(j_name==0.and.t_single)then
              write(file_ps,1041) trim(file_master),jfile,trim(pg_ext)
1041   		format(a,'_rdf','_',i2.2,a)      
              write(file_res,1042) trim(file_master),jfile
1042   		format(a,'_rdf','_',i2.2,'.txt')      
            else
              write(c_jfile,'("_",i2.2)') jfile
              if(nfile>nfile_min) then
                if(nfile_min<=9999) then
                  write(c_nfile_min,'(i4.4)') nfile_min
                elseif(nfile_min>=10000) then
                  write(c_nfile_min,'(i8)') nfile_min
                endif
                c_nfile_min = '_'//adjustl(c_nfile_min)
              else
                c_nfile_min = ''
              endif
    
              if(nfile<=9999) then
                write(c_nfile,'(i4.4)') nfile
              elseif(nfile>=10000) then
                write(c_nfile,'(i8)') nfile
              endif
              c_nfile = '_'//adjustl(c_nfile)
    
              file_res = trim(file_master)//'_pdf'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//'.txt'							
              file_ps  = trim(file_master)//'_pdf'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//trim(pg_ext)
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
							
					ier = PGOPEN(file_ps//"/"//trim(pg_out))
					IF (ier.LE.0) STOP
					CALL PGASK(.FALSE.)     ! would not ask for <RET>
 			  	CALL PGPAP(11.0,.6)     ! define the plot area as landscape
					CALL PGSUBP(1,1)				! PGPLOT window of 1x1 no of panes
					CALL PGSCRN(0, 'white', IER)	!plot on white background
					CALL PGSCRN(1, 'black', IER)

					plot_header = trim(subst_name)//'  '//file_ps//'  '//trim(at_weight_scheme(j_weight))//' weights  '//trim(masks)//'  '//trim(smooth)	
 
          CALL PGSCI (1)  !white
          CALL PGSCH(1.)
          CALL PGSLW(2)
          CALL PGSLS (1)  !full
          CALL PGENV(x_start,x_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes w/o grid
          CALL PGLAB(trim(x_label(j_pdf)), trim(y_label(j_pdf)),trim(plot_header))  !put the axis labels
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW(5)			!operates in steps of 5
				  CALL PGLINE(n_plot,x(i_start:i_end),rdf_p2_plot(n_atom+1,n_atom+1,i_start:i_end))  !plots the curve
					y_plot = .9*c_max
					CALL PGSTBG(0)																				 !erase graphics under text
C		  		CALL PGTEXT (x_plot,y_plot,trim(curve_label(j_pdf))//trim(subst_name)) !the x_plot,y_plot are in world (axis units) coordinates
          CALL PGTEXT (x_plot,y_plot,'Total     x 1.0')
          CALL PGSLW(2)

          do j=1,n_part
            if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
            write(plot_title,115) at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
          	y_plot = (.9-.06*j)*c_max
            CALL PGSCI (j+1)  !red-green-blue
            CALL PGLINE(n_plot,x(i_start:i_end),part_scale(j)*rdf_p2_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title)				!the x_plot,y_plot are in world (axis units) coordinates
          	CALL PGSLW(2)
          enddo    
					CALL PGCLOS
					write(*,*) ' Postscript output written to: ',file_ps	
					write(9,*)
					write(9,*) '  ',trim(int_mode),' integration',n_int,' pairs'
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
				    write(4,*) 'OMP processes  = ', proc_num_in									
            write(4,*) 'Integration ',trim(int_mode),'  ',n_int,'cell pairs, j_rand',j_rand
            write(4,*) 'Supercell size:',n_row																				
            write(4,*) 'Unit cell parameter:',a_par																				
            write(4,*) 'Atoms/unit_cell:',n_atom
            write(4,*) 	 'Atoms  :',(('    '//at_name_plot(i)),i=1,n_atom)
            write(4,'(1x,"Atom no. ",50i8)') (i,i=1,n_atom)
            write(4,125) (at_occup_r(i),i=1,n_atom)
            write(4,126) (at_weight(i),i=1,n_atom)
            write(4,127) (at_mask(i),i=1,n_atom)
            write(4,*)  'Smoothing FWHM [Å]',n_smooth_fwhm*pdf_step
125     format(1x,"Occup's: ",50f8.4)								
126     format(1x,'Weights: ',50f8.4)								
127     format(1x,'Masks	: ',50i8)								
            write(4,*)
            write(4,*) 'Output: ',pdf_out(j_pdf)

C *** write the PDFs
            if(j_out==1) then
              write(4,107) ((ii,jj,ii=1,jj),jj=1,n_atom_tot)
              write(4,106) trim(x_label(j_pdf)),((trim(at_name_plot(ii)),trim(at_name_plot(jj)),ii=1,jj),jj=1,n_atom_tot)
              do i=1,n_pdf
                write(4,105) r(i),rdf_p2_plot(n_atom+1,n_atom+1,i),((rdf_p2_plot(ii,jj,i),ii=1,jj),jj=1,n_atom_tot)
              enddo
            else
              write(4,107) (ind_part(1,j),ind_part(2,j),j=1,n_part)
              write(4,106) trim(x_label(j_pdf)),(trim(at_name_plot(ind_part(1,j))),trim(at_name_plot(ind_part(2,j))),j=1,n_part)
              do i=i_start,i_end
                write(4,105) r(i),rdf_p2_plot(n_atom+1,n_atom+1,i),(rdf_p2_plot(ind_part(1,j),ind_part(2,j),i),j=1,n_part)
              enddo
            endif
105   	format(1x,f7.3,1x,f8.3,32(2x,5f8.3))	
107   	format(1x,"Partials index  ",32(2x,5(2x,2i3)))	
C106   	format('    r[Å]  PDF_tot   ',32(2x,5(a,'_',a,3x)))	
106   	format('   ',a,'   Total     ',32(2x,5(a,'_',a,3x)))	
            close(4)
            write(*,*) ' Text output written to: ',file_res	  
            write(9,*) ' Text output written to: ',file_res	  

				  endif 																							!j_txt
				endif 																							!j_ps

C ***   all done, now decide what comes next
			
			 way_point: do
					write(*,*) 'Choose a PDF options (MODE is ',trim(pdf_out(j_pdf)),', FILE output is ',trim(ps_out(j_ps+1)),', SIZE is ',trim(size_out(j_out+1)),'):'
	        write(*,*) '       1   REPLOT the PDFs '
	        write(*,'("        2   select max ",i2," partial PDFs & replot")') n_part_max
	        write(*,*) '       3   adjust partial PDF scales & replot '
          write(*,*) '       4   modify atom WEIGHTS (',trim(at_weight_scheme(j_weight)),')'
          write(*,*) '       5   edit atom MASKS ',trim(masks(8:))
          write(*,'("        6   create/modify max ",i2," PSEUDO_ATOMS ")') n_pseudo_max
					write(*,*) '       7   toggle FILE output ',trim(ps_out(mod(j_ps+1,2)+1)),' (mind the J_TXT switch in .PAR)'
					write(*,*) '       8   toggle .TXT output SIZE to ',trim(size_out(mod(j_out+1,2)+1))
          write(*,*) '       9   RESTART with updated weights, masks & pseudo_atoms '
          write(*,*)        
          write(*,*) '       10  select the PDF MODE, actual: ',trim(pdf_out(j_pdf))
          
          write(*,*) '       0   EXIT'  

          read(*,*) jj
          
          select case(jj)
						case(1) 
							cycle plot_loop

						case(2) 
              write(*,*) '("Confirm/modify up to ", i2," pairs of partial/pseudo PDF indices (0 0 erase, -1 -1 skip the rest):")',n_part_max
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
                    write(*,*) trim(at_name_pseudo(i)),n_atom+i,'  masks:',ind_pseudo(:,n_atom+i)
                  enddo
                endif
                
                write(*,'("Modify (number), create (max+1) or exit (0): ")',advance = 'no')
                read(*,*) j
                
                if(j>n_pseudo_max) then
                  write(*,*) 'Maximum number of pseudo-atoms reached, consider modifying existing ones'
                  write(*,'("Modify (number): ")',advance = 'no')
                  read(*,*) j
                  j = j-n_atom
                endif
                  
                if(j<=0) then
                  exit pseudo_loop
                elseif(j==1) then
                  write(*,*) 'The TOT pseudo cannot be modified'
                  cycle pseudo_loop
                elseif(j>1.and.j<=n_pseudo) then
                  write(*,'("New indices: ")',advance = 'no')
                  read(*,*)ind_pseudo(:,j)
                  cycle pseudo_loop
                else
                  write(*,'("Type pseudo_atom name (<=4 char): ")',advance = 'no')
                  read(*,*) at_name_pseudo(n_pseudo+1)
                  write(*,'("Type pseudo_atom indices (n_atom): ")',advance = 'no')
                  read(*,*)ind_pseudo(:,n_pseudo+1)
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

						case(10) 
              write(*,*) 'Select the PDF MODE: ',(j,'  ',trim(pdf_out(j)),j=1,n_mode)
              read(*,*) j_pdf
              exit plot_loop

						case(0) 
              exit at_weights_loop
					end select

				enddo way_point
			
        enddo plot_loop

        deallocate(rdf_p2_plot,rdf_p2_ws,at_name_plot,rdf_fft)  			
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
      
C **********************************************************************************************************
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

C **** string conversion to all lower case
C          
 			subroutine down_case (string)

			character(*), intent(inout)	:: string
			integer											:: j, nc
			character(len=26), parameter	:: lower = 'abcdefghijklmnopqrstuvwxyz'
			character(len=26), parameter	:: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

			do j = 1, len(string)
				nc = index(upper, string(j:j))
				if (nc > 0) string(j:j) = lower(nc:nc)
			end do

			end subroutine down_case	     
     
