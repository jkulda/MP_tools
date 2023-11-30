  
program mp_pdf55

! *************************************************************************************
! *****
! ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! *****
! *****   calculates the pair distribution functions (PDF) for simulated supercell data
! *****
!**********					Copyright (C) 2023  Jiri Kulda, Grenoble/Prague          **********
!**	
!** This file is part MP_TOOLS developed and maintained by Jiri Kulda <jkulda@free.fr>
!**
!**	MP_TOOLS are free software: you can use it, redistribute it and/or modify it 
!**	under the terms of the GNU General Public License as published by the Free Software 
!**	Foundation, either version 3 of the License, or (at your option) any later version, 
!**	for details see <https://www.gnu.org/licenses/>
!**
!**	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!**	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!**	See the GNU General Public License for more details.
!**
! ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! *****
! ***** Ver. 1.1 The correlation is asymmetric: n_int random cell against the complete supercell
! ***** Ver. 1.2 the whole supercell is read into memory once for ever
! ***** Ver. 1.21 the correlation is symmetric n_int random cells against n_int random cells
! *****
! ***** Ver. 1.31 between two delayed snapshot frame uses the 1.21 symmetric algorithm
! ***** Ver. 1.31 now the delay may be 0, so 1.31 replaces also ver. 1.21 and will override it
! *****
! ***** Ver. 1.32 the IDOM dialog moved after reading first input file
! ***** Ver. 1.32 the TIME_LOOP replaces the MASTER_LOOP
! ***** Ver. 1.32 normalization of correlation functions by true number of active cells n_norm = 45**3
! ***** 
! ***** Ver. 1.33 the reference snapshot and the first one in the time-series defined independently
! ***** Ver. 1.33 hence instantaneous (lag 0) correlations can be calculated explicitly
! ***** 
! ***** Ver. 1.4 	- flexible cell structure
! ***** 				 	- use of auxiliary <file_title.par> file (correlation pairs etc.)
! ***** Ver. 1.41	- corrected inconsistency in output <rdf>.txt file name
! *****						- correct histogram identification by j_hist_type
! ***** 
! ***** Ver. 1.42	- lighter & parallelised version (no correlation histograms, no IDOM)
! *****						- integrates the PGPLOT graphics of mp_rplot.f
! ***** 					- only PDF
! ***** 					- for correlation histograms (refer to mp_corr)
! ***** 					- modified .PAR file format, consequent .LOG
! ***** 
! ***** Ver. 1.50	- start of new series (identical to mp_pdf 1.42) as MP_PDF
! ***** 
! ***** Ver. 1.51	- all data arrays allocatable, no predefined array size limits
! ***** 					- supercell format in .PAR
! *****						- number of correlation pairs read from .PAR (don't use ver. 1.50!)
! *****
! ***** Ver 1.53  - record length 4096, all data 32bit length
! *****						- at_mass is saved as at_veloc(4,j,k)
! *****						- at_charge is saved as at_pos(4,j,k)
! *****           - CELL and QUICK data type (BULK not eligible for PDF)
! *****						- secondary sampling volume clipped for highly anisotropic or 2D supercells
! *****
! ***** 
! ***** Ver 1.54  - NAMELIST I/O for .PAR files and .DAT headers implemented 
! *****           - PNG driver for PGPLOT included
! *****           
! ***** Ver 1.55  - completely rewritten MC accumulation part  
! *****           - consistent MC results between OMP and non-OMP runs
! *****           - estimate of MC statistical accuracy
! *****           - automatic accumulation of the complete correlation landscape
! *****           - pseudoatoms to represent group correlations
! *****           - overlay with adjustable partial PDF scales
! ***** 
! ***** Ver 1.56  - completely rewritten g(r)/S(Q) part  
! *****           - Q-space results consistent with MP_SQL6
! *****            
! ***** file name convention:
! *****     <master_file>_n<snapshot number>.dat
! *****     example: H05_10K_n0001.dat from H05_10K.txt
! *****
! ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
! ***** record numbers are directly related to cell positions and atom types
! *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
! *****						1 going for the header record
! *****
! ***** atom positions are converted from reduced lattice coordinates (x) to real space distances
! ***** 
  
  use omp_lib
  use singleton

! *** notes on FT & SINGLETON:
!                 - check whether automatic scaling by 1/sqrt(N) inside SINGLETON is blocked or not - now it is!
!                 - if YES  you need to do it explicitly within the code	
!                 - the scaled FT transform restitutes properly f(x) after going forth and back
!                 - NEVERTHELESS the FT value will depend on the sampling frequency as SQRT(N)
!                 - if the FT value itself is important and should be invariant w.r.t. the sampling frequency the sample sequence has to be scaled by 1/SQRT(N) as well
  		
!!  		use mp_nopgplot					! uncomment when not able to use PGPLOT, compile and link together with the module mp_nopgplot.f

  integer,parameter :: l_rec  =  1024		    !record length in real(4)
  integer,parameter :: n_mc  =  1e6		      !unit of MC events count
  real(8), parameter   :: pi=3.141592653589793238462643383279502884197
  real(8), parameter   :: twopi = 2.d0*pi
  character(4),allocatable  ::	at_name_par(:),at_label(:),at_name_ext(:),at_name_pseudo(:),at_name_plot(:),pdf_out(:)
  character(16),allocatable ::	curve_label(:),x_label(:),y_label(:)

  integer,allocatable ::  at_mask(:),ind_part(:,:),ind_ext(:),ind_pseudo(:,:),numbers(:),ind_hist(:,:,:),ind_pdf(:,:,:),ind_at1(:),l_rdf_p2(:,:,:,:),n_cut(:)

  real(8),allocatable ::  rdf_p2(:,:,:),rdf_p2_n(:,:,:)                                               !the histogram has to be real(8) (long mantissa) to avoid conversion errors from integer
  real,allocatable ::  ft_wind(:),rdf_p2_ws(:,:,:),rdf_p2_plot(:,:,:),rdf_err(:),rdf_fft(:),x_ffq_av_sq(:),x_ffq_sq_av(:),f_smooth(:)
  real,allocatable ::  r_cut(:),x_ext(:),y_ext(:,:),ext_scale(:),ext_dy(:)
  real,allocatable ::  at_base(:,:),b_coh(:),at_weight(:),x_ffpar(:,:),x_ffq(:,:),at_weight_matrix(:,:),at_mask_matrix(:,:),at_av_matrix(:,:),rdf_norm(:,:)
  
  real,allocatable,target ::  r(:),q(:)
  real,pointer ::      x(:)

  integer ::       i_rdf,n_pdf,n_pdf4,j_acc,n_int,n_pseudo,i_ref,j_rand,n_ind(3),n_skip0,n_skip1,n_skip2,n_skip2_tot
  real    ::       rdf_dist,pdf_range,pdf_range_ext,pdf_range_sq_ext,rdf_sum,rdf_tot_sum,at_weight_av_sq,at_weight_sq_av,at_pos(3),at_pos2(3),base(3),q_xff
  real 		:: 			 arg,c_min,c_max,c_max2,c1,c2,x_start,x_end,pdf_step,q_step,ro_0,c_smooth,rnd(5),data_in(128)
  
  character(4)   :: version,head,atom,ps_out(2),size_out(2)
  character(10)  :: at_weight_scheme(3),pg_out,pg_ext,section,c_date,c_time,c_zone,c_nfile_min,c_nfile,c_jfile
  character(16)  :: sim_type_par,string16,filter_name
  character(40)  :: subst_name,file_master,file_inp,file_out,time_stamp,int_mode,x_file_name,string,mp_tool
  character(60)  :: file_dat,file_dat_t0,file_res,file_ps,file_log,masks,smooth
  character(256) :: cwd_path,plot_header,plot_title,line
  character(l_rec):: header_record
  
  logical ::  nml_in,found,found_txt,found_ps,t_single,single_bulk

  integer(8) ::  n_norm,n_h_max                !these numbers may be biiiiig!   
  integer ::  j_proc,proc_num,proc_num_in,thread_num,thread_num_max,m,jm,j1m,j2m,j_mode,n_mode
  integer ::  i_start,i_end,n_step,nfile_step,n_h,j_head_in,hist_ind(3),j_verb,n_tot,n_head
  integer ::  i,j,k,ii,jj,ind,ind2,jat,i_rec,ind_rec,nrec,nfile,nfile_min,nfile_max,jfile,n_at1,n_at1_tot,n_at2,i_at2
  integer ::  j1,j2,j_plane,j_grid,j_logsc,j_ps,j_out,j_txt,i_seed,i_r1,i_r2,i_r11,i_r12,i_time(8)
  integer ::  at_ind(3),at_ind2(3),d_ind(3),d_ind_shift(3),j_dom,n_plot,n_smooth,n_part,n_part_max,n_part_ext,n_atom_tot,n_pseudo_max,n_pdf_grid(3)
  integer ::  j_atom,j_weight,j_xff,j_edit,j_mask,ifile,j_int,n_atom,sc_c2,sc_c1,ier,ios,n_pix_step,i_cut,n_ext_skip,n_x,j_x,j_y,n_line,j_ext 
  
  real :: t_dump,filter_fwhm,tt0,tt,b_sum,at_sum,rand1,rand2,sc_r,at_displ,p_size,rdf_tot_err,pdf_pix,pdf_pix_shift(3),a_par_grid(3),sum_sq,sum_sq1,sum_sq2,sum_sq3

  integer :: rand1_seed(8),rand2_seed(8),seed_size

  real :: pdf_grid_min(3),pdf_grid_max(3),diff_pos(3),diff_pos_norm,x_plot,y_plot,smooth_fwhm,tot_scale,d_pos(3)
  real :: a_cell(3,3),a_cell_inv(3,3),a_cell_grid(3,3),a_cell_grid_inv(3,3)

! **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
  character(4),allocatable :: at_name_out(:)
  integer(4),allocatable   :: at_ind_in(:),nsuper_r(:)	

  real(4),allocatable ::	at_pos_in(:),at_occup_r(:),at_occup_1(:),at_occup_2(:),at_occup_k1(:,:),at_occup_k2(:,:),r_min(:,:),part_scale(:),at_scf(:)

  character(16)  :: sim_type,dat_type,input_method,file_par,dat_source
  integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,idum,rand,mult,j_pgc
  real(4)        :: rec_zero(l_rec),t_ms,t_step,t0,t1,a_par(3),angle(3),a_par_pdf(3),temp

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
  namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
  namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
 
  namelist /mp_gen/ j_verb,j_proc       
  namelist /mp_out/ j_weight,j_logsc,j_ps,j_txt,p_size,j_grid,pg_out,j_out,j_pgc       
                      !general rule: namelists of tools should only contain their local parameters
                      !what is of global interest they should pass into data_header
  namelist /mp_pdf/ pdf_range,pdf_step,q_step,x_end,a_par_pdf,pdf_pix,pdf_pix_shift,j_rand,n_h,j_acc,j_mode,n_part_max,n_pseudo_max,n_cond,q_xff    
			
! **** PGPLOT stuff
  INTEGER :: PGOPEN,j_xserv, NXSUB, NYSUB
  REAL::     XTICK,YTICK
  CHARACTER(10) :: XOPT,YOPT

  j_xserv = 0
  
!!! ********************* Test of RANDOM_NUMBER generator (start) *******************************      
!!!        
!!!        - remove the CC to check possible local and processor-dependent effects
!!!        - default seed size 8, can be changed by
!!!        
!!!     CALL RANDOM_SEED(SIZE=K)               !Puts size of seed into K
!!!     CALL RANDOM_SEED(PUT = SEED (1 : K))   ! Define a new seed
!!!     CALL RANDOM_SEED(GET = old_seed)       ! Read the current seed
!!
!!      integer ::  old_seed(8),new_seed(8)		
!!      	
!!			old_seed = 0
!!			new_seed = (/1,2,3,4,5,6,7,8/)
!!			CALL RANDOM_SEED         !Processor initialization
!!      do i=1,4
!!      CALL RANDOM_SEED(GET = old_seed)  ! Read current seed
!!      call random_number(rand1)
!!      write(*,*) 'default',rand1,old_seed
!!      enddo
!!
!!      CALL RANDOM_SEED(PUT = new_seed) ! Define seed
!!      do i=1,4
!!      CALL RANDOM_SEED(GET = old_seed)  ! Read current seed
!!      call random_number(rand1)
!!      write(*,*) 'new',rand1,old_seed
!!      enddo
!!
!!      CALL RANDOM_SEED(PUT = new_seed) ! check the reproducibility
!!      do i=1,4
!!      CALL RANDOM_SEED(GET = old_seed)  ! Read current seed
!!      call random_number(rand1)
!!      write(*,*) 'new_again',rand1,old_seed
!!      enddo
!! 
!! C *** Output of test:
!! 
!! default  0.774232209      -328864355 -1812850565  1795335799  -423129627   258875061   746672867 -1013455889  1109822438
!! default  0.499739051     -1772520199  1041085153  -345388489  -275392182 -1581757171   246595309  -441046729 -1125993601
!! default  0.551364720       512484964  1716921112  -363648824 -1071663748 -1457211347  -619233579   394076467  -840249998
!! default  0.671572745      -592193939 -1415474521  -203200136 -1628592313  1380558838  1345367254  1710769802   823960254
!! new  0.471070886               1           2           3           4           5           6           7           8
!! new  0.117344737       -16064360  1080567055  2124207785  -289480891   105059873 -2032997881  -987088397  1102541325
!! new  0.357547939     -1494515791  2126487000  -379374510 -1783270911 -1145639077 -1358320318  -112593222  -960835836
!! new  0.318134785       957948723  -227147975 -1341437354  -415669684 -1377879149  2060778058   557735840  1455132063
!! new_again  0.471070886               1           2           3           4           5           6           7           8
!! new_again  0.117344737       -16064360  1080567055  2124207785  -289480891   105059873 -2032997881  -987088397  1102541325
!! new_again  0.357547939     -1494515791  2126487000  -379374510 -1783270911 -1145639077 -1358320318  -112593222  -960835836
!! new_again  0.318134785       957948723  -227147975 -1341437354  -415669684 -1377879149  2060778058   557735840  1455132063
!!     
!!! ********************* Test of RANDOM_NUMBER generator (end) *******************************      




! ********************* Initialization *******************************      
  version = '1.56'
  mp_tool = 'MP_PDF '//version

  write(*,*) '*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  write(*,*)
  
! ********************* Get a time stamp and open a .LOG file *******************************

  call getcwd(cwd_path)
  call date_and_time(c_date,c_time,c_zone,i_time)
  write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
  write(file_log,'("mp_tools_",a,".log")') trim(c_date)
  inquire(file=file_log,exist=found_txt)
  if(found_txt)then
    open (9,file=file_log,position='append')
  else
    open (9,file=file_log)
  endif

  write(9,*) 
  write(9,*) 
  write(9,*) trim(time_stamp),'  ',trim(mp_tool),'  ',trim(cwd_path)
  write(9,*) 
  
! *********************  OpenMP initialization start  *******************************      
!
!    cf. line ≈410 (using input from .PAR)
!
! ********************* OpenMP initialization end *******************************      

! *** other initialisations
  n_mode = 7
  allocate(pdf_out(n_mode),curve_label(n_mode),x_label(n_mode),y_label(n_mode))
  curve_label = (/'g_tot','G_tot','F_tot','S_tot','Z_tot','I_tot'/) 
  pdf_out = (/'g(r)','G(r)','F(Q)','S(Q)','Z(Q)','I(Q)','I(Q)'/)			
  y_label = pdf_out
  y_label(7) = 'I(Q) [unscaled]'
  x_label = (/'r [A]  ','r [A]  ','Q [A-1]','Q [A-1]','Q [A-1]','Q [A-1]','Q [A-1]'/)
  ps_out = (/'OFF ','ON  '/)			!PGPLOT
  size_out = (/'S   ','XXL '/)		!PGPLOT
  j_out = 1
  j_ext = 0
  
! *** Generate data file access
  write(*,*) 'Input data file_master: '
  read(*,*) file_master 
        
  if(j_verb==1)	then	
    write(*,*) 'Read data files number min, step, max (0 0 no numbers, single file): '
    read(*,*)   nfile_min,nfile_step,nfile_max
  else
    write(*,*) 'Read data files number min, max (0 0 no numbers, single file): '
    read(*,*)   nfile_min,nfile_max
    if(nfile_max<nfile_min) nfile_max = nfile_min
    nfile_step = 1
  endif
  t_single = nfile_min==0.and.nfile_max==0

  if(t_single)then
    nfile_min = 1		!for conformity with loop control conventions
    nfile_max = 1
    nfile_step = 1
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
  
! *** Read the header record 

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
  elseif(head.eq.'TIME'.or.head.eq.'STAT') then                                  !old w/o structure
    write(*,*)		'Input data:  ','old header format'
    nml_in = .false.
   read(1,rec=1) sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,n_cond,n_rec					
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
  allocate(at_base(n_atom,3),at_weight(n_atom),at_weight_matrix(n_atom,n_atom),at_mask(n_atom),at_mask_matrix(n_atom,n_atom),rdf_err(n_atom))
  allocate(b_coh(n_atom),x_ffpar(n_atom,9),SOURCE=.0)												!we need this to read .par
  at_mask = 1

  if(nml_in) then      !new structure with namelist
    i_rec = i_rec+1   							
    read(1,rec=i_rec) header_record
    read(header_record,nml=data_header_2)	
  else                                  !old w/o structure
    read(1,rec=1) sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell_out,n_cond,n_rec,n_tot,&			
 &    at_name_out,at_occup_r,nsuper_r																		
    if(head.eq.'TIME') sim_type = 'TIMESTEP'
    if (head.eq.'STAT') sim_type = 'STATIC'
  endif 
  
  single_bulk = n_row(1)==1.and.n_row(2)==1.and.n_row(3)==1

  if (n_head==4)then
    i_rec = i_rec+1   							
    read(1,rec=i_rec) header_record
    read(header_record,nml=data_header_3)	   !get a_cell, a_cell_inv
  else
    a_cell = .0
    a_cell_inv = .0
    do i=1,3                               !generate a_cell, a_cell_inv for orthogonal cases covered by earlier mp_tools versions
      if(single_bulk) then
        a_cell(i,i) = 1.
        a_cell_inv(i,i) = 1.
      else
        a_cell(i,i) = a_par(i)              ! for n_row==1 the cell is always cubic and the effective lattice parameter is 1Å
        a_cell_inv(i,i) = 1./a_par(i)
      endif
    enddo
  endif		  
  close(1)
  
! **** Read the auxiliary file <file_par.par> with structure parameters, atom names and further info
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
  
  at_weight_scheme(1) = 'Unit'
  at_weight_scheme(2) = 'Neutron'
  at_weight_scheme(3) = 'Xray'
!   at_weight_scheme(4) = 'FZ partial'

  n_h = 0
  j_acc = 1
  j_mode = 1
  x_end = 20.
  pdf_step = 0.02
  pdf_pix = .5
  a_par_pdf = 1.
  n_part_max = 4         !defaults for &mp_pdf
  n_pseudo_max = 0
  n_part_ext = 0
  j_rand = 1
  q_xff = .0
  rewind(4)
  read(4,nml=mp_pdf)
  
  allocate(ind_pseudo(n_atom,n_atom+n_pseudo_max+1),at_name_pseudo(n_pseudo_max+1))       !1st pseudo is TOT by default
  ind_pseudo = 0
  at_name_pseudo = ''
  at_name_pseudo(1) = 'TOT'
  ind_pseudo = 0
  do j=1,n_atom
    ind_pseudo(j,j) = 1
  enddo
  ind_pseudo(:,n_atom+1) = 1

  allocate(ind_part(2,n_part_max),ind_ext(n_part_max),part_scale(n_part_max),ext_scale(n_part_max),ext_dy(n_part_max),r_cut(n_part_max),n_cut(n_part_max))       !1st pseudo is TOT by default
  allocate(at_av_matrix(n_atom,n_atom))       !1st pseudo is TOT by default 
  ind_part = 0
  ind_ext = 0
  ext_scale = 1.
  ext_dy = .0
  
  if(j_acc==3) then
    write(*,*) 'Setting J_ACC = 2 to the recommended MC integration algorithm (check for other J_ACC choices in .PAR)'
    j_acc = 2
  endif
  
! *** Read the atom positions       
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

! *** Read the partial PDF definitions       
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
        read(4,*,iostat=ios) string,ind_part(:,j)      !,r_cut(j)	!indices of partial PDFs to be displayed reading until end of the list
        if(ios/=0) then
          n_part = j-1
          exit
        endif
        n_part = j
      enddo
    endif
  endif
  
  if(n_part==0) ind_part = 0
  part_scale = 1.
  tot_scale = 1.
      

  
! *** Read the pseudo atom definitions       
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
			  

! *** Check atom names against the .PAR input       
  do j=1,n_atom
    if(at_name_par(j)/=at_name_out(j)) then
      write(*,*) 'Not-matching  atom names in .PAR and .DAT: ',j,at_name_par(j),at_name_out(j)
      write(*,*) 'Prefer .DAT? (1/0)'
      read(*,*) ii
      if(ii==1) at_name_par = at_name_out
      exit 
    endif
  enddo			

! *** read neutron scattering lengths (always)
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
        
! *** read Xray formfactor parameters 
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
    write(*,*) 'check atom name spelling and the neutron_xs.txt table; use unit weights'
  enddo xff_loop
  close(4)


! *** write overview of atom data
  write(*,*) 'Input method:  ',input_method
  if(input_method=='CELL') at_occup_r = at_occup_r/sum(at_occup_r)

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

! **** use the .PAR value of a_par in case of input in LATTICE units      
  if(a_par(1)*a_par(2)*a_par(3)==1.and.a_par_pdf(1)*a_par_pdf(2)*a_par_pdf(3)/=0.) then
    a_par = a_par_pdf
    write(*,*) 'Setting a_par to',a_par
  endif

! **** setup the MC integration 
!
  if(n_h==0) then
    write(string,*) n_h_max
    write(*,*) 'MC sampling pairs per frame ([x 10^6]: '                !,trim(adjustl(string)),' max):'
    read(*,*)   n_h
  endif

! **** check the suitability of the chosen MC integration algorithm 
  if(n_cond==0.and.j_acc>0) then                         !box with periodic bound_cond
    write(*,*) 'ATTENTION: the simulated data are not subject to periodic boundary conditions!'
    write(*,*) 'ATTENTION: an adapted MC integration algorithm has to be adopted, setting *** j_acc = 0 ***'
    j_acc = 0
  endif

! **** check the PDF range and number of MC pairs n_h 

  if(j_acc>0) then                         !box with periodic bound_cond
    if(pdf_range>minval(a_par*n_row)) then
      write(*,*) 'PDF range',pdf_range,' exceeds box size',minval(a_par*n_row)
      write(*,*) 'Setting PDF range to',minval(a_par*n_row)
      pdf_range = minval(a_par*n_row)        
    endif
    n_h_max = real(n_tot)**2/(2.*n_mc)    !(1+n_tot/n_mc)*n_atom*product(n_row)/2 product(n_row)/2 is roughly the volume of sphere with radius of pdf_range_max !n_mc=1e6

    if(n_h.gt.n_h_max) then                   !n_h_max is in units of n_mc to avoid overflow for large boxes
      write(*,*) 'WARNING: n_h exceeds max number of 10^6 atom pairs ',n_h_max
      write(*,*) 'type in a reduced n_h (<',n_h_max,'):'
      read(*,*)   n_h
    endif

  else                                      !box with non-periodic bound_cond (j_acc==0)
    write(*,*) 'ATTENTION: simple MC w/o periodic boundary conditions!'
     write(*,*) 'a_par,pdf_range',a_par,pdf_range
    do
      if(pdf_range> .333*minval(a_par*n_row)) then
        write(*,*) 'Box with non-periodic boundary conditions:'
        write(*,*) 'PDF range',pdf_range,' exceeds 1/3 box size',.333*minval(a_par*n_row)
        pdf_range = int(.333*minval(a_par*n_row))        
        write(*,*) 'Setting PDF range to',pdf_range
      endif      

      n_at1 = n_tot*product(a_par-2.*pdf_range*(/1.,1.,1./))/product(a_par)    !this just a first estimate of n_at1
      n_at2 = n_tot*1.3333*pi*pdf_range**3/product(a_par)                      !this just a first estimate of n_at2_max
      n_h_max =  real(n_at1*n_at2)/n_mc               
      write(*,*) 'n_tot,n_at1,n_at2,n_h_max',n_tot,n_at1,n_at2,n_h_max

      if(n_h<=n_h_max) exit                       !n_h_max is in units of n_mc to avoid overflow for large boxes
      
      write(*,*) 'WARNING: n_h exceeds max number of 10^6 atom pairs!'
      write(*,*) 'Reduce pdf_range and/or n_h :',pdf_range,n_h_max
      read(*,*)   pdf_range,n_h
    enddo        
  endif
  pdf_range_sq_ext = (pdf_range+10*pdf_step)**2 ! square of pdf_range slightly extended to avoid cuttoff within the last pixel
  
  n_pdf = 4*int(pdf_range/(4.*pdf_step))+1
  pdf_range = (n_pdf-1)*pdf_step             !update pdf_range to be consistent with n_pdf
  n_int = n_h*n_mc                           !number of MC cycles per snapshot
  int_mode = 'Monte Carlo'
  write(*,*) trim(int_mode),' integration over',n_h,'*10^6 cell pairs'

 
! **** establish the overlay PDF_GRID 
  n_pdf_grid = nint(n_row/pdf_pix)     
  d_pos = pdf_pix                     !d_pos is overlay grid cell size in lattice units, has to be commensurate with the box size, should be commensurate with the cell size
  a_par_grid = a_par*pdf_pix           !a_par_grid is its lattice parameter in Å 
!      pdf_pix_shift                      !pdf_pix_shift specified in .PAR is offset of PDF overlay grid to lattice to get most atoms in centres of its cells
  a_cell = transpose(a_cell)          ! transpose a_cell to be able to do matmul(a_cell,lattice_vector) to transform lattice_vector to cartesian coordinates
  a_cell_grid = a_cell*pdf_pix
  a_cell_grid_inv = a_cell_inv/pdf_pix
  
  if(single_bulk) then
    a_cell_grid = .0
    a_cell_grid_inv = .0
    do i=1,3
      a_cell_grid(i,i) = a_par_grid(i)
      a_cell_grid_inv(i,i) = 1./a_par_grid(i)
    enddo
  endif
  
  if(j_verb==1) write(*,*) 'Overlay grid size, step and offset:',n_pdf_grid,d_pos,pdf_pix_shift
  if(j_verb==1) write(*,*) 'a_par_grid',a_par_grid

  n_pix_step = nint(norm2(a_par_grid)/pdf_step)       !we shall avoid correlations within the 1st pixel range
  if(j_verb==1) write(*,*) 'pdf_range,n_pdf,n_pix_step',pdf_range,n_pdf,n_pix_step
 
  
! *********************  OpenMP initialization start  *******************************      
!
  proc_num_in = j_proc
  thread_num_max = omp_get_max_threads( )			!this gives maximum number of threads available (limited by concurrent tasks??)
  proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
  if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
  call omp_set_num_threads(proc_num_in)
  thread_num = omp_get_num_threads( )				  !this gives threads available at this moment: here always = 1 as we are outside of a PARALLEL range


  if(proc_num_in.gt.1) then			
    write(*,*) 
    write(*,*) 'OMP threads max          = ', thread_num_max
    write(*,*) 'OMP processes requested  = ', proc_num_in										
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
!
! ********************* OpenMP initialization end *******************************      

! *** initialise the random_number generator
  call random_seed(size=seed_size)               !returns size of seed into seed_size
  if(j_verb==1) write(*,*) 'Random_numbers seed_size',seed_size
  allocate(numbers(seed_size))

  if(j_rand==0) then                     
    write(*,*) 'Random_number: j_rand =',j_rand,'  the system will supply unique, machine dependent seeds each time this code runs'
    call random_seed(get=numbers)               !Gets actual seeds 
    if(j_verb==1) then
      write(*,*) 'Random_number seed size:',seed_size
      write(*,*) 'Random_seed:',numbers
      write(*,*) 'Reference 1st 5 random numbers:',(rnd(i),i=1,5)
    endif
  elseif(j_rand==1) then                     !if j_rand==1 generate k-dependent standard seeds for each of the OMP threads
    write(*,*) 'Random_number: j_rand =',j_rand,'  the system will supply k-dependent standard seeds for each of the OMP threads (use only for testing the consistence of OMP_on/OMP_off results)'
  elseif(j_rand>1) then                     !if j_rand>1 generate a seed for later reference & numerical reproducibility checks
    write(*,*) 'Random_number: j_rand =',j_rand,'  this seeding reference can be used to exactly reproduce this MC-run later on'
    numbers(1) = rand(huge(1)/j_rand)   !this is to seed the RAND generator
    do i=1,seed_size
      numbers(i) = huge(1)*rand(0)    !use a trivial random number generator to produce the seeds (they could even be all the same small ones, but ...)
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

! *** Allocate and clear the histogram arrays for accumulation across several snapshots	       
  allocate (at_pos_in(4*n_tot),at_ind_in(4*n_tot),ind_at1(n_tot))
  allocate(ind_pdf(n_pdf_grid(1),n_pdf_grid(2),n_pdf_grid(3)))
  allocate(rdf_p2(n_atom,n_atom,n_pdf),rdf_p2_n(n_atom,n_atom,n_pdf),l_rdf_p2(n_atom,n_atom,n_pdf,n_h))
  allocate(at_occup_1(n_atom),at_occup_2(n_atom),at_occup_k1(n_atom,n_h),at_occup_k2(n_atom,n_h),r_min(n_atom,n_atom),rdf_norm(n_atom,n_atom))

  r_min = .0
  ind_pdf = 0
  rdf_p2 = .0
  l_rdf_p2 = 0

  call cpu_time(t0)
  
  if(t_single)then
    nfile_min = 1
    nfile_max = 1
    nfile_step = 1
  endif
  
! *************  the file loop: cycle over the tt files to augment statistics
  CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
  
  nfile = 0
  n_skip2_tot = 0
  n_cut = 0
  at_occup_1 = .0
  at_occup_2 = .0
  n_at1_tot = 0

  file_loop: do ifile=nfile_min,nfile_max,nfile_step									

    n_skip0 = 0
    n_skip1 = 0

! ***  open the t0 file (binary MD snapshot file)
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

!     CALL SYSTEM_CLOCK (COUNT = sc_c2)
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
     write(*,*) 'Read finished'


! **** shortlisting atom_1 candidates in a non-periodic box (j_acc==0)
    if(j_acc==0) then
      ind_at1 = 0
      n_at1 = 0

      if(single_bulk) then
        do i=1,n_tot
          if(sum(abs(at_pos_in(4*(i-1)+1:4*(i-1)+3)))==.0) cycle                     !skip void positions
          if(minval(.5*a_par-pdf_range-abs(at_pos_in(4*(i-1)+1:4*(i-1)+3)))>=.0) then
            n_at1 = n_at1+1
            ind_at1(n_at1) = i
          endif
        enddo
      else                                            ! i.e. n_row(:)==1
        do i=1,n_tot
          if(sum(abs(at_pos_in(4*(i-1)+1:4*(i-1)+3)))==.0) cycle                     !skip void positions
          if(minval(n_row/2.-pdf_range/a_par-abs(at_pos_in(4*(i-1)+1:4*(i-1)+3)))>=.0) then
            n_at1 = n_at1+1
            ind_at1(n_at1) = i
          endif
        enddo
      endif

      if(n_h<=nint((n_at1*n_at2)/real(n_mc))) then
        n_at2 = nint(n_h*n_mc/real(n_at1))
        n_int = n_at1*n_at2
      else
        write(*,*) 'WARNING: n_h exceeds max number of 10^6 atom pairs ',nint((n_at1*1000)/1.e6)
        write(*,*) 'restart with an adequate n_h and use a larger number of frames'
        stop
      endif
      write(*,*) 'Atom_1 & atom_2 pools:',n_at1,n_at2
      
      do i=1,n_at1
        at_occup_1(at_ind_in(4*(ind_at1(i)-1)+4)) = at_occup_1(at_ind_in(4*(ind_at1(i)-1)+4))+1
      enddo
      n_at1_tot = n_at1_tot+n_at1

    else
! *** reindexing atoms on a finer PDF overlay grid	
      ii =1		
      do i=1,n_tot
        if(maxval(abs(at_pos_in(4*(i-1)+1:4*(i-1)+3)))>.0) then
          if(single_bulk) then                                            ! i.e. n_row(:)==1  always orthorhombic
            d_ind = 1+nint((at_pos_in(4*(i-1)+1:4*(i-1)+3)/a_par+.5)/d_pos)   !+d_ind_shift
          else
            d_ind = 1+nint((at_pos_in(4*(i-1)+1:4*(i-1)+3)-pdf_pix_shift+n_row/2)/d_pos)   !+d_ind_shift
          endif

          do j=1,3                                                  !handle basis atoms displaced out of the nominal box
            if(d_ind(j)<=0) then
              d_ind(j) = d_ind(j)+n_pdf_grid(j)
              if(single_bulk)then
                at_pos_in(4*(i-1)+j) = at_pos_in(4*(i-1)+j)+a_par(j)
              else
                at_pos_in(4*(i-1)+j) = at_pos_in(4*(i-1)+j)+n_row(j)
              endif
            endif

            if(d_ind(j)>n_pdf_grid(j)) then
              d_ind(j) = d_ind(j)-n_pdf_grid(j)
               if(single_bulk)then
                at_pos_in(4*(i-1)+j) = at_pos_in(4*(i-1)+j)-a_par(j)
              else
                at_pos_in(4*(i-1)+j) = at_pos_in(4*(i-1)+j)-n_row(j)
              endif
            endif
          enddo

          do j=1,3                                                  !handle basis atoms displaced out of the nominal box
            if(d_ind(j)<=0.or.d_ind(j)>n_pdf_grid(j)) then
              write(*,*)'i,d_ind(:)',i,d_ind(:)
              read(*,*)
            endif
          enddo


         if(ind_pdf(d_ind(1),d_ind(2),d_ind(3))/= 0) then
            jj = ind_pdf(d_ind(1),d_ind(2),d_ind(3))
            write(*,*) 'Overlay grid cell already taken (you can accept a few by RETURN, else slightly decrease pdf_pix in .PAR):'
            write(*,*) 'at_ind_in,at_pos_in,d_ind',at_ind_in(4*(i-1)+1:4*(i-1)+4),'  ',at_pos_in(4*(i-1)+1:4*(i-1)+3),'  ',d_ind
            write(*,*) 'IN:jj,at_pos_in',jj,at_pos_in(4*(jj-1)+1:4*(jj-1)+3)
           read(*,*)
           n_skip0 = n_skip0+1
          endif
          ind_pdf(d_ind(1),d_ind(2),d_ind(3)) = i
          at_ind_in(4*(i-1)+1:4*(i-1)+3) = d_ind            !re-use at_ind_in(:,1:3) to store the new indices		
        endif	
      enddo
    endif
    
    CALL SYSTEM_CLOCK (COUNT = sc_c2)
    write(*,*) 
    if(j_verb==1) write(*,*) 'Input finished         ',(sc_c2-sc_c1)/sc_r,' sec'
    
    if(n_skip0>0) then
      write(*,*) 'Number of atoms skipped due to overlay grid cell double occupancy:', n_skip0
      if(n_skip0>10)then
        write(*,*) 'Modify .PAR pdf_pix (0), continue (1)?'
        read(*,*) jj
        if(jj==0) stop
      endif
    endif
  
! *************  the integration loop: cycle over the site pairs to accumulate the PDFs ******************

    write(*,*)    
    write(*,*) 'Accumulating the PDFs ...'     

! *** basic MC for small boxes without periodic boundary conditions (n_cond=0) - simply taking pairs of atoms
!
    if(j_acc==0) then 
      if(.not.single_bulk) at_pos_in = at_pos_in*a_par
      n_skip1 = 0
      n_skip2 = 0
      j_int = 0

      do k=1,n_at1                           ! has to run over ALL n_at1 - their order is not random a priori
        at_pos(:) = at_pos_in(4*(ind_at1(k)-1)+1:4*(ind_at1(k)-1)+3)       !take first atom from shortlist, the whole shortlist is to be taken as atom types may be ordered
        ii = at_ind_in(4*(ind_at1(k)-1)+4)

        i_at2 = 0
        do i=1,n_tot
          call random_number(rand1)
          i_r2 = n_tot*rand1+1
          at_pos2(:) = at_pos_in(4*(i_r2-1)+1:4*(i_r2-1)+3)       !randomly select second atom
          jj = at_ind_in(4*(i_r2-1)+4)

          if(sum((at_pos2-at_pos)*(at_pos2-at_pos))<=pdf_range_sq_ext) then
            i_rdf = nint(norm2(at_pos2-at_pos)/pdf_step)+1     
            if(i_rdf<=n_pdf) then 
              l_rdf_p2(ii,jj,i_rdf,1) = l_rdf_p2(ii,jj,i_rdf,1)+1 
              i_at2 = i_at2+1
              j_int = j_int+1
              at_occup_2(jj) = at_occup_2(jj)+1
              if(i_at2==n_at2) exit
            endif
          endif
        enddo    !n_at2
      enddo      !n_at1

     if(j_verb==1) write(*,*)'number of preset and really accumulated MC events',n_int,j_int
      rdf_p2 = rdf_p2+real(l_rdf_p2(:,:,:,1))
      l_rdf_p2 = .0
    endif     !j_acc==0


! *** MC generation of g(r) for large boxes with periodic boundary conditions (n_cond>=1), using the projective MC algorithm for J_ACC==2
!
    if(j_acc>=1) then
      n_at1_tot = n_h*n_mc/1000
      n_skip1 = 0
      n_skip2 = 0
      at_occup_k1 = .0
      at_occup_k2 = .0
      pdf_range_ext = pdf_range+1.4*norm2(a_par_grid)   !!! 1.4 is a security factor >1 to get the norm2() beyond rounding errors with some safety margin

! *** MonteCarlo integration

!$omp parallel shared(at_pos_in,at_ind_in,ind_pdf,l_rdf_p2,a_par,a_par_grid,a_cell,a_cell_inv,a_cell_grid,a_cell_grid_inv,&
!$omp& pdf_step,pdf_range_ext,n_tot,n_pdf_grid,n_row,n_atom,n_h,j_rand,n_skip1,n_skip2,n_at1_tot,at_occup_k1,at_occup_k2,single_bulk)&
!$omp& private(rnd,at_pos,at_pos2,diff_pos,diff_pos_norm,at_ind,at_ind2,d_ind,base,i_r1,i_rdf,idum,rand1,rand2,numbers,r,j_int,m,ii,jj)
!$omp do

    do k=1,n_h                           ! runs over 1000*1000 atom pairs for each k
      if(j_rand==1) then                 ! if j_rand=1 produce k-dependent standard seeds for each of the threads and cycles to test consistence of OMP_on and OMP_off results
        idum = j_rand+k
        numbers(1) = huge(1)*rand(idum)            !a dry call to initialise ran0
        do i=1,seed_size
          numbers(i) = huge(1)*rand(idum)          !use a trivial random number generator to produce the seeds (they could even be all the same small ones, but ...)
        enddo
        call random_seed(put=numbers)
        do i=1,5
          call random_number(rnd(i))
        enddo
        if(j_verb==1) then
          write(*,*) 'k =',k
          write(*,*) 'Random_seed:',numbers
          write(*,*) 'Reference 1st 5 random numbers:',(rnd(i),i=1,5)
        endif
      endif

      j_int = 0
      integration_loop_2: do              !runs over 1000 "first atoms"
        call random_number(rand1)			!the first cell
          i_r1 = n_tot*rand1+1	!i_r1 may become 1 to n_tot
          at_pos(:) = at_pos_in(4*(i_r1-1)+1:4*(i_r1-1)+3)       !randomly select first atom
          if(sum(abs(at_pos(:)))==.0) cycle integration_loop_2
          at_ind = at_ind_in(4*(i_r1-1)+1:4*(i_r1-1)+3)   !at_ind_in(:,1:3) has been re-used to store the overlay grid indices
          ii = at_ind_in(4*(i_r1-1)+4)
          at_occup_k1(ii,k) = at_occup_k1(ii,k)+1 

        m = 0
        m_cycle_2: do                 !runs over 1000 "second atoms"
          do
            do i=1,3
              call random_number(rand2)		!generate the position within a cube [-1,1]; only points within the unit sphere will be retained
              diff_pos(i) = 2.*rand2-1.
            enddo
            diff_pos_norm = sum(diff_pos*diff_pos)
            if(diff_pos_norm>.0.and.diff_pos_norm<=1.) exit     !only diff_pos within unit sphere is kept !throwing away also the unlikely event of diff_pos=.0, not d_ind=0!
          enddo



!!!!!!!!!!!!!!!!!! the whole difference between 'normal' and 'projective' MC sampling is here (plus normalisation ....) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(j_acc==2) then                                   
!               diff_pos = diff_pos/sqrt(diff_pos_norm)        !now we project random points onto spherical surface by d_ind/norm2(d_ind) ....
            call random_number(rand2)
            diff_pos = pdf_range_ext*rand2*diff_pos/sqrt(diff_pos_norm)      ! .... and select a random point on their radius (rand2), diff_pos is in Å now 
          else
            diff_pos = pdf_range_ext*diff_pos
          endif
!!!!!!!!!!!!!!!!!! the whole difference between 'normal' and 'projective' MC sampling is here (plus normalisation ....) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
!           d_ind = nint(diff_pos/a_par_grid)       ! .... convert the distance to overlay grid units  !use a_cell_inv to cover non-orthogonal cases
          d_ind = nint(matmul(a_cell_grid_inv,diff_pos))       ! .... convert the distance to overlay grid units  !use a_cell_inv to cover non-orthogonal cases
          
          if(sum(abs(d_ind))==0) cycle m_cycle_2
           
          at_ind2 = at_ind+d_ind          

! *** treat atoms out of the supercell by applying cyclic boundary conditions
          base = 0
          do j =1,3                                !n_pdf_grid is the supercell size on the overlay grid
            if(at_ind2(j)<1) then
              at_ind2(j) = at_ind2(j)+n_pdf_grid(j)
              if(single_bulk) then
                base(j) = -a_par(j)
              else
                base(j) = -n_row(j)
              endif
              if(at_ind2(j)<1) cycle m_cycle_2  !out-of-range due to security margin
            endif

            if(at_ind2(j)>n_pdf_grid(j)) then
              at_ind2(j) = at_ind2(j)-n_pdf_grid(j)
              if(single_bulk) then
                base(j) = a_par(j)
              else
                base(j) = n_row(j)
              endif
              if(at_ind2(j)>n_pdf_grid(j)) cycle m_cycle_2  !out-of-range due to security margin
            endif
          enddo

          do j =1,3                                !n_pdf_grid is the supercell size on the overlay grid
            if(at_ind2(j).lt.1.or.at_ind2(j).gt.n_pdf_grid(j)) then
              write(*,*)'m,at_ind,d_ind,at_ind2,diff_pos',m,at_ind,d_ind,at_ind2,diff_pos
              read(*,*)
            endif
          enddo

! *** retrieve the 2nd atom positions
          i = ind_pdf(at_ind2(1),at_ind2(2),at_ind2(3))
          if(i==0) then
            n_skip1 = n_skip1+1             
            cycle m_cycle_2                              !cycle over void indexing overlay grid
          endif
          at_pos2 = at_pos_in(4*(i-1)+1:4*(i-1)+3)             
          
          if(at_pos2(1)/=.0.and.at_pos2(2)/=.0.and.at_pos2(3)/=0.) at_pos2(:) = at_pos2(:)+base	!handle possible PBC effects
          jj = at_ind_in(4*(i-1)+4)

          at_pos2 = at_pos2-at_pos
          if(.not.single_bulk) at_pos2 = matmul(a_cell, at_pos2)      !use a_cell to cover non-orthogonal cases, here a_cell is already TRANSPOSED!

! *** accumulate all the PDF at once
       
          diff_pos_norm = sum(at_pos2*at_pos2)                  
          if(diff_pos_norm>=pdf_step.and.diff_pos_norm<=pdf_range_sq_ext) then         !contains a safety margin of .5A
            i_rdf = nint(sqrt(diff_pos_norm)/pdf_step)+1
            if(i_rdf<=n_pdf) then 
              l_rdf_p2(ii,jj,i_rdf,k) = l_rdf_p2(ii,jj,i_rdf,k)+1       
              m = m+1
              at_occup_k2(jj,k) = at_occup_k2(jj,k)+1
            endif
          else
            n_skip2 = n_skip2+1             
            cycle m_cycle_2
          endif
          if(m==1000) exit
        enddo m_cycle_2   !m
        j_int = j_int+1
        if(j_int==1000) exit integration_loop_2
      enddo integration_loop_2
    enddo		!k=1,n_h
!$omp end do
!$omp end parallel

      write(*,*)'MC hits to empty overlay grid cells [%]:',(100.*n_skip1)/(real(n_h*n_mc)+real(n_skip1))
      write(*,*)'MC hits out of PDF range [%]:',(100.*n_skip2)/(real(n_h*n_mc))
      n_skip2_tot = n_skip2_tot+n_skip2

! *** sum up contributions from individual cycles and estimate statistical error at the a_par(1) position

      do k=1,n_h  
        rdf_p2 = rdf_p2 + 1.*l_rdf_p2(:,:,:,k)
        at_occup_1 = at_occup_1+at_occup_k1(:,k)
        at_occup_2 = at_occup_2+at_occup_k2(:,k)
      enddo

      l_rdf_p2 = .0
      ind_pdf = 0
    endif  !j_acc=2
    
    if(j_verb==1) write(*,*) 'PDF size, MC events, file number, PDF histogram sum',n_pdf,n_int,ifile,sum(rdf_p2)       
            
  enddo file_loop

  
  CALL SYSTEM_CLOCK (COUNT = sc_c1)
  write(*,*) 
  write(*,*) 'Accumulation finished         ',(sc_c1-sc_c2)/sc_r,' sec'
  write(*,*) 

! *** estimate statistical error (Poisson statistics, sqrt(n)^-1) at the a_par(1) position
  n_pdf4 = (n_pdf-1)/4

  if(j_acc<=1) then
    rdf_tot_err = sqrt(real(n_pdf4-n_pix_step+1)/sum(rdf_p2(:,:,n_pix_step:n_pdf4)))         
  else
    rdf_tot_err = sqrt((n_pdf-1)*.5/sum(rdf_p2(:,:,(n_pdf-1)/2:n_pdf)))                   
  endif            
  rdf_err = rdf_tot_err/at_occup_r  


! *** normalise the accumulated PDF

  rdf_p2 = rdf_p2/real(nfile*n_int)
  
         
  at_occup_1 = at_occup_1/sum(at_occup_1)
  at_occup_2 = at_occup_2/sum(at_occup_2)


! *** Normalise the accumulated g_i(r)

  if(j_acc<=1) then
    arg = .0
    do i=2,n_pdf
      rdf_p2_n(:,:,i) = rdf_p2(:,:,i)/((i-1)**2)
      arg = arg+real(i-1)**2
    enddo
    rdf_p2_n = rdf_p2_n*arg
    
  elseif(j_acc==2) then            
    do i=2,n_pdf
      arg = 0.62035*a_par_grid(1)/((i-1)*pdf_step)              !correct for intercepted solid angle deviation from 1/R**2 !replacing pixel cube by sphere with equivalent volume r=a*(3/4/pi)**1/3
      arg = (atan(arg)/arg)**2     
      rdf_p2(:,:,i) = rdf_p2(:,:,i)*arg 
    enddo

    rdf_p2 = rdf_p2/sum(rdf_p2)
    
    do i=1,n_atom
      do j=1,n_atom
        i_cut = n_pix_step/2        !we are not going to care about what happens below, setting 0 there
        sum_sq1 = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut+1)/real(2*n_pdf4)) 
        i_cut = n_pix_step/2+1       !we are not going to care about what happens below, setting 0 there
        sum_sq2 = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut+1)/real(2*n_pdf4)) 
 
        if(sum_sq2<sum_sq1) then   
          do k=1,n_pdf/2
            sum_sq = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut+1)/real(2*n_pdf4))      

            if(sum_sq<sum_sq1) then
              sum_sq1 = sum_sq
              i_cut = i_cut+1
              cycle
            else
              i_cut = i_cut-1
              exit
            endif
            write(*,*) 'R_min cutoff not found, i_cut,sum_sq1',i_cut,sum_sq1
          enddo   !k
        else
          i_cut = n_pix_step/2-1
           do k=1,n_pdf/2
            sum_sq = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut+1)/real(2*n_pdf4))      

            if(sum_sq<sum_sq1) then
              sum_sq1 = sum_sq
              i_cut = i_cut-1
              cycle
            else
              i_cut = i_cut+1
              exit
            endif
            write(*,*) 'R_min cutoff not found, i_cut,sum_sq1',i_cut,sum_sq1
          enddo   !k
       endif
    
        sum_sq1 = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut)/real(2*n_pdf4))   !i_cut-1     !ex: n_pdf=1001,n_pdf4=250,pdf_range=(n_pdf-1)*pdf_step=1000*pdf_step
        sum_sq2 = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut+1)/real(2*n_pdf4)) !i_cut
        sum_sq3 = abs(at_occup_1(i)*at_occup_2(j)-sum(rdf_p2(i,j,2*n_pdf4+2:n_pdf))*(n_pdf-i_cut+2)/real(2*n_pdf4)) !i_cut+1   
        r_min(i,j) = i_cut-.5*(sum_sq3-sum_sq1)/(2.*sum_sq2-sum_sq1-sum_sq3)

        rdf_norm(i,j) = real(n_pdf)-r_min(i,j)+1.    
        r_min(i,j) = (r_min(i,j)-1.)*pdf_step
        
        rdf_p2_n(i,j,:) = rdf_p2(i,j,:)*rdf_norm(i,j)                    !multiplying by an "effective N_PDF", the rest will go into the 1st channel
        rdf_p2_n(i,j,1) = n_pdf*at_occup_1(i)*at_occup_2(j)-sum(rdf_p2_n(i,j,2:n_pdf))       !fill the first channel        
      enddo
    enddo

    if(j_verb==1) then
      write(*,*) 'PDF norm(i,j) [eff number of channels]'
      do ii=1,n_atom
          write(*,*) ii,rdf_norm(ii,:) 
      enddo      
      write(*,*) 

      write(*,*) 'R_min(i,j) [A]'
      do ii=1,n_atom
          write(*,*) ii,r_min(ii,:) 
      enddo      
      write(*,*)
    endif 

! *** analyse R_min cutoff effects
    if(j_verb==1) then    
      write(*,*) 'at1 at2    R_min    ' 
      do i=1,n_atom
        do j= 1,i
          write(*,*) at_name_par(i),at_name_par(j),.5*(r_min(i,j)+r_min(j,i))    
        enddo
      enddo
    endif

  endif     !(j_acc==2)

! *** Accumulated! set initial weights, plot range and partial PDFs

  allocate(r(n_pdf),q(n_pdf),ft_wind(n_pdf),x_ffq(n_atom,n_pdf),x_ffq_av_sq(n_pdf),x_ffq_sq_av(n_pdf),at_scf(n_atom))

  q_step = twopi/(2.*(n_pdf-1)*pdf_step)
  do i=1,n_pdf
    r(i) = (i-1)*pdf_step
    q(i) = (i-1)*q_step
  enddo
  q(1) = 1.e-8        !to avoid singularities 
  r(1) = 1.e-8        !to avoid singularities 
  
! *** generate the Xray form-factor table
    do i=1,n_pdf
      do j=1,n_atom
        x_ffq(j,i) = x_ffpar(j,9)
        do ii=1,4
          x_ffq(j,i) = x_ffq(j,i)+x_ffpar(j,2*ii-1)*exp(-.25*x_ffpar(j,2*ii)*(q(i)/twopi)**2)
        enddo
      enddo
      x_ffq_av_sq(i) = sum(at_occup_1(:)*x_ffq(:,i))*sum(at_occup_2(:)*x_ffq(:,i))
      x_ffq_sq_av(i) = sum(at_occup_1(:)*x_ffq(:,i)**2)      
    enddo

    j_xff = 1+nint(q_xff/q_step)            !reference to Q for effective XFF to be used in g(r) and G(r)

! *** set atom weights
  at_scf = 1.
  at_mask = 1
  if(j_weight==1) then
    at_weight(1:n_atom) = 1.
  elseif(j_weight==2) then
    at_weight(1:n_atom) = b_coh(1:n_atom)
  elseif(j_weight==3) then
    at_weight(1:n_atom) = x_ffq(1:n_atom,j_xff)
  endif


  at_weights_loop: do		

    ro_0 = sum(nsuper_r)/(product(a_par)*product(n_row))               !ro_0 is numeric density; N_TOT (number of lattice positions including the virtually unoccupied ones for mixed occup) can't be used here!
    
  
! *** extend the RDF matrix by pseudo_atom rows&columns 
    n_atom_tot = n_atom+n_pseudo             
    allocate(rdf_p2_ws(n_atom,n_atom,n_pdf),rdf_p2_plot(n_atom_tot,n_atom_tot,n_pdf),at_name_plot(n_atom_tot),rdf_fft(2*n_pdf-2)) 

    at_name_plot(1:n_atom) = at_name_par       
    if(n_pseudo>0) at_name_plot(n_atom+1:n_atom_tot) = at_name_pseudo(1:n_pseudo)

    write(*,*)
    write(9,*)
    write(*,'(1x,"Lattice parameter: ",3(1x,f8.4))')  a_par
    write(9,'(1x,"Lattice parameter: ",3(1x,f8.4))')  a_par
    write(*,'(1x,"Atoms:           ",51(1x,a8))')  (at_name_plot(i),i=1,n_atom+1)
    write(9,'(1x,"Atoms:           ",51(1x,a8))')  (at_name_plot(i),i=1,n_atom+1)
    write(*,'(1x,"Atoms no.:         ",51(1x,i8))') (i,i=1,n_atom+1)
    write(9,'(1x,"Atoms no.:         ",51(1x,i8))') (i,i=1,n_atom+1)
    write(*,'(1x,"Occupancy nominal: ",51(1x,f8.4))') (at_occup_r(i),i=1,n_atom),1.                      !1. for TOT
    write(9,'(1x,"Occupancy nominal: ",51(1x,f8.4))') (at_occup_r(i),i=1,n_atom),1.
    write(*,'(1x,"Occupancy real 1:  ",51(1x,f8.4))') (at_occup_1(i),i=1,n_atom),1.                      !1. for TOT
    write(9,'(1x,"Occupancy real 1:  ",51(1x,f8.4))') (at_occup_1(i),i=1,n_atom),1.
    write(*,'(1x,"Occupancy real 2:  ",51(1x,f8.4))') (at_occup_2(i),i=1,n_atom),1.                      !1. for TOT
    write(9,'(1x,"Occupancy real 2:  ",51(1x,f8.4))') (at_occup_2(i),i=1,n_atom),1.
    write(*,'(1x,a," weights:  ",50(1x,f8.4))') trim(at_weight_scheme(j_weight)),(sum(ind_pseudo(:,j)*at_weight),j=1,n_atom+1)
    write(9,'(1x,a," weights:  ",50(1x,f8.4))') trim(at_weight_scheme(j_weight)),(sum(ind_pseudo(:,j)*at_weight),j=1,n_atom+1)
    write(*,'(1x,"Mean rel_error:    ",50(1x,f8.4))') (rdf_err(i),i=1,n_atom),rdf_tot_err						
    write(*,*) 
    write(*,'(1x,"Atom density [Å-3] and masks",e14.6,4x,(50i3))')ro_0,(at_mask(i),i=1,n_atom)        
    write(9,'(1x,"Atom density [Å-3] and masks",e14.6,4x,(50i3))')ro_0,(at_mask(i),i=1,n_atom)        
    write(*,*) 
    
    if(minval(at_mask(1:n_atom))==1) then  !if all masks =1 don't put them into plot title
      masks = ''
    else
      write(masks,'("  Masks: ",50i1.1)') (at_mask(i),i=1,n_atom)			!to be used in plot titles
    endif

    at_weight_sq_av = sum(at_occup_r*(at_scf*at_weight)**2)       !mean squared amplitude per atom
    at_weight_av_sq = sum(at_occup_r*at_scf*at_weight)*sum(at_occup_r*at_scf*at_weight)       !mean amplitude per atom squared

    if(j_mode==5) then
      at_weight_matrix = 1.
      at_mask_matrix = 1.
    else
      do jj=1,n_atom
        do ii=1,n_atom
          at_weight_matrix(ii,jj) = at_weight(ii)*at_weight(jj)
          at_mask_matrix(ii,jj) = at_mask(ii)*at_scf(ii)*at_mask(jj)*at_scf(jj)
        enddo
      enddo
    endif

      if(j_verb==1) then
        write(*,*) 'at_weight',at_weight
        write(*,*) 'at_weight_sq_av,at_weight_av_sq',at_weight_sq_av,at_weight_av_sq
        write(*,*) 'at_weight matrix'
        do ii=1,n_atom
            write(*,*) ii,at_weight_matrix(ii,:) 
        enddo      
        write(*,*) 
      endif
           
! *** generate the smoothing profile (Gauss) and ...
!
    if(j_mode<=2) then
      write(*,'(a,f5.2,a)') 'Gaussian smooth FWHM in steps of',pdf_step,'Å (1 no smoothing, 4 - 20 useful)'
      read(*,*) smooth_fwhm
  
      if(nint(smooth_fwhm)==1) then
        n_smooth=1			               !no smoothing
      else
        n_smooth = nint(5*smooth_fwhm)/2+1
      endif
    else
      n_smooth=1			               !smoothing by FT will come later
    endif
    
    allocate(f_smooth(n_smooth))
    rdf_p2_ws = .0
    rdf_p2_plot = .0

    do j=1,n_smooth
      f_smooth(j) = 2.**(-((j-n_smooth/2-1)/(.5*smooth_fwhm))**2) 
    enddo
    f_smooth = f_smooth/sum(f_smooth(1:n_smooth))				!profile normalized to unit integral

! *** ... apply it to calculate g(r)
    if(n_smooth>1) write(*,*) 'Applying Gaussian smoothing with FWHM=',nint(smooth_fwhm),' steps of',pdf_step,'[A^]'
    write(smooth,'("Smooth FWHM",f5.2," [A]")') pdf_step*smooth_fwhm


! *** calculate g(r) components & total

    if(n_smooth==1.or.j_mode>=3) then
      do i = 1,n_pdf
        rdf_p2_ws(:,:,i) = rdf_p2_n(:,:,i)*at_weight_matrix*at_mask_matrix/at_weight_av_sq         
      enddo               
    else
      rdf_p2_ws(:,:,1) = rdf_p2_n(:,:,1)*at_weight_matrix*at_mask_matrix/at_weight_av_sq     !don't smooth the 1st channel nor take it into smoothing of the others
      do i = 2,n_pdf      
        do j = 1,n_smooth
          if (i-n_smooth/2+j-1.le.2) cycle  !put zeros out of range
          if (i-n_smooth/2+j-1.gt.n_pdf) then
            rdf_p2_ws(:,:,i)=rdf_p2_ws(:,:,i)+f_smooth(j)*rdf_p2_n(:,:,n_pdf)     !replace out-of-range points by the asymptote, placed there at the normalisation
          else
            rdf_p2_ws(:,:,i)=rdf_p2_ws(:,:,i)+rdf_p2_n(:,:,i-n_smooth/2+j-1)*f_smooth(j)
          endif
        enddo 
        rdf_p2_ws(:,:,i) = rdf_p2_ws(:,:,i)*at_mask_matrix*at_weight_matrix/at_weight_av_sq            !now we allow for b_coh in the g(r), it won't change the overall norm, just redistribute intensity between partials
      enddo
    endif
    deallocate(f_smooth)
 
!     arg = ro_0
!     write(*,*) 'Confirm/modify ro_0',arg
!     read(*,*) arg
!     if(arg/=ro_0) then
!       a_par = a_par*(ro_0/arg)**.3333333
!       write(*,*) 'Modified lattice parameter:',a_par
!       ro_0 = arg
!     endif
!     
! *** calculate G(r) components & total
    if(j_mode>=2) then
      do j=1,n_atom
        do i=1,n_atom  
          rdf_p2_ws(i,j,:) = 4.*pi*ro_0*r*(rdf_p2_ws(i,j,:)-at_occup_1(i)*at_occup_2(j)*at_weight_matrix(i,j)*at_mask_matrix(i,j)/at_weight_av_sq)    !for j_mode>=2 make it G(r)  
        enddo
      enddo
    endif

! *** calculate S(Q) components & total
    if(j_mode>=3) then
      write(*,*) 'S(Q) q_step, q_range:',q_step,(n_pdf-1)*q_step
      write(*,*) 'Hann smooth FWHM in Q_steps (<=1. no smoothing)'
      read(*,*) smooth_fwhm
    
      if(smooth_fwhm<=1.) then
        ft_wind = 1.
        smooth = ''
      else
        ft_wind = .0
        do i=1,min(n_pdf,nint(n_pdf/(smooth_fwhm-1.)))
          ft_wind(i) = .5*(1.+cos(pi*(smooth_fwhm-1.)*(i-1)/real(n_pdf-1)))       ! FT window == Hann profile
        enddo
        if(n_smooth>1) write(*,*) 'Applying Gaussian smoothing with FWHM=',smooth_fwhm,' steps of',q_step,'[A-1]'
        write(smooth,'("Smooth FWHM",f5.2," [A-1]")') q_step*smooth_fwhm
      endif
   
      do j=1,n_atom
        do i=1,n_atom 
          rdf_fft(1:n_pdf)= rdf_p2_ws(i,j,1:n_pdf)*ft_wind(1:n_pdf)
          do ii=1,n_pdf-2
            rdf_fft(2*n_pdf-2-ii+1)= -1.*rdf_p2_ws(i,j,ii+1)*ft_wind(ii+1)    !antisymmetric periodic continuation
          enddo
          rdf_fft = fft((1.,0.)*rdf_fft,inv=.false.)*pdf_step/(.0,-2.)          !pdf_step = pdf_range/(n_pdf-1) with 1/sqrt(n_pdf-1) coming once from FT and and once from sampling norm
          rdf_p2_ws(i,j,:) = rdf_fft(1:n_pdf)        
          rdf_p2_ws(i,j,2:n_pdf) =  rdf_p2_ws(i,j,2:n_pdf)/q(2:n_pdf)

          if(j_weight<=2) then
!           if(j_mode==3)   NOTHING TO DO
            if(j_mode==4) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)+at_occup_1(i)*at_occup_2(j)*at_weight_matrix(i,j)*at_mask_matrix(i,j)/at_weight_av_sq                 ! case of structure factor S(Q)
            if(j_mode==5) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)/(at_occup_2(j)*at_occup_1(i))         ! case of Faber-Ziman partials FZ(Q),at_weight_matrix=1     
            if(j_mode>=6) then
              rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)*at_weight_av_sq
              if(i==j) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)+at_occup_1(i)*at_mask(i)*at_weight_matrix(i,i)       ! case of total intensity I(Q)
              if(j_mode==6) rdf_p2_ws(i,j,2:n_pdf) = (rdf_p2_ws(i,j,2:n_pdf))/at_weight_sq_av       ! case of total intensity I(Q)
            endif
          elseif(j_weight==3) then
            rdf_p2_ws(i,j,:) = rdf_p2_ws(i,j,:)/(at_weight_matrix(i,j)/at_weight_av_sq)                   !suppress the effective constant Xray weight
            if(j_mode/=5) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)*x_ffq(i,2:n_pdf)*x_ffq(j,2:n_pdf)  !introduce Q-dependent Xray weight

            if(j_mode==3) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)/x_ffq_av_sq(2:n_pdf)                                                                                                         ! case of interference function F(Q)     
            if(j_mode==4) rdf_p2_ws(i,j,2:n_pdf) = (rdf_p2_ws(i,j,2:n_pdf)+at_occup_1(i)*at_occup_2(j)*at_mask(i)*at_mask(j)*x_ffq(i,2:n_pdf)*x_ffq(j,2:n_pdf))/x_ffq_av_sq(2:n_pdf) 
            if(j_mode==5) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)/(at_occup_2(j)*at_occup_1(i))           ! case of Faber-Ziman partials FZ(Q)- this is to avoid switching to j_weight=1 for the FZ partials 
            if(j_mode>=6.and.i==j) rdf_p2_ws(i,j,2:n_pdf) = (rdf_p2_ws(i,j,2:n_pdf)+at_occup_1(i)*at_mask(i)*x_ffq(i,2:n_pdf)**2)
            if(j_mode==6) rdf_p2_ws(i,j,2:n_pdf) = rdf_p2_ws(i,j,2:n_pdf)/x_ffq_sq_av(2:n_pdf)
          endif   
        enddo
      enddo
    endif

! *** generate rdf_p2_plot and sum up symmetric off-diagonal elements
         
    do i=1,n_pdf
      rdf_p2_plot(:,:,i) = matmul(matmul(transpose(ind_pseudo(:,1:n_atom_tot)),rdf_p2_ws(:,:,i)),ind_pseudo(:,1:n_atom_tot))
    enddo    

    if(j_mode/=5) then                 !sum up the off-diagonals to give the total cross-correlation (apart of the FZ case)
      do j=1,n_atom_tot
        do i=1,j  
          if(i/=j) rdf_p2_plot(i,j,:) = rdf_p2_plot(i,j,:)+rdf_p2_plot(j,i,:)
          if(i/=j) rdf_p2_plot(j,i,:) = rdf_p2_plot(i,j,:)
          if(i/=j.and.i>n_atom.and.j>n_atom) rdf_p2_plot(i,j,:) = 0.
         enddo
        tot_scale = 1.
      enddo
    else
        tot_scale = .0
        part_scale = 1.
    endif
    
    if(j_mode<=2) then
      x => r
    else        
      x => q
    endif

    x_start = x(1)        !origin of X_axis
    i_start = 2           ! 1st plotted point
    if(x_end>pdf_range) x_end = pdf_range
    if(j_mode>=3) then
      i_start = int(.4/q_step)     !for S(Q) avoid FT oscillations close to the origin
      if(x_end>(n_pdf-1)*q_step) x_end = (n_pdf-1)*q_step
    endif
    i_end = (x_end-x(1))/(x(2)-x(1))
    n_plot = i_end-i_start+1
     
! *** open the PGPLOT graphics window (X11)
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
  
! *** Set JK colors for line plots								  
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
    do j=1,n_part
      if(part_scale(j)==-1) then
        part_scale(j) = at_occup_1(ind_part(1,j))*at_occup_2(ind_part(2,j))*at_weight_matrix(ind_part(1,j),ind_part(2,j))/at_weight_av_sq
        if(ind_part(1,j)/=ind_part(2,j)) part_scale(j) = part_scale(j)+at_occup_2(ind_part(1,j))*at_occup_1(ind_part(2,j))*at_weight_matrix(ind_part(2,j),ind_part(1,j))/at_weight_av_sq
        part_scale(j) = 1./part_scale(j)
      endif
    enddo

    if(c_max==.0) then          !only do this at the real start of plot_loop or when switching to partials display
      c_min = minval(rdf_p2_plot(n_atom+1,n_atom+1,i_start+10:i_end))        !avoid maximum close to the origin
      c_max = maxval(rdf_p2_plot(n_atom+1,n_atom+1,i_start+10:i_end))
      c_max = .1*(int(10*c_max)+1)
    endif
          
    if(j_mode==5) then          
      c_min = minval(rdf_p2_plot(1:n_atom,1:n_atom,i_start+10:i_end))        !avoid maximum close to the origin
      c_min = .1*(int(10*c_min)-1)
      c_max = maxval(rdf_p2_plot(1:n_atom,1:n_atom,i_start+10:i_end))
      c_max = .1*(int(10*c_max)+1)
    endif
          
    write(*,*) 
    write(*,*) 'Vertical scale c_min,c_max',c_min,c_max
    
    write(plot_header,'(a,"    ",a," weights  ")') trim(file_dat_t0),trim(at_weight_scheme(j_weight))
    plot_header = trim(subst_name)//'  '//trim(pdf_out(j_mode))//'  '//trim(plot_header)//trim(masks)//'  '//trim(smooth)

    scale_loop: do

  call date_and_time(c_date,c_time,c_zone,i_time)
  write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
       
      CALL PGSLCT(j_xserv)
      CALL PGSCI (1)  !white
      CALL PGSCH(1.)
      CALL PGSLW(2)
      CALL PGSLS (1)  !full
      CALL PGENV(x_start,x_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
      CALL PGLAB(x_label(j_mode),y_label(j_mode),trim(plot_header))  !put the axis labels
      CALL PGSCI (1)  !white
      CALL PGSLW(5)			!operates in steps of 5

      x_plot = x_start+.75*(x_end-x_start)
      y_plot = c_max-.1*(c_max-c_min)

      if(j_ext==1) then
        CALL PGTEXT (x_plot,y_plot,trim(file_inp))
        do j=1,n_part_ext
          if(ext_scale(j)==0) cycle
          CALL PGSCI(15)  !grey
          CALL PGSLW(10)			!operates in steps of 5
          CALL PGLINE(n_x,x_ext(1:n_x),ext_scale(j)*(y_ext(1:n_x,j)+ext_dy(j)))  !plots the curve 
          CALL PGSTBG(0)																				 !erase graphics under text
          CALL PGSLW(5)			!operates in steps of 5
          y_plot = y_plot-.04*(c_max-c_min)
          write(plot_title,'(a,a,i1.1," x",2f7.3)') trim(pdf_out(j_mode)),'_ext',j,ext_scale(j),ext_dy(j)
          CALL PGTEXT (x_plot,y_plot,plot_title)
        enddo
      endif

      if(tot_scale/=0.and.j_mode/=5) then
        y_plot = y_plot-.06*(c_max-c_min)
        CALL PGSCI(1)  !white
        CALL PGSLW(5)			!operates in steps of 5
        CALL PGLINE(n_plot,x(i_start:i_end),tot_scale*rdf_p2_plot(n_atom+1,n_atom+1,i_start:i_end))  !plots the curve         
        CALL PGSTBG(0)																				 !erase graphics under text
        write(plot_title,'(a,a," x",f7.3)') 'Total','    ',tot_scale
        CALL PGTEXT (x_plot,y_plot,plot_title)
      endif

     do j=1,n_part
        if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
        if(part_scale(j)==0) cycle
        write(plot_title,'(a,a," x",f7.3)') at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
        y_plot = y_plot-.04*(c_max-c_min)
        CALL PGSCI (j+1)  !basic colors
        if (j_mode==5.or.tot_scale==0)then
          CALL PGSLW(5)
        else
          CALL PGSLW(2)
        endif
        CALL PGLINE(n_plot,x(i_start:i_end),part_scale(j)*rdf_p2_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
        CALL PGSTBG(0)																				 !erase graphics under text
        CALL PGSLW(5)			!operates in steps of 5
        CALL PGTEXT (x_plot,y_plot,plot_title)
      enddo

      x_plot = x_start+.8*(x_end-x_start)
      y_plot = c_min-.1*(c_max-c_min)
      CALL PGSCI (1)  !white needs to be reset after PGLAB
      CALL PGSTBG(0)																				 !erase graphics under text
      CALL PGSLW(2)			!operates in steps of 5
      CALL PGSCH(.6)
      CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp))
     
      write(*,*) 'Adjust vertical scale (min, max) (0 0 EXIT, -1 -1 to adjust plot range)'
      c1 = c_min
      c2 = c_max
      read(*,*) c1,c2
      if(c1==.0.and.c2==.0) then
        exit
      elseif(c1==-1.and.c2==-1) then              
        write(*,*) 'Confirm/adjust plot range, max =',(n_pdf-1)*(x(2)-x(1))
        write(*,'("x_start, x_end ",2f7.1,":  ")',advance='no') x_start,x_end 
        read(*,*) x_start,x_end 
        if(x_start<=.0) x_start = x(1)
        if(x_end>(n_pdf-1)*(x(2)-x(1))) x_end = (n_pdf-1)*(x(2)-x(1))
        i_start = (x_start-x(1))/(x(2)-x(1)) + 1
        i_end = nint((x_end-x(1))/(x(2)-x(1)))
        n_plot = i_end-i_start+1
      else
        c_min = c1
        c_max = c2
      endif
    enddo scale_loop
         
! **** Prepare and plot the same on .PS, look for existing output files in order not overwrite them				

    if(j_ps==1) then
      jfile = 1
      do						!look for existing .ps files to continue numbering
        if(t_single)then
          write(file_ps,'(a,"_pdf","_",i2.2,a)') trim(file_master),jfile,trim(pg_ext)
          write(file_res,'(a,"_pdf","_",i2.2,".txt")') trim(file_master),jfile
        else
          write(c_jfile,'("_",i2.2)') jfile
          if(nfile_min<=9999) then
            write(c_nfile_min,'(i4.4)') nfile_min
          elseif(nfile_min>=10000) then
            write(c_nfile_min,'(i8)') nfile_min
          endif
          c_nfile_min = '_'//adjustl(c_nfile_min)

          if(nfile>nfile_min) then
            if(nfile<=9999) then
              write(c_nfile,'(i4.4)') nfile
            elseif(nfile>=10000) then
              write(c_nfile,'(i8)') nfile
            endif
            c_nfile = '_'//adjustl(c_nfile)
          else
            c_nfile = ''
          endif

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

!       if(j_grid==1) then
!         CALL PGSCI (1)  
!         CALL PGSCH(1.)
!         CALL PGSLW(2)
!         CALL PGSLS (1)  !full
!         XOPT = 'ABCGN'				!draw grid lines - after the color map, otherwise obscured
!         YOPT = 'ABCGN'
!         XTICK = .0
!         YTICK = .0
!         NXSUB = 5
!         NYSUB = 5
!         CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
!       endif

      plot_header = trim(subst_name)//'  '//file_ps//'  '//trim(at_weight_scheme(j_weight))//' weights  '//trim(masks)//'  '//trim(smooth)	

      CALL PGSCI (1)  !white
      CALL PGSCH(1.)
      CALL PGSLW(2)
      CALL PGSLS (1)  !full
      CALL PGENV(x_start,x_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes w/o grid
      CALL PGLAB(trim(x_label(j_mode)), trim(y_label(j_mode)),trim(plot_header))  !put the axis labels
      CALL PGSCI (1)  !white needs to be reset after PGLAB
      CALL PGSLW(5)			!operates in steps of 5
      x_plot = x_start+.75*(x_end-x_start)
      y_plot = c_max-.1*(c_max-c_min)

      if(j_ext==1) then
        CALL PGTEXT (x_plot,y_plot,trim(file_inp))
        do j=1,n_part_ext
          CALL PGSCI(15)  !grey
          CALL PGSLW(10)			!operates in steps of 5
          CALL PGLINE(n_x,x_ext(1:n_x),ext_scale(j)*(y_ext(1:n_x,j)+ext_dy(j)))  !plots the curve 
          write(plot_title,'(a,a,i1.1," x",2f7.3)') trim(pdf_out(j_mode)),'_ext',j,ext_scale(j),ext_dy(j)
          y_plot = y_plot-.04*(c_max-c_min)
          CALL PGSTBG(0)																				 !erase graphics under text
          CALL PGSLW(5)			!operates in steps of 5
          CALL PGTEXT (x_plot,y_plot,plot_title)
        enddo
      endif

      if(tot_scale/=0.and.j_mode/=5) then
        y_plot = y_plot-.06*(c_max-c_min)
        CALL PGSCI(1)  !white
        CALL PGSLW(5)			!operates in steps of 5
        CALL PGLINE(n_plot,x(i_start:i_end),tot_scale*rdf_p2_plot(n_atom+1,n_atom+1,i_start:i_end))  !plots the curve         
        CALL PGSTBG(0)																				 !erase graphics under text
        write(plot_title,'(a,a," x",f7.3)') 'Total','    ',tot_scale
        CALL PGTEXT (x_plot,y_plot,plot_title)
      endif

      do j=1,n_part
        if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
        if(part_scale(j)==0) cycle
        CALL PGSCI (j+1)  !red-green-blue
        CALL PGSLW(2)
        CALL PGLINE(n_plot,x(i_start:i_end),part_scale(j)*rdf_p2_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
        write(plot_title,'(a,a," x",f7.3)') at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
        y_plot = y_plot-.04*(c_max-c_min)
        CALL PGSTBG(0)																				 !erase graphics under text
        CALL PGSLW(5)			!operates in steps of 5
        CALL PGTEXT (x_plot,y_plot,plot_title)				!the x_plot,y_plot are in world (axis units) coordinates
        CALL PGSLW(2)
      enddo    

! *** print footer with program version & date_and_time
      x_plot = x_start+.8*(x_end-x_start)
      y_plot = c_min-.1*(c_max-c_min)
      CALL PGSCI (1)  !white needs to be reset after PGLAB
      CALL PGSTBG(0)																				 !erase graphics under text
      CALL PGSLW(2)			!operates in steps of 5
      CALL PGSCH(.6)
      CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp))

      CALL PGCLOS
      write(*,*) ' Postscript output written to: ',file_ps	
      write(9,*)
      if(j_acc==2) then
        write(9,*) 'Integration ',trim(int_mode),'  ',n_int,'cell pairs, j_rand',j_rand
      else
        write(9,*) 'Integration ',trim(int_mode),'  ',n_int,'cell pairs'
      endif
      write(9,*) '  ','Masks:',(at_mask(i),i=1,n_atom)				
      if(n_smooth>1) write(9,*) '  ','Smoothing FWHM:',nint(smooth_fwhm)	
      write(9,*) ' Postscript output written to: ',file_ps

! *** save the PDF results into an ASCII file (each line corresponds to a distance point)
!     look for existing output files in order not overwrite them				
      if(j_txt==1) then
        open (4,file=file_res)

! *** write PDF file header
        write(4,*) 'Substance:',trim(subst_name),'       ',time_stamp
        write(4,*) 'Input files: ',trim(file_dat_t0),' to ',trim(file_dat),' step',nfile_step									   
        write(4,*) 'OMP processes  = ', proc_num_in									
        write(4,*) 'Integration ',trim(int_mode),'  ',n_int,'cell pairs, j_rand',j_rand
        write(4,*) 'Supercell size:',n_row																				
        write(4,*) 'Unit cell parameter:',a_par																				
        write(4,*) 'Atoms/unit_cell:',n_atom
        write(4,*) 	 'Atoms  :',(('    '//at_name_plot(i)),i=1,n_atom)
        write(4,'(1x,"Atom no. ",50i8)') (i,i=1,n_atom)
        write(4,'(1x,"Occups: ",50f8.4)') (at_occup_r(i),i=1,n_atom)
        write(4,'(1x,"Weights: ",a,5x,50f8.4)') at_weight_scheme(j_weight),(at_weight(i),i=1,n_atom)
        write(4,'(1x,"Masks	: ",50i8)') (at_mask(i),i=1,n_atom)
        if(j_mode<=2) write(4,*)  'Smoothing FWHM [Å]',smooth_fwhm*pdf_step
        if(j_mode>=3) write(4,*)  'Smoothing FWHM [Å-1]',smooth_fwhm*q_step
        write(4,*)
        write(4,*) 'Output: ',pdf_out(j_mode)

! *** write the PDFs
        if(j_out==1) then
          if(j_mode==5) then       
            write(4,'(1x,"Partials index",32(2x,5(2x,2i3)))') ((ii,jj,ii=1,jj),jj=1,n_atom)
            write(4,'(3x,a,1x,32(2x,5(a,"_",a,3x)))') trim(x_label(j_mode)),((trim(at_name_plot(ii)),trim(at_name_plot(jj)),ii=1,jj),jj=1,n_atom)
            do i=1,n_pdf
                write(4,'(1x,f7.3,32(2x,5f8.3))') x(i),((rdf_p2_plot(ii,jj,i),ii=1,jj),jj=1,n_atom)   !TOT doesn't make sense in the PARTIALS output
            enddo
          else
            write(4,'(1x,"Partials index  ",32(2x,5(2x,2i3)))') ((ii,jj,ii=1,jj),jj=1,n_atom)
            write(4,'(3x,a,"   Total     ",32(2x,5(a,"_",a,3x)))') trim(x_label(j_mode)),((trim(at_name_plot(ii)),trim(at_name_plot(jj)),ii=1,jj),jj=1,n_atom)
            do i=1,n_pdf
              write(4,'(1x,f7.3,1x,f8.3,32(2x,5f8.3))') x(i),rdf_p2_plot(n_atom+1,n_atom+1,i),((rdf_p2_plot(ii,jj,i),ii=1,jj),jj=1,n_atom)
            enddo
          endif
        else
          write(4,'(1x,"Partials index  ",32(2x,5(2x,2i3)))') (ind_part(1,j),ind_part(2,j),j=1,n_part)
          write(4,'(3x,a,"   Total     ",32(2x,5(a,"_",a,3x)))') trim(x_label(j_mode)),(trim(at_name_plot(ind_part(1,j))),trim(at_name_plot(ind_part(2,j))),j=1,n_part)
          do i=i_start,i_end
            write(4,'(1x,f7.3,1x,f8.3,32(2x,5f8.3))') x(i),rdf_p2_plot(n_atom+1,n_atom+1,i),(rdf_p2_plot(ind_part(1,j),ind_part(2,j),i),j=1,n_part)
          enddo
        endif
        close(4)
        write(*,*) ' Text output written to: ',file_res	  
        write(9,*) ' Text output written to: ',file_res	  

      endif 																							!j_txt
    endif 																							!j_ps

! ***   all done, now decide what comes next
  
   way_point: do
      write(*,*) 'Choose output options (MODE is ',trim(pdf_out(j_mode)),', FILE output is ',trim(ps_out(j_ps+1)),', SIZE is ',trim(size_out(j_out+1)),'):'
      write(*,*) '       1   REPLOT the last graph '
      write(*,'("        2   select max ",i2," partial PDFs & replot")') n_part_max
      write(*,*) '       3   adjust partials scales & replot '
      write(*,*) '       4   modify atom WEIGHTS (',trim(at_weight_scheme(j_weight)),')'
      write(*,*) '       5   edit atom MASKS ',trim(masks(8:))
      write(*,'("        6   create/modify max ",i2," PSEUDO_ATOMS ")') n_pseudo_max
      write(*,*) '       7   toggle FILE output ',trim(ps_out(mod(j_ps+1,2)+1)),' (mind the J_TXT switch in .PAR)'
      if(j_ext==0)then
        write(*,*) '       8   import external ',trim(pdf_out(j_mode)),' curve '
      else
        write(*,*) '       8   change/close external ',trim(pdf_out(j_mode)),' curve ',file_inp
      endif

!!!					write(*,*) '       8   toggle .TXT output SIZE to ',trim(size_out(mod(j_out+1,2)+1))
      write(*,*) '       9   RESTART with updated weights, masks & pseudo_atoms '
      write(*,*)        
      write(*,*) '       10  select the PDF MODE, actual: ',trim(pdf_out(j_mode))
      
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
          write(*,"('Confirm/modify total scale factor ',4f4.1)") tot_scale
          read(*,*) tot_scale
          write(*,'("Confirm/modify partial scale factors ",4f4.1,"(-1 unit asymptote)")') part_scale(1:n_part)
          read(*,*) part_scale(1:n_part)
          if(j_ext==1) then
            write(*,'("Confirm/modify external data scale factor",4f6.3)') ext_scale(1:n_part_ext)
            read(*,*) ext_scale(1:n_part_ext)
            write(*,'("Confirm/modify external data y-shift ",4f6.3)') ext_dy(1:n_part_ext)
            read(*,*) ext_dy(1:n_part_ext)
          endif
          cycle plot_loop

        case(4) 
          do
            write(*,*) 'Atom weights (1= unit, 2= neutron b_c, 3= Xray f(Q))'	
            read(*,*) j_weight
            if(j_weight<=3) exit
          enddo
          
          if(j_weight==1) then
            at_weight(1:n_atom) = 1.
          elseif(j_weight==2) then
            at_weight(1:n_atom) = b_coh(1:n_atom)
          elseif(j_weight==3) then
            at_weight(1:n_atom) = x_ffq(1:n_atom,j_xff)
          endif
          cycle way_point
              
        case(5) 
          write(*,'(" Actual masks:",(50i3))') (at_mask(i),i=1,n_atom)
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
          if(j_ext/=0) deallocate(x_ext,y_ext)
          write(*,*) 'Confirm or modify input file name ("=" close file):'
          read(*,*) file_inp
          if(index(file_inp,'=')/=0) then
            j_ext = 0
            cycle way_point
          else
            open(4,file=trim(file_inp),iostat=ios)
            if(ios.ne.0) write(*,*) 'File ',trim(file_inp),' not opened! IOS =',ios
            write(*,*) 'number of lines to skip, to read, column X, column Y(1) ... Y(4)? (pad by 0 for not used)'
            read(*,*) n_ext_skip,n_x,j_x,ind_ext

            n_part_ext = 0
            do j=1,n_part_max
              if(ind_ext(j)>0) n_part_ext = n_part_ext+1
            enddo
            if(n_part_ext==0) then
              write(*,*) 'Number of external data columns must be >0'
              cycle way_point
            endif
          
            allocate(x_ext(n_x),y_ext(n_x,n_part_ext))
            ext_scale = 1.
            if(n_ext_skip>=1) then
              do i=1,n_ext_skip
                read(4,*)
              enddo
            endif
        
            read(4,'(a)') line
            do i=1,128
              read(line,*,iostat=ios) data_in(1:i)
              if (ios/=0) exit
            enddo
            n_line = i-1
          
            do i=1,n_x
              read(4,*,iostat=ios) data_in(1:n_line)
              if(i==1) write(*,*) data_in(1:n_line)
              if(ios==-1) then
                n_x = i-1
                exit
              endif
              x_ext(i) = data_in(j_x)
              do j=1,n_part_ext
                y_ext(i,j) = data_in(ind_ext(j))
              enddo
            enddo
            close(4)
            write(*,*) n_x,'data points read in'
            j_ext = 1
           endif 
           cycle way_point

        case(9) 
          exit plot_loop

        case(10) 
          do
            write(*,*) 'Select the output MODE: ',(j,'  ',trim(pdf_out(j)),j=1,n_mode),'_unscaled'
            read(*,*) j_mode         
            if(j_mode>=1.and.j_mode<=n_mode) exit plot_loop
          enddo

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
      
! **** minimal random number generator (Numerical Recipees))
!           
!     Minimal” random number generator of Park and Miller. 
!     Returns a uniform random deviate between 0.0 and 1.0. 
!     Set or reset idum to any integer value (except the unlikely value MASK) to initialize the sequence; 
!     idum must not be altered between calls for successive deviates in a sequence.
!
FUNCTION rand(idum)
  INTEGER(4) idum,IA,IM,IQ,IR,MASK
  PARAMETER (IA=16807,IM=2147483647,IQ=127773,IR=2836,MASK=123459876)
  INTEGER(4) k,rand
!       write(*,*) 'R: idum0',idum
  idum=ieor(idum,MASK)
!      write(*,*) 'R: idum1',idum
  k=idum/IQ
!     write(*,*) 'R: k',k
  idum=IA*(idum-k*IQ)-IR*k
!     write(*,*) 'R: idum2',idum
  rand = idum-IM/2
!     write(*,*) 'R: rand',rand
  idum=ieor(idum,MASK)
!     write(*,*) 'R: idum2',idum
!     write(*,*) 
  return
END      
      
! **********************************************************************************************************
! **** string conversion to all upper case
!     
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

! **** string conversion to all lower case
!          
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
     
