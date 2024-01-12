
program mp_sql56

! *************************************************************************************
! *****
! *****  %%%%%%%%%%%%%%%%   		 program MP_SQL 1.56   		 %%%%%%%%%%%%%%%%%%%%%%
! *****
! *****   calculates and plots liquid scattering functions S(Q) from simulated data
! *****
!**********					Copyright (C) 2022  Jiri Kulda, Grenoble/Prague          **********
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
! *****  %%%%%%%%%%%%%%%%       program MP_SQL 1.55      %%%%%%%%%%%%%%%%%%%%%%%%
! ***** 
! ***** Ver 1.55  - forked from MP_SQOM55 on 01/08/2023
! *****           - uses spherical average of 3D NUFFT amplitudes/intensities
! ***** 
! ***** Ver 1.56  - optimised generation of diverse scattering functions
! *****           - true 3D_NUFFT replaced by the corresponding sequence of 2D_NUFFT slices
! *****           
! ***** 
! ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
! ***** record numbers are directly related to cell positions and atom types
! *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
! *****						1 going for the header record
! *****
! ***** atom positions are converted to and recorded in reduced lattice coordinates (x) 
! ***** 
  use omp_lib
			
  integer,parameter :: l_rec  =  1024		    !record length in real(4)
  integer,parameter :: n_mc  =  1e6		      !unit of MC events count
  real(8), parameter   :: pi=3.141592653589793238462643383279502884197
  real(8), parameter   :: twopi = 2.d0*pi
  real(4), parameter   :: twopi_s = 2.*pi
  real(4), parameter   :: pi_s = pi
  real(4), parameter   :: sqrt2 = sqrt(2.)
  real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
  real, parameter   :: h_bar = 6.350776 !DAPS*ps/2Pi Planck const
  real, parameter   :: h_daps = 39.9031	!DAPS*ps=DAPS/THz  
  real, parameter   :: n_mass = 1.008665	!atomic units  
  
  logical        :: found_txt,found_ps,found,t_single,nml_in,inv,read_ind,non_zero,inside
  character(1)   :: xyz(3)=(/'X','Y','Z'/)
  character(4)   :: atom,ps_out(2),version_t,head
  character(5)   :: pg_ext,c_dir(5)=(/'[00X]','[0X0]','[0XX]','[-XX]','[XXX]'/)
  character(10)	 :: prompt,space = '          '
  character(10)  :: string,section,c_date,c_time,c_zone,mode,ext,pg_out,c_nfile_min,c_nfile,c_jfile
  character(40)  :: subst_name,file_title,file_master,time_stamp,x_title,y_title,masks,at_weight_scheme(4),int_mode,mp_tool
  character(40)  :: file_dat,file_inp,file_dat_t0,file_res,file_ps,file_log,x_file_name,wedge_label,string_in,smooth
  character(256) :: line,plot_title1,plot_title2,plot_title,plot_header,scan_title,interp_title,data_path,cwd_path,rec_str
  character(l_rec):: header_record

  character(4),allocatable   :: at_label(:),at_name_par(:),at_name_ext(:),at_name_plot(:),at_name_pseudo(:),pdf_out(:)
  character(16),allocatable ::	curve_label(:),x_label(:),y_label(:)
  
  integer, allocatable :: i_site(:)		!index attributing site number to atom number
  integer,allocatable ::  ind_at(:),map_unit(:),at_mask(:),ind_part(:,:),ind_ext(:),numbers(:),ind_mo(:)              !,ind_pseudo(:,:),m_pseudo(:)

  real, allocatable :: at_base(:,:),at_weight(:),at_weight_matrix(:,:),at_mask_matrix(:,:),at_av_matrix(:,:),at_weight_av_sq,at_weight_sq_av,b_coh(:),x_ffpar(:,:),x_ffq(:,:),q_norm(:,:,:)    !,c_pseudo(:),c_pseudo_mean(:)
  real, allocatable :: sq_hist_tot(:),sq_plot_tot(:),sq_hist(:,:,:),sq_hist_k(:,:,:,:),sq_plot(:,:,:),x_ext(:),y_ext(:,:),ext_scale(:),ext_dy(:)			!s_base(:,:),q(:,:,:),at_vel_file(:,:)		!atom basis, site basis fractional coordinates
  real, allocatable :: csp(:),qq(:),ff(:),t_wind(:),f_smooth(:),at_scf(:) 								!q_at_pos(:),q_at_vel(:),at_pos_perp(:),rq_fft(:)
        
  real(8), allocatable :: xx(:,:)			
  real(8) ::     eps_fft
  integer(8) ::  n_tot,n_qq8(3),n_fft,n_mem
  integer(8), allocatable :: n_sup(:)		
!   integer(8), allocatable :: nul_opt		!nul_opt (since unallocated) used to pass a null pointer to FINUFFT...
  integer(8), pointer :: nul_opt => null()		!nul_opt (since unallocated) used to pass a null pointer to FINUFFT...
  complex(8), allocatable :: cf(:),ampl_tot(:),ampl_tot_3d(:,:,:)												! this is equivalent to complex*16

  complex, allocatable :: om_phase(:),ampl_atom_3d(:,:,:,:),ampl_q(:)

  integer :: j_head_in,ind_rec,hist_ind(3),m_dom1,m_dom2,j_verb,jm,j1m,j2m,j_bc,j_exit
  integer :: i,ii,iii,j,jj,jjj,i1,i2,i3,j1,j2,k,kk,kkk,l,ll,ios,ier,iflag,j_at,i_rec,nrec,nsuper,nlayer,n_r1,n_r2,n_row_save(3)
  integer :: j_basis,j_centred,j_test,nfile,nfile_min,nfile_max,ifile,jfile,nfile_step
  integer :: n_site,i_shift,i_time(8),j_weight,j_wind,j_qsq,n_smooth,n_smooth_fwhm
  integer :: n_qx,n_qy,n_qz,n_q1x,n_q2x,n_q1y,n_q2y,nq_tot,n_q1,n_q2,n_qxg,n_qyg,n_qxp,n_qyp,n_qq(3),n_hist,j_acc,j_rand,n_h,n_h_max,n_skip1,n_skip2

  integer :: n_plot,i_step,i_centre(3),n_ext_skip,n_x,j_x,j_y,j_ext,n_line
  integer :: sc_c1,sc_c2,sc_m,j_proc,proc_num,proc_num_in,thread_num,thread_num_max,j_mode,n_mode
  integer :: ind,j_site,j_q,n_int,n_frame,j_mask,j_logsc,j_grid,j_txt,j_ps,j_out,i_start,i_end,i_hist,n_pseudo_max,n_part,n_part_max,n_part_ext
  integer :: j_atom,n_atom,at_no,n_head,n_pseudo
  integer :: rand1_seed(8),rand2_seed(8),seed_size,d_ind(3)

  real :: at_mass,at_charge,bc_in,s_trig,q_sq,rn,bz_n,bz_nx,bz_ny,bz_nz,eps_x,rand1,rand2,rnd(5),rand0(5)
  real :: at_displ,at_pos(3),at_veloc(3),at_pos_min(3),at_pos_max(3),a_cell_par(3),a_cell(3,3),a_cell_inv(3,3),a_cell_1(3,3),pdf_pix,pdf_pix_shift(3)
  real :: t2,t3,t4,t_dump,t_step,t_step_in,t_tot,tt0,tt,f_width,t_width,temp_par
  real :: arg,q1_x,q1_y,q2_x,q2_y,qx_min,qx_max,qy_min,qy_max,q_min,q_max,qx_ming,qx_maxg,qq_plot,dq(3),dq_p(3),q_step,q_range,pdf_step,pdf_range,q_pos(3),q_pos_norm,q_xff
  real :: d_q(3),x_start,x_end,x_plot,y_plot
  real :: e1(3),e2(3),e3(3),e1_norm,e2_norm,e3_norm,e1_r(3),e2_r(3),e3_r(3),e1v(3),e2v(3),e3v(3),e1_r_norm,e2_r_norm,e3_r_norm,ep_angle,tpq_center(3),tpe1,tpe2,tpev,tpe1_shift,tpe2_shift
  real :: at_volume,cell_volume,cell_angle,angle_r(3),g_matrix(3,3),g_r_matrix(3,3)
  real :: q1(3),q2(3),q1_3d(3),q2_3d(3),q_v(3),q_nx,q_ny,wtime,data_in(128)
  real :: c1,c2,c_min,c_max,c_min_save,c_max_save,sc_r,part_scale(4),tot_scale	
  
! **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
!
  character(4)             :: version
  character(4),allocatable :: at_name_out(:)
  integer(4),allocatable ::   nsuper_r(:)
  integer(4),allocatable,target   :: at_ind_in(:)
  integer(4),pointer :: at_ind(:,:,:)
  
  real(4),allocatable ::	at_occup_r(:),cs_plot(:,:),cs_out(:,:),cs_pgplot(:,:)
  real(4),allocatable,target ::	at_pos_in(:,:),q_hist(:)
  real(4), pointer :: at_pos_file(:,:,:),x(:)
  real(4), pointer :: at_pos_file_bulk(:,:)

  character(16)  :: sim_type,dat_type,input_method,file_par,dat_source,string16,filter_name
  integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,n_tot_in,idum,j_pgc
  real(4)        :: rec_zero(l_rec),filter_fwhm,t_ms,t0,t1,a_par_pdf(3),a_par(3),angle(3),temp

! *** PGPLOT stuff
  integer :: PGOPEN,NC,plot_unit
  real :: TR(6),CONTRA,BRIGHT,p_size
  CHARACTER*4 XOPT, YOPT
  REAL XTICK, YTICK
  INTEGER NXSUB, NYSUB,j_xserv					

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
  namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
  namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
 
  namelist /mp_gen/ j_verb,j_proc       
  namelist /mp_out/ j_weight,j_logsc,j_txt,p_size,j_grid,pg_out,j_ps,j_out,j_pgc        
                      !general rule: namelists of tools should only contain their localisation parameters
                      !what is data-related should pass into data_header
  namelist /mp_pdf/ pdf_range,pdf_step,q_step,x_end,a_par_pdf,pdf_pix,pdf_pix_shift,j_rand,n_h,j_acc,j_mode,n_part_max,n_pseudo_max,n_cond,q_xff  
!
! ********************* Initialization *******************************      
  version_t = '1.56'
  prompt = 'MP_SQL>   '
  mp_tool = 'MP_SQL '//version_t

  print *,'*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  print *

! ********************* Get a time stamp and open a .LOG file *******************************

  call getcwd(cwd_path)
  call date_and_time(c_date,c_time,c_zone,i_time)

  write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
  write(file_log,'("mp_tools_",a,".log")') trim(c_date)
  inquire(file=file_log,exist=found_txt)
  if(found_txt)then
    open (9,file=file_log,position='append')
  else
    open (9,file=file_log,status ='new')
  endif

  write(9,*) 
  write(9,*) 
  write(9,*) trim(time_stamp),'   ',trim(mp_tool),'   ',trim(cwd_path)
  write(9,*) 

! *** other initialisations
  n_mode = 7
  allocate(pdf_out(n_mode),curve_label(n_mode),x_label(n_mode),y_label(n_mode))
!   curve_label = (/'g_tot','G_tot','F_tot','S_tot','I_tot','Z_tot'/)                                           !g(r), G(r) not used actually
  pdf_out = (/'g(r)','G(r)','F(Q)','S(Q)','Z(Q)','I(Q)','I(Q)'/)			
  y_label = pdf_out
  y_label(7) = 'I(Q) [unscaled]'
  x_label = (/'r [A]  ','r [A]  ','Q [A-1]','Q [A-1]','Q [A-1]','Q [A-1]','Q [A-1]'/)
  ps_out = (/'OFF ','ON  '/)			!PGPLOT
  plot_unit = 0
  p_size = 7.	                  !PGPLOT - screen plot size
  j_pgc = 6 										!PGPLOT - JK printer-friendly rainbow for colormaps (no black)
  eps_fft = 1.e-6						!finufft
  iflag = -1						!finufft			
  j_head_in = 0		! if header found will become 1
  at_weight_scheme(1) = 'Unit weights '
  at_weight_scheme(2) = 'Neutron weights'
  at_weight_scheme(3) = 'Xray weights'
  at_weight_scheme(4) = 'FZ partials '
  j_weight = 1
  j_wind = 0			!no FT window in space by default
  nfile_step = 1
  j_verb = 1
  j_test = 0
  s_trig = 0
  j_grid = 0
  j_xserv = 0
  j_ext = 0
  n_part_max = 4
  n_part_ext = 0
  q_step = .02  
  j_acc = 3     
  j_rand = 1
  n_h = 0

  
! *** Generate data file access
  print *,prompt, 'Data file_master & file numbers (min, max; 0 0 single file): '
  read(*,*) file_master,nfile_min,nfile_max
  nfile_step = 1

  t_single = nfile_max==0.and.nfile_min==0

  if(t_single.and.nfile_min==0)then
    nfile_min = 1		!for conformity with loop control conventions
    nfile_max = 1
    nfile_step = 1
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
    print *,space, 'File ',trim(file_dat_t0),' not found! Stop execution.'
    stop
  endif
  
! *** Read the header record 

  i_rec = 1   			
  read(1,rec=i_rec) string_in
    if(string_in(1:8)=='MP_TOOLS') then      !new structure with namelist
      nml_in = .true.
      read(1,rec=i_rec) dat_source,version,string16
      read(string16,*) n_head
    print *,space,		'Input data:  ',dat_source,version,n_head
    i_rec = i_rec+1   							
    read(1,rec=i_rec) header_record
    read(header_record,nml=data_header_1)	
    t0 = t_dump
    t_step_in = t_step
  elseif(head.eq.'TIME'.or.head.eq.'STAT') then                                  !old w/o structure
    print *,space,		'Input data:  ','old header format'
    nml_in = .false.
    read(1,rec=1) sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,n_cond,n_rec					
    n_head = 1
    call up_case(sim_type)
    input_method = 'CELL'
    if(index(sim_type,'BULK')/=0) input_method = 'BULK'
  else
    print *,space, 'header record wrong or missing'
    print *,space, trim(header_record)
    stop
  endif 

  allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom),n_sup(n_atom))
  allocate(at_label(n_atom),at_name_par(n_atom),at_name_plot(n_atom),at_name_ext(n_atom),at_base(n_atom,3),at_mask(n_atom))
  allocate(b_coh(n_atom),x_ffpar(n_atom,9),SOURCE=.0)												!we need this to read .par
  at_mask = 1

  if(nml_in) then      !new structure with namelist
    i_rec = i_rec+1   							
    read(1,rec=i_rec) header_record
    read(header_record,nml=data_header_2)
    if(n_head==4) then
      i_rec = i_rec+1   							
      read(1,rec=i_rec) header_record
      read(header_record,nml=data_header_3)
    endif			
  else                                  !old w/o structure
    read(1,rec=1) sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell_out,n_cond,n_rec,n_tot,at_name_out,at_occup_r,nsuper_r																		!char(16): sim_type,file_par, char(4): at_name
    if(head.eq.'TIME') sim_type = 'TIMESTEP'
    if (head.eq.'STAT') sim_type = 'STATIC'
  endif 
  close(1)
 
! *** cell and box geometry
  nlayer = n_row(1)*n_row(2) 
  nsuper = product(n_row) 
  n_sup(:) = nsuper_r(:)					!n_sup is a 64bit copy needed for FINUFFT, while storage is 32bit for space economy
  cell_volume = product(a_par)		!normalisation constant
  at_volume = product(n_row)*product(a_par)/n_tot
  at_occup_r = at_occup_r/(sum(at_occup_r))


  if(n_cond==0) then
    print *,space, 'Non-periodic boundary conditions!'
    print *,space, 'FT will use Hann window profile - results in 2x lower resolution in Q!'
    j_wind = 1
  endif
 
! **** Read the auxiliary file <file_title.par> with structure parameters, atom names and further info
  print *,prompt, 'Parameter file name (.par to be added) (confirm or type other name): ', file_par
  file_inp = trim(file_par)//'.par'

  open(4,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,space, 'File ',trim(file_inp),' not found! Stop execution.'
    stop
  endif

  write(9,*) 'Read the parameter file:  ',trim(file_inp)

  read(4,nml=mp_gen,iostat=ios)
    if(ios/=0) then
      print *,space, 'Error reading NML = mp_gen from ', trim(file_inp),', check that you have MP_TOOLS at least ',version
      stop
    endif
  rewind(4)

  read(4,nml=mp_out,iostat=ios)
    if(ios/=0) then
      print *,space, 'Error reading NML = mp_out from ', trim(file_inp),', check that you have MP_TOOLS at least ',version
      stop
    endif
    if(j_weight<1.or.j_weight>3) j_weight = 1
  
  call down_case(pg_out) 
  if(index(pg_out,'png')/=0) then
    pg_ext = '.png'
  else
    pg_ext = '.ps'
  endif


  rewind(4)
!   read(4,nml=mp_sql,iostat=ios)
  read(4,nml=mp_pdf,iostat=ios)
  if(ios/=0) then
    print *,space, 'Error',ios ,'reading NML = mp_pdf from ', trim(file_inp),', check that you have MP_TOOLS at least ',version
    stop
  endif
  if(j_mode<=2) then
    print *,space, 'Modes J_MODE=1..2 not available with MP_SQL, setting J_MODE=6 for I(Q)'
    j_mode = 6
  endif
  
!!  allocate(ind_pseudo(n_atom,n_atom+n_pseudo_max+1),at_name_pseudo(n_pseudo_max+1),c_pseudo(n_pseudo_max+1),c_pseudo_mean(n_pseudo_max+1),m_pseudo(n_pseudo_max+1))       !1st pseudo is TOT by default
!!  ind_pseudo = 0
!!  at_name_pseudo = ''
!!  at_name_pseudo(1) = 'TOT'
!!  ind_pseudo = 0
!!  do j=1,n_atom
!!    ind_pseudo(j,j) = 1
!!  enddo
!!  ind_pseudo(:,n_atom+1) = 1
!!  c_pseudo = .0
!!  m_pseudo = 0

  allocate(ind_part(2,n_part_max),ind_ext(n_part_max),ext_scale(n_part_max),ext_dy(n_part_max))       !1st pseudo is TOT by default 
  allocate(ind_mo(n_atom),at_av_matrix(n_atom,n_atom))       !1st pseudo is TOT by default 
  ind_part = 0
  ind_ext = 0
  ext_scale = 1.
  ext_dy = .0
  
  if(j_acc==0.or.j_acc==1) then
    print *,space, 'Setting J_ACC = 3 to the recommended Gauss integration algorithm (check for other J_ACC choices in .PAR)'
  endif

! *** Read the atom positions       
  rewind(4)
  section = 'atoms'
  do
    read(4,'(a)',iostat=ios) string
    if(ios/=0) then
      print *,space, 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
      stop
    endif
    if(string(1:5).eq.section) exit	!find the atoms part of the .par file
  enddo
  do j=1,n_atom
    read(4,*) at_label(j),at_name_par(j),arg,arg,arg,arg,at_name_ext(j)	!at_base & conc come from the header
    if(at_name_ext(j)=='n'.or.at_name_ext(j)=='N') at_name_ext(j)=''
  enddo

! check for partial/mixed site occupation
  if(sum(abs(at_base))==.0) then
    ind_mo = 1         !all atoms are non-distinguishable (FZ type)
  else
    ind_mo = 0          !atoms are distinguishable (distinct lattice positions) in principle, look for mixed occupations on a lattice (not single bulk)
    
    if(product(n_row)>1) then
      jj = 1
      do i=1,n_atom
        ii = 0
        do j=1,n_atom
          if(at_occup_r(j)>.0.and.at_occup_r(j)<1..and.ind_mo(j)==0) then        
            if(sum(abs(at_base(j,:)-at_base(i,:)))==.0) then
              ind_mo(j) = jj                     ! set to 1 for mixed occupation (undistinguishable atoms)
              ii = ii+1
            endif
          endif
        enddo
        if(ii==0) then
          cycle
        elseif(ii==1) then
          do k=1,n_atom
            if(ind_mo(k)==jj) ind_mo(k) = 0     !skip if for partial occupation, not mixture
          enddo
        elseif(ii>0) then
          jj = jj+1
        endif
      enddo
    endif   
  endif

  if(j_verb==1) print *,space, 'ind_mo',ind_mo
  
  at_av_matrix = .0
!   if(sum(ind_mo)==0.or.product(ind_mo)==1) then
  do j=1,n_atom
    do i=1,n_atom
      if(ind_mo(i)==ind_mo(j).and.(i==j.or.ind_mo(i)/=0)) at_av_matrix(i,j) = 1.
    enddo
  enddo

      if(j_verb==1) then
        print *,space, 'at_av matrix'
        do ii=1,n_atom
            print *,space, ii,at_av_matrix(ii,:) 
        enddo      
        print *
      endif


  
! *** Read the partial PDF definitions - used also for S(Q) partials   
!   if(n_part_max>0.and.product(ind_mo)==0) then
  if(n_part_max>0) then
    rewind(4)
    section = 'partial'
    do
      read(4,'(a)',iostat=ios) string
      if(ios/=0) then
        print *,space, 'Section title: PARTIAL_PDF  not found (can be added in dialogue)'    !n_part,n_pseudo
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
        if (ind_part(1,j)<ind_part(2,j)) ind_part(:,j) = cshift(ind_part(:,j),1)      !the 2nd one should be smaller  
        n_part = j
      enddo
    else
      n_part = 0
    endif
  endif
  

!!! *** Read the pseudo atom definitions       
!!  if(n_pseudo_max>0) then
!!    rewind(4)
!!    section = 'pseudo'
!!    do
!!      read(4,'(a)',iostat=ios) string
!!      if(ios/=0) then
!!        print *,'Section title: PSEUDO_ATOMS not found (can be added in dialogue)'
!!        n_pseudo = 1
!!        found = .false.
!!        exit
!!      endif
!!      found = (string(1:6)==section)
!!      if(found) exit	                  !found the 'pseudo_atom' part of the .par file
!!    enddo
!!    
!!    if(found) then
!!      do j=2,n_pseudo_max
!!        read(4,*,iostat=ios) at_name_pseudo(j),ind_pseudo(1:n_atom,n_atom+j)	              ! pseudo_atom name and indices
!!        if(ios/=0) then
!!          n_pseudo = j-1
!!          exit
!!        endif
!!        n_pseudo = j
!!      enddo
!!    endif
!!  endif
!!  close(4)
!!
!!  if(sum(ind_mo)>0) then      !add pseudo atoms for mixed occupances
!!    do j=1,n_atom 
!!      if(ind_mo(j)>0) then
!!        ind_pseudo(j,n_atom+n_pseudo+ind_mo(j)) = 1            !+1 is for TOT
!!        c_pseudo(ind_mo(j)) = c_pseudo(ind_mo(j))+at_occup_r(j)
!!        m_pseudo(ind_mo(j)) = m_pseudo(ind_mo(j))+1
!!        if(n_atom+n_pseudo+ind_mo(j)>n_atom+n_pseudo_max) then
!!          print *,'Increase N_PSEUDO_MAX in .PAR to >',n_atom+n_pseudo+maxval(ind_mo)
!!          stop
!!        endif
!!      endif
!!    enddo
!!  endif
!!  
!!  c_pseudo_mean = c_pseudo(1:maxval(ind_mo))/real(m_pseudo(1:maxval(ind_mo)))
!! 
!! if(j_verb==1) print *,'c_pseudo_mean,c_pseudo(1:maxval(ind_mo))',c_pseudo_mean,c_pseudo(1:maxval(ind_mo))
!! 
!!  do j=1,maxval(ind_mo)
!!    write(at_name_pseudo(n_pseudo+j),'("P",I1.1)')j   !put names P1,P2 to mixed atoms
!!  enddo
!!  n_pseudo = n_pseudo+maxval(ind_mo)
!!
!!  if(j_verb==1) then
!!    print *,'n_pseudo',n_pseudo
!!    do j=1,n_pseudo
!!      print *,'j,at_name_pseudo(j),ind_pseudo(:,j)',j,at_name_pseudo(j),ind_pseudo(:,n_atom+j)
!!    enddo
!!  endif
!!  

! *** Check atom names against the .PAR input       
  do j=1,n_atom
    if(at_name_par(j)/=at_name_out(j)) then
      print *,space, 'Not-matching  atom names in .PAR and .DAT: ',j,at_name_par(j),at_name_out(j)
      print *,prompt, 'Prefer .DAT? (1/0)'
      read(*,*) ii
      if(ii==1) at_name_par = at_name_out
      exit 
    endif
  enddo			

! *** read neutron scattering lengths
  open(4,file='neutron_xs.txt',action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    open(4,file='/usr/local/mp_tools/ref/neutron_xs.txt',action='read',status ='old',iostat=ios)
    if(ios.ne.0) then
      do
        print *,prompt, 'File neutron_xs.txt not found, type in valid access path/filename'
        read(*,'(a)') x_file_name
        open(4,file=trim(x_file_name),action='read',status ='old',iostat=ios)
        if(ios==0) exit
        print *,space, 'File',trim(x_file_name),' not found, try again ...'
      enddo
    endif
  endif
  search_loop: do j=1,n_atom
    rewind(4)
    do i=1,210
      read(4,*) atom
      if(atom==trim(at_label(j))) then
        backspace(4)
        read(4,*) atom,b_coh(j)
        cycle search_loop
      endif
    enddo
    print *,space, 'b_coh for ',trim(at_label(j)),' not found,'
    print *,space, 'check your spelling and the neutron_xs.txt table; use unit weights'
  enddo search_loop
  
  b_coh = .1*b_coh  !convert b_coh from FM to 10^12 cm
  print *

! *** read Xray formfactor parameters 
  open(4,file='xray_ff.txt',action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    open(4,file='/usr/local/mp_tools/ref/xray_ff.txt',action='read',status ='old',iostat=ios)
    if(ios.ne.0) then
      do
        print *,prompt, 'File xray_ff.txt not found, type valid access path/filename'
        read(*,'(a)') x_file_name
        open(4,file=trim(x_file_name),action='read',status ='old',iostat=ios)
        if(ios==0) exit
        print *,space, 'File',trim(x_file_name),' not found, try again ...'
      enddo
    endif
  endif
  atom_loop: do j=1,n_atom
    rewind(4)
    do i=1,210
      read(4,*) atom
      if(atom==trim(at_label(j))//trim(at_name_ext(j))) then
        backspace(4)
        read(4,*) atom,x_ffpar(j,1:9)
        cycle atom_loop
      endif
    enddo
    print *,space, 'Xray formfactor for ',trim(at_label(j))//trim(at_name_ext(j)),' not found,'
    print *,space, 'check your spelling and the neutron_xs.txt table; use unit weights'
  enddo atom_loop
  close(4)
  
  string_in = subst_name
  print *,prompt, 'Substance name (confirm, &append or replace): ', string_in
  read(*,*) string_in
  string = trim(adjustl(string_in))
  if(string_in/=subst_name) then
    if(string_in(1:1)=='&') then
      subst_name = trim(subst_name)//string_in(2:)
    else
      subst_name = string_in
    endif
  endif

  print *
  print *,space, 'Substance name: ', subst_name
  print *,space, 'Sim_type, dat_type, input method: ',sim_type,dat_type,input_method		

! *** write overview of atom data
  print *
  print *,space, 'Atoms from ',trim(file_inp), ' (name, occupancy,neutron b_coh, Xray formfactor parameters)'
  if(input_method== 'BULK') then
    do j=1,n_atom
      write(*,'(1x,a,5x,a4,2x,2(f8.4,2x),20f8.4)')	space,at_name_par(j),at_occup_r(j),b_coh(j),x_ffpar(j,:)
      write(9,'(5x,a4,2x,2(f8.4,2x),20f8.4)')	at_name_par(j),at_occup_r(j),b_coh(j),x_ffpar(j,:)
    enddo	
  else
    do j=1,n_atom
      write(*,'(1x,a,5x,a4,3f8.4,2x,2(f8.4,2x),20f8.4)')	space,at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j),x_ffpar(j,:)
      write(9,'(5x,a4,3f8.4,2x,2(f8.4,2x),20f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j),x_ffpar(j,:)
    enddo	
  endif							

! *** Create the a_cell_inv matrix for older data (orthorhombic only)

  if(n_head<4) then                                 !generate a_cell_inv for orthogonal lattices (older data without a_cell)
    a_cell_inv = .0
    
    if(product(n_row)>1) then
      do k=1,3
        a_cell_inv(k,k) = 1./a_par(k)
      enddo
    else
       do k=1,3
        a_cell_inv(k,k) = 1.
      enddo
   endif
  endif

!!!! *** generate triclinic a_cell for test purposes
!!!      if(j_test==1) then
!!!        n_head = 4
!!! 
!!! 				print *,
!!! 				print *,'Test of geometric relations in real & reciprocal space'
!!! 				print *,
!!!
!!!        a_cell = transpose(reshape((/
!!!     1   3.5,-1.,.5, 											!rows are elementary translations in real space
!!!     2  -1.2,4.5,.5, 
!!!     3   .5,.5,4.1/),(/3,3/)))
!!!
!!!				print *,'a_cell - rows are elementary translations in real space'
!!!				do i=1,3
!!!					print *,a_cell(i,:)
!!!				enddo
!!!        e1 = a_cell(1,:) 
!!!        e2 = a_cell(2,:) 
!!!        e3 = a_cell(3,:) 
!!!        e1_norm = norm2(e1)
!!!        e2_norm = norm2(e2)
!!!        e3_norm = norm2(e3)
!!!
!!!        do j=1,3
!!!          a_par(j) = norm2(a_cell(j,:))
!!!        enddo 
!!!			  print *,'Lattice parameters',a_par
!!!
!!!				print *,'Axis angles - real space:'
!!!				angle(1) = dot_product(a_cell(2,:),a_cell(3,:))/(a_par(2)*a_par(3))
!!!				angle(2) = dot_product(a_cell(1,:),a_cell(3,:))/(a_par(1)*a_par(3))
!!!				angle(3) = dot_product(a_cell(1,:),a_cell(2,:))/(a_par(1)*a_par(2))
!!!				angle = acos(angle)
!!!
!!!				print *,'angle_rad',angle
!!!				print *,'angle_deg',angle*180./pi
!!!			
!!!
!!!				print *,'Reciprocal space basis - from A_CELL inversion:'
!!!				a_cell_1 = a_cell
!!!				call gjinv(a_cell_1,3,3,a_cell_inv,3,ier)
!!!				if(ier==1) then
!!!					print *,'Singular cell vector matrix, check your HISTORY file!'
!!!					stop
!!!				endif
!!!				print *,'a_cell_inv  - columns are elementary translations in reciprocal space'
!!!				do i=1,3
!!!					print *,a_cell_inv(i,:)
!!!				enddo
!!!        
!!!
!!!        print *,'Verify A_CELL inversion'
!!!        print *,'check by matmul(a_cell,a_cell_inv)',matmul(a_cell,a_cell_inv)
!!!        print *,'check by matmul(a_cell_inv,a_cell)',matmul(a_cell_inv,a_cell)
!!!
!!!
!!!        cell_volume = dot_product(vector_product(e1,e2),e3)
!!!        print *,'Real space cell_volume by mixed product',cell_volume
!!!
!!!        cell_angle = sqrt(1.-cos(angle(1))**2-cos(angle(2))**2-cos(angle(3))**2+2.*cos(angle(1))*cos(angle(2))*cos(angle(3)))  !abc*sqrt(1−cos2α−cos2β−cos2γ+2cosαcosβcosγ) 
!!!        cell_volume = e1_norm*e2_norm*e3_norm*cell_angle
!!!       
!!!        print *,'Real space cell_volume by ""SQRT(COS)"" formula',cell_volume
!!!      
!!!				 print *
!!!				 print *,'Reciprocal space basis - vector_product formulas'
!!!				 e1_r = vector_product(e2,e3)/cell_volume
!!!				 e2_r = vector_product(e3,e1)/cell_volume
!!!				 e3_r = vector_product(e1,e2)/cell_volume
!!!			 
!!!				 a_cell_inv(:,1) = e1_r
!!!				 a_cell_inv(:,2) = e2_r
!!!				 a_cell_inv(:,3) = e3_r
!!!			 
!!!				 print *,'e1_r',e1_r
!!!				 print *,'e2_r',e2_r
!!!				 print *,'e3_r',e3_r
!!!				 print *,
!!!			 
!!!				 print *
!!!				 print *,'Axis angles - reciprocal space - dot_product'
!!!				 angle_r(1) = acos(dot_product(e2_r,e3_r)/(norm2(e2_r)*norm2(e3_r)))
!!!				 angle_r(2) = acos(dot_product(e3_r,e1_r)/(norm2(e3_r)*norm2(e1_r)))
!!!				 angle_r(3) = acos(dot_product(e1_r,e2_r)/(norm2(e1_r)*norm2(e2_r)))
!!!				 print *,'angle_r',angle_r
!!!       
!!!				 print *
!!!				 print *,'Angles - reciprocal space - COS formulas'
!!!				 angle_r(1) = acos((cos(angle(2))*cos(angle(3))-cos(angle(1)))/(sin(angle(2))*sin(angle(3))))
!!!				 angle_r(2) = acos((cos(angle(3))*cos(angle(1))-cos(angle(2)))/(sin(angle(3))*sin(angle(1))))
!!!				 angle_r(3) = acos((cos(angle(1))*cos(angle(2))-cos(angle(3)))/(sin(angle(1))*sin(angle(2))))
!!!				 print *,'angle_r',angle_r
!!!       
!!!       
!!!!      metrics matrices
!!!       
!!!					print *
!!!				  print *,'Metrics matrices'
!!!			    print *,'1 G_matrix - via basis vectors - real space metrics'
!!!				  forall (i=1:3,j=1:3)
!!!					  g_matrix(i,j) = dot_product(a_cell(i,:),a_cell(j,:))
!!!				  end forall
!!! 					do i=1,3
!!!						print *,g_matrix(i,:)
!!!					enddo
!!!      
!!!					print *,'G*_matrix via G^-1 - reciprocal space metrics'
!!!					a_cell_1 = g_matrix
!!!					call gjinv(a_cell_1,3,3,g_r_matrix,3,ier)
!!!					if(ier==1) then
!!!						print *,'Singular g_matrix, check ...'
!!!						stop
!!!					endif
!!!					do i=1,3
!!!						print *,g_r_matrix(i,:)
!!!					enddo
!!!        
!!!				  print *,
!!!				  print *,'G_matrix  - via angles - real space metrics'
!!!					g_matrix(1,1) = e1_norm**2  
!!!					g_matrix(2,2) = e2_norm**2  
!!!					g_matrix(3,3) = e3_norm**2  
!!!					g_matrix(2,1) =  e1_norm*e2_norm*cos(angle(3))
!!!					g_matrix(3,1) =  e1_norm*e3_norm*cos(angle(2))
!!!					g_matrix(3,2) =  e2_norm*e3_norm*cos(angle(1))
!!!					g_matrix(1,2) = g_matrix(2,1)
!!!					g_matrix(1,3) = g_matrix(3,1)
!!!					g_matrix(2,3) = g_matrix(3,2)
!!!					do i=1,3
!!!						print *,g_matrix(i,:)
!!!					enddo
!!!
!!! 			    print *,'2 G*_matrix  - via reciprocal basis vectors - reciprocal space metrics'
!!!					forall (i=1:3,j=1:3)
!!!						g_r_matrix(i,j) = dot_product(a_cell_inv(:,i),a_cell_inv(:,j))
!!!					end forall
!!!					do i=1,3
!!!						print *,g_r_matrix(i,:)
!!!					enddo
!!!
!!!! *** dot product
!!!					 print *,
!!!					 print *,'Dot product test'
!!!					 print *,'real space: (3,2,1).(1,2,3) explicit', dot_product(3.*e1+2.*e2+e3,1.*e1+2.*e2+3.*e3)
!!!					 print *,'real space: (3,2,1).(1,2,3) G-matrix', dot_product((/3.,2.,1./),matmul(g_matrix,(/1.,2.,3./)))
!!!
!!!					 print *,'reciprocal space: (3,2,1).(1,2,3) explicit', dot_product(3.*e1_r+2.*e2_r+e3_r,1.*e1_r+2.*e2_r+3.*e3_r)
!!!					 print *,'reciprocal space: (3,2,1).(1,2,3) G*-matrix',dot_product((/3.,2.,1./),matmul(g_r_matrix,(/1.,2.,3./)))
!!!					 print *,
!!!					 print *,
!!!
!!!         endif !j_test
!!!   

! *********************  OpenMP initialization start  *******************************      
!
  proc_num_in = j_proc
  if(j_verb.ge.1) print *,space, 'PAR proc_num_in          = ',proc_num_in
  thread_num_max = omp_get_max_threads( )								!this gives maximum number of threads available (limited by concurrent tasks??)
  proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
  if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
  call omp_set_num_threads(proc_num_in)
  
  if(j_verb.ge.1) then
    write (*,*) 'OMP processes available  = ', proc_num
    write (*,*) 'OMP threads maximum      = ', thread_num_max
    write (*,*) 'OMP processes requested  = ', proc_num_in
  endif


  if(proc_num_in.gt.1) then
    write(9,*) 'OMP threads maximum = ', thread_num_max
    write(9,*) 'OMP processes requested 	= ', proc_num_in
  else
    write(9,*) 'OpenMP not in use'
    print *,space, 'OpenMP not in use'
  endif
  write(9,*) 
!
! ********************* OpenMP initialization end *******************************      
!

!
! *** initialise the random_number generator for Mont-Carlo integration
!
  if(j_acc==2) then
    call random_seed(size=seed_size)               !Puts size of seed into seed_size
    allocate(numbers(seed_size))

    if(j_rand==0) then                     
      print *,space, 'Random_number: j_rand =',j_rand,'  the system will supply unique, machine dependent seeds each time this code runs'
      call random_seed(get=numbers)               !Gets actual seeds 
      if(j_verb==1) then
        print *,space, 'Random_number seed size:',seed_size
        print *,space, 'Random_seed:',numbers
        print *,space, 'Reference 1st 5 random numbers:',(rnd(i),i=1,5)
      endif
    elseif(j_rand==1) then                     !if j_rand>1 generate a seed for later reference & numerical reproducibility checks
      print *,space, 'Random_number: j_rand =',j_rand,'  the system will supply k-dependent standard seeds for each of the OMP threads (use only for testing the consistence of OMP_on/OMP_off results)'
    elseif(j_rand>1) then                     !if j_rand>1 generate a seed for later reference & numerical reproducibility checks
      print *,space, 'Random_number: j_rand =',j_rand,'  this seeding reference can be used to exactly reproduce this MC-run later on'
     idum = j_rand
      numbers(1) = huge(1)*rand(idum)            !a dry call to initialise ran0
      do i=1,seed_size
        numbers(i) = huge(1)*rand(idum)          !use a trivial random number generator to produce the seeds (they could even be all the same small ones, but ...)
      enddo
      call random_seed(put=numbers)       !Produce a seed to start, this permits to reproduce exactly the same results on the same system 
      do i=1,5
        call random_number(rnd(i))
      enddo
      if(j_verb==1) then
        print *,space, 'Random_seed:',numbers
        print *,space, 'Reference 1st 5 random numbers:',(rnd(i),i=1,5)
      endif
    endif
    print *
  endif

!
! **** cycle over snapshot files to accumulate input data
!
  print *,space, 'Input files:'

  CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
  sc_r = 1./sc_r

  allocate(at_ind_in(4*n_tot),ind_at(n_atom))			

  allocate(at_pos_in(4*n_tot,nfile))			

  ifile = 0
  file_loop: do jfile=nfile_min,nfile_max,nfile_step

! ***  open the binary MD snapshot file
    if(t_single)then
      write(file_dat,'("./data/",a,".dat")') trim(file_master)
    else
      if(jfile<=9999) then
        write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),jfile
      elseif(jfile>=10000) then
        write(string,'(i8)') jfile
        file_dat = './data/'//trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
      endif
    endif

    open(1,file=file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
    if(ios.ne.0) then
      print *,space, 'File ',trim(file_dat),' not opened! IOS =',ios
    print *,prompt, 'Skip this one (2), skip the rest (1), stop execution(0)?'
    read(*,*) jj
    if(jj==2) cycle file_loop
    if(jj==1) exit file_loop
    if(jj==0) stop
    endif
    ifile = ifile+1				
    if(jfile==nfile_min.or.ifile==100*(ifile/100)) print *,space, file_dat

    ind_at(1) = 0
    do j = 2,n_atom
      ind_at(j) = ind_at(j-1)+nsuper_r(j-1)
    enddo

    i_rec = n_head	
    do j=1,n_rec-1
      i_rec = i_rec+1
      if(jfile==nfile_min) read(1,rec=i_rec) at_ind_in((j-1)*l_rec+1:j*l_rec)			!only needed for 1_phonon with CELL 
    enddo	
    i_rec = i_rec+1
    if(jfile==nfile_min) read(1,rec=i_rec) at_ind_in((n_rec-1)*l_rec+1:4*n_tot)			

    do j=1,n_rec-1
      i_rec = i_rec+1
      read(1,rec=i_rec) at_pos_in((j-1)*l_rec+1:j*l_rec,ifile)			
    enddo	
    i_rec = i_rec+1
    read(1,rec=i_rec) at_pos_in((n_rec-1)*l_rec+1:4*n_tot,ifile)					
    close(1)
  enddo file_loop

  if(ifile.ne.nfile) then						!update nfile in case of shorter sequence
    nfile = ifile
  endif

  CALL SYSTEM_CLOCK (COUNT = sc_c2)
  print *,space, nfile,' files read in ',(sc_c2-sc_c1)*sc_r, ' sec SYS time'
 
  if(nsuper==1) then
    at_pos_min = -a_par/2.      !for n_row=1 BULK a_par represents the whole box
    at_pos_max = a_par/2.
  else
    at_pos_min = -n_row/2.
    at_pos_max = +n_row/2.         
  endif
  if(j_verb==1) print *,space, 'at_pos_min,at_pos_max',at_pos_min,at_pos_max


  jfile = 1				!now to be used as index for the successive output files			

  e1 = (/1,0,0/)
  e2 = (/0,1,0/)
  e3 = (/0,0,1/)


! **** map loop (goes till the very end)
!
  bz_loop: do

! *** define the momentum space range          
    print *,prompt, 'Q-range [Å-1] (0=END)'		
    read(*,*) q_range
    if(q_range==0.) exit bz_loop	

    bz_n = 2.*q_range
    
! **** prepare reciprocal lattice vectors for plotting etc.

    e1_r = a_cell_inv(:,1)					!elementary translation vectors spanning the reciprocal lattice, expressed in Euclidian RS coordinates, multiply by 2Pi is needed
    e2_r = a_cell_inv(:,2)
    e3_r = a_cell_inv(:,3)
    
    e1_r_norm = norm2(e1_r)			
    e2_r_norm = norm2(e2_r)
    e3_r_norm = norm2(e3_r)

      if(j_verb==1) then
        print *,'e1_r,2Pi*e1_r_norm',e1_r,twopi*e1_r_norm			!twopi*e1_r_norm is the length of e1_r in A-1 on the usual scale containing the 2Pi factor
        print *,'e2_r,2Pi*e2_r_norm',e2_r,twopi*e2_r_norm
        print *,'e3_r,2Pi*e3_r_norm',e3_r,twopi*e3_r_norm
      endif
    
    e3v = vector_product(e1_r,e2_r)
    e2v = vector_product(e3_r,e1_r)
    e1v = vector_product(e2_r,e3_r)
     
! *** to get the BZ_NX etc. range is BZ_N divided by the height of the parallepiped in the corresponding direction        
    bz_nx = bz_n/(dot_product(e1_r,e1v)/norm2(e1v))  !+.5			!+.5 is a safety margin to avoid rounding effects
    bz_ny = bz_n/(dot_product(e2_r,e2v)/norm2(e2v))  !+.5
    bz_nz = bz_n/(dot_product(e3_r,e3v)/norm2(e3v))  !+.5
    if(j_verb==1) print *,space, 'bz_n,bz_nx,bz_ny,bz_nz',bz_n,bz_nx,bz_ny,bz_nz

    d_q = twopi/(at_pos_max-at_pos_min)				!pseudo-cubic phase steps for NUFFT defined by the box size
    
    if(j_verb==1) print *,space, 'd_q,q_step,q_range',d_q,'   ',q_step,q_range

    n_qq(1) = 2*nint((.5*bz_nx)/d_q(1))+1       
    n_qq(2) = 2*nint((.5*bz_ny)/d_q(2))+1       
    n_qq(3) = 2*nint((.5*bz_nz)/d_q(3))+1   
    i_centre = n_qq/2+1
    if(j_verb==1) print *,space, 'n_qq,i_centre',n_qq,'   ',i_centre
    print *	
    
! **** test the size of Q-grid to fit into the 4Gb limit
!
    if(product(n_qq)*n_atom>huge(n_qq(1))) then
      t2 = product(real(n_qq))/(real(huge(n_qq(1))/real(n_atom)))
      print *,space, 'Q-range too large by a factor:',t2**.3333333	
      cycle bz_loop
    endif
  
    n_mem = 1+8.*(n_atom+4.)*(product(real(n_qq))/1.e9)					! 16*product(n_qq) is the size of NUFFT output (NUFFT needs 2x that), n_atom*8*product(n_qq) is needed for partial amplitudes; 1Gb is the rest (input etc.)
    write(*,'(1x,a,"NOTE: NUFFT_3D requires about",i4," Gb memory, in case of problems quit other applications or reduce Q-range!")') space,n_mem
    print *
    print *,space, 'Preparing FT ...'		


! *** generate look-up table for Q(i,j,k) norms (non-identical on (i,j,k) permutations in non-orthogonal systems)
    allocate(q_norm(n_qq(1),n_qq(2),n_qq(3)))				!,cross_sec
     
    do i=1,n_qq(1)
      do j=1,n_qq(2)
        do k=1,n_qq(3)
          q_norm(i,j,k) = norm2(matmul(a_cell_inv,((/i,j,k/)-i_centre)*d_q))
        enddo
      enddo
    enddo
    
! generate look-up table of Xray form factors
    n_hist = 1+q_range/q_step
    allocate(q_hist(n_hist),x_ffq(n_atom,n_hist))				!,cross_sec
  
    do i=1,n_hist                                           ! here q_hist contains the factor 2Pi
      q_hist(i) = (i-1)*q_step															! the form factor formula contains a*exp(-b*(q/(4*Pi))**2)				
      q_sq = .25*(q_hist(i)/twopi)**2																! the form factor formula contains a*exp(-b*(q/(4*Pi))**2)				
      do j_at=1,n_atom
        x_ffq(j_at,i) = x_ffpar(j_at,9)
        do ii=1,4
          x_ffq(j_at,i) = x_ffq(j_at,i)+x_ffpar(j_at,2*ii-1)*exp(-x_ffpar(j_at,2*ii)*q_sq)
        enddo
      enddo
    enddo       

! *** initialise the MC integration procedure
    if(j_acc==2) then
      n_h_max = int(.5*product(real(n_qq))/real(n_mc))   !(1+n_tot/n_mc)*n_atom*product(n_row)/2 product(n_row)/2 is roughly the volume of sphere with radius of pdf_range_max !n_mc=1e6
      write(string,'(i4)') n_h_max
      if(n_h==0) then
        print *,prompt, 'MC sampling pairs per frame ([x 10^6], ',trim(adjustl(string)),' max):'
        read(*,*)   n_h
      endif
      if(n_h.gt.n_h_max) then                   !n_h_max is in units of n_mc to avoid overflow for large boxes
        print *,space, 'WARNING: n_h exceeds max number of 10^6 atom pairs ',n_h_max
        print *,prompt, 'type in a reduced n_h (<',n_h_max,'):'
        read(*,*)   n_h
      endif
      n_int = n_h*n_mc                           !number of MC cycles per snapshot
      int_mode = 'Monte Carlo'
      print *,space, trim(int_mode),' integration over',n_h,'*10^6 cell pairs'
    endif
      

! **** Start the FT and integration cycle over model frames
!
    allocate(sq_hist_tot(n_hist),sq_plot_tot(n_hist),sq_hist(n_atom,n_atom,n_hist),sq_plot(n_atom,n_atom,n_hist))	
    allocate(ampl_atom_3d(n_qq(1),n_qq(2),n_qq(3),n_atom),at_weight(n_atom),at_scf(n_atom),ampl_q(n_atom))   

    n_qq8 = n_qq
    sq_hist_tot = .0
    sq_hist = .0
  
    print *,space, 'Doing FT ...'		

    call cpu_time(t0)				
    call cpu_time(t1)				
    CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
    sc_r = 1./sc_r

    frame_loop: do ifile=1,nfile

      if(input_method=='CELL') then
        at_pos_file(1:4,1:nsuper,1:n_atom) => at_pos_in(:,ifile)
        at_pos_file_bulk(1:4,1:1) => at_pos_in(1:4,ifile)       !this is pro forma because of OMP
      else
        at_pos_file(1:4,1:1,1:1) => at_pos_in(1:4,ifile)       !this is pro forma because of OMP
        at_pos_file_bulk(1:4,1:n_tot) => at_pos_in(:,ifile)
      endif				
 
      do j = 1,n_atom
        n_fft = n_sup(j)			!we need n_fft 64bit
        allocate(xx(n_sup(j),2),cf(n_sup(j)),ampl_tot(n_qq(1)*n_qq(2)))																	

        do k = 1,n_qq8(3)    !cycle over horizontal layers
          kk = k-(n_qq8(3)+1)/2  
          if(input_method=='CELL') then
            ii = 0
            do i=1,nsuper
              non_zero = (at_pos_file(1,i,j)/=.0.or.at_pos_file(2,i,j)/=.0.or.at_pos_file(3,i,j)/=.0)
              if(non_zero) then
                ii = ii+1				!avoiding "empty" positions
                xx(ii,1:2) = at_pos_file(1:2,i,j)*d_q(1:2)
                cf(ii) = cexp((0.,-1.)*kk*d_q(3)*at_pos_file(3,i,j))
              endif
            enddo              
          else					!BULK
            do i=1,n_sup(j)
              xx(i,1:2) = at_pos_file_bulk(1:2,ind_at(j)+i)*d_q(1:2)
              cf(i) = cexp((0.,-1.)*kk*d_q(3)*at_pos_file_bulk(3,ind_at(j)+i))
            enddo							
          endif		!input_method

          call finufft2d1(n_fft,xx(:,1),xx(:,2),cf,iflag,eps_fft,n_qq8(1),n_qq8(2),ampl_tot,nul_opt,ier)
         
         ampl_atom_3d(:,:,k,j) = reshape(source=ampl_tot,shape=[n_qq8(1),n_qq8(2)])     !/sqrt(real(n_fft))												
        enddo  !k

        deallocate(xx,cf,ampl_tot)						
      enddo  !j (n_atom)
 
      call cpu_time(t2)
      CALL SYSTEM_CLOCK (COUNT = sc_c2)
      print *	
      print *,space, 'FINUFFT on',ifile,'snapshot CPU_TIME',t2-t1,'  SYS_TIME',(sc_c2-sc_c1)*sc_r
      print *					
   
      call cpu_time(t1)				
    CALL SYSTEM_CLOCK (COUNT = sc_c1)		

    if(j_acc==2) then																	!do MC integration
        print *,space, 'Doing MC integration ...'
        allocate(sq_hist_k(n_atom,n_atom,n_hist,n_h))
        sq_hist_k = 0
        n_skip1 = 0
        n_skip2 = 0
        ampl_atom_3d = ampl_atom_3d/sqrt(real(sum(n_sup)))

!!!!CCC OMP continuation line must have a valid character in column 6  CCCCCCC
!$omp parallel shared(sq_hist_k,ampl_atom_3d,q_norm,q_step,q_range,n_h,n_atom,n_qq,i_centre,n_int,j_rand,n_skip1,n_skip2)&
!$omp& private(q_pos,q_pos_norm,d_ind,i_hist,idum,rand0,rand1,numbers,k,j,ii,jj,iii,jjj,kkk)
!$omp do

      do k=1,n_h                           
        if(j_rand==1) then                 ! if j_rand=1 produce k-dependent standard seeds for each of the threads and cycles to test consistence of OMP_on and OMP_off results
          idum = j_rand+k
          numbers(1) = huge(1)*rand(idum)            !a dry call to initialise ran0
          do i=1,seed_size
            numbers(i) = huge(1)*rand(idum)          !use a trivial random number generator to produce the seeds (they could even be all the same small ones, but ...)
          enddo
          call random_seed(put=numbers)
          do i=1,5
            call random_number(rand0(i))
          enddo
          if(j_verb==1) then
            print *,space, 'k =',k
            print *,space, 'Random_seed:',numbers
            print *,space, 'Reference 1st 5 random numbers:',(rand0(i),i=1,5)
          endif
        endif
          
        do j = 1,n_mc
          do
            do i=1,3
              call random_number(rand1)		!generate the position within a cube [-1,1]; only points within the unit sphere will be retained
              q_pos(i) = 2.*rand1-1.
            enddo
            q_pos_norm = norm2(q_pos)
            if(q_pos_norm>.0.and.q_pos_norm<=1.) exit     !only q_pos within unit sphere is kept !throwing away also the unlikely event of q_pos=.0, not d_ind=0!
          enddo

          q_pos = q_pos/q_pos_norm        !now we project random points onto spherical surface by d_ind/norm2(d_ind) ....
          call random_number(rand1)

          d_ind = nint(((n_qq-i_centre)*q_pos)*rand1)

          if(sum(abs(d_ind))==0) then
            n_skip1 = n_skip1+1
            cycle	
          endif		

          iii = i_centre(1)+d_ind(1)
          jjj = i_centre(2)+d_ind(2)
          kkk = i_centre(3)+d_ind(3)

          if(q_norm(iii,jjj,kkk)<=q_range) then
            i_hist = nint(q_norm(iii,jjj,kkk)/q_step)+1    !Q=0 is the 1st channel
            do ii=1,n_atom
              do jj = 1,n_atom
                sq_hist_k(ii,jj,i_hist,k) = sq_hist_k(ii,jj,i_hist,k)+ampl_atom_3d(iii,jjj,kkk,ii)*conjg(ampl_atom_3d(iii,jjj,kkk,jj))
              enddo
            enddo
          else
            n_skip2 = n_skip2+1             
            cycle 
          endif
        enddo   !j=1,n_mc
      enddo		!k=1,n_h
!$omp end do
!$omp end parallel
              
      if(j_verb==1) print *,space, 'MC integration: n_skip1,n_skip2',n_skip1,n_skip2
      
      do k=1,n_h
        sq_hist = sq_hist+sq_hist_k(:,:,:,k)
      enddo
      deallocate(sq_hist_k)
     
    else                     !doing Gauss integration
      print *,space, 'Doing Gauss integration ...'
      allocate(sq_hist_k(n_atom,n_atom,n_hist,n_qq(1)))
      sq_hist_k = 0

      ampl_atom_3d = ampl_atom_3d/sum(n_sup)
              
!$omp parallel shared(n_qq,n_atom,q_norm,q_range,q_step,sq_hist_k,ampl_atom_3d) private(i_hist)
!$omp do
        do i=1,n_qq(1)
          do j=1,n_qq(2)
            do k=1,n_qq(3)            
              if(q_norm(i,j,k)<=q_range) then
                i_hist = nint(q_norm(i,j,k)/q_step)+1    !Q=0 is the 1st channel
                do ii=1,n_atom
                  do jj = 1,n_atom
                    sq_hist_k(ii,jj,i_hist,i) = sq_hist_k(ii,jj,i_hist,i)+ampl_atom_3d(i,j,k,ii)*conjg(ampl_atom_3d(i,j,k,jj))
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel

        do k=1,n_qq(1)
          sq_hist = sq_hist+sq_hist_k(:,:,:,k)
        enddo
        deallocate(sq_hist_k)
      endif   !j_acc  					
    enddo frame_loop  
    

    deallocate(ampl_atom_3d)						

!
!! *** NOTE: by now the asymptotes of sq_hist are at_occup_r(ii) for the diagonal terms and 0 for the others 
!! *** NOTE: we shall remove at_occup_r(ii) now  
!
! *** normalise the S(Q) to the Faber-Ziman form
				
    if(j_acc==2) then																	!normalise the MC integral
    
      do ii=1,n_atom
        do jj = 1,n_atom
          sq_hist(ii,jj,2:n_hist) = sq_hist(ii,jj,2:n_hist)/(q_step*nfile)
          sq_hist(ii,jj,2:n_hist) = sq_hist(ii,jj,2:n_hist)*q_range/(n_int-n_skip1-n_skip2)          
        enddo
      enddo				
    else                                                !normalise the Gauss integral

      do ii=1,n_atom
        do jj = 1,n_atom
          sq_hist(ii,jj,2:n_hist) = sq_hist(ii,jj,2:n_hist)/(q_step*q_hist(2:n_hist)**2)
          sq_hist(ii,jj,2:n_hist) = sq_hist(ii,jj,2:n_hist)*2.*pi**2/(nfile*at_volume)          
        enddo
      enddo
    endif
  
    call cpu_time(t2)
    CALL SYSTEM_CLOCK (COUNT = sc_c2)
    print *,space, 'S(Q) histogram CPU_TIME',t2-t1,'  SYS_TIME',(sc_c2-sc_c1)*sc_r
    print *					


! **** start the weight loop to form the S(Q) output 
!
! *** set atom amplitudes and weights
!
    at_mask = 1
    at_scf = 1.
    if(n_part==0) ind_part = 0
    part_scale = 1.
    tot_scale = 1.
    
    if(j_mode<=5.and.product(ind_mo)/=1) then
      print *,space, 'Modes J_MODE=3..5 only available for fully disordered systems, setting J_MODE=6 for I(Q)'
      j_mode = 6
    endif
  

    
    allocate(at_weight_matrix(n_atom,n_atom),at_mask_matrix(n_atom,n_atom))
    
    at_weights_loop: do		
      
      do j_atom=1,n_atom
        if(j_weight<=1.) then
          at_weight(j_atom) = 1.
        elseif(j_weight==2) then					!neutrons
          at_weight(j_atom) = b_coh(j_atom)
        elseif(j_weight==3) then					!neutrons
          at_weight(j_atom) = x_ffq(j_atom,1)          !just for overall normalisation
        else
          at_weight(j_atom) = 1.
        endif
      enddo

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
!         print *,space, 'at_weight',at_weight
!         print *,space, 'at_weight_sq_av,at_weight_av_sq',at_weight_sq_av,at_weight_av_sq
        print *,space, 'at_weight matrix'
        do ii=1,n_atom
            print *,space, ii,at_weight_matrix(ii,:) 
        enddo      
        print *
      endif


! *** average weights for mixed occupancies

     do j=1,n_atom
       at_weight(j) = sum(at_weight*at_occup_r*at_av_matrix(j,:))/sum(at_occup_r*at_av_matrix(j,:))     
     enddo

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

! *** generate the smoothing profile (Gauss) and ...
!
      n_smooth_fwhm = 1
      write(*,"(1x,a,'Gaussian smooth FWHM in Q_steps of',f6.3,'[Å-1] (1 no smoothing = default)')") space,q_step
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

! *** ... apply it
!       if(n_smooth>1) print *,space, 'Applying Gaussian smoothing with FWHM=',n_smooth_fwhm,' steps of',q_step,'[A-1]'
      write(smooth,'("Smooth FWHM",f6.3," [A-1]")') q_step*n_smooth_fwhm
      
      if(n_smooth==1) then
        do i = int(.5/q_step),n_hist      
!           sq_plot(:,:,i) = sq_hist(:,:,i)*at_weight_matrix*at_mask_matrix
          sq_plot(:,:,i) = sq_hist(:,:,i)
        enddo
      else
        sq_plot_tot = .0
        sq_plot = .0        
        do i = int(.5/q_step),n_hist      
          do j = 1,n_smooth
            if (i-n_smooth/2+j-1.le.2) cycle  !leave zeros out of range
            if (i-n_smooth/2+j-1>=n_hist) then
               sq_plot(:,:,i)= sq_plot(:,:,i)+sq_hist(:,:,n_hist-j)*f_smooth(j)      ! sq_plot(:,:,n_hist) is halfway to 0, ≈anything else is good enough
            else
              sq_plot(:,:,i)= sq_plot(:,:,i)+sq_hist(:,:,i-n_smooth/2+j-1)*f_smooth(j)
            endif
          enddo 
        enddo
      endif
      deallocate(f_smooth)
     
      sq_plot(:,:,1:int(.5/q_step)) = .0

      if(j_weight<=2) then
        if(j_mode/=5) then
          do i=int(.5/q_step),n_hist 
            sq_plot(:,:,i) = sq_plot(:,:,i)*at_weight_matrix*at_mask_matrix
          enddo
        endif

        if(j_mode==6) sq_plot = sq_plot/at_weight_sq_av
!             if(j_mode==7) nothing to do
        
        if(product(ind_mo)==1.and.j_mode<5) then    !mixture of indistinguishable atoms on a single sublattice (or amorphous)
          do ii=1,n_atom
            sq_plot(ii,ii,2:n_hist) = sq_plot(ii,ii,2:n_hist)-at_occup_r(ii)*at_mask(ii)*(at_scf(ii)*at_weight(ii))**2  !remove the b^2 asymptote
          enddo
          do jj=1,n_atom
            do ii=1,n_atom 
!               if(j_mode==3)   !NOTHING ELSE TO DO for F(Q)
              if(j_mode==4) sq_plot(ii,jj,:) = sq_plot(ii,jj,:)+at_occup_r(ii)*at_occup_r(jj)*at_weight_matrix(ii,jj)*at_mask_matrix(ii,jj)      !add <b>^2 asymptote for structure factor S(Q)
              sq_plot(ii,jj,:) = sq_plot(ii,jj,:)/at_weight_av_sq
            enddo
          enddo

        elseif(product(ind_mo)==1.and.j_mode==5) then
          do jj=1,n_atom
            do ii=1,n_atom 
             if(ii==jj) sq_plot(ii,jj,:) = sq_plot(ii,jj,:)-at_occup_r(ii)
             sq_plot(ii,jj,:) = sq_plot(ii,jj,:)/(at_occup_r(ii)*at_occup_r(jj))           !remove the occupancy scaling for FZ(Q), weights already =1
            enddo
          enddo       
        endif
      endif


      if(j_weight==3) then 
        if(j_mode/=5) then 
          do i=int(.5/q_step),n_hist 
            at_weight_av_sq = sum(at_occup_r*at_scf*x_ffq(:,i))*sum(at_occup_r*at_scf*x_ffq(:,i))        !mean amplitude squared per atom
            at_weight_sq_av = sum(at_occup_r*(at_scf*x_ffq(:,i)**2))                                     !mean squared amplitude per atom
            do ii=1,n_atom
              do jj = 1,n_atom
                sq_plot(ii,jj,i) = sq_plot(ii,jj,i)*(at_mask(ii)*at_scf(ii)*x_ffq(ii,i)*at_mask(jj)*at_scf(jj)*x_ffq(jj,i))
              enddo
            enddo

            if(j_mode==6) sq_plot(:,:,i) = sq_plot(:,:,i)/at_weight_sq_av
!             if(j_mode==7) nothing to do

            if(product(ind_mo)==1.and.j_mode<5) then    !mixture of indistinguishable atoms on a single sublattice (or amorphous)
              do ii=1,n_atom
                do jj = 1,n_atom
                  if(ii==jj) sq_plot(ii,jj,i) = sq_plot(ii,jj,i)-at_occup_r(ii)*at_mask(ii)*(at_scf(ii)*x_ffq(ii,i))**2         !remove the b^2 asymptote
                  if(j_mode==4) sq_plot(ii,jj,i) = sq_plot(ii,jj,i)+at_occup_r(ii)*at_mask(ii)*at_scf(ii)*x_ffq(ii,i)*at_occup_r(jj)*at_mask(jj)*at_scf(jj)*x_ffq(jj,i)
!                 if(j_mode==3)   !NOTHING ELSE TO DO for F(Q)
                enddo
              enddo
              sq_plot(:,:,i) = sq_plot(:,:,i)/at_weight_av_sq
            endif
          enddo

        elseif(product(ind_mo)==1.and.j_mode==5) then
          do jj=1,n_atom
            do ii=1,n_atom 
              if(ii==jj) sq_plot(ii,jj,:) = sq_plot(ii,jj,:)-at_occup_r(ii)
              sq_plot(ii,jj,:) = sq_plot(ii,jj,:)/(at_occup_r(ii)*at_occup_r(jj))           !remove the occupancy scaling for FZ(Q), weights already =1
            enddo
          enddo       
        endif  !(j_mode<=5.and.product(ind_mo)==1)
      endif

      sq_plot_tot = .0

      if(j_mode/=5) then
        do ii=1,n_atom
          do jj = 1,n_atom
              sq_plot_tot(:) = sq_plot_tot(:)+sq_plot(ii,jj,:)
          enddo
        enddo
        do j=1,n_atom
          do i=1,j  
            if(i/=j) sq_plot(i,j,:) = sq_plot(i,j,:)+sq_plot(j,i,:)
            if(i/=j) sq_plot(j,i,:) = sq_plot(i,j,:)
          enddo
        enddo
        tot_scale = 1.
      else
        tot_scale = .0
        part_scale = 1.
      endif

! **** start the output part: plot & text
!
      at_name_plot(1:n_atom) = at_name_par       
                    
      write(9,*)
      write(*,'(1x,a,"Atoms:         ",51(1x,a8))')  space,(at_name_plot(i),i=1,n_atom)
      write(9,'(1x,"Atoms:         ",51(1x,a8))')  (at_name_plot(i),i=1,n_atom)
      write(*,'(1x,a,"Atoms no.:  ",51(1x,i8))') space,(i,i=1,n_atom)
      write(9,'(1x,"Atoms no.:  ",51(1x,i8))') (i,i=1,n_atom)
      write(*,'(1x,a,"Occupancy:       ",51(1x,f8.4))') space,(at_occup_r(i),i=1,n_atom)                    
      write(9,'(1x,"Occupancy:       ",51(1x,f8.4))') (at_occup_r(i),i=1,n_atom)
      write(*,'(1x,a,"Amplitudes:      ",51(1x,f8.4))') space,(at_weight(i),i=1,n_atom)                    
      write(9,'(1x,"Amplitudes:      ",51(1x,f8.4))') (at_weight(i),i=1,n_atom)
      write(*,'(1x,a,"Scale factors:   ",51(1x,f8.4))') space,(at_scf(i),i=1,n_atom)                    
      write(9,'(1x,"Scale factors:   ",51(1x,f8.4))') (at_scf(i),i=1,n_atom)

      print *
      print *,space, 'Actual masks:'
      write(*,'(11x,(50i3))') (at_mask(i),i=1,n_atom)        
      print *
    
      if(minval(at_mask(1:n_atom))==1) then  !if all masks =1 don't put them into plot title
        masks = ''
      else
        write(masks,'("  Masks: ",50i1.1)') (at_mask(i),i=1,n_atom)			!to be used in plot titles
      endif

      x => q_hist

      x_end = x(n_hist) 
      x_start = x(1)
      i_start = int(.5/q_step)+n_smooth/2
      i_end = (x_end-x(1))/(x(2)-x(1))
      n_plot = i_end-i_start+1
				
! *** open the PGPLOT graphics window (X11)
    if (j_xserv.LE.0) then          
      j_xserv = PGOPEN('/xserv')
    else
      call PGSLCT(j_xserv)
    endif   
    if (j_xserv.LE.0) then    
      print *,space, 'Could not open PGPLOT /xserv'
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

    if(c_max==.0) then          !only do this at the real start of plot_loop
      c_min = amin1(minval(sq_plot_tot(i_start:n_hist/2)),minval(sq_plot(:,:,i_start:n_hist/2)))        !avoid maximum close to the origin
      c_min = .1*(int(10*c_min)-1)
      c_max = amax1(maxval(sq_plot_tot(i_start:n_hist/2)),maxval(sq_plot(:,:,i_start:n_hist/2)))
      c_max = .1*(int(10*c_max)+1)
    endif
          
    print *,space, 'Vertical scale c_min,c_max',c_min,c_max
    
    if (j_mode==5)then
      write(plot_header,'(a,"    ",a)') trim(file_dat_t0),trim(at_weight_scheme(1))         !implicit unit weights for Z(Q)
    else    
      write(plot_header,'(a,"    ",a)') trim(file_dat_t0),trim(at_weight_scheme(j_weight))
    endif
    
    plot_header = trim(subst_name)//'  '//trim(pdf_out(j_mode))//'  '//trim(plot_header)//trim(masks)      //'  '//trim(smooth)

    scale_loop: do
       
  call date_and_time(c_date,c_time,c_zone,i_time)
  write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)

      CALL PGSLCT(j_xserv)
      CALL PGSCI (1)  !white
      CALL PGSCH(1.)
      CALL PGSLW(2)
      CALL PGSLS(1)  !full
      CALL PGENV(x_start,x_end,c_min,c_max,0,j_grid+1)       !j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
      CALL PGLAB(x_label(j_mode),y_label(j_mode),trim(plot_header))  !put the axis labels
      CALL PGSCI(1)  !white
      CALL PGSLW(5)			!operates in steps of 5

      x_plot = x_start+.75*(x_end-x_start)
      y_plot = c_max-.1*(c_max-c_min)
      
      if(j_ext==1) then
        CALL PGTEXT (x_plot,y_plot,trim(file_inp))
        do j=1,n_part_ext
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
        CALL PGLINE(n_plot,x(i_start:i_end),tot_scale*sq_plot_tot(i_start:i_end))  !plots the curve         
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
        if (j_mode==5)then
          CALL PGSLW(5)
        else
          CALL PGSLW(2)
        endif
        CALL PGLINE(n_plot,x(i_start:i_end),part_scale(j)*sq_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
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

       
      print *,prompt, 'Adjust vertical scale (min, max) (0 0 EXIT, -1 -1 to adjust plot range)'
      c1 = c_min
      c2 = c_max
      read(*,*) c1,c2
      if(c1==.0.and.c2==.0) then
        exit
      elseif(c1==-1.and.c2==-1) then              
        print *,prompt, 'Confirm/adjust plot range, max =',n_hist*(x(2)-x(1))
        write(*,'("x_start, x_end ",2f7.1,":  ")',advance='no') x_start,x_end 
        read(*,*) x_start,x_end 
        if(x_start.lt.q_step) x_start = x(1)
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
          write(file_ps,'(a,"_sql_",i2.2,a)') trim(file_master),jfile,trim(pg_ext)
          write(file_res,'(a,"_sql_",i2.2,".txt")') trim(file_master),jfile
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
          file_res = trim(file_master)//'_sql'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//'.txt'							
          file_ps  = trim(file_master)//'_sql'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//trim(pg_ext)
        endif

        inquire(file=file_ps,exist=found_ps)
        inquire(file=file_res,exist=found_txt)
        if(.not.found_txt.and.(.not.found_ps)) exit				
        jfile = jfile+1
        if(jfile==100) then
          print *,prompt,'Tidy up .txt/.ps files to restart count from 01 and type [RET]'
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

      plot_header = trim(subst_name)//'  '//file_ps//'  '//trim(at_weight_scheme(j_weight))//'  '//trim(masks)//'  '//trim(smooth)	

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
          write(plot_title,'(a,a,i1.1," x",2f7.3)') trim(pdf_out(j_mode)),'_ext',j,ext_scale(j),ext_dy(j)
          y_plot = y_plot-.04*(c_max-c_min)
          CALL PGSCI(15)  !grey
          CALL PGSLW(10)			!operates in steps of 5
          CALL PGLINE(n_x,x_ext(1:n_x),ext_scale(j)*(y_ext(1:n_x,j)+ext_dy(j)))  !plots the curve 
          CALL PGSTBG(0)																				 !erase graphics under text
          CALL PGSLW(5)			!operates in steps of 5
          CALL PGTEXT (x_plot,y_plot,plot_title)
        enddo
      endif

      if(tot_scale/=0) then
        y_plot = y_plot-.06*(c_max-c_min)
        CALL PGSCI(1)  !white
        CALL PGLINE(n_plot,x(i_start:i_end),tot_scale*sq_plot_tot(i_start:i_end))  !plots the curve         
        CALL PGSTBG(0)																				 !erase graphics under text
        write(plot_title,'(a,a," x",f7.3)') 'Total','    ',tot_scale
        CALL PGTEXT (x_plot,y_plot,plot_title)
      endif

      do j=1,n_part
        if(ind_part(1,j)==0.or.ind_part(2,j)==0) cycle
        if(part_scale(j)==0) cycle
        write(plot_title,'(a,a," x",f7.3)') at_name_plot(ind_part(1,j)),at_name_plot(ind_part(2,j)),part_scale(j)
        y_plot = y_plot-.04*(c_max-c_min)
        CALL PGSCI (j+1)  !red-green-blue
        if(j_mode==5) then
          CALL PGSLW(5)
        else
          CALL PGSLW(2)
        endif
        CALL PGLINE(n_plot,x(i_start:i_end),part_scale(j)*sq_plot(ind_part(1,j),ind_part(2,j),i_start:i_end))  !plots the curve
        CALL PGSTBG(0)																				 !erase graphics under text
        CALL PGSLW(5)			!operates in steps of 5
        CALL PGTEXT (x_plot,y_plot,plot_title)				!the x_plot,y_plot are in world (axis units) coordinates
        CALL PGSLW(2)
      enddo    

      x_plot = x_start+.8*(x_end-x_start)
      y_plot = c_min-.1*(c_max-c_min)
      CALL PGSCI (1)  !white needs to be reset after PGLAB
      CALL PGSTBG(0)																				 !erase graphics under text
      CALL PGSLW(2)			!operates in steps of 5
      CALL PGSCH(.6)
      CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp))

      CALL PGCLOS
      print *,space, ' Postscript output written to: ',file_ps	
      write(9,*)
      write(9,*) '  ',' S(Q) by a 3D_averaged NUFFT '
      write(9,*) '  ','Masks:',(at_mask(i),i=1,n_atom)				
      if(n_smooth>1) write(9,*) '  ','Smoothing FWHM:',n_smooth_fwhm	
      write(9,*) 'Postscript output written to: ',file_ps

! *** save the S(Q) results into an ASCII file (each line corresponds to a distance point)
!     look for existing output files in order not overwrite them				
      if(j_txt==1) then
        open (4,file=file_res)

    ! *** write the S(Q) file header
        write(4,*) 'Substance:   ',trim(subst_name),'       ',trim(mp_tool)//' '//trim(time_stamp)
        write(4,*) 'Input files: ',trim(file_dat_t0),' to ',trim(file_dat),' step',nfile_step									   
        write(4,*) 'OMP processes  = ', proc_num_in									
        write(4,*) 'Supercell size:',n_row																				
        write(4,*) 'Unit cell parameter:',a_par																				
        write(4,*) 'Atoms/unit_cell:',n_atom
        write(4,*) 'Atoms       :',(('    '//at_name_plot(i)),i=1,n_atom)
        write(4,'(1x,"Atom no.    :",50i8)') (i,i=1,n_atom)
        write(4,'(1x,"Occupancies: ",50f8.4)') (at_occup_r(i),i=1,n_atom)
        write(4,'(1x,"Amplitudes  : ",50f8.4)') (at_scf(i)*at_weight(i),i=1,n_atom)
        write(4,'(1x,"Scale factors:",51(1x,f8.4))') (at_scf(i),i=1,n_atom)                    
        write(4,'(1x,"Masks	      : ",50i8)') (at_mask(i),i=1,n_atom)
        write(4,*)  'Smoothing FWHM [Å]',n_smooth_fwhm*q_step
        write(4,*)
        write(4,*) 'Output: ',pdf_out(j_mode)

    ! *** write the S(Q)_l
        i_start = int(.5/q_step)+n_smooth/2
        if(j_mode==5) then       
          write(4,'(a,50(2i3,11x))') '    Q [A-1]          ',((j,k,k=1,j),j=1,n_atom)
          do i=i_start,n_hist
            write(4,*) q_hist(i),((sq_plot(j,k,i),k=1,j),j=1,n_atom)
          enddo
        else
          write(4,'(a,50(2i3,11x))') '    Q [A-1]          total          ',((j,k,k=1,j),j=1,n_atom)
          do i=i_start,n_hist
            write(4,*) q_hist(i),sq_plot_tot(i),((sq_plot(j,k,i),k=1,j),j=1,n_atom)
          enddo
        endif
        
        close(4)
        print *,space, ' Text output written to: ',file_res	  
        write(9,*) ' Text output written to: ',file_res	  
      endif 																							!j_txt
    endif 																							  !j_ps

! ***   all done, now decide what comes next
  
    way_point: do
      print *                 
      print *,prompt, 'Choose output options (MODE is ',trim(pdf_out(j_mode)),', FILE output is ',trim(ps_out(j_ps+1))       !,', SIZE is ',trim(size_out(j_out+1)),'):'
      print *,space, '       1   REPLOT the last graph'
      write(*,'(10x,"        2   select max ",i2," partial PDFs & replot")') n_part_max
      print *,space, '       3   adjust partials scales & replot '
      print *,space, '       4   modify atom WEIGHTS (',trim(at_weight_scheme(j_weight)) 
      print *,space, '       5   edit atom MASKS ',trim(masks(8:))
      print *,space, '       6   edit atom scale factors ',at_scf
      print *,space, '       7   toggle FILE output ',trim(ps_out(mod(j_ps+1,2)+1)),' (mind the J_TXT switch in .PAR)'
      if(j_ext==0)then
        print *,space, '       8   import external ',trim(pdf_out(j_mode)),' curve '
      else
        print *,space, '       8   change/close external ',trim(pdf_out(j_mode)),' curve ',file_inp
      endif
     
!!!					print *,space, '       8   toggle .TXT output SIZE to ',trim(size_out(mod(j_out+1,2)+1))
      print *,space, '       9   RESTART with updated weights & masks'
      if(product(ind_mo)==1) print *,space, '      10   select the OUTPUT MODE, actual: ',trim(pdf_out(j_mode))
      print *                 
      print *,space, '       0   EXIT'  

      read(*,*) jj
  
      select case(jj)
        case(1) 
          cycle plot_loop

        case(2) 
          print *,prompt, '("Confirm/modify up to ", i2," pairs of partial PDF indices (0 0 erase, -1 -1 skip the rest):")',n_part_max
          do j=1,n_part_max
            j1 = ind_part(1,j)
            j2 = ind_part(2,j) 
            write(*,'(i2,": [",i2,",",i2,"]  ")',advance='no') j,j1,j2
            read(*,*) j1,j2
            if(j1==-1.and.j2==-1) exit
            ind_part(1,j) = j1
            ind_part(2,j) = j2
            if (ind_part(1,j)<ind_part(2,j)) ind_part(:,j) = cshift(ind_part(:,j),1)      !the 2nd one should be smaller  
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
          write(*,"('Confirm/modify total scale factor ',f6.3)") tot_scale
          read(*,*) tot_scale
          write(*,"('Confirm/modify partial scale factors ',4f6.3)") part_scale(1:n_part)
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
            print *,prompt, 'Atom weights ( 1= uniform, 2= neutron b_c^2, 3= Xray f(Q))'	
            read(*,*) j_weight
            tot_scale = 1.
            if(j_weight>=1.and.j_weight<=3) exit
          enddo
          cycle way_point
      
        case(5) 
          write(*,'(" Actual masks:",(50i3))') (at_mask(i),i=1,n_atom)
          print *,prompt,'Type in new ones (0/1):'
          do
            read(*,*)(at_mask(i),i=1,n_atom)
            if(any(at_mask(1:n_atom).ge.0).and.any(at_mask(1:n_atom).le.1)) exit
            print *,space, 'Input out of range, repeat ...'
          enddo
          cycle way_point
          
        case(6) 
          print *,space, ' Actual values:',(at_scf(i),i=1,n_atom)
          print *,prompt, ' Confirm or type in new ones:'
          read(*,*) at_scf
          exit plot_loop

        case(7) 
          j_ps = j_ps+1
          j_ps = mod(j_ps,2)
          cycle way_point

        case(8)
          if(j_ext/=0) deallocate(x_ext,y_ext)
          ext_scale = 1.
          print *,prompt, 'Confirm or modify input file name ("=" close file):'
          read(*,*) file_inp
          if(index(file_inp,'=')/=0) then
            j_ext = 0
            cycle way_point
          else
            open(4,file=trim(file_inp),iostat=ios)
            if(ios.ne.0) print *,space, 'File ',trim(file_inp),' not opened! IOS =',ios
            print *,prompt, 'number of lines to skip, to read, column X, column Y(1)...Y(4) (0 if not used)?'
            read(*,*) n_ext_skip,n_x,j_x,ind_ext

            n_part_ext = 0
            do j=1,n_part_max
              if(ind_ext(j)>0) n_part_ext = n_part_ext+1
            enddo
            if(n_part_ext==0) then
              print *,space, 'Number of external data columns must be >0'
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
            print *,space, n_x,'data points read in'
            j_ext = 1
           endif 
         cycle way_point
             
        case(9) 
          exit plot_loop

        case(10) 
          do
            print *,prompt, 'Select the output MODE: ',(j,'  ',trim(pdf_out(j)),j=3,n_mode),'_unscaled'
            read(*,*) j_mode         
            if(j_mode>=3.and.j_mode<=n_mode) exit plot_loop
          enddo

        case(0) 
          exit bz_loop

      end select

    enddo way_point
  
    enddo plot_loop

  enddo at_weights_loop

  deallocate(x_ffq,q_norm,q_hist)		

  enddo bz_loop

  close(3)
  flush(9)
  close(9)
  CALL PGEND
              
  contains

  function vector_product(e1,e2)  
    real ::   e1(3),e2(3),vector_product(3)
  
      vector_product = (/e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)/) 
  
  end function     
!

!-----------------------------------------------------------------------
! gjinv - Invert a matrix, Gauss-Jordan algorithm
! A is destroyed.
!
!___Name_______Type_______________In/Out____Description_________________
!   A(LDA,N)   Real               In        An N by N matrix
!   LDA        Integer            In        Row bound of A
!   N          Integer            In        Order of matrix
!   B(LDB,N)   Real               Out       Inverse of A
!   LDB        Integer            In        Row bound of B
!   IERR       Integer            Out       0 = no errors
!                                           1 = singular matrix
!-----------------------------------------------------------------------
    SUBROUTINE GJINV (A, LDA, N, B, LDB, IERR)
      IMPLICIT NONE
      INTEGER LDA, N, LDB, IERR
      REAL A(LDA,N), B(LDB,N)

       REAL EPS                                  ! machine constant
       PARAMETER (EPS = 1.1920929E-07)
       INTEGER I, J, K, P                        ! local variables
       REAL F, TOL

!-----------------------------------------------------------------------
!             Begin.
!-----------------------------------------------------------------------
       IF (N < 1) THEN            ! Validate.
         IERR = -1
         RETURN
       ELSE IF (N > LDA .OR. N > LDB) THEN
         IERR = -2
         RETURN
       END IF
       IERR = 0

       F = 0.                     ! Frobenius norm of A
       DO J = 1, N
         DO I = 1, N
           F = F + A(I,J)**2
         END DO
       END DO
       F = SQRT(F)
       TOL = F * EPS

       DO J = 1, N                ! Set B to identity matrix.
         DO I = 1, N
           IF (I .EQ. J) THEN
             B(I,J) = 1.
           ELSE
             B(I,J) = 0.
           END IF
         END DO
       END DO

!             Main loop
       DO K = 1, N
         F = ABS(A(K,K))          ! Find pivot.
         P = K
         DO I = K+1, N
           IF (ABS(A(I,K)) > F) THEN
             F = ABS(A(I,K))
             P = I
           END IF
         END DO

         IF (F < TOL) THEN        ! Matrix is singular.
           IERR = 1
           RETURN
         END IF

         IF (P .NE. K) THEN       ! Swap rows.
           DO J = K, N
             F = A(K,J)
             A(K,J) = A(P,J)
             A(P,J) = F
           END DO
           DO J = 1, N
             F = B(K,J)
             B(K,J) = B(P,J)
             B(P,J) = F
           END DO
         END IF

         F = 1. / A(K,K)          ! Scale row so pivot is 1.
         DO J = K, N
           A(K,J) = A(K,J) * F
         END DO
         DO J = 1, N
           B(K,J) = B(K,J) * F
         END DO

         DO 10 I = 1, N           ! Subtract to get zeros.
           IF (I .EQ. K) GO TO 10
           F = A(I,K)
           DO J = K, N
             A(I,J) = A(I,J) - A(K,J) * F
           END DO
           DO J = 1, N
             B(I,J) = B(I,J) - B(K,J) * F
           END DO
  10     CONTINUE
       END DO

       RETURN
      END  ! of gjinv



  SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set the RAINBOW palette of colors to be used by PGIMAG.
!-----------------------------------------------------------------------
    INTEGER TYPE
    REAL CONTRA, BRIGHT
!
    REAL GL(2), GR(2), GG(2), GB(2)
    REAL RL(9), RR(9), RG(9), RB(9)
    REAL JKL(9), JKR(9), JKG(9), JKB(9)
    REAL HL(5), HR(5), HG(5), HB(5)
    REAL WL(10), WR(10), WG(10), WB(10)
    REAL AL(20), AR(20), AG(20), AB(20)
!
    DATA GL /0.0, 1.0/
    DATA GR /0.0, 1.0/
    DATA GG /0.0, 1.0/
    DATA GB /0.0, 1.0/
!
    DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
    DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
    DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
    DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
    DATA JKL / 0,   0.1,  0.2,  0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
    DATA JKR / 1.0, 0.5,  0.1,  0.1,  0.3,  1.0,  1.0, 1.0, 1.0/
    DATA JKG / 1.0, 0.7,  0.4,   .8,  1.,   1.0,  0.6, 0.0, 1.0/
    DATA JKB / 1.0, 1.0,  0.9,   .6,  0.3,  0.0,  0.0, 0.0, 1.0/
!
    DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
    DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
    DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
    DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
    DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
    DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
    DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
    DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
    DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,&
   &         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
    DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,&
   &         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
    DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,&
   &         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
    DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,&
   &         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
    IF (TYPE.EQ.1) THEN
!        -- gray scale
       CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
    ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
       CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
    ELSE IF (TYPE.EQ.3) THEN
!        -- heat
       CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
    ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
       CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
    ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
       CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
    ELSE IF (TYPE.EQ.6) THEN
!        -- JK rainbow (printer friendly)
       CALL PGCTAB(JKL, JKR, JKG, JKB, 9, CONTRA, BRIGHT)
    END IF
  END

end program mp_sql56


!! Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
!! Contact: greengard@cims.nyu.edu
!! 
!! This software is being released under a FreeBSD license
!! (see license.txt in this directory). 
!!
subroutine dirft2d1(nj,xj,yj,cj, iflag, ms,mt,fk)
  implicit none
  integer iflag
  integer(8) nj, ms, mt
  real(8) xj(nj), yj(nj)
  complex*16 cj(nj), fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2)
!     ------------------------------------------------------------------
!     direct computation of nonuniform FFT
!
!                   nj
!     fk(k1,k2) =  SUM cj(j) exp(+/-i k1 xj(j)) exp(+/-i k2 yj(j))
!                  j=1
!
!     for -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2
!
!     If (iflag .ge.0) the + sign is used in the exponential.
!     If (iflag .lt.0) the - sign is used in the exponential.
!
!***********************************************************************
  integer j, k1, k2
  complex*16 zf, cm1, z1n(-ms/2:(ms-1)/2)
!
  do k2 = -mt/2, (mt-1)/2
     do k1 = -ms/2, (ms-1)/2
        fk(k1,k2) = dcmplx(0d0,0d0)
     enddo
  enddo
!
  do j = 1, nj
!
!     ----------------------------------------------------------
!     Precompute exponential for exp(+/-i k1 xj)
!     ----------------------------------------------------------
!
     if (iflag .ge. 0) then
        zf = dcmplx(dcos(xj(j)),+dsin(xj(j)))
     else
        zf = dcmplx(dcos(xj(j)),-dsin(xj(j)))
     endif
     z1n(0) = (1d0,0d0)
     do k1 = 1, (ms-1)/2
        z1n(k1) = zf*z1n(k1-1)
        z1n(-k1)= dconjg(z1n(k1))
     enddo
     if (ms/2*2.eq.ms) z1n(-ms/2) = dconjg(zf*z1n(ms/2-1))
!
!     ----------------------------------------------------------
!     Loop over k2 for yj
!     ----------------------------------------------------------
     if (iflag .ge. 0) then
        zf = dcmplx(dcos(yj(j)),+dsin(yj(j)))
     else
        zf = dcmplx(dcos(yj(j)),-dsin(yj(j)))
     endif
!
     cm1 = cj(j)
     do k2 = 0, (mt-1)/2
        do k1 = -ms/2, (ms-1)/2
          fk(k1,k2) = fk(k1,k2) + cm1*z1n(k1)
        enddo
        cm1 = cm1*zf
     enddo
!
     zf = dconjg(zf)
     cm1 = cj(j)
     do k2 = -1, -mt/2, -1
        cm1 = cm1*zf
        do k1 = -ms/2, (ms-1)/2
          fk(k1,k2) = fk(k1,k2) + cm1*z1n(k1)
        enddo
     enddo
  enddo
end
!
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
 
