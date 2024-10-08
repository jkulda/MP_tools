  
program mp_dplot56

! *************************************************************************************
! *****
! ***** %%%%%%%%%%%%%%%%         program MP_DPLOT  1.56        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! *****
! *****   plots diverse real-space variables from simulated data files
! *****
!**********          Copyright (C) 2023  Jiri Kulda, Grenoble/Prague          **********
!**  
!** This file is part MP_TOOLS developed and maintained by Jiri Kulda <jkulda@free.fr>
!**
!**  MP_TOOLS are free software: you can use it, redistribute it and/or modify it 
!**  under the terms of the GNU General Public License as published by the Free Software 
!**  Foundation, either version 3 of the License, or (at your option) any later version, 
!**  for details see <https://www.gnu.org/licenses/>
!**
!**  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!**  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!**  See the GNU General Public License for more details.
!**
! ***** Ver. 1.54 - NAMELIST I/O for .PAR files and .DAT headers implemented 
! *****           - PNG driver for PGPLOT included
! *****           
! ***** Ver. 1.55 - bug fixes 
! ***** 
! ***** Ver. 1.56 - bug fixes 
! ***** 
! ***** file name convention:
! *****     <master_file>_n<snapshot number>.dat
! *****     example: H05_10K_n0001.dat from H05_10K.txt
! *****
! ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
! ***** record numbers are directly related to cell positions and atom types
! *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
! *****            1 going for the header record
! *****
! ***** atom positions are converted from reduced lattice coordinates (x) to real space distances
! ***** 

  integer,parameter :: l_rec  =  1024        !record length in real(4)
  real,parameter    :: pi = 3.14159
  real,parameter    :: e_C = 16.02176634       ! e/Å^2 to C/m^2
  character(4),allocatable  ::  at_name_par(:),at_label(:),file_ext
  character(10),allocatable ::  dom_ind(:)

  integer,allocatable ::  i_dom(:,:,:,:,:),i_dom_p(:,:,:,:),i_dom_out(:,:,:),ind(:,:,:),mask(:)

  real,allocatable :: displ_field(:,:,:,:,:,:),displ_field1(:,:,:,:,:,:),displ_norm(:,:,:,:,:),displ_norm_tot(:,:)
  real,allocatable :: vel_field(:,:,:,:,:,:),vel_norm(:,:,:,:,:),vel_field1(:,:,:,:,:,:),vel_norm1(:,:,:,:,:)
  real,allocatable :: polar(:,:,:,:,:),polar_tot(:,:),polar_tot_abs(:,:),pol_norm(:,:,:,:),pol_norm_tot(:),pol_norm_tot2(:)
  real,allocatable :: at_base(:,:),at_charge(:),at_charge1(:),displ_plot(:,:,:,:),displ_vect(:,:,:)

  real     ::  c_min,c_max,c_min_save,c_max_save,filter_fwhm,p_size,res(3),pol(3),res2,charge_mom_c(3),charge_mom_s(3),charge_mom_cell(3)
  real     ::  e1_norm(3),e2_norm(3),ev_norm(3),ed_norm(3),e1p(3),e2p(3),evp(3),e1p_norm,e2p_norm,evp_norm,ep_angle,x1,y1,x2,y2,x_plot,y_plot
  
  character(4)   :: c_int(2),c_fil(2),version,head,atom
  character(10)	 :: prompt,space = '          ',cs_string(3)
  character(10)  :: pg_out,string,section,c_date,c_time,c_zone,c_jt,c_slice,c_mode,c_jfile,at_name,dom_name,polar_name
  character(16)  :: sim_type_par,data_type,string16,wedge_label,filter_name,c_e1(3),c_e2(3),c_x,c_y,header_1
  character(40)  :: subst_name,file_master,file_inp,file_out,time_stamp,int_mode,x_file_name,mp_tool,z_str,ev_str
  character(60)  :: file_dat,file_dat_t0,file_res,file_ps,file_log,line,header_2
  character(128) :: cwd_path,plot_title
  character(l_rec):: header_record
  
  logical ::  nml_in,found_txt,found_ps,t_single,mode_s,mode_pix,mode_pts,mode_vect

  integer(8) ::  n_in,n_out,n_inbox,n_outbox,n_pdf_range(3),n_mc,n_mc_max,n_norm,j_ox
  integer ::  j_proc,proc_num,proc_num_in,j_name,jm,j1m,j2m,j_mode,mode
  integer ::  i_start,i_end,j_part(4),n_step,nfile_step,n_hsum,n_h,j_head_in,j_verb,n_tot,n_head,j_auto,j_frame,n_slice,j_atom,n_hist
  integer ::  i,j,k,m,i2,i3,ii,jj,at,i_rec,ind_rec,nrec,ncell,nrow,nlayer,nsuper,nfile,nfile_t0,nfile_min,nfile_max,jfile
  integer ::  j1,j2,j_plane,j_grid,j_logsc,j_ps,j_op,j_tot,j_txt,i_time(8),j_cycle,i_shift,j_adv,n_shift,j_shift,j_slice,j_sign,jt,jt0,jt_max
  integer ::  at_no,at_ind(3),at_ind2(3),d_ind(3),n_dom,j_base(3),n_plot,e1(3),e2(3),ev(3),e_slice(3),ind_c,ind_h(3),n_x,n_y
  integer ::  jat,j_edit,j_mask,ifile,ifile0,jfil_step,jint,jint2,n_atom,n_corr,sc_c2,sc_c1,ier,ios,j_pol,j_cs
  integer ::  j_weight,j_xray
  
! **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
  character(4),allocatable :: at_name_out(:),corr_name(:)
  integer(4),allocatable   :: nsuper_r(:),at_ind_in(:),at_ind_file(:,:,:,:,:)

  real(4),allocatable ::  at_occup_r(:),at_pos_hist(:,:,:,:,:,:),hist_plot(:,:,:),at_pos_ref(:,:,:,:,:),at_pos_ref1(:,:,:,:,:)
  real(4),allocatable,target ::  at_pos_in(:),at_pos1_in(:),at_vel_in(:),at_vel1_in(:)
  real,pointer        :: at_pos_file(:,:,:,:,:),at_pos1_file(:,:,:,:,:),at_vel_file(:,:,:,:,:),at_vel1_file(:,:,:,:,:)
  
  character(16)  :: sim_type,dat_type,input_method,file_par,dat_source
  integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,n_tot_in,idum,j_out,j_pgc,j_shell_plot
  real(4)        :: rec_zero(l_rec),t_ms,t0,t_dump,t_step,a_par(3),angle(3),a_cell(3,3),a_cell_inv(3,3),temp,temp_cs,hist_step(3)

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,temp_cs,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
  namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
  namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
 
  namelist /mp_gen/ j_verb,j_proc       
  namelist /mp_out/ j_weight,j_logsc,j_txt,p_size,j_grid,pg_out,j_ps,j_out,j_pgc        
                  !general rule: namelists of tools should only contain their local parameters
                      !what is of global interest they should pass into data_header
  
! *** PGPLOT variables 
  INTEGER   PGBEG,PGOPEN,map_unit(3)
  INTEGER   MXI, MXJ
  PARAMETER (MXI=21, MXJ=21)
  INTEGER   L, C1, C2,CI, NC
  REAL      F(MXI,MXJ),CR,CG,CB,CH,CL,CS
  REAL      FMIN,FMAX,TR(6), CONTRA, BRIGHT, C, S, ALEV(1),BLANK,SCALE,SIZE
  CHARACTER*16 VAL 
  CHARACTER*4 XOPT, YOPT
  REAL XTICK, YTICK
  INTEGER NXSUB, NYSUB					


! ********************* Initialization *******************************      
  version = '1.56'
  prompt = 'MP_DPLOT> '
  mp_tool = 'MP_DPLOT '//version

  print *,'*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  print *,'     *** for the moment only orthogonal output is available *** '
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
  write(9,*) trim(time_stamp),'  ',mp_tool,'  ',trim(cwd_path) 
  write(9,*) 
  
! *** Generate data file access
  print *,prompt, 'Data file_master & file numbers (min, max; 0 0 single file): '
  read(*,*) file_master,nfile_min,nfile_max
  nfile_step = 1
  
  print *,prompt, 'Display shells instead of cores? (1/0) '
  read(*,*) j_shell_plot
  
        
  t_single = nfile_min==0.and.nfile_max==0  
  j_pgc = 6       !PG_PLOT palette JK
  map_unit = 0
  n_hist = 101
  hist_step = .0025  
  j_ox = 3     !index of the 1st oxygen atom (assuming they are consecutive)
  polar_name = 'POLAR'
  n_corr = 6
  j_cs = 1
  cs_string = (/'cores ','shells','CSdist'/)
  mode_s = .false.
  mode_pix = .true.
  mode_pts = .false.
  mode_vect = .false.
  
  if(t_single)then
    nfile_min = 1    !for conformity with loop control conventions
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
  read(1,rec=i_rec) header_record
  head = header_record(1:4)
  call up_case(head)
  if(head == 'MP_T') then      !new structure with namelist
    nml_in = .true.
    read(1,rec=i_rec) dat_source,version,string16
    read(string16,*) n_head
    print *,space,    'Input data:  ',dat_source,version,n_head
    i_rec = i_rec+1                 
    read(1,rec=i_rec) header_record
    read(header_record,nml=data_header_1)  
    t0 = t_dump
!!        write(*,nml=data_header_1)
  else
    print *,space, 'header record wrong or missing'
    print *,space, trim(header_record)
    stop
  endif 
  
  if(trim(input_method)/='CELL') then
    print *,space, 'ERROR: the present data were produced by the ',trim(input_method),' method,'
    print *,space, 'while direct space mapping works with the CELL data only.'
    stop
  endif

  allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom))
  allocate(at_label(n_atom),at_name_par(n_atom),at_base(n_atom,3),at_charge(n_atom),at_charge1(n_atom),corr_name(n_corr-n_atom))
  


  i_rec = i_rec+1                 
  read(1,rec=i_rec) header_record
  read(header_record,nml=data_header_2)
  i_rec = i_rec+1                 
  read(1,rec=i_rec) header_record
  read(header_record,nml=data_header_3)
  close(1)

! **** Read the auxiliary file <file_par.par> with structure parameters, atom names and further info
!     print *,'Parameter file name (.par to be added) (confirm or type other name): ', file_par
!     read(*,*) file_par

  file_inp = trim(file_par)//'.par'
  print *,space, 'Parameter file name: ', file_inp

  open(4,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,space, 'File ',trim(file_inp),' not found! Stop execution.'
    stop
  endif

  write(9,*) 'Read the parameter file:  ',trim(file_inp)

  read(4,nml=mp_gen)
  rewind(4)
  read(4,nml=mp_out)
  rewind(4)

  j_adv = 2
  n_shift = 10
  at_name_par = at_name_out  
  close(4)

! *** write overview of atom data
  print *
  print *,space, 'Substance name: ',subst_name    
  print *,space, 'Atoms from ',trim(file_inp)
  do j=1,n_atom
    write(*,'(15x,i2,5x,a4,3f8.4,2x,f8.4)')  j,at_name_par(j),at_base(j,:),at_occup_r(j)
    write(9,'(5x,a4,3f8.4,2x,f8.4)')  at_name_par(j),at_base(j,:),at_occup_r(j)
  enddo
  
  corr_name(1) = trim(at_name_par(1))//trim(at_name_par(2))
  
!   print *,space,'Translation vector to correct at_base origin: (0 0 0 skip)'
!   read (*,*) res
!   do j=1,n_atom
!     at_base(j,:) = at_base(j,:)+res
!   enddo
!   
!   if(sum(abs(res))/=.0) then
!     do j=1,n_atom
!       write(*,*) 'Corrected values:'
!       write(9,*) 'Corrected values:'
!       write(*,'(15x,i2,5x,a4,3f8.4,2x,f8.4)')  j,at_name_par(j),at_base(j,:),at_occup_r(j)
!       write(9,'(5x,a4,3f8.4,2x,f8.4)')  at_name_par(j),at_base(j,:),at_occup_r(j)
!     enddo
!   endif

! *** define references for lattice search  
!     nrow = n_row(1)
  nsuper = n_row(1)*n_row(2)*n_row(3)


! *** Allocate and clear the histogram arrays for accumulation across several snapshots         
!     allocate (at_ind_file(4,n_row(1),n_row(2),n_row(3),jat),at_ind_in(4*n_tot))

  allocate(at_pos_in(4*n_tot),at_vel_in(4*n_tot))
  allocate(displ_field(3,n_row(1),n_row(2),n_row(3),n_atom,nfile),displ_norm(n_row(1),n_row(2),n_row(3),n_corr,nfile),displ_norm_tot(n_corr,nfile))
  allocate(at_pos_ref(3,n_row(1),n_row(2),n_row(3),n_corr))
  allocate(vel_field(3,n_row(1),n_row(2),n_row(3),n_atom,nfile),vel_norm(n_row(1),n_row(2),n_row(3),n_atom,nfile))
  allocate(i_dom(n_row(1),n_row(2),n_row(3),n_atom,nfile),i_dom_p(n_row(1),n_row(2),n_row(3),nfile))
  allocate(polar(3,n_row(1),n_row(2),n_row(3),nfile),pol_norm(n_row(1),n_row(2),n_row(3),nfile),polar_tot(3,nfile),polar_tot_abs(3,nfile),pol_norm_tot(nfile),pol_norm_tot2(nfile))
  if(j_shell_out==1) then 
    allocate(at_pos1_in(4*n_tot),at_vel1_in(4*n_tot),displ_field1(3,n_row(1),n_row(2),n_row(3),n_atom,nfile),at_pos_ref1(3,n_row(1),n_row(2),n_row(3),n_corr))
    allocate(vel_field1(3,n_row(1),n_row(2),n_row(3),n_atom,nfile),vel_norm1(n_row(1),n_row(2),n_row(3),n_atom,nfile))
  endif
  allocate (at_pos_hist(n_hist,n_hist,n_hist,n_corr,1+2*j_shell_out,nfile))

  print *,prompt, 'Displacement domain type (0 = NONE, 1=[100], 2=[110], 3= 111])?'
  read(*,*) n_dom		!displacement domain type (1=[100], 2=[110], 3= 111])
  if(n_dom==0) i_dom = 1
  
  print *,space, '          (reading input files ...)'
  print *
 			
!!! *** calculate the domain segment number - new convention different from MD_TOOLS
!
!			i_dom = 0																		!if n_dom = 0 nothing happens
!      do ii = 1,3
!      	pos_inp(ii) = at_pos_c(i,ii)-x_pos(jat,ii)
!      	pos_inp(ii) = pos_inp(ii)-nint(pos_inp(ii))
!      enddo
!			
!			if(n_dom.eq.1) then
!				i_dom = maxloc((abs(pos_inp)),dim=1)
!!					print *,'maxloc',i_dom
!				if(pos_inp(i_dom).gt.0.) i_dom = i_dom+3
!			else if (n_dom.eq.2) then
!				ii = minloc((abs(pos_inp)),dim=1)
!!					print *,'minloc',ii
!				i2 = ii+1
!				if(i2.gt.3) i2 = i2-3
!				i3 = ii+2
!				if(i3.gt.3) i3= i3-3
!				i_dom = mod(ii,3)*4+(sign(1.,pos_inp(i2))+1)+(sign(1.,pos_inp(i3))+1)/2 +1
!			else if (n_dom.eq.3) then
!				i_dom = 1
!				do ii=1,3
!					i_dom = i_dom+2**(ii-1)*(sign(1.0,pos_inp(ii))+1.)/2.
!				enddo 
!			endif

  jt0 = 1
  jt = 0
  at_pos_hist = .0

  file_loop: do ifile=nfile_min,nfile_max,nfile_step                  
    jt = jt+1
! ***  open the t0 file (binary MD snapshot file)
    if(t_single)then
      write(file_dat,'(a,".dat")') trim(file_master)
    else
      if(nfile_min<=9999) then
        write(file_dat,'(a,"_n",i4.4,".dat")') trim(file_master),nfile_min+(ifile-nfile_min)*nfile_step
      elseif(nfile_min>=10000) then
        write(string,'(i8)') nfile_min+(ifile-nfile_min)*nfile_step
        file_dat = trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
      endif
    endif

!       print *
!       print *'input: ',file_dat
    write(9,*)'input: ',file_dat

    open(1,file='./data/'//file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
    if(ios.ne.0) then
      print *,space, 'File ',trim(file_dat),' not opened! IOS =',ios
      print *,prompt, 'Skip(1), stop execution(0)?'
      read(*,*) jj
      if(jj==1) exit file_loop
      if(jj==0) stop
    endif

    i_rec = n_head+n_rec      !skipping at_ind
  
    do j=1,n_rec-1                    ! read CORE data
      i_rec = i_rec+1
      read(1,rec=i_rec) at_pos_in((j-1)*l_rec+1:j*l_rec)      
    enddo  
    i_rec = i_rec+1
    read(1,rec=i_rec) at_pos_in((n_rec-1)*l_rec+1:4*n_tot)  
    
    if(n_traj==1) then
      do j=1,n_rec-1                    ! read CORE data
        i_rec = i_rec+1
        read(1,rec=i_rec) at_vel_in((j-1)*l_rec+1:j*l_rec)      
      enddo  
      i_rec = i_rec+1
      read(1,rec=i_rec) at_vel_in((n_rec-1)*l_rec+1:4*n_tot) 
    endif 
    
    if(j_shell_out==1) then           ! read SHELL data
      i_rec = i_rec+(n_traj-1)*n_rec      !skip possible velocities & forces
      do j=1,n_rec-1
        i_rec = i_rec+1
        read(1,rec=i_rec) at_pos1_in((j-1)*l_rec+1:j*l_rec)      
      enddo  
      i_rec = i_rec+1
      read(1,rec=i_rec) at_pos1_in((n_rec-1)*l_rec+1:4*n_tot)  

      if(n_traj==1) then
        do j=1,n_rec-1                    ! read CORE data
          i_rec = i_rec+1
          read(1,rec=i_rec) at_vel1_in((j-1)*l_rec+1:j*l_rec)      
        enddo  
        i_rec = i_rec+1
        read(1,rec=i_rec) at_vel1_in((n_rec-1)*l_rec+1:4*n_tot) 
      endif 
    endif        
    
    close(1)

    if(j_shell_plot==0)then
      at_pos_file(1:4,1:n_row(1),1:n_row(2),1:n_row(3),1:n_atom) => at_pos_in(:)
      at_vel_file(1:4,1:n_row(1),1:n_row(2),1:n_row(3),1:n_atom) => at_vel_in(:)
      header_1 = 'cores'
    else
      at_pos_file(1:4,1:n_row(1),1:n_row(2),1:n_row(3),1:n_atom) => at_pos1_in(:)
      at_vel_file(1:4,1:n_row(1),1:n_row(2),1:n_row(3),1:n_atom) => at_vel1_in(:)
      header_1 = 'shells'
    endif


    vel_field(:,:,:,:,:,jt) = at_vel_file(1:3,:,:,:,:)
    do jj = 1,n_atom
      vel_norm(:,:,:,jj,jt) = sqrt(at_vel_file(1,:,:,:,jj)**2+at_vel_file(2,:,:,:,jj)**2+at_vel_file(3,:,:,:,jj)**2)
    enddo

    if(j_shell_out==1) then
      at_pos1_file(1:4,1:n_row(1),1:n_row(2),1:n_row(3),1:n_atom) => at_pos1_in(:)
      at_vel1_file(1:4,1:n_row(1),1:n_row(2),1:n_row(3),1:n_atom) => at_vel1_in(:)
      vel_field1(:,:,:,:,:,jt) = at_vel1_file(1:3,:,:,:,:)              
      do jj = 1,n_atom
        vel_norm1(:,:,:,jj,jt) = sqrt(at_vel1_file(1,:,:,:,jj)**2+at_vel1_file(2,:,:,:,jj)**2+at_vel1_file(3,:,:,:,jj)**2)
      enddo
    endif

    i3 = 0
    do jat=1,n_atom     !correct shells torn away from cores and bring back down atoms shifted to top
      do k=1,n_row(3)
        do j=1,n_row(2)
          do i=1,n_row(1)
            if(at_pos_file(1,1,j,k,jat)>.0) at_pos_file(1,1,j,k,jat) = at_pos_file(1,1,j,k,jat)-n_row(1)
            if(at_pos_file(2,i,1,k,jat)>.0) at_pos_file(2,i,1,k,jat) = at_pos_file(2,i,1,k,jat)-n_row(2)
            if(at_pos_file(3,i,j,1,jat)>.0) at_pos_file(3,i,j,1,jat) = at_pos_file(3,i,j,1,jat)-n_row(3)

            if(j_shell_out==1) then
              if(at_pos1_file(1,1,j,k,jat)>.0) at_pos1_file(1,1,j,k,jat) = at_pos1_file(1,1,j,k,jat)-n_row(1)
              if(at_pos1_file(2,i,1,k,jat)>.0) at_pos1_file(2,i,1,k,jat) = at_pos1_file(2,i,1,k,jat)-n_row(2)
              if(at_pos1_file(3,i,j,1,jat)>.0) at_pos1_file(3,i,j,1,jat) = at_pos1_file(3,i,j,1,jat)-n_row(3)
            endif
          enddo
        enddo
      enddo
    enddo

! ***** calculate reference positions for NN distances
! only valid for ABO3 cubic/orthorhombic perovskites
! all lattice positions are defined by neighbouring oxygen atoms
! cycling through the three oxygens and shifting their indexes provides all the necessary input

   at_pos_ref = .0

! *** the A-atom
   at_pos_ref(:,:,:,:,1) = at_pos_ref(:,:,:,:,1)+at_pos_file(1:3,:,:,:,j_ox)+at_pos_file(1:3,:,:,:,j_ox+1)+at_pos_file(1:3,:,:,:,j_ox+2)
   at_pos_ref(:,:,:,:,1) = at_pos_ref(:,:,:,:,1)+cshift(at_pos_file(1:3,:,:,:,j_ox),-1,2)+cshift(at_pos_file(1:3,:,:,:,j_ox),-1,3)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox),-1,2),-1,3)
   at_pos_ref(:,:,:,:,1) = at_pos_ref(:,:,:,:,1)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),-1,3)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),-1,4)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox+1),-1,3),-1,4)
   at_pos_ref(:,:,:,:,1) = at_pos_ref(:,:,:,:,1)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),-1,2)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),-1,4)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox+2),-1,2),-1,4)
   at_pos_ref(1,1,:,:,1) = at_pos_ref(1,1,:,:,1)-4*n_row(1)
   at_pos_ref(2,:,1,:,1) = at_pos_ref(2,:,1,:,1)-4*n_row(2)
   at_pos_ref(3,:,:,1,1) = at_pos_ref(3,:,:,1,1)-4*n_row(3)
   at_pos_ref(:,:,:,:,1) = at_pos_ref(:,:,:,:,1)/12.

! *** the B-atom
   at_pos_ref(:,:,:,:,2) = at_pos_ref(:,:,:,:,2)+at_pos_file(1:3,:,:,:,j_ox)+at_pos_file(1:3,:,:,:,j_ox+1)+at_pos_file(1:3,:,:,:,j_ox+2)
   at_pos_ref(:,:,:,:,2) = at_pos_ref(:,:,:,:,2)+cshift(at_pos_file(1:3,:,:,:,j_ox),1,4)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),1,2)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),1,3)
   at_pos_ref(1,n_row(1),:,:,2) = at_pos_ref(1,n_row(1),:,:,2)+n_row(1)
   at_pos_ref(2,:,n_row(2),:,2) = at_pos_ref(2,:,n_row(2),:,2)+n_row(2)
   at_pos_ref(3,:,:,n_row(3),2) = at_pos_ref(3,:,:,n_row(3),2)+n_row(3)
   at_pos_ref(:,:,:,:,2) = at_pos_ref(:,:,:,:,2)/6.


! *** the O1 [.5 .5 0] atom
   at_pos_ref(:,:,:,:,3) = at_pos_ref(:,:,:,:,3)+at_pos_file(1:3,:,:,:,j_ox+1)+at_pos_file(1:3,:,:,:,j_ox+2)
   at_pos_ref(:,:,:,:,3) = at_pos_ref(:,:,:,:,3)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),1,2)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),1,3)
   at_pos_ref(:,:,:,:,3) = at_pos_ref(:,:,:,:,3)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),-1,4)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),-1,4)
   at_pos_ref(:,:,:,:,3) = at_pos_ref(:,:,:,:,3)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox+1),1,2),-1,4)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox+2),1,3),-1,4)
   at_pos_ref(1,n_row(1),:,:,3) = at_pos_ref(1,n_row(1),:,:,3)+2*n_row(1)
   at_pos_ref(2,:,n_row(2),:,3) = at_pos_ref(2,:,n_row(2),:,3)+2*n_row(2)
   at_pos_ref(3,:,:,1,3) = at_pos_ref(3,:,:,1,3)-4*n_row(3)
   at_pos_ref(:,:,:,:,3) = at_pos_ref(:,:,:,:,3)/8.

! *** the O2 [0 .5 .5] atom
   at_pos_ref(:,:,:,:,4) = at_pos_ref(:,:,:,:,4)+at_pos_file(1:3,:,:,:,j_ox)+at_pos_file(1:3,:,:,:,j_ox+2)
   at_pos_ref(:,:,:,:,4) = at_pos_ref(:,:,:,:,4)+cshift(at_pos_file(1:3,:,:,:,j_ox),1,4)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),1,3)
   at_pos_ref(:,:,:,:,4) = at_pos_ref(:,:,:,:,4)+cshift(at_pos_file(1:3,:,:,:,j_ox),-1,2)+cshift(at_pos_file(1:3,:,:,:,j_ox+2),-1,2)
   at_pos_ref(:,:,:,:,4) = at_pos_ref(:,:,:,:,4)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox),1,4),-1,2)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox+2),1,3),-1,2)
   at_pos_ref(3,:,:,n_row(3),4) = at_pos_ref(3,:,:,n_row(3),4)+2*n_row(3)
   at_pos_ref(2,:,n_row(2),:,4) = at_pos_ref(2,:,n_row(2),:,4)+2*n_row(2)
   at_pos_ref(1,1,:,:,4) = at_pos_ref(1,1,:,:,4)-4*n_row(1)
   at_pos_ref(:,:,:,:,4) = at_pos_ref(:,:,:,:,4)/8.


! *** the O3 [.5 0 .5] atom
   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)+at_pos_file(1:3,:,:,:,j_ox)+at_pos_file(1:3,:,:,:,j_ox+1)
   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)+cshift(at_pos_file(1:3,:,:,:,j_ox),1,4)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),1,2)
   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)+cshift(at_pos_file(1:3,:,:,:,j_ox),-1,3)+cshift(at_pos_file(1:3,:,:,:,j_ox+1),-1,3)
   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox),1,4),-1,3)+cshift(cshift(at_pos_file(1:3,:,:,:,j_ox+1),1,2),-1,3)
   at_pos_ref(1,n_row(1),:,:,5) = at_pos_ref(1,n_row(1),:,:,5)+2*n_row(1)
   at_pos_ref(3,:,:,n_row(3),5) = at_pos_ref(3,:,:,n_row(3),5)+2*n_row(3)
   at_pos_ref(2,:,1,:,5) = at_pos_ref(2,:,1,:,5)-4*n_row(2)
   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)/8.

!!! *** the B-atom from B-neighbours
!!   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)+cshift(at_pos_file(1:3,:,:,:,2),1,4)+cshift(at_pos_file(1:3,:,:,:,2),1,2)+cshift(at_pos_file(1:3,:,:,:,2),1,3)
!!   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)+cshift(at_pos_file(1:3,:,:,:,2),-1,4)+cshift(at_pos_file(1:3,:,:,:,2),-1,2)+cshift(at_pos_file(1:3,:,:,:,2),-1,3)
!!!    at_pos_ref(1,n_row(1),:,:,2) = at_pos_ref(1,n_row(1),:,:,2)+n_row(1)
!!!    at_pos_ref(2,:,n_row(2),:,2) = at_pos_ref(2,:,n_row(2),:,2)+n_row(2)
!!!    at_pos_ref(3,:,:,n_row(3),2) = at_pos_ref(3,:,:,n_row(3),2)+n_row(3)
!!   at_pos_ref(:,:,:,:,5) = at_pos_ref(:,:,:,:,5)/6.

! *** the A-atom from B-atoms (or vice-versa)
   at_pos_ref(:,:,:,:,6) = at_pos_ref(:,:,:,:,6)+at_pos_file(1:3,:,:,:,2)+cshift(at_pos_file(1:3,:,:,:,2),-1,2)+cshift(at_pos_file(1:3,:,:,:,2),-1,3)+cshift(cshift(at_pos_file(1:3,:,:,:,2),-1,2),-1,3)
   at_pos_ref(:,:,:,:,6) = at_pos_ref(:,:,:,:,6)+cshift(at_pos_file(1:3,:,:,:,2),-1,4)+cshift(cshift(at_pos_file(1:3,:,:,:,2),-1,2),-1,4)
   at_pos_ref(:,:,:,:,6) = at_pos_ref(:,:,:,:,6)+cshift(cshift(at_pos_file(1:3,:,:,:,2),-1,3),-1,4)+cshift(cshift(cshift(at_pos_file(1:3,:,:,:,2),-1,2),-1,3),-1,4)
   at_pos_ref(1,1,:,:,6) = at_pos_ref(1,1,:,:,6)-4*n_row(1)
   at_pos_ref(2,:,1,:,6) = at_pos_ref(2,:,1,:,6)-4*n_row(2)
   at_pos_ref(3,:,:,1,6) = at_pos_ref(3,:,:,1,6)-4*n_row(3)
   at_pos_ref(:,:,:,:,6) = at_pos_ref(:,:,:,:,6)/8.

   if(j_shell_out==1) then
     at_pos_ref1 = .0

! *** the A-atom
    at_pos_ref1(:,:,:,:,1) = at_pos_ref1(:,:,:,:,1)+at_pos1_file(1:3,:,:,:,j_ox)+at_pos1_file(1:3,:,:,:,j_ox+1)+at_pos1_file(1:3,:,:,:,j_ox+2)
    at_pos_ref1(:,:,:,:,1) = at_pos_ref1(:,:,:,:,1)+cshift(at_pos1_file(1:3,:,:,:,j_ox),-1,2)+cshift(at_pos1_file(1:3,:,:,:,j_ox),-1,3)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox),-1,2),-1,3)
    at_pos_ref1(:,:,:,:,1) = at_pos_ref1(:,:,:,:,1)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),-1,3)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),-1,4)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox+1),-1,3),-1,4)
    at_pos_ref1(:,:,:,:,1) = at_pos_ref1(:,:,:,:,1)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),-1,2)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),-1,4)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox+2),-1,2),-1,4)
    at_pos_ref1(1,1,:,:,1) = at_pos_ref1(1,1,:,:,1)-4*n_row(1)
    at_pos_ref1(2,:,1,:,1) = at_pos_ref1(2,:,1,:,1)-4*n_row(2)
    at_pos_ref1(3,:,:,1,1) = at_pos_ref1(3,:,:,1,1)-4*n_row(3)
    at_pos_ref1(:,:,:,:,1) = at_pos_ref1(:,:,:,:,1)/12.

! *** the B-atom
     at_pos_ref1(:,:,:,:,2) = at_pos_ref1(:,:,:,:,2)+at_pos1_file(1:3,:,:,:,j_ox)+at_pos1_file(1:3,:,:,:,j_ox+1)+at_pos1_file(1:3,:,:,:,j_ox+2)
     at_pos_ref1(:,:,:,:,2) = at_pos_ref1(:,:,:,:,2)+cshift(at_pos1_file(1:3,:,:,:,j_ox),1,4)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),1,2)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),1,3)
     at_pos_ref1(1,n_row(1),:,:,2) = at_pos_ref1(1,n_row(1),:,:,2)+n_row(1)
     at_pos_ref1(2,:,n_row(2),:,2) = at_pos_ref1(2,:,n_row(2),:,2)+n_row(2)
     at_pos_ref1(3,:,:,n_row(3),2) = at_pos_ref1(3,:,:,n_row(3),2)+n_row(3)
     at_pos_ref1(:,:,:,:,2) = at_pos_ref1(:,:,:,:,2)/6.

! *** the O1 [.5 .5 0] atom
     at_pos_ref1(:,:,:,:,3) = at_pos_ref1(:,:,:,:,3)+at_pos1_file(1:3,:,:,:,j_ox+1)+at_pos1_file(1:3,:,:,:,j_ox+2)
     at_pos_ref1(:,:,:,:,3) = at_pos_ref1(:,:,:,:,3)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),1,2)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),1,3)
     at_pos_ref1(:,:,:,:,3) = at_pos_ref1(:,:,:,:,3)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),-1,4)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),-1,4)
     at_pos_ref1(:,:,:,:,3) = at_pos_ref1(:,:,:,:,3)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox+1),1,2),-1,4)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox+2),1,3),-1,4)
     at_pos_ref1(1,n_row(1),:,:,3) = at_pos_ref1(1,n_row(1),:,:,3)+2*n_row(1)
     at_pos_ref1(2,:,n_row(2),:,3) = at_pos_ref1(2,:,n_row(2),:,3)+2*n_row(2)
     at_pos_ref1(3,:,:,1,3) = at_pos_ref1(3,:,:,1,3)-4*n_row(3)
     at_pos_ref1(:,:,:,:,3) = at_pos_ref1(:,:,:,:,3)/8.

! *** the O2 [0 .5 .5] atom
     at_pos_ref1(:,:,:,:,4) = at_pos_ref1(:,:,:,:,4)+at_pos1_file(1:3,:,:,:,j_ox)+at_pos1_file(1:3,:,:,:,j_ox+2)
     at_pos_ref1(:,:,:,:,4) = at_pos_ref1(:,:,:,:,4)+cshift(at_pos1_file(1:3,:,:,:,j_ox),1,4)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),1,3)
     at_pos_ref1(:,:,:,:,4) = at_pos_ref1(:,:,:,:,4)+cshift(at_pos1_file(1:3,:,:,:,j_ox),-1,2)+cshift(at_pos1_file(1:3,:,:,:,j_ox+2),-1,2)
     at_pos_ref1(:,:,:,:,4) = at_pos_ref1(:,:,:,:,4)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox),1,4),-1,2)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox+2),1,3),-1,2)
     at_pos_ref1(3,:,:,n_row(3),4) = at_pos_ref1(3,:,:,n_row(3),4)+2*n_row(3)
     at_pos_ref1(2,:,n_row(2),:,4) = at_pos_ref1(2,:,n_row(2),:,4)+2*n_row(2)
     at_pos_ref1(1,1,:,:,4) = at_pos_ref1(1,1,:,:,4)-4*n_row(1)
     at_pos_ref1(:,:,:,:,4) = at_pos_ref1(:,:,:,:,4)/8.

! *** the O3 [.5 0 .5] atom
     at_pos_ref1(:,:,:,:,5) = at_pos_ref1(:,:,:,:,5)+at_pos1_file(1:3,:,:,:,j_ox)+at_pos1_file(1:3,:,:,:,j_ox+1)
     at_pos_ref1(:,:,:,:,5) = at_pos_ref1(:,:,:,:,5)+cshift(at_pos1_file(1:3,:,:,:,j_ox),1,4)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),1,2)
     at_pos_ref1(:,:,:,:,5) = at_pos_ref1(:,:,:,:,5)+cshift(at_pos1_file(1:3,:,:,:,j_ox),-1,3)+cshift(at_pos1_file(1:3,:,:,:,j_ox+1),-1,3)
     at_pos_ref1(:,:,:,:,5) = at_pos_ref1(:,:,:,:,5)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox),1,4),-1,3)+cshift(cshift(at_pos1_file(1:3,:,:,:,j_ox+1),1,2),-1,3)
     at_pos_ref1(1,n_row(1),:,:,5) = at_pos_ref1(1,n_row(1),:,:,5)+2*n_row(1)
     at_pos_ref1(3,:,:,n_row(3),5) = at_pos_ref1(3,:,:,n_row(3),5)+2*n_row(3)
     at_pos_ref1(2,:,1,:,5) = at_pos_ref1(2,:,1,:,5)-4*n_row(2)
     at_pos_ref1(:,:,:,:,5) = at_pos_ref1(:,:,:,:,5)/8.

! *** the A-atom from B-atoms (or vice-versa)
     at_pos_ref1(:,:,:,:,6) = at_pos_ref1(:,:,:,:,6)+at_pos1_file(1:3,:,:,:,2)+cshift(at_pos1_file(1:3,:,:,:,2),-1,2)+cshift(at_pos1_file(1:3,:,:,:,2),-1,3)+cshift(cshift(at_pos1_file(1:3,:,:,:,2),-1,2),-1,3)
     at_pos_ref1(:,:,:,:,6) = at_pos_ref1(:,:,:,:,6)+cshift(at_pos1_file(1:3,:,:,:,2),-1,4)+cshift(cshift(at_pos1_file(1:3,:,:,:,2),-1,2),-1,4)
     at_pos_ref1(:,:,:,:,6) = at_pos_ref1(:,:,:,:,6)+cshift(cshift(at_pos1_file(1:3,:,:,:,2),-1,3),-1,4)+cshift(cshift(cshift(at_pos1_file(1:3,:,:,:,2),-1,2),-1,3),-1,4)
     at_pos_ref1(1,1,:,:,6) = at_pos_ref1(1,1,:,:,6)-4*n_row(1)
     at_pos_ref1(2,:,1,:,6) = at_pos_ref1(2,:,1,:,6)-4*n_row(2)
     at_pos_ref1(3,:,:,1,6) = at_pos_ref1(3,:,:,1,6)-4*n_row(3)
     at_pos_ref1(:,:,:,:,6) = at_pos_ref1(:,:,:,:,6)/8.
   endif

! ***** calculate NN distances      
    do jat=1,n_corr
      do k=1,n_row(3)    
        do j=1,n_row(2)    
          do i=1,n_row(1)
            jj = 1 
            res2 = .0
            if(jat<=n_atom) then
              res = at_pos_file(1:3,i,j,k,jat)-at_pos_ref(1:3,i,j,k,jat)
              displ_field(:,i,j,k,jat,jt) = res
              if(norm2(res)>1.) then
                print *,'i,j,k,jat,at_pos_file(1:3,i,j,k,jat),at_pos_ref(1:3,i,j,k,jat)',i,j,k,jat,at_pos_file(1:3,i,j,k,jat),at_pos_ref(1:3,i,j,k,jat)
                read(*,*)
              endif
            else
              res = at_pos_file(1:3,i,j,k,1)-at_pos_ref(1:3,i,j,k,jat)      ! "1"   ???????
            endif

            ind_h = n_hist/2+1+nint(res/hist_step)
            if(minval(ind_h)<1.or.maxval(ind_h)>n_hist) cycle
            at_pos_hist(ind_h(1),ind_h(2),ind_h(3),jat,1,jt) = at_pos_hist(ind_h(1),ind_h(2),ind_h(3),jat,1,jt)+1
            res2 = dot_product(res,res)
            displ_norm(i,j,k,jat,jt) = sqrt(res2)   !displacement magnitude is calculated for CORES only
  
            if(jat<=n_atom) then
              if(n_dom==1) then
                jj = maxloc((abs(res)),dim=1)
              else if (n_dom==2) then
                ii = minloc((abs(res)),dim=1)
                i2 = ii+1
                if(i2.gt.3) i2 = i2-3
                i3 = ii+2
                if(i3.gt.3) i3= i3-3
                jj = mod(ii,3)*2+abs(sign(1.,res(i2))+sign(1.,res(i3)))*.5+1
              else if (n_dom==3) then
                jj = 0
                do ii=1,3
                  jj = jj+(2**(ii-1))*(sign(1.0,res(ii))+1.)*.5
                enddo 
                if(jj>=0.and.jj<=3) then
                  jj = 4-jj
                else
                  jj = jj-3
                endif
              endif
              i_dom(i,j,k,jat,jt) = jj
            endif
            
!               print *,'n_dom,res,jj,ii,i2,i3',n_dom,res,jj,ii,i2,i3
!               read(*,*)

            if(j_shell_out==1) then
              if(jat<=n_atom) then
                res = at_pos1_file(1:3,i,j,k,jat)-at_pos_ref1(1:3,i,j,k,jat)
                if(norm2(res)>1.) then
                  print *,'i,j,k,jat,at_pos1_file(1:3,i,j,k,jat),at_pos_ref1(1:3,i,j,k,jat)',i,j,k,jat,at_pos1_file(1:3,i,j,k,jat),at_pos_ref1(1:3,i,j,k,jat)
                  read(*,*)
                endif
                displ_field1(:,i,j,k,jat,jt) = res
              else
                res = at_pos1_file(1:3,i,j,k,1)-at_pos_ref1(1:3,i,j,k,jat)     ! "1"   ???????
              endif
              ind_h = n_hist/2+1+nint(res/hist_step)
              if(minval(ind_h)>=1.and.maxval(ind_h)<=n_hist) at_pos_hist(ind_h(1),ind_h(2),ind_h(3),jat,2,jt) = at_pos_hist(ind_h(1),ind_h(2),ind_h(3),jat,2,jt)+1

              if(jat<=n_atom) res = at_pos1_file(1:3,i,j,k,jat)-at_pos_file(1:3,i,j,k,jat)     !CS_displacement
              ind_h = n_hist/2+1+nint(res/hist_step)
              if(minval(ind_h)>=1.and.maxval(ind_h)<=n_hist) at_pos_hist(ind_h(1),ind_h(2),ind_h(3),jat,3,jt) = at_pos_hist(ind_h(1),ind_h(2),ind_h(3),jat,3,jt)+1
            endif                   
          enddo
        enddo   
      enddo
      if(jat>n_atom) at_pos_hist(:,:,:,jat,3,jt) = 0
    enddo 
    
! ***  calculate nominal unit cell electric polarisation for regularly occupied lattice
    if(ifile==nfile_min.and.sum(at_occup_r)/n_atom==1.) then
      at_charge = at_pos_file(4,1,1,1,:)          
      if(j_verb==1) then  
        print *,prompt, 'Atom (core) charges',at_charge, 'confirm or type new values'
        read(*,*) at_charge
      endif

      if(j_shell_out==1) then
        at_charge1 = at_pos1_file(4,1,1,1,:) 
        if(j_verb==1) then  
          print *,prompt, 'Shell charges',at_charge1, 'confirm or type new values'
          read(*,*) at_charge1
        endif
      endif 

      charge_mom_c = 0.
      charge_mom_s = 0.
      do jat=1,n_atom                       
        charge_mom_c = charge_mom_c + at_charge(jat)*at_base(jat,:)          !core   
        if(j_shell_out==1) charge_mom_s = charge_mom_s + at_charge1(jat)*at_base(jat,:)         !shell
      enddo
      charge_mom_cell = charge_mom_c+charge_mom_s 
      if(j_shell_out==0) then
        print *,space, 'Static unit cell dipole moment:',charge_mom_c     
      else
        print *,space, 'Static unit cell dipole moment: cores ',charge_mom_c
        print *,space, '                                shells',charge_mom_s
        print *,space, '                                total ',charge_mom_cell
      endif
      print *
    endif 

! ***  calculate dynamic electric polarisation of each unit cell
    do k=1,n_row(3)    
      do j=1,n_row(2)    
        do i=1,n_row(1)
          charge_mom_c = 0.
          charge_mom_s = 0.
          do jat=1,n_atom                       
              charge_mom_c = charge_mom_c + at_pos_file(4,i,j,k,jat)*(at_pos_file(1:3,i,j,k,jat)-at_pos_file(1:3,i,j,k,1))          !core   !use true local charges
              if(j_shell_out==1) charge_mom_s = charge_mom_s + at_pos1_file(4,i,j,k,jat)*(at_pos1_file(1:3,i,j,k,jat)-at_pos_file(1:3,i,j,k,1))        !shell
          enddo

          pol = charge_mom_c+charge_mom_s
          polar(:,i,j,k,jt) = pol
          pol_norm(i,j,k,jt) = norm2(pol)
! *** calculate the domain (quadrant) number for polar(ization)
          if(n_dom==1) then
            jj = maxloc((abs(pol)),dim=1)
          else if (n_dom==2) then
            ii = minloc((abs(pol)),dim=1)
            i2 = ii+1
            if(i2.gt.3) i2 = i2-3
            i3 = ii+2
            if(i3.gt.3) i3= i3-3
            jj = mod(ii,3)*2+(sign(1.,pol(i3))+1)/2
          else if (n_dom==3) then
            jj = 1
            do ii=1,3
              jj = jj+2**(ii-1)*(sign(1.0,pol(ii))+1.)/2.-4
              if(jj<1) jj = 1-jj
            enddo 
          endif
          i_dom_p(i,j,k,jt) = jj        
        enddo
      enddo   
    enddo  

    
    do i=1,3
      polar_tot(i,jt) = sum(polar(i,:,:,:,jt))
    enddo
    polar_tot(:,jt) = polar_tot(:,jt)*a_par*e_c/(nsuper*product(a_par))
    pol_norm_tot(jt) = norm2(polar_tot(:,jt))
    
    pol_norm_tot2(jt) = .0
    do k=1,n_row(3)
      do j=1,n_row(2)
        do i=1,n_row(1)
          pol_norm_tot2(jt) = pol_norm_tot2(jt)+norm2(polar(:,i,j,k,jt)*a_par)
        enddo
      enddo
    enddo
    pol_norm_tot2(jt) = pol_norm_tot2(jt)*e_c/(nsuper*product(a_par))
    
    if(j_verb==1) then
      print *,space,'Total polarisation vector [C/m^2]',jt,polar_tot(:,jt)
      print *,space,'Total polarisation norm, sum_norm [C/m^2]',jt,pol_norm_tot(jt),pol_norm_tot2(jt)
      print *
    endif
    
  enddo file_loop

  jt_max = jt   ! jt_max = 1 means a single explicit JT while t_single means absence of JT numbers
  j_ps = 0
  j_sign = 1
  

!!! *** Polarisation domains:
!!         type {1 0 0}: 1 = [-1 0 0], 2 = [0 -1 0], 3 = [0 0 -1], 4 = [1 0 0], 5 = [0 1 0], 6 = [0 0 1]
!!
!!         type {1 1 0}: 1 = [-1 -1 0],   2 = [-1 1 0],  3 = [1 -1 0],  4 = [1 1 0], 
!!                       5 = [0 -1 -1],   6 = [0 -1 1],  7 = [0 1 -1],  8 = [0 1 1], 
!!                       9 = [-1  0 -1], 10 = [1 0 -1], 11 = [-1 0 1], 12 = [1 0 1]
!!
!!         type {1 1 1}: 1 = [-1 -1 -1], 2 = [1 -1 -1], 3 = [-1 1 -1], 4 = [1 1 -1],        !binary code from left: take the -1 for 0 & add +1
!!                       5 = [-1 -1  1], 6 = [1 -1  1], 7 = [-1 1  1], 8 = [1 1  1]

!     select case(n_dom)      !n_dom becomes the number of possible domains
!       case(0)
!         n_dom = 1
!       case(1)
!         n_dom = 6
!       case(2)
!         n_dom = 12
!       case(3)
!         n_dom = 8
!     end select
!
!     allocate(dom_ind(n_dom),mask(n_dom))
!
!     select case(n_dom)      !n_dom becomes the number of possible domains
!       case(1)
!         dom_ind = ''
!       case(6)
!         dom_ind = (/'[-1 0 0]','[ 0-1 0]','[ 0 0-1]','[ 1 0 0]','[ 0 1 0]','[ 0 0 1]'/)
!       case(12)
!         dom_ind = (/'[-1-1 0]','[-1 1 0]','[ 1-1 0]','[ 1 1 0]','[ 0-1-1]','[ 0-1 1]','[ 0 1-1]','[ 0 1 1]','[-1 0-1]','[ 1 0-1]','[-1 0 1]','[ 1 0 1]'/)
!       case(8)
!         dom_ind = (/'[-1-1-1]','[ 1-1-1]','[-1 1-1]','[ 1 1-1]','[-1-1 1]','[ 1-1 1]','[-1 1 1]','[ 1 1 1]'/)
!     end select
! 
  select case(n_dom)      !n_dom becomes the number of possible domains
    case(0)
      n_dom = 1
    case(1)
      n_dom = 3
    case(2)
      n_dom = 6
    case(3)
      n_dom = 4
  end select

  allocate(dom_ind(n_dom),mask(n_dom))

  select case(n_dom)      !n_dom becomes the number of possible domains
    case(1)
      dom_ind = ''
    case(3)
      dom_ind = (/'[ 1 0 0]','[ 0 1 0]','[ 0 0 1]'/)
    case(6)
      dom_ind = (/'[ 1-1 0]','[ 1 1 0]','[ 0 1-1]','[ 0 1 1]','[-1 0 1]','[ 1 0 1]'/)
    case(4)
      dom_ind = (/'[-1-1 1]','[ 1-1 1]','[-1 1 1]','[ 1 1 1]'/)
  end select
  


! **** choose displayed variable and range
  plot_loop: do      
    scale = 1
            
!         print *,'Display: 1 displacement domains, n_dom(6[100], 12[110], 8[111]) '
!         print *,'         2 displacement magnitude, scale '
!         print *,'         3 displacement out-of-plane, scale '
!         print *,'         4 in-plane displacement vectors, scale (~50)'
!         print *,'         5 polarisation domains, n_dom(6[100], 12[110], 8[111]) '
!         print *,'         6 polarisation magnitude, scale '
!         print *,'         7 polarisation out-of-plane, scale '
!         print *,'         8 in-plane polarisation vectors, scale (~50)'
!         print *,'         0 exit (0 0)'

    print *,prompt, 'Display: 1 displacement in-plane (vector), scale (~20) '
    print *,space, '         2 displacement in direction (value), scale (~5) '
    print *,space, '         3 velocity in-plane (vector), scale (~1) '
    print *,space, '         4 velocity in direction (value), scale (~1)'
    print *,space, '         5 polarisation in-plane (vector), scale (~5)  '
    print *,space, '         6 polarisation in direction (value), scale (~1) '
    print *,space, '         7 displacement histogram, scale (= c_max) '
    print *,space, '         0 exit (0 0)'
    
!     if(j_verb==1) then
      print *
      print *,space, '         10 displacement domains, 1 '
      print *,space, '         11 atom mass, .01 '          
      print *,space, '         12 atom charge, .01 '          
      print *,space, '         13 displacement magnitude, scale (~20) '
      print *,space, '         14 displacement component, scale (~20) '
      print *,space, '         15 velocity magnitude, scale (~1) '
      print *,space, '         16 velocity component, scale (~1) '
      print *,space, '         17 polarisation magnitude, scale(~5)  '
      print *,space, '         18 polarisation component, scale (~5) '
      print *,space, '         19 polarisation domains, 1 '          
      print *,space, '         20 bond length, (~5) '          
!     endif

    fmin = 0.
    read(*,*) mode,scale
    if(mode==2.or.mode==4.or.mode==6) then
      print *,prompt, 'Direction vector components:'
      read(*,*) ed_norm
      ed_norm = ed_norm/norm2(ed_norm)
    endif                             

    if(mode==7) then
      print *,prompt, 'displacements: 1=cores, 2=shells, 3=CS_shift'
      read(*,*) j_cs
    endif                             

    if(mode==0) exit plot_loop

    j_cycle = 1
    j_auto = 0     !start by cycling layers and advancing manually to set the right scale etc.

! *** reset masks          
    mask = 1

! *** define plot geometry
    if (mode==14.or.mode==16.or.mode==18) then
      print *,prompt, 'Display plane ((0=general, 1=(hk0), 2=(hhl)), component index [1-3]'
      read(*,*) j_plane,ind_c
    else
      print *,prompt, 'Display plane ((0=general, 1=(hk0), 2=(hhl)):'
      read(*,*) j_plane
    endif

    if(j_plane==1)then
      e1 = (/1,0,0/)
      e2 = (/0,1,0/)
    elseif(j_plane==2)then
      e1 = (/1,1,0/)
      e2 = (/0,0,1/)
    else						
      do
        print *,prompt, 'Input perpendicular vectors (integer) to define the Q-plane e1(3),e2(3):'
        read(*,*) e1,e2
        if(dot_product(e1,e1).ne.0..and.dot_product(e2,e2).ne.0..and.dot_product(e1,e2)==0.) then
          exit
        else
          cycle
        endif
      enddo
    endif

    ev = (/e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)/)
    write(ev_str,'(3i2)') ev
    e1_norm = e1/sqrt(dot_product(1.*e1,1.*e1))
    e2_norm = e2/sqrt(dot_product(1.*e2,1.*e2))
    ev_norm = ev/sqrt(dot_product(1.*ev,1.*ev))
    
    ii = 0
    jj = 100
    do i=1,3
      if(abs(ev(i))==0) cycle
      if(abs(ev(i))<jj) then
        jj = ev(i)
        ii = i
      endif
    enddo
    
    e_slice = 0
    e_slice(ii) = sign(1,jj)
    n_slice = n_row(ii)           !number of horizontal slices corresponding to EV
    if(mode==7) n_slice = n_hist

!     define n_x,n_y according to geometry
    jj = maxloc(e1,1)
    n_x = n_row(jj)
    jj = maxloc(e2,1)
    n_y = n_row(jj)
    
    if(mode==7)  n_x = n_hist
    if(mode==7)  n_y = n_hist
    
!       print *,space, 'n_x,n_y',n_x,n_y
!       print *,space, 'e1,e2,ev',e1,e2,ev
!       print *,space, 'e_slice',e_slice
    
    e1p = matmul(a_cell,real(e1))			
    e2p = matmul(a_cell,real(e2))
    evp = matmul(a_cell,real(ev))
    e1p_norm = norm2(e1p)			
    e2p_norm = norm2(e2p)
    evp_norm = norm2(evp)
    ep_angle = dot_product(e1p,e2p)/(e1p_norm*e1p_norm)

    if(j_verb==1) then
      print *,space, 'e1p,e1p_norm',e1p,e1p_norm
      print *,space, 'e2p,e2p_norm',e2p,e2p_norm
      print *,space, 'evp,evp_norm',evp,evp_norm
      print *,space, 'ep_angle',ep_angle 		
    endif

    
    if(.not.allocated(displ_plot)) allocate(hist_plot(n_x,n_y,n_corr),displ_plot(2,n_x,n_y,n_corr),displ_vect(2,n_x,n_y),i_dom_out(n_x,n_y,n_atom+1),ind(3,n_x,n_y))
    
    jat = 0             !display all atoms initially
    j_atom = 0
    mode_s = .false.
    c_max_save = .0

    atom_loop: do
      if(mode==20) then
        print *,prompt, 'Atom pair numbers: (0 0 =END)'
        read(*,*) jat,j_atom
        if(jat==0.or.j_atom==0) exit atom_loop          

      elseif(mode==7) then
        jat = 0         
        j_atom = 0

      else
        print *,prompt, 'Atom number for domain & bond reference (0=NONE, -1=POLARISATION, 99=END)'
        read(*,*) j_atom
        if(j_atom==99) exit atom_loop
        if(j_atom<-1.or.j_atom>n_atom) then
          print *,space, 'wrong atom number(s), retype ...'
          cycle atom_loop
        endif
        if(j_atom==0.and.mode==20) then
          print *,space, 'Option 20: reference atom is obligatory, retype ...'
          cycle atom_loop
        endif
        if(j_atom==-1) then
          dom_name = 'POLAR'
          j_atom = n_atom+1
        elseif(j_atom==0) then
          dom_name = ''
        else
          dom_name = at_name_par(j_atom)
        endif
      endif
      
! **** collect map data, start from the beginning
      j_frame = 0
      j_auto = 0
      jt = jt0 
      j_shift = 1
    
      if(mode==7) then
        j_slice = n_hist/2+1        
        x1 = -hist_step(1)*(n_x/2)          ! plot frame
        x2 = hist_step(1)*(n_x/2)
        y1 = -hist_step(2)*(n_y/2)
        y2 = hist_step(2)*(n_y/2)
      else
        j_slice = 1
        x1 = 1.          ! plot frame
        x2 = n_x
        y1 = 1.
        y2 = n_y
      endif

      slice_loop: do    !cycle over t-snapshots or supercell layers          
        k = j_slice
        write(z_str,'(i3)') k   !for MODE=7 will be overwritten later on
        
        i_dom_out = 0  

! *** generate index table and domain reference for the actual slice
        do i=1,n_x         
          do j=1,n_y
            ind(:,i,j) = i*e1+j*e2+k*e_slice 
            if(mode/=7) then
              do ii=1,3
                ind(ii,i,j) = mod(ind(ii,i,j)-1,n_row(ii))+1
                if(ind(ii,i,j)<=0) ind(ii,i,j) = ind(ii,i,j)+n_row(ii)
              enddo
            else
              do ii=1,3
                ind(ii,i,j) = mod(ind(ii,i,j)-1,n_hist)+1
                if(ind(ii,i,j)<=0) ind(ii,i,j) = ind(ii,i,j)+n_hist
              enddo
            endif

            if(j_atom==0) i_dom_out(i,j,:) = 1
            
            if((j_atom>=1).and.(j_atom<=n_atom)) then
              do jj=1,n_corr     
!                 i_dom_out(i,j,jj) = i_dom(ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt)
                i_dom_out(i,j,jj) = i_dom(ind(1,i,j),ind(2,i,j),ind(3,i,j),j_atom,jt)
              enddo
            endif

            if(j_atom==n_atom+1) then
              do jj=1,n_corr      
                i_dom_out(i,j,jj) = i_dom_p(ind(1,i,j),ind(2,i,j),ind(3,i,j),jt)
              enddo
            endif
            
          enddo
        enddo

!     print *,prompt, 'Display: 1 displacement in-plane (vector), scale (~20) '
!     print *,space, '         2 displacement in direction (value), scale (~5) '
!     print *,space, '         3 velocity in-plane (vector), scale (~1) '
!     print *,space, '         4 velocity in direction (value), scale (~1)'
!     print *,space, '         5 polarisation in-plane (vector), scale (~5)  '
!     print *,space, '         6 polarisation in direction (value), scale (~1) '
!     print *,space, '         7 displacement histogram, scale (~1) '
!     print *,space, '         0 exit (0 0)'
!     
! ***         define displ_plot(1,n_x,n_y,jat) according to mode
       mode_pix = .false.
       mode_pts = .false.
       mode_vect = .false.
       
        select case(mode)             !scalar plots
          case(1)
            header_2 = 'displacement'
            mode_vect = .true.
            do i=1,n_x         
              do j=1,n_y
                do jj=1,n_atom
                  displ_plot(1,i,j,jj) = dot_product(displ_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt),e1_norm)
                  displ_plot(2,i,j,jj) = dot_product(displ_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt),e2_norm)
                enddo
                displ_plot(1,i,j,n_atom+1) = dot_product(polar(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt),e1_norm)
                displ_plot(2,i,j,n_atom+1) = dot_product(polar(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt),e2_norm)
              enddo
            enddo
          case(2)         
            header_2 = 'displ_along'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                do jj=1,n_atom
                  displ_plot(1,i,j,jj) = dot_product(displ_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt),ed_norm)
                enddo
                displ_plot(1,i,j,n_atom+1) = dot_product(polar(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt),ed_norm)
              enddo
            enddo
          case(3)
            header_2 = 'velocity'
            mode_vect = .true.
            do i=1,n_x         
              do j=1,n_y
                do jj=1,n_atom
                  displ_plot(1,i,j,jj) = dot_product(vel_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt),e1_norm)
                  displ_plot(2,i,j,jj) = dot_product(vel_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt),e2_norm)
                enddo
              enddo
            enddo
          case(4)         
            header_2 = 'vel_along'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                do jj=1,n_atom
                  displ_plot(1,i,j,jj) = dot_product(vel_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt),ed_norm)
                enddo
              enddo
            enddo
          case(5)         
            header_2 = 'polarisation'
            mode_vect = .true.
            jat = n_atom+1
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,jat) = dot_product(polar(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt),e1_norm)
                displ_plot(2,i,j,jat) = dot_product(polar(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt),e2_norm)
              enddo
            enddo
          case(6)         
            header_2 = 'pol_along'
            mode_pts = .true.
            jat = n_atom+1
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,jat) = dot_product(polar(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt),ed_norm)
              enddo
            enddo
          case(7) 
            header_2 = 'displ_hist'
            mode_pix = .true.
            c_max = scale
!             c_max = .0
            write(z_str,'("z = ",f7.3)') (k-(n_hist/2+1))*hist_step(3)
            do i=1,n_x         
              do j=1,n_y
!                 hist_plot(i,j,:) = at_pos_hist(ind(1,i,j),ind(2,i,j),ind(3,i,j),:,1,jt)
!                 do jj=1,n_atom
                do jj=1,n_corr
                  displ_plot(1,i,j,jj) = sum(at_pos_hist(ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,j_cs,:))
                enddo
              enddo
            enddo
!     if(j_verb==1) then
!       print *
!       print *,space, '         10 displacement domains, 1 '
!       print *,space, '         11 atom mass, .01 '          
!       print *,space, '         12 atom charge, .01 '          
!       print *,space, '         13 displacement magnitude, scale (~20) '
!       print *,space, '         14 displacement component, scale (~20) '
!       print *,space, '         15 velocity magnitude, scale (~1) '
!       print *,space, '         16 velocity component, scale (~1) '
!       print *,space, '         17 polarisation magnitude, scale(~5)  '
!       print *,space, '         18 polarisation component, scale (~5) '
!       print *,space, '         19 polarisation domains, 1 '          
!       print *,space, '         20 bond length, (~5) '          
!     endif
! 
          case(10)         
            header_2 = 'displ_domains'
            mode_pix = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = i_dom(ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom,jt)
              enddo
            enddo
          case(11)         
            header_2 = 'atom_mass'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = at_vel_file(4,ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom)   !take masses from the last input frame
              enddo
            enddo
          case(12)         
            header_2 = 'atom_charget'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = at_pos_file(4,ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom)   !take charges from the last input frame
              enddo
            enddo
          case(13)                
            header_2 = 'displ_magnitude'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = displ_norm(ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom,jt)
              enddo
            enddo
          case(14)         
            header_2 = 'displ_component'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = displ_field(ind_c,ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom,jt)
              enddo
            enddo
          case(15)                
            header_2 = 'vel_magnitude'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = vel_norm(ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom,jt)
              enddo
            enddo
          case(16)         
            header_2 = 'vel_component'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,1:n_atom) = vel_field(ind_c,ind(1,i,j),ind(2,i,j),ind(3,i,j),1:n_atom,jt)
              enddo
            enddo
          case(17)         
            header_2 = 'polar_magnitude'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,n_atom+1) = pol_norm(ind(1,i,j),ind(2,i,j),ind(3,i,j),jt)
              enddo
            enddo
          case(18)         
            header_2 = 'polar_component'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,n_atom+1) = polar(ind_c,ind(1,i,j),ind(2,i,j),ind(3,i,j),jt)
              enddo
            enddo
          case(19)         
            header_2 = 'polar_domains'
            mode_pix = .true.
            do i=1,n_x         
              do j=1,n_y
                displ_plot(1,i,j,j_atom) = i_dom_p(ind(1,i,j),ind(2,i,j),ind(3,i,j),jt)
              enddo
            enddo
          case(20)         
            header_2 = 'bond_length'
            mode_pts = .true.
            do i=1,n_x         
              do j=1,n_y
                do jj=1,n_atom
                  displ_plot(1,i,j,jj) = norm2(displ_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),jj,jt)-&
                &      displ_field(1:3,ind(1,i,j),ind(2,i,j),ind(3,i,j),j_atom,jt))
                enddo
              enddo
            enddo
        end select

      if(j_verb==1) then
        print *,space, 'Test output vs. input (last frame): atom, slice',jat,k
        do
          print *,prompt, 'Position: i,j? (0 0 exit)'
          read(*,*) i,j
          if(i==0.or.j==0) exit
          print *,space, 'Core:  at_pos_file ', at_pos_file(:,i,j,k,jat)
          print *,space, 'Core:  displ_field ', displ_field(:,i,j,k,jat,jt),i_dom(i,j,k,jat,jt)
          if(j_shell_out==1) then
            print *,space, 'Shell: at_pos1_file', at_pos1_file(:,i,j,k,jat)
            print *,space, 'Shell: displ_field1', displ_field1(:,i,j,k,jat,jt)
          endif
          print *,space, 'Plot:  displ_plot  ', displ_plot(1,i,j,jat)
        enddo
      endif

!! *** PGPLOT: set the coordinate transformation matrix: world coordinate = pixel number.
!!
!!    The transformation matrix TR is used to calculate the world
!!    coordinates of the center of the "cell" that represents each
!!    array element. The world coordinates of the center of the cell
!!    corresponding to array element A(I,J) are given by:
!!             X = TR(1) + TR(2)*I + TR(3)*J
!!             Y = TR(4) + TR(5)*I + TR(6)*J
!!    Usually TR(3) and TR(5) are zero -- unless the coordinate
!!    transformation involves a rotation or shear.  The corners of the
!!    quadrilateral region that is shaded by PGIMAG are given by
!!    applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).
!!

      if(mode==7) then
          TR(1) = -hist_step(1)*(n_hist/2+1)
          TR(2) = hist_step(1)
          TR(3) = 0.0
          TR(4) = -hist_step(2)*(n_hist/2+1)
          TR(5) = 0.0
          TR(6) = hist_step(2)
      else
          TR(1) = 0.0
          TR(2) = 1.0
          TR(3) = 0.0
          TR(4) = 0.0
          TR(5) = 0.0
          TR(6) = 1.0
      endif

      do i =1,3
        write(c_e1(i),*) e1(i)
        write(c_e2(i),*) e2(i)
      enddo

      if(mode==7) then
        c_x = 'X ['//trim(adjustl(c_e1(1)))//' '//trim(adjustl(c_e1(2)))//' '//trim(adjustl(c_e1(3)))//']'
        c_y = 'Y ['//trim(adjustl(c_e2(1)))//' '//trim(adjustl(c_e2(2)))//' '//trim(adjustl(c_e2(3)))//']'
      else
        c_x = 'N_x ['//trim(adjustl(c_e1(1)))//' '//trim(adjustl(c_e1(2)))//' '//trim(adjustl(c_e1(3)))//']'
        c_y = 'N_y ['//trim(adjustl(c_e2(1)))//' '//trim(adjustl(c_e2(2)))//' '//trim(adjustl(c_e2(3)))//']'
      endif

      if(nfile_min<=9999) then
        write(file_dat,'(a,"_n",i4.4,".dat")') trim(file_master),jt
      elseif(nfile_min>=10000) then
        write(string,'(i8)') jt
        file_dat = trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
      endif

      if(jt<=9999) then
        write(line,'(" frame = ",i4.4," layer = ",i2," mode = ",i2)') jt,j_slice,mode
      else
        write(line,'(" frame = *",i4.4," layer = ",i2," mode = ",i2)') mod(jt,10000),j_slice,mode
      endif     
      plot_title = trim(file_dat)//'  '//trim(at_name)//'/'//trim(dom_name)//trim(line)

      header_2 = trim(line)//' '//trim(header_2)

      map_unit(1) = PGOPEN('/xserv')   !open the PGPLOT window, 1,1 no of panes
      CALL PGASK(.FALSE.)     ! would not ask for <RET>


!       CALL PGSHLS (21,.0,.4,.7)     !my blue
!       CALL PGSHLS (22,120.,.5,1.)   !my red
!       CALL PGSHLS (23,240.,.35,.8)  !my green
!       CALL PGSHLS (24,60.,.4,.9)    !my violet
!       CALL PGSHLS (25,170.,.5,.9)   !my yellow
!       CALL PGSHLS (26,320.,.4,.9)   !my turquoise
!       CALL PGSHLS (27,.0,.7,.0)     !light grey
!       CALL PGSHLS (28,30.,.5,1.)    !my other
!       CALL PGSHLS (29,150.,.5,1.)   !my other
!       CALL PGSHLS (30,270.,.5,1.)   !my other
!       CALL PGSHLS (31,300.,.5,1.0)  !my other
!       CALL PGSHLS (32,.0,.3,.0)     !dark grey  !set my colors


      if(mode_s) then
        CALL PGPAP(p_size,1.0)     ! define the plot area as a 1x1 square
        CALL PGSUBP (1,1)
      else
        CALL PGPAP(p_size,1.5)     ! define the plot area as a 3x2 rectangle
        CALL PGSUBP (2,3)
      endif

      CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
      CALL PGSCRN(1, 'black', IER)
        
! *** plot the histograms
!       if(mode==7) then
      if(mode_pix) then       !here come pixel maps
        c_min = .0
        if(c_max_save==.0) then
!           c_max = maxval(displ_plot(1,1:n_x,1:n_y,:))
          c_max = scale
        else
          c_max = c_max_save
        endif
 !        write(*,*) 'c_min,c_max',c_min,c_max
        print *,space,'Plot sum, j_slice, z',j_slice,'  ',trim(z_str),(sum(displ_plot(1,:,:,jj)),jj=1,n_corr)

        CALL PGQCIR(C1, C2)
        NC = MAX(0, C2-C1+1)
        BRIGHT = 0.5
        CONTRA  = 1.0
        CALL PALETT(j_pgc, CONTRA, BRIGHT)    !default j_pgc=6 - JK rainbow
        CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
        CALL PGSCRN(1, 'black', IER)
        CALL PGSCH(1.)					!set character height					
        CALL PGSLS (1)  !FULL
        CALL PGSLW(2)
        CALL PGSCI (1)  !black(0 = white)
  
        if(mode_s) then
          plot_title = trim(at_name_par(jat))//' '//trim(cs_string(j_cs))//' ('//trim(ev_str)//')  '//trim(z_str)//' r.l.u.'		
          CALL PGENV(x1-.5*hist_step(1),x2+.5*hist_step(1),y1-.5*hist_step(2),y2+.5*hist_step(2),0,1) !draw the axes
          CALL PGSCH(1.2)					!set character height
          CALL PGLAB('X','Y',plot_title)  					!put the axis labels
          CALL PGIMAG(displ_plot(1,1:n_x,1:n_y,jat),n_x,n_y,1,n_x,1,n_y,.0,c_max,TR)    
          CALL PGSCH(1.)					!set character height					
          CALL PGWEDG('RI', 1., 3.,.0, c_max,'')           ! R is right (else L,T,B), I is PGIMAG (G is PGGRAY)
         
          call date_and_time(c_date,c_time,c_zone,i_time)
          write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
          x_plot = x1+.75*(x2-x1)
          y_plot = y2-.10*(y2-y1)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSCH(.7)
          CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'  '//trim(file_master)//' '//trim(header_1)//'   '//trim(header_2))
          CALL PGSCH(1.)
        
          if(j_grid==1) then
            XOPT = 'G'				!draw grid lines - after the color map, otherwise obscured
            YOPT = 'G'
            XTICK = 1.0
            YTICK = 1.0
            NXSUB = 5
            NYSUB = 5
            CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
          endif  
        else
          do jj=1,min(n_corr,6)
            CALL PGENV(x1-.5*hist_step(1),x2+.5*hist_step(1),y1-.5*hist_step(2),y2+.5*hist_step(2),0,1) !draw the axes
            CALL PGLAB('X','Y','')  					!put the axis labels
            CALL PGIMAG(displ_plot(1,1:n_x,1:n_y,jj),n_x,n_y,1,n_x,1,n_y,.0,c_max,TR)    
            CALL PGSCH(1.2)					!set character height
            if(jj<=n_atom) then
              plot_title = trim(at_name_par(jj))//' '//trim(cs_string(j_cs))//' ('//trim(ev_str)//')  '//trim(z_str)//' Å'		
            else
              plot_title = trim(corr_name(jj-n_atom))//' '//trim(cs_string(j_cs))//' ('//trim(ev_str)//')  '//trim(z_str)//' Å'		
            endif
            CALL PGTEXT(x1+.1*(x2-x1),y2-.1*(y2-y1),plot_title)
            CALL PGSCH(1.)					!set character height					
            CALL PGWEDG('RI', 1., 3.,.0, c_max,'')           ! R is right (else L,T,B), I is PGIMAG (G is PGGRAY)
           
            if(jj==1) then          ! *** print header with program version & date_and_time with the 1st pane
              call date_and_time(c_date,c_time,c_zone,i_time)
              write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
              x_plot = x1         !+.75*(x2-x1)
              y_plot = y2+.05*(y2-y1)
              CALL PGSCI (1)  !white needs to be reset after PGLAB
              CALL PGSCH(1.4)
              CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'  '//trim(file_master)//' '//trim(header_1)//'   '//trim(header_2))
              CALL PGSCH(1.)
            endif
          
            if(j_grid==1) then
              XOPT = 'G'				!draw grid lines - after the color map, otherwise obscured
              YOPT = 'G'
              XTICK = 1.0
              YTICK = 1.0
              NXSUB = 5
              NYSUB = 5
              CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
            endif
          enddo
        endif

! *******************************************

! *** plot the scalar fields
      elseif(mode_pts) then             !here come points
!         print *,'mode_pts'
        displ_plot = displ_plot*scale   ! for vectors the SCALE is embedded in PGVECT
        displ_plot(:,:,:,n_atom+1) = displ_plot(:,:,:,n_atom+1)*.2  !reduce the scale for POLAR

        if(mode_s) then
          if(jat<=n_atom) then
            plot_title = trim(file_dat)//'  '//trim(at_name_par(jat))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          else
             plot_title = trim(file_dat)//'  '//trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          endif
          CALL PGSLS (1)  !FULL
          CALL PGSLW(2)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSCH(.7)
          CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
          CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
          CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
          CALL PGQCIR(C1, C2)
          NC = MAX(0, C2-C1+1)
          BRIGHT = 0.5
          CONTRA  = 1.0
          CALL PALETT(6, CONTRA, BRIGHT)
          CALL PGSLW (2)          ! *** write the domain legend
          CALL PGSCI(1)
          CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
          CALL PGSCH(.8)
          do k=1,n_dom
            if(mask(k)==0) cycle
!             CALL PGSCI(20+k)
            CALL PGSCI(k)
            CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
          enddo
          do i=1,n_x         
            do j=1,n_y
              if(mask(i_dom_out(i,j,jat))==0) cycle
!               CALL PGSCI(20+i_dom_out(i,j,jat))
              CALL PGSCI(i_dom_out(i,j,jat))
              if(n_dom==1) CALL PGSCI(24)
              if(displ_plot(1,i,j,jat)>=.0)then
                CALL PGSFS (1)            !full circles
              else
                CALL PGSFS (2)            !hollow circles
              endif              
              CALL PGCIRC (real(i),real(j),sqrt(abs(displ_plot(1,i,j,jat))))
            enddo
          enddo

          call date_and_time(c_date,c_time,c_zone,i_time)
          write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
          x_plot = x1+.15*(x2-x1)
          y_plot = y1-.1*(y2-y1)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW (1)          ! *** write the domain legend
          CALL PGSCH(.5)
          CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//' '//trim(header_1)//'   '//trim(header_2))
          CALL PGSCH(1.)

        else
          do jj=1,min(n_corr,6)
            if(jj<=n_atom) then
              plot_title = trim(at_name_par(jj))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            else
              plot_title = trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            endif
            CALL PGSLS (1)  !FULL
            CALL PGSLW(2)
            CALL PGSCI (1)  !white needs to be reset after PGLAB
            CALL PGSCH(1.)
            CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
            CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
            CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
            CALL PGQCIR(C1, C2)
            NC = MAX(0, C2-C1+1)
            BRIGHT = 0.5
            CONTRA  = 1.0
            CALL PALETT(6, CONTRA, BRIGHT)
            CALL PGSLW (2)          ! *** write the domain legend
            CALL PGSCI(1)
            CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
            CALL PGSCH(.8)
            do k=1,n_dom
              if(mask(k)==0) cycle
  !             CALL PGSCI(20+k)
              CALL PGSCI(k)
              CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
            enddo
            do i=1,n_x         
              do j=1,n_y
                if(mask(i_dom_out(i,j,jj))==0) cycle
 !                CALL PGSCI(20+i_dom_out(i,j,jj))
                CALL PGSCI(i_dom_out(i,j,jj))
                if(n_dom==1) CALL PGSCI(8)
                if(displ_plot(1,i,j,jj)>=.0)then
                  CALL PGSFS (1)            !full circles
                else
                  CALL PGSFS (2)            !hollow circles
                endif              
                CALL PGCIRC (real(i),real(j),sqrt(abs(displ_plot(1,i,j,jj))))
              enddo
            enddo
            if(jj==1) then          ! *** print header with program version & date_and_time with the 1st pane
              call date_and_time(c_date,c_time,c_zone,i_time)
              write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
              x_plot = x1 +.15*(x2-x1)
              y_plot = y2+.12*(y2-y1)
              CALL PGSCI (1)  !white needs to be reset after PGLAB
              CALL PGSLW (1)          ! *** write the domain legend
              CALL PGSCH(.9)
              CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'  '//trim(file_master)//' '//trim(header_1)//'   '//trim(header_2))
              CALL PGSCH(1.)
            endif
          enddo
        endif

! *** plot the vector fields
      elseif(mode_vect) then             !here go vector fields
        BLANK = .0    ! blanking
        NC = 0        !arrows centred
        CALL PGSCH(.4)

        if(mode_s) then
          if(jat<=n_atom) then
            plot_title = trim(file_dat)//'  '//trim(at_name_par(jat))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          else
             plot_title = trim(file_dat)//'  '//trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          endif
          CALL PGSLS (1)  !FULL
          CALL PGSLW(2)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSCH(.7)
          CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
          CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
          CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
          CALL PGQCIR(C1, C2)
          NC = MAX(0, C2-C1+1)
          BRIGHT = 0.5
          CONTRA  = 1.0
          CALL PALETT(6, CONTRA, BRIGHT)
          CALL PGSLW (2)          ! *** write the domain legend
          CALL PGSCI(1)
          CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
          CALL PGSCH(.8)
          do k=1,n_dom
            if(mask(k)==0) cycle
!             CALL PGSCI(20+k)
            CALL PGSCI(k)
            CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
          enddo
          do j=1,n_dom
            if(mask(j)==0)cycle
            displ_vect = .0
            where(i_dom_out(:,:,jat)==j)
              displ_vect(1,:,:) = displ_plot(1,:,:,jat)
              displ_vect(2,:,:) = displ_plot(2,:,:,jat)
            end where
            CALL PGSCI(j)
!            CALL PGSCI(20+j)
            if(jat<=n_atom) then
              CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,scale,nc,TR,blank)   
            else
              CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,.2*scale,nc,TR,blank)   
            endif           
          enddo

          call date_and_time(c_date,c_time,c_zone,i_time)
          write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
          x_plot = x1+.75*(x2-x1)
          y_plot = y1-.1*(y2-y1)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW (1)          ! *** write the domain legend
          CALL PGSCH(.5)
          CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp))
          CALL PGSCH(1.)

        else
          do jj=1,min(n_corr,6)
            if(jj<=n_atom) then
              plot_title = trim(at_name_par(jj))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            else
              plot_title = trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            endif
            CALL PGSLS (1)  !FULL
            CALL PGSLW(2)
            CALL PGSCI (1)  !white needs to be reset after PGLAB
            CALL PGSCH(1.)
            CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
            CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
            CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
            CALL PGQCIR(C1, C2)
            NC = MAX(0, C2-C1+1)
            BRIGHT = 0.5
            CONTRA  = 1.0
            CALL PALETT(6, CONTRA, BRIGHT)
            CALL PGSLW (2)          ! *** write the domain legend
            CALL PGSCI(1)
            CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
            CALL PGSCH(.8)
            do k=1,n_dom
              if(mask(k)==0) cycle
  !             CALL PGSCI(20+k)
              CALL PGSCI(k)
              CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
            enddo
            do j=1,n_dom
              if(mask(j)==0)cycle
              displ_vect = .0
              where(i_dom_out(:,:,jj)==j)
                displ_vect(1,:,:) = displ_plot(1,:,:,jj)
                displ_vect(2,:,:) = displ_plot(2,:,:,jj)
              end where
!               CALL PGSCI(20+j)
              CALL PGSCI(j)
              if(jj<=n_atom) then
                CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,scale,nc,TR,blank)
              else
                CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,.2*scale,nc,TR,blank)
             endif              
            enddo
            if(jj==1) then          ! *** print header with program version & date_and_time with the 1st pane
              call date_and_time(c_date,c_time,c_zone,i_time)
              write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
              x_plot = x1 +.7*(x2-x1)
              y_plot = y2+.12*(y2-y1)
              CALL PGSCI (1)  !white needs to be reset after PGLAB
              CALL PGSLW (1)          ! *** write the domain legend
              CALL PGSCH(.9)
              CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'       '//trim(file_master))
              CALL PGSCH(1.)
            endif
          enddo
        endif
        
      else
        print *,'ERR: no plot option selected'
      endif
    call PGCLOS

    if(j_auto.eq.0) then
      print *,prompt, 'Increment (OPTIONS = 0, END = 99): [1]'
      read(*,*) j_shift
      if(j_shift==99) exit atom_loop

      if(j_shift.eq.0) then
        options_loop: do
          print *,prompt, 'Options (type input after option no.):' 
          if(j_auto==0) print *,space, '   1   toggle MAN to AUTO advance'          !j_auto
          if(j_auto==1) print *,space, '   1   toggle AUTO to MAN advance'          !j_auto
          if(jt_max>1.and.j_frame==0)print *,space, '   2   toggle SLICE to FRAME advance' !j_frame
          if(jt_max>1.and.j_frame==1)print *,space, '   2   toggle FRAME to SLICE advance' !j_frame
          write(*,"(10x,'    3   ZOOM (actual frame:',4f6.1,')')") x1,y1,x2,y2
          if(mode==7) then
            print *,space, '   4   modify C_MAX (actual',c_max,')'
          else
            print *,space, '   4   modify scale (actual',scale,', -1 invert sign)'
          endif
          write(*,"(10x,'    5   mask domains (actual :',12i2,')')",advance='no') mask
                                                           write(*,"(')')")
          print *,space, '   6   record graphics output'
          print *,space, '   7   display atom no. (0=ALL, 6=POLAR) '
          print *,space, '   0   EXIT'
!           read(*,*) j_op

          read(*,'(a)') string
          read(string,*) j_op

          select case (j_op)
            case(1)
              j_auto = j_auto+1
              j_auto = mod(j_auto,2)
              if(j_auto==1) then
                print *,prompt, '       increment & time delay (≥1 sec): '
                read(*,*) j_shift,j_adv
                i_shift = 0
              endif

            case(2)
              j_frame = j_frame+1
              j_frame = mod(j_frame,2)
            case(3)
              read(string,*,iostat=ios) j_op,x1,y1,x2,y2
              if(ios/=0) then     
                print *,prompt, '       corner indices (bottom left & top right, confirm/re-type):'
                read(*,*) x1,y1,x2,y2
              endif

            case(4)
              read(string,*,iostat=ios) j_op,res2
              if(ios/=0) then     
                print *,prompt, '         new scale factor (-1 invert sign)'
                read(*,*) res2
              endif
              if(res2==-1) then
                scale = -scale
              else
                scale = res2
              endif
              if (mode==7) c_max_save = abs(res2)
              exit options_loop
            case(5)
              read(string,*,iostat=ios) j_op,mask
              if(ios/=0) then              
                write(*,"('        confirm or retype masks: ')",advance='no')
                read(*,*) mask
              endif
            case(6)
              j_ps = 1
              exit options_loop
            case(7)
              read(string,*,iostat=ios) j_op,jat
              if(ios/=0) then     
                print *,space,'Atom number (0=ALL): '
                read(*,*) jat
              endif
              mode_s = .true.
              if(jat==0) mode_s = .false.
              cycle slice_loop
            case(0)
              j_shift = 1
              exit options_loop
          end select
        enddo options_loop
      endif
    else
      call sleep(j_adv)
    endif
          
!! ***  save GRAPHICS 

    if(j_ps.eq.1) then  
      if(jt<=9999) then
        write(c_jt,'(i4.4)') jt
      elseif(jt>=10000) then
        write(c_jt,'(i8)') jt
      endif
      c_jt = '_'//adjustl(c_jt)      

      write(c_mode,'("_",i2.2)') mode
      if(j_slice<=99) then
        write(c_slice,'(i2.2)') j_slice
      elseif(j_slice>=100) then
        write(c_slice,'(i8)') j_slice
      endif
      c_slice = '_'//adjustl(c_slice)      

      file_ext = 'png'
      pg_out = '/PNG'

      jfile = 0
      do						!look for existing .png files to continue numbering
        write(c_jfile,'("_",i2.2)') jfile
        if(jfile==0) c_jfile=''             !don't write file count, if only a single copy


        if(jt==0)then
          file_ps  = trim(file_master)//trim(c_slice)//trim(c_mode)//trim(c_jfile)//'_dp.png'
        else			      
          file_ps  = trim(file_master)//trim(c_jt)//trim(c_slice)//trim(c_mode)//trim(c_jfile)//'_dp.png'
        endif
        
        inquire(file=file_ps,exist=found_ps)
        if(.not.found_ps) exit				

        jfile = jfile+1
        if(jfile==100) then
          print *,prompt,'Tidy up .txt/.ps/.png files to restart count from 01 and type [RET]'
          read(*,*)
          jfile = 1
        endif							
      enddo

      j_ps = PGOPEN(trim(file_ps)//trim(pg_out))
      if(j_ps.LE.0) then
        print *,space, 'The .PS file could not be opened: check your PGPLOT installation (drivers)!'
        cycle slice_loop
      endif
      print *,space, 'Saving PS file: ',file_ps

!!      CALL PGPAP(8.0,1.0)     ! define the plot area
!!      CALL PGSLW (2)
!!      CALL PGSCH(0.8)
!!      CALL PGENV(x1-.5,x2-.5,y1+.5,y2+.5,0,1) !draw the axes
!!      CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
!!      CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels
!!
!!
!!      CALL PGSHLS (20,.0,.3,.0)     !dark grey  !set my colors
!!      CALL PGSHLS (21,.0,.4,.7)     !my blue
!!      CALL PGSHLS (22,120.,.5,1.)   !my red
!!      CALL PGSHLS (23,240.,.35,.8)  !my green
!!      CALL PGSHLS (24,60.,.4,.9)    !my violet
!!      CALL PGSHLS (25,170.,.5,.9)   !my yellow
!!      CALL PGSHLS (26,320.,.4,.9)   !my turquoise
!!      CALL PGSHLS (27,.0,.7,.0)     !light grey
!!      CALL PGSHLS (28,30.,.5,1.)    !my other
!!      CALL PGSHLS (29,150.,.5,1.)   !my other
!!      CALL PGSHLS (30,270.,.5,1.)   !my other
!!      CALL PGSHLS (31,300.,.5,1.0)  !my other
!!
!!
!!! *** write the domain legend
!!      CALL PGSLW (2)
!!      CALL PGSCH(.5)
!!      do k=1,n_dom
!!        if(mask(k)==0) cycle
!!        CALL PGSCI(19+k)
!!        CALL PGMTXT ( 'RV',.5,1.-.03*k, .0,trim(dom_ind(k)))
!!!                  PGMTXT (SIDE, DISP, COORD, FJUST, TEXT)
!!      enddo
!!        
!!! *** plot the scalar fields
!!        if(mode/=1.and.mode/=3.and.mode/=5) then
!!          displ_plot = displ_plot*scale   ! for vectors the SCALE is embedded in PGVECT
!!
!!          do i=1,n_x         
!!            do j=1,n_y
!!              if(mask(i_dom_out(i,j,jat))==0) cycle
!!              CALL PGSCI(19+i_dom_out(i,j,jat))
!!              if(displ_plot(1,i,j,jat)>=.0)then
!!                CALL PGSFS (1)            !full circles
!!              else
!!                CALL PGSFS (2)            !hollow circles
!!              endif              
!!              CALL PGCIRC (real(i),real(j),sqrt(abs(displ_plot(1,i,j,jat))))
!!            enddo
!!          enddo
!!
!!! *** plot the vector fields
!!        else
!!          BLANK = .0    ! blanking
!!          NC = 0        !arrows centred
!!          CALL PGSCH(.4)
!!
!!          do j=1,n_dom
!!            if(mask(j)==0)cycle
!!            displ_vect = .0
!!            where(i_dom_out(:,:,jat)==j)
!!              displ_vect(1,:,:) = displ_plot(1,:,:,jat)
!!              displ_vect(2,:,:) = displ_plot(2,:,:,jat)
!!            end where
!!            CALL PGSCI(19+j)
!!            CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,scale,nc,TR,blank)              
!!          enddo
!!        endif

      if(mode_s) then
        CALL PGPAP(p_size,1.0)     ! define the plot area as a 1x1 square
        CALL PGSUBP (1,1)
      else
        CALL PGPAP(p_size,1.5)     ! define the plot area as a 3x2 rectangle
        CALL PGSUBP (2,3)
      endif

      CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
      CALL PGSCRN(1, 'black', IER)
        
! *** plot the histograms
!       if(mode==7) then
      if(mode_pix) then       !here come pixel maps
        c_min = .0
        if(c_max_save==.0) then
          c_max = maxval(displ_plot(1,1:n_x,1:n_y,:))
        else
          c_max = c_max_save
        endif
        write(*,*) 'c_min,c_max,c_max_save',c_min,c_max,c_max_save
        print *,space,'Plot sum, j_slice, z',j_slice,'  ',trim(z_str),(sum(displ_plot(1,:,:,jj)),jj=1,n_corr)

        CALL PGQCIR(C1, C2)
        NC = MAX(0, C2-C1+1)
        BRIGHT = 0.5
        CONTRA  = 1.0
        CALL PALETT(j_pgc, CONTRA, BRIGHT)    !default j_pgc=6 - JK rainbow
        CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
        CALL PGSCRN(1, 'black', IER)
        CALL PGSCH(1.)					!set character height					
        CALL PGSLS (1)  !FULL
        CALL PGSLW(2)
        CALL PGSCI (1)  !black(0 = white)
  
        if(mode_s) then
          plot_title = trim(at_name_par(jat))//' '//trim(cs_string(j_cs))//' ('//trim(ev_str)//')  '//trim(z_str)//' r.l.u.'		
          CALL PGENV(x1-.5*hist_step(1),x2+.5*hist_step(1),y1-.5*hist_step(2),y2+.5*hist_step(2),0,1) !draw the axes
          CALL PGSCH(1.2)					!set character height
          CALL PGLAB('X','Y',plot_title)  					!put the axis labels
          CALL PGIMAG(displ_plot(1,1:n_x,1:n_y,jat),n_x,n_y,1,n_x,1,n_y,.0,c_max,TR)    
          CALL PGSCH(1.)					!set character height					
          CALL PGWEDG('RI', 1., 3.,.0, c_max,'')           ! R is right (else L,T,B), I is PGIMAG (G is PGGRAY)
         
          call date_and_time(c_date,c_time,c_zone,i_time)
          write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
          x_plot = x1+.75*(x2-x1)
          y_plot = y2-.10*(y2-y1)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSCH(.7)
          CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'       '//trim(file_master))
          CALL PGSCH(1.)
        
          if(j_grid==1) then
            XOPT = 'G'				!draw grid lines - after the color map, otherwise obscured
            YOPT = 'G'
            XTICK = 1.0
            YTICK = 1.0
            NXSUB = 5
            NYSUB = 5
            CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
          endif  
        else
          do jj=1,min(n_corr,6)
            CALL PGENV(x1-.5*hist_step(1),x2+.5*hist_step(1),y1-.5*hist_step(2),y2+.5*hist_step(2),0,1) !draw the axes
            CALL PGLAB('X','Y','')  					!put the axis labels
            CALL PGIMAG(displ_plot(1,1:n_x,1:n_y,jj),n_x,n_y,1,n_x,1,n_y,.0,c_max,TR)    
            CALL PGSCH(1.2)					!set character height
            if(jj<=n_atom) then
              plot_title = trim(at_name_par(jj))//' '//trim(cs_string(j_cs))//' ('//trim(ev_str)//')  '//trim(z_str)//' Å'		
            else
              plot_title = trim(corr_name(jj-n_atom))//' '//trim(cs_string(j_cs))//' ('//trim(ev_str)//')  '//trim(z_str)//' Å'		
            endif
            CALL PGTEXT(x1+.1*(x2-x1),y2-.1*(y2-y1),plot_title)
            CALL PGSCH(1.)					!set character height					
            CALL PGWEDG('RI', 1., 3.,.0, c_max,'')           ! R is right (else L,T,B), I is PGIMAG (G is PGGRAY)
           
            if(jj==1) then          ! *** print header with program version & date_and_time with the 1st pane
              call date_and_time(c_date,c_time,c_zone,i_time)
              write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
              x_plot = x1         !+.75*(x2-x1)
              y_plot = y2+.05*(y2-y1)
              CALL PGSCI (1)  !white needs to be reset after PGLAB
              CALL PGSCH(1.4)
              CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'       '//trim(file_master))
              CALL PGSCH(1.)
            endif
          
            if(j_grid==1) then
              XOPT = 'G'				!draw grid lines - after the color map, otherwise obscured
              YOPT = 'G'
              XTICK = 1.0
              YTICK = 1.0
              NXSUB = 5
              NYSUB = 5
              CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
            endif
          enddo
        endif

! *******************************************

! *** plot the scalar fields
      elseif(mode_pts) then             !here come points
!         print *,'mode_pts'
        displ_plot = displ_plot*scale   ! for vectors the SCALE is embedded in PGVECT
        displ_plot(:,:,:,n_atom+1) = displ_plot(:,:,:,n_atom+1)*.2  !reduce the scale for POLAR

        if(mode_s) then
          if(jat<=n_atom) then
            plot_title = trim(file_dat)//'  '//trim(at_name_par(jat))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          else
             plot_title = trim(file_dat)//'  '//trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          endif
          CALL PGSLS (1)  !FULL
          CALL PGSLW(2)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSCH(.7)
          CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
          CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
          CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
          CALL PGQCIR(C1, C2)
          NC = MAX(0, C2-C1+1)
          BRIGHT = 0.5
          CONTRA  = 1.0
          CALL PALETT(6, CONTRA, BRIGHT)
          CALL PGSLW (2)          ! *** write the domain legend
          CALL PGSCI(1)
          CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
          CALL PGSCH(.8)
          do k=1,n_dom
            if(mask(k)==0) cycle
!             CALL PGSCI(20+k)
            CALL PGSCI(k)
            CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
          enddo
          do i=1,n_x         
            do j=1,n_y
              if(mask(i_dom_out(i,j,jat))==0) cycle
!               CALL PGSCI(20+i_dom_out(i,j,jat))
              CALL PGSCI(i_dom_out(i,j,jat))
              if(n_dom==1) CALL PGSCI(4)
              if(displ_plot(1,i,j,jat)>=.0)then
                CALL PGSFS (1)            !full circles
              else
                CALL PGSFS (2)            !hollow circles
              endif              
              CALL PGCIRC (real(i),real(j),sqrt(abs(displ_plot(1,i,j,jat))))
            enddo
          enddo

          call date_and_time(c_date,c_time,c_zone,i_time)
          write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
          x_plot = x1+.75*(x2-x1)
          y_plot = y1-.1*(y2-y1)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW (1)          ! *** write the domain legend
          CALL PGSCH(.5)
          CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp))
          CALL PGSCH(1.)

        else
          do jj=1,min(n_corr,6)
            if(jj<=n_atom) then
              plot_title = trim(at_name_par(jj))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            else
              plot_title = trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            endif
            CALL PGSLS (1)  !FULL
            CALL PGSLW(2)
            CALL PGSCI (1)  !white needs to be reset after PGLAB
            CALL PGSCH(1.)
            CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
            CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
            CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
            CALL PGQCIR(C1, C2)
            NC = MAX(0, C2-C1+1)
            BRIGHT = 0.5
            CONTRA  = 1.0
            CALL PALETT(6, CONTRA, BRIGHT)
            CALL PGSLW (2)          ! *** write the domain legend
            CALL PGSCI(1)
            CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
            CALL PGSCH(.8)
            do k=1,n_dom
              if(mask(k)==0) cycle
  !             CALL PGSCI(20+k)
              CALL PGSCI(k)
              CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
            enddo
            do i=1,n_x         
              do j=1,n_y
                if(mask(i_dom_out(i,j,jj))==0) cycle
 !                CALL PGSCI(20+i_dom_out(i,j,jj))
                CALL PGSCI(i_dom_out(i,j,jj))
                if(n_dom==1) CALL PGSCI(4)
                if(displ_plot(1,i,j,jj)>=.0)then
                  CALL PGSFS (1)            !full circles
                else
                  CALL PGSFS (2)            !hollow circles
                endif              
                CALL PGCIRC (real(i),real(j),sqrt(abs(displ_plot(1,i,j,jj))))
              enddo
            enddo
            if(jj==1) then          ! *** print header with program version & date_and_time with the 1st pane
              call date_and_time(c_date,c_time,c_zone,i_time)
              write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
              x_plot = x1 +.7*(x2-x1)
              y_plot = y2+.12*(y2-y1)
              CALL PGSCI (1)  !white needs to be reset after PGLAB
              CALL PGSLW (1)          ! *** write the domain legend
              CALL PGSCH(.9)
              CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'       '//trim(file_master))
              CALL PGSCH(1.)
            endif
          enddo
        endif

! *** plot the vector fields
      elseif(mode_vect) then             !here go vector fields
        BLANK = .0    ! blanking
        NC = 0        !arrows centred
        CALL PGSCH(.4)

        if(mode_s) then
          if(jat<=n_atom) then
            plot_title = trim(file_dat)//'  '//trim(at_name_par(jat))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          else
             plot_title = trim(file_dat)//'  '//trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
          endif
          CALL PGSLS (1)  !FULL
          CALL PGSLW(2)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSCH(.7)
          CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
          CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
          CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
          CALL PGQCIR(C1, C2)
          NC = MAX(0, C2-C1+1)
          BRIGHT = 0.5
          CONTRA  = 1.0
          CALL PALETT(6, CONTRA, BRIGHT)
          CALL PGSLW (2)          ! *** write the domain legend
          CALL PGSCI(1)
          CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
          CALL PGSCH(.8)
          do k=1,n_dom
            if(mask(k)==0) cycle
!             CALL PGSCI(20+k)
            CALL PGSCI(k)
            CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
          enddo
          do j=1,n_dom
            if(mask(j)==0)cycle
            displ_vect = .0
            where(i_dom_out(:,:,jat)==j)
              displ_vect(1,:,:) = displ_plot(1,:,:,jat)
              displ_vect(2,:,:) = displ_plot(2,:,:,jat)
            end where
            CALL PGSCI(j)
!            CALL PGSCI(20+j)
            if(jat<=n_atom) then
              CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,scale,nc,TR,blank)   
            else
              CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,.2*scale,nc,TR,blank)   
            endif           
          enddo

          call date_and_time(c_date,c_time,c_zone,i_time)
          write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
          x_plot = x1+.75*(x2-x1)
          y_plot = y1-.1*(y2-y1)
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW (1)          ! *** write the domain legend
          CALL PGSCH(.5)
          CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp))
          CALL PGSCH(1.)

        else
          do jj=1,min(n_corr,6)
            if(jj<=n_atom) then
              plot_title = trim(at_name_par(jj))//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            else
              plot_title = trim(polar_name)//' ('//trim(ev_str)//') N_z ='//trim(z_str)	
            endif
            CALL PGSLS (1)  !FULL
            CALL PGSLW(2)
            CALL PGSCI (1)  !white needs to be reset after PGLAB
            CALL PGSCH(1.)
            CALL PGENV(x1-.5,x2+.5,y1-.5,y2+.5,0,1) !draw the axes
            CALL PGSTBG (0)     !puts opaque bgr to text (erase previous one)
            CALL PGLAB(c_x, c_y,plot_title)  !put the axis labels          
            CALL PGQCIR(C1, C2)
            NC = MAX(0, C2-C1+1)
            BRIGHT = 0.5
            CONTRA  = 1.0
            CALL PALETT(6, CONTRA, BRIGHT)
            CALL PGSLW (2)          ! *** write the domain legend
            CALL PGSCI(1)
            CALL PGMTXT ( 'RV',.5,1.-.05, .0,trim(dom_name))
            CALL PGSCH(.8)
            do k=1,n_dom
              if(mask(k)==0) cycle
  !             CALL PGSCI(20+k)
              CALL PGSCI(k)
              CALL PGMTXT ( 'RV',.5,1.-.05*(k+1), .0,trim(dom_ind(k)))
            enddo
            do j=1,n_dom
              if(mask(j)==0)cycle
              displ_vect = .0
              where(i_dom_out(:,:,jj)==j)
                displ_vect(1,:,:) = displ_plot(1,:,:,jj)
                displ_vect(2,:,:) = displ_plot(2,:,:,jj)
              end where
!               CALL PGSCI(20+j)
              CALL PGSCI(j)
              if(jj<=n_atom) then
                CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,scale,nc,TR,blank)
              else
                CALL PGVECT(displ_vect(1,:,:),displ_vect(2,:,:),n_x,n_y,1,n_x,1,n_y,.2*scale,nc,TR,blank)
             endif              
            enddo
            if(jj==1) then          ! *** print header with program version & date_and_time with the 1st pane
              call date_and_time(c_date,c_time,c_zone,i_time)
              write(time_stamp,'(a,"/",a,"/",a," ",a,":",a,":",a)') c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
              x_plot = x1 +.7*(x2-x1)
              y_plot = y2+.12*(y2-y1)
              CALL PGSCI (1)  !white needs to be reset after PGLAB
              CALL PGSLW (1)          ! *** write the domain legend
              CALL PGSCH(.9)
              CALL PGTEXT (x_plot,y_plot,trim(mp_tool)//' '//trim(time_stamp)//'       '//trim(file_master))
              CALL PGSCH(1.)
            endif
          enddo
        endif
        
      else
        print *,'ERR: no plot option selected'
      endif


        call PGCLOS
        j_ps = 0      !make just one .PS
        cycle slice_loop
      endif     !make .PS (j_ps.eq.1)

! *** stop the movie?
      if(j_auto==1.and.i_shift==n_shift) then !stop/cont automatic show
        print *,prompt, 'Continue/stop? (1/0)'
        read(*,*) j_auto
        i_shift = 0
        if(j_auto==0) cycle slice_loop
      endif                


! *** what next?
      if(j_frame==0)then                !cycle layer/snapshot (1/2)
        j_slice=j_slice+j_shift
      else
        jt = jt+j_shift
        print *,'jt=',jt
      endif                                  
      i_shift = i_shift+1

      if(j_slice<1) j_slice = 1
      if(j_slice>n_slice) j_slice = 1

      if(jt<jt0) jt = jt0
      if(jt>jt_max) jt = jt0

      enddo slice_loop

      if((mode>4.and.mode<10).or.mode>16) exit atom_loop
    enddo atom_loop

    deallocate(hist_plot,displ_plot,displ_vect,i_dom_out,ind)
  enddo plot_loop
!
! Close the device and exit.
!
888      CALL PGEND
  deallocate(displ_field,displ_norm,i_dom)

end program mp_dplot56

  
! **** string conversion to all upper case
!     
subroutine up_case (string)

  character(*), intent(inout)  :: string
  integer                      :: j, nc
  character(len=26), parameter  :: lower = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter  :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  do j = 1, len(string)
    nc = index(lower, string(j:j))
    if (nc > 0) string(j:j) = upper(nc:nc)
  end do

end subroutine up_case       


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
!        -- JK rainbow
         CALL PGCTAB(JKL, JKR, JKG, JKB, 9, CONTRA, BRIGHT)
      END IF
      END
      
   
  
     
