      
program mp_lbin56
 
! ************************************************************************************* 
! ***** 
! *****  %%%%%%%%%%%%%%%%   		  program MP_LBIN 1.56   		 %%%%%%%%%%%%%%%%%%%%%%%%%% 
! ***** 
! *****   converts MD trajectory data (LAMMPS or equivalent) to the MP_TOOLS binary form 
! ***** 
!**********					Copyright (C) 2022  Jiri Kulda, Grenoble/Prague          ********** 
! ************************************************************************************* 
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
! *****	Reads LAMMPS dump files and records binary data files (based on MB_LBINARY44.F) 
! ***** 
! ***** Ver. 1.44 - identifies atoms according to the fractional positions in the PAR file  
! *****						- abandons the "cryst_structure" TYPE 
! *****						- distinguishes low-T and hi-T behaviour of cores and shells:  
! *****									low-T strongly bound, k_B*T/2 = E_kin = E_kin_C + E_kin_S 
! *****									hi-T independent movement, k_B*T/2 = E_kin = E_kin_C = E_kin_S 
! ***** Ver. 1.52 - adaptation to new .PAR file  
! *****						- atoms are identified according to their atom_name (2nd column in .PAR) 
! ***** Ver. 1.53 - uses new semi-sequential binary file format  
! *****						- record length 4096, all data 32bit length 
! *****						- at_mass_in is saved as at_veloc(4,j,k) 
! *****						- at_charge_in is saved as at_pos(4,j,k) 
! *****           - CELL, FAST and BULK data type 
! ***** Ver. 1.54 - CELL and BULK data type revised (FAST discontinued, BULK may now use N_ROW) 
! *****           - using NAMELIST for parameters and headers 
! ***** 
! ***** Ver. 1.55 - bug fixes 
! ***** 
! ***** Ver. 1.56 - handles triclinic unit cell geometry
! ***** 
! ***** reads ASCII text output from LAMMPS(Sandia) MD simmulations (trajectory file) 
! ***** writes binary direct-access files for each snapshot 
! ***** depending on settings both CORE and SHELL records are read and recorded 
! ***** file name convention: 
! *****     <master_file>_n<snapshot number>.dat 
! *****     example: H05_10K_n0001.dat from H05_10K.txt 
! ***** 
! ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer) 
! ***** record numbers are directly related to cell positions and atom types 
! *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1) 
! *****						1 going for the header record 
! ***** 
! ***** atom positions are converted to and recorded in reduced lattice coordinates (x)  
! ***** 

  real,parameter    :: pi = 3.14159
  real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
  integer,parameter :: l_rec  =  1024		    !record length in real(4)

  logical :: found,found_txt,t_single,at_list
  character(4) :: at_name_in,col(32),version
  character(10)	 :: prompt,space = '          '
  character(10) :: sim_style,c_date,c_time,c_zone,ext,number,items(4),pg_out
  character(40)  :: subst_name,string,string2,section,mp_tool
  character(128) :: line,cwd_path,data_path,time_stamp,rec_str
  character(128) :: file_dat,file_trajectory,file_master,file_inp,file_log
  character(l_rec) :: header_record
  
  integer ::  at_no,at_ind_in(3),at_ind_base(3),i_dom,n_dom,i_dom_rec,n_at_cell
  integer ::  j_at,j_yes,j_shell,nskip,nfile_min,nfile_max,nfile_step,i_time(8),j_proc,j_mult,j_test
  integer ::  ios,ierr,i,j,k,ii,jj,jl,jat,jhead,j_step,j_label,n_label,n_site,j_eq,n_tstep0,n_tstep,n_items
  integer ::  inrec,jrec,nrec,i_rec,l_rec4,ifile,n_tot,n_tot_in,nsuper,nrow,nlayer,ind,ind_rec
  integer ::  sc_c1,sc_c2,sc_m,j_data,j_struct,j_verb,j_shrec,j_ext,j_time,j_basis,j_centred
  integer ::  ind_id,ind_type,ind_atom,ind_pos,ind_vel,ind_mass,ind_charge,ind_force,n_col
  integer ::  nt_min,nt_step,nt_max,n_save_min,i_save,i_traj		

  real :: at_mass_in,at_mass_in2,at_charge_in,at_charge_in2,item_value(4),sc_r
  real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force_in(3),at_force_in2(3),at_base_shift(3),a_scale
  real ::	at_pos_centre(3),a_cell(3,3),a_cell_1(3,3),a_cell_inv(3,3),a_cell_par(9),a_cell_lo(3),a_cell_hi(3),a_cell_half(3),cell_par,atp,b_coh
  real :: t_dump0,t0,t1,t2,dt,t_step,pos_inp(3),xy(3)
  real :: temp_par,temp_r_c,temp_r_s,eps_x
  
  character(4),allocatable :: at_name(:),at_label(:)
  character(16),allocatable :: data_line(:)
  integer,allocatable ::  ind_l(:),i_site(:,:),ind_at(:)
  integer,allocatable,target :: i_series(:),at_ind_out(:)
  integer,pointer ::  jr(:)
  real,allocatable ::  e_kin(:,:),e_kin_s(:,:),at_base_in(:,:),at_base(:,:),at_occup(:),at_mass_in_c(:),at_mass_in_s(:)

! **** the following variables MUST have the following 32bit sizes or multiples because of alignement in the binary output file 
!
  character(4),allocatable :: at_name_par(:),at_name_out(:)
  integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)

  real(4),allocatable ::	at_pos_c(:,:),at_veloc_c(:,:),at_force_c(:,:),at_occup_r(:)
  real(4),allocatable ::	at_pos_s(:,:),at_veloc_s(:,:),at_force_s(:,:)

  character(16)  :: sim_type,input_method,dat_type,dat_origin,dat_source,file_par,filter_name,pos_units
  integer(4)     :: n_row(3),n_atom,n_eq,j_shell_out,n_traj,n_cond,idum,n_rec,n_head,n_head_in1,n_head_in2
  real(4)        :: rec_zero(l_rec),t_ms,t_dump,filter_fwhm,a_par(3),angle(3),temp

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
  namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
  namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
 
  namelist /mp_gen/ j_verb,j_proc       
                      !general rule: namelists of tools should only contain their local parameters
                      !what is of global interest they should pass into data_header
  namelist /mp_bin/ subst_name,sim_type,dat_type,input_method,pos_units,data_path,ext,rec_str,&
 &					j_mult,n_head_in1,n_head_in2,n_tot_in,n_atom,n_row,j_basis,j_centred,j_test,j_shrec,&
 &          a_cell_par,a_scale,at_base_shift,eps_x,temp_par,t_step

!     namelist /mp_sqom/ n_int,s_trig,j_oneph,j_qsq
!     namelist /mp_pdf/ n_pdf,pdf_step,j_gauss,n_h,j_weight,j_smooth,n_corr

  data rec_zero/l_rec*.0/

! 
! ********************* Initialization *******************************      
  version = '1.56'
  prompt = 'MP_LBIN>  '
  mp_tool = 'MP_LBIN '//version

  print *,'*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  print *
  dat_source = 'MP_TOOLS'

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
  write(9,*) trim(time_stamp),'  ',mp_tool,'  ',trim(cwd_path) 
  write(9,*) 

! *** diverse initialisations
  l_rec4 = l_rec/4
  j_test = 1
  n_head = 3
  n_head_in1 = 0
  n_head_in2 = 0
  rec_str = ''
  idum = 0
  j_label = 0
  n_site = 0
  n_eq = 1								!later introduce sites, basis etc.
  i_dom = 0
  n_cond = 2               !orthorhombic periodic bound_cond by default
  n_tot_in = 0
  at_list = .true.
  t_ms = .0
  t_step = .0
  t_dump = .0
  t_dump0 = .0
  temp_par = .0
  at_base_shift = .0
  a_scale = 1.
  filter_name = ''
  filter_fwhm = .0
  pos_units = 'ANGSTROM'
  
! *** read auxiliary file <file_par.par> with structure parameters, atom names and further info
  print *,prompt, 'Parameter file name (.par will be added)'
  read(*,*) file_par
  file_inp = trim(file_par)//'.par'

  open(4,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,space, 'File ',trim(file_inp),' not found! Stop execution.'
    stop
  endif

  write(9,*) 'Reading parameter file:  ',trim(file_inp)

! *** Read the mp_gen and mp_bin namelists       
  read(4,nml=mp_gen) 
  rewind(4)

  data_path = './data/'     !the default location

  read(4,nml=mp_bin) 
  call up_case(sim_type)
  call up_case(dat_type)
  call up_case(input_method)
  call up_case(pos_units)
  call down_case(rec_str)

  if(ext=='ext'.or.ext=='EXT') ext=''
  if(ext/=''.and.index(ext,'.')==0) ext='.'//ext
  
  if(dat_type=='LAMMPS') then
    dat_origin = 'LAMMPS'
    if(pos_units/='ANGSTROM') then
      print *,space, 'Position units ', trim(pos_units),' not implemented with LAMMPS data type, ANGSTROM needed or use GENERAL data type! '
      stop
    endif
  elseif(dat_type=='DL_POLY') then
    print *,space, 'For input of DL_POLY data use MP_DBIN'
    stop
  elseif(dat_type=='GENERAL') then
    print *,space, 'Position units:', trim(pos_units)
  endif

  
  if(dat_type=='GENERAL') then    !it is assumed that here the box is only defined at the very beginning of the simulations and is given at the start in the .PAR file
    a_cell = transpose(reshape(a_cell_par,(/3,3/)))
  endif 
  
  allocate(at_name_par(n_atom),at_label(n_atom),at_base_in(n_atom,3),at_base(n_atom,3))
  allocate(ind_l(n_atom),i_site(n_atom,n_atom),at_occup(n_atom))
  allocate(at_mass_in_c(n_atom),at_mass_in_s(n_atom))

  ind_l = 0
  nsuper = n_row(1)*n_row(2)*n_row(3)
  nlayer = n_row(1)*n_row(2)						!only to be used for record number calculation for CELL data

! *** Read the atom positions       
  section = 'atoms'
  rewind(4)
  do
    read(4,'(a)',iostat=ios) string 
    if(ios/=0) then
      print *,space, 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
      stop
    endif
    if(string(1:5).eq.section) exit	
  enddo

  do j=1,n_atom
    read(4,*) string,at_name_par(j),at_base_in(j,:),at_occup(j)
    found = .false.
    do i=1,j_label
      if(trim(string).eq.at_label(i)) then
        found = .true.
        ind_l(i) = ind_l(i)+1
        i_site(i,ind_l(i)) = j
        exit
      endif
    enddo
    if(found) cycle
    j_label = j_label+1
    at_label(j_label) = trim(string)
    ind_l(i) = 1
    i_site(j_label,ind_l(i)) = j
  enddo
  close(4)		
  
  do j=1,n_atom
    at_base(j,:) = at_base_in(j,:)+at_base_shift
  enddo

  string = subst_name
  print *,prompt, 'Substance name (confirm, &append or replace): ', string
  read(*,*) string
  string = trim(adjustl(string))
  if(string/=subst_name) then
    if(string(1:1)=='&') then
      subst_name = trim(subst_name)//string(2:)
    else
      subst_name = string
    endif
  endif

  print *
  print *,space, 'Substance name: ', subst_name
  print *,space, 'Sim_type, dat_type, input method: ',sim_type,dat_type,input_method		
			
  n_label = j_label
  print *, space, 'Atom labels (chemical species) ',n_label
  do j=1,n_label
    print *,space, ind_l(j),at_label(j),(i_site(j,i),i=1,ind_l(j))
  enddo	

  print *,space, 'Atom types: '	  
  do j=1,n_atom
    if(input_method=='CELL') then
      print *,space, j,at_name_par(j),at_base_in(j,:)
    else
      print *,space, j,at_name_par(j)
    endif
  enddo

  at_charge_in = .0
  j_shell = 0
  jhead = 0
  nskip = 20 ! number of records to test
  j_time = 0				!1 indicates presence of TIME
  n_items = 0
  file_master = '' 
! ***********************************
! n_row_par 32 32 32 & bounds:   standard box  cell & bulk periodic         !bounds may follow on N_ROW line in .PAR
! n_row_par  1  1  1 & bounds:   standard box  bulk aperiodic              
! n_row_par 32 32 32 & no_bounds: standard box  cell & bulk periodic    r.l.u.
! n_row_par  1  1  1 & no_bounds: min_max box  bulk aperiodic          r.l.u. 
! **********************************


  print *,prompt, 'Read snapshots number: min, max'
  read(*,*) nfile_min, nfile_max 
  nfile_step = 1

  write(9,*) 'Read snapshots number: min, step, max',nfile_min, nfile_step,nfile_max
  print *,prompt, 'Saved snapshot numbers start:'
  read(*,*) n_save_min 
  write(9,*) 'Saved snapshot numbers start:',n_save_min   
   
! *** input cycle over MD trajectory files
  CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)			
  i_save = n_save_min

  print *,prompt, 'Master for MD input filename: '
  read(*,*) file_master 
!!			trajectory_inp = (file_master/='') 

!!			if(trajectory_inp) then 
!!				i = len_trim(file_master) 
!!				if(file_master(i:i)=='_') then 
!!					file_master_out = file_master(1:i-1) 
!!				else 
!!					file_master_out = file_master(1:i) 
!!				endif

  print *,prompt, 'Read input files number (0 0 0 no numbers): min,step,max'
  read(*,*) nt_min,nt_step,nt_max
  t_single = (nt_min==0.and.nt_max==0)
  if(nt_step==0) nt_step=1

!! *** treat 1st frame of 1st history file to get info
  i_traj = nt_min

  if(t_single) then
    file_trajectory = file_master
  else 
!!					if(i_traj>=1.and.i_traj<=9)    write(number,'(i1.1)') i_traj 
!!					if(i_traj>=10.and.i_traj<=99)  write(number,'(i2.2)') i_traj 
!!					if(i_traj>=100.and.i_traj<=999)write(number,'(i3.3)') i_traj 
!!					if(i_traj>=1000.and.i_traj<=9999)write(number,'(i4.4)') i_traj
    write(number,'(i8)') i_traj
    number = trim(adjustl(number))
  endif


!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  GENERAL  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  if(dat_type=='GENERAL') then      !parse input record structure and position for data input in 1st file to test for SHELL presence
    if(t_single) then
      file_inp = trim(data_path)//trim(adjustl(file_master))//trim(ext)	
    else
      file_trajectory = trim(file_master)//trim(number)
      file_inp = trim(data_path)//trim(adjustl(file_trajectory))//trim(ext)
    endif						
    print *,space, '1st input file to read:',file_inp
    write(9,*) '1st input file to read:',file_inp

! *** parse the data record description 
! *** covered items:     id type element mass q x y z vx vy vz fx fy fz nn 
! *** obligatory items:  x y z  
! *** NN is a placeholder for input positions to be skipped, hence no ind_nn 
! *** GENERAL method needs al least ELEMENT X Y Z to be defined

    call down_case(rec_str)		!for LAMMPS rec_str comes from data, for GENERAL from .PAR
    n_col = 1
    do 
      read(rec_str,*,iostat=ios) (col(i),i=1,n_col)
      if(ios/=0) exit
      n_col = n_col+1
    enddo      
    n_col = n_col-1
   
    i = 0
    ind_id = i
    ind_type = i
    ind_atom = i
    ind_mass = i
    ind_charge = i
    ind_pos = i
    ind_vel = i			
    ind_force = i	

    do i=1,n_col
      select case(col(i))
        case('id') 
          ind_id = i
        case('type') 
          ind_type = i
        case('elem') 
          ind_atom = i
        case('mass') 
          ind_mass = i
        case('q') 
          ind_charge = i
        case('x') 
          ind_pos = i
        case('vx') 
          ind_vel = i			
        case('fx') 
          ind_force = i			
      end select
    enddo  
    if(ind_pos==0) then
      print *,space, 'The keys X Y Z are missing in the data.'
      stop
    endif
    
    at_list = ind_atom>0.or.ind_type>0     !at_list=.false. means no explicit atom identification

    n_traj = 0
    if(ind_vel/=0) n_traj = 1
    if(ind_force/=0) n_traj = 2
  
    allocate(data_line(n_col))
     
! *** open 1st frame and position to data read to test for shell presence  !SHELLs only with LAMMPS data type 
!			open (1,file=file_inp,action='read',status ='old') 
!			do i=1,n_head_in1 
!				read(1,*)                                             
!			enddo 
!			do i=1,n_head_in2 
!				read(1,*) 
!			enddo 
!
endif !'GENERAL'


!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  LAMMPS  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! *** look for the first input file with MD trajectory for data type LAMMPS   
  if(dat_type=='LAMMPS') then			!goes down to line 699 !retrieve the header info from 1st file and position for reading 1st data line

    file_inp = trim(data_path)//trim(adjustl(file_master))//'.'//trim(adjustl(number))							
    open (1,file=file_inp,action='read',status ='old')
    write(9,*) '1st snapshot file to read:',file_inp

! *** look for keyword 'TIMESTEP'   
    do i=1,nskip
      read(1,'(a)') line
      call up_case(line)
      if(index(line,'TIMESTEP')/=0) then
        sim_type = 'TIMESTEP'
        read(1,*) n_tstep0
        exit
      endif
      if(i==nskip) then
        print *,space, 'Mandatory TIMESTEP info missing! Check SIM_TYPE in .PAR'
        stop
      endif
    enddo
    rewind(1)

    do i=1,nskip
      read(1,'(a)') line
      call up_case(line)
      if((index(line,'TIME')/=0.and.index(line,'TIMESTEP')==0)) then
        read(1,*) t_dump0
        t_ms = t_dump0/n_tstep0
        n_items = n_items+1
        j_time = 1
        exit
      else
        print *,space,'Input data do not contain explicit TIME information!'
        print *,prompt,'Type in the microstep T_MS [ps] (frequent value = .0002)'
        read(*,*) t_ms
        t_dump0 = t_ms*n_tstep0
        exit
      endif
    enddo
    t_step = t_ms*nt_step
    print *,space,'The frequency scale will rely on T_MS=',t_ms
    print *,space,'If strange, modify this value and restart the whole treatment.'
    write(9,*)'Input data do not contain explicit TIME information!'
    write(9,*)'The frequency scale will rely on T_MS=',t_ms

    rewind(1)
    
    do i=1,nskip
      read(1,'(a)') line
      call up_case(line)
     
      if(index(line,'NUMBER')/=0) then
        read(1,*) n_tot_in
        n_items = n_items+1
        exit
      endif
    enddo
    rewind(1)

!      a_cell = .0

    do i=1,nskip
      read(1,'(a)') line
      call up_case(line)
      if(index(line,'BOUNDS')/=0) then
        n_head_in1 = i-1
        if(index(line,'PP')/=0) then              !periodic boundary conditions
          n_cond = 1
        else
          n_cond = 0                      !aperiodic boundary conditions
        endif


!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!           
!  contrary to the GENERAL data for the LAMMPS case the cell parameters will be retrieved on the fly 
!  the following code is just for initial INFO 

       a_cell = .0
       if(index(line,'XY')/=0) then       !the box is triclinic 
          do k=1,3 
            read(1,*) a_cell_lo(k),a_cell_hi(k),xy(k)   !these are bounding box dimensions 
          enddo 
           
          a_cell(1,1) = a_cell_hi(1)-max(0.0,xy(1),xy(2),xy(1)+xy(2))-(a_cell_lo(1)-min(0.0,xy(1),xy(2),xy(1)+xy(2))) 
          a_cell(2,1) = xy(1) 
          a_cell(2,2) = a_cell_hi(2)-max(0.0,xy(3))-(a_cell_lo(2)-min(0.0,xy(3))) 
          a_cell(3,1) = xy(2) 
          a_cell(3,2) = xy(3) 
          a_cell(3,3) = a_cell_hi(3)-a_cell_lo(3) 

        else                               !the box is orthorombic 
          do k=1,3 
            read(1,*) a_cell_lo(k),a_cell_hi(k) 
            a_cell(k,k) = a_cell_hi(k)-a_cell_lo(k) 
          enddo 
        endif 
        n_items = n_items+1
        exit
      endif
    enddo

    rewind(1)
    do i=1,nskip
      read(1,'(a)') line
      call up_case(line)
      if(index(line,': ATOMS').gt.0) exit
      if(i.eq.nskip) then
        print *,space, 'ATOMS: data not found'
        stop
      endif
    enddo
    rec_str = line(index(line,'S')+1:) 
   
! *** parse the data record description 
! *** covered items:     id type element mass q x y z vx vy vz fx fy fz nn 
! *** obligatory items:  (type OR element) x y z  
! *** NN is a placeholder for input positions to be skipped, hence no ind_nn 
! *** GENERAL method needs ELEMENT X Y Z to be defined

    call down_case(rec_str)		!for LAMMPS rec_str comes from data, for GENERAL from .PAR
    n_col = 1
    do 
      read(rec_str,*,iostat=ios) (col(i),i=1,n_col)
      if(ios/=0) exit
      n_col = n_col+1
    enddo      
    n_col = n_col-1
   
    i = 0
    ind_id = i
    ind_type = i
    ind_atom = i
    ind_mass = i
    ind_charge = i
    ind_pos = i
    ind_vel = i			
    ind_force = i	

    do i=1,n_col
      select case(col(i))
        case('id') 
          ind_id = i
        case('type') 
          ind_type = i
        case('elem') 
          ind_atom = i
        case('mass') 
          ind_mass = i
        case('q') 
          ind_charge = i
        case('x') 
          ind_pos = i
        case('vx') 
          ind_vel = i			
        case('fx') 
          ind_force = i			
      end select
    enddo 

    if(ind_pos==0) then
      print *,space, 'The keys X Y Z are missing in the data.'
      stop
    endif

    n_traj = 0
    if(ind_vel/=0) n_traj = 1
    if(ind_force/=0) n_traj = 2
  
    allocate(data_line(n_col))

  if(n_tot_in==0) then
    print *,prompt, 'Info on total number of atoms (input lines) missing, type in N_TOT: '
    read(*,*) n_tot_in
  endif

   
! *** Test possible core/shell presence, assuming all atoms would be the same

!  typical data: 
!  1 1 Na1_c 20.9 0.167166 0.254627 0.115284 3.148 -3.06563 0.629693  
!  2 9 Na1_s 2.09 0.165161 0.250872 0.113053 2.9716 -2.99374 0.691807 
 
  if(ind_atom/=0) then
    read(1,*)  data_line  
    string = 	data_line(ind_atom)
    read(1,*)  data_line  
    string2 = 	data_line(ind_atom)
    if((index(string,'_c').gt.1.or.index(string,'_C').gt.1).and.(index(string2,'_s').gt.1.or.index(string2,'_S').gt.1))	j_shell = 1		
  else
    j_shell = 0
  endif
  
  close(1)
  
  if(j_shell==1) then
    print *,space,'Core & shell data found'
    if(j_shrec==0) print *,space,'Shell data NOT to be recorded (change this in .PAR)'
    if(j_shrec==1) print *,space,'Shell data WILL be recorded'
    if(j_shrec==1) write(9,*)'Shell data WILL be recorded'
  endif
  j_shell_out = j_shell*j_shrec

     
!!			if(j_shell.eq.0) then 
!!				sim_style = 'atom' 
!!			elseif(j_shell.eq.1) then 
!!				sim_style = 'core/shell' 
!!			else 
!!				sim_style = 'undefined!' 
!!			endif

  endif !'LAMMPS'

!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCC  END LAMMPS  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! *** produce informative values on A_PAR and ANGLES, for GENERAL data they will serve throughout ....
! *** calculate angles of unit cell and the inverse of a_cell                    
    do j=1,3
      a_cell(j,:) = a_cell(j,:)/n_row(j)        !for n_row>1 a_cell becomes unit cell
      a_par(j) = norm2(a_cell(j,:))
      if(nsuper==1) a_cell(j,:) = a_cell(j,:)/a_par(j)         !effective a_par=1 and Q will be in [A-1]
    enddo 

    angle(1) = dot_product(a_cell(2,:),a_cell(3,:))/(a_par(2)*a_par(3))
    angle(2) = dot_product(a_cell(1,:),a_cell(3,:))/(a_par(1)*a_par(3))
    angle(3) = dot_product(a_cell(1,:),a_cell(2,:))/(a_par(1)*a_par(2))
    angle = acos(angle)
    
    if(j_verb==1) then
      print *,space, 'angle_rad',angle
      print *,space, 'angle_deg',angle*180./pi
    endif

    a_cell_1 = a_cell
    call gjinv(a_cell_1,3,3,a_cell_inv,3,ierr)       !!! a_cell would get destroyed here !!!
    if(ierr==1) then
      print *,space, 'Singular cell vector matrix, check your HISTORY file!'
      stop
    endif
     
  print *,space, 'Simulation type = ',sim_type
  print *,space, 'Trajectory recording mode =', n_traj
  if(sim_type == 'TIMESTEP') print *,space, 'Trajectory time start, step [ps]:',t_dump0,t_step
  print *,space, 'Simulation box:',n_row
  print *,space, 'Lattice parameter:',a_par
  print *,space, 'Angles:           ',angle

  write(9,*) 'sim_type = ',sim_type
  write(9,*) 'Substance:',subst_name,'n_atom, n_eq,n_row:',n_atom,n_eq,n_row
  if(sim_type == 'TIMESTEP') write(9,*) 'Trajectory time start, step [ps]:',t_dump0,t_step
  write(9,*) 'Simulation box:',n_row
  write(9,*) 'Lattice parameter:',a_par
  write(9,*) 'Angles:           ',angle
  
  if(at_list.eqv..false.) then
    print *,space, 'Atom species not identified explicitly, using their order in the list'
    write(9,*) 'Atom species not identified explicitly, using their order in the list'
  endif

  if(input_method=='CELL'.and.j_basis==0) then               
    print *,space, 'Atoms types not identified explicitly,'
    print *,space, 'trying to identify them by fractional positions'
    print *,space, '      tolerance range:',eps_x,' check whether the whole box isn''t shifted!'
    write(9,*) 'Atoms types not identified explicitly,trying to identify them by fractional positions,'
    write(9,*) '      tolerance range:',eps_x
  endif
   
!!      print *,space, 'domain analysis:',n_dom,' 0=no, 1=[100], 2=[110], 3=[111]'
  print *,prompt, 'go ahead? (1/0)'
  read(*,*) j_yes
  if (j_yes.eq.0) stop

  print *,space, '  (working)'
   
! *** handle the supercell parameters and the origin of the supercell coordinate system
  if(input_method=='CELL'.or.input_method=='FAST') then
    if(nsuper==1) then
      print *,space, 'CELL and FAST input methods not compatible with N_ROW = 1; use BULK'
      stop
    endif
    n_tot = n_atom*nsuper
  else						!input_method=='BULK'
    n_tot = n_tot_in
  endif
  
  at_ind_base = n_row/2+1

!     if(input_method=='CELL') then          
!			n_tot = n_atom*n_row(1)*n_row(2)*n_row(3) 
!    	  at_ind_base = n_row/2+1 
!		else						!input_method=='BULK' 
!			n_tot = n_tot_in 
!		endif

  allocate(at_name(n_tot),SOURCE='    ')
  allocate(e_kin(n_atom,3),at_occup_r(n_atom))
  allocate(nsuper_r(n_atom),SOURCE=0) 
  if(j_shell.eq.1)allocate(e_kin_s(n_atom,3),SOURCE=0.0)
  allocate(at_ind(4,n_tot),i_series(n_tot))
  i_series = (/ (i, i = 1, n_tot) /)


!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  INPUT  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! *** Now ready to cycle over trajectory files, each snapshot to be saved in a separate binary file

  CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)

! *** prepare atom indexing for list-directed data files
  if(.not.at_list) then
      if(sum(int(at_occup))/=n_tot_in) then
      print *,space, 'Atom occupancies not consistent with total number of atoms in input (N_TOT_IN)'
      print *,space, 'Check your data and .PAR file'
      stop
    endif
    ii = 0
    if(j_basis==1) then
      do j_at=1,n_atom
        do i=1,int(at_occup(j_at))
          ii = ii+1
          at_ind(1,ii) = ii     
          at_ind(4,ii) = j_at      !prefill atom_type j into at_ind
        enddo
      enddo
    else
      jl = 0
      do j_label=1,n_label
        do j=1,ind_l(j_label)
          jl = jl+1
          do i=1,int(at_occup(jl))
            ii = ii+1
            at_ind(1,ii) = ii     
            at_ind(4,ii) = j_label      !prefill atom_label jl into at_ind
          enddo
        enddo
      enddo
    endif
  endif

  ifile = nfile_min
  nsuper_r = 0        !n_super will be accumulated during the read of the 1st snapshot file

  trajectory_loop: do i_traj=nt_min,nt_max,nt_step 
    if(t_single) then
      file_trajectory = file_master
    else 
!!					if(i_traj>=1.and.i_traj<=9)    write(number,'(i1.1)') i_traj 
!!					if(i_traj>=10.and.i_traj<=99)  write(number,'(i2.2)') i_traj 
!!					if(i_traj>=100.and.i_traj<=999)write(number,'(i3.3)') i_traj 
!!					if(i_traj>=1000.and.i_traj<=9999)write(number,'(i4.4)') i_traj 
!!					if(i_traj>=10000.and.i_traj<=99999)write(number,'(i5.5)') i_traj
      write(number,'(i8)') i_traj
      number = trim(adjustl(number))
    endif

    if(dat_type=='LAMMPS') then
      file_inp = trim(data_path)//trim(adjustl(file_master))//'.'//trim(adjustl(number))
    else            !GENERAL
      if(t_single) then
        file_inp = trim(data_path)//trim(adjustl(file_master))//trim(ext)	
      else
        file_trajectory = trim(file_master)//trim(number)
        file_inp = trim(data_path)//trim(adjustl(file_trajectory))//trim(ext)
      endif						
    endif
  
    if(j_verb==1.or.i_traj==1) print *,space, 'Reading trajectory file(s):  ',trim(file_inp)

!!				t1 = .0 
!!				do  
!!					inquire (file=file_inp,exist=found) 
!!					if(.not.found) then 
!!						if(t1==0.) print *,'File ',trim(file_inp),' not found! Waiting for it ...'			!if(t1==0.)  
!!						CALL SLEEP(10) 
!!						t1 = t1+10. 
!!					else 
!!						exit 
!!					endif 
!!				enddo

    open (1,file=file_inp,action='read',status ='old',iostat=ios)
    if(ios.ne.0) then
      print *,space, "Can't open file ",trim(file_inp),'!'
      stop
    endif				

! *** open the input file and position it to data read the 1st frame 
    open (1,file=file_inp,action='read',status ='old')
    do i=1,n_head_in1
      read(1,*)
    enddo

! *** take fresh lattice parameters for each LAMMPS snapshot, GENERAL uses the a_cell input from .PAR
    if(dat_type=='LAMMPS') then
      read(1,'(a)') line
      call up_case(line)
      if(index(line,'BOUNDS')==0) then
        print *,space, 'Inconsistent data files: BOUNDS not in place!' 
        stop           
      endif
    
      a_cell = .0
      if(index(line,'XY')/=0) then       !the box is triclinic
        do k=1,3
          read(1,*) a_cell_lo(k),a_cell_hi(k),xy(k)   !these are bounding box dimensions
        enddo
        a_cell(1,1) = a_cell_hi(1)-max(0.0,xy(1),xy(2),xy(1)+xy(2))-(a_cell_lo(1)-min(0.0,xy(1),xy(2),xy(1)+xy(2)))
        a_cell(2,1) = xy(1)
        a_cell(2,2) = a_cell_hi(2)-max(0.0,xy(3))-(a_cell_lo(2)-min(0.0,xy(3)))
        a_cell(3,1) = xy(2)
        a_cell(3,2) = xy(3)
        a_cell(3,3) = a_cell_hi(3)-a_cell_lo(3)
       else                               !the box is orthorombic
        do k=1,3
          read(1,*) a_cell_lo(k),a_cell_hi(k)
          a_cell(k,k) = a_cell_hi(k)-a_cell_lo(k)
        enddo
      endif
      read(1,*)           !line     skip the header line with ATOMS       
   
! *** calculate angles of unit cell and the inverse of a_cell                    
      do j=1,3
        a_cell_half(j) = .5*sum(a_cell(:,j))
        a_cell(j,:) = a_cell(j,:)/n_row(j)        !for n_row>1 a_cell becomes unit cell
        a_par(j) = norm2(a_cell(j,:))
        if(nsuper==1) a_cell(j,:) = a_cell(j,:)/a_par(j)         !effective a_par=1 and Q will be in [A-1]
      enddo 

      angle(1) = dot_product(a_cell(2,:),a_cell(3,:))/(a_par(2)*a_par(3))
      angle(2) = dot_product(a_cell(1,:),a_cell(3,:))/(a_par(1)*a_par(3))
      angle(3) = dot_product(a_cell(1,:),a_cell(2,:))/(a_par(1)*a_par(2))
      angle = acos(angle)
    
      if(j_verb==1) then
        print *,space, 'angle_rad',angle
        print *,space, 'angle_deg',angle*180./pi
      endif

      a_cell_1 = a_cell
      call gjinv(a_cell_1,3,3,a_cell_inv,3,ierr)       !!! a_cell would get destroyed here !!!
      if(ierr==1) then
        print *,space, 'Singular cell vector matrix, check your HISTORY file!'
        stop
      endif
  
      if(j_verb==1) then
        print *,space, 'a_cell'
        do k=1,3
          print *,space, a_cell(k,:)
        enddo
        print *,space, 'a_cell_inv'
        do k=1,3
          print *,space, a_cell_inv(k,:)
        enddo
      endif	
    endif      !LAMMPS
	
! *** cycle over snapshots, each to be saved in a separate binary file 

  frame_loop: do 				!we have to go through all the snapshots in the .TXT input

    if(ifile>nfile_max) exit trajectory_loop

! *** test end of trajectory file and position for data input
    if(n_head_in2>0) then
      do i=1,n_head_in2
        read(1,*,iostat=ios)
        if(ios<0) exit frame_loop
      enddo
    else
      read(1,*,iostat=ios) data_line
      if(ios<0) exit frame_loop
      backspace(1)
    endif

    string = ''
    allocate(at_pos_c(4,n_tot))
    if(n_traj>=1) allocate(at_veloc_c(4,n_tot),source=.0)
    if(n_traj==2) allocate(at_force_c(4,n_tot),source=.0)
    if(j_shell_out==1) then
      allocate(at_pos_s(4,n_tot))
      if(n_traj>=1) allocate(at_veloc_s(4,n_tot),source=.0)
      if(n_traj==2) allocate(at_force_s(4,n_tot),source=.0)
    endif

    e_kin = .0
    if(j_shell.eq.1) e_kin_s = .0

! ***** now read the data
    read_loop: do i=1,n_tot

      read(1,*,iostat=ios)  data_line   		
      if(ios<0)then
          exit read_loop
      elseif(ios>0) then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
      endif
      
      if(ind_id/=0)then
        read(data_line(ind_id),*,iostat=ios) inrec
        if(ios/=0)then
          print *,space, 'Data read problem: ios =',ios
          print *,prompt, 'record:',i,' data:',data_line
          read(*,*)
          cycle read_loop
        endif
      else
        inrec = i
      endif
      if(ind_type/=0)read(data_line(ind_type),*,iostat=ios) jl
      
      if(ind_atom/=0)string = data_line(ind_atom)
      if(ios/=0)then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
        cycle read_loop
      endif

      if(ind_mass/=0)read(data_line(ind_mass),*,iostat=ios) at_mass_in
      if(ios/=0)then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
      endif

      if(ind_charge/=0)read(data_line(ind_charge),*,iostat=ios) at_charge_in
      if(ios/=0)then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
        cycle read_loop
      endif

      read(data_line(ind_pos:ind_pos+2),*,iostat=ios) at_pos_in
      if(ios/=0)then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
        cycle read_loop
      endif
      at_pos_in = at_pos_in/a_scale								!apply a scaling to convert from exotic position units; normally a_scale = 1.
                 
      if(ind_vel/=0)read(data_line(ind_vel:ind_vel+2),*,iostat=ios) at_veloc_in
      if(ios/=0)then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
        cycle read_loop
      endif

      if(ind_force/=0)read(data_line(ind_force:ind_force+2),*,iostat=ios) at_force_in
      if(ios/=0)then
        print *,space, 'Data read problem: ios =',ios
        print *,prompt, 'record:',i,' data:',data_line
        read(*,*)
         cycle read_loop
      endif

        if(j_shell.eq.1) then
          ii = max(index(string,'_c'),index(string,'_C'))  					!supposed core comes first
          at_name_in = string(1:ii-1)
        else
          at_name_in = trim(string)
        endif 

       
! *** read also the SHELL line if needed 
      if(at_list) then            !for at_list=.false. the shell-test has no sense, possibly list shells as atoms ....
        if(j_shell.eq.1) then     !this is based on testing the presence of '_c' and '_s' in subsequent lines; 
          read(1,*)  data_line   		
          if(ios.ne.0.and.i.lt.n_tot) then
            print *,space, 'EOF or EOR at',i,jat,string 
            stop
          endif

          if(ind_mass/=0)read(data_line(ind_mass),*) at_mass_in2
          if(ind_charge/=0)read(data_line(ind_charge),*) at_charge_in2 
          read(data_line(ind_pos:ind_pos+2),*) at_pos_in2
          if(ind_vel/=0)read(data_line(ind_vel:ind_vel+2),*) at_veloc_in2
          if(ind_force/=0)read(data_line(ind_force:ind_force+2),*) at_force_in2

          string2 = data_line(ind_atom)
          ii = max(index(string2,'_s'),index(string2,'_S'))  
           
          if(string2(1:ii-1)/=trim(at_name_in)) then
            print *,space, 'inconsistent core/shell names',i,string,string2
            write(9,*) 'inconsistent core/shell names',i,string,string2,  'STOP!'
            stop
          endif
        endif   !j_shell

        if(ind_type==0) then            ! for LAMMPS itself should never occur - TYPE is one of the default identifiers
          if(j_basis==0) then
            do j=1,n_atom
              if (at_name_in.eq.at_label(j)) then     !look for label
                jl=j			!atom label found
                exit
              endif
            enddo
            if(jl.ne.j) then			!after "normal" loop exit
              print *,space, 'atom ',at_name_in,' not found in .PAR'
              stop
            endif
          elseif(j_basis==1)  then
!            jl = jl      !nothing to do
          elseif(j_basis==2)  then
            do j=1,n_atom
              if (at_name_in.eq.at_name_par(j)) then    !look for type
                jl=j			!atom label found
                exit
              endif
            enddo
            if(jl/=j) then			!after "normal" loop exit
              print *,space, 'atom ',at_name_in,' not found in .PAR'
              stop
            endif
          else 
            print *,space, 'ERROR: wrong J_BASIS',j_basis
            stop           
          endif
        endif  !ind_type
      else
        jl = at_ind(4,inrec)    ! nothing to do: at_inds are already defined, labels or types
      endif  !at_list

! *** transform the atom positions to pseudo-cubic r.l.u. coordinates 
!
      if(pos_units=='ANGSTROM') then
        if(j_centred==0) then
          at_pos_centre = a_cell_half
        else
          at_pos_centre = .0
        endif
        at_pos_in = at_pos_in-at_pos_centre
        at_pos_in = matmul(a_cell_inv,at_pos_in)
        if(j_shell_out==1) then                     !shells only with LAMMPS which is ANGSTROM
          if(j_centred==0) at_pos_in2 = at_pos_in2-at_pos_centre
          at_pos_in2 = matmul(a_cell_inv,at_pos_in2)
        endif       				
      elseif(pos_units=='LATTICE') then
        if(j_centred==0) then
          at_pos_centre = n_row/2
        else
          at_pos_centre = .0
        endif
        at_pos_in = at_pos_in-at_pos_centre   			  
      elseif(pos_units=='BOX') then
        if(j_centred==0) then
          at_pos_centre = .5
        else
          at_pos_centre = .0
        endif
        at_pos_in = at_pos_in-at_pos_centre   			  
        if(nsuper/=1) at_pos_in = at_pos_in*n_row      !bring at_pos_in from the (-.5,.5) range to the r.l.u. box range
      endif

! *** index BULK data 
      if(input_method=='BULK') then            !jl identifes label (chem species), jat basis position in CELL
        jat = jl				          
        at_ind(1,inrec) = inrec
        at_ind(2:3,inrec) = 0
        at_ind(4,inrec) = jat
        jrec = inrec

! *** index CELL data 
      elseif(input_method=='CELL') then
        if(j_basis==1) then
          jat = jl
        else		                    !trying to identify atoms by their position
          if(ind_l(jl).eq.1) then
            jat = i_site(jl,1)	    !atom name and fractional position is identified
          else
            do ii=1,ind_l(jl)		    !multiple basis positions correspond to the at_label - LAMMPS: identification by position needed in each frame 
              pos_inp = at_pos_in-at_base_in(i_site(jl,ii),:)
              jat = i_site(jl,ii)
              if(maxval(abs(pos_inp-anint(pos_inp))).le.eps_x) exit !atom found
            enddo
            if(maxval(abs(pos_inp-anint(pos_inp))).gt.eps_x) then
              print *,space, 'Identification by position not succeeded: frame, record, atom ',ifile,inrec,at_label(jl)
              print *,space, 'Input position ',at_pos_in
              print *,space, 'Candidates (JAT, ATOM, BASIS POSITION, MAX_DIFF,EPS_X): '
              do ii=1,ind_l(jl)
                pos_inp = at_pos_in-at_base_in(i_site(jl,ii),:)
                print *,space, i_site(jl,ii),at_name_par(i_site(jl,ii)),at_base_in(i_site(jl,ii),:),maxval(abs(pos_inp-anint(pos_inp))),eps_x
              enddo
              print *,prompt, 'Type JAT and confirm/modify AT_POS_IN:'
              read(*,*) jat,at_pos_in
              print *,space, 'Other possible solutions:'
              print *,space, '  1/check the ATOMS basis, 2/ try to slightly increase EPS_X, 3/ use the BULK input method, 4/edit the trajectory file '
              print *,space, '      (working ...) '
            endif
          endif						
        endif

        if(jat<1.or.jat>n_atom) then
          print *,space, 'Atom ID wrong:',jat,' not compatible with N_ATOM =',n_atom
        endif
        
!! *** calculate cell indices ix,iy,iz from atomic positions shifted to cell origin

        at_ind_in = anint(at_pos_in-at_base_in(jat,:))+at_ind_base		!they will serve as pointers to the right order of atom records
        at_pos_in = at_pos_in+at_base_shift			!now the supercell will be centred & basis positions origin at 0
        if(j_shell_out.eq.1) at_pos_in2 = at_pos_in2+at_base_shift

        do k=1,3
          if(at_ind_in(k)==0) then
            at_ind_in(k) = at_ind_in(k)+n_row(k)
          endif
          if(at_ind_in(k)==n_row(k)+1) then
            at_ind_in(k) = at_ind_in(k)-n_row(k)
          endif
        enddo
      
        jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1)

        if(jrec>n_tot_in) then
          print *,space, 'JREC wrong:',jrec,' > ',n_tot,' check n_tot_in, n_row and j_centred in the .PAR'
          print *,space, 'jrec,at_ind_in',jrec,at_ind_in
          print *,space, 'at_pos_in',at_pos_in,at_pos_in2
          print *,space, 'n_tot,nsuper,nlayer,n_row',n_tot,nsuper,nlayer,n_row
          stop
        endif

        at_ind(1:3,jrec) = at_ind_in
        at_ind(4,jrec) = jat					
      endif		!input_method CELL

! *** enforce atom positions within the box limits by periodic boundary conditions (in case non-periodic case these atoms have no weight due to the FT window) 
!
      if(nsuper/=1.or.(nsuper==1.and.pos_units=='BOX')) then
        do k=1,3
          if(at_pos_in(k)<-n_row(k)/2.) at_pos_in(k) = at_pos_in(k)+n_row(k)   
          if(at_pos_in(k)>=n_row(k)/2.) at_pos_in(k) = at_pos_in(k)-n_row(k)
          if(j_shell_out==1) then
              if(at_pos_in2(k)<-n_row(k)/2.) at_pos_in2(k) = at_pos_in2(k)+n_row(k)   
              if(at_pos_in2(k)>=n_row(k)/2.) at_pos_in2(k) = at_pos_in2(k)-n_row(k)
          endif
        enddo
      endif

! *** prepare the output data 
!	
      if(nsuper==1.and.pos_units=='BOX') then
        at_pos_in = at_pos_in*a_par      !bring at_pos_in from the (-.5,.5) range back to the  box range 
        if(j_shell_out==1) at_pos_in2 = at_pos_in2*a_par
      endif														 

      at_pos_c(1:3,jrec) = at_pos_in
      at_pos_c(4,jrec) = at_charge_in															 

      if(n_traj>=1) then
        at_veloc_c(1:3,jrec) = at_veloc_in
        at_veloc_c(4,jrec) = at_mass_in	
      endif						
      if(n_traj==2) at_force_c(1:3,jrec) = at_force_in

      if(j_shell_out==1) then				!only if the shell data are going to be recorded
        at_pos_s(1:3,jrec) = at_pos_in2
        at_pos_s(4,jrec) = at_charge_in2							
        if(n_traj>=1) then
          at_veloc_s(1:3,jrec) = at_veloc_in2
          at_veloc_s(4,jrec) = at_mass_in2							
        endif
        if(n_traj==2) at_force_s(1:3,jrec) = at_force_in2
      endif

! *** accumulate the occupation number and the kinetic energy to refine the real temperature
      if(ifile==nfile_min) nsuper_r(jat) = nsuper_r(jat)+1
      if(n_traj>=1) then
        do k = 1,3
          e_kin(jat,k) = e_kin(jat,k) + at_mass_in*at_veloc_in(k)**2			!at_mass_c(i)=at_veloc_c(4,i)
          if(j_shell==1) e_kin_s(jat,k) = e_kin_s(jat,k) + at_mass_in2*at_veloc_in2(k)**2
        enddo
      endif !input method
    enddo read_loop

! *** indexing for the output
    allocate(at_ind_out(n_tot),at_name_out(n_atom),ind_at(n_atom))
    if(input_method=='BULK') then
      ind_at(1) = 0
      do jat = 2,n_atom
        ind_at(jat) = ind_at(jat-1)+nsuper_r(jat-1)
      enddo
    
      do i=1,n_tot
        jat = at_ind(4,i)
        ind_at(jat) = ind_at(jat)+1
        at_ind_out(ind_at(jat)) = at_ind(1,i)
      enddo
      jr => at_ind_out(1:n_tot)
    else
      jr => i_series
    endif			

    call cpu_time(t2)
    if(i_traj==nt_min.and.ifile==nfile_min) then												!analyze in detail the 1st snapshot
      print *,space, '1st snapshot: total of',jrec,' atoms read in',t2-t1,' sec'
    endif
                     
! *** normalize the kinetic energy and get the true temperature
  if(n_traj>=1) then
    do j=1,n_atom*n_eq
      e_kin(j,:) = e_kin(j,:)/nsuper_r(j)
      if(j_shell.eq.1) e_kin_s(j,:) = e_kin_s(j,:)/nsuper_r(j)
    enddo	 

    if(j_verb==1.and.ifile==nfile_max) print *
    if(j_verb==1.and.ifile==nfile_max) print *,space, 'Last frame no.',ifile

    if(j_verb==1.and.(ifile==nfile_min.or.ifile==nfile_max)) then
      print *
      print *,space, 'Cores E_kin(jat,:)',(e_kin(jat,:),jat=1,n_atom)
      if(j_shell.eq.1) print *,space, 'Shells E_kin(jat,:)',(e_kin_s(jat,:),jat=1,n_atom)
    endif

    temp_r_c = sum(e_kin(1:n_atom,:))/(n_atom*3*k_B) !the true core temperature

    if(j_shell.eq.1) then
      temp_r_s = sum(e_kin_s(1:n_atom,:))/(n_atom*3*k_B) !the true shell temperature
      if (abs(temp_r_c-temp_r_s).le..1*temp_r_c) then
        temp = (temp_r_c+temp_r_s)*.5						!hi-T limit: independent C and S vibrations
        if((ifile==nfile_min.or.ifile==nfile_max)) print *,space, 'Hi-T limit: independent C and S vibrations'
      else
        temp = temp_r_c+temp_r_s	!low-T limit: strongly bound C and S vibrations
        if((ifile==nfile_min.or.ifile==nfile_max)) print *,space, 'Low-T limit: strongly bound C and S vibrations'
      endif        
      if(ifile==nfile_min.or.ifile==nfile_max) print *,space, 'Real temperature: core/shell/total ',temp_r_c,temp_r_s,temp
      if(ifile==nfile_min.or.ifile==nfile_max) write(9,*) 'Real temperature: core/shell/total ',temp_r_c,temp_r_s,temp
    else
      temp = temp_r_c
      temp_r_s = .0
      if(ifile==nfile_min.or.ifile==nfile_max) print *,space, 'Real temperature: cores only ',temp
      if(ifile==nfile_min.or.ifile==nfile_max) write(9,*) 'Real temperature: cores only ',temp
    endif
  else
    temp = temp_par
    if(n_traj==1) then
      if(ifile==nfile_min.or.ifile==nfile_max) print *,space, 'Using nominal temperature [K] ',temp
      if(ifile==nfile_min.or.ifile==nfile_max) write(9,*) 'Using nominal temperature [K] ',temp
    endif
  endif
  
! *** First snapshot in the series only  
! *** get the atom position occupation numbers and compare them with those from the .par file
  at_name_out = at_name_par(1:n_atom)
  if(input_method=='CELL') then
    at_occup_r = (1.*nsuper_r)/nsuper 
!!        if(i_traj==nt_min.and.ifile==nfile_min) then
    if(ifile==nfile_min) then
      print *,space, 'Occupancies: nominal 		real'
      do ii=1,n_atom
        print *,space, '     ',at_name_out(ii),at_occup(ii),at_occup_r(ii)
      enddo
    endif
  else 
    at_occup_r = nsuper_r/real(n_tot) 
!!        if(i_traj==nt_min.and.ifile==nfile_min) then
    if(ifile==nfile_min) then
      print *,space, 'Bulk concentrations:'
      do ii=1,n_atom
        print *,space, '     ',at_name_out(ii),at_occup_r(ii)
      enddo
    endif
  endif 		!i_traj==nt_min.and.ifile==nfile_min

  t_dump = i_save*t_step


! *** define the record structure
  n_rec = (n_tot/l_rec4)														!for each position there are 4 components
  if(mod(n_tot,l_rec4)/=0) n_rec = n_rec+1

! *** generate output filename
  if(i_save<=9999) then
    write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_par),i_save
  elseif(i_save>=10000) then
    write(string,'(i8)') i_save
    file_dat = './data/'//trim(file_par)//'_n'//trim(adjustl(string))//'.dat'
  endif

  if(i_save==n_save_min.or.i_save==10*(i_save/10)) print *,space,trim(file_dat)
  i_save = i_save+1       

  call cpu_time(t1)
  open(2,file=file_dat,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4

! *** write the header record
  n_head = 4                        !number of header lines for including full cell info for non-orthogonal lattices
  write(string,*) n_head
  if(j_verb==1.and.ifile==nfile_min) print *,space, 'n_head',n_head

  header_record = dat_source//version//string//'   '//trim(time_stamp)
  i_rec = 1
  write(2,rec=i_rec) header_record

  write(header_record,nml=data_header_1)	
  i_rec = i_rec+1
  write(2,rec=i_rec) header_record

  write(header_record,nml=data_header_2)	
  i_rec = i_rec+1
  write(2,rec=i_rec) header_record

  write(header_record,nml=data_header_3)	
  i_rec = i_rec+1
  write(2,rec=i_rec) header_record

! *** do the rest	
  do i=1,n_rec-1
    i_rec = i_rec+1
    write(2,rec=i_rec) (at_ind(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
  enddo
  i = n_rec
  i_rec = i_rec+1
  write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
  write(2,rec=i_rec) (at_ind(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot) 
    
  do i=1,n_rec-1
    i_rec = i_rec+1
    write(2,rec=i_rec) (at_pos_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
  enddo
  i = n_rec
  i_rec = i_rec+1
  write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
  write(2,rec=i_rec) (at_pos_c(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)

  if(n_traj>=1) then
    do i=1,n_rec-1
      i_rec = i_rec+1
      write(2,rec=i_rec) (at_veloc_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
    enddo
    i = n_rec
    i_rec = i_rec+1
    write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
    write(2,rec=i_rec) (at_veloc_c(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)
  endif

  if(n_traj==2) then
    do i=1,n_rec-1
      i_rec = i_rec+1
      write(2,rec=i_rec) (at_force_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
    enddo
    i = n_rec
    i_rec = i_rec+1
    write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
    write(2,rec=i_rec) (at_force_c(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)
  endif

  if(j_shell_out==1) then
    do i=1,n_rec-1
      i_rec = i_rec+1
      write(2,rec=i_rec) (at_pos_s(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
    enddo
    i = n_rec
    i_rec = i_rec+1
    write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
    write(2,rec=i_rec) (at_pos_s(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)

    if(n_traj>=1) then
      do i=1,n_rec-1
        i_rec = i_rec+1
        write(2,rec=i_rec) (at_veloc_s(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
      enddo
      i = n_rec
      i_rec = i_rec+1
      write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
      write(2,rec=i_rec) (at_veloc_s(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)
    endif

    if(n_traj==2) then
      do i=1,n_rec-1
        i_rec = i_rec+1
        write(2,rec=i_rec) (at_force_s(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
      enddo
      i = n_rec
      i_rec = i_rec+1
      write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
      write(2,rec=i_rec) (at_force_s(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)
    endif
  endif 

  close(2)

  ifile = ifile+1

  deallocate(at_pos_c,ind_at,at_ind_out,at_name_out)
  if(n_traj>=1) deallocate(at_veloc_c)
  if(n_traj==2) deallocate(at_force_c)
  if(j_shell_out==1) deallocate(at_pos_s)
  if(j_shell_out==1.and.n_traj>=1) deallocate(at_veloc_s)
  if(j_shell_out==1.and.n_traj==2) deallocate(at_force_s)


  enddo frame_loop 
  close(1)
  
  enddo trajectory_loop 

  deallocate(at_name,at_ind,e_kin,at_occup_r,nsuper_r,i_series)
  if(j_shell.eq.1) deallocate(e_kin_s)

  deallocate(at_name_par,at_label,at_base_in,at_base,ind_l,i_site,at_occup,at_mass_in_c,at_mass_in_s)

  CALL SYSTEM_CLOCK (COUNT = sc_c2)
  print *,space, 'Trajectory finished: ',ifile-nfile_min,' .dat files written in SYS time',(sc_c2-sc_c1)/sc_r,' sec'
  write(9,*) 'Trajectory finished: ',ifile-nfile_min,' .dat files written in SYS time',(sc_c2-sc_c1)/sc_r,' sec'

  stop
end program mp_lbin56

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

