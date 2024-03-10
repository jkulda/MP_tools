      
program mp_dbin56
 
! ************************************************************************************* 
! ***** 
! *****  %%%%%%%%%%%%%%%%   		  program MP_BIN 1.56   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! ***** 
! *****   converts MD trajectory data (DL_POLY or equivalent) to the MP_TOOLS binary form 
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
! *****   %%%%%%%%%%%%%%%%   			program MP_BIN 1.56  				 %%%%%%%%%%%%%%%%%%%%%%%% 
! ***** 
! *****		!!!! incompatible with the old MD 1.4 file format because of N_ROW(3) !!!! 
! ***** 
! ***** 
! ***** Ver. 1.1 - reads properly the header information from the ASCII input file 
! ***** Ver. 1.2 - looks for the header line 'timestep' instead of asking lines to skip 
! ***** Ver. 1.3 - revised: 
! *****		- only the "cryst_structure" type is kept, no substance-related data are left in the code 
! *****		- looks for auxiliary file with structure info <title.par> 
! *****		- the crystal structure is flexible (up to 20 atoms), the 3 oxygens in perovskite are still decodable 
! *****		- cubic supercell is assumed, nrow etc. is then calculated, only n_atom to be supplied 
! *****		- the header now contains n_atom, j_force and temp (on last position for compatibility) 
! *****		- when atom forces (j_force = 1) are present in the trajectory, they can be included at the end of the binary records 
! ***** Ver. 1.4 - dialogue revision 
! *****		- ./data subdirectory used for trajectory and primary .bin files 
! ***** Ver. 1.41 - cell parameter updated for each snapshot 
! *****           - "history" replaced by "trajectory" 
! ***** Ver. 1.42 - skips first nfile_min-1 snapshots and reads the sequence with a step nfile_step until nfile_max 
! ***** Ver. 1.43 - refines the real temperature from kinetic energies of cores 
! ***** Ver. 1.44 - identifies atoms according to the fractional positions in the PAR file  
! *****						- abandons the "cryst_structure" TYPE 
! *****						- distinguishes low-T and hi-T behaviour of cores and shells:  
! *****															low-T strongly bound, k_B*T/2 = E_kin = E_kin_C + E_kin_S 
! *****															hi-T independent movement, k_B*T/2 = E_kin = E_kin_C = E_kin_S 
! ***** Ver. 1.45 - reads-in the whole snapshot in one go  
! ***** 					- uses orthorhombic supercell with N_ROW(3)  
! ***** 					- detects the NROW, if cubic, and refines the lattice parameters  
! *****           - at_ind(i,1:3) contains cell indices, at_ind(i,4) atom identifier (1:n_atom) 
! *****           - writes a terminal record with fictive atom name 'End' 
! ***** Ver. 1.46 - corrected: missing/wrong velocities input, test of j_force on output, atom names 
! *****						- correct handling of partial and mixed occupancies  
! ***** 
! ***** Ver. 1.50 - takes over ver. 1.46 with minor bug fixes 
! *****           - allows for orthorombic (non-cubic) supercells by n_row(3) 
! ***** Ver. 1.51 - restructured data access: possibility to point towards ASCII data on an external volume 
! ***** 					- possibility to specify trajectory filename extension 
! ***** 					- all data arrays allocatable, no predefined array size limits 
! ***** 					- supercell format in .PAR 
! ***** Ver. 1.52 - cycle over a series of trajectory files to create a contiguous series of binary snapshot files 
! *****						- waits until a subsequent trajectory file appears in its directory 
! ***** Ver. 1.53 - new binary data format  
! *****						- record length 4096, all data 32bit length 
! *****						- at_mass is saved as at_veloc(4,j,k) 
! *****						- at_charge is saved as at_pos(4,j,k) 
! *****						- at_displ is not supported anymore 
! *****						- cycle over trajectory files renewed 
! *****           - CELL, FAST and BULK input methods 
! *****           - BULK input uses simulation box sizes 
! *****           - j_force replaced by n_traj 
! *****           - j_centred in .PAR indicates coordinate centring in box (1/0) 
! ***** 
! *****   FROZEN on 05/09/2022 10:45 and forked as mp_bin54 
! ***** 
! *****  Ver. 1.54 - uses NAMESLIST to read .PAR files and save headers 
! ***** 
! *****  Ver. 1.55 - bug fixes 
! ***** 
! *****  Ver. 1.56 - handles triclinic unit cell geometry
! *****            - for core-shell model records the centre-of-mass positions&velocities for J_SHREC=0
! *****
! ***** reads ASCII text output from MD simmulations (trajectory file) 
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

  logical :: found,found_txt,t_single
  character(4)   :: at_name_in,at_name_in2,at_name_in3,head,sim_type_lc,pg_out,version
  character(10)	 :: prompt,space = '          '
  character(10)  :: c_date,c_time,c_zone,ext,number
  character(16)  :: filter_name,pos_units
  character(40)  :: subst_name,string,section,mp_tool
  character(128) :: line,cwd_path,data_path,time_stamp,rec_str
  character(128) :: file_master,file_master_out,file_dat,file_trajectory,file_inp,file_log
  character(4*l_rec) :: header_record

  integer ::  i_dom,i_dom_rec,n_dom,n_tot_in,j_struct,n_tot,i_save
  integer ::  at_no,at_ind_base(3),at_ind_in(3)
  integer ::  j_yes,nskip,nfile_min,nfile_max,nfile_step,i_time(8),n_save_min,izero,indzero(3)
  integer ::  ios,ios_t,ierr,i,j,k,ii,i2,i3,jl,jat,j_step,j_label,n_label,j_first,j_read,j_verb,j_proc,j_shrec,j_test
  integer ::  inrec,jrec,nrec,i_rec,l_rec4,ifile,ncell,nsuper,nrow,nlayer,n_site,j_shell
  integer ::  sc_c1,sc_c2,sc_m,nt_min,nt_max,nt_step,i_traj,j_mult,j_basis,j_centred

  real :: at_mass_in,at_mass_in2,at_charge_in,at_charge_in2,at_displ_in,sc_r
  real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force_in(3),at_force_in2(3),at_base_shift(3)
  real ::	dummy,at_pos2(3),at_pos3(3),at_veloc2(3),a_cell(3,3),a_cell_1(3,3),a_cell_inv(3,3),a_cell_par(3),a_cell_half(3),at_pos_centre(3)
  real :: t1,t2,filter_fwhm,t_step,zero,at_zero(3),pos_inp(3)
  real :: temp_par,temp_r_c,temp_r_s,temp_r_cg,temp_cs,eps_x

  character(4),allocatable :: at_name(:),at_label(:)
  integer,allocatable ::  ind_l(:),i_site(:,:),ind_at(:)
  integer,allocatable,target :: i_series(:),at_ind_out(:)
  integer,pointer ::  jr(:)
  real,allocatable ::  sum_pos(:,:,:),ord(:),e_kin(:,:),e_kin_s(:,:),e_kin_cg(:,:),at_base_in(:,:),at_base(:,:),at_occup(:)

! **** the following variables MUST have the following 32bit sizes or multiples because of alignement in the binary output file 
!
  character(4),allocatable :: at_name_par(:),at_name_out(:)
  integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)

  real(4),allocatable ::	at_pos_c(:,:),at_veloc_c(:,:),at_force_c(:,:),at_occup_r(:)
  real(4),allocatable ::	at_pos_s(:,:),at_veloc_s(:,:),at_force_s(:,:)

  character(16)  :: sim_type,input_method,dat_type,dat_source,file_par
  integer(4)     :: n_row(3),n_atom,n_eq,j_shell_out,n_traj,n_cond,idum,n_rec,n_head,n_head_in1,n_head_in2
  real(4)        :: rec_zero(l_rec),t_ms,t_dump,a_par(3),angle(3),temp

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,temp_cs,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
  namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
  namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
  namelist /mp_gen/ j_verb,j_proc       
                      !general rule: namelists of tools should only contain their local parameters
                      !what is of global interest they should pass into data_header
  namelist /mp_bin/ subst_name,sim_type,dat_type,input_method,pos_units,data_path,ext,rec_str,&
 &					j_mult,n_head_in1,n_head_in2,n_tot_in,n_atom,n_row,j_basis,j_centred,j_test,j_shrec,&
 &          a_cell_par,at_base_shift,eps_x,temp_par,t_step

!     namelist /mp_sqom/ n_int,s_trig,j_oneph,j_qsq
!     namelist /mp_pdf/ n_pdf,pdf_step,j_gauss,n_h,j_weight,j_smooth,n_corr

  data rec_zero/l_rec*.0/

! 
! ********************* Initialization *******************************      
  version = '1.56'
  prompt = 'MP_DBIN>  '
  mp_tool = 'MP_DBIN '//version

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
  t_single = .true.
  nt_min = 1
  nt_max = 1
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
  j_shrec = 0
  t_ms = .0
  temp_par = .0
  temp_cs = .0
  at_base_shift = .0
  filter_name = 'nn'
  filter_fwhm = .0
  pos_units = 'ANGSTROM'
   
! *** read auxiliary file <file_par.par> with structure parameters, atom names and further info
  print *,prompt, 'Parameter file name (.par will be added)'
  do
    read(*,*) string
    if(len(trim(string))<=16) exit
    print *,space,'Name exceeds 16 character, modify & retype ...'
  enddo
  file_par = trim(string)
  file_inp = trim(file_par)//'.par'

  open(4,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,space, 'File ',trim(file_inp),' not found! Stop execution.'
    stop
  endif

  write(9,*) 'Read parameter file:  ',trim(file_inp)

  read(4,nml=mp_gen)

  j_test = 1
  data_path = './data/'
  ext = ''
  rewind(4)
  read(4,nml=mp_bin)
  call up_case(sim_type)
  call up_case(dat_type)
  call up_case(input_method)
  call up_case(pos_units)
  
  if(dat_type/='DL_POLY') then
    print *,space, 'Data type in .PAR ', trim(pos_units),' not corresponding to DL_POLY, continue (1/0)?'
    read(*,*) j_yes
    if(j_yes==0) stop    
  endif

  if(pos_units/='ANGSTROM') then
    print *,space, 'Position units ', trim(pos_units),' not implemented with DL_POLY data type, ANGSTROM needed!'
    stop
  endif
  
  if(j_centred/=1) then
    print *,space, 'Atom coordinates not centred (J_CENTRED/=1, unusual with DL_POLY data type), continue (1/0)?'
    read(*,*) j_yes
    if(j_yes==0) stop    
  endif
  
  if(ext=='ext'.or.ext=='EXT') ext=''
  if(ext/=''.and.index(ext,'.')==0) ext='.'//ext

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
  print *,space, 'Sim_type, dat_type, input method, position units: ',sim_type,dat_type,input_method,pos_units		
  
  allocate(at_name_par(n_atom),at_label(n_atom),at_base_in(n_atom,3),at_base(n_atom,3))   !at_base would include at_base_shift & saved in data file header
  allocate(ind_l(n_atom),i_site(n_atom,n_atom),at_occup(n_atom))

  nsuper = n_row(1)*n_row(2)*n_row(3)
  nlayer = n_row(1)*n_row(2)						!to be used for record number calculation for CELL data
   
! *** Read the atom positions       
  section = 'atoms'
  rewind(4)
  do
    read(4,'(a)',iostat=ios) string
    if(ios/=0) then
      print *,space, 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
      stop
    endif
    if(string(1:6).eq.section) exit	!find the mp_simple part of the .par file
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
  
  n_label = j_label
  print *,space, 'Atom labels ',n_label
  do j=1,n_label
    print *,space, ind_l(j),at_label(j),(i_site(j,i),i=1,ind_l(j))
  enddo	

  print *,space, trim(subst_name),' structure info (atoms): '	  
  do j=1,n_atom
    if(input_method=='CELL'.or.input_method=='FAST') then
      print *,space, j,at_name_par(j),at_base_in(j,:)
    else
      print *,space, j,at_name_par(j)
    endif
  enddo 
			
  angle = 90.				!assuming orthogonal lattice

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
  i_traj = 0

  print *,prompt, 'Master for MD trajectory filename: '
  read(*,*) file_master

  print *,prompt, 'Read input files number (0 0 no numbers): min,max'
  read(*,*) nt_min,nt_max
  nt_step = 1
  t_single = (nt_min==0.and.nt_max==0)

!! *** tread 1st frame of 1st history file to get info
  i_traj = nt_min

  if(t_single) then
    file_trajectory = file_master
  else
    if(i_traj>=1.and.i_traj<=9)    write(number,'(i1.1)') i_traj
    if(i_traj>=10.and.i_traj<=99)  write(number,'(i2.2)') i_traj
    if(i_traj>=100.and.i_traj<=999)write(number,'(i3.3)') i_traj
    if(i_traj>=1000.and.i_traj<=9999)write(number,'(i4.4)') i_traj
    if(i_traj>=1000)then 
      write(number,'(i8)') i_traj
      number = trim(adjustl(number))
    endif					
    file_trajectory = trim(file_master)//trim(number)
  endif

  file_inp = trim(data_path)//trim(adjustl(file_trajectory))//trim(ext)

  print *,space,'Input trajectory file:  ',trim(file_inp)

  open (1,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,space, "Can't open the file ",trim(file_inp),'!'
    stop
  endif

  write(9,*) 'Reading MD trajectory file:  ',trim(file_inp)
             
! ***  skip the first <nskip> lines until 'timestep' !normally they are two
   nskip = 20 ! number of records to test
  call down_case_2(sim_type(1:4),sim_type_lc) 
  do i=1,nskip
    read(1,*) line
    if(index(line,sim_type(1:4)).ne.0) then
      head = sim_type(1:4)
      exit
    elseif(index(line,sim_type_lc).ne.0) then
      head = sim_type_lc
      exit
    endif        
    if(i.eq.nskip) then
      print *,space, 'MP_BIN: trajectory header not found'
      stop
    endif
  enddo 
  
  backspace(1)
     
! *** read the header of the first snapshot  
! 
! *** typical header of a snapshot (4 lines): 
! timestep      750125      552960           1           3  3.99999990E-04   270.04999     
!      194.2110000000        0.0000000000        0.0000000000           
!        0.0000000000      194.2110000000        0.0000000000           
!        0.0000000000        0.0000000000      194.2110000000          

  read(1,*) string,j_step,n_tot_in,n_traj,n_cond,t_ms,t_dump		!t_ms is MD microstep, t_dump is the snapshot time; n_cond=0 non-periodic bound_cond (all others are periodic)
                                                                  !n_tot_in total number of atoms in cells n_atom*nrow**3
                                                                     
  if(n_cond==0) then
    print *,prompt, 'Non-periodic boundary conditions! Continue? (1/0)'
    read(*,*) j
    if(j==0) stop
    print *,space, 'This may have adverse effects on further treatment ...'		  
  endif
                                                                  
  do j=1,3
    read(1,*) a_cell(j,:)
  enddo
  
  do j=1,3
    a_cell(j,:) = a_cell(j,:)/n_row(j)        !a_cell becomes unit cell size
    a_par(j) = norm2(a_cell(j,:))
  enddo			

  angle(1) = dot_product(a_cell(2,:),a_cell(3,:))/(a_par(2)*a_par(3))
  angle(2) = dot_product(a_cell(1,:),a_cell(3,:))/(a_par(1)*a_par(3))
  angle(3) = dot_product(a_cell(1,:),a_cell(2,:))/(a_par(1)*a_par(2))
  angle = acos(angle)

!   if(j_verb==1) then
!     print *,'angle_rad',angle
!     print *,'angle_deg',angle*180./pi
!   endif

  a_cell_1 = a_cell
  call gjinv(a_cell_1,3,3,a_cell_inv,3,ierr)
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

! *** read 1st atom record         
  read(1,*) at_name_in,at_no,at_mass_in,at_charge_in,at_displ_in
    read(1,*) at_pos_in
    if(n_traj>=1) read(1,*) at_veloc_in
    if(n_traj==2) read(1,*) at_force_in

! *** read 2nd atom record         
  read(1,*) at_name_in2,at_no,at_mass_in,at_charge_in,at_displ_in
    read(1,*) at_pos2
    if(n_traj>=1) read(1,*) at_veloc2
    if(n_traj==2) read(1,*) at_force_in

! *** read 3rd atom record         
  read(1,*) at_name_in3,at_no,at_mass_in,at_charge_in,at_displ_in
    read(1,*) at_pos3
    if(n_traj>=1) read(1,*) at_veloc_in
    if(n_traj==2) read(1,*) at_force_in
  close(1)

! *** analyse the input
  if(at_name_in.eq.at_name_in2) then
    j_shell = 0				
    print *,space, 'No shells'               !', supercell is:',nrow,'^3',' j_shell =,',j_shell
  else if(trim(at_name_in)//'s'.eq.at_name_in2.and.at_name_in.eq.at_name_in3) then
    j_shell = 1
    n_tot_in = n_tot_in/2
    print *,space, 'Found a shell candidate   ',at_name_in2		!,' j_shell =',j_shell    !', supercell is:',nrow,'^3'
  else
    print *,space, 'Strange primary data, please check them!'
    print *,space, at_name_in,trim(at_name_in)//'s'
    print *,space, at_name_in2
    print *,space, at_name_in3
    stop
  endif
  
  if(n_traj==0) then
    print *,prompt, 'Type in nominal temperature [K]:'
    read(*,*) temp_par
  endif
  
  if(n_traj>=1.and.at_veloc_in(1)==at_veloc2(1).and.at_veloc_in(2)==at_veloc2(2).and.at_veloc_in(3)==at_veloc2(3)) then
    print *,space, 'Strange velocities:'
    print *,space, at_veloc_in   
    print *,space, at_veloc2          
    print *,prompt, 'Type in nominal temperature [K]:'
    read(*,*) temp_par
  endif
  
  if(j_shell==1) then
    print *,space,'Core & shell data found'
    if(j_shrec==0) print *,space,'Core-shell centre-of-mass data will be recorded (default, change in .PAR with caution)'
    if(j_shrec==1) print *,space,'Core-shell data WILL be recorded individually'
    if(j_shrec==1) write(9,*)'Core-shell data WILL be recorded individually'
  endif
  j_shell_out = j_shell*j_shrec

  print *,space, 'Simulation type = ',sim_type
  print *,space, 'Trajectory recording mode =', n_traj
  print *,space, 'Boundary conditions ', n_cond
  print *,space, 'Trajectory time start [ps]:',t_dump
  if(n_traj==0) write(9,*) 'Using nominal temperature [K] ',temp_par

  write(9,*) 'Simulation type = ',sim_type
  write(9,*) 'Trajectory recording mode =', n_traj
  write(9,*) 'Boundary conditions ', n_cond
  write(9,*) 'Trajectory time start, step [ps]:',t_dump,t_ms
  if(n_traj==0) write(9,*) 'Using nominal temperature [K] ',temp_par

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

  allocate (at_ind(4,n_tot),SOURCE=0)
  allocate(at_name(n_tot),SOURCE='    ')
  allocate(e_kin(n_atom,3),at_occup_r(n_atom))
  allocate(nsuper_r(n_atom),SOURCE=0) !atom number jat is the first of the four indices
  if(j_shell.eq.1)allocate(e_kin_s(n_atom,3),SOURCE=0.0)
  if(j_shell.eq.1)allocate(e_kin_cg(n_atom,3),SOURCE=0.0)
  allocate(i_series(n_tot))
  i_series = (/ (i, i = 1, n_tot) /)

! *** Now ready to cycle over trajectory files, each snapshot to be saved in a separate binary file

  CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)

  nsuper_r = 0
  ifile = 1

  trajectory_loop: do i_traj=nt_min,nt_max
    print *
    if(t_single) then
      file_trajectory = file_master
    else
      if(i_traj>=1.and.i_traj<=9)    write(number,'(i1.1)') i_traj
      if(i_traj>=10.and.i_traj<=99)  write(number,'(i2.2)') i_traj
      if(i_traj>=100.and.i_traj<=999)write(number,'(i3.3)') i_traj
      if(i_traj>=1000.and.i_traj<=9999)write(number,'(i4.4)') i_traj
      if(i_traj>=1000)then 
        write(number,'(i8)') i_traj
        number = trim(adjustl(number))
      endif					
      file_trajectory = trim(file_master)//trim(number)
    endif

    file_inp = trim(data_path)//trim(adjustl(file_trajectory))//trim(ext)

    if(.not.t_single) print *,space, 'Input trajectory file:  ',trim(file_inp)

    t1 = .0
    do 
      inquire (file=file_inp,exist=found)
      if(.not.found) then
        if(t1==0.) print *,space, 'File ',trim(file_inp),' not found! Waiting for it ...'			!if(t1==0.) 
        CALL SLEEP(10)
        t1 = t1+10.
      else
        exit
      endif
    enddo

  open (1,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,space, "Can't open the file ",trim(file_inp),'!'
    stop
  endif
    
               
! *** cycle over snapshots, each to be saved in a separate binary file
  j_read = 0

  frame_loop: do 				!we have to go through all the snapshots in the .TXT input

    do 						
      read(1,*,iostat=ios) line
      if(ios<0) then   !end of file 
!!						print *,'End of trajectory file: ',ios
        exit frame_loop 
      endif 
!!          if(ios>0) cycle     !corrupt input data, go on
      if(ios>0) then
        print *,space, 'Line after:',line
        print *,space, 'Input error, corrupt trajectory data?'
        stop
      endif
      if(index(line,head).ne.0) exit
    enddo 
    

! *** skip the snapshots that are not requested to be read
    if((ifile.lt.nfile_min).or.(mod(ifile-nfile_min,nfile_step).ne.0)) then  !skip the snapshot
      ifile = ifile+1
      cycle frame_loop         ! cycle the frame_loop
    else
      backspace(1)	
      j_read = j_read+1
    endif                    

! *** no need to re-read the first header line but 
! *** take fresh lattice parameters for each snapshot - they may evolve
    read(1,*) string,j_step,n_tot_in,n_traj,n_cond,t_ms,t_dump

    do j=1,3
      read(1,*) a_cell(j,:)
    enddo 

    do j=1,3
      a_cell_half(j) = .5*sum(a_cell(:,j))
      if(j_centred==0) then
        at_pos_centre(j) = a_cell_half(j)   
      else
        at_pos_centre(j) = .0
      endif
      a_cell(j,:) = a_cell(j,:)/n_row(j)        !for n_row>1 a_cell becomes unit cell
      a_par(j) = norm2(a_cell(j,:))
      if(nsuper==1) a_cell(j,:) = a_cell(j,:)/a_par(j)         !effective a_par=1 and Q will be in [A-1]
    enddo 
        
    if(i_traj==nt_min.and.ifile==nfile_min) then
      if(nsuper/=1) then
        print *,space, 'Lattice parameter estimate =',a_par 
      else
        print *,space, 'Simulation box size =',a_par
      endif			
    endif
    
    angle(1) = dot_product(a_cell(2,:),a_cell(3,:))/(a_par(2)*a_par(3))
    angle(2) = dot_product(a_cell(1,:),a_cell(3,:))/(a_par(1)*a_par(3))
    angle(3) = dot_product(a_cell(1,:),a_cell(2,:))/(a_par(1)*a_par(2))
    angle = acos(angle)

    a_cell_1 = a_cell
    call gjinv(a_cell_1,3,3,a_cell_inv,3,ierr)       !!! a_cell would get destroyed here !!!
    if(ierr==1) then
      print *,space, 'Singular cell vector matrix, check your HISTORY file!'
      stop
    endif 

! *** get the actual time step of the sequence
    if(ifile==nfile_min) then
      t_step = t_dump
    elseif(ifile==nfile_min+1) then
      t_step = t_dump-t_step
      t_dump = t_dump-t_step !get the old t_dump for the old frame           
! *** correct t_step in the first snapshot				  
      open(2,file=file_dat,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4
      write(header_record,nml=data_header_1)	
      write(2,rec=2) header_record
      close(2)
      t_dump = t_dump+t_step !get the right t_dump for the present frame          
    endif
     
! ***  read-in the text of a snapshot in one go
    if(ifile==nfile_min) print *,space, 'reading the 1st snapshot (takes a few seconds) ...'

    allocate(at_pos_c(4,n_tot),SOURCE=0.0)
    if(n_traj>=1) allocate(at_veloc_c(4,n_tot),SOURCE=0.0)
    if(n_traj==2) allocate (at_force_c(4,n_tot),SOURCE=0.0)
    if(j_shell_out.eq.1) then
      allocate(at_pos_s(4,n_tot),SOURCE=0.0)
      if(n_traj>=1) allocate(at_veloc_s(4,n_tot),SOURCE=0.0)
      if(n_traj==2) allocate (at_force_s(4,n_tot),SOURCE=0.0)
    endif
  
    e_kin = .0
    if(j_shell.eq.1) e_kin_s = .0
    if(j_shell.eq.1) e_kin_cg = .0

    call cpu_time(t1)				

    read_loop: do  inrec=1,n_tot_in        !swallow the snapshot

! *** first read the CORE data  

      read(1,*,iostat=ios_t) at_name_in,at_no,at_mass_in,at_charge_in,at_displ_in

      if(ios_t<0.or.at_name_in.eq.head) then 
        exit read_loop    !end of snapshot, end of file
      elseif(ios_t>0) then
        print *,space, 'Error on ASCII input IOS: ',ios_t
      endif

      read(1,*) at_pos_in
      if(j_centred==0) at_pos_in = at_pos_in-at_pos_centre


      if(n_traj>=1) read(1,*,iostat=ios_t) at_veloc_in
      if(ios_t/=0) then
        print *,space, 'Input problem: ios_t,inrec,ifile,at_name_in,head',ios_t,inrec,ifile,at_name_in,head
      endif
      if(n_traj==2) read(1,*) at_force_in

! 
! *** now read the SHELL data if needed
      if(j_shell.eq.1) then
        read(1,*) at_name_in2,at_no,at_mass_in2,at_charge_in2,at_displ_in
        if(trim(at_name_in)//'s'.ne.at_name_in2) then
          print *,space,'wrong core-shell sequence',inrec,at_name(inrec),at_name_in
          stop
        endif
        read(1,*) at_pos_in2
        if(j_centred==0) at_pos_in2 = at_pos_in2-at_pos_centre
        if(n_traj>=1) read(1,*) at_veloc_in2
        if(n_traj==2) read(1,*) at_force_in2
        at_no = at_no/2
      endif		!j_shell

! ***  treat the CORE data and get the right labels & positions 
      jl = 0				
      if(j_basis==0) then
        do j=1,n_atom
          if (at_name_in.eq.at_label(j)) then
            jl=j			!atom label found
            exit
          endif
        enddo
        if(jl.ne.j) then			!after "normal" loop exit
          print *,space, 'atom ',at_name_in,' not found in .PAR'
          stop
        endif
      else
         do j=1,n_atom
          if (at_name_in.eq.at_name_par(j)) then
            jl=j			!atom label found
            exit
          endif
        enddo
        if(jl.ne.j) then			!after "normal" loop exit
          print *,space, 'atom ',at_name_in,' not found in .PAR'
          stop
        endif
      endif

! *** transform the atom positions to pseudo-cubic r.l.u. coordinates 
! *** DL_POLY data are always CENTRED and in ANGSTROM units  
!
        at_pos_in = matmul(a_cell_inv,at_pos_in)    			! at_pos_in are in Eucleidian direct space coordinates, we need to get them into the lattice coordinates to be able to use periodic boundary contions for FT
                                                          ! for nsuper=1 nothing happens as a_cell and a_cell_inv are unit matrices
        if(j_shell==1) at_pos_in2 = matmul(a_cell_inv,at_pos_in2)     
        
! *** index BULK data 
      if(input_method=='BULK') then         !jl identifes label (chem species), jat basis position in CELL
        jat = jl	
        at_ind(1,inrec) = inrec
        at_ind(2:3,inrec) = 0
        at_ind(4,inrec) = jat
        jrec = inrec

! *** index CELL data 
      elseif(input_method=='CELL'.or.input_method=='FAST')then		
        if(j_basis==1) then
          jat = jl        !atom site was identified by name
        else
          if(input_method=='FAST') then
            if(at_no<=2*nsuper) jat = (at_no-1)/nsuper+1
            if(at_no>2*nsuper) jat = mod((at_no-2*nsuper-1),3)+3
          else													!input_method=='CELL'
            if(ind_l(jl).eq.1) then
              jat = i_site(jl,1)	!identify atom site by fractional position 
            else
              do ii=1,ind_l(jl)
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
        endif     !j_basis

!! *** calculate cell indices ix,iy,iz from atomic positions shifted to cell origin 

        at_ind_in = anint(at_pos_in-at_base_in(jat,:))+at_ind_base  !they will serve as pointers to the right order of atom records
        at_pos_in = at_pos_in+at_base_shift			!now the supercell will be centred & basis origin in at_base_shift
        if(j_shell.eq.1) at_pos_in2 = at_pos_in2+at_base_shift

        do k=1,3
          if(at_ind_in(k)==0) then
            at_ind_in(k) = at_ind_in(k)+n_row(k)
          elseif(at_ind_in(k)==n_row(k)+1) then
            at_ind_in(k) = at_ind_in(k)-n_row(k)
          endif
        enddo
   
        jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1)
      
        if(jrec<1.or.jrec>n_tot.or.jat<1.or.jat>n_atom) then
          print *,space, 'JREC wrong:',jrec,' > ',n_tot,' check n_tot_in, n_row and j_centred in the .PAR'
          print *,space, 'jrec,at_ind_in',jrec,at_ind_in
          print *,space, 'at_pos_in',at_pos_in,at_pos_in2
          print *,space, 'n_tot,nsuper,nlayer,n_row',n_tot,nsuper,nlayer,n_row
          stop
        endif

        at_ind(1:3,jrec) = at_ind_in
        at_ind(4,jrec) = jat
      endif		!input_method CELL
     
! *** enforce atom positions within the box limits by periodic boundary conditions (in case non-periodic case these atoms have no weight due to the FT window; for nsuper=1 no need - it's guaranteed) 
!
      if(nsuper/=1) then
        do k=1,3
          if(at_pos_in(k)<n_row(k)/2.) at_pos_in(k) = at_pos_in(k)+n_row(k)   
          if(at_pos_in(k)>=n_row(k)/2.) at_pos_in(k) = at_pos_in(k)-n_row(k)
          if(j_shell==1) then
              if(at_pos_in2(k)<n_row(k)/2.) at_pos_in2(k) = at_pos_in2(k)+n_row(k)   
              if(at_pos_in2(k)>=n_row(k)/2.) at_pos_in2(k) = at_pos_in2(k)-n_row(k)
          endif
        enddo
      endif

! *** prepare the output data     ! at_pos_in for BULK with n_row=1 will be recorded in ANGSTROM units 
!

      if(j_shell==0) then				    !no shells at all
        at_pos_c(1:3,jrec) = at_pos_in
        at_pos_c(4,jrec) = at_charge_in
        if(n_traj>=1) then
          at_veloc_c(1:3,jrec) = at_veloc_in
          at_veloc_c(4,jrec) = at_mass_in	
        endif						
        if(n_traj==2) at_force_c(1:3,jrec) = at_force_in
      elseif(j_shell==1.and.j_shell_out==0) then				!calculate core-shell centre of mass									 
        at_pos_c(1:3,jrec) = (at_mass_in*at_pos_in+at_mass_in2*at_pos_in2)/(at_mass_in+at_mass_in2)
        at_pos_c(4,jrec) = at_charge_in+at_charge_in2
        if(n_traj>=1) then
          at_veloc_c(1:3,jrec) = (at_mass_in*at_veloc_in+at_mass_in2*at_veloc_in2)/(at_mass_in+at_mass_in2)
          at_veloc_c(4,jrec) = at_mass_in+at_mass_in2	
        endif						
        if(n_traj==2) at_force_c(1:3,jrec) = at_force_in+at_force_in2
      elseif(j_shell_out==1) then				!only if the shell data are going to be recorded
        at_pos_c(1:3,jrec) = at_pos_in
        at_pos_c(4,jrec) = at_charge_in
        at_pos_s(1:3,jrec) = at_pos_in2
        at_pos_s(4,jrec) = at_charge_in2							
        if(n_traj>=1) then
          at_veloc_c(1:3,jrec) = at_veloc_in
          at_veloc_c(4,jrec) = at_mass_in	
          at_veloc_s(1:3,jrec) = at_veloc_in2
          at_veloc_s(4,jrec) = at_mass_in2							
        endif
        if(n_traj==2) at_force_s(1:3,jrec) = at_force_in2
      else
        print *,space,'ERR: j_shell_out out of range 0..1 '
        stop
      endif

! *** accumulate the occupation number and the kinetic energy to refine the real temperature
      if(ifile == nfile_min) nsuper_r(jat) = nsuper_r(jat)+1
      if(n_traj>=1) then
        do k = 1,3
          e_kin(jat,k) = e_kin(jat,k) + at_mass_in*at_veloc_in(k)**2			!at_mass_c(i)=at_veloc_c(4,i)
          if(j_shell==1) then
            e_kin_s(jat,k) = e_kin_s(jat,k) + at_mass_in2*at_veloc_in2(k)**2        !we drop the 1/2 factor everywhere as we will do later with k_B
            e_kin_cg(jat,k) = e_kin_cg(jat,k) + (at_mass_in*at_veloc_in(k)+at_mass_in2*at_veloc_in2(k))**2/(at_mass_in+at_mass_in2)
          endif
        enddo
      endif
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
      do j=1,n_atom
        e_kin(j,:) = e_kin(j,:)/(nsuper_r(j))
        if(j_shell.eq.1) e_kin_s(j,:) = e_kin_s(j,:)/(nsuper_r(j))
        if(j_shell.eq.1) e_kin_cg(j,:) = e_kin_cg(j,:)/nsuper_r(j)
      enddo

      if(j_verb==1.and.i_traj==nt_min.and.ifile==nfile_min) then
        print *
        print *,space, 'Cores E_kin(jat,:)',(.5*e_kin(jat,:),jat=1,n_atom)
        if(j_shell.eq.1) print *,space, 'Shells E_kin(jat,:)',(.5*e_kin_s(jat,:),jat=1,n_atom)
      endif

! *** get the true temperature
    temp_r_c = sum(e_kin(1:n_atom,:))/(n_atom*3*k_B) !the true core temperature          !we drop the 1/2 factor with k_B as we did before with m*v**2
    temp = temp_r_c
    
    if(j_shell.eq.1) then
      temp_r_s = sum(e_kin_s(1:n_atom,:))/(n_atom*3*k_B) !the true shell temperature
      temp_r_cg = sum(e_kin_cg(1:n_atom,:))/(n_atom*3*k_B) !the true LAMMPS temperature
      
      temp = temp_r_cg                      !temperature of the core-shell centre-of-mass
      temp_cs = temp_r_c+temp_r_s-temp_r_cg !internal temperature of the core-shell pair
      
      if (abs(temp_r_c-temp_r_s).le..1*temp_r_c) then
        if((ifile==nfile_min.or.ifile==nfile_max)) print *,space, 'Hi-T limit: independent C and S vibrations'
      else
        if((ifile==nfile_min.or.ifile==nfile_max)) print *,space, 'Low-T limit: strongly bound C and S vibrations'
      endif        
      if(ifile==nfile_min.or.ifile==nfile_max) print *,space, 'Real temperature: core, shell      ',temp_r_c,temp_r_s
      if(ifile==nfile_min.or.ifile==nfile_max) write(9,*) 'Real temperature: core, shell      ',temp_r_c,temp_r_s
      if(ifile==nfile_min.or.ifile==nfile_max) print *,space, 'Real temperature: CG, CS           ',temp_r_cg,temp_cs
      if(ifile==nfile_min.or.ifile==nfile_max) write(9,*) 'Real temperature: CG, CS           ',temp_r_cg,temp_cs
    else
      temp_r_s = .0
      temp = temp_r_c
      if(ifile==nfile_min.or.ifile==nfile_max) print *,space, 'Real temperature: atoms/cores_only ',temp
      if(ifile==nfile_min.or.ifile==nfile_max) write(9,*) 'Real temperature: atoms/cores_only ',temp
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
                  !analyze in detail the 1st snapshot
    at_name_out = at_name_par(1:n_atom)
    if(input_method=='CELL'.or.input_method=='FAST') then
      at_occup_r = (1.*nsuper_r)/nsuper
      if(i_traj==nt_min.and.ifile==nfile_min) then
        print *,space, 'Occupancies: nominal 		real'
        do ii=1,n_atom
          print *,space, '     ',at_name_out(ii),at_occup(ii),at_occup_r(ii)
        enddo
      endif
    else 
      at_occup_r = nsuper_r/real(n_tot)
      if(i_traj==nt_min.and.ifile==nfile_min) then
        print *,space, 'Bulk concentrations:'
        do ii=1,n_atom
          print *,space, '     ',at_name_out(ii),at_occup_r(ii)
        enddo
      endif
    endif 		!i_traj==nt_min.and.ifile==nfile_min   

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
    if(angle(1)==.5*pi.and.angle(2)==.5*pi.and.angle(3)==.5*pi) then
      n_head = 3                        !number of header lines
    else
      n_head = 4                        !number of header lines for non-orthogonal lattices
    endif
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
   
   if(n_head==4) then
     write(header_record,nml=data_header_3)	
     i_rec = i_rec+1
     write(2,rec=i_rec) header_record
   endif

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

    if(j_verb==1.and.ifile==nfile_min)         call cpu_time(t2)
    if(j_verb==1.and.ifile==nfile_min) print *,space,'New binary output',t2-t1,' sec'
  
      backspace(1) 

      deallocate(at_pos_c,ind_at,at_ind_out,at_name_out)
      if(n_traj>=1) deallocate(at_veloc_c)
      if(j_shell_out==1) deallocate(at_pos_s)
      if(j_shell_out==1.and.n_traj>=1) deallocate(at_veloc_s)
      if(n_traj==2) deallocate(at_force_c)
      if(n_traj==2.and.j_shell_out==1) deallocate(at_force_s)

      if(ifile==nfile_max) exit trajectory_loop				
      if(ios_t<0) exit     !end of trajectory file
      ifile = ifile+1

    enddo frame_loop
    close(1) 

  enddo trajectory_loop
  deallocate(at_name,at_ind,e_kin,at_occup_r,nsuper_r,i_series)
    if(j_shell.eq.1) deallocate(e_kin_s)
  
  CALL SYSTEM_CLOCK (COUNT = sc_c2)
  print *,space, 'Trajectory files finished: ',i_save-n_save_min,' .dat files written in',(sc_c2-sc_c1)/sc_r,' sec (SYS)'
  write(9,*) 'Trajectory files finished: ',i_save-n_save_min,' .dat files written in',(sc_c2-sc_c1)/sc_r,' sec (SYS)'
  write(9,*) 

  stop
end program mp_dbin56




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

!     
subroutine up_case_2 (string_in,string_out)

  character(*), intent(in)	:: string_in
  character(*), intent(out)	:: string_out
  integer											:: j, nc
  character(len=26), parameter	:: lower = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter	:: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  do j = 1, len(string_in)
    nc = index(lower, string_in(j:j))
    if(nc > 0)then
       string_out(j:j) = upper(nc:nc)
    else
       string_out(j:j) = string_in(j:j)
    endif
  end do

end subroutine up_case_2	     

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
  
!          
subroutine down_case_2 (string_in,string_out)

  character(*), intent(in)	:: string_in
  character(*), intent(out)	:: string_out
  integer											:: j, nc
  character(len=26), parameter	:: lower = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter	:: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  do j = 1, len(string_in)
    nc = index(upper, string_in(j:j))
    if(nc > 0)then
       string_out(j:j) = lower(nc:nc)
    else
       string_out(j:j) = string_in(j:j)
    endif
  end do

end subroutine down_case_2	      
! 
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
   