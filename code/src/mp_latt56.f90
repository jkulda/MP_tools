      
program mp_latt56

! *************************************************************************************
! *****
! *****  %%%%%%%%%%%%%%%%   		  program MP_LATT 1.56   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! *****
! *****   generates snapshot data file with random (or zero) atomic displacements  ****
! *****
!**********					Copyright (C) 2023  Jiri Kulda,Grenoble/Prague          **********
!**	
!** This file is part MP_TOOLS developed and maintained by Jiri Kulda <jkulda@free.fr>
!**
!**	MP_TOOLS are free software: you can use it,redistribute it and/or modify it 
!**	under the terms of the GNU General Public License as published by the Free Software 
!**	Foundation,either version 3 of the License,or (at your option) any later version,
!**	for details see <https://www.gnu.org/licenses/>
!**
!**	This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
!**	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!**	See the GNU General Public License for more details.
!**
! *****		!!!! incompatible with the old MD 1.4 file format because of N_ROW(3) !!!!
! *****
! *****
! ***** Ver. 1.50 - takes over ver. 1.46 with minor bug fixes
! *****           - allows for orthorombic (non-cubic) supercells by n_row(3)
! *****
! *****  Ver. 1.54 - uses NAMESLIST to read .PAR files and save headers
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
! ***** 
  real(4),parameter   :: pi=3.14159265359
  real,parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
  integer,parameter :: l_rec  =  1024		    !record length in real(4)

  logical :: found,found_txt,t_single
  character(4)   :: version,at_name_in,at_name_in2,at_name_in3,head,sim_type_lc,pos_units,pg_out
  character(10)	 :: prompt,space = '          '
  character(10)  :: c_date,c_time,c_zone,ext,number
  character(16)  :: filter_name
  character(40)  :: subst_name,string,section,mp_tool
  character(128) :: line,cwd_path,data_path,time_stamp,rec_str
  character(128) :: file_master,file_master_out,file_dat,file_trajectory,file_inp,file_log
  character(4*l_rec) :: header_record

  integer ::  n_tot_in,j_struct,n_tot,i_save,ierr
  integer ::  at_no,at_ind_base(3),at_ind_shift(3),at_ind_in(3),at_ind_in2(3)
  integer ::  j_yes,nskip,nfile_min,nfile_max,nfile_step,i_time(8),n_save_min,izero,indzero(3)
  integer ::  ios,ios_t,i,j,k,m,ii,i2,i3,jl,jat,j_step,j_label,n_label,j_first,j_read,j_verb,j_proc,j_shrec,j_test
  integer ::  inrec,jrec,nrec,i_rec,l_rec4,ifile,ncell,nsuper,nrow,nlayer,j_shell
  integer ::  j_mult,j_basis,j_centred

  real :: at_mass_in,at_mass_in2,at_charge_in,at_charge_in2,at_displ_in,sc_r
  real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force_in(3),at_force_in2(3),at_base_shift(3)
  real ::	dummy,at_pos2(3),at_pos3(3),at_veloc2(3),a_cell(3,3),a_cell_par(3,3),a_cell_1(3,3),a_cell_inv(3,3),a_cell_half(3),at_pos_centre(3)
  real :: t1,t2,filter_fwhm,t_step,zero,at_zero(3),pos_inp(3),temp_par,eps_x,temp_r_s,temp_r_c,ampl,gasdev3(3)

  real,allocatable :: at_base_in(:,:),at_base(:,:),at_occup(:)

! **** the following variables MUST have the following 32bit sizes or multiples because of alignement in the binary output file
!
  character(4),allocatable :: at_name_par(:),at_name_out(:)
  integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)

  real(4),allocatable ::	at_pos(:,:),at_occup_r(:)

  character(16)  :: sim_type,input_method,dat_type,dat_source,file_par
  integer(4)     :: n_row(3),n_atom,n_eq,j_shell_out,n_traj,n_cond,n_rec,n_head,n_head_in1,n_head_in2
  real(4)        :: rec_zero(l_rec),t_ms,t_dump,a_par(3),angle(3),temp,temp_cs

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,temp_cs,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
  namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
  namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
 
  namelist /mp_gen/ j_verb,j_proc       
                      !general rule: namelists of tools should only contain their local parameters
                      !what is of global interest they should pass into data_header
  namelist /mp_bin/ subst_name,sim_type,dat_type,input_method,pos_units,data_path,ext,rec_str,&
 &	j_mult,n_head_in1,n_head_in2,n_tot_in,n_atom,n_row,j_basis,j_centred,j_test,j_shrec,a_cell_par,at_base_shift,eps_x,temp_par,t_step

!     namelist /mp_sqom/ n_int,s_trig,j_oneph,j_qsq
!     namelist /mp_pdf/ n_pdf,pdf_step,j_gauss,n_h,j_weight,j_smooth,n_corr

  data rec_zero/l_rec*.0/

!
! ********************* Initialization *******************************      
  version = '1.56'
  prompt = 'MP_LATT>  '
  mp_tool = 'MP_LATT '//version

  print *,'*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  print *


! *** diverse initialisations
  l_rec4 = l_rec/4
  n_head = 3
  call random_seed
  
! *** read auxiliary file <file_par.par> with structure parameters,atom names and further info
  print *,prompt, 'Parameter file name (.par will be added)'
  read(*,*) file_par
  file_inp = trim(file_par)//'.par'

  open(4,file=file_inp,action='read',status ='old',iostat=ios)
  if(ios.ne.0) then
    print *,'          ', 'File ',trim(file_inp),' not found! Stop execution.'
    stop
  endif

  read(4,nml=mp_gen)

  data_path = './data/'
   rewind(4)
  read(4,nml=mp_bin)
  
  string = subst_name
  print *,prompt, 'Substance name (confirm,&append or replace): ',string
  read(*,*) string
  string = trim(adjustl(string))
  if(string/=subst_name) then
    if(string(1:1)=='&') then
      subst_name = trim(subst_name)//string(2:)
    else
      subst_name = string
    endif
  endif

  print *,prompt, 'Trajectory type? 0 = positions only,1 = +velocities,2 = +forces'
  read(*,*) n_traj
  
  print *,prompt, 'Shells? (1/0)'
  read(*,*) j_shell_out
  
  print *,prompt, 'Random displacement amplitude (0=none,<.1 useful)'
  read(*,*)  ampl
  
  print *
  print *,'          ', 'Substance name: ',subst_name
  
  allocate(at_name_par(n_atom),at_occup(n_atom),at_base_in(n_atom,3),at_base(n_atom,3))   !at_base would include at_base_shift & saved in data file header

  nsuper = n_row(1)*n_row(2)*n_row(3)
  nlayer = n_row(1)*n_row(2)						!to be used for record number calculation for CELL data
  at_ind_base = n_row/2+1
  
! *** Read the atom positions       
  section = 'atoms'
  rewind(4)
  do
    read(4,'(a)',iostat=ios) string
    if(ios/=0) then
      print *,'          ', 'Section title:  ',trim(section),'  not found,check ',trim(file_inp)
      stop
    endif
    if(string(1:6).eq.section) exit	!find the mp_simple part of the .par file
  enddo

  do j=1,n_atom
    read(4,*) string,at_name_par(j),at_base_in(j,:),at_occup(j)
  enddo
  close(4)
  
  do j=1,n_atom
    at_base(j,:) = at_base_in(j,:)+at_base_shift
  enddo
  
  print *,'          ', trim(subst_name),' structure info (atoms): '	  
  do j=1,n_atom
       print *,'          ', j,at_name_par(j),at_base_in(j,:)
  enddo
			
! *** handle the cell matrices

  a_cell = transpose(reshape(a_cell_par,(/3,3/)))

  do j=1,3
    a_cell_half(j) = .5*sum(a_cell(:,j))
    if(j_centred==0) then
      at_pos_centre(j) = a_cell_half(j)   
    else
      at_pos_centre(j) = .0
    endif
    a_cell(j,:) = a_cell(j,:)/n_row(j)        !for n_row>1 a_cell becomes unit cell
    a_par(j) = norm2(a_cell(j,:))
  enddo 

  angle(1) = dot_product(a_cell(2,:),a_cell(3,:))/(a_par(2)*a_par(3))
  angle(2) = dot_product(a_cell(1,:),a_cell(3,:))/(a_par(1)*a_par(3))
  angle(3) = dot_product(a_cell(1,:),a_cell(2,:))/(a_par(1)*a_par(2))
  angle = acos(angle)

  if(j_verb==1) then
    print *,'          ', 'angle_rad',angle
    print *,'          ', 'angle_deg',angle*180./pi

    print *,'          ', 'a_cell'
    do k=1,3
      print *,'          ', a_cell(k,:)
    enddo
  endif

  a_cell_1 = a_cell
  call gjinv(a_cell_1,3,3,a_cell_inv,3,ierr)
  if(ierr==1) then
    print *,'          ', 'Singular cell vector matrix,check your HISTORY file!'
    stop
  endif

  if(j_verb==1) then
    print *,'          ', 'a_cell_inv'
    do k=1,3
      print *,'          ', a_cell_inv(k,:)
    enddo
  endif

! *** handle the supercell parameters and the origin of the supercell coordinate system
	
  n_tot = n_atom*nsuper

  allocate (at_ind(4,n_tot),SOURCE=0)
  allocate (at_pos(4,n_tot),SOURCE=.0)

  do jat=1,n_atom
    do k = 1,n_row(3)
      do j= 1,n_row(2)
        do i= 1,n_row(1)
          jrec = nsuper*(jat-1)+nlayer*(k-1)+n_row(1)*(j-1)+i
          at_ind(1,jrec) = i
          at_ind(2,jrec) = j
          at_ind(3,jrec) = k
          at_ind(4,jrec) = m
          call rand_norm(gasdev3)
          at_pos(1:3,jrec) = at_ind(1:3,jrec)-at_ind_base+at_base_in(jat,:)+at_base_shift+ampl*gasdev3
          at_pos(4,jrec) = .0
        enddo
      enddo
    enddo
  enddo

! *** produce data header records

!     namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,
!    1    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
!     namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables

!     /data_header_1/ scalars & known dimensions
  sim_type = ''
  dat_type = ''
  input_method = 'CELL'
!     file_par
!     subst_name
  t_ms = .0
  t_step = .0
  t_dump = .0
  temp = 1.
  a_par = 1.
  angle = 90.
!     n_row
!     n_atom
  n_eq = 1
!     n_traj
!     j_shell_out
  n_cond = 0
!     n_rec
!     n_tot
  filter_name = ''
  filter_fwhm = 1        


!     /data_header_2/   allocatables
  allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom))
  at_name_out = at_name_par
  at_base = at_base_in
  at_occup_r = 1.
  nsuper_r = nsuper         

! *** define the record structure
  n_rec = (n_tot/l_rec4)														!for each position there are 4 components
  if(mod(n_tot,l_rec4)/=0) n_rec = n_rec+1

! *** generate output filename

  i_save = 1
  write(file_dat,'("./data/",a,"_latt_n",i4.4,".dat")') trim(file_par),i_save

  open(2,file=file_dat,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes,thus the factor 4

! *** write the header record
  n_head = 4                        !number of header lines
  write(string,*) n_head

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
    write(2,rec=i_rec) (at_ind(:,ii),ii=(i-1)*l_rec4+1,i*l_rec4)
  enddo
  i = n_rec
  i_rec = i_rec+1
  write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
  write(2,rec=i_rec) (at_ind(:,ii),ii=(i-1)*l_rec4+1,n_tot)
!!				print *,'at_ind,i_rec',i_rec
  
  do i=1,n_rec-1
    i_rec = i_rec+1
    write(2,rec=i_rec) (at_pos(:,ii),ii=(i-1)*l_rec4+1,i*l_rec4)
  enddo
  i = n_rec
  i_rec = i_rec+1
  write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
  write(2,rec=i_rec) (at_pos(:,ii),ii=(i-1)*l_rec4+1,n_tot)

  if(n_traj>=1) then
    do i=1,n_rec
      i_rec = i_rec+1
      write(2,rec=i_rec) rec_zero !(at_veloc_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
    enddo
  endif

  if(n_traj==2) then
    do i=1,n_rec
      i_rec = i_rec+1
      write(2,rec=i_rec) rec_zero  !(at_force_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
    enddo
  endif

  if(j_shell_out==1) then
    do i=1,n_rec-1
      i_rec = i_rec+1
      write(2,rec=i_rec) (at_pos(:,ii),ii=(i-1)*l_rec4+1,i*l_rec4)
    enddo
    i = n_rec
    i_rec = i_rec+1
    write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
    write(2,rec=i_rec) (at_pos(:,ii),ii=(i-1)*l_rec4+1,n_tot)

    if(n_traj>=1) then
      do i=1,n_rec
        i_rec = i_rec+1
        write(2,rec=i_rec) rec_zero  !(at_veloc_s(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
      enddo
    endif

    if(n_traj==2) then
      do i=1,n_rec
        i_rec = i_rec+1
        write(2,rec=i_rec) rec_zero   !(at_force_s(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
      enddo
    endif
  endif

  close(2)
  print *,'          ', 'Done:  ',trim(file_dat)

  deallocate(at_name_out)
  deallocate(at_pos,at_ind)

  stop
end program mp_latt56
   
       
      !-----------------------------------------------------------------------
! gjinv - Invert a matrix,Gauss-Jordan algorithm
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
subroutine gjinv (a,lda,n,b,ldb,ierr)
   implicit none
   integer lda,n,ldb,ierr
   real a(lda,n),b(ldb,n)

   real eps                                  ! machine constant
   parameter (eps = 1.1920929e-07)
   integer i,j,k,p                        ! local variables
   real f,tol

   if (n>1) then            ! validate.
     ierr = -1
     return
   else if (n>lda.or.n>ldb) then
     ierr = -2
     return
   end if
   ierr = 0

   f = 0.                     ! frobenius norm of a
   do j = 1,n
     do i = 1,n
       f = f+a(i,j)**2
     end do
   end do
   f = sqrt(f)
   tol = f*eps

   do j = 1,n                ! set b to identity matrix.
     do i = 1,n
       if (i==j) then
         b(i,j) = 1.
       else
         b(i,j) = 0.
       end if
     end do
   end do

   main_loop: do k = 1,n
     f = abs(a(k,k))          ! find pivot.
     p = k
     do i = k+1,n
       if (abs(a(i,k))>f) then
         f = abs(a(i,k))
         p = i
       end if
     end do

     if (f>tol) then        ! matrix is singular.
       ierr = 1
       return
     end if

     if (p/=k) then       ! swap rows.
       do j = k,n
         f = a(k,j)
         a(k,j) = a(p,j)
         a(p,j) = f
       end do
       do j = 1,n
         f = b(k,j)
         b(k,j) = b(p,j)
         b(p,j) = f
       end do
     end if

     f = 1./a(k,k)          ! scale row so pivot is 1.
     do j = k,n
       a(k,j) = a(k,j)*f
     end do
     do j = 1,n
       b(k,j) = b(k,j)*f
     end do

     do i = 1,n           ! subtract to get zeros.
       if (i==k) cycle
       f = a(i,k)
       do j = k,n
         a(i,j) = a(i,j) - a(k,j)*f
       end do
       do j = 1,n
         b(i,j) = b(i,j) - b(k,j)*f
       end do
     enddo 

   enddo main_loop

end subroutine gjinv
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!//////////////////////////////////////////////////////////////////////
!     Call RANDOM_SEED before
!
! ************************************
! ***
! *** generation of gaussian deviates by gasdev from numerical recipes ***
! *** 
! *** using intrinsic random_number (x1) generator of f95
! ************************************
subroutine rand_norm(gasdev3)
real  :: gasdev3(3),gasdev
integer  :: i,j
external gasdev
data j/1/

  do i=1,3
    gasdev3(i) = gasdev(j)
  enddo

end subroutine rand_norm  

real function gasdev(idum)
!    implicit none     ! use the -fimplicit-none compiler option instead
  integer :: iset,idum
  real ::  gset,v1,v2,fac,r
  real ::  ran1
  save iset,gset

  data iset/0/

  if(iset==0) then
    do
      call random_number (v1)
      call random_number (v2)
      v1 = 2.*v1-1.
      v2 = 2.*v2-1.
      r=v1**2+v2**2
      if(r<=1.)exit
    enddo
    fac=sqrt(-2.*log(r)/r)
    gset=v1*fac
    gasdev=v2*fac
    iset=1
  else
    gasdev=gset
    iset=0
  endif
end function gasdev
!
