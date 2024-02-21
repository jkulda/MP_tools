  
program mp_lc56

! *************************************************************************************
! *****
! ***** %%%%%%%%%%%%%%%%   		  program MP_LC  1.56   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! *****
! *****               produces a linear combination of snapshots 
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
! ***** Ver. 1.55 - bug fixes 
! ***** 
! ***** Ver. 1.56 - handles triclinic unit cell geometry
! ***** 
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
  integer,parameter :: l_rec  =  1024		    !record length in real(4)

  character(4)   :: c_int(2),c_fil(2),version,head,atom
  character(10)	 :: prompt,space = '          '
  character(10)  :: pg_out,section,c_date,c_time,c_zone
  character(16)  :: sim_type_par,data_type,string16,f_name,filter_name
  character(40)  :: subst_name,file_master_a,file_master_b,file_master_out,time_stamp,int_mode,string,mp_tool
  character(60)  :: file_dat_a,file_dat_b,file_dat_out,file_log,line
  character(128) :: cwd_path
  character(4*l_rec):: header_record
  character(4*l_rec),allocatable:: header(:),dat_char(:)
  
  logical ::  nml_in,found,found_txt,t_single

  integer ::  j_head_in,j_verb,j_over,n_tot,n_head
  integer ::  i,j,k,ii,j0,jj,ind,ind2,i_rec,n_rec_0,n_rec_tot,i_time(8)
  integer ::  ifile,jfile,jfile_a,jfile_b,ifile0,ier,ios
  integer ::  nfile,nfile_a,nfile_b,nfile_min,nfile_min_out,nfile_min_a,nfile_step_a,nfile_max_a,nfile_min_b,nfile_step_b,nfile_max_b
  
  real    ::  coeff_a,coeff_b
  
  character(16)  :: sim_type,dat_type,input_method,file_par,dat_source
  integer(4)     :: n_row(3),n_atom,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,n_tot_in,idum
  real(4)        :: rec_zero(l_rec),t_ms,t_step,t_step_in,t_dump,t_dump_in,t0,t1,a_par(3),angle(3),a_par_pdf(3),temp,temp_cs,f_fwhm,filter_fwhm

  real(4),allocatable:: dat_file_a(:,:),dat_file_b(:,:)

  namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,temp_cs,a_par,angle,&
 &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
!!      namelist /data_header_2/at_name_out,at_base,at_occup_r,at_base,nsuper_r           !allocatables
 
!!      namelist /mp_gen/ j_verb,j_proc       
!!      namelist /mp_out/ j_weight,j_xray,j_logsc,j_txt,j_grid,pg_out       
!!      										!general rule: namelists of tools should only contain their local parameters
!!                          !what is of global interest they should pass into data_header
!!			namelist /mp_pdf/ n_pdf,pdf_step,a_par_pdf,n_h,j_acc,j_smooth,n_corr
!!			

! ********************* Initialization *******************************      
  version = '1.56'
  prompt = 'MP_LC>    '
  mp_tool = 'MP_LC '//version

  print *,'*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  print *,'***           linear combination of snapshots             ***'			
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
    open (9,file=file_log)
  endif

  write(9,*) 
  write(9,*) 
  write(9,*) trim(time_stamp),'  ',mp_tool,'  ',trim(cwd_path) 
  write(9,*) 
  
! *** Generate data file access
  print *,prompt, 'Data sequence A - file_master, file number min, step, max:  '
  read(*,*) file_master_a,nfile_min_a,nfile_step_a,nfile_max_a
        
!     print *,'Sequence A: data files number min, step, max: '
!     read(*,*)   nfile_min_a,nfile_step_a,nfile_max_a
  
  nfile_a = ((nfile_max_a-nfile_min_a)/nfile_step_a)+1
  
!		print *,'Sequence B: input data file_master: '
!		read(*,*) file_master_b 
!					
!     print *,'Sequence B: data files number min, step, max: '
!     read(*,*)   nfile_min_b,nfile_step_b,nfile_max_b
  
  print *,prompt, 'Data sequence B - file_master, file number min, step, max:  '
  read(*,*) file_master_b,nfile_min_b,nfile_step_b,nfile_max_b
  nfile_b = ((nfile_max_b-nfile_min_b)/nfile_step_b)+1

!      print *,'nfile_a,nfile_b',nfile_a,nfile_b
        
  if(nfile_a>1.and.nfile_b>1) then
    if(nfile_a/=nfile_b) then
      print *,space, 'Unequal length of file sequences: use the shorter one (1) or stop(0)? '
      read(*,*) jj      
      if(jj==1) then
        nfile = min(nfile_a,nfile_b)
      else
        stop
      endif
    else
      nfile = nfile_a      
    endif
  else
    nfile = max(nfile_a,nfile_b)      
  endif
!        print *,'nfile',nfile
      
! ***  Edit output file name
  do 
    print *,prompt, 'Output data file_master, file number min: '
    read(*,*) string,nfile_min_out
    file_master_out = trim(adjustl(string))
     
    jfile = nfile_min_out
    write(file_dat_out,'("./data/",a,"_n",i4.4,".dat")') trim(file_master_out),jfile
    inquire(file=file_dat_out,exist=found)
    
    if(.not.found) exit
    print *,space, 'File ',trim(file_dat_out),' already exists, overwrite(1/0)?'
    read(*,*) j_over
    if(j_over==1) exit
  enddo 
     
  print *,prompt, 'Linear combination coeffs for sequences A & B:'
  read(*,*) coeff_a,coeff_b
  
        call cpu_time(t0)

  j_over = 0            !do not automatically overwrite binary data files
  file_loop: do ifile=1,nfile

    if(ifile==1.or.nfile_a>1) then
      jfile_a = nfile_min_a+(ifile-1)*nfile_step_a
      if(jfile_a<=9999) then
        write(file_dat_a,'("./data/",a,"_n",i4.4,".dat")') trim(file_master_a),jfile_a
      elseif(jfile_a>=10000) then
        write(string,'(i8)') jfile_a
        file_dat_a = './data/'//trim(file_master_a)//'_n'//trim(adjustl(string))//'.dat'
      endif

      open(1,file=file_dat_a,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
      if(ios.ne.0) then
        print *,space, 'File A ',trim(file_dat_a),' not opened! IOS =',ios
        print *,space, 'Stop!'
        stop
      else
        if(ifile==1) print *,space, 'File A ',trim(file_dat_a)
      endif
    endif

    if(ifile==1.or.nfile_b>1) then
      jfile_b = nfile_min_b+(ifile-1)*nfile_step_b
      if(jfile_b<=9999) then
        write(file_dat_b,'("./data/",a,"_n",i4.4,".dat")') trim(file_master_b),jfile_b
      elseif(jfile_b>=10000) then
        write(string,'(i8)') jfile_b
        file_dat_b = './data/'//trim(file_master_b)//'_n'//trim(adjustl(string))//'.dat'
      endif

      open(2,file=file_dat_b,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
      if(ios.ne.0) then
        print *,space, 'File B ',trim(file_dat_b),' not opened! IOS =',ios
        print *,space, 'Stop!'
        stop
      else
        if(ifile==1) print *,space, 'File B ',trim(file_dat_b)
      endif
    endif

!! *** Read the header record 
    if(ifile==1) then
      i_rec = 1
      read(1,rec=i_rec) header_record
      head = header_record(1:4)
      call up_case(head)
      if(head == 'MP_T') then      !new structure with namelist
        nml_in = .true.
        read(1,rec=i_rec) dat_source,version,string16
        read(string16,*) n_head
!           print *,		'Input data:  ',dat_source,version,n_head
    
        allocate(header(n_head))
        header(i_rec) = header_record

        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_1)	
      else
        print *,space, 'header record old, wrong or missing'
        print *,space, trim(header_record)
        stop
      endif 

      if(trim(input_method)/='CELL') then
        print *,space, 'WARNING: the input data were not produced by the CELL input method,'
        print *,space, 'the algorithm as such may work, BUT as there is no guaranteed relationship'
        print *,space, 'between the data position in the .BIN file and the position of the related ,'
        print *,space, 'atom in the simulation box, the results when combining several snapshots'
        print *,space, 'may become COMPLETELY IMPREVISIBLE!'
        print *
        print *,prompt, 'Do you wish to continue anyways, at YOUR risks & perils? (1/0)'
        read(*,*) jj
        if(jj/=1) stop
      endif

!! *** Get the total record number
      i_rec = n_head	
      do
        i_rec = i_rec+1
        read(1,rec=i_rec,iostat=ios)  
        if(ios/=0) exit    
      enddo 
      n_rec_tot = i_rec-1
      n_rec_0 = n_head+n_rec

      jj = 0
      do i=n_head+1,n_rec_tot
        read(2,rec=i,iostat=ios) 
        if(ios/=0) then
          print *,space, 'File B: end of file at i_rec=',i
          print *,prompt, 'Continue with this record number (1) or stop(0)?'
          read(*,*) jj
          if(jj==0) then
            stop
          else
            n_rec_tot = i-1
            exit
          endif
         endif   			 
      enddo	

!         print *,'n_rec_0,n_rec_tot',n_rec_0,n_rec_tot

      allocate(dat_char(n_rec_0),dat_file_a(n_rec_tot-n_rec_0,l_rec),dat_file_b(n_rec_tot-n_rec_0,l_rec))
    endif
    
! *** read files A & B
    if(ifile==1.or.nfile_a>1) then
      do i=1,n_rec_0
        read(1,rec=i,iostat=ios) dat_char(i)
        if(ios/=0) then
          print *,space, 'File A: end of file',ifile,'at i_rec=',i
          print *,space, 'Stop!'
          stop
        endif   			 
      enddo
    endif

    do i=1,n_rec_tot-n_rec_0
      read(1,rec=n_rec_0+i,iostat=ios) dat_file_a(i,:)
      if(ios/=0) then
        print *,space, 'File A: end of file',ifile,'at i_rec=',n_rec_0+i
        print *,space, 'Stop!'
        stop
      endif   			 
    enddo	

    if(ifile==1.or.nfile_b>1) then
      do i=1,n_rec_tot-n_rec_0
        read(2,rec=n_rec_0+i,iostat=ios) dat_file_b(i,:)
        if(ios/=0) then
          print *,space, 'File B: end of file',ifile,'at i_rec=',n_rec_0+i
          print *,space, 'Stop!'
          stop
        endif   			 
      enddo	
    endif

    jfile = nfile_min_out+ifile-1
    if(jfile<=9999) then
      write(file_dat_out,'("./data/",a,"_n",i4.4,".dat")') trim(file_master_out),jfile
    elseif(ifile>=10000) then
      write(string,'(i8)') jfile
      file_dat_out = './data/'//trim(file_master_out)//'_n'//trim(adjustl(string))//'.dat'
    endif
        
!       inquire(file=file_dat_out,exist=found)
!       
!       if(found.and.j_over/=1) then
!         print *,'File ',trim(file_dat_out),' already exists, overwrite(1/0)?'
!         read(*,*) j_over
!         if(j_over==0) stop
!       endif          

    open(3,file=file_dat_out,access='direct',form='unformatted',recl=4*l_rec)

    if(ifile==1) then 
       print *,'                 (working ... may take a short while)'
    endif
  

! *** write the header record & at_ind

    do i=1,n_rec_0
      write(3,rec=i) dat_char(i)     !header & at_ind are copied from file A       
    enddo

! *** do the rest	
    do i=1,n_rec_tot-n_rec_0
      write(3,rec=n_rec_0+i) coeff_a*dat_file_a(i,:)+coeff_b*dat_file_b(i,:)	
    enddo
    close(3)

    if(ifile==1) print *,space,'A: ',trim(file_dat_a),'   ','B: ',trim(file_dat_b),'   ','Out: ',trim(file_dat_out)
    if(ifile<100.and.ifile==10*(ifile/10)) print *,space, '   ',trim(file_dat_a),'   ','   ',trim(file_dat_b),'   ','     ',trim(file_dat_out)
    if(ifile>=100.and.ifile==100*(ifile/100)) print *,space, '   ',trim(file_dat_a),'   ','   ',trim(file_dat_b),'   ','     ',trim(file_dat_out)

  enddo file_loop	
  
  print *,space, 'Finished: ',file_dat_out
  print *,space, 'Total of',nfile,' files written'  
                    
end program mp_lc56
   
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
 

