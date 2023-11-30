   
program mp_insp55

! *************************************************************************************
! *****
! *****  %%%%%%%%%%%%%%%%   		  program MP_INSP 1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! *****
! *****   inspects single atom properties in a binary MP_TOOLS data file  
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
! ***** %%%%%%%%%%%%%%%%           program MP_INSP 1.55   					 %%%%%%%%%%%%%%%%%%%%%%%%
! *****
! ***** Ver. 1.1 - original version
! ***** Ver. 1.2  
! ***** Ver. 1.3 - reads the extended header with n_atom, j_force and temp (on end positions for compatibility)
! *****	     		 - if atom forces are present, they are printed as well
! ***** Ver. 1.4 - dialogue revision
! *****					 - ./data subdirectory used for HIST and primary .bin files
! ***** Ver. 1.41 - minor bug corrections
! *****
! ***** Ver. 1.51 - non-cubic supercell with n_row(3) in header
! ***** Ver. 1.52 - provision for filenames without number
! *****
! ***** Ver. 1.53 - uses new semi-sequential binary file format 
! *****						- record length 4096, all data 32bit length
! *****						- at_mass is saved as at_veloc(4,j,k)
! *****						- at_charge is saved as at_pos(4,j,k)
! ***** Ver. 1.535 - uses n_traj of DL_POLY instead of j_force 
! *****
! *****
! *****   FROZEN on 02/10/2022 00:01 and committed as MP_INSP54
! *****
! ***** Ver. 1.55 - bug fixes 
! ***** 
! ***** Ver. 1.56 - accepts and displays triclinic unit cell geometry
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
! ***** atom positions are converted to and recorded in reduced lattice coordinates (x) 
! ***** 
    integer,parameter :: l_rec  =  1024		    !record length in real(4)
    
    character(4),allocatable      :: at_name_out(:)
    character(4)			:: version
    character(16)     :: sim_type,file_par,dat_type,dat_source,input_method,string,filter_name
    character(40)     :: subst_name,file_master,string_in,mp_tool
    character(60) 		:: file_dat,file_out
    character(1024)    :: header_record

    logical           ::  nml_in
    integer(4),allocatable ::  at_ind(:),ind_at(:),nsuper_r(:)    
    real(4),allocatable ::  at_pos_c(:),at_veloc_c(:),at_force_c(:),at_occup_r(:),at_base(:,:)     
    real(4),allocatable ::  at_pos_s(:),at_veloc_s(:),at_force_s(:)      
    
    integer(4) ::  at_no,n_traj,n_cond,n_atom,n_eq,j_shell_out,jpos_in,jrec_in,jrec_save,n_rec,n2_rec,l_rec4,i_rec
    integer(4) ::  ios,idum,i,j,k,jj,jat,ind,nfile,ncell,n_tot,nsuper,nlayer,nrow,n_row(3),jt,j_layer,j_row,j_pos
    
    integer(4) :: j_name,n_head

    real(4) :: a_par(3),angle(3),a_cell(3,3),a_cell_inv(3,3),t_ms,t_step,t_dump,filter_fwhm,temp

    namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,&
  &    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
    namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
    namelist /data_header_3/ a_cell,a_cell_inv                                !optional header containing non-orthogonal cell description
   
! ********************* Initialization *******************************      
  version = '1.56'
  mp_tool = 'MP_INSP '//version

  write(*,*) '*** Program ',trim(mp_tool),' ** Copyright (C) Jiri Kulda (2023) ***'
  write(*,*)
! *******************************************************************      

    n_traj = 0
    l_rec4 = l_rec/4
    dat_source = 'nn'
    version = 'nn'
    sim_type = 'nn'
    dat_type = 'nn'
    input_method = 'nn'
    t_ms = .0
    t_step = .0
    filter_name = 'nn'
    filter_fwhm = .0
    
    write(*,*) 'Master filename: '
    read(*,*) file_master 

! **** open a t-snapshot file and read its header 

    file_loop: do 

      jrec_save = 0

      write(*,*) 'snapshot file number (0 for no number):'
      read(*,*) jt

      if(jt==0)then
        file_dat = './data/'//trim(file_master)//'.dat'
      else
        if(jt<=9999) then
          write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),jt
        elseif(jt>=10000) then
          write(string,'(i8)') jt
          file_dat = './data/'//trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
        endif
      endif 			

      write(*,*)'Reading ',file_dat
      write(*,*) 
      open (1,file=file_dat,action='read',status ='old',access='direct',form='unformatted',recl=4*l_rec)

      i_rec = 1   			
      read(1,rec=i_rec) string_in
      if(string_in(1:8)=='MP_TOOLS') then      !new structure with namelist
        nml_in = .true.
        read(1,rec=i_rec) dat_source,version,string
        read(string,*) n_head
        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_1)			
        allocate(at_name_out(n_atom*n_eq),nsuper_r(n_atom*n_eq))
        allocate(at_occup_r(n_atom*n_eq),at_base(n_atom,3),SOURCE=.0)
        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_2)			!read the allocatables
        i_rec = i_rec+1   							
        if(n_head==4) then
          read(1,rec=i_rec) header_record
          read(header_record,nml=data_header_3)			!read the cell parameters
        endif
      else                                  !old w/o structure
        nml_in = .false.
        read(1,rec=i_rec) sim_type,file_par,t_ms,t_dump,temp,a_par,angle,n_row,n_atom,n_eq				
        allocate(at_name_out(n_atom*n_eq),at_occup_r(n_atom*n_eq),nsuper_r(n_atom*n_eq))
        read(1,rec=i_rec)& 
    &		sim_type,file_par,t_ms,t_dump,temp,a_par,angle,n_row,n_atom,n_eq,n_traj,j_shell_out,&
    &   n_cond,n_rec,n_tot,at_name_out,at_occup_r,nsuper_r																		!char(16): sim_type,file_par, char(4): at_name_out
        call up_case(sim_type)
        n_head = 1
        input_method = 'CELL'
        if(index(sim_type,'BULK')/=0) input_method = 'BULK'
      endif

      if(nml_in) then
      write(*,*) 'Code & version, no. of header lines:   ',dat_source,version,n_head	
      else
        write(*,*) 'Data source, version & number of header lines:  ','old header format'	
      endif
                    
      write(*,*) 'Substance name:                        ',subst_name
      write(*,*) 'Data & simulation type, input method:',dat_type,sim_type,input_method
      write(*,*) 'Time structure t_ms,t_step,t_dump:  ',t_ms,t_step,t_dump
      write(*,*) 'Shells, trajectory type, boundary cnd: ',j_shell_out,'  ',n_traj,'  ',n_cond
      if(filter_fwhm/=0.) write(*,*) 'Time filter name, fwhm:                ',trim(filter_name),filter_fwhm
      write(*,*) 'Supercell & temperature:            ',n_row,temp
      write(*,*) 'Unit cell parameter(3), angle(3):   ',a_par,'    ',angle
      if(n_head==4) then
      write(*,*) 
        write(*,*) '  a_cell                                                  a_cell_inv'
        do k=1,3
          write(*,*) a_cell(k,:),'    ',a_cell_inv(k,:)
        enddo
      endif

      write(*,*) 
      write(*,*) 'Atom numbers:                       ',n_atom,nsuper_r,n_tot
      write(*,*) 'Atoms & occupancies: '
      do j=1,n_atom
        write(*,*) '                   ',at_name_out(j),at_base(j,:),at_occup_r(j)
      enddo
      write(*,*) 
      
      nrow = n_row(1)
      nlayer = n_row(1)*n_row(2)						!only to be used for record number calculation
      nsuper = n_row(1)*n_row(2)*n_row(3)

      if(input_method=='FAST') input_method = 'CELL'

      allocate(at_ind(l_rec),at_pos_c(l_rec),at_veloc_c(l_rec),at_force_c(l_rec))
      if(j_shell_out==1) allocate(at_pos_s(l_rec),at_veloc_s(l_rec),at_force_s(l_rec))
      allocate(ind_at(n_atom))
      
      ind_at(1) = 0
      do jat = 2,n_atom
        ind_at(jat) = ind_at(jat-1)+nsuper_r(jat-1)
      enddo

      master_loop: do 	
        if(input_method=='CELL'.or.input_method=='FAST')	then
          write(*,*) 'jat, j_pos,j_row,j_layer (0 0 0 0 = new file number, 9 9 9 9 = END): '
          read(*,*) jat,i,j,k
          if (jat==0) exit master_loop
          if (jat==9.and.i==9.and.j==9.and.k==9) exit file_loop
          if(jat>n_atom) cycle

          ind = nsuper*(jat-1)+nlayer*(k-1)+n_row(1)*(j-1)+i
        else !input_method == 'BULK'
          write(*,*) 'jat, relative record_index (0 = new file number, -1 = END): '
          read(*,*) jat,ind
          if (jat==0) exit master_loop
          if (jat==-1) exit file_loop	
          if(jat>n_atom) cycle
          if(ind>nsuper_r(jat))	then
            write(*,*) 'IND out of range!'
            cycle master_loop
          endif			
          ind = ind+ind_at(jat)
        endif
        
        jrec_in = ind/(l_rec4)
        jpos_in = modulo(ind,l_rec4)
        if(jpos_in==0) then
          jpos_in = jpos_in+l_rec4
        else
          jrec_in = jrec_in+1
        endif
      
        if(jrec_in/=jrec_save) then			!read from disk only when necessary
          i_rec = jrec_in+n_head
          read(1,rec=i_rec) at_ind
        
          i_rec = jrec_in+n_rec+n_head
          read(1,rec=i_rec,iostat=ios) at_pos_c
          if(ios/=0) write(*,*)'at_pos_c: ios,i_rec',ios,i_rec
        
          if(index(sim_type,'static')==0.and.index(sim_type,'Static')==0.and.index(sim_type,'STATIC')==0) then
            i_rec = jrec_in+2*n_rec+n_head
            if(n_traj>=1) read(1,rec=i_rec) at_veloc_c
            if(n_traj==2) read(1,rec=jrec_in+3*n_rec+1) at_force_c  
      
            n2_rec = (2+n_traj)*n_rec+n_head
            if(j_shell_out==1) then
              read(1,rec=jrec_in+n2_rec) at_pos_s
              if(n_traj>=1) read(1,rec=jrec_in+n2_rec+n_rec) at_veloc_s
              if(n_traj==2) read(1,rec=jrec_in+n2_rec+2*n_rec) at_force_s  					
            endif
          endif
          jrec_save = jrec_in
        endif

        jj = 4*(jpos_in-1)
        write(*,*) 'Atom, type, mass, charge    ',at_name_out(jat),at_ind(jj+4),at_veloc_c(jj+4),at_pos_c(jj+4)   !at_ind(jj+4)=0 for unoccupied
        if(input_method=='CELL'.or.input_method=='FAST')	then
          write(*,*) 'Cell index            ',at_ind(jj+1:jj+3)
        else
          write(*,*) 'Atom index            ',at_ind(jj+1)					
        endif
        
        write(*,*) 'Atom position            ',at_pos_c(jj+1:jj+3)
        if(index(sim_type,'static')==0.and.index(sim_type,'Static')==0.and.index(sim_type,'STATIC')==0) then
          if(n_traj>=1) write(*,*) 'Atom velocity             ',at_veloc_c(jj+1:jj+3)
          if(n_traj==2) write(*,*) 'Atom force                ',at_force_c(jj+1:jj+3)
          write(*,*)
          if(j_shell_out==1) then
            write(*,*) trim(at_name_out(jat))//'_s',at_veloc_s(jj+4),at_pos_s(jj+4)
            write(*,*) 'Shell position                ',at_pos_s(jj+1:jj+3)
            if(n_traj>=1) write(*,*) 'Shell velocity                ',at_veloc_s(jj+1:jj+3)
            if(n_traj==2) write(*,*) 'Shell force               ',at_force_s(jj+1:jj+3)
          endif
        endif				
        write(*,*)         
      enddo master_loop
      close(1)

      deallocate(at_name_out,at_base,at_occup_r,nsuper_r,ind_at,at_ind,at_pos_c,at_veloc_c,at_force_c)
      if(j_shell_out==1) deallocate(at_pos_s,at_veloc_s,at_force_s)
    enddo file_loop

    stop
  end program mp_insp55
    
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

   
    