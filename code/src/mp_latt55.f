      
      program mp_latt55

C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		  program MP_LATT 1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   generates a lattice snapshot data file with zero atomic displacements  ****
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
C *****   %%%%%%%%%%%%%%%%   			program MP_LATT 1.55  				 %%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****		!!!! incompatible with the old MD 1.4 file format because of N_ROW(3) !!!!
C *****
C *****
C ***** Ver. 1.50 - takes over ver. 1.46 with minor bug fixes
C *****           - allows for orthorombic (non-cubic) supercells by n_row(3)
C *****
C *****  Ver. 1.54 - uses NAMESLIST to read .PAR files and save headers
C *****
C ***** reads ASCII text output from MD simmulations (trajectory file)
C ***** writes binary direct-access files for each snapshot
C ***** depending on settings both CORE and SHELL records are read and recorded
C ***** file name convention:
C *****     <master_file>_n<snapshot number>.dat
C *****     example: H05_10K_n0001.dat from H05_10K.txt
C *****
C ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
C ***** record numbers are directly related to cell positions and atom types
C *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
C *****						1 going for the header record
C *****
C ***** atom positions are converted to and recorded in reduced lattice coordinates (x) 
C ***** 
      real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
      integer,parameter :: l_rec  =  1024		    !record length in real(4)

			logical :: found,found_txt,t_single
      character(4)   :: at_name_in,at_name_in2,at_name_in3,head,sim_type_lc,pos_units,pg_out
      character(10)  :: c_date,c_time,c_zone,ext,number
      character(16)  :: filter_name
      character(40)  :: subst_name,string,section
      character(128) :: line,cwd_path,data_path,time_stamp,rec_str
      character(128) :: file_master,file_master_out,file_dat,file_trajectory,file_inp,file_log
      character(4*l_rec) :: header_record

      integer ::  n_tot_in,j_struct,n_tot,i_save
      integer ::  at_no,at_ind_base(3),at_ind_shift(3),at_ind_in(3),at_ind_in2(3)
      integer ::  j_yes,nskip,nfile_min,nfile_max,nfile_step,i_time(8),n_save_min,izero,indzero(3)
      integer ::  ios,ios_t,i,j,k,m,ii,i2,i3,jl,jat,j_step,j_label,n_label,j_first,j_read,j_verb,j_proc,j_shrec,j_test
      integer ::  inrec,jrec,nrec,i_rec,l_rec4,ifile,ncell,nsuper,nrow,nlayer,j_shell
      integer ::  j_mult,j_basis,j_centred

      real :: at_mass_in,at_mass_in2,at_charge_in,at_charge_in2,at_displ_in,sc_r
      real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force_in(3),at_force_in2(3),at_base_shift(3)
      real ::	dummy,at_pos2(3),at_pos3(3),at_veloc2(3),a_cell(3,3),a_cell_par(3),a_cell_half(3),at_posentre(3)
      real :: t1,t2,filter_fwhm,t_step,zero,at_zero(3),pos_inp(3),temp_par,eps_x,temp_r_s,temp_r_c

      real,allocatable :: at_base_in(:,:),at_base(:,:),at_occup(:)

C **** the following variables MUST have the following 32bit sizes or multiples because of alignement in the binary output file
C
      character(4),allocatable :: at_name_par(:),at_name_out(:),version
      integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)

      real(4),allocatable ::	at_pos(:,:),at_occup_r(:)

      character(16)  :: sim_type,input_method,dat_type,dat_source,file_par
      integer(4)     :: n_row(3),n_atom,n_eq,j_shell_out,n_traj,n_cond,n_rec,n_head,n_head_in1,n_head_in2
      real(4)        :: rec_zero(l_rec),t_ms,t_dump,a_par(3),angle(3),temp

      namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,
     1    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
      namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables
     
      namelist /mp_gen/ j_verb,j_proc       
      										!general rule: namelists of tools should only contain their local parameters
                          !what is of global interest they should pass into data_header
			namelist /mp_bin/ subst_name,sim_type,dat_type,input_method,pos_units,data_path,ext,rec_str,
     1					j_mult,n_head_in1,n_head_in2,n_tot_in,n_atom,n_row,j_basis,j_centred,j_test,j_shrec,
     2          a_cell_par,at_base_shift,eps_x,temp_par,t_step

!     namelist /mp_sqom/ n_int,s_trig,j_oneph,j_qsq
!     namelist /mp_pdf/ n_pdf,pdf_step,j_gauss,n_h,j_weight,j_smooth,n_corr

      data rec_zero/l_rec*.0/

C
C *****************************************************************************************
C
			write(*,*) '*** Program MP_LATT 1.55 ** Copyright (C) Jiri Kulda (2023) ***'
      write(*,*)
      dat_source = 'MP_TOOLS'
      version = '1.55'


C *** diverse initialisations
			l_rec4 = l_rec/4
			n_head = 3
      
C *** read auxiliary file <file_par.par> with structure parameters, atom names and further info
      write(*,*) 'Parameter file name (.par will be added)'
      read(*,*) file_par
      file_inp = trim(file_par)//'.par'

			open(4,file=file_inp,action='read',status ='old',iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_inp),' not found! Stop execution.'
				stop
			endif

      read(4,nml=mp_gen)

      data_path = './data/'
       rewind(4)
      read(4,nml=mp_bin)
      
      string = subst_name
      write(*,*) 'Substance name (confirm, &append or replace): ', string
      read(*,*) string
      string = trim(adjustl(string))
      if(string/=subst_name) then
        if(string(1:1)=='&') then
          subst_name = trim(subst_name)//string(2:)
        else
          subst_name = string
        endif
      endif

      write(*,*) 'Trajectory type? 0 = positions only, 1 = +velocities, 2 = +forces'
      read(*,*) n_traj
      
      write(*,*) 'Shells? (1/0)'
      read(*,*) j_shell_out
      
      write(*,*) 
      write(*,*) 'Substance name: ', subst_name
			
			allocate(at_name_par(n_atom),at_occup(n_atom),at_base_in(n_atom,3),at_base(n_atom,3))   !at_base would include at_base_shift & saved in data file header

			nsuper = n_row(1)*n_row(2)*n_row(3)
			nlayer = n_row(1)*n_row(2)						!to be used for record number calculation for CELL data
			at_ind_base = n_row/2+1
			
C *** Read the atom positions       
      section = 'atoms'
      rewind(4)
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
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
      
      write(*,*) trim(subst_name),' structure info (atoms): '	  
      do j=1,n_atom
           write(*,*) j,at_name_par(j),at_base_in(j,:)
      enddo
			
C *** handle the supercell parameters and the origin of the supercell coordinate system
	
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
              at_pos(1:3,jrec) = at_ind(1:3,jrec)-at_ind_base+at_base_in(jat,:)+at_base_shift
              at_pos(4,jrec) = .0
            enddo
          enddo
        enddo
      enddo

C *** produce data header records

C     namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_step,t_dump,temp,a_par,angle,
C    1    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot,filter_name,filter_fwhm             !scalars & known dimensions
C     namelist /data_header_2/at_name_out,at_base,at_occup_r,nsuper_r           !allocatables

C     /data_header_1/ scalars & known dimensions
      sim_type = ''
      dat_type = ''
      input_method = 'CELL'
C     file_par
C     subst_name
      t_ms = .0
      t_step = .0
      t_dump = .0
      temp = 1.
      a_par = 1.
      angle = 90.
C     n_row
C     n_atom
      n_eq = 1
C     n_traj
C     j_shell_out
      n_cond = 0
C     n_rec
C     n_tot
      filter_name = ''
      filter_fwhm = 1        


C     /data_header_2/   allocatables
			allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom))
      at_name_out = at_name_par
      at_base = at_base_in
      at_occup_r = 1.
      nsuper_r = nsuper         

C *** define the record structure
				n_rec = (n_tot/l_rec4)														!for each position there are 4 components
				if(mod(n_tot,l_rec4)/=0) n_rec = n_rec+1

C *** generate output filename

  			i_save = 1
				write(file_dat,103) trim(file_par),i_save
103     format('./data/',a,'_latt_n',i4.4,'.dat')

				open(2,file=file_dat,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4

C *** write the header record
				n_head = 3                        !number of header lines
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

C *** do the rest	
				do i=1,n_rec-1
					i_rec = i_rec+1
					write(2,rec=i_rec) (at_ind(:,ii),ii=(i-1)*l_rec4+1,i*l_rec4)
				enddo
				i = n_rec
				i_rec = i_rec+1
				write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
				write(2,rec=i_rec) (at_ind(:,ii),ii=(i-1)*l_rec4+1,n_tot)
CC				write(*,*)'at_ind,i_rec',i_rec
				
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
      write(*,*) 'Done:  ',trim(file_dat)

			deallocate(at_name_out)
			deallocate(at_pos,at_ind)
			
      stop
      end program mp_latt55
   
   
