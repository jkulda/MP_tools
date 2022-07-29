      
      program mp_dbin53

C *************************************************************************************
C *****
C ***** %%%%%%%%%%%%%%%%   		  program MP_DBIN 1.53   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   converts ASCII input data (DISCUS or equivalent) to the MP_TOOLS binary form
C *****
C**********					Copyright (C) 2022  Jiri Kulda, Grenoble/Prague          **********
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
C ***** %%%%%%%%%%%%%%%%   			program MP_DBIN 1.53  				 %%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****		read ASCII trajectories from DISCUS output files
C *****
C *****
C ***** Ver. 1.50 - based on mp_bin ver. 1.46 with minor bug fixes
C ***** Ver. 1.51 - restructured data access: possibility to point towards ASCII data on an external volume
C ***** 					- possibility to specify trajectory filename extension
C ***** Ver. 1.53 - cycle over a series of trajectory files to create a contiguous series of binary snapshot files
C *****						- record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****           - CELL, QUICK and BULK data type
C *****
C ***** reads ASCII output from diverse simulations and writes binary direct-access files for each snapshot
C ***** file name convention:
C *****     <master_file>_n<snapshot number>.dat
C *****     example: H05_10K_n0001.dat from H05_10K.txt
C *****
C ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
C ***** record numbers are directly related to cell positions and atom types
C *****    jrec = nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
C *****
C ***** atom positions are converted to and recorded in reduced (x) coordinates
C ***** segment (octant) of atom displacement is identified, numbered (idom) and recorded
C ***** at_pos_c(i,j,k):  idom = 4*(sign(i)+1.)/2 + 2*(sign(j)+1.)/2 + (sign(k)+1.)/2 + 1
C ***** example: at_pos_c(i,j,k) = (-.13 -.05 -.07) >> idom = 1 (≈ -1 -1 -1)
C ***** example: at_pos_c(i,j,k) = ( .13  .05  .07) >> idom = 8 (≈  1  1  1)
C *****

      real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
      integer,parameter :: l_rec  =  1024		    !record length in real(4)

			logical :: found,found_txt,t_single,fixed_form,short_form
      character(4) :: at_name_in,at_name_in2,at_name_in3
      character(10) :: struct_name,c_date,c_time,c_zone,ext,number,data_type
      character(16) :: string,section
      character(128) :: line,cwd_path,data_path,file_dat_c,file_dat_s,file_trajectory,f_master,file_inp,file_log,time_stamp
      character(4),allocatable :: at_label(:)

      integer ::  at_no,at_ind_base(3),at_ind_in(3),izero,indzero(3),n_at_cell,j_struct,n_tot,i_atom,i_save
      integer ::  j_yes,j_ext,j_shell,nskip,nfile_min,nfile_max,nfile_step,i_time(8),n_save_min,n_row_save(3)
      integer ::  ios,i,j,k,ii,i2,i3,jl,jat,jhead,jrec,j_step,j_label,n_label,j_first,j_read,j_verb,j_shrec,j_test,j_sim
      integer ::  nrec,ifile,ncell,nsuper,nrow,nlayer,n_row_eff(3),n_site,l_rec4,i_rec,n_items
      integer ::  sc_c1,sc_c2,sc_m,nt_min,nt_max,it

      integer,allocatable ::  ind_l(:),i_check(:),i_site(:,:),ind_at(:)
      integer,allocatable,target :: i_series(:),at_ind_out(:)
      integer,pointer ::  jr(:)
      real,allocatable ::  sum_pos(:,:,:),ord(:),e_kin(:,:),e_kin_s(:,:),a_par_at(:,:),x_pos(:,:),at_occup(:)
      real,allocatable ::	at_mass_c(:),at_charge_c(:),at_displ_c(:)

      real :: at_mass_inp,at_charge_inp,at_displ_inp,sc_r
      real ::	at_pos_in(3),at_pos_min(3),at_pos_max(3),dummy
      real :: t1,t2,zero,at_zero(3),pos_inp(3),eps_x

C **** the following variables MUST have the following dimensions (4) or multiples because of alignement in the binary output file
C
      character(4),allocatable :: at_name_par(:),at_name_out(:)
      integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)

      real(4),allocatable ::	at_occup_r(:)
      real(4),allocatable ::	at_pos_c(:,:)

      character(16)  :: sim_type,file_title
      integer(4)     :: n_row(3),n_atom,n_eq,j_force,j_shell_out,n_traj,n_cond,idum,n_rec
      real(4)        :: rec_zero(l_rec),t_ms,t_step,a_par(3),angle(3),temp

      data rec_zero/l_rec*.0/
      data zero/.0/,at_zero/3*.0/,izero/0/,indzero/3*0/

C
C *****************************************************************************************
C
			write(*,*) '*** Program MP_DBIN 1.53 ** Copyright (C) Jiri Kulda 2022 ***'
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
			write(9,*) trim(time_stamp),'  MP_DBIN 1.53  ',trim(cwd_path)
			write(9,*) 
      
C *** read auxiliary file <file_title.par> with structure parameters, atom names and further info
      write(*,*) 'Parameter file name (.par will be added)'
      read(*,*) file_title
      file_inp = trim(file_title)//'.par'

			open(4,file=file_inp,action='read',status ='old',iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_inp),' not found! Stop execution.'
				stop
			endif

      write(9,*) 'Read parameter file:  ',trim(file_inp)
      write(*,*) 'Read parameter file:  ',trim(file_inp)

      section = 'mp_gen'
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
          stop
        endif
        if(string(1:6).eq.section) exit	!find the mp_gen part of the .par file
      enddo
			read(4,*) 		j_verb		! (0/1) verbose command line output
			read(4,*) 		!j_proc		!j_proc requested number of OpenMP parallel processes, 0 = automatic maximum
			read(4,*) 		!j_shrec=0: SHELLS data won't be recorded even if exist, (1 = SHELL data, if exist, will be recorded)
			read(4,*) 		!j_grid		! (0/1) grid overlay on PG_PLOT graphs 
			read(4,*) 		sim_type	! 'timestep' for MD, 'static' for static MD or DISCUS, must contain '_bulk' for BULK data type 
			rewind(4)

			data_type = 'cell'
			if(index(sim_type,'bulk')/=0) data_type = 'bulk'
			if(index(sim_type,'quick')/=0) data_type = 'quick'
			write(*,*) 'Data type: ',data_type

      section = 'mp_bin'
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
          stop
        endif
        if(string(1:6).eq.section) exit	!find the mp_bin part of the .par file
      enddo
102   format(a)

			read(4,*) struct_name
			read(4,*) n_atom,idum,		!number of atoms, displacement domain type (1=[100], 2=[110], 3= 111])
     1									eps_x 	!fractional coordinate tolerance range
			read(4,*) n_row_save 
			allocate(at_name_par(n_atom),at_label(n_atom))
			allocate(x_pos(n_atom,3),i_site(n_atom,n_atom),at_occup(n_atom))
			allocate(a_par_at(n_atom,3),at_occup_r(n_atom))
			allocate(nsuper_r(n_atom),ind_l(n_atom))
			n_eq = 1
			ind_l= 0
			j_label = 0
			n_site = 0
			
      if(data_type=='cell'.or.data_type=='quick')	then
				do j=1,n_atom
					read(4,*) string,at_name_par(j),x_pos(j,:),at_occup(j)
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
					at_label(j_label) = string
					ind_l(i) = 1
					i_site(j_label,ind_l(i)) = j
				enddo
				n_label = j_label
				write(*,*) 'Atom labels ',n_label
				do j=1,n_label
					write(*,*) ind_l(j),at_label(j),(i_site(j,i),i=1,ind_l(j))
				enddo	
			
				write(*,*) trim(struct_name),' structure info (atoms): '	  
				do j=1,n_atom
					write(*,*) j,at_name_par(j),x_pos(j,:)
				enddo
			else							! data_type = 'bulk'
				do j=1,n_atom
					read(4,*) at_label(j)
					ind_l(j) = j
				enddo
				n_label = n_atom
				write(*,*) 'Atom labels ',n_label
				do j=1,n_label
					write(*,*) ind_l(j),at_label(j)
				enddo	
			endif

			read(4,*) j_test		! number of snapshots to test the lattice parameter (1 minimum, 3 prudent; more superfluous) 
      read(4,'(a)',iostat=ios) data_path
      ii = index(data_path,'#')
      if(ii/=0) data_path = trim(data_path(1:ii-1))
      if(data_path(1:1)/='.'.and.data_path(1:1)/='/') data_path = 'data/'
			if(j_verb==1) write(*,*) 'Input data path:  ',trim(data_path)
      read(4,*,iostat=ios) ext
			j_ext = 1
      if(ios/=0.or.ext(1:1)=='#'.or.ext(1:3)=='ext'.or.ext(1:3)=='EXT') j_ext=0			!no file extension for ASCII input
			if(trim(ext)=='dat'.or.trim(ext)=='DAT') then
				write(*,*) 'PAR_file: the extension .DAT is reserved for the binary output ...'
				stop
			endif
			close(4)
			

C *** initialise variables from standard MP_TOOLS file header, that are not in use here
			l_rec4 = l_rec/4
			nskip = 50
			n_traj = 0		
			n_cond = 0
			j_force = 0
			j_shell_out = 0
			t_ms = .0
			t_step = .0
			temp = .0

C *** look for input file with MD trajectory
	  	master_loop: do
				write(*,*) 'Filename master for input data (e = end = END): '
				read(*,*) f_master
				if(trim(f_master)=='e'.or.trim(f_master)=='end'.or.trim(f_master)=='END') exit master_loop
				write(*,*) 'Read data files number min, max (step=1, 0 0 no numbers): '
				read(*,*) nt_min,nt_max 
				t_single = nt_min==0.and.nt_max==0

				if(t_single)then			
					nt_min = 1				!the message has passed, now make them usable to control the loop
					nt_max = 1
				endif
CC				write(*,*) 'Type in nominal temperature [K]:'
CC				read(*,*) temp

C *** input cycle over MD trajectory files
				j_read = 0

				file_loop: do it=nt_min,nt_max
					if(t_single) then
						file_trajectory = f_master
						if(j_ext==0) then
							file_inp = trim(data_path)//trim(adjustl(file_trajectory))
						else
							file_inp = trim(data_path)//trim(adjustl(file_trajectory))//'.'//trim(ext)
						endif
						inquire (file=file_inp,exist=found)
						if(.not.found) then
							write(*,*) 'Input data file ',trim(file_inp),' not found ...'
							cycle master_loop
						endif
					else
						write(number,'(i4.4)') it		!test a fixed_format counter (f_master0001)
						file_trajectory = trim(f_master)//trim(number)
						if(j_ext==0) then
							file_inp = trim(data_path)//trim(adjustl(file_trajectory))
						else
							file_inp = trim(data_path)//trim(adjustl(file_trajectory))//'.'//trim(ext)
						endif
						inquire (file=file_inp,exist=found)
						if(found) fixed_form = found
					endif

					if(.not.found) then
						if(it>=1.and.it<=9)       write(number,'(i1.1)') it			!test the short_format
						if(it>=10.and.it<=99)     write(number,'(i2.2)') it
						if(it>=100.and.it<=999)   write(number,'(i3.3)') it
						if(it>=1000.and.it<=9999) write(number,'(i4.4)') it
						file_trajectory = trim(f_master)//trim(number)
						if(j_ext==0) then
							file_inp = trim(data_path)//trim(adjustl(file_trajectory))
						else
							file_inp = trim(data_path)//trim(adjustl(file_trajectory))//'.'//trim(ext)
						endif
						inquire (file=file_inp,exist=found)
						if(found) short_form = found
					endif

					if(.not.found) then
						write(*,*) 'Input data file not found ...'
						cycle master_loop
					endif

					write(*,*) 'Input data file:  ',trim(file_inp)
					write(9,*) 'Input data file:  ',trim(file_inp)

					open (1,file=file_inp,action='read',status ='old',iostat=ios)
					if(ios.ne.0) then
						write(*,*) "Can't open the file ",trim(file_inp),'!'
						stop
					endif

					write(9,*) 'Reading input file:  ',trim(file_inp)
CC					if(j_verb==1) write(*,*) 'Reading input file:1  ',trim(file_inp)

C ***  read lines in file header     
				do i=1,nskip
					read(1,99) line
99      format(a80)
					if(index(line,'#')==1) cycle			
					if(index(line,'TITLE')/=0.or.index(line,'Title')/=0.or.index(line,'title')/=0) then
						backspace(1)
						read(1,*) line,file_title
						exit
					endif
				enddo
				rewind(1)
           						
			n_items = 0
			do i=1,nskip
				read(1,99) line
				if(index(line,'#')==1) cycle			
				if(index(line,'CELL')/=0.or.index(line,'Cell')/=0.or.index(line,'cell')/=0) then
					backspace(1)
					read(1,*) line,a_par,angle
					n_items = n_items+1
					exit
				endif
			enddo
			if(n_items==0) then
				write(*,*) 'The CELL input missing, STOP!'
				stop
			endif
			rewind(1)
      
			n_items = 0
			do i=1,nskip
				read(1,99) line
				if(index(line,'#')==1) cycle			
				if(index(line,'NCELL')/=0.or.index(line,'Ncell')/=0.or.index(line,'ncell')/=0) then
					backspace(1)
					if(data_type=='cell'.or.data_type=='quick') then
						read(1,*) line,n_row,n_at_cell
						if(n_row(1)/=n_row_save(1).or.n_row(2)/=n_row_save(2).or.n_row(3)/=n_row_save(3)) then
							write(*,*) 'Supercell dimensions differs in .PAR and in data:',n_row_save,n_row
							write(*,*) 'Type the correct N_ROW value:'
							read(*,*) n_row
						endif						
						nlayer = n_row(1)*n_row(2)						!only to be used for record number calculation
						nsuper = n_row(1)*n_row(2)*n_row(3)

						if(n_at_cell/=n_atom.and.data_type=='cell') then
							write(*,*) 'Number of atoms in unit cell differs in .PAR and in data:',n_atom,n_at_cell
							write(*,*) 'Type the correct N_ATOM value:'
							read(*,*) n_atom
						endif
						n_tot = n_atom*n_row(1)*n_row(2)*n_row(3)
						n_items = n_items+1
						exit
					else			!data_type == 'bulk'
						read(1,*) line,n_row,idum,n_tot,n_at_cell
						if(n_row(1)/=n_row_save(1).or.n_row(2)/=n_row_save(2).or.n_row(3)/=n_row_save(3)) then
							write(*,*) 'Supercell dimensions differs in .PAR and in data:',n_row_save,n_row
							write(*,*) 'Type the correct N_ROW value:'
							read(*,*) n_row
						endif						
						nlayer = 0						
						nsuper = 0

						if(n_at_cell/=n_atom) then
							write(*,*) 'Number of atoms in unit cell differs in .PAR and in data:',n_atom,n_at_cell
							write(*,*) 'Type the correct N_ATOM value:'
							read(*,*) n_atom
						endif
						n_items = n_items+1
						exit
					endif
				endif
			enddo
			if(n_items==0) then
				write(*,*) 'The NCELLS input is missing, STOP'
				stop
			endif

			nsuper_r = 0					!needed for accumulation of atom counts for occupancies
			write(*,*) 'n_atom,n_at_cell,n_tot,nsuper,nlayer',n_atom,n_at_cell,n_tot,nsuper,nlayer,n_row(1)
			rewind(1)
 

C ***  Prepare to read-in the text of a snapshot in one go

				allocate(at_pos_c(4,n_tot),source=.0)
				allocate(at_ind(4,n_tot),source=0)					!,at_displ_c(n_tot),i_check(n_tot)) !atom number jat is the 4th of the four indices
				allocate(at_ind_out(n_tot),i_series(n_tot))
				i_series = (/ (i, i = 1, n_tot) /)

				if(data_type/='bulk') then		!prefill the atom numbers so that empty positions conform
					do j=1,n_atom
						do i=1,nsuper
							at_ind(4,(j-1)*nsuper+i) = j
						enddo
					enddo
				endif

				if(it==nt_min) write(*,*) 'reading the 1st input file (takes a few seconds) ...'

        call cpu_time(t1)

				do i=1,nskip
					read(1,99) line
					if(index(line,'ATOMS')/=0.or.index(line,'Atoms')/=0.or.index(line,'atoms')/=0) exit
				enddo
				i_atom = 0
				at_pos_min = .0
				at_pos_max = .0

        read_loop: do          !swallow the snapshot

					if(data_type=='quick')then
          	read(1,*,iostat=ios) at_name_in,jat,at_pos_in
					else
	          read(1,*,iostat=ios) at_name_in,at_pos_in
					endif

					if(ios>0) then 
						write(*,*) 'Error on ASCII input IOS: ',ios
						backspace(1)
						read(1,'(a)')line
						write(*,*) line
						stop
          elseif(ios<0.or.at_name_in.eq.'end') then
CC 						if(j_verb==1) write(*,*) 'End of file IOS: ',ios
          	exit read_loop    !end of snapshot, end of file
          endif
					i_atom = i_atom+1
					
CC					if(i_atom==10000*(i_atom/10000)) write(*,*) i_atom

					jl = 0
					do j=1,n_atom
						if (at_name_in.eq.at_label(j)) then
							jl=j			!atom label found
							exit
						endif
					enddo

					if(jl.ne.j) then			!after "normal" loop exit
						write(*,*) 'atom ',at_name_in,' not found in .PAR'
						stop
					endif

					if(data_type=='bulk') then
						jat = jl				
						do k=1,3
							if(at_pos_in(k)<at_pos_min(k)) at_pos_min(k) = at_pos_in(k)
							if(at_pos_in(k)>at_pos_max(k)) at_pos_max(k) = at_pos_in(k)
							if(at_pos_in(k)<.0) at_pos_in(k) = at_pos_in(k)+n_row(k)					!bring the position into the box
							if(at_pos_in(k)>real(n_row(k))) at_pos_in(k) = at_pos_in(k)-n_row(k)							
						enddo
						jrec = i_atom
						at_pos_c(1:3,jrec) = at_pos_in-n_row/2
						at_pos_c(4,jrec) = .0							!could be atom charge if individualised

						at_ind(1,jrec) = jrec
						at_ind(2:3,jrec) = 0
						at_ind(4,jrec) = jat

					elseif(data_type=='cell') then
						if(ind_l(jl).eq.1) then
							jat = i_site(jl,1)	!atom name and fractional position is identified
						else
							do ii=1,ind_l(jl)
								pos_inp = at_pos_in-x_pos(i_site(jl,ii),:)
								jat = i_site(jl,ii)
								if(maxval(abs(pos_inp-anint(pos_inp))).lt.eps_x) exit !atom found
							enddo
						endif

						at_ind_in = anint(at_pos_in-x_pos(jat,:))+1				!they will serve as pointers to the right order of atom records
						at_pos_in = at_pos_in-n_row/2						        !now the supercell will be centred as from DL_POLY

						do k=1,3
							if(at_ind_in(k).eq.0) then
								at_ind_in(k) = at_ind_in(k)+n_row(k)
								at_pos_in = at_pos_in+n_row(k)
CC								write(*,*)'i,k,at_name,at_pos,at_ind',i,k,at_name(i),at_ind(i,:),at_pos_c(i,:)
							endif
							if(at_ind_in(k).eq.n_row(k)+1) then
								at_ind_in(k) = at_ind_in(k)-n_row(k)
								at_pos_in = at_pos_in-n_row(k)
CC								write(*,*)'i,k,at_name,at_pos,at_ind',i,k,at_name(i),at_ind(i,:),at_pos_c(i,:)
							endif
						enddo

						jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1) 					         
						at_ind(1:3,jrec) = at_ind_in
						at_ind(4,jrec) = jat
						at_pos_c(1:3,jrec) = at_pos_in
						at_pos_c(4,jrec) = .0

					else	!data_type=='quick'				         
CC            jat = mod((i_atom-1),n_atom)+1 !jat comes on input
            ii = (i_atom-1)/n_at_cell+1
						at_ind_in(3) = (ii-1)/nlayer+1
						at_ind_in(2) = mod((ii-1),nlayer)/n_row(1)+1
						at_ind_in(1) = mod(ii-1,n_row(1))+1

						jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1) 	
						at_ind(1:3,jrec) = at_ind_in
						at_ind(4,jrec) = jat

						do k=1,3
							if(at_pos_in(k)<0.) at_pos_in(k) = at_pos_in(k)+n_row(k)					!bring the position into the box
							if(at_pos_in(k)>real(n_row(k))) at_pos_in(k) = at_pos_in(k)-n_row(k)							
						enddo
						at_pos_c(1:3,jrec) = at_pos_in-n_row/2   !now the supercell will be centred as from DL_POLY
						at_pos_c(4,jrec) = .0
 					endif		!data_type

        	nsuper_r(jat) = nsuper_r(jat)+1        
				enddo read_loop

				n_tot = jrec
				j_read = j_read+1
CC				write(*,*) 'i_atom,jat,ii,at_ind_in',i_atom,jat,ii,at_ind_in
CC				write(*,*) at_pos_c(:,1)
				
C *** for 'BULK' the atom cloud is centred in the supercell box defined by n_row 
C *** now order the AT_IND and AT_POS arrays by atoms
				allocate(at_name_out(n_atom),ind_at(n_atom))


				if(data_type=='bulk') then
CC					do i=1,n_tot
CC						at_pos_c(1:3,i) =	at_pos_c(1:3,i) - .5*(at_pos_max-at_pos_min)
CC					enddo
CC					n_row_eff = int(at_pos_max-at_pos_min)+1
				
					write(*,*) 'Position min',at_pos_min			
					write(*,*) 'Position max',at_pos_max			
CC					write(*,*) 'Effective n_row',n_row_eff			
					write(*,*) 'Nominal n_row  ',n_row

					ind_at(1) = 0
					do jat = 2,n_atom
						ind_at(jat) = ind_at(jat-1)+nsuper_r(jat-1)
					enddo
				
					do i=1,n_tot
						jat = at_ind(4,i)
						ind_at(jat) = ind_at(jat)+1
						at_ind_out(ind_at(jat)) = at_ind(1,i)
					enddo
					jr => at_ind_out
				else
					jr => i_series
				endif			
				
        call cpu_time(t2)
				if(it==nt_min) write(*,*) '        total of',i_atom,' atoms read in',t2-t1,' sec'
					
C *** get the atom position occupation numbers and compare them with those from the .par file
				if(it==nt_min) then												!analyze in detail the 1st snapshot
					if(data_type=='cell'.or.data_type=='quick') then
						at_occup_r = (1.*nsuper_r)/nsuper
						at_name_out = at_name_par(1:n_atom)
						write(*,*) 'Occupancies: nominal 		real'
						do ii=1,n_atom
							write(*,*) '     ',at_name_out(ii),at_occup(ii),at_occup_r(ii)
						enddo
					else
						at_occup_r = nsuper_r/real(n_tot)
						at_name_out = at_label(1:n_atom)
						write(*,*) 'Bulk concentrations:'
						do ii=1,n_atom
							write(*,*) '     ',at_name_out(ii),at_occup_r(ii)
						enddo
					endif
				endif

				
C *** define the record structure
				n_rec = (n_tot/l_rec4)														!for each position there are 4 components
				if(mod(n_tot,l_rec4)/=0) n_rec = n_rec+1
																				
C *** generate output filename
				if(t_single)then
					file_dat_c = './data/'//trim(f_master)//'.dat' 	  	
				else
					write(number,'(i4.4)') it
					if(scan(f_master,'_')==len_trim(f_master))then
						file_dat_c = './data/'//trim(f_master)//trim(number)//'_n.dat'
					else
						file_dat_c = './data/'//trim(f_master)//'_n'//trim(number)//'.dat'
					endif 	  	
				endif
				write(*,*)trim(file_dat_c)


        call cpu_time(t1)
				open(2,file=file_dat_c,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4
				
C **** write the header
				i_rec = 1
				write(2,rec=i_rec)rec_zero				!fill zeros first
				write(2,rec=i_rec) 
     1		sim_type,file_title,t_ms,t_step,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,
     2    n_cond,n_rec,n_tot,at_name_out,at_occup_r(1:n_atom),nsuper_r(1:n_atom)	

C *** do the rest

				do i=1,n_rec-1
					i_rec = i_rec+1
					write(2,rec=i_rec) (at_ind(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
				enddo
				i = n_rec
				i_rec = i_rec+1
				write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
				write(2,rec=i_rec) (at_ind(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
				
				do i=1,n_rec-1
					i_rec = i_rec+1
					write(2,rec=i_rec) (at_pos_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
				enddo
				i = n_rec
				i_rec = i_rec+1
				write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
				write(2,rec=i_rec) (at_pos_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)

				close(2)
				close(1)

        call cpu_time(t2)
				if(it==nt_min) write(*,*)'        binary output file written in',t2-t1,' sec'

				deallocate(at_pos_c,at_ind,at_ind_out,at_name_out,ind_at,i_series) 
			enddo file_loop
			
      write(*,*) 'Trajectory files finished: ',j_read,' .dat files written'
      write(9,*) 'Trajectory files finished: ',j_read,' .dat files written'
      write(9,*) 
CC				CALL SYSTEM_CLOCK (COUNT = sc_c2)
CC				write(*,*) '  SYS_TIME',(sc_c2-sc_c1)/sc_r

			enddo master_loop
			
      stop
      end program mp_dbin53
      