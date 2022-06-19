      
      program mp_bin53

C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		  program MP_BIN 1.53   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   converts MD trajectory data (DL_POLY or equivalent) to the MP_TOOLS binary form
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
C *****   %%%%%%%%%%%%%%%%   			program MP_BIN 1.53  				 %%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****		!!!! incompatible with the old MD 1.4 file format because of N_ROW(3) !!!!
C *****
C *****
C ***** Ver. 1.1 - reads properly the header information from the ASCII input file
C ***** Ver. 1.2 - looks for the header line 'timestep' instead of asking lines to skip
C ***** Ver. 1.3 - revised:
C *****		- only the "cryst_structure" type is kept, no substance-related data are left in the code
C *****		- looks for auxiliary file with structure info <title.par>
C *****		- the crystal structure is flexible (up to 20 atoms), the 3 oxygens in perovskite are still decodable
C *****		- cubic supercell is assumed, nrow etc. is then calculated, only n_atom to be supplied
C *****		- the header now contains n_atom, j_force and temp (on last position for compatibility)
C *****		- when atom forces (j_force = 1) are present in the trajectory, they can be included at the end of the binary records
C ***** Ver. 1.4 - dialogue revision
C *****		- ./data subdirectory used for trajectory and primary .bin files
C ***** Ver. 1.41 - cell parameter updated for each snapshot
C *****           - "history" replaced by "trajectory"
C ***** Ver. 1.42 - skips first nfile_min-1 snapshots and reads the sequence with a step nfile_step until nfile_max
C ***** Ver. 1.43 - refines the real temperature from kinetic energies of cores
C ***** Ver. 1.44 - identifies atoms according to the fractional positions in the PAR file 
C *****						- abandons the "cryst_structure" TYPE
C *****						- distinguishes low-T and hi-T behaviour of cores and shells: 
C *****															low-T strongly bound, k_B*T/2 = E_kin = E_kin_C + E_kin_S
C *****															hi-T independent movement, k_B*T/2 = E_kin = E_kin_C = E_kin_S
C ***** Ver. 1.45 - reads-in the whole snapshot in one go 
C ***** 					- uses orthorhombic supercell with N_ROW(3) 
C ***** 					- detects the NROW, if cubic, and refines the lattice parameters 
C *****           - at_ind(i,1:3) contains cell indices, at_ind(i,4) atom identifier (1:n_atom)
C *****           - writes a terminal record with fictive atom name 'End'
C ***** Ver. 1.46 - corrected: missing/wrong velocities input, test of j_force on output, atom names
C *****						- correct handling of partial and mixed occupancies 
C *****
C ***** Ver. 1.50 - takes over ver. 1.46 with minor bug fixes
C *****           - allows for orthorombic (non-cubic) supercells by n_row(3)
C ***** Ver. 1.51 - restructured data access: possibility to point towards ASCII data on an external volume
C ***** 					- possibility to specify trajectory filename extension
C ***** 					- all data arrays allocatable, no predefined array size limits
C ***** 					- supercell format in .PAR
C ***** Ver. 1.52 - cycle over a series of trajectory files to create a contiguous series of binary snapshot files
C *****						- waits until a subsequent trajectory file appears in its directory
C ***** Ver. 1.53 - new binary data format 
C *****						- record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****						- at_displ is not supported anymore
C *****						- cycle over trajectory files renewed
C *****           - CELL, QUICK and BULK data type
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
      character(4)   :: at_name_in,at_name_in2,at_name_in3
      character(10)  :: struct_name,c_date,c_time,c_zone,ext,number,data_type
      character(16)  :: string
      character(128) :: line,cwd_path,data_path,file_dat_c,file_dat_s,file_trajectory,file_inp,file_log,time_stamp
      character(128) :: file_dat_new,file_master,file_master_out

      integer ::  at_no,at_ind_base(3),at_ind_in(3),i_dom,i_dom_rec,n_dom,n_at_cell,j_struct,n_tot,i_atom,i_save
      integer ::  j_yes,j_ext,nskip,nfile_min,nfile_max,nfile_step,i_time(8),n_save_min,izero,indzero(3)
      integer ::  ios,ios_t,i,j,k,ii,i2,i3,jl,jat,jhead,j_step,j_label,n_label,j_first,j_read,j_corr,j_verb,j_shrec,j_test
      integer ::  i_rec,jrec,nrec,l_rec4,ifile,ncell,nsuper,nrow,nlayer,n_site,j_shell,n_row_save(3),n_row_eff(3)
      integer ::  sc_c1,sc_c2,sc_m,nt_min,nt_max,i_traj,j_mult,j_end

      real :: at_mass_in,at_mass_in2,at_charge_in,at_charge_in2,at_displ_in,sc_r
      real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force_in(3),at_force_in2(3)
      real ::	at_pos_min(3),at_pos_max(3),dummy,at_pos2(3),at_pos3(3),at_veloc2(3),a_cell(3,3)
      real :: t1,t2,zero,at_zero(3),pos_inp(3),temp,eps_x,temp_r_s,temp_r_c

      character(4),allocatable :: at_name(:),at_label(:)
      integer,allocatable ::  ind_l(:),i_check(:),nsuper_r(:),i_site(:,:),ind_at(:)
      integer,allocatable,target :: i_series(:),at_ind_out(:)
      integer,pointer ::  jr(:)
      real,allocatable ::  sum_pos(:,:,:),ord(:),e_kin(:,:),e_kin_s(:,:),a_par_at(:,:),x_pos(:,:),at_occup(:)

C **** the following variables MUST have the following 32bit sizes or multiples because of alignement in the binary output file
C
      character(4),allocatable :: at_name_par(:),at_name_out(:)
      integer(4),allocatable   :: at_ind(:,:)

      real(4),allocatable ::	at_pos_c(:,:),at_veloc_c(:,:),at_force_c(:,:),at_occup_r(:)
      real(4),allocatable ::	at_pos_s(:,:),at_veloc_s(:,:),at_force_s(:,:)

      character(16)  :: sim_type,file_title
      integer(4)     :: n_row(3),n_atom,n_eq,j_force,j_shell_out,n_traj,n_cond,idum,n_rec
      real(4)        :: rec_zero(l_rec),t_ms,t_step,a_par(3),angle(3),temp_r

      data rec_zero/l_rec*.0/

C
C *****************************************************************************************
C
			write(*,*) '*** Program MP_BIN 1.53 ** Copyright (C) Jiri Kulda (2019,2021,2022) ***'
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
			write(9,*) 
			write(9,*) trim(time_stamp),'  MP_BIN 1.53  ',trim(cwd_path)
      
C *** read auxiliary file <file_title.par> with structure parameters, atom names and further info
      write(*,*) 'Parameter file name (.par will be added)'
      read(*,*) file_title
      write(file_inp,101) trim(file_title)
101   format(a,'.par')

			open(4,file=file_inp,action='read',status ='old',iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_inp),' not found! Stop execution.'
				stop
			endif

      write(9,*) 'Read parameter file:  ',trim(file_inp)

      do
        read(4,'(a)') string
        if(string(1:6).eq.'mp_gen') exit	!find the mp_gen part of the .par file
      enddo
			read(4,*) 		j_verb		! (0/1) verbose command line output
			read(4,*) 		!j_proc		!j_proc requested number of OpenMP parallel processes, 0 = automatic maximum
			read(4,*) 		j_shrec		! j_shrec=0: SHELLS data won't be recorded even if exist, (1 = SHELL data, if exist, will be recorded)
			read(4,*) 		!j_grid		! (0/1) grid overlay on PG_PLOT graphs 
			read(4,*) 		sim_type	! 'timestep' for MD, 'static' for DISCUS, must contain '_bulk' for bulk data type
			rewind(4)

			data_type = 'cell'
			if(index(sim_type,'bulk')/=0) data_type = 'bulk'
CC			if(index(sim_type,'quick')/=0) data_type = 'quick'
			write(*,*) 'Data type: ',data_type
 
      do
        read(4,102) string
        if(string(1:6).eq.'mp_bin') exit	!find the mp_bin part of the .par file
      enddo
102   format(a)

			read(4,*) struct_name
			read(4,*) n_atom,n_dom,		!number of atoms, displacement domain type (1=[100], 2=[110], 3= 111])
     1									eps_x 	!fractional coordinate tolerance range
			read(4,*) n_row_save						!nrow(3) supercell format
			allocate(at_name_par(n_atom),at_label(n_atom),x_pos(n_atom,3),ind_l(n_atom),i_site(n_atom,n_atom),at_occup(n_atom))

			l_rec4 = l_rec/4
			t_single = .true.
			nt_min = 1
			nt_max = 1
			ind_l = 0
			j_ext = 1
			n_eq = 1									!for future use with multiplicities
			j_label = 0
			n_site = 0
			j_read = 0
			j_corr = 0
			n_row = n_row_save
			
      if(data_type=='cell')	then			
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
					at_label(j_label) = trim(string)
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
			else										! data_type = 'bulk'
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
			read(4,*) j_mult			!1 read multiple trajectories in a row, 0 just a single (long) one	
      if(j_mult/=0)	t_single = .false.
      close(4)

			angle = 90.				!assuming orthogonal lattice
			temp = .0


      write(*,*) 'Read snapshots number: min, step, max'
      read(*,*) nfile_min, nfile_step,nfile_max 
      write(9,*) 'Read snapshots number: min, step, max',nfile_min, nfile_step,nfile_max
      write(*,*) 'Saved snapshot numbers start:'
      read(*,*) n_save_min 
      write(9,*) 'Saved snapshot numbers start:',n_save_min   

C *** input cycle over MD trajectory files
			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)			
			nt_min = 1
			i_traj = 0
			ifile = 1
			i_save = n_save_min
			j_end = 0		

	  
			write(*,*) 'Master for MD trajectory filename: '
			read(*,*) file_master

      i = len_trim(file_master)
      if(file_master(i:i)=='_') then
				file_master_out = file_master(1:i-1)
      else
				file_master_out = file_master(1:i)
      endif

			if(.not.t_single) then 
				write(*,*) 'Read trajectories number (step=1): min, max'
				read(*,*) nt_min,nt_max 
			endif
			i_traj = nt_min-1

CC			trajectory_loop: do i_traj=nt_min,nt_max
			trajectory_loop: do 

				i_traj = i_traj+1
				if(i_traj>nt_max) exit trajectory_loop

				if(t_single) then
					file_trajectory = file_master
				else
					if(i_traj>=1.and.i_traj<=9)    write(number,'(i1.1)') i_traj
					if(i_traj>=10.and.i_traj<=99)  write(number,'(i2.2)') i_traj
					if(i_traj>=100.and.i_traj<=999)write(number,'(i3.3)') i_traj
					file_trajectory = trim(file_master)//trim(number)
				endif

				if(j_ext==0) then
					file_inp = trim(data_path)//trim(adjustl(file_trajectory))
				else
					file_inp = trim(data_path)//trim(adjustl(file_trajectory))//'.'//trim(ext)
				endif
				write(*,*) 'Input trajectory file:  ',trim(file_inp)

				t1 = .0
				do 
					inquire (file=file_inp,exist=found)
					if(.not.found) then
						if(t1==0.) write(*,*) 'File ',trim(file_inp),' not found! Waiting for it ...'			!if(t1==0.) 
						CALL SLEEP(10)
						t1 = t1+10.
					else
						exit
					endif
				enddo

				open (1,file=file_inp,action='read',status ='old',iostat=ios)
				if(ios.ne.0) then
					write(*,*) "Can't open the file ",trim(file_inp),'!'
					stop
				endif

				write(9,*) 'Reading MD trajectory file:  ',trim(file_inp)
					      
C ***  skip the first <nskip> lines until 'timestep' !normally they are two
				jhead = 0
				nskip = 20 ! number of records to test
				do i=1,nskip
					read(1,*) line
C        write(*,*) i,line
					jhead = index(line,sim_type(1:4))
					if(jhead.ne.0) exit
					if(i.eq.nskip) then
						write(*,*) 'MP_BIN: trajectory header not found'
						stop
					endif
				enddo 
				backspace(1)
    
C *** read the header of the first snapshot 
C
C *** typical header of a snapshot (4 lines):
C timestep      750125      552960           1           3  3.99999990E-04   270.04999    
C      194.2110000000        0.0000000000        0.0000000000          
C        0.0000000000      194.2110000000        0.0000000000          
C        0.0000000000        0.0000000000      194.2110000000          

      read(1,*) string,j_step,n_at_cell,n_traj,n_cond,t_ms,t_step		!t_ms is MD microstep, t_step is the snapshot time
																																			!n_at_cell total number of atoms in cells n_atom*nrow**3
      j_force = 0
      if(n_traj.eq.2) j_force = 1

			do j=1,3
				read(1,*) a_cell(j,:)
			enddo
			
			do j=1,3
			 a_par(j) = a_cell(j,j)/n_row(j)
			enddo

			if(i_traj==nt_min) then
C *** read 1st atom record         
			read(1,*) at_name_in,at_no,at_mass_in,at_charge_in,at_displ_in
				read(1,*) at_pos_in
				read(1,*) at_veloc_in
				if(j_force.eq.1) read(1,*) at_force_in

C *** read 2nd atom record         
			read(1,*) at_name_in2,at_no,at_mass_in,at_charge_in,at_displ_in
				read(1,*) at_pos2
				read(1,*) at_veloc2
				if(j_force.eq.1) read(1,*) at_force_in

C *** read 3rd atom record         
			read(1,*) at_name_in3,at_no,at_mass_in,at_charge_in,at_displ_in
				read(1,*) at_pos3
				read(1,*) at_veloc_in
				if(j_force.eq.1) read(1,*) at_force_in

C *** analyse the input
			if(at_name_in.eq.at_name_in2) then
				j_shell = 0				
				write(*,*) 'No shells'               !', supercell is:',nrow,'^3',' j_shell =,',j_shell
			else if(trim(at_name_in)//'s'.eq.at_name_in2.and.at_name_in.eq.at_name_in3) then
				j_shell = 1
				n_at_cell = n_at_cell/2
				write(*,*) 'Found a shell candidate   ',at_name_in2		!,' j_shell =',j_shell    !', supercell is:',nrow,'^3'
			else
				write(*,*) 'Strange primary data, please check them!'
				write(*,*) at_name_in,trim(at_name_in)//'s'
				write(*,*) at_name_in2
				write(*,*) at_name_in3
				stop
			endif

			if(at_veloc_in(1)==at_veloc2(1).and.at_veloc_in(2)==at_veloc2(2).and.at_veloc_in(3)==at_veloc2(3)) then
				write(*,*) 'Strange velocities:'
				write(*,*) at_veloc_in   
				write(*,*) at_veloc2          
				write(*,*) 'Type in nominal temperature [K]:'
				read(*,*) temp
			endif
			
      if(j_shell==1) then
				if(j_shrec==0) write(*,*)'Shell data NOT to be recorded (change this in .PAR)'
				if(j_shrec==1) write(*,*)'Shell data files WILL be recorded'
				if(j_shrec==1) write(9,*)'Shell data files WILL be recorded'
      endif
      j_shell_out = j_shell*j_shrec

      if(j_verb==1)  then
				write(*,*) 'j_shell =',j_shell
				write(*,*) 'sim_type = ',sim_type
				write(*,*) 'trajectory recording mode =', n_traj
				write(*,*) 'boundary conditions ', n_cond
      endif

				write(9,*) 'sim_type = ',sim_type
				write(9,*) 'trajectory recording mode =', n_traj
				write(9,*) 'boundary conditions ', n_cond
				write(9,*) 'Trajectory time start, step [ps]:',t_step,t_ms
				if(temp.ne.0.) write(9,*) 'Using nominal temperature [K] ',temp
			endif		!i_traj==nt_min

				write(*,*) 'Trajectory time start, step [ps]:',t_step,t_ms
			
      if(data_type=='cell') then
      	at_ind_base = n_row/2
				nsuper = n_row(1)*n_row(2)*n_row(3)
				nlayer = n_row(1)*n_row(2)						!only to be used for record number calculation
				n_tot = n_atom*n_row(1)*n_row(2)*n_row(3)
CC				if(n_tot/=n_at_cell) then
CC					write(*,*) 'Total number of atoms in .PAR and in DATA not coherent:',n_tot,n_at_cell
CC					stop 
CC				endif			
			else						!data_type=='bulk'
				n_tot = n_at_cell
			endif
      
			allocate(at_name(n_tot),SOURCE='    ')
			allocate(e_kin(n_atom,3),a_par_at(n_atom,3),at_occup_r(n_atom))
			allocate(i_check(n_tot),nsuper_r(n_atom),SOURCE=0) !atom number jat is the first of the four indices
			if(j_shell.eq.1)allocate(e_kin_s(n_atom,3),SOURCE=0.0)
			allocate(i_series(n_tot))
			i_series = (/ (i, i = 1, n_tot) /)
			j_read = 0
			j_corr = 0
			
			rewind(1)
     				
C *** cycle over snapshots, each to be saved in a separate binary file

      file_loop: do 				!we have to go through all the snapshots in the .TXT input

				if(ifile>nfile_max) exit trajectory_loop
				
				if(j_end/=0)then
					deallocate(at_pos_c,at_veloc_c,at_ind)
					if(j_shell==1) deallocate(at_pos_s,at_veloc_s)
					if(j_force==1) deallocate(at_force_c)
					if(j_force==1.and.j_shell==1) deallocate(at_force_s)
				endif
				j_end = 1			!at normal end of file_loop j_end will become 0

				allocate(at_ind(4,n_tot),SOURCE=0)
				allocate(at_pos_c(4,n_tot),at_veloc_c(4,n_tot),SOURCE=0.0)
				if(j_force.eq.1) allocate (at_force_c(4,n_tot),SOURCE=0.0)
				if(j_shell.eq.1) then
					allocate(at_pos_s(4,n_tot),at_veloc_s(4,n_tot),SOURCE=0.0)
					if(j_force.eq.1) allocate (at_force_s(4,n_tot),SOURCE=0.0)
				endif
			
				e_kin = .0
				if(j_shell.eq.1) e_kin_s = .0
				nsuper_r = 0
				
        do 						!i=1,nskip
          read(1,*,iostat=ios) line
          if(ios<0) then   !end of file
						write(*,*) 'End of trajectory file: ',ios
						exit file_loop 
					endif
CC          if(ios>0) cycle     !corrupt input data, go on
          if(ios>0) then
          	write(*,*) 'Input error, corrupt trajectory data?'
          	stop
          endif
          jhead = index(line,sim_type(1:4))
          if(jhead.ne.0) exit
        enddo 
        
 
C *** skip the snapshots that are not requested to be read
        if((ifile.lt.nfile_min).or.(mod(ifile-nfile_min,nfile_step).ne.0)) then  !skip the snapshot
        	ifile = ifile+1
          cycle file_loop         ! cycle the file_loop
        else
	        backspace(1)	
	        j_read = j_read+1
        endif                    
 
C *** no need to re-read the first header line but
C *** take fresh lattice parameters for each snapshot - they may evolve
        read(1,*) string,j_step,n_at_cell,n_traj,n_cond,t_ms,t_step
        do j=1,3
          read(1,*) a_cell(j,:)
 		 			a_par(j) = a_cell(j,j)/n_row(j)
       enddo  
        
      
    
				if(i_traj==nt_min.and.ifile==nfile_min) then
					write(*,*) 'a_par estimate =',a_par
					write(*,*) 'OK? (1/0)'
					read(*,*)	j_yes

					if(j_yes.ne.1) then
						write(*,*) 'input a better A_PAR(1:3)'
						read(*,*)  a_par
					endif
				endif

C ***  read-in the text of a snapshot in one go
				if(ifile==nfile_min) write(*,*) 'reading the 1st snapshot (takes a few seconds) ...'

        call cpu_time(t1)
				
				i_atom = 0
				at_pos_min = .0
				at_pos_max = .0

        read_loop: do          !swallow the snapshot

C *** first read the CORE data  

          read(1,*,iostat=ios_t) at_name_in,at_no,at_mass_in,at_charge_in,at_displ_in
CC						if(ios/=0) write(*,*) 'Error on ASCII input IOS: ',ios
CC						if(ios>0) cycle file_loop     !corrupt input data, jump out to new snapshot

          if(ios_t<0.or.at_name_in.eq.sim_type(1:4)) then
CC          	write(*,*) 'End of snapshot'
CC          	backspace(1)
          	exit read_loop    !end of snapshot, end of file
          endif
					i_atom = i_atom+1

          read(1,*) at_pos_in
          read(1,*) at_veloc_in
          if(j_force.eq.1) read(1,*) at_force_in

					at_pos_in = at_pos_in/a_par

C
C *** now read the SHELL data if needed
          if(j_shell.eq.1) then
          	read(1,*) at_name_in2,at_no,at_mass_in2,at_charge_in2,at_displ_in
						if(trim(at_name_in)//'s'.ne.at_name_in2) then
							write(*,*)'wrong core-shell sequence',i_atom,at_name(i_atom),at_name_in
							stop
						endif
						read(1,*) at_pos_in2
						read(1,*) at_veloc_in2
						if(j_force.eq.1) read(1,*) at_force_in2
						at_pos_in2 = at_pos_in2/a_par
					endif		!j_shell


C ***  treat the CORE data and get the right labels & positions 
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
						enddo
						jrec = i_atom

						at_ind(1,jrec) = jrec
						at_ind(2:3,jrec) = 0
						at_ind(4,jrec) = jat

					else			!data_type=='cell'
					
						if(ind_l(jl).eq.1) then
							jat = i_site(jl,1)	!atom name and fractional position is identified
						else
							do ii=1,ind_l(jl)
								pos_inp = at_pos_in-x_pos(i_site(jl,ii),:)
								jat = i_site(jl,ii)
								if(maxval(abs(pos_inp-anint(pos_inp))).lt.eps_x) exit !atom found
							enddo
						endif
          
                 
CC *** calculate cell indices ix,iy,iz from atomic positions shifted to cell origin       

						at_ind_in = anint(at_pos_in-x_pos(jat,:))+at_ind_base+1	!they will serve as pointers to the right order of atom records
				
						do k=1,3
							if(at_ind_in(k)==0) then
								at_ind_in(k) = at_ind_in(k)+n_row(k)
								at_pos_in(k) = at_pos_in(k)+n_row(k)
 								if(j_shell_out.eq.1) at_pos_in2(k) = at_pos_in2(k)+n_row(k)
							endif
							if(at_ind_in(k)==n_row(k)+1) then
								at_ind_in(k) = at_ind_in(k)-n_row(k)
								at_pos_in(k) = at_pos_in(k)-n_row(k)
 								if(j_shell_out.eq.1) at_pos_in2(k) = at_pos_in2(k)-n_row(k)
							endif
						enddo
					
						jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1)
						
						if(jrec>n_tot) then
						  write(*,*) 'jrec,at_ind_in',jrec,at_ind_in
						  write(*,*) 'at_pos_in',at_pos_in,at_pos_in2
						  write(*,*) 'n_tot,nsuper,nlayer,n_row',n_tot,nsuper,nlayer,n_row
						endif

						at_ind(1:3,jrec) = at_ind_in
						at_ind(4,jrec) = jat					
					endif		!data_type
				
					at_pos_c(1:3,jrec) = at_pos_in
					at_pos_c(4,jrec) = at_charge_in															 
					at_veloc_c(1:3,jrec) = at_veloc_in
					at_veloc_c(4,jrec) = at_mass_in							
					if(j_force.eq.1) at_force_c(1:3,jrec) = at_force_in

					if(j_shell_out==1) then				!only if the shell data are going to be recorded
						at_pos_s(1:3,jrec) = at_pos_in2
						at_pos_s(4,jrec) = at_charge_in2							
						at_veloc_s(1:3,jrec) = at_veloc_in2
						at_veloc_s(4,jrec) = at_mass_in2							
						if(j_force.eq.1) at_force_s(1:3,jrec) = at_force_in2
					endif

C *** accumulate the occupation number and the kinetic energy to refine the real temperature
					nsuper_r(jat) = nsuper_r(jat)+1
					do k = 1,3
						e_kin(jat,k) = e_kin(jat,k) + at_mass_in*at_veloc_in(k)**2			!at_mass_c(i)=at_veloc_c(4,i)
						if(j_shell==1) e_kin_s(jat,k) = e_kin_s(jat,k) + at_mass_in2*at_veloc_in2(k)**2
					enddo
				enddo read_loop
	
	
				allocate(at_ind_out(n_tot),at_name_out(n_atom),ind_at(n_atom))
			
				if(data_type=='bulk') then
					n_row_eff = int(at_pos_max-at_pos_min)+1
				
					if(i_traj==nt_min.and.ifile==nfile_min) then												!analyze in detail the 1st snapshot
						write(*,*) 'Position min',at_pos_min			
						write(*,*) 'Position max',at_pos_max			
						write(*,*) 'Effective n_row',n_row_eff			
						write(*,*) 'Nominal n_row  ',n_row
					endif
CC					write(*,*) 'n_atom,n_tot',n_atom,n_tot
					
					ind_at(1) = 0
					do jat = 2,n_atom
						ind_at(jat) = ind_at(jat-1)+nsuper_r(jat-1)
					enddo
CC				write(*,*) 'ind_at',ind_at
				
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
					write(*,*) '1st snapshot: total of',i_atom,' atoms read in',t2-t1,' sec'
				endif
 				
  				
				if(temp==.0) then
					do j=1,n_atom
						e_kin(j,:) = e_kin(j,:)/(nsuper_r(j))
						if(j_shell.eq.1) e_kin_s(j,:) = e_kin_s(j,:)/(nsuper_r(j))
					enddo

					if(j_verb==1.and.i_traj==nt_min.and.ifile==nfile_min) then
						write(*,*)
						write(*,*) 'Cores E_kin(jat,:)',(e_kin(jat,:),jat=1,n_atom)
						if(j_shell.eq.1) write(*,*) 'Shells E_kin(jat,:)',(e_kin_s(jat,:),jat=1,n_atom)
					endif

C *** get the true temperature
					temp_r_c = sum(e_kin(1:n_atom,:))/(n_atom*3*k_B) !the true core temperature

					if(j_shell.eq.1) then
						temp_r_s = sum(e_kin_s(1:n_atom,:))/(n_atom*3*k_B) !the true temperature
						if (abs(temp_r_c-temp_r_s).le..1*temp_r_c) then
							temp_r = (temp_r_c+temp_r_s)*.5						!hi-T limit: independent C and S vibrations
							if(j_read.le.j_test) write(*,*) 'Hi-T limit: independent C and S vibrations'
						else
							temp_r = temp_r_c+temp_r_s	!low-T limit: strongly bound C and S vibrations
							if(j_read.le.j_test) write(*,*) 'Low-T limit: strongly bound C and S vibrations'
						endif        
						if(j_read.le.j_test) write(*,*) 'Real temperature: core/shell/total ',temp_r_c,temp_r_s,temp_r
						if(j_read.le.j_test) write(9,*) 'Real temperature: core/shell/total ',temp_r_c,temp_r_s,temp_r
					else
						temp_r = temp_r_c
						temp_r_s = .0
						if(j_read.le.j_test) write(*,*) 'Real temperature: cores only ',temp_r
						if(j_read.le.j_test) write(9,*) 'Real temperature: cores only ',temp_r
					endif
				else
					temp_r = temp
					write(*,*) 'Using nominal temperature [K] ',temp
				endif
  				
      
C *** First snapshot in the series only
C *** get the atom position occupation numbers and compare them with those from the .par file
											!analyze in detail the 1st snapshot
				if(data_type=='cell') then
					at_occup_r = (1.*nsuper_r)/nsuper
					at_name_out = at_name_par(1:n_atom)
					if(i_traj==nt_min.and.ifile==nfile_min) then
						write(*,*) 'Occupancies: nominal 		real'
						do ii=1,n_atom
							write(*,*) '     ',at_name_out(ii),at_occup(ii),at_occup_r(ii)
						enddo
					endif
				else 
					at_occup_r = nsuper_r/real(n_tot)
					at_name_out = at_label(1:n_atom)
					if(i_traj==nt_min.and.ifile==nfile_min) then
						write(*,*) 'Bulk concentrations:'
						do ii=1,n_atom
							write(*,*) '     ',at_name_out(ii),at_occup_r(ii)
						enddo
					endif
				endif 		!i_traj==nt_min.and.ifile==nfile_min

 
C *** define the record structure
				n_rec = (n_tot/l_rec4)														!for each position there are 4 components
				if(mod(n_tot,l_rec4)/=0) n_rec = n_rec+1
																
CC				if(j_verb==1.and.ifile==nfile_min) write(*,*) 'n_tot,l_rec4,n_rec',n_tot,l_rec4,n_rec			

C *** generate output filename
				write(file_dat_c,103) trim(file_title),i_save
103     format('./data/',a,'_n',i4.4,'.dat')
				write(*,*)trim(file_dat_c)

        call cpu_time(t1)
				open(2,file=file_dat_c,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4

C *** write the header record
				i_rec = 1
				write(2,rec=i_rec)rec_zero				!fill zeros first
				write(2,rec=i_rec) 
     1		sim_type,file_title,t_ms,t_step,temp_r,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,
     2    n_cond,n_rec,n_tot,at_name_out,at_occup_r(1:n_atom),nsuper_r(1:n_atom)		

C *** do the rest	
				do i=1,n_rec-1
					i_rec = i_rec+1
					write(2,rec=i_rec) (at_ind(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
				enddo
				i = n_rec
				i_rec = i_rec+1
				write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
				write(2,rec=i_rec) (at_ind(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)
CC				write(*,*)'at_ind,i_rec',i_rec
				
				do i=1,n_rec-1
					i_rec = i_rec+1
					write(2,rec=i_rec) (at_pos_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
				enddo
				i = n_rec
				i_rec = i_rec+1
				write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
				write(2,rec=i_rec) (at_pos_c(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)

				do i=1,n_rec-1
					i_rec = i_rec+1
					write(2,rec=i_rec) (at_veloc_c(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
				enddo
				i = n_rec
				i_rec = i_rec+1
				write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
				write(2,rec=i_rec) (at_veloc_c(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)

				if(j_force==1) then
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

					do i=1,n_rec-1
						i_rec = i_rec+1
						write(2,rec=i_rec) (at_veloc_s(:,jr(ii)),ii=(i-1)*l_rec4+1,i*l_rec4)
					enddo
					i = n_rec
					i_rec = i_rec+1
					write(2,rec=i_rec)rec_zero										!the last one has to be padded by 0
					write(2,rec=i_rec) (at_veloc_s(:,jr(ii)),ii=(i-1)*l_rec4+1,n_tot)

					if(j_force==1) then
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
CC  				write(*,*)'at_shells,i_rec',i_rec

  			close(2)

  			i_save = i_save+1
				ifile = ifile+1
 
CC				if(j_verb==1.and.ifile==nfile_min)         call cpu_time(t2)
CC				if(j_verb==1.and.ifile==nfile_min) write(*,*)'New binary output',t2-t1,' sec'
			
        	backspace(1) 

					deallocate(at_pos_c,at_veloc_c,at_ind,ind_at,at_ind_out,at_name_out)
					if(j_shell==1) deallocate(at_pos_s,at_veloc_s)
					if(j_force==1) deallocate(at_force_c)
					if(j_force==1.and.j_shell==1) deallocate(at_force_s)

        	j_end = 0 
        	if(ios_t<0) exit     !end of trajectory file
				enddo file_loop
				close(1)
CC				write(*,*) 'End file_loop, i_traj',i_traj
				deallocate(at_name,e_kin,a_par_at,at_occup_r,i_check,nsuper_r,i_series)
				if(j_shell.eq.1) deallocate(e_kin_s)

			enddo trajectory_loop
			
			CALL SYSTEM_CLOCK (COUNT = sc_c2)
CC			write(*,*) '  SYS_TIME',(sc_c2-sc_c1)/sc_r
      write(*,*) 'Trajectory files finished: ',i_save-1,' .dat files written in',(sc_c2-sc_c1)/sc_r,' sec (SYS)'
      write(9,*) 'Trajectory files finished: ',i_save-1,' .dat files written in',(sc_c2-sc_c1)/sc_r,' sec (SYS)'
      write(9,*) 

      stop
      end program mp_bin53
      