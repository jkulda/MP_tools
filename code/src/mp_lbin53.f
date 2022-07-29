      
      program mp_lbin53

C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		  program MP_LBIN 1.53   		 %%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   converts MD trajectory data (LAMMPS or equivalent) to the MP_TOOLS binary form
C *****
C**********					Copyright (C) 2022  Jiri Kulda, Grenoble/Prague          **********
C *************************************************************************************
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
C ***** %%%%%%%%%%%%%%%%   			program MP_LBIN 1.53  				 %%%%%%%%%%%%%%%%%%%%%%%%
C *****	Reads LAMMPS dump files and records binary data files (based on MB_LBINARY44.F)
C *****
C ***** Ver. 1.44 - identifies atoms according to the fractional positions in the PAR file 
C *****						- abandons the "cryst_structure" TYPE
C *****						- distinguishes low-T and hi-T behaviour of cores and shells: 
C *****									low-T strongly bound, k_B*T/2 = E_kin = E_kin_C + E_kin_S
C *****									hi-T independent movement, k_B*T/2 = E_kin = E_kin_C = E_kin_S
C ***** Ver. 1.52 - adaptation to new .PAR file 
C *****						- atoms are identified according to their atom_name (2nd column in .PAR)
C ***** Ver. 1.53 - uses new semi-sequential binary file format 
C *****						- record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****           - CELL, QUICK and BULK data type
C *****
C ***** reads ASCII text output from LAMMPS(Sandia) MD simmulations (trajectory file)
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

      integer,parameter :: n_dim  =  3		      !number of position dimensions
      integer,parameter :: l_rec  =  1024		    !record length in real(4)

      real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   

			logical :: found
      character(4) :: at_name_inp,col(32)
      character(10) :: sim_style,subst_name,c_date,c_time,c_zone,ext,number,data_type,items(4)
      character(16) :: string,string2,section
      character(128) :: line,cwd_path,data_path,time_stamp
      character(128) :: file_dat,file_trajectory,file_master,file_inp,file_log
      
      integer ::  at_no,at_ind_in(3),at_ind_s(3),at_ind_base(3),i_dom,n_dom,i_dom_rec,n_at_cell
      integer ::  j_at,j_yes,j_shell,nskip,nfile_min,nfile_max,nfile_step,i_time(8)
      integer ::  ios,i,j,k,ii,jl,jat,jhead,j_step,j_label,n_label,n_site,j_eq,j_head_in,n_tstep0,n_tstep,n_items
      integer ::  i_rec,jrec,nrec,l_rec4,ifile,ncell,n_tot,nsuper,nrow,nlayer,ind,ind_rec,n_row_eff(3)
      integer ::  sc_c1,sc_c2,sc_m,j_data,j_struct,j_verb,j_shrec,j_ext,j_time,j_basis
      integer ::  ind_id,ind_type,ind_atom,ind_pos,ind_vel,ind_mass,ind_charge,ind_force,n_col			

      character(4),allocatable :: at_name(:),at_label(:)
      character(16),allocatable :: data_line(:)
      integer,allocatable ::  ind_l(:),i_site(:,:),ind_at(:)
      real,allocatable ::  e_kin(:,:),e_kin_s(:,:),x_pos(:,:),at_occup(:),at_mass_c(:),at_mass_s(:)


CC      real, allocatable :: eq_pos(:,:),at_base(:,:),conc(:)

      real :: at_mass,at_mass2,at_charge,at_charge2,item_value(4),sc_r
      real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force(3),at_force2(3)
      real ::	at_pos_0(3),at_pos_min(3),at_pos_max(3),a_cell(3,3),cell_par,atp,b_coh
      real :: t_step0,t0,t1,t2,dt,pos_inp(3)
      real :: temp,temp_r_c,temp_r_s,eps_x,x_lo,x_hi,y_lo,y_hi,z_lo,z_hi
      
C **** the following variables MUST have the following dimensions (4) or multiples because of alignement in the binary output file
C
      character(4),allocatable :: at_name_par(:),at_name_out(:)
      integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)
      integer,allocatable,target :: i_series(:),at_ind_out(:)
      integer,pointer ::  jr(:)

      real(4),allocatable ::	at_pos_c(:,:),at_veloc_c(:,:),at_force_c(:,:),at_occup_r(:)
      real(4),allocatable ::	at_pos_s(:,:),at_veloc_s(:,:),at_force_s(:,:)

      character(16)  :: sim_type,file_title
      integer(4)     :: n_row(3),n_atom,n_eq,j_force,j_shell_out,n_traj,n_cond,idum,n_rec
      real(4)        :: rec_zero(l_rec),t_ms,t_step,a_par(3),angle(3),temp_r

      data rec_zero/l_rec*.0/

C
C *****************************************************************************************
C
			write(*,*) '*** Program MP_LBIN 1.53 ** Copyright (C) Jiri Kulda (2019,2021,2022) ***'
      write(*,*)

C ********************* Get a time stamp and open a .LOG file *******************************
      call getcwd(cwd_path)
	    call date_and_time(c_date,c_time,c_zone,i_time)
	    write(time_stamp,110) c_date(1:4),c_date(5:6),c_date(7:8),c_time(1:2),c_time(3:4),c_time(5:6)
110		format(a,'/',a,'/',a,' ',a,':',a,':',a)	    
	    write(file_log,111) trim(c_date)
111   format('mp_tools_',a,'.log')
			inquire(file=file_log,exist=found)
			if(found)then
		    open (9,file=file_log,position='append')
		  else
		    open (9,file=file_log)
		  endif

			write(9,*) 
			write(9,*) trim(time_stamp),'  MP_LBIN 1.53  ',trim(cwd_path)
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
			read(4,*) 		j_shrec		! j_shrec=0: SHELLS data won't be recorded even if exist, (1 = SHELL data, if exist, will be recorded)
			read(4,*) 		!j_grid		! (0/1) grid overlay on PG_PLOT graphs 
			read(4,*) 		sim_type	! 'timestep' for MD, 'static' for static MD or DISCUS, must contain '_bulk' for BULK data type 
			rewind(4)

			data_type = 'cell'
			if(index(sim_type,'bulk')/=0) data_type = 'bulk'
			if(index(sim_type,'quick')/=0) then
				write(*,*) 'LAMMPS: the data_type QUICK is too dangerous, use CELL or BULK'
				stop
			endif
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

			read(4,*) subst_name
			read(4,*) n_atom,n_dom,		!number of atoms, displacement domain type (1=[100], 2=[110], 3= 111])
     1									eps_x 	!fractional coordinate tolerance range
			read(4,*) n_row						!nrow(3) supercell format
			allocate(at_name_par(n_atom),at_label(n_atom),x_pos(n_atom,3),e_kin(n_atom,3),e_kin_s(n_atom,3))
			allocate(ind_l(n_atom),i_site(n_atom,n_atom),nsuper_r(n_atom),at_occup(n_atom),at_occup_r(n_atom))
			allocate(at_mass_c(n_atom),at_mass_s(n_atom))

			idum = 0
			l_rec4 = l_rec/4
			ind_l = 0
			j_label = 0
			n_site = 0
			n_eq = 1								!later introduce sites, basis etc.
			j_force = 0
			nsuper = n_row(1)*n_row(2)*n_row(3)
			nlayer = n_row(1)*n_row(2)						!only to be used for record number calculation
			i_dom = 0
			n_cond = 0
			
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
					at_label(j_label) = string
					ind_l(i) = 1
					i_site(j_label,ind_l(i)) = j
				enddo
				n_label = j_label
				write(*,*) 'Atom labels ',n_label
				do j=1,n_label
					write(*,*) ind_l(j),at_label(j),(i_site(j,i),i=1,ind_l(j))
				enddo	
			
				write(*,*) trim(subst_name),' structure info (atoms): '	  
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

			read(4,*) !j_test		! number of snapshots to test the lattice parameter (1 minimum, 3 prudent; more superfluous) 
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

			rewind(4)
      section = 'mp_lammps'
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
          stop
        endif
        if(string(1:9).eq.section) exit	!find the mp_lammps part of the .par file
      enddo
			read(4,*) j_basis  	!atom type numbers (2nd column in data) corr. to  0 = chemical species, 1 = basis positions (PREFERABLE!
			read(4,*) t_ms			!MD microstep in time
			close(4)
			
			angle = 90.				!unit cell angles
			at_charge = .0
			temp = .0

C *** look for input file with MD trajectory
	  
      write(*,*) 'snapshot filename master (w/o .1000): '
      read(*,*) file_master
      
      write(*,*) 'snapshots to read: nfile_min, nfile_step, nfile_max'
      read(*,*) nfile_min, nfile_step,nfile_max 
      write(9,*) 'Read snapshots number: min, step, max',nfile_min, nfile_step,nfile_max
      write(*,*) 'Saved snapshot numbers start:',nfile_min/nfile_step
      write(9,*) 'Saved snapshot numbers start:',nfile_min/nfile_step
      write(9,*)
            
			write(string,*) nfile_min
      file_inp = trim(data_path)//trim(adjustl(file_master))//'.'//trim(adjustl(string))
      open (1,file=file_inp,action='read',status ='old')

      write(9,*) '1st snapshot file to read:',file_inp

      j_shell = 0
      jhead = 0
      nskip = 20 ! number of records to test
			j_time = 0				!1 indicates presence of TIME
      n_items = 0
      
C ***  skip the first <nskip> lines until 'TIMESTEP' !normally there are none     
			do i=1,nskip
				read(1,99) line
99      format(a80)
				if(index(line,'TIMESTEP')/=0) then
					read(1,*) n_tstep0
					n_items = n_items+1
					exit
				endif
			enddo
			rewind(1)

			do i=1,nskip
				read(1,99) line
				if(index(line,'TIME')/=0.and.index(line,'TIMESTEP')==0) then
					read(1,*) t_step0
					n_items = n_items+1
					j_time = 1
					exit
				else
					write(*,*)'Input data do not contain explicit TIME information!'
					write(*,*)'The frequency scale will rely on T_MS from the .PAR file!'
					write(*,*)'Continue? (1/0)'
					read(*,*) ii
					if(ii==0) stop
					write(9,*)'Input data do not contain explicit TIME information!'
					write(9,*)'The frequency scale will rely on the .PAR file T_MS=',t_ms
					exit
				endif
			enddo
			rewind(1)
				
			do i=1,nskip
				read(1,99) line
				if(index(line,'NUMBER')/=0) then
					read(1,*) n_tot
					n_items = n_items+1
					exit
				endif
			enddo
			rewind(1)

			do i=1,nskip
				read(1,99) line
				if(index(line,'BOUNDS')/=0) then
					read(1,*) x_lo,x_hi
					read(1,*) y_lo,y_hi
					read(1,*) z_lo,z_hi
					n_items = n_items+1
					exit
				endif
			enddo
			rewind(1)

			if(n_items<3) then
				write(*,*) 'Snapshot header info not complete:'
				write(*,*) 't_step,n_tstep,n_tot,x_lo,x_hi',t_step,n_tstep,n_tot,x_lo,x_hi
				stop
			endif		 
				
			      do i=1,nskip
        read(1,'(a)') line
CC99      format(a80)
        if(index(line,'ITEM: ATOMS').gt.0) exit
        if(i.eq.nskip) then
          write(*,*) 'ATOMS: data not found'
          stop
        endif
      enddo 
      
C *** parser of the data record description
C
			n_col = 1
			do 
				read(line(index(line,'S')+1:),*,iostat=ios) (col(i),i=1,n_col)
				if(ios/=0) exit
				n_col = n_col+1
			enddo      
      n_col = n_col-1
      
CC      write(*,*) 'line',line
CC      write(*,*) 'n_col,col', n_col,(col(i),i=1,n_col)

C *** typical case:  id type element mass x y z vx vy vz

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
      
      allocate(data_line(n_col))
      
C *** look for possible core/shell presence, assuming all atoms would be the same
C  typical data:
C  1 1 Na1_c 20.9 0.167166 0.254627 0.115284 3.148 -3.06563 0.629693 
C  2 9 Na1_s 2.09 0.165161 0.250872 0.113053 2.9716 -2.99374 0.691807 
  	
CC    	read(1,*,iostat=ios) ii,jat,string				
CC    	read(1,*,iostat=ios) ii,jat,string2

			read(1,*)  data_line  
   		string = 	data_line(ind_atom)
			read(1,*)  data_line  
 			string2 = 	data_line(ind_atom)

    	if((index(string,'_c').gt.1.or.index(string,'_C').gt.1).and.
     1    	(index(string2,'_s').gt.1.or.index(string2,'_S').gt.1))	j_shell = 1		
      close(1)
      
      if(j_shell==1) then
       	write(*,*)'Core & shell data found'
				if(j_shrec==0) write(*,*)'Shell data NOT to be recorded (change this in .PAR)'
				if(j_shrec==1) write(*,*)'Shell data WILL be recorded'
				if(j_shrec==1) write(9,*)'Shell data WILL be recorded'
      endif
      j_shell_out = j_shell*j_shrec

        
C *** initialise the process from data in the first frame
			if(j_shell.eq.1) n_tot = n_tot/2		!shell model counts both cores and shells as atoms
CC			write(*,*) 'n_tot',n_tot
			nsuper = n_tot/(n_atom*n_eq)
      nrow = anint(real(nsuper)**(1./3.))
      nlayer = nrow**2  ! 48 gives   2304
			n_row = nrow					!so far assuming cubic supercells

      at_ind_base = 1  
      
      a_par(1) = (x_hi-x_lo)/nrow							!assuming LAMMPS system created by LATTICE command
      a_par(2) = (y_hi-y_lo)/nrow
      a_par(3) = (z_hi-z_lo)/nrow
      
      at_pos_0 = (/x_lo,y_lo,z_lo/)

			if(j_shell.eq.0) then
				sim_style = 'atom'
			elseif(j_shell.eq.1) then
				sim_style = 'core/shell'
			else
				sim_style = 'undefined!'
			endif
			
CC			sim_type = 'timestep'
			write(9,*) 'sim_type, sim_style = ',sim_type,sim_style
			if(j_verb==1) then
	      write(*,*) 'sim_type, sim_style = ',sim_type,sim_style
CC      write(*,*) 'trajectory recording mode =', n_traj
CC      write(*,*) 'boundary conditions ', n_cond
CC      write(*,*) 't_step,t_ms:',t_step,t_ms
				write(*,*) 'nominal temperature [K] ',temp
				write(*,*) 'substance:  ',subst_name,'n_atom, n_eq, nrow:',n_atom,n_eq,nrow
			endif
			write(*,*) 'a_par:',a_par,' cubic supercell assumed!'
      write(9,*) 'nominal temperature [K] ',temp
      write(9,*) 'substance:',subst_name,'n_atom, n_eq,nrow:',n_atom,n_eq,nrow
      write(9,*) 'a_par:',a_par
CC      write(*,*) 'fractional position tolerance range:',eps_x
CC      write(*,*) 'domain analysis:',n_dom,' 0=no, 1=[100], 2=[110], 3=[111]'
      write(*,*) 'go ahead? (1/0)'
      read(*,*) j_yes
		  if (j_yes.eq.0) stop

      write(9,*)			
      write(9,*) 'Trajectory time start, step [ps]:',t_step0,t_ms

			allocate(i_series(n_tot))
			i_series = (/ (i, i = 1, n_tot) /)

C **** the file_loop 

			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
			sc_r = 1./sc_r
			call cpu_time(t0)

      file_loop: do ifile=nfile_min,nfile_max,nfile_step
				allocate(at_pos_c(4,n_tot),at_veloc_c(4,n_tot),at_ind(4,n_tot))
				if(j_shell_out==1) allocate(at_pos_s(4,n_tot),at_veloc_s(4,n_tot))
				if(j_force==1) allocate(at_force_c(4,n_tot),source=.0)
				if(j_force==1.and.j_shell_out==1) allocate(at_force_s(4,n_tot),source=.0)

				e_kin = .0
				e_kin_s = .0
				 
				nsuper_r = 0
				j_head_in = 1

				write(string,*) ifile
      	file_inp = trim(data_path)//trim(adjustl(file_master))//'.'//trim(adjustl(string))
				open (1,file=file_inp,action='read',status ='old',iostat=ios)
				if(ifile==nfile_min.or.ifile==10*nfile_step*(ifile/(10*nfile_step))) write(*,*)trim(file_inp),'   '
				if(ifile==nfile_min) write(*,*) '(this may take time)'
				if(ifile==nfile_max) write(9,*) 'last snapshot file to read:',file_inp
			     
C ***  start reading the text and recording the binary output

				do i=1,nskip
					read(1,99) line
					if(j_time==0)then
						if(index(line,'TIMESTEP')/=0) then 
							read(1,*) n_tstep
							t_step = n_tstep*t_ms
							exit
						endif
					else
						if(index(line,'TIME')/=0.and.index(line,'TIMESTEP')==0) then
							read(1,*) t_step
							exit
						endif
					endif
				enddo 

C *** check the t_ms value

				if(j_time==1.and.ifile==nfile_min+nfile_step) then
					dt = (t_step-t_step0)/real(nfile_step)
					if(dt/=t_ms) then
						write(*,*) 'Correcting t_ms(PAR)',t_ms,' to',dt,' from the snapshot data'
						write(9,*) 'Correcting t_ms(PAR)',t_ms,' to',dt,' from the snapshot data'
						t_ms = dt
					endif					
				endif

			
				do i=1,nskip
					read(1,99) line
					jhead = index(line,'ITEM: ATOMS')
					if(jhead.ne.0) exit
				enddo 

				at_pos_min = .0
				at_pos_max = .0							

        read_loop: do i = 1,n_tot				!LAMMPS dump is potentially disordered, but the true N_TOT is there

					read(1,*)  data_line   		
          if(ios.ne.0.and.i.lt.n_tot) then
						write(*,*) 'EOF or EOR at',ii,jat,string 
						stop
          endif
c   				write(*,*) data_line
   				read(data_line(ind_id),*) jrec
   				if(ind_type/=0)read(data_line(ind_type),*) jl
					
   				string = data_line(ind_atom)
   				if(ind_mass/=0)read(data_line(ind_mass),*) at_mass
   				if(ind_charge/=0)read(data_line(ind_charge),*) at_charge
   				read(data_line(ind_pos:ind_pos+2),*) at_pos_in
   				if(ind_vel/=0)read(data_line(ind_vel:ind_vel+2),*) at_veloc_in
   				if(ind_force/=0)read(data_line(ind_force:ind_force+2),*) at_force
   				
C *** read also the SHELL line if needed  
					if(j_shell.eq.1) then
						read(1,*)  data_line   		
						if(ios.ne.0.and.i.lt.n_tot) then
							write(*,*) 'EOF or EOR at',ii,jat,string 
							stop
						endif
c   				write(*,*) data_line
						string2 = data_line(ind_atom)
						if(ind_mass/=0)read(data_line(ind_mass),*) at_mass2
						if(ind_charge/=0)read(data_line(ind_charge),*) at_charge2
						read(data_line(ind_pos:ind_pos+2),*) at_pos_in2
						if(ind_vel/=0)read(data_line(ind_vel:ind_vel+2),*) at_veloc_in2
						if(ind_force/=0)read(data_line(ind_force:ind_force+2),*) at_force2
					
						if(string2(1:j_eq-1)/=string(1:j_eq-1).or.jl-n_atom*n_eq/=jat) then
							write(*,*) 'inconsistent core/shell data',ii,jat,jl,string,string2
							write(9,*) 'inconsistent core/shell data',ii,jat,jl,string,string2,  'STOP!'
							stop
						endif
					endif
  
          if(j_shell.eq.1) then
          	ii = index(string,'_c')						!supposed core comes first
          	at_name_inp = string(1:ii-1)
          else
          	at_name_inp = trim(string)
          endif 

   				at_pos_in = (at_pos_in-at_pos_0)/a_par
  				if(j_shell_out.eq.1) at_pos_in2 = (at_pos_in2-at_pos_0)/a_par

					if(data_type=='bulk') then
						jat = jl				
						do k=1,3
							if(at_pos_in(k)<at_pos_min(k)) at_pos_min(k) = at_pos_in(k)
							if(at_pos_in(k)>at_pos_max(k)) at_pos_max(k) = at_pos_in(k)
							if(at_pos_in(k)<.0) at_pos_in(k) = at_pos_in(k)+n_row(k)					!bring the position into the box
							if(at_pos_in(k)>real(n_row(k))) at_pos_in(k) = at_pos_in(k)-n_row(k)							
						enddo
						jrec = i
						at_pos_c(1:3,jrec) = at_pos_in-n_row/2
						at_pos_c(4,jrec) = .0							!could be atom charge if individualised

						if(j_shell_out==1) then
							do k=1,3
								if(at_pos_in2(k)<0.) at_pos_in2(k) = at_pos_in2(k)+n_row(k)					!bring the position into the box
								if(at_pos_in2(k)>real(n_row(k))) at_pos_in2(k) = at_pos_in2(k)-n_row(k)							
							enddo
							at_pos_s(1:3,jrec) = at_pos_in2-n_row/2
							at_pos_s(4,jrec) = .0							!could be atom charge if individualised
						endif
						
						at_ind(1,jrec) = jrec
						at_ind(2:3,jrec) = 0
						at_ind(4,jrec) = jat

					else			!data_type=='cell'
						if(j_basis==1) then
							jat = jl
						elseif(at_name_inp.eq.trim(at_label(jl))) then		
							if(ind_l(jl).eq.1) then
								jat = i_site(jl,1)	!atom name and fractional position is identified
							else
								do ii=1,ind_l(jl)		!multiple basis positions correspond to the at_label
									pos_inp = at_pos_in-x_pos(i_site(jl,ii),:)
									jat = i_site(jl,ii)
									if(maxval(abs(pos_inp-anint(pos_inp))).lt.eps_x) exit !atom found
								enddo
							endif						
						else																					!if still not found GIVE UP!
							write(*,*) 'Atom names do not match:',jrec,jl,at_pos_in,at_name_inp,'   ',at_label(jl)
							write(*,*) 'Check your data and the MP_LAMMPS settings in your .PAR file (j_basis)!'
							stop
						endif

						at_ind_in = anint(at_pos_in-x_pos(jat,:))+at_ind_base		!they will serve as pointers to the right order of atom records
						at_pos_in = at_pos_in-n_row/2						!now the supercell will be centred as from DL_POLY
  					if(j_shell_out.eq.1) at_pos_in2 = at_pos_in2-n_row/2
					
C *** check supercell boundary wrapping
						do k=1,3
							if(at_ind_in(k)==0) then
								at_ind_in(k)=at_ind_in(k)+n_row(k)
								at_pos_in(k)=at_pos_in(k)+n_row(k)
 								if(j_shell_out.eq.1) at_pos_in2(k) = at_pos_in2(k)+n_row(k)
							endif
							if(at_ind_in(k)==nrow+1) then
								at_ind_in(k)=at_ind_in(k)-n_row(k)
								at_pos_in(k)=at_pos_in(k)-n_row(k)
 								if(j_shell_out.eq.1) at_pos_in2(k) = at_pos_in2(k)-n_row(k)
							endif
						enddo
 
						jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1)
						at_ind(1:3,jrec) = at_ind_in
						at_ind(4,jrec) = jat
					endif		!data_type
          
 					at_pos_c(1:3,jrec) = at_pos_in
 					at_pos_c(4,jrec) = at_charge							
          at_veloc_c(1:3,jrec) = at_veloc_in
          at_veloc_c(4,jrec) = at_mass	
          if(j_force==1) at_force_c(1:3,jrec) = at_force

 					if(j_shell_out.eq.1) then
						at_pos_s(1:3,jrec) = at_pos_in2
						at_pos_s(4,jrec) = at_charge2							
						at_veloc_s(1:3,jrec) = at_veloc_in2
						at_veloc_s(4,jrec) = at_mass2	
						if(j_force==1) at_force_s(1:3,jrec) = at_force2
 					endif

C *** accumulate the occupation number and the kinetic energy to refine the real temperature
          nsuper_r(jat) = nsuper_r(jat)+1
					do k = 1,3
						e_kin(jat,k) = e_kin(jat,k) + at_mass*at_veloc_in(k)**2
						if(j_shell==1) e_kin_s(jat,k) = e_kin_s(jat,k) + at_mass2*at_veloc_in2(k)**2
					enddo
	      enddo read_loop

C *** for 'BULK' centre the atom cloud, produce virtual n_row and order the AT_IND and AT_POS arrays by atoms
				allocate(at_ind_out(n_tot),at_name_out(n_atom),ind_at(n_atom))


				if(data_type=='bulk') then
CC					do i=1,n_tot
CC						at_pos_c(1:3,i) =	at_pos_c(1:3,i) - .5*(at_pos_max-at_pos_min)
CC					enddo
CC
CC        	if(j_shell_out.eq.1) then
CC						do i=1,n_tot
CC							at_pos_s(1:3,i) =	at_pos_s(1:3,i) - .5*(at_pos_max-at_pos_min)
CC						enddo
CC        	endif

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
				if(ifile==nfile_min) write(*,*) '        total of',n_tot,' atoms read in',t2-t0,' sec'


C *** get the atom position occupation numbers and compare them with those from the .par file
				if(ifile.eq.nfile_min) then
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

C *** normalize the kinetic energy
        do j=1,n_atom*n_eq
	        e_kin(j,:) = e_kin(j,:)/nsuper
	        if(j_shell.eq.1) e_kin_s(j,:) = e_kin_s(j,:)/nsuper
        enddo
	
C *** retrieve masses for display
				do j=1,n_atom
					do i=1,nlayer
						if(at_veloc_c(4,(j-1)*nsuper+i)/=0.)exit
					enddo
					at_mass_c(j) = at_veloc_c(4,(j-1)*nsuper+i)
					if(j_shell.eq.1) at_mass_s(j) = at_veloc_s(4,(j-1)*nsuper+i)
				enddo
			
C *** get the true temperature
				if(j_verb==1.and.ifile==nfile_min) then
					write(*,*)
					write(*,*) 'nominal k_B*T/2 [DAPS]',.5*k_B*temp
					write(*,*) 'at_mass_c:',(at_mass_c(jat),jat=1,n_atom*n_eq)
					write(*,*) 'e_kin_c(jat,:)',(e_kin(jat,:),jat=1,n_atom*n_eq)
				endif
				temp_r = sum(e_kin(1:n_atom*n_eq,:))/(n_atom*n_eq*n_dim*k_B) !the true temperature

				if(j_shell.eq.1) then
					if(j_verb==1.and.ifile==nfile_min) then
						write(*,*) 'at_mass_s:',(at_mass_s(jat),jat=1,n_atom)
						write(*,*) 'e_kin_s(jat,:)',(e_kin_s(jat,:),jat=1,n_atom)
					endif
				
					temp_r_s = sum(e_kin_s(1:n_atom*n_eq,:))/(n_atom*n_eq*n_dim*k_B)  !the true temperature
					temp_r_c = sum(e_kin(1:n_atom*n_eq,:))/(n_atom*n_eq*n_dim*k_B) !the true temperature

					if (abs(temp_r_c-temp_r_s).le..1*temp_r_c) then
						temp_r = (temp_r_c+temp_r_s)*.5						!hi-T limit: independent C and S vibrations
					else
						temp_r = temp_r_c+temp_r_s	!low-T limit: strongly bound C and S vibrations
					endif
				endif

C **** use the new binary format *** do all done by array(1:4,i) to be compatible with the record length ***
C *** generate output filename

				write(file_dat,103) trim(file_master),ifile/nfile_step
103     format('./data/',a,'_n',i4.4,'.dat')
CC				if(ifile==nfile_min.or.ifile==10*nfile_step*(ifile/(10*nfile_step))) write(*,*)trim(file_dat)

				if(j_verb==1) then
					if(j_shell==0) then
						write(*,*) trim(file_dat),'  T_tot = ',temp_r
						write(9,*) trim(file_dat),'  T_tot = ',temp_r
					elseif(j_shell==1) then
						write(*,*) trim(file_dat),'  T_core/shell/total = ',temp_r_c,temp_r_s,temp_r
						write(9,*) trim(file_dat),'  T_core/shell/total = ',temp_r_c,temp_r_s,temp_r
					endif
				endif
        

				open(2,file=file_dat,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4

				n_rec = (n_tot/l_rec4)														!for each position there are 4 components
				if(mod(n_tot,l_rec4)/=0) n_rec = n_rec+1																
CC				if(ifile==nfile_min) write(*,*) 'n_tot,l_rec4,n_rec',n_tot,l_rec,n_rec			
				
C *** write the header record
				i_rec = 1
				write(2,rec=i_rec)rec_zero				!fill zeros first
				write(2,rec=i_rec) 
     1		sim_type,file_title,t_ms,t_step,temp_r,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,
     2    n_cond,n_rec,n_tot,at_name_par,at_occup_r(1:n_atom),nsuper_r(1:n_atom)																		!char(16): sim_type,file_title, char(4): at_name

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
				close(1)
				deallocate(at_pos_c,at_veloc_c,at_ind,ind_at,at_ind_out,at_name_out)
				if(j_shell_out==1) deallocate(at_pos_s,at_veloc_s)
				if(j_force==1) deallocate(at_force_c)
				if(j_force==1.and.j_shell_out==1) deallocate(at_force_s)

      enddo file_loop 

			CALL SYSTEM_CLOCK (COUNT = sc_c2)
CC			write(*,*) nfile,' files read in ', t2-t1,' SYS time',(sc_c2-sc_c1)*sc_r

			ifile = 1+(nfile_max-nfile_min)/nfile_step
      write(*,*) 'Trajectory finished: ',ifile,' .dat files written in SYS time',(sc_c2-sc_c1)*sc_r,' sec'
      write(9,*) 'Trajectory finished: ',ifile,' .dat files written in SYS time',(sc_c2-sc_c1)*sc_r,' sec'
      stop

      end program mp_lbin53
      