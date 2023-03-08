      
      program mp_dbin55

C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		  program MP_BIN 1.55   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
C *****   %%%%%%%%%%%%%%%%   			program MP_BIN 1.55  				 %%%%%%%%%%%%%%%%%%%%%%%%
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
C *****           - CELL, FAST and BULK input methods
C *****           - BULK input uses simulation box sizes
C *****           - j_force replaced by n_traj
C *****           - j_centred in .PAR indicates coordinate centring in box (1/0)
C *****
C *****   FROZEN on 05/09/2022 10:45 and forked as mp_bin54
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

      integer ::  i_dom,i_dom_rec,n_dom,n_tot_in,j_struct,n_tot,i_save
      integer ::  at_no,at_ind_base(3),at_ind_shift(3),at_ind_in(3),at_ind_in2(3)
      integer ::  j_yes,nskip,nfile_min,nfile_max,nfile_step,i_time(8),n_save_min,izero,indzero(3)
      integer ::  ios,ios_t,i,j,k,ii,i2,i3,jl,jat,j_step,j_label,n_label,j_first,j_read,j_verb,j_proc,j_shrec,j_test
      integer ::  inrec,jrec,nrec,i_rec,l_rec4,ifile,ncell,nsuper,nrow,nlayer,n_site,j_shell
      integer ::  sc_c1,sc_c2,sc_m,nt_min,nt_max,nt_step,i_traj,j_mult,j_basis,j_centred

      real :: at_mass_in,at_mass_in2,at_charge_in,at_charge_in2,at_displ_in,sc_r
      real ::	at_pos_in(3),at_pos_in2(3),at_veloc_in(3),at_veloc_in2(3),at_force_in(3),at_force_in2(3),at_base_shift(3)
      real ::	dummy,at_pos2(3),at_pos3(3),at_veloc2(3),a_cell(3,3),a_cell_par(3),a_cell_half(3),at_pos_centre(3)
      real :: t1,t2,filter_fwhm,t_step,zero,at_zero(3),pos_inp(3),temp_par,eps_x,temp_r_s,temp_r_c

      character(4),allocatable :: at_name(:),at_label(:)
      integer,allocatable ::  ind_l(:),i_site(:,:),ind_at(:)
      integer,allocatable,target :: i_series(:),at_ind_out(:)
      integer,pointer ::  jr(:)
      real,allocatable ::  sum_pos(:,:,:),ord(:),e_kin(:,:),e_kin_s(:,:),at_base_in(:,:),at_base(:,:),at_occup(:)

C **** the following variables MUST have the following 32bit sizes or multiples because of alignement in the binary output file
C
      character(4),allocatable :: at_name_par(:),at_name_out(:),version
      integer(4),allocatable   :: at_ind(:,:),nsuper_r(:)

      real(4),allocatable ::	at_pos_c(:,:),at_veloc_c(:,:),at_force_c(:,:),at_occup_r(:)
      real(4),allocatable ::	at_pos_s(:,:),at_veloc_s(:,:),at_force_s(:,:)

      character(16)  :: sim_type,input_method,dat_type,dat_source,file_par
      integer(4)     :: n_row(3),n_atom,n_eq,j_shell_out,n_traj,n_cond,idum,n_rec,n_head,n_head_in1,n_head_in2
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
			write(*,*) '*** Program MP_BIN 1.55 ** Copyright (C) Jiri Kulda (2019,2021,2022) ***'
      write(*,*)
      dat_source = 'MP_TOOLS'
      version = '1.55'

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
			write(9,*) trim(time_stamp),'  MP_BIN 1.55  ',trim(cwd_path)
C		  write(*,*) trim(time_stamp),'  MP_BIN 1.55  ',trim(cwd_path)
			write(9,*) 

C *** diverse initialisations
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
			n_cond = 0
			n_tot_in = 0
      t_ms = .0
      temp_par = .0
      at_base_shift = .0
			filter_name = 'nn'
			filter_fwhm = .0

      
C *** read auxiliary file <file_par.par> with structure parameters, atom names and further info
      write(*,*) 'Parameter file name (.par will be added)'
      read(*,*) file_par
      file_inp = trim(file_par)//'.par'

			open(4,file=file_inp,action='read',status ='old',iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_inp),' not found! Stop execution.'
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
      
      if(ext=='ext'.or.ext=='EXT') ext=''
      if(ext/=''.and.index(ext,'.')==0) ext='.'//ext

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

      write(*,*) 
      write(*,*) 'Substance name: ', subst_name
			write(*,*) 'Sim_type, dat_type, input method: ',sim_type,dat_type,input_method		
			
			allocate(at_name_par(n_atom),at_label(n_atom),at_base_in(n_atom,3),at_base(n_atom,3))   !at_base would include at_base_shift & saved in data file header
			allocate(ind_l(n_atom),i_site(n_atom,n_atom),at_occup(n_atom))

			nsuper = n_row(1)*n_row(2)*n_row(3)
			nlayer = n_row(1)*n_row(2)						!to be used for record number calculation for CELL data
			
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
      write(*,*) 'Atom labels ',n_label
      do j=1,n_label
        write(*,*) ind_l(j),at_label(j),(i_site(j,i),i=1,ind_l(j))
      enddo	
    
      write(*,*) trim(subst_name),' structure info (atoms): '	  
      do j=1,n_atom
        if(input_method=='CELL'.or.input_method=='FAST') then
          write(*,*) j,at_name_par(j),at_base_in(j,:)
        else
          write(*,*) j,at_name_par(j)
        endif
      enddo
CC			write(*,*) 'dat_type:  ',dat_type
CC			write(*,*) 'ext:  ',ext
			
			angle = 90.				!assuming orthogonal lattice

      write(*,*) 'Read snapshots number: min, max'
      read(*,*) nfile_min, nfile_max 
			nfile_step = 1

      write(9,*) 'Read snapshots number: min, step, max',nfile_min, nfile_step,nfile_max
      write(*,*) 'Saved snapshot numbers start:'
      read(*,*) n_save_min 
      write(9,*) 'Saved snapshot numbers start:',n_save_min   
      
C *** input cycle over MD trajectory files
			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)			
			i_save = n_save_min
			i_traj = 0
	  
			write(*,*) 'Master for MD trajectory filename: '
			read(*,*) file_master

CC      i = len_trim(file_master)
CC      if(file_master(i:i)=='_') then
CC				file_master_out = file_master(1:i-1)
CC      else
CC				file_master_out = file_master(1:i)
CC      endif

			write(*,*) 'Read input files number (0 0 no numbers): min,max'
			read(*,*) nt_min,nt_max
			nt_step = 1
			t_single = (nt_min==0.and.nt_max==0)

CC *** tread 1st frame of 1st history file to get info
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

      write(*,*) 'Input trajectory file:  ',trim(file_inp)

      open (1,file=file_inp,action='read',status ='old',iostat=ios)
      if(ios.ne.0) then
        write(*,*) "Can't open the file ",trim(file_inp),'!'
        stop
      endif

      write(9,*) 'Reading MD trajectory file:  ',trim(file_inp)
					      
C ***  skip the first <nskip> lines until 'timestep' !normally they are two
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

      read(1,*) string,j_step,n_tot_in,n_traj,n_cond,t_ms,t_dump		!t_ms is MD microstep, t_dump is the snapshot time
																																			!n_tot_in total number of atoms in cells n_atom*nrow**3
			do j=1,3
				read(1,*) a_cell(j,:)
			enddo
			
			do j=1,3
			 a_par(j) = a_cell(j,j)/n_row(j)
			enddo

C *** read 1st atom record         
			read(1,*) at_name_in,at_no,at_mass_in,at_charge_in,at_displ_in
				read(1,*) at_pos_in
				if(n_traj>=1) read(1,*) at_veloc_in
				if(n_traj==2) read(1,*) at_force_in

C *** read 2nd atom record         
			read(1,*) at_name_in2,at_no,at_mass_in,at_charge_in,at_displ_in
				read(1,*) at_pos2
				if(n_traj>=1) read(1,*) at_veloc2
				if(n_traj==2) read(1,*) at_force_in

C *** read 3rd atom record         
			read(1,*) at_name_in3,at_no,at_mass_in,at_charge_in,at_displ_in
				read(1,*) at_pos3
				if(n_traj>=1) read(1,*) at_veloc_in
				if(n_traj==2) read(1,*) at_force_in
			close(1)

C *** analyse the input
			if(at_name_in.eq.at_name_in2) then
				j_shell = 0				
				write(*,*) 'No shells'               !', supercell is:',nrow,'^3',' j_shell =,',j_shell
			else if(trim(at_name_in)//'s'.eq.at_name_in2.and.at_name_in.eq.at_name_in3) then
				j_shell = 1
				n_tot_in = n_tot_in/2
				write(*,*) 'Found a shell candidate   ',at_name_in2		!,' j_shell =',j_shell    !', supercell is:',nrow,'^3'
			else
				write(*,*) 'Strange primary data, please check them!'
				write(*,*) at_name_in,trim(at_name_in)//'s'
				write(*,*) at_name_in2
				write(*,*) at_name_in3
				stop
			endif
			
			if(n_traj==0) then
				write(*,*) 'Type in nominal temperature [K]:'
				read(*,*) temp_par
			endif
			
			if(n_traj>=1.and.at_veloc_in(1)==at_veloc2(1).and.at_veloc_in(2)==at_veloc2(2).and.at_veloc_in(3)==at_veloc2(3)) then
				write(*,*) 'Strange velocities:'
				write(*,*) at_veloc_in   
				write(*,*) at_veloc2          
				write(*,*) 'Type in nominal temperature [K]:'
				read(*,*) temp_par
			endif
			
      if(j_shell==1) then
       	write(*,*)'Core & shell data found'
				if(j_shrec==0) write(*,*)'Shell data NOT to be recorded (change this in .PAR)'
				if(j_shrec==1) write(*,*)'Shell data WILL be recorded'
				if(j_shrec==1) write(9,*)'Shell data WILL be recorded'
      endif
      j_shell_out = j_shell*j_shrec

			write(*,*) 'Simulation type = ',sim_type
			write(*,*) 'Trajectory recording mode =', n_traj
			write(*,*) 'Boundary conditions ', n_cond
      write(*,*) 'Trajectory time start [ps]:',t_dump
      if(n_traj==0) write(9,*) 'Using nominal temperature [K] ',temp_par

      write(9,*) 'Simulation type = ',sim_type
      write(9,*) 'Trajectory recording mode =', n_traj
      write(9,*) 'Boundary conditions ', n_cond
      write(9,*) 'Trajectory time start, step [ps]:',t_dump,t_ms
      if(n_traj==0) write(9,*) 'Using nominal temperature [K] ',temp_par

C *** handle the supercell parameters and the origin of the supercell coordinate system
      if(input_method=='CELL'.or.input_method=='FAST') then
				if(nsuper==1) then
					write(*,*) 'CELL and FAST methods not compatible with N_ROW = 1'
					stop
				endif
				n_tot = n_atom*nsuper
			else						!input_method=='BULK'
				n_tot = n_tot_in
			endif
			
			if(nsuper/=1)then
				if(j_centred==1) then
						at_ind_shift = 0
				else
						at_ind_shift = n_row/2
				endif
				at_ind_base = n_row/2+1-at_ind_shift
			endif

      allocate (at_ind(4,n_tot),SOURCE=0)
			allocate(at_name(n_tot),SOURCE='    ')
			allocate(e_kin(n_atom,3),at_occup_r(n_atom))
			allocate(nsuper_r(n_atom),SOURCE=0) !atom number jat is the first of the four indices
			if(j_shell.eq.1)allocate(e_kin_s(n_atom,3),SOURCE=0.0)
			allocate(i_series(n_tot))
			i_series = (/ (i, i = 1, n_tot) /)

C *** Now ready to cycle over trajectory files, each snapshot to be saved in a separate binary file

			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)

			nsuper_r = 0
			ifile = 1

			trajectory_loop: do i_traj=nt_min,nt_max
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
				
			     				
C *** cycle over snapshots, each to be saved in a separate binary file
			j_read = 0

      frame_loop: do 				!we have to go through all the snapshots in the .TXT input

        do 						
          read(1,*,iostat=ios) line
          if(ios<0) then   !end of file
CC						write(*,*) 'End of trajectory file: ',ios
						exit frame_loop 
					endif
CC          if(ios>0) cycle     !corrupt input data, go on
          if(ios>0) then
          	write(*,*) 'Line after:',line
          	write(*,*) 'Input error, corrupt trajectory data?'
          	stop
          endif
          if(index(line,head).ne.0) exit
        enddo 
        
 
C *** skip the snapshots that are not requested to be read
        if((ifile.lt.nfile_min).or.(mod(ifile-nfile_min,nfile_step).ne.0)) then  !skip the snapshot
        	ifile = ifile+1
          cycle frame_loop         ! cycle the frame_loop
        else
	        backspace(1)	
	        j_read = j_read+1
        endif                    
 
C *** no need to re-read the first header line but
C *** take fresh lattice parameters for each snapshot - they may evolve
        read(1,*) string,j_step,n_tot_in,n_traj,n_cond,t_ms,t_dump
        do j=1,3
          read(1,*) a_cell(j,:)
          if(input_method=='BULK') then        
            a_cell_half(j) = .5*a_cell(j,j)
            if(j_centred==0) then
              at_pos_centre(j) = a_cell_half(j)
            else
              at_pos_centre(j) = .0
            endif
		 			endif
  		 		a_par(j) = a_cell(j,j)/n_row(j)
        enddo 
            
				if(i_traj==nt_min.and.ifile==nfile_min) then
				  if(nsuper/=1) then
            write(*,*) 'Lattice parameter estimate =',a_par
            write(*,*) 'OK? (1/0)'
            read(*,*)	j_yes

            if(j_yes.ne.1) then
              write(*,*) 'input a better A_PAR(1:3)'
              read(*,*)  a_par
            endif
				  else
					  write(*,*) 'Simulation box size =',a_par
					endif			
				endif

C *** get the actual time step of the sequence
				if(ifile==nfile_min) then
				  t_step = t_dump
				elseif(ifile==nfile_min+1) then
				  t_step = t_dump-t_step
				  t_dump = t_dump-t_step !get the old t_dump for the old frame          
C *** correct t_step in the first snapshot				  
				  open(2,file=file_dat,access='direct',form='unformatted',recl=4*l_rec)		! l_rec is in 32 bit words = 4 bytes, thus the factor 4
          write(header_record,nml=data_header_1)	
          write(2,rec=2) header_record
          close(2)
				  t_dump = t_dump+t_step !get the right t_dump for the present frame          
				endif
				
C ***  read-in the text of a snapshot in one go
				if(ifile==nfile_min) write(*,*) 'reading the 1st snapshot (takes a few seconds) ...'

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
				
        call cpu_time(t1)				

        read_loop: do  inrec=1,n_tot_in        !swallow the snapshot

C *** first read the CORE data  

          read(1,*,iostat=ios_t) at_name_in,at_no,at_mass_in,at_charge_in,at_displ_in
CC						if(ios/=0) write(*,*) 'Error on ASCII input IOS: ',ios
CC						if(ios>0) cycle frame_loop     !corrupt input data, jump out to new snapshot

          if(ios_t<0.or.at_name_in.eq.head) then
CC          	write(*,*) 'End of snapshot'
CC          	backspace(1)
          	exit read_loop    !end of snapshot, end of file
          endif

          read(1,*) at_pos_in
          if(n_traj>=1) read(1,*,iostat=ios_t) at_veloc_in
          if(ios_t/=0) then
          	write(*,*) 'Input problem: ios_t,inrec,ifile,at_name_in,head',ios_t,inrec,ifile,at_name_in,head
					endif
          if(n_traj==2) read(1,*) at_force_in

C
C *** now read the SHELL data if needed
          if(j_shell.eq.1) then
          	read(1,*) at_name_in2,at_no,at_mass_in2,at_charge_in2,at_displ_in
						if(trim(at_name_in)//'s'.ne.at_name_in2) then
							write(*,*)'wrong core-shell sequence',inrec,at_name(inrec),at_name_in
							stop
						endif
						read(1,*) at_pos_in2
						if(n_traj>=1) read(1,*) at_veloc_in2
						if(n_traj==2) read(1,*) at_force_in2
						at_no = at_no/2
					endif		!j_shell

C ***  treat the CORE data and get the right labels & positions 
					jl = 0				
					if(j_basis==0) then
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
          else
             do j=1,n_atom
              if (at_name_in.eq.at_name_par(j)) then
                jl=j			!atom label found
                exit
              endif
            enddo
            if(jl.ne.j) then			!after "normal" loop exit
              write(*,*) 'atom ',at_name_in,' not found in .PAR'
              stop
            endif
          endif

C *** treat BULK data 
          if(input_method=='BULK') then         !jl identifes label (chem species), jat basis position in CELL
            jat = jl	
            at_ind(1,inrec) = inrec
            at_ind(2:3,inrec) = 0
            at_ind(4,inrec) = jat

            at_pos_in = at_pos_in-at_pos_centre
            do k=1,3
              if(at_pos_in(k)<-a_cell_half(k)) at_pos_in(k) = at_pos_in(k)+a_cell(k,k)   
              if(at_pos_in(k)>=a_cell_half(k)) at_pos_in(k) = at_pos_in(k)-a_cell(k,k)
            enddo
            if(nsuper/=1) at_pos_in = at_pos_in/a_par   !otherwise at_pos stay in Å

            if(j_shell_out==1) then
              at_pos_in2 = at_pos_in2-at_pos_centre
              do k=1,3
                if(at_pos_in2(k)< -a_cell_half(k)) at_pos_in2(k) = at_pos_in2(k)+a_cell(k,k)
                if(at_pos_in2(k)>= a_cell_half(k)) at_pos_in2(k) = at_pos_in2(k)-a_cell(k,k)
              enddo
              if(nsuper/=1) at_pos_in2 = at_pos_in2/a_par
            endif
		
C *** treat CELL data 
          elseif(input_method=='CELL'.or.input_method=='FAST')then		
            at_pos_in = at_pos_in/a_par
            at_pos_in2 = at_pos_in2/a_par
      
            if(j_basis==1) then
              jat = jl        !atom site was identified by name
            else
              if(input_method=='FAST') then
                if(at_no<=2*nsuper) jat = (at_no-1)/nsuper+1
                if(at_no>2*nsuper) jat = mod((at_no-2*nsuper-1),3)+3
              else													!'CELL'
                if(ind_l(jl).eq.1) then
                  jat = i_site(jl,1)	!identify atom site by fractional position 
                else
                  do ii=1,ind_l(jl)
                    pos_inp = at_pos_in-at_base_in(i_site(jl,ii),:)
                    jat = i_site(jl,ii)
                    if(maxval(abs(pos_inp-anint(pos_inp))).lt.eps_x) exit !atom found
                  enddo
                  if(maxval(abs(pos_inp-anint(pos_inp))).gt.eps_x) then
                    write(*,*) 'Identification by position not succeeded: atom, record ',at_label(jl),inrec
                    write(*,*) 'Possible solutions (modify the .PAR file):'
                    write(*,*) '  1/check the ATOMS basis, 2/ try to slightly increase EPS, 3/ use the BULK input method '
                    stop
                  endif
                endif
              endif     
            endif     !j_basis

CC *** calculate cell indices ix,iy,iz from atomic positions shifted to cell origin

            at_ind_in = anint(at_pos_in-at_base_in(jat,:))+at_ind_base  !they will serve as pointers to the right order of atom records
            at_pos_in = at_pos_in-at_ind_shift+at_base_shift			!now the supercell will be centred & basis origin in at_base_shift
C           if(j_shell_out.eq.1) at_pos_in2 = at_pos_in2-at_ind_shift+at_base_shift
                                                                   !both at_ind_base and at_ind_shift are now based on j_centred (line 417) 
            do k=1,3
              if(at_ind_in(k)==0) then
                at_ind_in(k) = at_ind_in(k)+n_row(k)
                at_pos_in(k) = at_pos_in(k)+n_row(k)
C               if(j_shell_out.eq.1) at_pos_in2(k) = at_pos_in2(k)+n_row(k)
C             endif
              elseif(at_ind_in(k)==n_row(k)+1) then
                at_ind_in(k) = at_ind_in(k)-n_row(k)
                at_pos_in(k) = at_pos_in(k)-n_row(k)
C               if(j_shell_out.eq.1) at_pos_in2(k) = at_pos_in2(k)-n_row(k)
              endif
            enddo

            if(j_shell_out.eq.1) then
              at_ind_in2 = anint(at_pos_in2-at_base_in(jat,:))+at_ind_base  !they will serve as pointers to the right order of atom records
              at_pos_in2 = at_pos_in2-at_ind_shift+at_base_shift

              do k=1,3
                if(at_ind_in2(k)==0) then
                  at_ind_in2(k) = at_ind_in2(k)+n_row(k)
                  at_pos_in2(k) = at_pos_in2(k)+n_row(k)
                elseif(at_ind_in2(k)==n_row(k)+1) then
                  at_ind_in2(k) = at_ind_in2(k)-n_row(k)
                  at_pos_in2(k) = at_pos_in2(k)-n_row(k)
                endif
              enddo
            endif     !'SHELL'
       
            jrec = nsuper*(jat-1)+nlayer*(at_ind_in(3)-1)+n_row(1)*(at_ind_in(2)-1)+at_ind_in(1)
          
            if(jrec>n_tot.or.jat<1.or.jat>n_atom) then
						  write(*,*) 'JREC wrong:',jrec,' > ',n_tot,' check n_tot_in, n_row and j_centred in the .PAR'
              write(*,*) 'jrec,at_ind_in',jrec,at_ind_in
              write(*,*) 'at_pos_in',at_pos_in,at_pos_in2
              write(*,*) 'n_tot,nsuper,nlayer,n_row',n_tot,nsuper,nlayer,n_row
						  stop
            endif

            at_ind(1:3,jrec) = at_ind_in
            at_ind(4,jrec) = jat
          endif		!input_method CELL
				
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

C *** accumulate the occupation number and the kinetic energy to refine the real temperature
					if(ifile == nfile_min) nsuper_r(jat) = nsuper_r(jat)+1
          if(n_traj>=1) then
            do k = 1,3
              e_kin(jat,k) = e_kin(jat,k) + at_mass_in*at_veloc_in(k)**2			!at_mass_c(i)=at_veloc_c(4,i)
              if(j_shell==1) e_kin_s(jat,k) = e_kin_s(jat,k) + at_mass_in2*at_veloc_in2(k)**2
            enddo
          endif
				enddo read_loop
				
C *** indexing for the output
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
					write(*,*) '1st snapshot: total of',jrec,' atoms read in',t2-t1,' sec'
				endif
 				
C *** normalize the kinetic energy and get the true temperature  				
				if(n_traj>=1) then
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
							temp = (temp_r_c+temp_r_s)*.5						!hi-T limit: independent C and S vibrations
							if(j_read.le.j_test) write(*,*) 'Hi-T limit: independent C and S vibrations'
						else
							temp = temp_r_c+temp_r_s	!low-T limit: strongly bound C and S vibrations
							if(j_read.le.j_test) write(*,*) 'Low-T limit: strongly bound C and S vibrations'
						endif        
						if(j_read.le.j_test) write(*,*) 'Real temperature: core/shell/total ',temp_r_c,temp_r_s,temp
						if(j_read.le.j_test) write(9,*) 'Real temperature: core/shell/total ',temp_r_c,temp_r_s,temp
					else
						temp = temp_r_c
						temp_r_s = .0
						if(j_read.le.j_test) write(*,*) 'Real temperature: cores only ',temp
						if(j_read.le.j_test) write(9,*) 'Real temperature: cores only ',temp
					endif
				else
					temp = temp_par
					if(j_read.le.j_test) write(*,*) 'Using nominal temperature [K] ',temp
				endif
  				
      
C *** First snapshot in the series only
C *** get the atom position occupation numbers and compare them with those from the .par file
											!analyze in detail the 1st snapshot
				at_name_out = at_name_par(1:n_atom)
				if(input_method=='CELL'.or.input_method=='FAST') then
					at_occup_r = (1.*nsuper_r)/nsuper
					if(i_traj==nt_min.and.ifile==nfile_min) then
						write(*,*) 'Occupancies: nominal 		real'
						do ii=1,n_atom
							write(*,*) '     ',at_name_out(ii),at_occup(ii),at_occup_r(ii)
						enddo
					endif
				else 
					at_occup_r = nsuper_r/real(n_tot)
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
        if(i_save<=9999) then
				  write(file_dat,103) trim(file_par),i_save
				elseif(i_save>=10000) then
				  write(string,'(i8)') i_save
				  file_dat = './data/'//trim(file_par)//'_n'//trim(adjustl(string))//'.dat'
				endif
103     format('./data/',a,'_n',i4.4,'.dat')

				if(i_save==n_save_min.or.i_save==10*(i_save/10)) write(*,*)trim(file_dat)
  			i_save = i_save+1

        call cpu_time(t1)
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
CC  				write(*,*)'at_shells,i_rec',i_rec

  			close(2)
 
				if(j_verb==1.and.ifile==nfile_min)         call cpu_time(t2)
				if(j_verb==1.and.ifile==nfile_min) write(*,*)'New binary output',t2-t1,' sec'
			
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
CC				write(*,*) 'End frame_loop, i_traj',i_traj

			enddo trajectory_loop
			deallocate(at_name,at_ind,e_kin,at_occup_r,nsuper_r,i_series)
				if(j_shell.eq.1) deallocate(e_kin_s)
			
			CALL SYSTEM_CLOCK (COUNT = sc_c2)
      write(*,*) 'Trajectory files finished: ',i_save-n_save_min,' .dat files written in',(sc_c2-sc_c1)/sc_r,' sec (SYS)'
      write(9,*) 'Trajectory files finished: ',i_save-n_save_min,' .dat files written in',(sc_c2-sc_c1)/sc_r,' sec (SYS)'
      write(9,*) 

      stop
      end program mp_dbin55
   
   
   

C **** string conversion to all upper case
C     
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

C     
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

C **** string conversion to all lower case
C          
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
     
C          
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

     
C          
   