
      program mp_dos53
			
C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		  program MP_DOS  1.53   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C ***** vibrational density of states (DOS) from the velocity autocorrelation function
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
C**
C *****  %%%%%%%%%%%%%%%%        program MP_DOS 1.53        %%%%%%%%%%%%%%%%%%%%%%%%
C ***** 
C ***** Ver. 1.50 - start of a new series (identical to 1.45)
C ***** 
C ***** Ver. 1.52 - all data arrays allocatable, no predefined array size limits
C ***** 					- compatible with supercell format in .PAR
C ***** 
C ***** Ver. 1.53 - introduces semi-sequential binary data format 
C *****						- record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****           - CELL, QUICK and BULK data type
C ***** 
C ***** file name convention:
C *****     <master_file>_n<snapshot number>.dat
C *****     example: H05_10K_n0001.dat from H05_10K.txt
C *****
C ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
C ***** record numbers are directly related to cell positions and atom types
C *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
C *****						1 going for the header record
C *****
C ***** velocities are given in Å/ps (angstrom/picosecond), frequency in THz
C ***** atom positions are converted to and recorded in reduced lattice coordinates (x) 
C ***** 

  		use omp_lib
			use singleton

      integer,parameter :: l_rec  =  1024		    !record length in real(4)
      real*8, parameter   :: pi=3.141592653589793238462643383279502884197
      real*8, parameter   :: twopi = 2.d0*pi
      real*4, parameter   :: twopi_s = 2.*pi
      real*4, parameter   :: pi_s = pi
      real*4, parameter   :: sqrt2 = sqrt(2.)
      real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
   		real, parameter   :: h_bar = 6.350776 !DAPS*ps/2Pi Planck const
   		real, parameter   :: h_daps = 39.9031	!DAPS*ps=DAPS/THz  
   		real, parameter   :: n_mass = 1.008665	!atomic units  
			
			logical        :: found_txt,found_ps	
      character(4),allocatable   :: at_label(:),at_name_par(:)
      character(10)  :: subst_name,string,ax_label,c_date,c_time,c_zone
      character(16)  :: sim_type_par,data_type
      character(40)  :: file_master,file_inp,time_stamp,x_title,y_title,at_weight_scheme(2)
      character(40)  :: file_dat,file_dat_t0,file_res,file_ps,file_log
      character(128) :: header_record,plot_title,plot_title_2,scan_title,cwd_path
      character(256) :: header_record_2
      
			integer,allocatable :: i_site(:),at_mask(:)
      real, allocatable :: at_base(:,:),b_coh(:),at_weight(:),conc(:)
      real, allocatable :: cs_plot(:,:),cs_out(:,:),wind(:),ff(:)		
      complex, allocatable :: cs3(:,:,:),cs_inp(:)

      integer ::  j_head_in,ind_rec,hist_ind(3),m_dom1,m_dom2,j_verb,j_name,jm,j1m,j2m,j_bc
      integer ::  i,ii,j,jj,k,l,ios,ier,iflag,j_at,j_rec,nrec,n_tot,nsuper,nrow,nlayer,j_t
      integer ::  nfile,nfile_min,nfile_max,ifile,jfile,n_site,j_dir,i_shift,i_time(8)

      integer :: n_qx,n_qy,n_qz,n_qx_plot,n_qy_plot,j_bin,j_bragg,n_xsym,n_ysym,j_sym,n_q_plot,n_qdisp,n_q1x,n_q2x,n_q1y,n_q2y
      integer :: j_proc,proc_num,proc_num_in,thread_num,thread_num_max,sc_c2,sc_c1
      integer :: ind,j_site,j_q,j_en,n_en,n_int,n_frame,j_disp,j_log,j_txt,j_grid,j_ps,j_freq,n_freq_plot,n_freq_min
      integer ::  j_atom,n_atom,j_at_weight
      real :: at_mass,at_charge,sc_r
      real :: t2,t_step,t_tot,at_weight_sum,wind_sum
      real :: f_plot,ff_plot,f_max_plot,freq_step,f_min,ff_min,f_max
      real :: c_min,c_max,c_min_save,c_max_save,x_plot,y_plot		

C **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
      character(4),allocatable :: at_name(:)
      integer,allocatable ::  nsuper_r(:)
      integer(4),allocatable,target   :: at_ind_in(:,:)
      integer(4), pointer :: at_ind(:,:,:)

      real(4),allocatable ::	at_occup(:)
      real(4),allocatable,target ::	at_vel_in(:,:)
      real(4),pointer :: at_vel_file(:,:,:)
      
      character(16)  :: sim_type,file_title
      integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell,n_traj,n_cond,n_rec,idum
      real(4)        :: rec_zero(l_rec),t_ms,t0,t1,a_par(3),angle(3),temp

C **** PGPLOT stuff
			integer :: PGOPEN,j_xserv,C1,C2,NC
			real :: TR(6),CONTRA,BRIGHT
      CHARACTER*4 XOPT, YOPT
      REAL XTICK, YTICK
      INTEGER NXSUB, NYSUB					

C
C ********************* Initialization *******************************      

			write(*,*) '*** Program MP_DOS 1.53 ** Copyright (C) Jiri Kulda (2022) ***'
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
		    open (9,file=file_log,status ='old',position='append')
		  else
		    open (9,file=file_log,status ='new')
		  endif
                 
			write(9,*)
			write(9,*)
			write(9,*) trim(time_stamp),'  MP_DOS 1.53  ',trim(cwd_path)
			write(9,*) 

C **** read auxiliary file <file_title.par> with structure parameters, atom names and further info
      write(*,*) 'parameter file name (.par will be added)'
      write(*,*) 'Check partial occupation factors - if relevant!'
      read(*,*) file_title
      write(file_inp,101) trim(file_title)
100   format(a)
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
			read(4,*) 		j_proc		!j_proc requested number of OpenMP parallel processes, 0 = automatic maximum
			read(4,*) 		j_name		! j_name: 0 master_string (_ & 4 digits counter), 1 "old" way (_C_ & 4 digits counter, "cores" only)
			read(4,*) 		j_grid		! (0/1) grid overlay on PG_PLOT graphs 
			read(4,*) sim_type_par	!'timestep' for MD, 'static' for DISCUS, phase field
			rewind(4)

      do
        read(4,100) string
        if(string(1:6).eq.'mp_bin') exit	
      enddo
			read(4,*) subst_name
			read(4,*) n_atom
			read(4,*) 				!n_row						!nrow(3) supercell format
			allocate(at_label(n_atom),at_name_par(n_atom))
			allocate(at_base(n_atom,3),b_coh(n_atom),at_mask(n_atom),at_weight(n_atom),conc(n_atom))		!we need this to read .par

			do j=1,n_atom
        read(4,*) at_label(j),at_name_par(j),at_base(j,:),conc(j),b_coh(j)
CC        write(*,*) at_label(j),at_name_par(j),at_base(j,:),conc(j),b_coh(j)
			enddo

			close(4)

C *********************  OpenMP initialization start  *******************************      
C
			proc_num_in = j_proc
			if(j_verb.ge.1)write (*,*) 'PAR proc_num_in          = ',proc_num_in
			thread_num_max = omp_get_max_threads( )								!this gives maximum number of threads available (limited by concurrent tasks??)
			proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
			if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
			call omp_set_num_threads(proc_num_in)
			
			if(j_verb.ge.1) then
CC				write (*,*) 'OMP processors available = ', proc_num
				write (*,*) 'OMP threads max          = ', thread_num_max
				write (*,*) 'OMP processes requested  = ', proc_num_in
			endif


 			if(proc_num_in.gt.1) then
				write(9,*) 'OMP threads max    = ', thread_num_max
				write(9,*) 'OMP processes requested 	= ', proc_num_in
CC				write(9,*) 'OMP processes allocated 		= ', proc_num
CC				write(*,*) 'OMP threads allocated    = ', thread_num
			else
				write(9,*) 'OpenMP not in use'
				write(*,*) 'OpenMP not in use'
			endif
			write(9,*) 
C
C ********************* OpenMP initialization end *******************************      

C *** other initialisations
			j_xserv = 0			
      j_head_in = 0		! if header found will become 1
  	  j_shell = 0
  	  j_log = 0
			j_sym = 0													!don't symmetrise (1/0) the plot
			
      
C *** Generate data file access
			write(*,*) 'Input data filename master (w/o _n): '
			read(*,*)   file_master 
      write(*,*) 'File number: min, max'
      read(*,*)   nfile_min,nfile_max								!here j_step=1
      nfile = nfile_max-nfile_min+1

C *** check for header record      
			write(file_dat_t0,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min
	    open (1,file=file_dat_t0,status ='old',access='direct',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File_1 ',trim(file_dat_t0),' not found! Stop execution.'
				stop
			endif
			
      read(1,rec=1) header_record
      if(header_record(1:4).eq.'time') then
				read(1,rec=1) 
     1		sim_type,file_title,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell,n_cond,n_rec					!char(16): sim_type,file_title, char(4): at_name
      else
      	write(*,*) 'Simulation type: ',sim_type,' not suitable for DOS calculation'
      	stop
      endif 
      close(1)
      
      if(sim_type/=sim_type_par) then
      	write(*,*) 'Simulation type differs in .PAR and in .DAT: ',trim(sim_type),'  ',trim(sim_type_par)
      	write(*,*) 'Continue (1/0)?'
      	read(*,*) j
      	if(j==0) stop
      endif

			data_type = 'cell'
			if(index(sim_type,'bulk')/=0) data_type = 'bulk'
			if(index(sim_type,'quick')/=0) data_type = 'quick'
			write(*,*) 'Data type: ',data_type
 
      if(n_at*n_eq/=n_atom) then
        write(*,*) 'Header:',sim_type,file_title,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq 
      	write(*,*) 'Unit cell number of atoms in .PAR/.DAT not compatible:',n_atom,n_at
      	stop
      endif

      allocate(at_name(n_atom),at_occup(n_atom),nsuper_r(n_atom))

C **** open a second data file to see the time step between recorded snapshots
			write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min+1
	    open (1,file=file_dat,status ='old',access='direct',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File_2 ',trim(file_dat),' not found! Stop execution.'
				stop
			endif

				read(1,rec=1) 
     1		sim_type,file_title,t_ms,t1,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell,
     2    n_cond,n_rec,n_tot,at_name,at_occup,nsuper_r																		!char(16): sim_type,file_title, char(4): at_name
      close(1)

			t_step = t1-t0						!this is the "macroscopic" time step between recorded snapshots
			t_tot = (nfile-1)*t_step

C *********************      Initial info output   ************************

			if(data_type=='cell'.or.data_type=='quick') then
				do j=1,n_at
					if(at_name(j)/=at_name_par(j)) then
						write(*,*) 'Atom names in .PAR and .DAT do not match: ',j,at_name_par(j),at_name(j)
						write(*,*) 'Prefer .PAR? (1/0)'
						read(*,*) ii
						if(ii==1) at_name(j) = at_name_par(j)
					endif
				enddo

				do j=1,n_at
					if(at_occup(j)/=conc(j)) then
						write(*,*) 'Atom occupation factors in .PAR and .DAT do not match: ',j,conc(j),at_occup(j)
						write(*,*) 'Prefer .DAT? (1/0)'
						read(*,*) ii
						if(ii==1)conc(j) = at_occup(j)
					endif
				enddo
			else		!data_type=='bulk'
				do j=1,n_at
					if(at_name(j)/=at_label(j)) then
						write(*,*) 'Atom names in .PAR and .DAT do not match: ',j,at_label(j),at_name(j)
						write(*,*) 'Prefer .PAR? (1/0)'
						read(*,*) ii
						if(ii==1) at_name(j) = at_label(j)
					endif
				enddo
			endif

			write(*,*)
      write(*,*) 'Substance name: ',subst_name	  
			write(*,*) 'Atoms from ',trim(file_inp)
      do j=1,n_atom
        write(*,'(5x,a4,3f8.4,2x,2f8.4)')	at_name(j),at_base(j,:),conc(j),b_coh(j)
        write(9,'(5x,a4,3f8.4,2x,2f8.4)')	at_name(j),at_base(j,:),conc(j),b_coh(j)
      enddo
			
			allocate(cs_inp(nfile),cs3(nfile,n_tot,3),wind(nfile))
			allocate(at_ind_in(4*n_tot,nfile),at_vel_in(4*n_tot,nfile))
 			
C *** cycle over snapshot files to accumulate input data
C
			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
			sc_r = 1./sc_r
			ifile = 0				
      write(*,*) 'Input files:'

      file_loop: do jfile=nfile_min,nfile_min+nfile-1

C ***  open the binary MD snapshot file
				write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),jfile
CC        write(*,*) file_dat
        open (1,file=file_dat,action='read',status ='old',access='direct',form='unformatted',recl=4*l_rec,iostat=ios)
        if(ios.ne.0) then
        	write(*,*) 'File ',trim(file_dat),' not opened! Exit reading / Stop execution (1/0)?'
        	read(*,*) jj
        	if(jj==1) exit file_loop
        	if(jj==0) stop
        endif
 
 				ifile = ifile+1				       	
        if(jfile==1.or.jfile==10*(jfile/10)) write(*,*) file_dat

C *** for the whole sequence the same header status is assumed 
C *** read the input .dat file as a whole, only at_name, at_pos and idom are really needed

				j_rec = 1				!1 goes for the header, which is not read
				do j=1,n_rec-1
					j_rec = j_rec+1
					read(1,rec=j_rec) at_ind_in((j-1)*l_rec+1:j*l_rec,ifile)			
				enddo	
				j_rec = j_rec+1
				read(1,rec=j_rec) at_ind_in((n_rec-1)*l_rec+1:4*n_tot,ifile)			
				
				j_rec = j_rec+n_rec			!skip the at_pos records
				
				do j=1,n_rec-1
					j_rec = j_rec+1
					read(1,rec=j_rec) at_vel_in((j-1)*l_rec+1:j*l_rec,ifile)			
				enddo	
				j_rec = j_rec+1
				read(1,rec=j_rec) at_vel_in((n_rec-1)*l_rec+1:4*n_tot,ifile)	

				close(1)
			enddo file_loop

				at_ind(1:4,1:n_tot,1:nfile) => at_ind_in				
				at_vel_file(1:4,1:n_tot,1:nfile) => at_vel_in

			
			if(ifile.ne.nfile) nfile = ifile

			CALL SYSTEM_CLOCK (COUNT = sc_c2)
			write(*,*) nfile,' files read in',(sc_c2-sc_c1)*sc_r,' sec SYS time'
			write(*,*) 

			jfile = 1				!now to be used as index for the successive output files
						
			call cpu_time(t1)

			f_max_plot = .5/t_step
      freq_step = 1./t_tot
      n_freq_plot = .5*t_tot/t_step
      n_freq_plot = n_freq_plot+1
      write(*,*) 'n_freq_plot',n_freq_plot

			write(*,*) 't_range [ps] =',t_tot,'t_step [ps] =',t_step
			write(*,*) 'freq_max [THz] =',f_max_plot,'freq_step_min [THz] =',1./t_tot
			write(*,*) 

					write(9,*) 
					write(9,*) 'Input files:',nfile,' starting from ',trim(file_dat_t0)
					write(9,*) '  t0,t_step,n_row,n_atom:',t0,t_step,n_row,n_atom
					write(9,*) '  temperature [K]',temp
        	write(9,121) (at_name(j),j=1,n_atom)
121				format(' Atoms:        ',20(a,3x))
        	write(9,122) (conc(j),j=1,n_atom)
122				format(' Occupations:',20f7.3)
        	write(9,123) (b_coh(j),j=1,n_atom)
123				format(' b_coh:      ',20f7.3)
					write(9,*) 

			allocate(cs_plot(n_freq_plot,n_atom+1),cs_out(n_freq_plot,n_atom+1),ff(n_freq_plot))
			
			do i=1,nfile
				wind(i) = .5*(1.-cos(twopi*i/real(nfile+1)))				!use the Hann window
			enddo
			wind_sum = sum(wind)

			write(*,*) 'OpenMP FFT of velocities'
			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				
!$omp parallel shared(at_vel_file,cs3,wind,n_tot,nfile) private(cs_inp)
CC				if(j_verb==1)	write(*,*) 'OMP threads available:', omp_get_num_threads( )	
!$omp do
				do ind=1,n_tot
					do j=1,3
						cs_inp(:) = at_vel_file(j,ind,:)*wind
						cs3(:,ind,j) = fft(cs_inp,inv=.false.)/sqrt(1.*nfile)
					enddo
				enddo
!$omp end do
!$omp end parallel

CC				CALL SYSTEM_CLOCK (COUNT = sc_c2)
CC				if(j_verb.ge.1) write(*,*) (sc_c2-sc_c1)/sc_r

			cs_plot = .0
			CALL SYSTEM_CLOCK (COUNT = sc_c1)				

!$omp parallel shared(cs_plot,at_ind,cs3,n_tot,n_freq_plot) 
!$omp do
				do i=1,n_freq_plot
					do ind=1,n_tot
						cs_plot(i,at_ind(4,ind,1)) = cs_plot(i,at_ind(4,ind,1))+sum(cs3(i,ind,:)*conjg(cs3(i,ind,:)))		!1: indices in all files must be the same
					enddo
				enddo
!$omp end do
!$omp end parallel

CC				CALL SYSTEM_CLOCK (COUNT = sc_c2)
CC				if(j_verb.ge.1) write(*,*) 'OpenMP time integration',(sc_c2-sc_c1)/sc_r	

				n_int = 1
				n_frame = nfile-n_int+1

C *** cycle over DOS with different atom at_weights

			at_weight_scheme(1) = 'Uniform weights'
			at_weight_scheme(2) = 'Neutron weights'
			at_mask = 1

								
			at_weights_loop: do	

CC				write(*,*) 'Atom weights: 1= uniform, 2= neutron coherent (b_c^2); 0= EXIT, negative = EDIT'	
				write(*,*) 'Atom weights: 1= uniform, 2= neutron coherent (b_c^2); 0= EXIT'	
				do
					read(*,*) jj
					if(jj==0) exit at_weights_loop
					if(abs(jj)==1.or.abs(jj)==2) exit
					write(*,*) 'Input out of range, repeat ...'
				enddo

				j_at_weight = abs(jj)
				at_weight = .0
				if(j_at_weight==1) then
					at_weight(1:n_atom) = 1.
				elseif(j_at_weight==2) then
					at_weight(1:n_atom) = b_coh(1:n_atom)*b_coh(1:n_atom)
				endif
				write(9,*)
				write(*,'(1x,"Atoms:         ",50a8)')  (at_name(i),i=1,n_atom)
				write(9,'(1x,"Atoms:         ",50a8)')  (at_name(i),i=1,n_atom)
				write(*,'(1x,"Atoms no.:  ",50i8)') (i,i=1,n_atom)
				write(9,'(1x,"Atoms no.:  ",50i8)') (i,i=1,n_atom)
				write(*,113) trim(at_weight_scheme(j_at_weight)),(at_weight(i),i=1,n_atom)
				write(9,113) trim(at_weight_scheme(j_at_weight)),(at_weight(i),i=1,n_atom)
113     format(1x,a,' at_weights: ',50f8.4)								
				write(*,*) 'Actual masks:'
				write(*,103) (at_mask(i),i=1,n_atom)
103     format((50i3))
				if(jj<0) then
					write(*,*)'Type in new ones:'
					do
						read(*,*)(at_mask(i),i=1,n_atom)
						if(any(at_mask(1:n_atom).ge.0).and.any(at_mask(1:n_atom).le.1)) exit
						write(*,*) 'Input out of range, repeat ...'
					enddo
				endif

				at_weight = at_weight*at_mask
				at_weight_sum = sum(conc(1:n_atom)*at_weight(1:n_atom)*nsuper_r(1:n_atom))*wind_sum**2			
C
C **** generate the plot data
C	
				do i=1,n_freq_plot
					cs_out(i,1:n_atom) = cs_plot(i,1:n_atom)*at_weight
					cs_out(i,n_atom+1) = sum(cs_out(i,1:n_atom))
				enddo

				cs_out = cs_out*nfile/at_weight_sum

				f_min = 0.
				f_max = f_max_plot
								
C
C **** generate the plot
C	
				if(j_xserv==0) then
					j_xserv = PGOPEN('/xserv')
					if(j_verb.ge.1) write(*,*) 'j_xserv',j_xserv
					CALL PGASK(.FALSE.)     ! would not ask for <RET>
					CALL PGPAP(7.0,1.)     ! define the plot areaCC						CALL PGERAS
					CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
					CALL PGSCRN(1, 'black', IER)
				else
					call PGSLCT(j_xserv)
				endif
C
C **** do the plots
C

				do j=1,n_freq_plot
					ff(j) = f_min+(j-1)*(f_max-f_min)/real(n_freq_plot-1)
				enddo
				c_min = .0
				c_max = maxval(cs_out)

				write(plot_title,106) trim(subst_name),temp,trim(at_weight_scheme(j_at_weight))
106     format(' Vibrational DOS:  ',a,'  T =',f6.1,' [K]  ',a)   
				write(*,*) plot_title
				do
					CALL PGSCH(1.)					!set character height					
					CALL PGSCI (1)  !white
					CALL PGSLS (1)  !full
					CALL PGENV(f_min,f_max,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes w/o grid
					CALL PGLAB('E [THz]', 'DOS [arb. units]',trim(plot_title))  !put the axis labels
          CALL PGSCI (1)  !white
          CALL PGSLW(5)			!operates in steps of 5
					CALL PGLINE(n_freq_plot,ff,cs_out(1:n_freq_plot,n_atom+1))  !plots the total DOS
					 write(plot_title_2,114) trim(subst_name)
114        format(a)   
					x_plot = f_min+.85*(f_max-f_min)
					y_plot = .9*c_max
					CALL PGSTBG(0)																				 !erase graphics under text
					CALL PGTEXT (x_plot,y_plot,plot_title_2)
          CALL PGSLW(2)
					do j=1,n_atom
					 write(plot_title_2,115) at_name(j)
115        format(a)   
					 y_plot = (.9-.06*j)*c_max
						CALL PGSCI (j+1)  !red-green-blue
						CALL PGLINE(n_freq_plot,ff,cs_out(1:n_freq_plot,j))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title_2)
          	CALL PGSLW(2)			
					enddo
					c_min_save = c_min
					c_max_save = c_max
					write(*,*) 'c_min,c_max (0 0 exit, 1 1 exit & save (PS, txt))',c_min,c_max
					read(*,*)c_min,c_max								
					if(c_min==c_max) then
						j_ps=0
						if(c_min==1) j_ps=1
						exit
					endif
				enddo

C **** make an optional hardcopy and .txt output
				if (j_ps.eq.1) then				!j_ps	
					jfile = 1
					do						!look for existing .txt files to continue numbering
						write(file_res,108) trim(file_master),nfile_min,nfile,jfile
108   			format(a,'_',i4.4,'_',i4.4,'_dos_',i2.2,'.txt')
						inquire(file=file_res,exist=found_txt)
						write(file_ps,117) trim(file_master),nfile_min,nfile,jfile	
117   			format(a,'_',i4.4,'_',i4.4,'_dos_',i2.2,'.ps')
						inquire(file=file_ps,exist=found_ps)

						if(.not.found_txt.and.(.not.found_ps)) exit				
						jfile = jfile+1
						if(jfile==100) then
							write(*,*)'Tidy up .txt/.ps files to restart count from 01 and type [RET]'
							read(*,*)
							jfile = 1
						endif							
					enddo

					write(9,*) 'Output files:  ',trim(file_ps),'   ',trim(file_res)
					write(9,*)					

					write(*,*)
					write(*,*) 'Text output:  ',trim(file_res)
					write(*,*) 'PS output  :  ',trim(file_ps)
					write(*,*)
CC					plot_title = trim(plot_title)//('  '//trim(file_ps))	
								
C **** Output the DOS to a text file (linear scale)			
					open (3,file=file_res)																		!open the output file
					write(3,*) '*****    MP_DOS53: total and partial densities of vibrational states     *****'
					write(3,*) 
					write(3,*) 'Input files:  ',trim(file_dat_t0),' to ',trim(file_dat)
					write(3,*) '  t0,t_step,n_row,n_atom:',t0,t_step,n_row,n_atom
					write(3,*) '  temperature [K]',temp
					write(3,*) 
					write(3,121) (at_name(j),j=1,n_atom)
					write(3,122) (conc(j),j=1,n_atom)
					write(3,123) (at_weight(j),j=1,n_atom)
					write(3,*) 
					write(3,124) (at_name(j),j=1,n_atom)
					do i=1,n_freq_plot
						write(3,109) ff(i),cs_out(i,n_atom+1),(cs_out(i,j),j=1,n_atom) 
					enddo
					close(3)
109   		format(1x,f6.3,2x,f8.3,2x,1000(1x,f8.3))			
124				format(' f [THz] ','   Total  ',20(a9))

C **** Prepare and plot the same on .PS		
CC					write(file_ps,118) trim(file_master),nfile_min,nfile,jfile
CC118   		format(a,'_dos_',i4.4,'_',i4.4,'_',i2.2,'.ps/VCPS')
					IF (PGOPEN(file_ps//'/VCPS').LE.0) STOP

					CALL PGASK(.FALSE.)     ! would not ask for <RET>
					CALL PGPAP(7.0,1.)     ! define the plot areaCC						CALL PGERAS
					CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
					CALL PGSCRN(1, 'black', IER)

					c_min = c_min_save 
					c_max = c_max_save

					CALL PGSCH(1.)					!set character height					
					CALL PGSCI (1)  !white
					CALL PGSLS (1)  !full
					CALL PGENV(f_min,f_max,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
CC					CALL PGLAB('E [THz]', 'DOS [arb. units]',trim(plot_title))  !put the axis labels
					CALL PGLAB('E [THz]', 'DOS [arb. units]',' ')  !put the axis labels
					CALL PGSCH(.7)					!set character height					
					CALL PGMTXT ('T', 3., .5, .5, trim(plot_title))			!put plot title on 2 lines
					CALL PGMTXT ('T', 1., .5, .5, trim(file_ps))
					CALL PGSCH(1.0)					!set character height					
          CALL PGSCI (1)  !white
          CALL PGSLW(5)			!operates in steps of 5
					CALL PGLINE(n_freq_plot,ff,cs_out(1:n_freq_plot,n_atom+1))  !plots the total DOS
					 write(plot_title_2,114) trim(subst_name)
CC114        format('DOS_total ',a)   
					x_plot = f_min+.75*(f_max-f_min)
					y_plot = .9*c_max
					CALL PGSTBG(0)																				 !erase graphics under text
					CALL PGTEXT (x_plot,y_plot,plot_title_2)
          CALL PGSLW(2)
					do j=1,n_atom
					 write(plot_title_2,115) at_name(j)
CC115        format('DOS_part  ',a)   
					 y_plot = (.9-.06*j)*c_max
						CALL PGSCI (j+1)  !red-green-blue
						CALL PGLINE(n_freq_plot,ff,cs_out(1:n_freq_plot,j))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title_2)
          	CALL PGSLW(2)			
					enddo

				CALL PGCLOS
				endif		!j_ps = 1
			enddo at_weights_loop
			CALL PGEND
C
C **** deallocate and cycle loops
C				
				deallocate(cs_plot,cs_out)				!wind,
				deallocate(cs3,ff)

				flush(9)
CC			deallocate(cs)				!,cross_sec

			close(3)
			close(9)
            
      end program mp_dos53
      
      