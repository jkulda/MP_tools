
      program mp_sqom53

C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		 program MP_SQOM 1.53   		 %%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   calculates and plots scattering functions (S(Q), S(Q,w)) from simulated data
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
C *****  %%%%%%%%%%%%%%%%       program MP_SQOM 1.53      %%%%%%%%%%%%%%%%%%%%%%%%
C ***** 
C ***** Ver 1.51  - originates from mp_sqom49.f of 2022-01-16 16:57 
C ***** 					- optimized algorithm using coherence profile integration w/o OMP 
C ***** 					- other parts parallelised via OMP without touching FINUFFT (OMP already in)
C ***** 
C ***** Ver 1.52  - major changes against v. 1.51 
C ***** 					- binary data input follows the standard MP_TOOLS way 
C ***** 					- the NUFFT in Q is atom-resolved to enable playing with atom masks
C ***** 					- new waypoint towards various plotting options
C ***** 					- linear plots skipped from plotting options (remaining a secret one)
C ***** 					- gives possibility to use the Hann window for space FT
C ***** 					- employs a speckle filter
C ***** 
C ***** Ver 1.53  - introduces semi-sequential binary data format 
C *****						- record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****           - CELL, QUICK and BULK data type
C ***** 
C ***** indexing convention: supercell(u,v,w) = at_ind(j,k,l) = supercell(j_pos,j_row,j_layer)
C ***** record numbers are directly related to cell positions and atom types
C *****    jrec = 1+nsuper*(jat-1)+nlayer*(at_ind(3)-1)+nrow*(at_ind(2)-1)+at_ind(1)
C *****						1 going for the header record
C *****
C ***** atom positions are converted to and recorded in reduced lattice coordinates (x) 
C ***** 
  		use omp_lib
			use singleton
CC  		use mp_nopgplot					! uncomment when not able to use PGPLOT, compile and link together with the module mp_nopgplot.f
			
      integer,parameter :: l_rec  =  1024		    !record length in real(4)
      real(8), parameter   :: pi=3.141592653589793238462643383279502884197
      real(8), parameter   :: twopi = 2.d0*pi
      real(4), parameter   :: twopi_s = 2.*pi
      real(4), parameter   :: pi_s = pi
      real(4), parameter   :: sqrt2 = sqrt(2.)
      real, parameter   :: k_B = .831444 !DAPS/K Boltzmann's constant 0.08617333262145 meV/K   
   		real, parameter   :: h_bar = 6.350776 !DAPS*ps/2Pi Planck const
   		real, parameter   :: h_daps = 39.9031	!DAPS*ps=DAPS/THz  
   		real, parameter   :: n_mass = 1.008665	!atomic units  
			
			logical        :: found_txt,found_ps,t_single,inv,read_ind,non_zero,inside
      character(1)   :: xyz(3)=(/'X','Y','Z'/)
      character(4)   :: atom,ps_out(2)
      character(5)   :: c_dir(5)=(/'[00X]','[0X0]','[0XX]','[-XX]','[XXX]'/)
      character(10)  :: subst_name,string,ax_label,c_date,c_time,c_zone,mode,data_type
      character(40)  :: file_master,file_inp,time_stamp,x_title,y_title,masks,at_weight_scheme(3)
      character(40)  :: file_dat,file_dat_t0,file_res,file_ps,file_log,x_file_name
      character(128) :: header_record,plot_title1,plot_title2,plot_title,scan_title,cwd_path
      character(256) :: header_record_2

      character(4),allocatable   :: at_label(:),at_name_par(:),at_name_ext(:)
      character(10),allocatable ::	corr_name(:)
      
      integer, allocatable :: i_site(:)		!index attributing site number to atom number
      integer,allocatable ::  ind_at(:),at_mask(:),map_unit(:) 
 
      real, allocatable :: at_base(:,:),conc(:),b_coh(:),b_coh_save(:),x_ffpar(:,:),x_ffq(:,:,:)			!s_base(:,:),q(:,:,:),at_vel_file(:,:)		!atom basis, site basis fractional coordinates
 			real, allocatable :: cs_plot(:,:),cs_out(:,:),csp(:),qq(:),ff(:),t_wind(:) 								!q_at_pos(:),q_at_vel(:),at_pos_perp(:),rq_fft(:)
			complex, allocatable :: u_x(:,:),u_y(:,:),u_qx(:,:),u_qy(:,:)
			 			
 			real(8), allocatable :: xf(:),yf(:)			
      real(8) ::     eps,twopi_nr(3)
      integer(8) ::  n_tot,nsuper8,n_qx8,n_qy8,n_fft
      integer(8), allocatable :: nul_opt,n_sup(:)		!nul_opt (since unallocated) used to pass a null pointer to FINUFFT...
      complex(8), allocatable :: cf(:),ampl_tot(:),ampl_tot_2d(:,:)												! this is equivalent to complex*16

      complex :: 			bragg,cs1,cs2_mean,c_f,phase_bz
      complex, allocatable :: om_phase(:),cs2(:),cs3(:),cs_mean(:,:,:),cs5(:,:),cs6(:,:),cs7(:,:),cs4(:,:)
      complex, allocatable, target :: cs_atom(:,:,:,:),cs(:,:,:)
      complex, pointer     :: cs_p1(:),cs_p2(:)

      integer ::  j_head_in,ind_rec,hist_ind(3),m_dom1,m_dom2,j_verb,jm,j1m,j2m,j_bc,j_exit,j_xray
      integer ::  i,ii,iii,j,jj,i1,i2,i3,k,kk,l,ll,ios,ier,iflag,j_at,j_rec,nrec,nsuper,nrow,nlayer,n_r1,n_r2,n_row_save(3)
      integer ::  nfile,nfile_0,nfile_min,nfile_max,ifile,jfile,nfile_step,nqx_min,nqy_min
      integer ::  n_site,j_dir,i_shift,i_time(8),j_weight
      integer :: n_qx,n_qy,n_qz,n_qdisp,n_q1x,n_q2x,n_q1y,n_q2y,nq_tot,n_bz,n_freq_plot,n_freq_min,n_freq_max
      integer :: e1(3),e2(3),ev(3),j_scan,j_freq,i_plot,j_plot,n_plot,i_step
      integer :: sc_c1,sc_c2,sc_m,j_proc,proc_num,proc_num_in,thread_num,thread_num_max
      integer :: ind,j_site,j_q,n_int,n_frame,j_disp,j_mask,j_logscale,j_grid,j_txt,j_ps
      integer :: j_name,j_atom,n_atom,at_no,j_atc1,j_atc2,j_ft,j_coh,j_zero,j_oneph,j_fft

      real :: at_mass,at_charge,bc_in,s_trig,q_sq,rn
      real :: at_displ,at_pos(3),at_veloc(3),at_pos_min(3),at_pos_max(3),at_pos_centre(3)
      real :: t2,t3,t4,t_sum,t_step,t_step0,t_tot,tt0,tt,temp,xx,f_width,t_width
      real :: f_plot,ff_plot,f_max_plot,freq_step,f_min,ff_min,f_max,ff_max,f_map_min,f_map_max,f_map_step,f_phase
      real :: arg,q1_x,q1_y,q2_x,q2_y,qx_min,qx_max,qy_min,qy_max,q_min,q_max,qq_min,qq_max,qq_plot,dq(3),dq_p(3),q_step
      real :: tr_shift_x,tr_shift_y,d_x,d_y,q_x,q_y
      real :: e1_norm,e2_norm,ev_norm,tpq_center(3),tpe1,tpe2,tpev
      real :: q1(3),q2(3),q1_3d(3),q2_3d(3),q_v(3),q_center(3),q_nx,q_ny,wtime
      real :: c_min,c_max,c_min_save,c_max_save,sc_r		!,tau(3),c_vert(3),b_perp(3),q_tot(3),q_unit(3),q_norm

C **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
C
      character(4),allocatable :: at_name_inp(:)
      integer(4),allocatable ::   nsuper_r(:)
      integer(4),allocatable,target   :: at_ind_in(:)
      integer(4),pointer :: at_ind(:,:,:)
      
      real(4),allocatable ::	at_occup(:)
      real(4),allocatable,target ::	at_pos_in(:,:)
 			real(4), pointer :: at_pos_file(:,:,:,:)
 			real(4), pointer :: at_pos_file_bulk(:,:,:)

      character(16)  :: sim_type,file_title
      integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell,n_traj,n_cond,n_rec,n_tot_in,idum
      real(4)        :: rec_zero(l_rec),t_ms,t0,t1,a_par(3),angle(3),temp_r

C *** PGPLOT stuff
			integer :: PGOPEN,C1,C2,NC,plot_unit,map_unit_1,map_unit_2
			real :: TR(6),CONTRA,BRIGHT,p_size
      CHARACTER*4 XOPT, YOPT
      REAL XTICK, YTICK
      INTEGER NXSUB, NYSUB					

C
C ********************* Initialization *******************************      

			write(*,*) '*** Program MP_SQOM 1.53 ** Copyright (C) Jiri Kulda (2022) ***'
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
		    open (9,file=file_log,status ='new')
		  endif

			write(9,*) 
			write(9,*) 
			write(9,*) trim(time_stamp),'  MP_SQOM 1.53  ',trim(cwd_path)
			write(9,*) 

C **** read auxiliary file <file_title.par> with structure parameters, atom names and further info
      write(*,*) 'parameter file name (.par will be added)'
      write(*,*) 'Check partial occupation factors - if relevant!'
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
			read(4,*) 		j_proc		!j_proc requested number of OpenMP parallel processes, 0 = automatic maximum
			read(4,*) 		j_name		! j_name: 0 master_string (_ & 4 digits counter), 1 "old" way (_C_ & 4 digits counter, "cores" only)
			read(4,*) 		j_grid		! (0/1) grid overlay on PG_PLOT graphs 
			read(4,'(a)') sim_type	!'timestep' for MD, 'static' for DISCUS, phase field etc.
			rewind(4)
			
      do
        read(4,'(a)') string
        if(string(1:6).eq.'mp_bin') exit	
      enddo
      
			read(4,*) subst_name
			read(4,*) n_atom
			read(4,*) n_row 
			allocate(at_label(n_atom),at_name_ext(n_atom),at_name_par(n_atom),at_base(n_atom,3),at_mask(n_atom))
			allocate(conc(n_atom),b_coh(n_atom),x_ffpar(n_atom,9),SOURCE=.0)												!we need this to read .par
			do j=1,n_atom
        read(4,*) at_label(j),at_name_par(j),at_base(j,:),conc(j),b_coh(j),at_name_ext(j)	!for BULK the at_base and conc are not significant
        if(at_name_ext(j)=='n'.or.at_name_ext(j)=='N') at_name_ext(j)=''
			enddo
			n_row_save = n_row
			b_coh_save = b_coh						!keep the original b_coh values for reference
			rewind(4)

      do
        read(4,'(a)') string
        if(string(1:7).eq.'mp_sqom') exit	
      enddo
			read(4,*) n_int					!n_int time-integration range (number of snapshots): 0 = automatic choice of Nfile/2
			read(4,*) j_zero				!j_zero=1 eliminate points with an exactly zero coordinate (0=OFF) [in SQOM521] otherwise s_trig=2  the triger level of the speckle filter
			s_trig = .0
			read(4,*) j_weight			!j_weight=1 uniform weighting, 2 neutron b_coh; 3 Xray form-factors 
			read(4,*) 							!j_xray				!j_xray=1	load Xray formfactor table
			read(4,*) j_logscale		!j_logscale=1  maps log scale
			read(4,*) j_txt					!j_txt=1 	accompany PS maps output by intensities in .txt file
			read(4,*) j_oneph				!start in NUFFT (0), 1phonon (1), runtime toggle 0,1
			close(4)

			j_xray = 0
			j_coh = 1

C *********************  OpenMP initialization start  *******************************      
C
			proc_num_in = j_proc
			if(j_verb.ge.1)write (*,*) 'PAR proc_num_in          = ',proc_num_in
			thread_num_max = omp_get_max_threads( )								!this gives maximum number of threads available (limited by concurrent tasks??)
			proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
			if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
			call omp_set_num_threads(proc_num_in)
CC			thread_num = omp_get_num_threads( )				!this gives threads available at this moment: here always = 1 as we are outside of a PARALLEL range
			
			if(j_verb.ge.1) then
				write (*,*) 'OMP processes available  = ', proc_num
				write (*,*) 'OMP threads maximum      = ', thread_num_max
				write (*,*) 'OMP processes requested  = ', proc_num_in
			endif


 			if(proc_num_in.gt.1) then
				write(9,*) 'OMP threads maximum = ', thread_num_max
				write(9,*) 'OMP processes requested 	= ', proc_num_in
			else
				write(9,*) 'OpenMP not in use'
				write(*,*) 'OpenMP not in use'
			endif
			write(9,*) 
C
C ********************* OpenMP initialization end *******************************      

C *** other initialisations
			ps_out = (/'OFF ','ON  '/)			!PGPLOT
			plot_unit = 0
			map_unit_1 = -1	
			map_unit_2 = -1
			p_size = 7.											!PGPLOT
			eps = 1.e-6						!finufft
			iflag = -1						!finufft			
      j_head_in = 0		! if header found will become 1
			at_weight_scheme(1) = 'Unit weights '
			at_weight_scheme(2) = 'Neutron weights'
			at_weight_scheme(3) = 'Xray weights'
			at_mask = 1
			f_max = .0


C *********************      Initial info output   ************************
			write(*,*)
      write(*,*) 'Substance name: ',subst_name	  
			write(*,*) 'Atoms from ',trim(file_inp)
      do j=1,n_atom
        write(*,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),conc(j),b_coh(j)
        write(9,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),conc(j),b_coh(j)
      enddo

C *** Generate data file access
			if(j_name==0)then			!THE NEW WAY
				write(*,*) 'Input data file_master: '
			else
				write(*,*) 'Input data file_master (_C_ will be added): '
			endif
			read(*,*) file_master 
			if(j_name==1) file_master = trim(file_master)//'_c'
						
      write(*,*) 'Read data files number min, max (0 0 no numbers, single file): '
      read(*,*)   nfile_min,nfile_max
      nfile_step = 1
			t_single = nfile_min==0.and.nfile_max==0

			if(t_single)then
				nfile_min = 1		!for conformity with loop control conventions
				nfile_max = 1
			endif
			
      nfile = ((nfile_max-nfile_min)/nfile_step)+1
      
			if(j_name==0.and.t_single)then
				write(file_dat_t0,'("./data/",a,".dat")') trim(file_master)
			else
				write(file_dat_t0,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min
			endif

	    open (1,file=file_dat_t0,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_dat_t0),' not found! Stop execution.'
				stop
			endif
			
C *** check for header record      
      read(1,rec=1) header_record
      if(header_record(1:4).eq.sim_type(1:4)) then

				read(1,rec=1) 
     1		sim_type,file_title,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell,n_cond,n_rec		
				
        j_head_in = 1
      else
      	write(*,*) 'header record wrong or missing'
      	write(*,*) trim(header_record)
      	stop
      endif 
      
			data_type = 'cell'
			if(index(sim_type,'bulk')/=0) data_type = 'bulk'
			if(index(sim_type,'quick')/=0) data_type = 'quick'
			write(*,*) 'Data type: ',data_type
      
      if(n_at/=n_atom.and.(data_type=='cell'.or.data_type=='quick')) then
        write(*,*) sim_type,file_title,n_traj,n_cond,t0,t_step0,a_par,n_row,n_at,j_force,temp  
      	write(*,*) 'Unit cell number of atoms in .PAR/.DAT not compatible:',n_atom,n_at
      	stop
      else
        n_atom = n_at
      endif
      allocate(at_name_inp(n_atom),at_occup(n_atom),nsuper_r(n_atom),n_sup(n_atom),ind_at(n_atom))

			read(1,rec=1) 
     1		sim_type,file_title,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell,n_cond,n_rec,n_tot_in,			
     2    at_name_inp,at_occup,nsuper_r																		!char(16): sim_type,file_title, char(4): at_name
      close(1) 
 
C **** for an MD snapshot sequence open a second data file to see the time step between recorded snapshots
			if(sim_type(1:4).eq.'time'.and.(.not.t_single)) then
				write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min+nfile_step
				open (1,file=file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
				if(ios.ne.0) then
					write(*,*) 'File ',trim(file_dat),' not found! Stop execution.'
					stop
				endif

				read(1,rec=1) 
     1		sim_type,file_title,t_ms,t1,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell,n_cond,n_rec,jj,			
     2    at_name_inp,at_occup,nsuper_r																		!char(16): sim_type,file_title, char(4): at_name
				close(1)

				t_step = t1-t0						!this is the "macroscopic" time step between recorded snapshots
				t_tot = (nfile-1)*t_step
			else
				t_step = .0						!for total scattering (unrelated snapshots)
				t_tot = .0
			endif

			do j=1,n_at
				if(at_name_inp(j)/=at_name_par(j)) then
					write(*,*) 'Atom names in .PAR and .DAT do not match: ',j,at_name_par(j),at_name_inp(j)
					write(*,*) 'Prefer .DAT? (1/0)'
					read(*,*) ii
					if(ii==1) at_name_par(j) = at_name_inp(j) 
				endif
			enddo
			
 			nrow = n_row(1)
      nlayer = n_row(1)*n_row(2) 
      nsuper = n_row(1)*n_row(2)*n_row(3) 
			n_sup(:) = nsuper_r(:)					!n_sup is a 64bit copy needed for FINUFFT, while storage is 32bit for space economy

      if(data_type=='cell'.or.data_type=='quick') then
				n_tot = n_atom*nsuper		! for allocation purposes, allowing for vacancies; the header contains counted n_tot_in
			else
				n_tot = n_tot_in			
			endif

C *** cycle over snapshot files to accumulate input data
C

			allocate(at_pos_in(4*n_tot,nfile),at_ind_in(4*n_tot))			

      write(*,*) 'Input files:'

			CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
			sc_r = 1./sc_r

			ifile = 0
      file_loop: do jfile=nfile_min,nfile_max,nfile_step

C ***  open the binary MD snapshot file
				if(j_name==0.and.t_single)then
					write(file_dat,'("./data/",a,".dat")') trim(file_master)
				else
					write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),jfile
				endif

        open(1,file=file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
        if(ios.ne.0) then
        	write(*,*) 'File ',trim(file_dat),' not opened! IOS =',ios
        	write(*,*) 'Skip(1), stop execution(0)?'
        	read(*,*) jj
        	if(jj==1) exit file_loop
        	if(jj==0) stop
        endif
				ifile = ifile+1				
        if(jfile==nfile_min.or.ifile==100*(ifile/100)) write(*,*) file_dat
CC       if(jfile==nfile_min.or.ifile==10*(ifile/10)) then
CC 					CALL SYSTEM_CLOCK (COUNT = sc_c2)				!, COUNT MAX = sc_m)
CC          write(*,*) (sc_c2-sc_c1)*sc_r,file_dat				!,jfile,ifile
CC        endif

				ind_at(1) = 0
				do j = 2,n_atom
					ind_at(j) = ind_at(j-1)+nsuper_r(j-1)
				enddo


				j_rec = 1	
				do j=1,n_rec-1
					j_rec = j_rec+1
					if(jfile==nfile_min) read(1,rec=j_rec) at_ind_in((j-1)*l_rec+1:j*l_rec)			!only needed for 1_phonon with CELL and QUICK
				enddo	
				j_rec = j_rec+1
				if(jfile==nfile_min) read(1,rec=j_rec) at_ind_in((n_rec-1)*l_rec+1:4*n_tot)			
				
				do j=1,n_rec-1
					j_rec = j_rec+1
					read(1,rec=j_rec) at_pos_in((j-1)*l_rec+1:j*l_rec,ifile)			
				enddo	
				j_rec = j_rec+1
				read(1,rec=j_rec) at_pos_in((n_rec-1)*l_rec+1:4*n_tot,ifile)					
				close(1)
						
			enddo file_loop

			if(data_type=='cell'.or.data_type=='quick') then
				at_ind(1:4,1:nsuper,1:n_atom) => at_ind_in
				at_pos_file(1:4,1:nsuper,1:n_atom,1:nfile) => at_pos_in
			else
				at_pos_file_bulk(1:4,1:n_tot,1:nfile) => at_pos_in
			endif				

			if(ifile.ne.nfile) then						!update nfile in case of shorter sequence
				nfile = ifile
			endif
			nfile_0 = nfile				!for later memory
			
				CALL SYSTEM_CLOCK (COUNT = sc_c2)
				write(*,*) nfile,' files read in ',(sc_c2-sc_c1)*sc_r, ' sec SYS time'


C *** set the time and frequency constants and the time-integration window

			if((nfile-1)*t_step/=t_tot) then				!nfile might have changed ...
				write(*,*)'Updating t_tot from',t_tot,' to',(nfile-1)*t_step
				t_tot = (nfile-1)*t_step
			endif
			
			if(t_step/=.0) then
				if(n_int==0) n_int = nfile/2
				if(2*(n_int/2)==n_int) n_int = n_int+1				!make it odd				
				if(j_verb==1) write(*,*) 'Time integration window set to nfile/2:',n_int

				f_max_plot = .5/t_step
				freq_step = 1./((n_int-1)*t_step)
				n_freq_max = (n_int-1)/2+1					!f_max_plot/freq_step			!also n_freq_max= .5*t_int/t_step
				write(*,*) 'MD time sequence:  t_total [ps]',t_tot,'t_step [ps]',t_step
				write(*,'("Maximum energy [THz]            ",f6.2,"      energy step [THz]",f6.2)') f_max_plot,freq_step
		
				t_width = .5*(n_int-1)*t_step		!effective width
				f_width = 1./t_width				!resolution= 4*f_max_plot/n_int
				
				write(*,'("Time-integration (BN) window FWHM [ps]",f8.2)') t_width
				write(*,'("Energy resolution FWHM [THz]            ",f8.2)') f_width
				write(*,*) 
			endif

				allocate(t_wind(n_int),SOURCE=0.0)
				allocate(om_phase(n_int),SOURCE=(.0,.0)) 

			at_pos_min = -n_row/2.
			at_pos_max = n_row/2.-1.
			at_pos_centre = .0
			twopi_nr = twopi/(1.0d0*n_row)
			 			
			
			j_fft = 1
			jfile = 1				!now to be used as index for the successive output files			
			j_plot = 0			!j_plot=0 is initialised and identifies 1st pass: we start with a total scattering map

			map_loop: do
				do
					if(j_plot==-6)	then
						do i=1,3
							write(*,*) xyz(i),' min, max',at_pos_min(i),at_pos_max(i),' new values: (,, = NO CHANGE)'
							read(*,*) at_pos_min(i),at_pos_max(i)
							if(at_pos_min(i)<-n_row(i)/2.) at_pos_min(i) = -n_row(i)/2.
							if(at_pos_max(i) > n_row(i)/2.) at_pos_max(i) = n_row(i)/2.
							write(*,*) xyz(i),' min, max',at_pos_min(i),at_pos_max(i)
							write(9,*) xyz(i),' min, max',at_pos_min(i),at_pos_max(i)
						enddo 
						nrow = int(at_pos_max(1)-at_pos_min(1)+1)
						at_pos_centre = ceiling(.5*(at_pos_max+at_pos_min))
					endif			

CC					write(*,*) 		'nrow,at_pos_centre',nrow,at_pos_centre		
					j_plot = 0			!j_plot=0 back to ZERO: we start again with a total scattering map

					write(*,*) 'Number of Brillouin zones (1 ... , 0=END), Q-plane (0=general, 1=(hk0), 2=(hhl))'		
					read(*,*) n_bz,j_dir
					if(n_bz==0.) exit map_loop	
					
C *** read Xray formfactor parameters if needed (j_weight=3, x_ffpar(1,1)=.0 )
					if((j_xray>=1.or.j_weight==3).and.x_ffpar(1,1)==.0) then
						open(4,file='xray_ff.txt',action='read',status ='old',iostat=ios)
						if(ios.ne.0) then
							open(4,file='/usr/local/mp_tools/ref/xray_ff.txt',action='read',status ='old',iostat=ios)
							if(ios.ne.0) then
								do
									write(*,*) 'File xray_ff.txt not found, type valid access path/filename'
									read(*,'(a)') x_file_name
									open(4,file=trim(x_file_name),action='read',status ='old',iostat=ios)
									if(ios==0) exit
									write(*,*) 'File',trim(x_file_name),' not found, try again ...'
								enddo
							endif
						endif
						atom_loop: do j=1,n_atom
							rewind(4)
							do i=1,210
								read(4,*) atom
								if(atom==trim(at_label(j))//trim(at_name_ext(j))) then
									backspace(4)
									read(4,*) atom,x_ffpar(j,1:9)
									cycle atom_loop
								endif
							enddo
							write(*,*) 'not found'
						enddo atom_loop
						j_xray = 1
						write(*,*)
						write(*,*) 'Xray form-factors'
						do j=1,n_atom
							write(*,*) at_label(j),x_ffpar(j,:)
						enddo			
					endif

					if(j_dir==1)then
						e1 = (/1,0,0/)
						e2 = (/0,1,0/)
						exit
					elseif(j_dir==2)then
						e1 = (/1,1,0/)
						e2 = (/0,0,1/)
						exit
					else						
						do
							write(*,*) 'Input perpendicular vectors to define the Q-plane e1(3),e2(3):'
							read(*,*) e1,e2
							if(dot_product(e1,e1).ne.0..and.dot_product(e2,e2).ne.0..and.dot_product(e1,e2)==0.) then
								exit
							else
								cycle
							endif
						enddo
						exit
					endif
				enddo

				ev = (/e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)/)
				e1_norm = sqrt(dot_product(1.*e1,1.*e1))
				e2_norm = sqrt(dot_product(1.*e2,1.*e2))
				ev_norm = sqrt(dot_product(1.*ev,1.*ev))

				tpe1 = twopi/dot_product(e1,anint(at_pos_max-at_pos_min+1)/e1_norm**2)
				tpe2 = twopi/dot_product(e2,anint(at_pos_max-at_pos_min+1)/e2_norm**2)
				tpev = twopi/dot_product(ev,(at_pos_max-at_pos_min+1)/ev_norm**2)
				n_qx = sum(abs(e1)*(at_pos_max-at_pos_min+1))*n_bz/e1_norm**2			    
				n_qy = sum(abs(e2)*(at_pos_max-at_pos_min+1))*n_bz/e2_norm**2

				n_qx8 = n_qx
				n_qy8 = n_qy
				nsuper8 = nsuper

				write(*,*) 'Vertical axis:', ev

				allocate(x_ffq(n_atom,n_qx,n_qy))

				c_min_save = .0
				c_max_save = .0
				
C **** cycle over q_center(3)-values
			
			bz_loop: do				
				write(*,*) 'Q_center [hkl] (9 9 9 END):'
				read(*,*) q_center

				if(q_center(1)==9.) exit bz_loop
				
				qx_min = dot_product(e1,q_center)/e1_norm**2-.5*n_bz
				qx_max = dot_product(e1,q_center)/e1_norm**2+.5*n_bz
				qy_min = dot_product(e2,q_center)/e2_norm**2-.5*n_bz
				qy_max = dot_product(e2,q_center)/e2_norm**2+.5*n_bz
				tpq_center = twopi_s*q_center
				
				q_v = (dot_product(ev,q_center)/ev_norm)*ev

C *** generate the Xray form-factor tables
			  if(j_weight==3) then
					do i=1,n_qx
						do j=1,n_qy
							q1 = (qx_min+(i-1)/n_qx)*(e1/a_par)
							q2 = (qy_min+(j-1)/n_qy)*(e2/a_par)
							q_sq = dot_product(q1,q1)+dot_product(q2,q2)
							do j_at=1,n_atom
								x_ffq(j_at,i,j) = x_ffpar(j_at,9)
								do ii=1,4
									x_ffq(j_at,i,j) = x_ffq(j_at,i,j)+x_ffpar(j_at,2*ii-1)*exp(-.25*x_ffpar(j_at,2*ii)*q_sq)
								enddo
							enddo
						enddo
					enddo
			 endif

C *** generate the Fourier time window table				
			if(n_int<=1) then
				t_wind = 1.
			else
				do i=1,n_int
CC					t_wind(i) = .5*(1.-cos(twopi*(i-1)/real(n_int-1)))												!max is at n_int/2+1		Hann window
					t_wind(i) = 0.355768-0.487396*cos(twopi*(i)/real(n_int+1))							!max is at n_int/2+1		Nuttall window
					t_wind(i) = t_wind(i)+0.144232*cos(2.*twopi*(i)/real(n_int+1))						
					t_wind(i) = t_wind(i)-0.012604*cos(3.*twopi*(i)/real(n_int+1))							
				enddo
			endif
			
        if(j_oneph==1) then
          write(*,*) 'Single-phonon FFT start'
        else
          if(j_fft==1) write(*,*) 'FINUFFT start'
          if(j_fft==0) write(*,*) 'Simple FT start'
        endif
				  
				t_sum = .0
				call cpu_time(t1)

				CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
				sc_r = 1./sc_r

				allocate(cs_atom(nfile,n_atom,n_qx,n_qy))
							 
				frame_loop: do ifile=1,nfile
				
				  if(j_oneph==1) then      !do one-phonon approximation
				    allocate(u_x(nrow,nrow),u_y(nrow,nrow),u_qx(nrow,nrow),u_qy(nrow,nrow))
				    nqx_min = nrow*(qx_min-floor(qx_min))
				    nqy_min = nrow*(qy_min-floor(qy_min))

  					do j = 1,n_atom
  					  u_x = (.0,.0)
  					  u_y = (.0,.0)

              if(j_dir==1)then											!(hk0) plane - vertical [0 0 1]
							ii = 0
              do i=1,nsuper
              	if(at_pos_file(1,i,j,ifile)/=.0.or.at_pos_file(2,i,j,ifile)/=.0.or.at_pos_file(3,i,j,ifile)/=.0) then
 									k = at_ind(1,i,j)
									l = at_ind(2,i,j)
									d_x = at_pos_file(1,i,j,ifile)-at_base(j,1)-at_ind(1,i,j)+n_row(1)/2
									d_y = at_pos_file(2,i,j,ifile)-at_base(j,2)-at_ind(2,i,j)+n_row(2)/2
									u_x(k,l) = u_x(k,l)+d_x			!*c_f     !staying in the plane
									u_y(k,l) = u_y(k,l)+d_y			!*c_f
 								endif
              enddo

              elseif(j_dir==2) then     ! the (h h l) plane

              do i=1,nsuper
              	if(at_pos_file(1,i,j,ifile)/=.0.or.at_pos_file(2,i,j,ifile)/=.0.or.at_pos_file(3,i,j,ifile)/=.0) then
                  k = modulo(at_ind(1,i,j)-1+at_ind(2,i,j)-1,nrow)+1
                  l = at_ind(3,i,j)
                  d_x = at_pos_file(1,i,j,ifile)-at_base(j,1)-at_ind(1,i,j)+n_row(1)/2+1
                  d_x = d_x+at_pos_file(2,i,j,ifile)-at_base(j,2)-at_ind(2,i,j)+n_row(2)/2+1
                  d_x = d_x/sqrt2
                  d_y = at_pos_file(3,i,j,ifile)-at_base(j,3)-at_ind(3,i,j)+n_row(3)/2+1
                  u_x(k,l) = u_x(k,l)+d_x   !staying in the plane
                  u_y(k,l) = u_y(k,l)+d_y
 								endif
              enddo

            endif		!j_dir				

              u_qx = fft(u_x,inv=.false.)/(1.*nrow)     !in fact 1/sqrt(nrow**2)
              u_qy = fft(u_y,inv=.false.)/(1.*nrow)     !FT in 0th Bz possibly with a q_z component

              do ii=1,n_bz      
                do jj=1,n_bz              
                   do l=1,nrow
                    do k=1,nrow
                      q_x = (qx_min+(ii-1)+(k-1)/(n_row(1)*1.))            !*(e1/a_par)
                      q_y = (qy_min+(jj-1)+(l-1)/(n_row(2)*1.))            !*(e2/a_par)
                      kk = modulo(nqx_min+k,nrow)
                      ll = modulo(nqy_min+l,nrow)                 
                      if(kk==0) kk = nrow                 
                      if(ll==0) ll = nrow 
                      if(j_dir==1) then               
                        phase_bz = cexp((0.,-1.)*twopi_s*(q_x*at_base(j,1)+q_y*at_base(j,2)))
                      elseif(j_dir==2) then               
                        phase_bz = cexp((0.,-1.)*twopi_s*(q_x*(at_base(j,1)+at_base(j,2))+q_y*at_base(j,3)))
                      endif
                      cs_atom(ifile,j,(ii-1)*nrow+k,(jj-1)*nrow+l) = phase_bz*(q_x*u_qx(kk,ll)+q_y*u_qy(kk,ll))
                    enddo
                  enddo
                enddo
              enddo
            enddo !j
            deallocate(u_x,u_y,u_qx,u_qy)

				  else                    !j_oneph/=1 do the complete NUFFT thing
 						do j = 1,n_atom
							allocate(xf(n_sup(j)),yf(n_sup(j)),cf(n_sup(j)))																	
							allocate(ampl_tot(n_qx*n_qy),ampl_tot_2d(n_qx,n_qy))
 
   						if(data_type=='cell'.or.data_type=='quick') then
								ii = 0
								do i=1,nsuper
									non_zero = (at_pos_file(1,i,j,ifile)/=.0.or.at_pos_file(2,i,j,ifile)/=.0.or.at_pos_file(3,i,j,ifile)/=.0)
									inside = .true.
									do iii=1,3
										inside = inside.and.((at_ind(iii,i,j)-n_row(iii)/2-1)>=at_pos_min(iii).and.
     1											(at_ind(iii,i,j)-n_row(iii)/2-1)<=at_pos_max(iii))
									enddo

									if(non_zero.and.inside) then
										ii = ii+1
										xf(ii) = dot_product(e1,(at_pos_file(1:3,i,j,ifile)-at_pos_centre))*tpe1
										yf(ii) = dot_product(e2,(at_pos_file(1:3,i,j,ifile)-at_pos_centre))*tpe2
										cf(ii) = cexp((0.,-1.)*dot_product(tpq_center,at_pos_file(1:3,i,j,ifile)))
									endif
								enddo              
							else					!bulk
								ii = 0
								do i=1,n_sup(j)
									inside = .true.
									do iii=1,3
										inside = inside.and.(at_pos_file_bulk(iii,ind_at(j)+i,ifile)>=at_pos_min(iii).and.
     1										at_pos_file_bulk(iii,ind_at(j)+i,ifile)<=at_pos_max(iii)+1)
									enddo
									if(inside)then
										ii = ii+1
										xf(ii) = dot_product(e1,(at_pos_file_bulk(1:3,ind_at(j)+i,ifile)-at_pos_centre))*tpe1
										yf(ii) = dot_product(e2,(at_pos_file_bulk(1:3,ind_at(j)+i,ifile)-at_pos_centre))*tpe2
										cf(ii) = cexp((0.,-1.)*dot_product(tpq_center,at_pos_file_bulk(1:3,ind_at(j)+i,ifile)))
               		endif
								enddo							
							endif		!data_type
              n_fft = ii			!we need n_fft 64bit
 						if(j_fft==1)then
             	call finufft2d1(n_fft,xf,yf,cf,iflag,eps,n_qx8,n_qy8,ampl_tot,nul_opt,ier)
            else
            	write(*,*) 'Normal FT can be very loooong!'
							call dirft2d1(n_fft,xf,yf,cf,iflag,n_qx8,n_qy8,ampl_tot)	
						endif							!direct FT for check
              ampl_tot_2d = reshape(source=ampl_tot,shape=[n_qx,n_qy])
              cs_atom(ifile,j,:,:) = ampl_tot_2d/sqrt(real(n_sup(j)))

 				 			deallocate(xf,yf,cf,ampl_tot_2d,ampl_tot)					!,om_phase(nfile))	
 						enddo  !j (n_atom)
         endif ! j_oneph
          if(mod(ifile,200)==0) write(*,*) 'FFT done frame',ifile
				enddo frame_loop  !ifile		!finish all file_related normalisation here



				call cpu_time(t2)
				CALL SYSTEM_CLOCK (COUNT = sc_c2)
				if(j_verb==1) write(*,*) 'FINUFFT on',nfile,'snapshots CPU_TIME',t2-t1,'  SYS_TIME',(sc_c2-sc_c1)*sc_r
					
C **** calculate the mean value over the whole trajectory, in fact the 0th component of timeFT = elastic scattering

				allocate(cs_mean(n_atom,n_qx,n_qy),source=(.0,.0))
CC				cs_mean = (.0,.0)
				do j=1,n_atom
					do ii=1,n_qx
						do jj=1,n_qy
							cs_mean(j,ii,jj) = sum(cs_atom(:,j,ii,jj))/real(nfile)
						enddo
					enddo
				enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC *** Frequency loop
CC						write(*,*) 'Choose a plot option (PS/TXT file output is ',trim(ps_out(j_ps+1)),'):'
CC							write(*,*) '       0  EXIT'
CC							write(*,*) '       1  explore total scattering (-1 edit atom masks)'
CC						if(t_step>.0) then			
CC							write(*,*) '       2  explore E=const maps         (-2 edit atom masks)'
CC							write(*,*) '       3  make a stack of E=const maps (-3 reset atom masks)'
CC							write(*,*) '       4  explore E(Q) maps            (-4 edit atom masks)'
CC							write(*,*) '       5  make a stack of E(Q) maps    (-5 reset atom masks)'
CC						endif
CC							write(*,*) '       6  change the HKL plane, range and centre (BZ)'
CC							write(*,*) '                                       (-6 change the real space range)'
CC							if(data_type=='cell'.or.data_type=='quick') then
CC								if(j_oneph==0) write(*,*) '       7  toggle the NU_FFT mode to ONE_PHONON (go on via 6)'
CC								if(j_oneph==1) write(*,*) '       7  toggle the ONE_PHONON mode to NU_FFT (go on via 6)'
CC							endif
CC							if(j_ps==0) write(*,*) '       8  toggle the PS/TXT output ON (mind the TXT switch in PAR)'
CC							if(j_ps==1) write(*,*) '       8  toggle the PS/TXT output OFF (mind the TXT switch in PAR)'
CC							write(*,*) '       9  change the time integration window width, weighting etc.'		!include here the straight FT option
CC							read(*,*) j_plot


CC			j_plot = 0									!j_plot=0 identifies 1st pass: we start with a total scattering map
			j_ps = 0										!no file output for the start
			dq_p = .0

			freq_loop: do								!goes almost till the end

				allocate(cs(nfile,n_qx,n_qy),SOURCE=(.0,.0))

				if(j_plot==0.or.j_plot==-3.or.j_plot==-5) then			!at the beginning of plotting (0) or upon request (-3,-5) reset masks
					at_mask = 1
					j_atc1 = 0
					j_atc2 = 0
					write(*,'("Atom masks reset:",50i3)')(at_mask(i),i=1,n_atom),j_atc1,j_atc2
				else
					write(*,'("Actual atom masks:",50i3)')(at_mask(i),i=1,n_atom),j_atc1,j_atc2
				endif

				if(j_plot==-1.or.j_plot==-2.or.j_plot==-3.or.j_plot==-5)	then	!edit atom masks
					write(*,*)'Type in atom masks:'
					do
						read(*,*)(at_mask(i),i=1,n_atom)
						if(all(at_mask(1:n_atom).ge.0).and.all(at_mask(1:n_atom).le.1)) then 
							j_atc1 = 0
							j_atc2 = 0	
							exit					
						elseif(at_mask(1).eq.2) then
							j_atc1 = at_mask(2)
							j_atc2 = at_mask(3)
							exit
						else
							write(*,*) 'Input out of range, repeat ...'
						endif
					enddo
				endif

				write(masks,'(50i1.1)') (at_mask(i),i=1,n_atom)			!to be used in plot titles
				
				if(abs(j_plot)<=1)then			!t_step=0. identifies total scattering
					j_disp = 0								!E = const map
					f_plot = 0.
					n_plot = 1
					i_step = 1
					n_freq_plot = 1
					n_freq_min  = n_freq_plot
					write(*,*) 'Total scattering (instant integration range)'

				elseif(abs(j_plot)==2) then
					write(*,*) 'E=const map E_plot [THz]:'
					read(*,*) f_map_min
						j_disp = 0									! E = const map
						n_freq_min  = nint(f_map_min/freq_step)+1
						n_plot =1

					elseif(abs(j_plot)==3) then
					write(*,*) 'E=const map stack:  E_min, E_step [THz], n_step (≤8 for PGPLOT screen):'
					read(*,*) f_map_min,f_map_step,n_plot
					do
						j_disp = 0									! E = const map
						j_scan = 0
						j_exit = 0
						n_freq_min  = nint(f_map_min/freq_step)+1
						i_step = nint(f_map_step/freq_step)
						if(i_step>=1) exit
						write(*, '("E_step too small, setting it to ",f5.2,"THz")') freq_step
						f_map_step = freq_step
					enddo

				elseif(abs(j_plot)==4.or.abs(j_plot)==5) then
					do
						write(*,*) 'E(Q) map: initial and final points Q1, Q2 [hkl]:'
						read(*,*) q1_3d,q2_3d					
						q1 =q1_3d-q_v
						q2 =q2_3d-q_v	
						if(dot_product(q1,ev).ne.0.)then
							write(*,*) 'Q1 not in the Q-plane:',q1
							cycle
						endif
						if(dot_product(q2,ev).ne.0.)then
							write(*,*) 'Q2 not in the Q-plane:',q2
							cycle
						endif
						exit
					enddo

					j_disp = 1									! E(Q) "dispersion" map
					n_freq_plot = n_freq_max
					n_freq_min  = 1
					n_plot = 1
					i_step = 1
					dq_p = .0
					

					if(abs(j_plot)==5) then
						write(*,*) 'perpendicular Q_step length [rlu]:'
						read(*,*) q_step
						n_plot = 7
						dq = q2-q1
						dq_p = (/dq(2)*ev(3)-dq(3)*ev(2),dq(3)*ev(1)-dq(1)*ev(3),dq(1)*ev(2)-dq(2)*ev(1)/) !dq_p is perpendicular in plane to dq
						dq_p = -1.*dq_p/sqrt(dot_product(dq_p,dq_p))													
						dq_p = q_step*dq_p													
					endif
				endif

				allocate(map_unit(n_plot),SOURCE=0)

C **** plot loop

			plot_loop: do i_plot=1,n_plot,i_step		!executed just once for j_plot= 1,2,4
			
				if(j_verb==1.and.(j_plot==3.or.j_plot==5)) write(*,*) 'Plot no.:',i_plot

C **** now AT_MASK are (re)set, we can make the right replica of CS_ATOM
				cs = (.0,.0)

				if(abs(j_plot)<=1) then
!$omp parallel shared(cs,cs_atom,at_mask,nfile,n_atom) 
!$omp do
					do ifile=1,nfile				!total scattering
						do j=1,n_atom
							if(at_mask(j)==1)then					
								if(j_weight==1)then
									cs(ifile,:,:) = cs(ifile,:,:)+cs_atom(ifile,j,:,:)
								elseif(j_weight==2)then
									cs(ifile,:,:) = cs(ifile,:,:)+b_coh(j)*cs_atom(ifile,j,:,:)
								elseif(j_weight==3)then
									cs(ifile,:,:) = cs(ifile,:,:)+x_ffq(j,:,:)*cs_atom(ifile,j,:,:)
								endif
							endif
						enddo
					enddo
!$omp end do
!$omp end parallel

				elseif(abs(j_plot)>1) then
!$omp parallel shared(cs,cs_atom,at_mask,nfile,n_atom) 
!$omp do
					do ifile=1,nfile
						do j=1,n_atom
							if(at_mask(j)==1)then
								if(j_weight==1)then
									cs(ifile,:,:) = cs(ifile,:,:)+(cs_atom(ifile,j,:,:)-cs_mean(j,:,:))
								elseif(j_weight==2)then
									cs(ifile,:,:) = cs(ifile,:,:)+b_coh(j)*(cs_atom(ifile,j,:,:)-cs_mean(j,:,:))
								elseif(j_weight==3)then
									cs(ifile,:,:) = cs(ifile,:,:)+x_ffq(j,:,:)*(cs_atom(ifile,j,:,:)-cs_mean(j,:,:))
								endif
							endif
						enddo
					enddo
!$omp end do
!$omp end parallel
				endif
CC				write(*,*) 'cs',cs(:,5,5)

C **** now generate the plot data according to the particular case specs

			if(abs(j_plot)<=1) then															!total scattering
				allocate(cs_plot(n_qx,n_qy),cs_out(n_qx,n_qy), SOURCE=0.0)					
				do ifile=1,nfile
					cs_plot = cs_plot+cs(ifile,:,:)*conjg(cs(ifile,:,:))
				enddo
				cs_plot = cs_plot/real(nfile)
			endif																						!total scattering

			if(abs(j_plot)>=2) then													!E-resolved cases (≈ 130 lines)
				if(j_disp>0) then
					q1_x = dot_product(e1,q1+(i_plot-n_plot/2-1)*dq_p)/e1_norm**2
					q1_y = dot_product(e2,q1+(i_plot-n_plot/2-1)*dq_p)/e2_norm**2				
					q2_x = dot_product(e1,q2+(i_plot-n_plot/2-1)*dq_p)/e1_norm**2
					q2_y = dot_product(e2,q2+(i_plot-n_plot/2-1)*dq_p)/e2_norm**2
					if(q1_x.lt.qx_min.or.q1_x.gt.qx_max.or.q1_y.lt.qy_min.or.q1_y.gt.qy_max)then
						write(*,*) 'Q1 out of the map range:',q1_x,q1_y
						cycle
					endif
					if(q2_x.lt.qx_min.or.q2_x.gt.qx_max.or.q2_y.lt.qy_min.or.q2_y.gt.qy_max)then
						write(*,*) 'Q2 out of the map range:',q2_x,q2_y
						cycle
					endif
								
					n_q1x = abs(q1_x-qx_min)*n_row(1)*e1_norm+1			!find indices of the first and last Q_points of the dispersion plot
					n_q2x = abs(q2_x-qx_min)*n_row(1)*e1_norm+1
					n_q1y = abs(q1_y-qy_min)*n_row(2)*e2_norm+1
					n_q2y = abs(q2_y-qy_min)*n_row(2)*e2_norm+1

					if(abs(n_q2x-n_q1x).ge.abs(n_q2y-n_q1y))then
						n_qdisp = abs(n_q2x-n_q1x)
						q_min = q1_x
						q_max = q2_x
						ax_label = 'Q_x [rlu]'
CC						if(j_verb==1) write(*,*) ax_label,n_q1x,n_q2x,n_qdisp
					else
						n_qdisp = abs(n_q2y-n_q1y)
						q_min = q1_y
						q_max = q2_y
						ax_label = 'Q_y [rlu]'
CC						if(j_verb==1) write(*,*) ax_label,n_q1y,n_q2y,n_qdisp
					endif
					q_nx = (n_q2x-n_q1x)/(1.*n_qdisp)			
					q_ny = (n_q2y-n_q1y)/(1.*n_qdisp)			
CC					write(*,*) 'n_qdisp,q_nx,q_ny 1',n_qdisp,q_nx,q_ny,q_min,q_max

				endif !if(j_disp>0)

				CALL SYSTEM_CLOCK (COUNT = sc_c1)
				call cpu_time(t1)

				if(j_disp==0) then	
					f_plot = f_map_min+(i_plot-1)*freq_step
					j_freq = n_freq_min+i_plot-1										!for map
					allocate(cs_plot(n_qx,n_qy),cs_out(n_qx,n_qy),cs3(nfile-n_int+1))					
					cs_plot = .0
						do i=1,n_int
							om_phase(i) = cexp((0.,-1.)*twopi_s*(j_freq-1)*(i-n_int/2-1)/real(n_int))		!n_int defines the effective "measurement" length
						enddo
						om_phase = om_phase*t_wind																					!window will be included in om_phase

!$omp parallel shared(cs_plot,om_phase,nfile,n_int,n_qx,n_qy) private (cs3,cs_p1,cs_p2) 
!$omp do
						do ii=1,n_qx
							do jj=1,n_qy
								cs3 = .0
								if(j_atc1== 0) then
									cs_p1 => cs(n_int/2+1:nfile-n_int/2,ii,jj)					!point to a subsection with negative shift in time
									do i=1,n_int
										cs_p2 => cs(i:nfile-n_int+i,ii,jj)								!point to a subsection with positive shift in time
										cs3 = cs3+om_phase(i)*cs_p2												!Time window to reduce cut-off effects is contained in the phase factor
									enddo									
									cs3 = cs3*(conjg(cs_p1))
								else
									cs_p1 => cs_atom(n_int/2+1:nfile-n_int/2,j_atc1,ii,jj)					!point to a subsection with negative shift in time
									do i=1,n_int
										cs_p2 => cs_atom(i:nfile-n_int+i,j_atc2,ii,jj)								!point to a subsection with positive shift in time
										cs3 = cs3+om_phase(i)*cs_p2												!Time window to reduce cut-off effects is contained in the phase factor
									enddo									
									cs3 = cs3*(conjg(cs_p1))								
								endif											

								do i=1,nfile-n_int+1
									cs_plot(ii,jj) = cs_plot(ii,jj)+cs3(i)
								enddo
							enddo												
						enddo																		
!$omp end do
!$omp end parallel

				deallocate(cs3)
				
				
				elseif(j_disp==1) then							!for dispersion
					allocate(cs_plot(n_freq_plot,n_qdisp),cs_out(n_freq_plot,n_qdisp),qq(n_qdisp),ff(n_freq_plot))
					allocate(cs2(n_int),cs3(n_int),cs4(n_int,n_qdisp))

					cs4 = (.0,.0)
!$omp parallel shared(cs,cs4,t_wind,nfile,n_int,n_qdisp,n_qx,n_qy,n_q1x,n_q1y,q_nx,q_ny) private (cs1,cs2,cs3) 
!$omp do
					do ifile = 1,nfile-n_int+1
						do ii = 1,n_qdisp
							if(j_atc1==0) then
								cs1 = cs(ifile+n_int/2,n_q1x+nint((ii-1)*q_nx),n_q1y+nint((ii-1)*q_ny))	!here pointers would be difficult, but the arrays are smaller
								do i=1,n_int																																		!our FT center is on n/2+1
										cs2(i) = cs(ifile+i-1,n_q1x+nint((ii-1)*q_nx),n_q1y+nint((ii-1)*q_ny))*t_wind(i)
								enddo		
								cs3 = fft(cs2,inv=.false.)/sqrt(1.*n_int)
							else
								cs1 =cs_atom(ifile+n_int/2,j_atc1,n_q1x+nint((ii-1)*q_nx),n_q1y+nint((ii-1)*q_ny))	!here pointers would be difficult, but the arrays are smaller
								do i=1,n_int																																		!our FT center is on n/2+1
									cs2(i) =cs_atom(ifile+i-1,j_atc2,n_q1x+nint((ii-1)*q_nx),n_q1y+nint((ii-1)*q_ny))*t_wind(i)
								enddo
								cs3 = fft(cs2,inv=.false.)/sqrt(1.*n_int)
							endif	
								cs4(:,ii) = cs4(:,ii)+cs3*conjg(cs1)			
						enddo
					enddo
!$omp end do
!$omp end parallel
					cs_plot(:,:) = abs(cs4(1:n_freq_plot,:))
					deallocate(cs2,cs3,cs4)
				endif		!j_disp

				CALL SYSTEM_CLOCK (COUNT = sc_c2)
				call cpu_time(t2)
				if(j_verb==1.and.(j_plot==2.or.j_plot==4))
     1			write(*,*) 'Time FT (normal, optimised & OMP):  PROC time ', t2-t1,' SYS time',(sc_c2-sc_c1)*sc_r

				n_frame = nfile-n_int+1
				
 	     	if(j_disp.ge.1) cs_plot = transpose(cs_plot)			!from now on the shape is cs_plot(n_qdisp,n_freq_plot)
				cs_plot = cs_plot/real(.5*n_frame*n_int)					!.5 comes for the integral of the time window profile

			endif !(abs(j_plot)>1)													!E-resolved cases

C *** apply the speckle filter
				if(j_plot>0.and.s_trig>0.) then				!the speckle filter can be turned off effectively by setting s_trig=0
					if(j_verb==1) write(*,*) 'Speckle filter s_trig:',s_trig
					if(j_disp==0) then
						nq_tot = n_qx*n_qy
					elseif(j_disp==1) then
						nq_tot = n_qdisp*n_freq_plot					
					endif
					allocate(csp(nq_tot))
					if(j_disp==1) cs_plot = transpose(cs_plot)
					csp = reshape(source=cs_plot,shape=[nq_tot])
					do i2=1,nq_tot
						i1 = i2-1
						if(i1<1) i1 = i1+nq_tot
						i3 = i2+1
						if(i3>nq_tot) i3 = i3-nq_tot
						if(csp(i2)<min(csp(i1),csp(i3))/s_trig) csp(i2)=min(csp(i1),csp(i3))
						if(csp(i2)>s_trig*max(csp(i1),csp(i3))) csp(i2)=max(csp(i1),csp(i3))
					enddo
					if(j_disp==0) then
						cs_plot = reshape(source=csp,shape=[n_qx,n_qy])
					elseif(j_disp==1) then
						cs_plot = reshape(source=csp,shape=[n_freq_plot,n_qdisp])
					endif

					if(j_disp==0) then				!2nd round perpendicular
						cs_plot = transpose(cs_plot)										!repeat filter in the perpendicular direction
						csp = reshape(source=cs_plot,shape=[nq_tot])
						do i2=1,nq_tot
							i1 = i2-1
							if(i1<1) i1 = i1+nq_tot
							i3 = i2+1
							if(i3>nq_tot) i3 = i3-nq_tot
							if(csp(i2)<min(csp(i1),csp(i3))/s_trig) csp(i2)=min(csp(i1),csp(i3))
							if(csp(i2)>s_trig*max(csp(i1),csp(i3))) csp(i2)=max(csp(i1),csp(i3))
						enddo
						cs_plot = reshape(source=csp,shape=[n_qy,n_qx])
					endif
					cs_plot = transpose(cs_plot)										
					deallocate(csp)
				endif				!s_trig>0.

C *** make a copy of CS_PLOT for .TXT output and linear plots
				cs_out = cs_plot										

				if(j_logscale==1) cs_plot = log10(cs_plot)	
				
C **** Draw the plot  
				f_min = 0.
				if(f_max==.0) f_max = f_max_plot
				
				if(c_min_save==0..and.c_max_save==0.) then
					c_min = anint(minval(cs_plot)-1.)
					c_max = anint(maxval(cs_plot)+1.)
				endif
				if(j_logscale==1.and.c_min<-8.) c_min = -8.
				write(*,*) 'c_min,c_max',c_min,c_max
			
C *** PGPLOT: set the coordinate transformation matrix: world coordinate = pixel number.
				if(j_disp==0) then
					tr_shift_x = .5
					tr_shift_y = .5
					if(n_qx==2*(n_qx/2)) tr_shift_x = 1.
					if(n_qy==2*(n_qy/2)) tr_shift_y = 1.
					TR(2) = real(n_bz)/n_qx
					TR(1) = qx_min-tr_shift_x*TR(2)		
					TR(3) = 0.0
					TR(5) = 0.0
					TR(6) = real(n_bz)/n_qy
					TR(4) = qy_min-tr_shift_y*TR(6)
				elseif(j_disp==1)then
					if(n_qdisp==2*(n_qdisp/2)) tr_shift_x = 1.
					TR(2) = (q_max-q_min)/(n_qdisp)
					TR(1) = q_min-tr_shift_x*TR(2)
					TR(3) = 0.0
					TR(5) = 0.0
					TR(6) = f_max/(n_freq_plot-1.)
					TR(4) = -TR(6)
				endif

	
				scale_loop: do
						if(map_unit(i_plot).le.0)then    !open the PGPLOT window, 1,1 no of panes
							map_unit(i_plot) = PGOPEN('/xserv')   !open the PGPLOT window, 1,1 no of panes
							CALL PGQCIR(C1, C2)
							NC = MAX(0, C2-C1+1)
							CALL PGASK(.FALSE.)     ! would not ask for <RET>
							BRIGHT = 0.5
							CONTRA  = 1.0
							CALL PALETT(6, CONTRA, BRIGHT)
							CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
							CALL PGSCRN(1, 'black', IER)
							CALL PGSCH(1.)					!set character height					
							CALL PGSLS (1)  !FULL
							CALL PGSLW(2)
							CALL PGSCI (1)  !black(0 = white)
						else
							CALL PGSLCT(map_unit(i_plot))
						endif

					plot_title1 = trim(subst_name)//'_'//trim(masks)
						if(j_oneph<=0)then
							mode = 'NU_FFT'
						else
							mode = 'Single_ph'
						endif

					
					if(j_disp==0) then
						if(abs(j_plot)<=1) then
							write(plot_title2,113) q_center,trim(at_weight_scheme(j_weight)),trim(mode)
						else 
							write(plot_title2,114) q_center,f_plot,trim(at_weight_scheme(j_weight)),trim(mode)
						endif
CC112   			format('  Q = [',3f6.2,']   Elastic scattering  ',a,'  ',a)   
113   			format('  Q = [',3f6.2,']   Total scattering  ',a,'  ',a)   
114   			format('  Q = [',3f6.2,']   Energy =',f6.3,' [THz]  ',a,'  ',a)   
						write(x_title,131)e1
						write(y_title,131)e2
131					format('[',i1,2i2,']',' [r.l.u.]')
						call PGSLCT(map_unit(i_plot))
						CALL PGPAP(p_size,e2_norm/e1_norm)     ! define the plot area as a square
						CALL PGENV(qx_min,qx_max,qy_min,qy_max,0,2) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
						CALL PGLAB(trim(x_title),trim(y_title),' ')  					!put the axis labels
						CALL PGMTXT ('T', 3., .5, .5, trim(plot_title1))			!put plot title on 2 lines
						CALL PGSCH(.7)					!set character height					
						CALL PGMTXT ('T', 1., .5, .5, trim(plot_title2))
						CALL PGSCH(1.)					!set character height					
						CALL PGIMAG(cs_plot(1:n_qx,1:n_qy),n_qx,n_qy,1,n_qx,1,n_qy,c_min,c_max,TR)
					elseif(j_disp==1) then
						write(plot_title2,115) q1_3d+(i_plot-n_plot/2-1)*dq_p,q2_3d+(i_plot-n_plot/2-1)*dq_p,
     1						trim(at_weight_scheme(j_weight)),trim(mode)
115   			format('  Q_1 = [',3f6.2,']','  Q_2 = [',3f6.2,'] ',a,' ',a)   
						call PGSLCT(map_unit(i_plot))
						CALL PGPAP(p_size,1.)     ! define the plot area as a rectangle
						CALL PGENV(q_min,q_max,f_min,f_max,0,2) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
						CALL PGLAB(trim(ax_label),'f [THz]',' ')  !put the axis labels
						CALL PGMTXT ('T', 3., .5, .5, trim(plot_title1))			!put plot title on 2 lines
						CALL PGSCH(.7)					!set character height					
						CALL PGMTXT ('T', 1., .5, .5, trim(plot_title2))
						CALL PGSCH(1.)					!set character height					
						CALL PGIMAG(cs_plot(1:n_qdisp,1:n_freq_plot),n_qdisp,n_freq_plot,1,n_qdisp,1,n_freq_plot,c_min,c_max,TR)
					endif
				
					if(j_grid==1) then
						XOPT = 'G'				!draw grid lines - after the color map, otherwise obscured
						YOPT = 'G'
						XTICK = 1.0
						YTICK = 1.0
						NXSUB = 5
						NYSUB = 5
						CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
					endif
						c_min_save = c_min
						c_max_save = c_max

					if(abs(j_plot)==3.or.abs(j_plot)==5)	exit scale_loop	!the stacks go straight out of the scale_loop

						if(abs(j_plot)<=1) then
							write(*,'("Modify intensity scale: min,max (0 0 or 9 9 EXIT)",2f8.1)') c_min,c_max
						elseif(abs(j_plot)==2) then
							write(*,'("Modify intensity scale: min,max (0 0 change E_plot, 9 9 EXIT)",2f8.1)') c_min,c_max
						elseif(abs(j_plot)==4) then
							write(*,'("Modify intensity scale: min,max (0 0 change Q1,Q2, 9 9 EXIT)",2f8.1)') c_min,c_max
						endif

						read(*,*)c_min,c_max								!,j_grid
						if(c_min==c_max) exit scale_loop

					enddo scale_loop
					
					if(c_min==c_max) then
CC						j_ps = 0						!j_ps is set in the options list under 8
						j_scan = 0
						j_exit = 0
CC						if(c_min==1) j_ps=1
						if(c_min==2) j_scan=1			!go for linear scans
						if(c_min==9.or.(c_min==0.and.abs(j_plot)<=1)) j_exit=1		!will ask for new options at the end of the frequency_loop
					endif

					c_min = c_min_save
					c_max = c_max_save
					
C **** make an optional hardcopy and .txt output
				if (j_ps.eq.1) then				!j_ps				
					jfile = 1
					do						!look for existing .txt files to continue numbering
			if(j_name==0.and.t_single)then
						write(file_res,107) trim(file_master),jfile
						write(file_ps,116) trim(file_master),jfile	
			else
						write(file_res,108) trim(file_master),nfile_min,nfile,jfile
						write(file_ps,117) trim(file_master),nfile_min,nfile,jfile	
			endif
107   			format(a,'_sq_',i2.2,'.txt')
108   			format(a,'_sq_',i4.4,'_',i4.4,'_',i2.2,'.txt')
116   			format(a,'_sq_',i2.2,'.ps')
117   			format(a,'_sq_',i4.4,'_',i4.4,'_',i2.2,'.ps')
						inquire(file=file_res,exist=found_txt)
						inquire(file=file_ps,exist=found_ps)

						if(.not.found_txt.and.(.not.found_ps)) exit				
						jfile = jfile+1
						if(jfile==100) then
							write(*,*)'Tidy up .txt/.ps files to restart count from 01 and type [RET]'
							read(*,*)
							jfile = 1
						endif							
					enddo

					write(9,*) 
					write(9,*) 'Input files:  ',trim(file_dat_t0),' to ',trim(file_dat)
					write(9,*) '  t0,t_step,n_row(3),n_atom:',t0,t_step,n_row,n_atom
					write(9,*) '  temperature [K]',temp
        	write(9,121) (at_name_par(j),j=1,n_atom)
121				format(' Atoms:        ',20(a,3x))
        	write(9,122) (conc(j),j=1,n_atom)
122				format(' Occupations:',20f7.3)
        	write(9,123) trim(at_weight_scheme(j_weight)),(b_coh(j),j=1,n_atom)
123				format(a,' b_coh:      ',20f7.3)
        	write(9,125) (at_mask(j),j=1,n_atom)
125				format(' Atom masks: ',i3,19i7)
        	write(9,126) j_atc1,j_atc2
126				format(' Atom pair correlation: ',2i3)
					write(9,*) 'n_bz:',n_bz
					write(9,124) e1,e2
124				format(' Q-plane: e1 = [',3i2,'], e2 = [',3i2,']')
					if(j_disp==0) then
						write(9,*) 'Q center, Energy transfer [THz]:',q_center,f_plot
						if(abs(j_plot)<=1) then
							write(9,*) 'Total scattering'
							write(9,*) 'Energy resolution [THz]     INF'
						else
							write(9,*) 'Time-integration (BN window FWHM) [ps]',t_width
							write(9,*) 'Energy resolution [THz]',f_width
						endif
					elseif(j_disp==1) then
						write(9,*) 'Dispersion range [rlu]:',q1,'  ',q2
						write(9,*) 'Time-integration (BN window FWHM) [ps]',t_width
						write(9,*) 'Energy resolution [THz]',f_width
					endif
					if(j_txt==0)then
						write(9,*) 'Output file:   ',trim(file_ps)
					else
						write(9,*) 'Output files:  ',trim(file_ps),'   ',trim(file_res)
					endif
					write(9,*)					
					
				
C **** Output the intensity map to a text file (linear scale)			
					if(j_txt==1) then			
						write(*,*) file_res
						open (3,file=file_res)																		!open the output file
						write(3,*) '*****    MP_SQ49: total scattering cross-section     *****'
						write(3,*) 
						write(3,*) 'Input files:  ',trim(file_dat_t0),' to ',trim(file_dat)
						write(3,*) '  t0,t_step,n_row(3),n_atom:',t0,t_step,n_row,n_atom
						write(3,*) '  temperature [K]',temp
						write(3,121) (at_name_par(j),j=1,n_atom)
						write(3,122) (conc(j),j=1,n_atom)
        		write(3,123) trim(at_weight_scheme(j_weight)),(b_coh(j),j=1,n_atom)
        		write(3,125) (at_mask(j),j=1,n_atom)
        		write(3,126) j_atc1,j_atc2
						write(3,124) e1,e2
						write(3,*) 'Time-integration (BN window FWHM) [ps]:',t_width

						if(j_disp==0) then
							if(abs(j_plot)==1) then
								write(3,*) 'Total scattering (infinite energy width)', 0.,9999.
							else
								write(3,*) 'Energy transfer & resolution [THz]:',f_plot,f_width
							endif
							write(3,*) 'Q center & map range (B_zones):',q_center,n_bz
							write(3,*) 'Plot size X (lines):  ',n_qx,qx_min,qx_min+(n_qx-1.)*(qx_max-qx_min)/(n_qx)
							write(3,*) 'Plot size Y (columns):',n_qy,qy_min,qy_min+(n_qy-1.)*(qy_max-qy_min)/(n_qy)
							write(3,*)
							do i=1,n_qx
								write(3,109) (cs_out(i,j),j=1,n_qy)
							enddo
						elseif(j_disp==1) then
							write(3,*) 'Energy resolution [THz]',f_width
							write(3,*) 'Dispersion range [rlu]:',q1,'  ',q2
							if((n_q2x-n_q1x).ge.(n_q2y-n_q1y)) then
								write(3,*) 'Plot size X (lines, Q):  ',n_qdisp,q1_x,q1_x+(n_qdisp-1.)*(q2_x-q1_x)/(n_qdisp)
							else
								write(3,*) 'Plot size X (lines, Q):  ',n_qdisp,q1_y,q1_y+(n_qdisp-1.)*(q2_y-q1_y)/(n_qdisp)
							endif
							write(3,*) 	 'Plot size Y (columns, f):',n_freq_plot,f_min,f_max
							write(3,*)
							do i=1,n_qdisp
								write(3,109) (cs_out(i,j),j=1,n_freq_plot)
							enddo
						endif
						close(3)
					endif		!j_txt
109   		format(200(1x,5e10.3))			

C **** Prepare and plot the same on .PS		
					write(*,*) file_ps
					IF (PGOPEN(file_ps//'/VCPS').LE.0) STOP
					CALL PGASK(.FALSE.)     ! would not ask for <RET>
					if(j_disp==0)then
						CALL PGPAP(7.0,e2_norm/e1_norm)     ! define the plot area as a square
					else
						CALL PGPAP(7.0,1./e1_norm)     ! define the plot area as a rectangle
					endif
					CALL PGQCIR(C1, C2)
					NC = MAX(0, C2-C1+1)
					CALL PALETT(6, CONTRA, BRIGHT)

					plot_title = trim(plot_title1)//('  '//trim(file_ps))	
       
					if(j_disp==0) then
						CALL PGENV(qx_min,qx_max,qy_min,qy_max,0,2) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
						CALL PGLAB(trim(x_title),trim(y_title),' ')  !put the axis labels
						CALL PGSCH(.7)					!set character height					
						CALL PGMTXT ('T', 3., .5, .5, trim(plot_title))			!put plot title on 2 lines
						CALL PGMTXT ('T', 1., .5, .5, trim(plot_title2))
						CALL PGSCH(1.0)					!set character height					
						CALL PGIMAG(cs_plot(1:n_qx,1:n_qy),n_qx,n_qy,1,n_qx,1,n_qy,c_min,c_max,TR)
					elseif(j_disp==1) then
						CALL PGENV(q_min,q_max,f_min,f_max,0,2) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
						CALL PGLAB(trim(ax_label),'f [THz]',' ')  !put the axis labels
						CALL PGSCH(.7)					!set character height					
						CALL PGMTXT ('T', 3., .5, .5, trim(plot_title))			!put plot title on 2 lines
						CALL PGMTXT ('T', 1., .5, .5, trim(plot_title2))
						CALL PGSCH(1.)					!set character height					
						CALL PGIMAG(cs_plot(1:n_qdisp,1:n_freq_plot),n_qdisp,n_freq_plot,1,n_qdisp,1,n_freq_plot,c_min,c_max,TR)
					endif

					if(j_grid==1) then
						XOPT = 'G'				!draw grid lines - after the color map, otherwise obscured
						YOPT = 'G'
						XTICK = 1.0
						YTICK = 1.0
						NXSUB = 5
						NYSUB = 5
						CALL PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
					endif
					CALL PGCLOS
				endif		!j_ps = 1
							
C
C **** generate linear scan profiles
C				
				if(j_scan.ge.1)then
					scan_loop: do
						if(plot_unit==0) then
							plot_unit = PGOPEN('/xserv')
							write(*,*) 'plot_unit',plot_unit
							CALL PGASK(.FALSE.)     ! would not ask for <RET>
							CALL PGPAP(7.0,1.5)     ! define the plot areaCC						CALL PGERAS
							call PGSUBP(1,4)
							CALL PGSCRN(0, 'white', IER)	!sets the color index of background to WHITE
							CALL PGSCRN(1, 'black', IER)
						else
							call PGSLCT(plot_unit)
						endif

						page_loop:do ii=1,4
							write(*,*) 'Q_plot, E_plot [THz]? (0 0 = END; Q = -9: E = const; Q=const: E = -9 ) '
							read(*,*) qq_plot,ff_plot

C
C **** check the limits
C
							if(qq_plot==.0.and.ff_plot==0.) then
								exit scan_loop																							! EXIT linear plots
							endif
							
							if(qq_plot.eq.-9.and.(ff_plot.lt.0..or.ff_plot.gt.f_max))then !E_const: E out of range
								cycle page_loop	
							endif	

							if(ff_plot.eq.-9.and.(qq_plot.lt.min(q_min,q_max).or.qq_plot.gt.max(q_min,q_max))) then
							 cycle page_loop 																							!Q=const: Q out of range
							endif

C
C **** do the plots
C
						if(qq_plot.eq.-9.)then														!do E_const over full Q_range
								j_scan = 2				
								do j=1,n_qdisp
									qq(j) = q_min+(j-1)*(q_max-q_min)/real(n_qdisp)
								enddo
								j_freq = nint(ff_plot/freq_step)+1
								if(ii==1)then
									c_min = .0
									c_max = maxval(cs_out(1:n_qdisp,j_freq))
									if(c_max.le.1.) then
										c_max = 1.
									else
										c_max = 10.*(int(c_max/10.)+1)
									endif
								else
									c_min = c_min_save
									c_max = c_max_save
								endif
							write(scan_title,105) trim(plot_title1),ff_plot
105       		format(a,'   E =',f5.2,' [THz]')   
							do
								CALL PGSCH(2.)					!set character height					
								CALL PGSCI (1)  !white
								CALL PGSLS (1)  !full
								CALL PGENV(qq(1),qq(n_qdisp),c_min,c_max,0,1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
								CALL PGLAB(trim(ax_label), 'cts',trim(scan_title))  !put the axis labels
								CALL PGSCI (ii+1)  !red-green-blue
								CALL PGLINE(n_qdisp,qq,cs_out(1:n_qdisp,j_freq))  !plots the curve
								c_min_save = c_min
								c_max_save = c_max
								write(*,*) 'ii',ii
								if(ii.eq.1)then
									write(*,*) 'No c_min/c_max adjustment in 1st panel'
									c_min = .0
									c_max = .0
								else
									write(*,*) 'c_min, c_max',c_min, c_max,' better values? (0 0 = END)'
									read(*,*) c_min, c_max
								endif
								if(c_min==c_max.and.c_min==.0) then
									c_min = c_min_save
									c_max = c_max_save
									cycle page_loop
								endif
								CALL PGERAS
								CALL PGPANL(1,ii-1)
							enddo
							
						elseif(ff_plot.eq.-9.)then																		!do Q_const over E_range < f_max
							j_q = nint((qq_plot-q_min)/(q_max-q_min)*n_qdisp)+1
							ff_min = 0.
							ff_max = f_max
							j_freq = nint(ff_max/freq_step)+1
							do j=1,n_freq_plot
								ff(j) = ff_min+(j-1)*(ff_max-ff_min)/real(n_freq_plot-1)
							enddo
							if(ii==1)then
								c_min = .0
								c_max = 10.
							else
								c_min = c_min_save
								c_max = c_max_save
							endif
							write(scan_title,106) trim(plot_title1),qq_plot
106       		format(a,'   q = ',f5.2,' [r.l.u.]')   
							do
								CALL PGSCH(2.)					!set character height					
								CALL PGSCI (1)  !white
								CALL PGSLS (1)  !full
								CALL PGENV(ff_min,ff_max,c_min,c_max,0,1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
								CALL PGLAB('E [THz]', 'cts',trim(scan_title))  !put the axis labels
								CALL PGSCI (ii+1)  !red-green-blue
								CALL PGLINE(n_freq_plot,ff,cs_out(j_q,1:j_freq))  !plots the curve
								c_min_save = c_min
								c_max_save = c_max
								if(ii.eq.1)then
									write(*,*) 'No c_min/c_max adjustment in 1st panel'
									c_min = .0
									c_max = .0
								else
									write(*,*) 'c_min, c_max',c_min, c_max,' better values? (0 0 = END)'
									read(*,*) c_min, c_max
								endif
								if(c_min==c_max.and.c_min==.0) then
									c_min = c_min_save
									c_max = c_max_save
									cycle page_loop
								endif
								CALL PGERAS
								CALL PGPANL(1,ii-1)
							enddo
						endif
					enddo	page_loop
					write(*,*) 'Make a copy of the graphics before you close it [RET]'
					read(*,*)
				enddo scan_loop				
			endif		!j_scan = 1

C
C **** deallocate and cycle loops
C				
				deallocate(cs_plot,cs_out)				!,cross_sec
				if(abs(j_plot)>1.and.j_disp.ge.1) deallocate(qq,ff)

			enddo plot_loop
			CALL PGEND
			deallocate (map_unit)
			deallocate(cs)			!it can reallocate on entering the loop with the new nfile from 9

				do
					if(j_plot==0.or.abs(j_plot)==3.or.abs(j_plot)==5.or.j_exit==1) then
						write(*,*) 'Choose a plot option (PS/TXT file output is ',trim(ps_out(j_ps+1)),'):'
							write(*,*) '       0  EXIT'
							write(*,*) '       1  explore total scattering (-1 edit atom masks)'
						if(t_step>.0) then			
							write(*,*) '       2  explore E=const maps         (-2 edit atom masks)'
							write(*,*) '       3  make a stack of E=const maps (-3 reset atom masks)'
							write(*,*) '       4  explore E(Q) maps            (-4 edit atom masks)'
							write(*,*) '       5  make a stack of E(Q) maps    (-5 reset atom masks)'
						endif
							write(*,*) '       6  change the HKL plane, range and centre (BZ)'
							write(*,*) '                                       (-6 change the real space range)'
							if(data_type=='cell'.or.data_type=='quick') then
								if(j_oneph==0) write(*,*) '       7  toggle the NU_FFT mode to ONE_PHONON (go on via 6)'
								if(j_oneph==1) write(*,*) '       7  toggle the ONE_PHONON mode to NU_FFT (go on via 6)'
							endif
							if(j_ps==0) write(*,*) '       8  toggle the PS/TXT output ON (mind the TXT switch in PAR)'
							if(j_ps==1) write(*,*) '       8  toggle the PS/TXT output OFF (mind the TXT switch in PAR)'
							write(*,*) '       9  change the time integration window width, weighting etc.'		!include here the straight FT option
							read(*,*) j_plot
					endif

					if(j_plot==0) exit map_loop
					if(abs(j_plot)>=1.and.abs(j_plot)<=5) cycle freq_loop
					if(abs(j_plot)==6) exit bz_loop

					if(j_plot==7) then
						j_oneph = j_oneph+1
						j_oneph = mod(j_oneph,2)
						if(data_type=='bulk') j_oneph = 0
						exit bz_loop              !new FT(Q) is needed
					endif

					if(j_plot==8) then
						j_ps = j_ps+1
						j_ps = mod(j_ps,2)
						cycle              
					endif

					if(j_plot==9) then
						write(*,*) 'Present values: j_weight,nfile,n_int,f_max,plot_size,j_fft:',j_weight,nfile,n_int,f_max,p_size,j_fft
						write(*,*) 'New values(after j_weight change go on via 6):'
						read(*,*) j_weight,nfile,n_int,f_max,p_size,j_fft
CC						if(j_weight==3.and.j_xray==0) then
CC							write(*,*) 'Xray form-factor data are missing:'
CC							write(*,*) '  locate first the XRAY_FF.TXT file, set J_XRAY=1 in the .PAR file and restart'
CC						endif
						if(2*(n_int/2)==n_int) n_int = n_int+1				!make it odd				
						n_freq_max = (n_int-1)/2+1					!f_max_plot/freq_step			!also n_freq_max= .5*t_int/t_step
						t_width = .5*n_int*t_step
						f_width = 1./t_width
						freq_step = 1./(n_int*t_step)
						write(*,*) 'Time-integration (BN window FWHM) [ps]',t_width
						write(*,*) 'Energy resolution [THz]',f_width

						deallocate(t_wind,om_phase)
						allocate(t_wind(n_int),SOURCE=0.0)
						allocate(om_phase(n_int),SOURCE=(.0,.0)) 
						do i=1,n_int
							t_wind(i) = .5*(1.-cos(twopi*(i)/real(n_int+1)))												!max is at n_int/2+1
						enddo
						cycle
					endif
				enddo

				flush(9)

			enddo freq_loop

			enddo bz_loop
			deallocate(cs_atom,cs_mean,x_ffq)				!,cross_sec
			enddo map_loop

			close(3)
			close(9)
			CALL PGEND
                  
      contains
      
C
      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set the RAINBOW palette of colors to be used by PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL JKL(9), JKR(9), JKG(9), JKB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA JKL / 0,   0.1,  0.2,  0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA JKR / 1.0, 0.5,  0.1,  0.1,  0.3,  1.0,  1.0, 1.0, 1.0/
      DATA JKG / 1.0, 0.7,  0.4,   .8,  1.,   1.0,  0.6, 0.0, 1.0/
      DATA JKB / 1.0, 1.0,  0.9,   .6,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.6) THEN
C        -- JK rainbow
         CALL PGCTAB(JKL, JKR, JKG, JKB, 9, CONTRA, BRIGHT)
      END IF
      END

      end program mp_sqom53


cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a FreeBSD license
cc (see license.txt in this directory). 
cc
      subroutine dirft2d1(nj,xj,yj,cj, iflag, ms,mt,fk)
      implicit none
      integer iflag
      integer(8) nj, ms, mt
      real(8) xj(nj), yj(nj)
      complex*16 cj(nj), fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2)
c     ------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c                   nj
c     fk(k1,k2) =  SUM cj(j) exp(+/-i k1 xj(j)) exp(+/-i k2 yj(j))
c                  j=1
c
c     for -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c
c***********************************************************************
      integer j, k1, k2
      complex*16 zf, cm1, z1n(-ms/2:(ms-1)/2)
c
      do k2 = -mt/2, (mt-1)/2
         do k1 = -ms/2, (ms-1)/2
            fk(k1,k2) = dcmplx(0d0,0d0)
         enddo
      enddo
c
      do j = 1, nj
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k1 xj)
c     ----------------------------------------------------------
c
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(xj(j)),+dsin(xj(j)))
         else
            zf = dcmplx(dcos(xj(j)),-dsin(xj(j)))
         endif
         z1n(0) = (1d0,0d0)
         do k1 = 1, (ms-1)/2
            z1n(k1) = zf*z1n(k1-1)
            z1n(-k1)= dconjg(z1n(k1))
         enddo
         if (ms/2*2.eq.ms) z1n(-ms/2) = dconjg(zf*z1n(ms/2-1))
c
c     ----------------------------------------------------------
c     Loop over k2 for yj
c     ----------------------------------------------------------
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(yj(j)),+dsin(yj(j)))
         else
            zf = dcmplx(dcos(yj(j)),-dsin(yj(j)))
         endif
c
         cm1 = cj(j)
         do k2 = 0, (mt-1)/2
            do k1 = -ms/2, (ms-1)/2
              fk(k1,k2) = fk(k1,k2) + cm1*z1n(k1)
            enddo
            cm1 = cm1*zf
         enddo
c
         zf = dconjg(zf)
         cm1 = cj(j)
         do k2 = -1, -mt/2, -1
            cm1 = cm1*zf
            do k1 = -ms/2, (ms-1)/2
              fk(k1,k2) = fk(k1,k2) + cm1*z1n(k1)
            enddo
         enddo
      enddo
      end
c
c
c