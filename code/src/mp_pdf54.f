  
      program mp_pdf54

C *************************************************************************************
C *****
C ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.54   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   calculates the pair distribution functions (PDF) for simulated supercell data
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
C ***** %%%%%%%%%%%%%%%%   		  program MP_PDF  1.54   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C ***** Ver. 1.1 The correlation is asymmetric: N_HSUM random cell against the complete supercell
C ***** Ver. 1.2 the whole supercell is read into memory once for ever
C ***** Ver. 1.21 the correlation is symmetric N_HSUM random cells against N_HSUM random cells
C *****
C ***** Ver. 1.31 between two delayed snapshot frame uses the 1.21 symmetric algorithm
C ***** Ver. 1.31 now the delay may be 0, so 1.31 replaces also ver. 1.21 and will override it
C *****
C ***** Ver. 1.32 the IDOM dialog moved after reading first input file
C ***** Ver. 1.32 the TIME_LOOP replaces the MASTER_LOOP
C ***** Ver. 1.32 normalization of correlation functions by true number of active cells n_norm = 45**3
C ***** 
C ***** Ver. 1.33 the reference snapshot and the first one in the time-series defined independently
C ***** Ver. 1.33 hence instantaneous (lag 0) correlations can be calculated explicitly
C ***** 
C ***** Ver. 1.4 	- flexible cell structure
C ***** 				 	- use of auxiliary <file_title.par> file (correlation pairs etc.)
C ***** Ver. 1.41	- corrected inconsistency in output <rdf>.txt file name
C *****						- correct histogram identification by j_hist_type
C ***** 
C ***** Ver. 1.42	- lighter & parallelied version (no correlation histograms, no IDOM)
C *****						- integrates the PGPLOT graphics of mp_rplot.f
C ***** 					- only PDF
C ***** 					- for correlation histograms (refer to mp_corr)
C ***** 					- modified .PAR file format, consequent .LOG
C ***** 
C ***** Ver. 1.50	- start of new series (identical to mp_pdf 1.42) as MP_PDF
C ***** 
C ***** Ver. 1.51	- all data arrays allocatable, no predefined array size limits
C ***** 					- supercell format in .PAR
C *****						- number of correlation pairs read from .PAR (don't use ver. 1.50!)
C *****
C ***** Ver 1.53  - record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****           - CELL and QUICK data type (BULK not eligible for PDF)
C *****						- secondary sampling volume clipped for highly anisotropic or 2D supercells
C *****
C ***** 
C ***** Ver 1.54  - NAMELIST I/O for .PAR files and .DAT headers implemented 
C *****           - PNG driver for PGPLOT included
C *****           
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
C ***** atom positions are converted from reduced lattice coordinates (x) to real space distances
C ***** 
  
  		use omp_lib
  		
CC  		use mp_nopgplot					! uncomment when not able to use PGPLOT, compile and link together with the module mp_nopgplot.f

      integer,parameter :: l_rec  =  1024		    !record length in real(4)
      real,parameter    :: pi = 3.14159

      character(4),allocatable  ::	at_name_par(:),at_label(:)
      character(10),allocatable ::	corr_name(:)

      integer,allocatable ::  at_ind_file(:,:,:),at_mask(:),j_corr(:,:),mult_corr(:) 

			real,allocatable ::  at_pos_file(:,:,:),rdf_part(:,:),rdf_tot(:),r(:)
      real,allocatable ::  r_plot(:),rdf_tot_plot(:),rdf_part_plot(:,:)
      real,allocatable ::  b_coh(:),at_weight(:),at_base(:,:),at_occup_par(:),at_pos(:,:),at_pos2(:,:)
	  
			integer ::       i_rdf,n_pdf,j_acc,n_int,n_int2,n_corr
      real    ::       rdf_dist,rdf_range,rdf_sum,rdf_tot_sum,rdf_norm,rand_rdf(3),at_weight_sum
      real 		:: 			 c_min,c_max,c_max2,c1,c2,r_start,r_end,pdf_step,c_smooth,f_smooth(51)
      
      character(4)   :: c_int(2),c_fil(2),version,head,atom
      character(10)  :: at_weight_scheme(2),pg_out,string,section,c_date,c_time,c_zone,c_nfile_min,c_nfile,c_jfile
      character(16)  :: sim_type_par,data_type,string16
      character(40)  :: file_master,file_inp,file_out,time_stamp,int_mode,x_file_name
      character(60)  :: file_dat,file_dat_t0,file_res,file_ps,file_log,plot_title,line
      character(128) :: cwd_path,plot_header
      character(l_rec):: header_record
      
			logical ::  nml_in,found_txt,found_ps,t_single

      integer(8) ::  n_in,n_out,n_inbox,n_outbox,n_pdf_range(3),n_mc,n_mc_max,n_norm
      integer ::  j_proc,proc_num,proc_num_in,thread_num,thread_num_max,m,j_name,jm,j1m,j2m,j_mode
      integer ::  i_start,i_end,j_part(4),n_step,nfile_step,n_hsum,n_h,j_head_in,hist_ind(3),j_verb,n_tot,n_head
      integer ::  i,j,k,ii,jj,ind,ind2,jat,i_rec,ind_rec,nrec,ncell,nrow,nlayer,nsuper,nfile,nfile_t0,nfile_min,nfile_max,jfile
      integer ::  j1,j2,j_plane,j_grid,j_logsc,j_ps,j_tot,j_txt,i_seed,i_r1,i_r11,i_r12,i_time(8),j_sr,j_pb,j_sr2,j_pb2
      integer ::  at_no,at_ind(3),at_ind2(3),d_ind(3),j_dom,j_base(3),n_plot,n_smooth,n_smooth_fwhm
      integer ::  j_atom,j_weight,j_xray,j_edit,j_mask,j_smooth,ifile,ifile0,jfil_step,jint,jint2,n_atom,sc_c2,sc_c1,ier,ios
      
      real :: t_dump,tt0,tt,b_sum,rand1,rand2,sc_r,at_displ

      real :: diff_pos(3),x_plot,y_plot

C **** the following variables MUST have the following type(4) or multiples because of alignement in the binary output file
      character(4),allocatable :: at_name_out(:)
      integer(4),allocatable   :: at_ind_in(:),nsuper_r(:)	

      real(4),allocatable ::	at_pos_in(:),at_occup_r(:)

      character(16)  :: sim_type,dat_type,input_method,file_par,subst_name,dat_source
      integer(4)     :: n_row(3),n_at,n_eq,j_force,j_shell_out,n_traj,n_cond,n_rec,n_tot_in,idum
      real(4)        :: rec_zero(l_rec),t_ms,t0,t1,a_par(3),angle(3),a_par_pdf(3),temp

      namelist /data_header_1/sim_type,dat_type,input_method,file_par,subst_name,t_ms,t_dump,temp,a_par,angle,
     1    n_row,n_atom,n_eq,n_traj,j_shell_out,n_cond,n_rec,n_tot                         !scalars & known dimensions
      namelist /data_header_2/at_name_out,at_occup_r,nsuper_r           !allocatables
     
      namelist /mp_gen/ j_verb,j_proc       
      namelist /mp_out/ j_weight,j_xray,j_logsc,j_txt,j_grid,pg_out       
      										!general rule: namelists of tools should only contain their local parameters
                          !what is of global interest they should pass into data_header
			namelist /mp_pdf/ n_pdf,pdf_step,a_par_pdf,n_h,j_acc,j_smooth,n_corr
			
C **** PGPLOT stuff
      INTEGER :: PGOPEN,j_xserv


C ********************* Initialization *******************************      

			write(*,*) '*** Program MP_PDF 1.54 ** Copyright (C) Jiri Kulda (2019,2021,2022) ***'
			write(*,*)'			*** for the moment only MonteCarlo integration is available *** '
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
			write(9,*) trim(time_stamp),'  MP_PDF 1.54  ',trim(cwd_path)
			write(9,*) 
      
CCCCCC !!!!!!!!! for the moment OpenMP is BLOCKED! !!!!!!!!!!!!!!!!!!!
CCCCC *********************  OpenMP initialization start  *******************************      
CCCCC
CCCC			proc_num_in = 1
CCCCCC			proc_num_in = j_proc
CCCC			if(j_verb.ge.1)write (*,*) 'PAR proc_num_in          = ',proc_num_in
CCCC			thread_num_max = omp_get_max_threads( )								!this gives maximum number of threads available (limited by concurrent tasks??)
CCCC			proc_num = omp_get_num_procs( )							!this should give number of processors, but in reality gives threads (cf. other unix/linux process enquiries)
CCCC			if(proc_num_in==0) proc_num_in = proc_num/2 !ask just for one thread per core	
CCCC			call omp_set_num_threads(proc_num_in)
CCCCCC			thread_num = omp_get_num_threads( )				!this gives threads available at this moment: here always = 1 as we are outside of a PARALLEL range
CCCC			
CCCC			if(j_verb.ge.1) then
CCCCCC				write (*,*) 'OMP processors available = ', proc_num
CCCC				write (*,*) 'OMP threads max          = ', thread_num_max
CCCC				write (*,*) 'OMP processes requested  = ', proc_num_in
CCCC										
CCCC			endif
CCCC
CCCC
CCCC 			if(proc_num_in.gt.1) then
CCCC				write(9,*) 'OMP threads max    = ', thread_num_max
CCCC				write(9,*) 'OMP processes requested 	= ', proc_num_in
CCCCCC				write(9,*) 'OMP processes allocated 		= ', proc_num
CCCCCC				write(*,*) 'OMP threads allocated    = ', thread_num
CCCC			else
CCCC				write(9,*) 'OpenMP not in use'
CCCC				write(*,*) 'OpenMP not in use'
CCCC			endif
CCCC			write(9,*) 
CCCCC
CCCCC ********************* OpenMP initialization end *******************************      

C *** other initialisations
			n_mc = 1000000					!unit for counting MC cycles (n_h)
			c_int = (/'GS','MC'/)
			c_fil = (/'g','m'/)

C *** Generate data file access
			write(*,*) 'Input data file_master: '
			read(*,*) file_master 
						
      write(*,*) 'Read data files number min, max (0 0 no numbers, single file): '
      read(*,*)   nfile_min,nfile_max
      nfile_step = 1
			t_single = nfile_min==0.and.nfile_max==0

			if(t_single)then
				nfile_min = 1		!for conformity with loop control conventions
				nfile_max = 1
			endif
			
      nfile = ((nfile_max-nfile_min)/nfile_step)+1
      
			if(t_single)then
				write(file_dat_t0,'("./data/",a,".dat")') trim(file_master)
			else
        if(nfile_min<=9999) then
          write(file_dat_t0,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min
        elseif(nfile_min>=10000) then
          write(string,'(i8)') nfile_min
          file_dat_t0 = './data/'//trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
        endif
			endif

	    open (1,file=file_dat_t0,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_dat_t0),' not found! Stop execution.'
				stop
			endif
			
C *** Read the header record 

      i_rec = 1   			
      read(1,rec=i_rec) header_record
      head = header_record(1:4)
      call up_case(head)
      if(head == 'MP_T') then      !new structure with namelist
        nml_in = .true.
				read(1,rec=i_rec) dat_source,version,string16
				read(string16,*) n_head
        write(*,*)		'Input data:  ',dat_source,version,n_head
        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_1)	
        t0 = t_dump
CC        write(*,nml=data_header_1)
      elseif(head.eq.'TIME'.or.head.eq.'STAT') then                                  !old w/o structure
        write(*,*)		'Input data:  ','old header format'
        nml_in = .false.
       read(1,rec=1) 
     1		  sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell_out,n_cond,n_rec					
         n_head = 1
         call up_case(sim_type)
         input_method = 'CELL'
				 if(index(sim_type,'BULK')/=0) input_method = 'BULK'
      else
      	write(*,*) 'header record wrong or missing'
      	write(*,*) trim(header_record)
      	stop
      endif 
      
      if(input_method == 'BULK') then
        write(*,*) 'BULK data are not suitable for this PDF algorithm, works with CELL data only.'
        stop
      endif

      allocate(at_name_out(n_atom),at_occup_r(n_atom),nsuper_r(n_atom))
			allocate(at_label(n_atom),at_name_par(n_atom),at_base(n_atom,3),at_weight(n_atom),at_mask(n_atom))
			allocate(at_occup_par(n_atom),b_coh(n_atom),SOURCE=.0)												!we need this to read .par
      at_mask = 1

      if(nml_in) then      !new structure with namelist
        i_rec = i_rec+1   							
        read(1,rec=i_rec) header_record
        read(header_record,nml=data_header_2)			
      else                                  !old w/o structure
			  read(1,rec=1) 
     1		sim_type,file_par,t_ms,t0,temp,a_par,angle,n_row,n_at,n_eq,j_force,j_shell_out,n_cond,n_rec,n_tot,			
     2    at_name_out,at_occup_r,nsuper_r																		
     		if(head.eq.'TIME') sim_type = 'TIMESTEP'
     		if (head.eq.'STAT') sim_type = 'STATIC'
      endif 
      close(1)
      
      
C **** Read the auxiliary file <file_par.par> with structure parameters, atom names and further info
      write(*,*) 'Parameter file name (.par to be added) (confirm or type other name): ', file_par
      read(*,*) file_par
      file_inp = trim(file_par)//'.par'

      open(4,file=file_inp,action='read',status ='old',iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_inp),' not found! Stop execution.'
				stop
			endif

      write(9,*) 'Read the parameter file:  ',trim(file_inp)

      read(4,nml=mp_gen)
      rewind(4)
      read(4,nml=mp_out)
      rewind(4)
      read(4,nml=mp_pdf)

CC			write(*,*) 'Sim_type, dat_type, input method: ',sim_type,dat_type,input_method		

C *** Read the atom positions       
      rewind(4)
      section = 'atoms'
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
          stop
        endif
        if(string(1:5).eq.section) exit	!find the mp_simple part of the .par file
      enddo
			do j=1,n_atom
        read(4,*) at_label(j),at_name_par(j),at_base(j,:),at_occup_par(j)	!for BULK the at_base and at_occup_par are not significant
			enddo
			do j=1,n_atom
				if(at_name_out(j)/=at_name_par(j)) then
					write(*,*) 'Atom names in .DAT and .PAR do not match: ',j,at_name_out(j),at_name_par(j)
					stop
				endif
CC				if(at_occup_r(j)/=at_occup_par(j)) then
CC					write(*,*) 'Atom occupancies in .DAT and .PAR do not match: ',j,at_occup_r(j),at_occup_par(j)
CC					stop
CC				endif
			enddo			

C *** Read the atom positions       
      rewind(4)
      section = 'pdf_pairs'
      do
        read(4,'(a)',iostat=ios) string
        if(ios/=0) then
          write(*,*) 'Section title:  ',trim(section),'  not found, check ', trim(file_inp)
          stop
        endif
        if(string(1:9).eq.section) exit	!find the mp_simple part of the .par file
      enddo

			allocate(corr_name(n_corr),j_corr(n_corr,2),mult_corr(n_corr))
			do i =1,n_corr		!get info on correlation pairs
				read(4,*,iostat=ios) corr_name(i)
				if(corr_name(i).eq.'#'.or.corr_name(i).eq.' '.or.ios/=0) exit
				backspace(4)
				read(4,*) corr_name(i),j_corr(i,:),mult_corr(i)
			enddo
			close(4)

C *** read neutron scattering lengths (always)
      open(4,file='neutron_xs.txt',action='read',status ='old',iostat=ios)
      if(ios.ne.0) then
        open(4,file='/usr/local/mp_tools/ref/neutron_xs.txt',action='read',status ='old',iostat=ios)
        if(ios.ne.0) then
          do
            write(*,*) 'File neutron_xs.txt not found, type in valid access path/filename'
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
          if(atom==trim(at_label(j))) then
            backspace(4)
            read(4,*) atom,b_coh(j)
            cycle atom_loop
          endif
        enddo
        write(*,*) 'b_coh for ',trim(at_label(j)),' not found,'
        write(*,*) 'check your spelling and the neutron_xs.txt table; use unit weights'
      enddo atom_loop

C *** write overview of atom data
			write(*,*)
      write(*,*) 'Substance name: ',subst_name	  
			write(*,*) 'Atoms from ',trim(file_inp)
      do j=1,n_atom
        write(*,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j)
        write(9,'(5x,a4,3f8.4,2x,2f8.4)')	at_name_par(j),at_base(j,:),at_occup_r(j),b_coh(j)
      enddo
      b_sum = sum(at_occup_r(1:n_atom)*b_coh(1:n_atom))

			write(*,*) 'partial correlations'
			do i =1,n_corr
				write(*,'(5x,i4,2x,a8,2x,2i4,i4)') i,trim(corr_name(i)),j_corr(i,:),mult_corr(i)
			enddo

C **** use the .PAR value of a_par in case of input in LATTICE units      
      if(a_par(1)*a_par(2)*a_par(3)==1.and.a_par_pdf(1)*a_par_pdf(2)*a_par_pdf(3)/=0.) a_par = a_par_pdf
      write(*,*) 'Setting a_par to',a_par
      


C *** define references for lattice search	
			nrow = n_row(1)
      nlayer = n_row(1)*n_row(2) 						! 48**2 gives   2304
      nsuper = n_row(1)*n_row(2)*n_row(3) 	! 48**3 gives 110592  
      n_tot = n_atom*nsuper   
CC    	write(*,*) 'n_row',n_row
      rdf_range = n_pdf*pdf_step
      do j=1,3
      	n_pdf_range(j) = min(nint(rdf_range/a_par(j)),n_row(j))		!clip the integration range for strongly anisotropic supercells
      enddo
      n_norm = 8*n_pdf_range(1)*n_pdf_range(2)*n_pdf_range(3)		!number of unit cells within the PDF reach from R1
      n_mc_max = nsuper*(n_norm-1) !number of independent site correlation pairs
CC      write(*,*) 'n_pdf_range,n_norm,n_mc_max',n_pdf_range,n_norm,n_mc_max

C *** Allocate and clear the histogram arrays for accumulation across several snapshots	       
			allocate (at_pos_file(4,nsuper,n_atom),at_ind_file(4,nsuper,n_atom),at_pos_in(4*n_tot),at_ind_in(4*n_tot))
			allocate(at_pos(n_atom,3),at_pos2(n_atom,3),r(n_pdf),rdf_tot(n_pdf),rdf_part(n_corr,n_pdf))
			allocate(r_plot(n_pdf),rdf_tot_plot(n_pdf),rdf_part_plot(n_corr,n_pdf))

C *********************      Initial info output   ************************

      
CC *** for the moment the MC integration strategy is imposed!

			j_acc = 2
			
			n_h = n_mc_max/1000000
			write(string,"(i8)") n_h
			write(*,*) 'MC sampling pairs per frame ([x 10^6], ',trim(adjustl(string)),' max)'
			read(*,*)   n_h
			n_hsum = n_h*n_mc !number of MC cycles per snapshot
			if(n_hsum.gt.n_mc_max) then
				write(*,*) 'WARNING: n_hsum exceeds max number of cell pairs ',n_mc_max
				write(*,*) 'input reduced n_h (<',n_mc_max/n_mc,'):'
				read(*,*)   n_h
				n_hsum = n_h*n_mc !number of MC cycles per snapshot
			endif
			n_int = n_hsum	! for MC cycle over n_hsum randomly chosen cell pairs
			n_int2 = 1      ! each MC step deals with a single R1-R2 site pair
			call random_seed	!initiate the random number sequence for MC sampling
			int_mode = 'Monte Carlo'
			write(*,*) 'Monte Carlo integration over',n_int*n_int2,'cell pairs'
			
C *** cycle over PDFs with different atom at_weights

			at_weight_scheme(1) = 'Uniform'
			at_weight_scheme(2) = 'Neutron'
			at_mask = 1
			
			at_weights_loop: do			
				
CC				write(*,*) 'at_weighting scheme (1 uniform, 2 neutrons)? Edit atom masks (1/0)? (0 0 END)'
				write(*,*) 'Atom weights: 1= uniform, 2= neutron coherent (b_c^2); 0= EXIT, negative = EDIT'	
				do
					read(*,*) jj
					if(jj==0) exit at_weights_loop
					if(abs(jj)==1.or.abs(jj)==2) exit
					write(*,*) 'Input out of range, repeat ...'
				enddo

				j_weight = abs(jj)
				at_weight = .0
				if(j_weight==1) then
					at_weight(1:n_atom) = 1.
				elseif(j_weight==2) then
					at_weight(1:n_atom) = b_coh(1:n_atom)
				endif
				write(9,*)
				write(*,'(1x,"Atoms:         ",50a8)')  (at_name_par(i),i=1,n_atom)
				write(9,'(1x,"Atoms:         ",50a8)')  (at_name_par(i),i=1,n_atom)
				write(*,'(1x,"Atoms no.:  ",50i8)') (i,i=1,n_atom)
				write(9,'(1x,"Atoms no.:  ",50i8)') (i,i=1,n_atom)
				write(*,113) trim(at_weight_scheme(j_weight)),(at_weight(i),i=1,n_atom)
				write(9,113) trim(at_weight_scheme(j_weight)),(at_weight(i),i=1,n_atom)
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
				at_weight_sum = sum(at_weight)
				
C *************  the file loop: cycle over the tt files to augment statistics

      rdf_part = .0
      rdf_tot = .0

      call cpu_time(t0)
      
      if(t_single)then
      	nfile_min = 1
      	nfile_max = 1
      	nfile_step = 1
      endif
      
      file_loop: do ifile=nfile_min,nfile_max,nfile_step									

C ***  open the t0 file (binary MD snapshot file)
			if(t_single)then
				write(file_dat,'("./data/",a,".dat")') trim(file_master)
			else
        if(nfile_min<=9999) then
          write(file_dat,'("./data/",a,"_n",i4.4,".dat")') trim(file_master),nfile_min+(ifile-nfile_min)*nfile_step
        elseif(nfile_min>=10000) then
          write(string,'(i8)') nfile_min+(ifile-nfile_min)*nfile_step
          file_dat = './data/'//trim(file_master)//'_n'//trim(adjustl(string))//'.dat'
        endif
			endif

			write(*,*)
			write(*,*)'input: ',file_dat
			write(9,*)'input: ',file_dat

		  open(1,file=file_dat,status ='old',access='direct',action='read',form='unformatted',recl=4*l_rec,iostat=ios)
			if(ios.ne.0) then
				write(*,*) 'File ',trim(file_dat),' not opened! IOS =',ios
				write(*,*) 'Skip(1), stop execution(0)?'
				read(*,*) jj
				if(jj==1) exit file_loop
				if(jj==0) stop
			endif

			i_rec = n_head	
			do j=1,n_rec-1
				i_rec = i_rec+1
				read(1,rec=i_rec) at_ind_in((j-1)*l_rec+1:j*l_rec)			
			enddo	
			i_rec = i_rec+1
			read(1,rec=i_rec) at_ind_in((n_rec-1)*l_rec+1:4*n_tot)			
			
			do j=1,n_rec-1
				i_rec = i_rec+1
				read(1,rec=i_rec) at_pos_in((j-1)*l_rec+1:j*l_rec)			
			enddo	
			i_rec = i_rec+1
			read(1,rec=i_rec) at_pos_in((n_rec-1)*l_rec+1:4*n_tot)	
			
			at_ind_file(:,:,:) = reshape(at_ind_in,(/4,nsuper,n_atom/))
			
			at_pos_file(:,:,:) = reshape(at_pos_in,(/4,nsuper,n_atom/))
			
			close(1)
			
			write(*,*) 'Accumulating the PDFs ...'     

C *************  the integration loop: cycle over the site pairs to accumulate the PDFs ******************

			n_inbox = 0
			n_outbox = 0
			n_in = 0
			n_out = 0

C *** MonteCarlo integration
CC			elseif(j_acc==2) then
				write(*,*)
				CALL SYSTEM_CLOCK (COUNT = sc_c1, COUNT_RATE = sc_r)				!, COUNT MAX = sc_m)
				j_sr = 0
				j_pb = 0
				j_sr2 = 0
				j_pb2 = 0
				do k=1,n_h
	
CC	      integration_loop_2: do jint = 1,n_int
					integration_loop_2: do jint = 1,n_mc/1000
						call random_number(rand1)			!the first cell
						i_r1 = nsuper*rand1+1	!i_r1 may become 1 to nsuper

						do j = 1,n_atom
							at_pos(j,:) = at_pos_file(1:3,i_r1,j)
						enddo

						at_ind(3) = (i_r1-1)/nlayer+1
						at_ind(2) = mod((i_r1-1),nlayer)/n_row(1)+1
						at_ind(1) = mod(i_r1-1,n_row(1))+1

						ii=1
							if(at_pos(ii,1).eq.0..and.at_pos(ii,2).eq.0..and.at_pos(ii,3).eq.0.) j_sr=j_sr+1
						ii=2
							if(at_pos(ii,1).eq.0..and.at_pos(ii,2).eq.0..and.at_pos(ii,3).eq.0.) j_pb=j_pb+1
						
					
					do m=1,1000
					do i=1,3
						call random_number(rand2)		!find the second cell within a cube of 2*n_pdf_range around the first one
						d_ind(i) = anint(2.*n_pdf_range(i)*(rand2-.5))
						if(n_pdf_range(i)==1) d_ind(i) = 0
					enddo

C *** skip selfcorrelations
					if(d_ind(1)==0.and.d_ind(2)==0.and.d_ind(3)==0) then
						n_outbox = n_outbox+25
						n_out = n_out+13
						cycle
					endif

          at_ind2 = at_ind+d_ind
					
C *** treat atoms out of the supercell by applying cyclic boundary conditions
        	j_base = 0
          do ii =1,3
            if(at_ind2(ii).lt.1) then
              at_ind2(ii) = at_ind2(ii)+n_row(ii)
              j_base(ii) = -n_row(ii)
            endif

            if(at_ind2(ii).gt.n_row(ii)) then
              at_ind2(ii) = at_ind2(ii)-n_row(ii)
              j_base(ii) = n_row(ii)
            endif
          enddo

C *** retrieve the atom positions
					ind2 = nlayer*(at_ind2(3)-1)+n_row(1)*(at_ind2(2)-1)+at_ind2(1)
					if(ind2.lt.1.or.ind2.gt.nsuper) then
						write(*,*) 'out of range: at_ind2,ind2',at_ind2,ind2,nsuper
          endif
          
          do j = 1,n_atom
          	at_pos2(j,:) = at_pos_file(1:3,ind2,j)
						if(at_pos2(j,1)/=.0.and.at_pos2(j,2)/=.0.and.at_pos2(j,3)/=0.) at_pos2(j,:) = at_pos2(j,:)+j_base	
          enddo

CC								write(*,*) 'at_pos2',(at_pos2(j,:),j=1,n_atom)

						ii=2
							if(at_pos2(ii,1).eq.0..and.at_pos2(ii,2).eq.0..and.at_pos2(ii,3).eq.0.) j_pb2=j_pb2+1
						ii=1
							if(at_pos2(ii,1).eq.0..and.at_pos2(ii,2).eq.0..and.at_pos2(ii,3).eq.0.) then
								j_sr2=j_sr2+1
							endif

	
C *** accumulate the total PDF first
					do ii=1,n_atom
						do jj=1,n_atom
C *** handle unoccupied positions 
							if(at_pos(ii,1).eq.0..and.at_pos(ii,2).eq.0..and.at_pos(ii,3).eq.0.) cycle
							if(at_pos2(jj,1).eq.0..and.at_pos2(jj,2).eq.0..and.at_pos2(jj,3).eq.0.) cycle
							diff_pos = at_pos(ii,:)-at_pos2(jj,:) 
							i_rdf = anint(norm2(diff_pos*a_par)/pdf_step)
						  
							if(i_rdf.ge.1.and.i_rdf.le.n_pdf) then
								rdf_tot(i_rdf) = rdf_tot(i_rdf)+at_weight(ii)*at_weight(jj)
								n_inbox = n_inbox+1
							else
								n_outbox = n_outbox+1
							endif
						enddo	!jj
					enddo		!ii

C *** accumulate the partials
					if(j_corr(1,1).ne.0.and.j_corr(1,2).ne.0) then
						do j=1,n_corr
							j1 = j_corr(j,1)
							j2 = j_corr(j,2)

C *** handle unoccupied positions 
							if(at_pos(j1,1).eq.0..and.at_pos(j1,2).eq.0..and.at_pos(j1,3).eq.0.) cycle
							if(at_pos2(j1,1).eq.0..and.at_pos2(j1,2).eq.0..and.at_pos2(j1,3).eq.0.) cycle
							if(at_pos(j2,1).eq.0..and.at_pos(j2,2).eq.0..and.at_pos(j2,3).eq.0.) cycle
							if(at_pos2(j2,1).eq.0..and.at_pos2(j2,2).eq.0..and.at_pos2(j2,3).eq.0.) cycle

							diff_pos = at_pos2(j2,:)-at_pos(j1,:)
							i_rdf = anint(norm2(diff_pos*a_par)/pdf_step)
							if(i_rdf.ge.1.and.i_rdf.le.n_pdf) then
									rdf_part(j,i_rdf) = rdf_part(j,i_rdf)+1.
									n_in = n_in+1
							else
									n_out = n_out+1
							endif

							if(j1/=j2)then													!for the off-diagonals do also the other way round
								diff_pos = at_pos2(j1,:)-at_pos(j2,:)									
								i_rdf = anint(norm2(diff_pos*a_par)/pdf_step)
								if(i_rdf.ge.1.and.i_rdf.le.n_pdf) then
									rdf_part(j,i_rdf) = rdf_part(j,i_rdf)+1.
									n_in = n_in+1
								else
									n_out = n_out+1
								endif
							endif
						enddo 
					endif
					enddo !m

					enddo integration_loop_2
					CALL SYSTEM_CLOCK (COUNT = sc_c2)
					if(k==10*(k/10).or.k==n_h) write(*,*) (sc_c2-sc_c1)/sc_r,k*n_mc
				enddo		!k=1,n_h
      enddo file_loop

C *** pre-scale PDF and accumulate rdf_tot, rdf_sum     
      rdf_tot_sum = 0.		!this and all following Cn's would serve "empirical" PDF normalisation
      rdf_sum = 0.

				do i=1,n_pdf
					rdf_tot(i) = rdf_tot(i)/i**2
					rdf_tot_sum=rdf_tot_sum+rdf_tot(i)
				enddo

      do j=1,n_corr
         do i=1,n_pdf
          rdf_part(j,i) = at_weight(j_corr(j,1))*at_weight(j_corr(j,2))*mult_corr(j)*rdf_part(j,i)/i**2
CC          if(i.gt.n_pdf/2) rdf_sum=rdf_sum+rdf_part(j,i)
          rdf_sum=rdf_sum+rdf_part(j,i)
        enddo
      enddo

      rdf_norm = 4.*pi*pdf_step**3*at_weight_sum**2*n_int*nfile/(n_norm*a_par(1)*a_par(2)*a_par(3)) ! PDF norm from "first principles"      
      rdf_tot = rdf_tot/rdf_norm
      rdf_part = rdf_part/rdf_norm
 
C *** generate the abscissa values for later output
      do i=1,n_pdf
				r(i) = i*pdf_step
			enddo

C *** open the PGPLOT graphics window (X11)
      j_xserv = PGOPEN('/xserv')   
      if (j_xserv.LE.0) then    
      	write(*,*) 'Could not open PGPLOT /xserv'
      	STOP
      endif
CC      CALL PGPAP(7.0,1.5)     ! define the plot area as portrait		!version with 4 subpanes
CC	  	call PGSUBP(1,4)				! PGPLOT window of 1x4 no of panes
      CALL PGPAP(10.0,.6)     ! define the plot area as landscape
	  	call PGSUBP(1,1)				! PGPLOT window of 1x1 no of panes
      CALL PGASK(.FALSE.)     ! would not ask for <RET>
	  	CALL PGSCRN(0, 'white', IER)	!plot on white background
	  	CALL PGSCRN(1, 'black', IER)
	  	
	  	
C *** Prepare the data to be plotted and smooth the PDF
     	plot_loop: do
				write(*,*) 'Plot range r_1, r_2, max =',n_pdf*pdf_step,' r_2=0 END'
				read(*,*) r_start,r_end 
				if(r_start.lt.pdf_step) r_start = pdf_step
				if(r_end.eq.0.) exit plot_loop  
				i_start = (r_start-r(1))/pdf_step + 1
				i_end = (r_end-r(1))/pdf_step+1
				n_plot = i_end-i_start+1

				if(j_smooth==0) then
					write(*,*) 'Gaussian smooth FWHM npts (1 no smoothing, 4 - 10 useful)'
					read(*,*) n_smooth_fwhm
				else
					n_smooth_fwhm = j_smooth			!value from the .PAR file
				endif
			
				n_smooth = 5*n_smooth_fwhm/2+1
				if(n_smooth_fwhm==1) n_smooth=1			!no smoothing

				do j=1,n_smooth
					f_smooth(j) = 2.**(-((j-n_smooth/2-1)/(.5*n_smooth_fwhm))**2) 
				enddo
				f_smooth = f_smooth/sum(f_smooth(1:n_smooth))				!profile normalized to unit integral

				if(n_smooth>1) write(*,*) 'Applying Gaussian smoothing with FWHM=',n_smooth_fwhm,' steps of',pdf_step,'[A^]'
CC				write(*,*) 'f_smooth', (f_smooth(j),j=1,n_smooth)		

				rdf_tot_plot = .0
				rdf_part_plot = .0
				
				if(n_smooth==1) then
					do i = 1,n_plot
						r_plot(i) = r(i_start+i-1)
						rdf_tot_plot(i)=rdf_tot(i_start+i-1)
						rdf_part_plot(:,i)=rdf_part(:,i_start+i-1)
					enddo               
				else
					do i = 1,n_pdf
						do j = 1,n_smooth
							if (i-n_smooth/2+j-1.le.0.or.i-n_smooth/2+j-1.gt.n_pdf) cycle  !zeros out of range
							rdf_tot_plot(i)=rdf_tot_plot(i)+rdf_tot(i-n_smooth/2+j-1)*f_smooth(j)
							rdf_part_plot(:,i)=rdf_part_plot(:,i)+rdf_part(:,i-n_smooth/2+j-1)*f_smooth(j)
						enddo 
					enddo
	      endif

C *** Plot the results

      c_min = minval(rdf_tot_plot)
      c_min = min(c_min,-.1)
      c_max = maxval(rdf_tot_plot)
      c_max = .1*(int(10*c_max)+1)
      c_max2 = maxval(rdf_part_plot)
      c_max2 = .1*(int(10*c_max2)+1)

      partials_loop: do
        write(*,*)
        write(*,*) 'Plotting ',trim(subst_name),
     1        ':  1/0 for PDF_tot & 4 indices of partial PDFs to be added (0 = skip, -1 = end)'
        read(*,*) j_tot,j_part
        if(j_part(1).gt.n_corr.or.j_part(2).gt.n_corr.or.j_part(3).gt.n_corr.or.j_part(4).gt.n_corr) cycle
        if (j_part(1).eq.-1.or.j_tot==-1) exit partials_loop
        write(plot_header,'("PDF plot ",a,"    ",a,"weights")') trim(file_dat_t0),at_weight_scheme(j_weight)
				do
          write(plot_title,114) trim(subst_name)
114       format(a)   
					CALL PGSLCT(j_xserv)
          CALL PGSCI (1)  !white
          CALL PGSCH(1.)
          CALL PGSLW(2)
          CALL PGSLS (1)  !full
          CALL PGENV(r_start,r_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
          CALL PGLAB('r(A)', 'g(r)',trim(plot_header))  !put the axis labels
          CALL PGSCI (1)  !white
          CALL PGSLW(5)			!operates in steps of 5
						x_plot = r_start+.75*(r_end-r_start)
						y_plot = .9*c_max
          if(j_tot==1) then
						CALL PGLINE(n_plot,r_plot,rdf_tot_plot)  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
						CALL PGTEXT (x_plot,y_plot,plot_title)
						CALL PGSLW(2)
					endif

          do j=1,4
          	if(j_part(j)==0) cycle
          	 write(plot_title,115) at_name_par(j_corr(j_part(j),1)),at_name_par(j_corr(j_part(j),2))
          	 y_plot = (.9-.06*j)*c_max
115       	format(a,a)   
CC            CALL PGSCI (1)  !white		!version with 4 subpanes
CCc            CALL PGSCI (0)  !black
CC            CALL PGSLS (1)  !full
CC            CALL PGENV(r_start,r_end,c_min,c_max2,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes
CC            CALL PGLAB('r(A)', 'g(r)',plot_title)  !put the axis labels
            CALL PGSCI (j+1)  !red-green-blue
            CALL PGLINE(n_plot,r_plot,rdf_part_plot(j_part(j),:))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title)
          	CALL PGSLW(2)			
          enddo
          write(*,*) 'Modify vertical scale (min,max; 0,0 EXIT)'
          read(*,*) c1,c2
          if(c1==.0.and.c2==.0) exit
          c_min = c1
          c_max = c2
        enddo

          write(*,*) ' 1   edit partials '
          write(*,*) ' 2   save graph (.PS) & edit partials  '
          write(*,*) ' 3   save graph (.PS) & PDF (.TXT) '
          write(*,*) ' 4   change plot range '
          write(*,*) ' 5   edit at_weights (restart integration) '
          write(*,*) ' 0   EXIT'  

          read(*,*) jj
					j_ps = 0
          j_txt = 0
          if(jj==0) exit at_weights_loop
          if(jj==1) cycle partials_loop
          if(jj==2) j_ps = 1
          if(jj==3) then
          	j_ps = 1
          	j_txt = 1
          endif
          if(jj==4) cycle plot_loop
          if(jj<0.or.jj>4) cycle at_weights_loop
          
          
C **** Prepare and plot the same on .PS, look for existing output files in order not overwrite them				
				if(j_ps==1) then
					jfile = 1
					do						!look for existing .ps files to continue numbering

			if(j_name==0.and.t_single)then
						write(file_ps,1041) trim(file_master),trim(c_fil(j_acc)),jfile
1041   		format(a,'_rdf',a,'_',i2.2,'.ps')      
						write(file_res,1042) trim(file_master),trim(c_fil(j_acc)),jfile
1042   		format(a,'_rdf',a,'_',i2.2,'.txt')      
			else
            write(c_jfile,'("_",i2.2)') jfile
            if(nfile_min<=9999) then
              write(c_nfile_min,'(i4.4)') nfile_min
            elseif(nfile_min>=10000) then
              write(c_nfile_min,'(i8)') nfile_min
            endif
            c_nfile_min = '_'//adjustl(c_nfile_min)
    
            if(nfile<=9999) then
              write(c_nfile,'(i4.4)') nfile
            elseif(nfile>=10000) then
              write(c_nfile,'(i8)') nfile
            endif
            c_nfile = '_'//adjustl(c_nfile)
    
            file_res = trim(file_master)//'_pdf'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//'.txt'							
            file_ps  = trim(file_master)//'_pdf'//trim(c_nfile_min)//trim(c_nfile)//trim(c_jfile)//'.ps'
			endif
						inquire(file=file_ps,exist=found_ps)
						inquire(file=file_res,exist=found_txt)
						if(.not.found_txt.and.(.not.found_ps)) exit				
						jfile = jfile+1
						if(jfile==100) then
							write(*,*)'Tidy up .txt/.ps files to restart count from 01 and type [RET]'
							read(*,*)
							jfile = 1
						endif	
					enddo						

					ier = PGOPEN(file_ps//'/CP')
					IF (ier.LE.0) STOP
					CALL PGASK(.FALSE.)     ! would not ask for <RET>
					CALL PGPAP(11.0,.7)     ! define the plot area as landscape
					CALL PGSUBP(1,1)				! PGPLOT window of 1x1 no of panes
					CALL PGSCRN(0, 'white', IER)	!plot on white background
					CALL PGSCRN(1, 'black', IER)

					plot_title = file_ps	
 
          CALL PGSCI (1)  !white
          CALL PGSCH(1.)
          CALL PGSLW(2)
          CALL PGSLS (1)  !full
          CALL PGENV(r_start,r_end,c_min,c_max,0,j_grid+1) !PGENV(xmin,xmax,ymin,ymax,0,1) - draw the axes w/o grid
          CALL PGLAB('r(A)', 'g(r)',trim(plot_title))  !put the axis labels
          CALL PGSCI (1)  !white needs to be reset after PGLAB
          CALL PGSLW(5)			!operates in steps of 5
          CALL PGLINE(n_plot,r_plot,rdf_tot_plot)  !plots the curve
					y_plot = .9*c_max
					CALL PGSTBG(0)																				 !erase graphics under text
					CALL PGTEXT (x_plot,y_plot,'PDF_tot '//trim(subst_name)) !the x_plot,y_plot are in world (axis units) coordinates
          CALL PGSLW(2)

          do j=1,4
          	 write(plot_title,115) at_name_par(j_corr(j_part(j),1)),at_name_par(j_corr(j_part(j),2))
          	 y_plot = (.9-.06*j)*c_max
            CALL PGSCI (j+1)  !red-green-blue
            CALL PGLINE(n_plot,r_plot,rdf_part_plot(j_part(j),:))  !plots the curve
						CALL PGSTBG(0)																				 !erase graphics under text
          	CALL PGSLW(5)			!operates in steps of 5
						CALL PGTEXT (x_plot,y_plot,plot_title)				!the x_plot,y_plot are in world (axis units) coordinates
          	CALL PGSLW(2)
          enddo    
					CALL PGCLOS
					write(*,*) ' Postscript output written to: ',file_ps	
					write(9,*)
					write(9,*) '  ',trim(int_mode),' integration',n_int*n_int2,' pairs'
					write(9,*) '  ','Masks:',(at_mask(i),i=1,n_atom)				
					if(n_smooth>1) write(9,*) '  ','Smoothing FWHM:',n_smooth_fwhm	
					write(9,*) ' Postscript output written to: ',file_ps
				endif 																							!j_ps
				if(j_txt==1) exit plot_loop		!save PDF to .TXT
        enddo partials_loop

        enddo plot_loop
 
C *** save the PDF results into an ASCII file (each line corresponds to a distance point)
C     look for existing output files in order not overwrite them				
			if(j_txt==1) then
				open (4,file=file_res)

C *** write PDF file header
				write(4,*) 'Substance:',trim(subst_name),'       ',time_stamp
				write(4,*) 'Input files: ',trim(file_dat_t0),' to ',trim(file_dat),' step',nfile_step									   
				write(4,*) 'Supercell size:',nrow																				
				write(4,*) 'Integration ',trim(int_mode),'  ',n_int*n_int2,'cell pairs'
				write(4,*) 'Atoms/unit_cell:',n_atom
				write(4,*) 	 'Atoms  :',(('    '//at_name_par(i)),i=1,n_atom)
				write(4,'(1x,"Atom no. ",50i8)') (i,i=1,n_atom)
				write(4,125) (at_occup_r(i),i=1,n_atom)
				write(4,126) (at_weight(i),i=1,n_atom)
				write(4,127) (at_mask(i),i=1,n_atom)
125     format(1x,"Occup's: ",50f8.4)								
126     format(1x,'at_weights: ',50f8.4)								
127     format(1x,'Masks	: ',50i8)								
				write(4,*)	
				write(4,*)	'Atom pair partials:'
				write(4,106) (j_corr(i,1),i=1,n_corr)
				write(4,106) (j_corr(i,2),i=1,n_corr)
				write(4,*)  '  r(A^)   PDF_tot   PDF partials (cf. above)'

C *** write the PDFs
				do i=1,n_pdf
					write(4,105) r(i),rdf_tot(i),(rdf_part(j,i),j=1,n_corr)
				enddo
105   	format(1x,f7.3,1x,f8.3,(1x,10f8.3))	
106   	format(15x,(1x,10i8))	
				close(4)
				write(*,*) ' Text output written to: ',file_res	  
				write(9,*) ' Text output written to: ',file_res	  

				write(*,*)
				write(*,*) ' 0   EXIT'  
				write(*,*) ' 1   edit at_weights (restart integration) '
				read(*,*) jj
				if(jj==0) exit at_weights_loop
			endif	
      enddo at_weights_loop

			flush(9)
			close(9)
			            
      end program mp_pdf54
  
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


   
  
     