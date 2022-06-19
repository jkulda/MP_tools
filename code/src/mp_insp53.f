     
      program mp_insp53

C *************************************************************************************
C *****
C *****  %%%%%%%%%%%%%%%%   		  program MP_INSP 1.53   		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C *****
C *****   inspects single atom properties in a binary MP_TOOLS data file  
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
C ***** %%%%%%%%%%%%%%%%           program MP_INSP 1.53   					 %%%%%%%%%%%%%%%%%%%%%%%%
C *****
C ***** Ver. 1.1 - original version
C ***** Ver. 1.2  
C ***** Ver. 1.3 - reads the extended header with n_atom, j_force and temp (on end positions for compatibility)
C *****	     		 - if atom forces are present, they are printed as well
C ***** Ver. 1.4 - dialogue revision
C *****					 - ./data subdirectory used for HIST and primary .bin files
C ***** Ver. 1.41 - minor bug corrections
C *****
C ***** Ver. 1.51 - non-cubic supercell with n_row(3) in header
C ***** Ver. 1.52 - provision for filenames without number
C *****
C ***** Ver. 1.53 - uses new semi-sequential binary file format 
C *****						- record length 4096, all data 32bit length
C *****						- at_mass is saved as at_veloc(4,j,k)
C *****						- at_charge is saved as at_pos(4,j,k)
C *****
C ***** inspects single atom properties in a binary .dat file
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
C ***** atom positions are converted to and recorded in reduced lattice coordinates (x) 
C ***** 
      integer,parameter :: l_rec  =  1024		    !record length in real(4)
      
      character(4),allocatable      :: at_name(:)

      character(16)     :: sim_type,file_title,data_type
      character(40)     :: file_master
      character(60) 		:: file_dat,file_out
      character(128)    :: header_record

      integer(4),allocatable ::  at_ind(:),ind_at(:),nsuper_r(:)    
      real(4),allocatable ::  at_pos_c(:),at_veloc_c(:),at_force_c(:),at_occup(:)      
      real(4),allocatable ::  at_pos_s(:),at_veloc_s(:),at_force_s(:)      
      
      integer(4) ::  at_no,n_traj,n_cond,n_atom,n_eq,j_force,j_shell,jpos_in,jrec_in,jrec_save,n_rec,n2_rec,l_rec4,i_rec
      integer(4) ::  ios,idum,i,j,k,jj,jat,ind,nfile,ncell,n_tot,nsuper,nlayer,nrow,n_row(3),jt,j_layer,j_row,j_pos
      
      integer(4) :: j_name

      real :: a_par(3),angle(3),t_ms,t_step,temp
     
C ********************* Initialization *******************************      

			write(*,*) '*** Program MP_INSP 1.53 ** Copyright (C) Jiri Kulda (2022) ***'

			j_force = 0
			l_rec4 = l_rec/4

      write(*,*) 'Master filename: '
      read(*,*) file_master 

C **** open a t-snapshot file and read its header 
 
      file_loop: do 

				jrec_save = 0

				write(*,*) 'snapshot file number (0 for no number):'
				read(*,*) jt

				if(jt==0)then
					file_dat = './data/'//trim(file_master)//'.dat'
				else
					write(file_dat,101) trim(file_master),jt
				endif 			
101   	format('./data/',a,'_n',i4.4,'.dat')

				write(*,*)'Reading ',file_dat
				write(*,*) 
				open (1,file=file_dat,action='read',status ='old',access='direct',form='unformatted',recl=4*l_rec)

   			i_rec = 1
				read(1,rec=i_rec) sim_type,file_title,t_ms,t_step,temp,a_par,angle,n_row,n_atom,n_eq				

				allocate(at_name(n_atom*n_eq),at_occup(n_atom*n_eq),nsuper_r(n_atom*n_eq))

				read(1,rec=i_rec) 
     1		sim_type,file_title,t_ms,t_step,temp,a_par,angle,n_row,n_atom,n_eq,j_force,j_shell,
     2    n_cond,n_rec,n_tot,at_name,at_occup,nsuper_r																		!char(16): sim_type,file_title, char(4): at_name

				data_type = 'cell'
				if(index(sim_type,'bulk')/=0) data_type = 'bulk'
			  if(index(sim_type,'quick')/=0) data_type = 'quick'
								
				write(*,*) 'Simulation & data type:              ',sim_type,'    ',data_type
				write(*,*) 'Time structure t_ms,t_step:       ',t_ms,t_step
				write(*,*) 'Supercell & temperature:    ',n_row,temp
				write(*,*) 'Unit cell parameter(3), angle(3): ',a_par,angle
				write(*,*) 'Atom numbers:              ',n_atom,nsuper_r,n_tot
				write(*,*) 'Atoms & occupancies:                 ',at_name,at_occup
				write(*,*) 
				
				nrow = n_row(1)
				nlayer = n_row(1)*n_row(2)						!only to be used for record number calculation
				nsuper = n_row(1)*n_row(2)*n_row(3)
CC				n_tot = n_atom*nsuper
				
				allocate(at_ind(l_rec),at_pos_c(l_rec),at_veloc_c(l_rec),at_force_c(l_rec))
				if(j_shell==1) allocate(at_pos_s(l_rec),at_veloc_s(l_rec),at_force_s(l_rec))
				allocate(ind_at(n_atom))
				
				ind_at(1) = 0
				do jat = 2,n_atom
					ind_at(jat) = ind_at(jat-1)+nsuper_r(jat-1)
				enddo

				master_loop: do 	
					if(data_type=='cell'.or.data_type=='quick')	then
						write(*,*) 'jat, j_pos,j_row,j_layer (0 0 0 0 = new file number, 9 9 9 9 = END): '
						read(*,*) jat,i,j,k
						if (jat==0) exit master_loop
						if (jat==9.and.i==9.and.j==9.and.k==9) exit file_loop
						if(jat>n_atom) cycle

						ind = nsuper*(jat-1)+nlayer*(k-1)+n_row(1)*(j-1)+i
					else !data_type == 'bulk'
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
						i_rec = jrec_in+1
						read(1,rec=i_rec) at_ind
					
						i_rec = jrec_in+n_rec+1
						read(1,rec=i_rec,iostat=ios) at_pos_c
						if(ios/=0) write(*,*)'at_pos_c: ios,i_rec',ios,i_rec
					
						if(index(sim_type,'static')==0.and.index(sim_type,'Static')==0.and.index(sim_type,'STATIC')==0) then
							i_rec = jrec_in+2*n_rec+1
							read(1,rec=i_rec) at_veloc_c

							if(j_force==1) read(1,rec=jrec_in+3*n_rec+1) at_force_c  
				
							n2_rec = (3+j_force)*n_rec+1
							if(j_shell==1) then
								read(1,rec=jrec_in+n2_rec) at_pos_s
								read(1,rec=jrec_in+n2_rec+n_rec) at_veloc_s
								if(j_force==1) read(1,rec=jrec_in+n2_rec+2*n_rec) at_force_s  					
							endif
						endif
						jrec_save = jrec_in
					endif

					jj = 4*(jpos_in-1)
					write(*,*) 'Atom, type, mass, charge    ',at_name(at_ind(jj+4)),at_ind(jj+4),at_veloc_c(jj+4),at_pos_c(jj+4)
					if(data_type=='cell'.or.data_type=='quick')	then
						write(*,*) 'Cell index            ',at_ind(jj+1:jj+3)
					else
						write(*,*) 'Atom index            ',at_ind(jj+1)					
					endif
					
					write(*,*) 'Atom position            ',at_pos_c(jj+1:jj+3)
					if(index(sim_type,'static')==0.and.index(sim_type,'Static')==0.and.index(sim_type,'STATIC')==0) then
						write(*,*) 'Atom velocity             ',at_veloc_c(jj+1:jj+3)
						if(j_force.eq.1) write(*,*) 'Atom force                ',at_force_c(jj+1:jj+3)
						write(*,*)
						if(j_shell==1) then
							write(*,*) trim(at_name(at_ind(jj+4)))//'_s',at_veloc_s(jj+4),at_pos_s(jj+4)
							write(*,*) 'Shell position                ',at_pos_s(jj+1:jj+3)
							write(*,*) 'Shell velocity                ',at_veloc_s(jj+1:jj+3)
							if(j_force.eq.1) write(*,*) 'Shell force               ',at_force_s(jj+1:jj+3)
						endif
					endif				
					write(*,*)         
				enddo master_loop
				close(1)
	
				deallocate(at_name,at_occup,nsuper_r,ind_at,at_ind,at_pos_c,at_veloc_c,at_force_c)
				if(j_shell==1) deallocate(at_pos_s,at_veloc_s,at_force_s)
      enddo file_loop

      stop
      end program mp_insp53
      
     
      