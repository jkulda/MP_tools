      program genlat
!
!***********************************************************************
!
!     program to generate a perfect lattice with general unit cell
!
!     copyright daresbury laboratory
!     author - w.smith oct. 1992
!
!     revised for core-shell models by Jiri Kulda 2024
!***********************************************************************
!
  integer,parameter :: mxatom = 1000
  integer :: i,j,k,m,n,j_shell,nx,ny,nz,natms
  real :: cell(9),cprp(10),xb(mxatom),yb(mxatom),zb(mxatom),xs,ys,zs,xx,yy,zz,wdth
  character*8 :: name(mxatom),name_s(mxatom)
  character*80 :: title


  write(*,*)'enter suitable one-line title for lattice file'
  read(*,'(a80)')title
  write(*,*)'enter number of basis atoms (max 1000)'
  read(*,*)natms
  write(*,*)'define unit cell vector A(1-3)'
  read(*,*)cell(1),cell(2),cell(3)
  write(*,*)'define unit cell vector B(1-3)'
  read(*,*)cell(4),cell(5),cell(6)
  write(*,*)'define unit cell vector C(1-3)'
  read(*,*)cell(7),cell(8),cell(9)
  write(*,*)'define 3 unit cell multiplications LxA,MxB,NxC'
  read(*,*)nx,ny,nz
  write(*,*)'core-shell system? (1/0)'
  read(*,*) j_shell
  if(j_shell==0) then
    write(*,*)'enter ion name (<=8 chars) and fractional coordinates'
    do i=1,natms
       read(*,*)name(i),xb(i),yb(i),zb(i)
    enddo
  else
    write(*,*)'enter core & shell names (8 chars) and fractional coordinates'
    do i=1,natms
       read(*,*)name(i),name_s(i),xb(i),yb(i),zb(i)
    enddo
  endif

  open(10,file='LATTICE')
  write(10,'(a80)')title
  write(10,'(2i10)')0,3
  write(10,'(3f20.8)')dble(nx)*cell(1),dble(nx)*cell(2),dble(nx)*cell(3),dble(ny)*cell(4),dble(ny)*cell(5),dble(ny)*cell(6),&
&                     dble(nz)*cell(7),dble(nz)*cell(8),dble(nz)*cell(9)
!
!     set up lattice
  m=0
  do  n=1,natms
    do  k=1,nz
       do  j=1,ny
          do  i=1,nx
            m=m+1
            xs=dble(i-1)+xb(n)-0.5d0*dble(nx)
            ys=dble(j-1)+yb(n)-0.5d0*dble(ny)
            zs=dble(k-1)+zb(n)-0.5d0*dble(nz)
            xx=cell(1)*xs+cell(4)*ys+cell(7)*zs
            yy=cell(2)*xs+cell(5)*ys+cell(8)*zs
            zz=cell(3)*xs+cell(6)*ys+cell(9)*zs
            write(10,'(a8,i10,/,3f20.6)')name(n),m,xx,yy,zz
            if(j_shell==1) then
              m=m+1
              write(10,'(a8,i10,/,3f20.6)')name_s(n),m,xx,yy,zz
            endif
           enddo
        enddo
     enddo
  enddo
  call dcell(cell,cprp)
  wdth=0.5d0*min(dble(nx)*cprp(7),dble(ny)*cprp(8),dble(nz)*cprp(9))
  write(*,*)'number of ions in system = ',m
  write(*,*)'maximum radius of cutoff = ',wdth
  stop
  end




  subroutine dcell(aaa,bbb)
!
!***********************************************************************
!
!     dl_poly subroutine to calculate the dimensional properies of
!     a simulation cell specified by the input matrix aaa.
!     the results are returned in the array bbb, with :
!
!     bbb(1 to 3) - lengths of cell vectors
!     bbb(4 to 6) - cosines of cell angles
!     bbb(7 to 9) - perpendicular cell widths
!     bbb(10)     - cell volume
!
!     copyright daresbury laboratory
!     author - w. smith july 1992
!
!***********************************************************************
!

      real :: aaa(9),bbb(10),axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3
!
!     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))
!
!     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))
!
!     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)
!
!     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)
!
!     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end
