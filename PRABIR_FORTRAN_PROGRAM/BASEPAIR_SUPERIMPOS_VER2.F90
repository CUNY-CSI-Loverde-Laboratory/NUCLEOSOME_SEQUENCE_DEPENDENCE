!program for calculating RMSD by quarterion method

!The code is written by Prabir Khatua
!date: 21st November, 2014
!Modified on: 5th October, 2021

!------------------------------------------------------------------------------
program rmsd
real,dimension(:),allocatable:: x,y,z,xr,yr,zr,xm
integer,dimension(:),allocatable:: atp,phos_atom
integer:: dummyi,ounit,punit,base_pair,base_pair_2,t1unit,t2unit
character(len=100),dimension(:),allocatable:: dcdfile
character(len=100):: pdbfile,tail1,tail2
character(len=4),dimension(:),allocatable:: chr2,chr3
character(len=1),dimension(:),allocatable:: chr4
logical:: there,ANS
real(kind=8):: dx,dy,dz,bx,by,bz,cx,cy,cz
data inunit,punit,t1unit,t2unit /10,11,12,13/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of output file for tail1'
read(*,'(a)')tail1
write(*,*)tail1
open(t1unit,file=tail1,status='unknown',form='formatted',iostat=ios)

write(*,'(1x,a)')'Enter the name of output file for tail2'
read(*,'(a)')tail2
write(*,*)tail2
open(t2unit,file=tail2,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the no of total atoms'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter the skip value'
read(*,*)nskip
write(*,*)'nskip======>',nskip
write(*,'(1x,a)')'Enter the first and last atom nno u want to analyse'
read(*,*)ifirst,ilast
write(*,*)'ifirst,ilast=====>',ifirst,ilast

!Memory allocation

allocate(x(natom),y(natom),z(natom),xr(natom),yr(natom),zr(natom))
allocate(atp(natom),xm(natom),phos_atom(natom))
allocate(chr2(natom),chr3(natom),chr4(natom))
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of pdbfile'
read(*,'(a)')pdbfile
write(*,*)pdbfile
open(punit,file=pdbfile,status='old',form='formatted',iostat=ios)

np=0
base_pair=0
do i=1,ilast ! natom
read(punit,111)chr2(i),chr3(i),nres,xr(i),yr(i),zr(i),chr4(i)
111   format(12x,a4,1x,a4,1x,i4,4x,3f8.3,22x,a2)

chr4(i) =  adjustl(chr2(i))

if(chr4(i) == 'H')xm(i)=1.0080
if(chr4(i) == 'C')xm(i)=12.0110
if(chr4(i) == 'N')xm(i)=14.0070
if(chr4(i) == 'O')xm(i)=15.9990
if(chr4(i) == 'S')xm(i)=32.0600
if(chr4(i) == 'P')xm(i)=30.9740

if(chr4(i) == 'H')cycle

!if(chr3(i) == 'DA5'.or.chr3(i) == 'DA3'.or.&
!        &chr3(i) == 'DT5'.or.chr3(i) == 'DT3'.or.&
!        &chr3(i) == 'DG5'.or.chr3(i) == 'DG3'.or.&
!        &chr3(i) == 'DC5'.or.chr3(i) == 'DC3')cycle

if(i >= ifirst.and.i <= ilast)then
np=np+1
atp(np)=i

if(chr2(i) ==' P  ')then
        base_pair=base_pair+1
        phos_atom(base_pair)=i
endif
endif
enddo

write(*,*)'# OF SELECTED HEAVY ATOM====>',np
write(*,*)'# OF BASES====>',base_pair

base_pair_2=base_pair/2

call CENTER(natom,np,atp,xr,yr,zr,xm)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the number of trajectory files you want to analyze.'
read(*,*) ninput
write(*,*)'ninput==========>',ninput

allocate(dcdfile(ninput))

do inp=1,ninput
write(*,*)'Enter the name of trajectory file',inp
read(*,'(a)')dcdfile(inp)
write(*,*)dcdfile(inp)
inquire(file=dcdfile(inp),exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: CONFIGURATION FILE NOT FOUND'
stop
endif
enddo
!------------------------------------------------------------------------------
nused=0
nread=0
do inp=1,ninput
write(*,'(1x,a,i5,a)')'Enter name of trajectory file',inp
write(*,*)dcdfile(inp)
open(inunit,file=dcdfile(inp),status='old',form='unformatted')
read(inunit)dummyc,nframes,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
read(inunit)dummyi,dummyr
write(*,*)"dummyi,dummyr====>",dummyi,dummyr
read(inunit)natom
write(*,*)"natom====>",natom
write(*,*)"nframes===>",nframes
do iframe=1,nframes
read(inunit)dx,bx,dy,by,bz,dz
read(inunit)(x(j),j=1,natom)
read(inunit)(y(j),j=1,natom)
read(inunit)(z(j),j=1,natom)
nread=nread+1
if(mod(nread,nskip) /= 0)cycle
nused=nused+1
!------------------------------------------------------------------------------
call CENTER(natom,np,atp,x,y,z,xm)
call ALIGN(natom,np,atp,x,y,z,xr,yr,zr)

do i=1,base_pair_2
j=base_pair+1-i
ii=phos_atom(i)
jj=phos_atom(j)


!Self-checking

!if(nused == 1)then
!if((chr3(ii) == ' DA '.and.chr3(jj) /= ' DT ').or.&
!        (chr3(ii) == ' DT '.and.chr3(jj) /= ' DA ').or.&
!       &(chr3(ii) == ' DG '.and.chr3(jj) /= ' DC ').or.&
!       &(chr3(ii) == ' DC '.and.chr3(jj) /= ' DG '))then
!       write(*,*)'BASE PAIRING IS NOT RIGHT!'
!       write(*,*)'SOMETHING IS WRONG'
!       write(*,*)i,ii,jj,chr3(ii),chr3(jj)
!       stop

!endif
!endif

cx=(x(ii)+x(jj))/2.0
cy=(y(ii)+y(jj))/2.0
cz=(z(ii)+z(jj))/2.0

if(i <= base_pair_2/2)then
        write(t1unit,'(3f8.3)')cx,cy,cz
else
        write(t2unit,'(3f8.3)')cx,cy,cz
endif

enddo

write(t1unit,*)
write(t2unit,*)

enddo
enddo

end program RMSD
!------------------------------------------------------------------------------
subroutine ALIGN(npr,np,atp,xa,ya,za,xr,yr,zr)
real,dimension(*),intent(inout):: xa,ya,za
real,dimension(*),intent(in):: xr,yr,zr
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,npr
real,dimension(4):: q
real,dimension(4,4):: m

m(1,1)=0.0
m(1,2)=0.0
m(1,3)=0.0
m(1,4)=0.0
m(2,2)=0.0
m(2,3)=0.0
m(2,4)=0.0
m(3,3)=0.0
m(3,4)=0.0
m(4,4)=0.0
do j=1,np
i=atp(j)
m(1,1)=m(1,1)+xa(i)*xr(i)+ya(i)*yr(i)+za(i)*zr(i)
m(2,2)=m(2,2)+xa(i)*xr(i)-ya(i)*yr(i)-za(i)*zr(i)
m(3,3)=m(3,3)+ya(i)*yr(i)-za(i)*zr(i)-xa(i)*xr(i)
m(4,4)=m(4,4)+za(i)*zr(i)-xa(i)*xr(i)-ya(i)*yr(i)
m(1,2)=m(1,2)+ya(i)*zr(i)-za(i)*yr(i)
m(1,3)=m(1,3)+za(i)*xr(i)-xa(i)*zr(i)
m(1,4)=m(1,4)+xa(i)*yr(i)-ya(i)*xr(i)
m(2,3)=m(2,3)+xa(i)*yr(i)+ya(i)*xr(i)
m(2,4)=m(2,4)+za(i)*xr(i)+xa(i)*zr(i)
m(3,4)=m(3,4)+ya(i)*zr(i)+za(i)*yr(i)
end do
m(2,1)=m(1,2)
m(3,1)=m(1,3)
m(4,1)=m(1,4)
m(3,2)=m(2,3)
m(4,2)=m(2,4)
m(4,3)=m(3,4)

call JACOBI(q,m,4)
call ROTN(xa,ya,za,q,npr)

end subroutine ALIGN
!------------------------------------------------------------------------------
subroutine ROTN(xa,ya,za,q,n)
integer,intent(in):: n
real,dimension(*),intent(inout):: xa,ya,za
real,dimension(*),intent(in):: q
real,dimension(n):: a,b,c
real:: q0,q2,qd,qv,q1

do i=1,n
a(i)=xa(i)
b(i)=ya(i)
c(i)=za(i)
end do
q0=q(1)
q1=2*q0
q0=q0**2
q2=q(2)**2+q(3)**2+q(4)**2
qd=q0-q2
do i=1,n
qv=q(2)*a(i)+q(3)*b(i)+q(4)*c(i)
qv=2*qv
xa(i)=qd*a(i)+qv*q(2)+q1*(q(3)*c(i)-q(4)*b(i))
ya(i)=qd*b(i)+qv*q(3)+q1*(q(4)*a(i)-q(2)*c(i))
za(i)=qd*c(i)+qv*q(4)+q1*(q(2)*b(i)-q(3)*a(i))
end do

end subroutine ROTN
!------------------------------------------------------------------------------
SUBROUTINE JACOBI(Q,A,N)
INTEGER::I,J,IP,IQ,M,N
INTEGER::NROT,SWEEP
REAL,DIMENSION(1:4,1:4),INTENT(INOUT)::A
REAL,DIMENSION(1:4,1:4)::V
REAL,DIMENSION(1:4),INTENT(OUT)::Q
REAL,DIMENSION(1:4)::D(4),B(4),Z(4)
REAL::THETA,TAU,TRESH,C,S,T,G,H,SM,DMX
!
!     ----Jacobi Method---
!     ----Ref: Numerical Recipe in Fortran; page 460 
!     ----Initialization of eigenvector martix as unit matrix---
DO IP=1,N
DO IQ=1,N
V(IP,IQ)=0.0
END DO
V(IP,IP)=1.0
END DO
!     ----Initialization of eigenvector as the diagonal elements---
DO IP=1,N
B(IP)=A(IP,IP)
D(IP)=B(IP)
Z(IP)=0.0
END DO
!     ----Iteration---
NROT=0
SWEEP=0
112  SWEEP=SWEEP+1
IF(SWEEP>150)THEN
PRINT*,"PROGRAM STOP: TOO MANY ITERATION IN THE CALCULATION."
STOP
END IF
SM=0.0                                     !Sum of absolute value
DO IP=1,N-1                                !of super off-diagonal
DO IQ=IP+1,N                             !elements. Program will
SM=SM+ABS(A(IP,IQ))                      !stop if the sum become
END DO                                   !zero.
END DO
IF(SM==0.0)GOTO 132
!     -----------------------------------------------------------------
IF(SWEEP<4)THEN                         !Read the paragraph in
TRESH=0.2*SM/N**2                        !Numerical Recipes in
ELSE                                     !Fortan, on the page 
TRESH=0.0                                !No. 459
END IF
!     -----------------------------------------------------------------
DO IP=1,N-1
DO IQ=IP+1,N
!        ----------calculation in a sweep-----
G=100.0*ABS(A(IP,IQ))
!After 4 step skip the rotation if of diagonal element is small
IF((SWEEP>4).AND.(ABS(D(IP))+G==ABS(D(IP))).AND.(ABS(D(IQ))+G==ABS(D(IQ))))THEN
A(IP,IQ)=0.0
ELSE IF(ABS(A(IP,IQ))>TRESH)THEN    
H=D(IQ)-D(IP)                        
IF(ABS(H)+G==ABS(H))THEN           
T=A(IP,IQ)/H
ELSE                                     
THETA=0.5*H/A(IP,IQ)               !Rotation parameters calculation
T=1.0/(ABS(THETA)+SQRT(1.0+ THETA**2))  
IF(THETA<0.0)T=-T               
END IF                               
C=1.0/SQRT(1.0+T**2)                 
S=T*C                                
TAU=S/(1.0+C)
!     ----operation corresponding rotation-----------------------------
H=T*A(IP,IQ)
Z(IP)=Z(IP)-H
Z(IQ)=Z(IQ)+H
D(IP)=D(IP)-H
D(IQ)=D(IQ)+H
A(IP,IQ)=0.0
!          -----------------------
DO J=1,IP-1                          !Changes in ip coulumn
G=A(J,IP)                          !and iq coulumn upto
H=A(J,IQ)                          !ip row. Stop before
A(J,IP)=G-S*(H+G*TAU)              !ip row. 
A(J,IQ)=H+S*(G-H*TAU)              !
END DO                               !
!           -----------------------
DO J=IP+1,IQ-1                       !Change in ip row from
G=A(IP,J)                          !ip coulumn to iq coulu
H=A(J,IQ)                          !mn and change in iq co
A(IP,J)=G-S*(H+G*TAU)              !ulumn from ip row to 
A(J,IQ)=H+S*(G-H*TAU)              !iq row.
END DO
!           -----------------------
DO J=IQ+1,N                          !change in ip row from
G=A(IP,J)                          !iq coulumn to end coul
H=A(IQ,J)                          !umn and change in iq 
A(IP,J)=G-S*(H+G*TAU)              !row from iq coulumn 
A(IQ,J)=H+S*(G-H*TAU)              !to end coulumn.
END DO                               !
!           -----------------------
DO J=1,N                             !
G=V(J,IP)                          !Update of the eigen-
H=V(J,IQ)                          !vector matrix.
V(J,IP)=G-S*(H+G*TAU)              !
V(J,IQ)=H+S*(G-H*TAU)              !
END DO                               !
!           -----------------------
NROT=NROT+1
END IF
END DO
END DO
!     ----------------
DO IP=1,N
B(IP)=B(IP)+Z(IP)                          !Update of d with sum 
D(IP)=B(IP)                                !of ta(p,q) and reinit
Z(IP)=0.0                                  !ialization of z
end do
!     ----------------
SWEEP=SWEEP+1
GOTO 112
132   DMX=D(1)
M=1
DO IP=1,N
IF(D(IP)>DMX)THEN
DMX=D(IP)
M=IP
END IF
END DO
DO IP=1,N
Q(IP)=V(IP,M)
END DO
END SUBROUTINE JACOBI
!------------------------------------------------------------------------------
subroutine CENTER(npr,np,atp,x,y,z,xm)
real,dimension(*),intent(inout):: x,y,z
real,dimension(*),intent(in):: xm
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,npr

cx=0.0
cy=0.0
cz=0.0
tm=0.0
do i=1,np
ii=atp(i)
cx=cx+xm(ii)*x(ii)
cy=cy+xm(ii)*y(ii)
cz=cz+xm(ii)*z(ii)
tm=tm+xm(ii)
enddo

cx=cx/tm
cy=cy/tm
cz=cz/tm

do i=1,npr
x(i)=x(i)-cx
y(i)=y(i)-cy
z(i)=z(i)-cz
enddo

end subroutine CENTER
!------------------------------------------------------------------------------
