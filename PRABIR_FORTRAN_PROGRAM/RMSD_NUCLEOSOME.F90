!program for calculating RMSD by quarterion method

!The code is written by Prabir Khatua
!date: 21st November, 2014
!Modified on: 5th October, 2021

!------------------------------------------------------------------------------
program rmsd
real,dimension(:),allocatable:: x,y,z,xr,yr,zr,xm
integer,dimension(:),allocatable:: atp
integer:: dummyi,ounit,punit
character(len=100),dimension(:),allocatable:: dcdfile
character(len=100):: pdbfile,outfile
character(len=4):: chr2,chr3
character(len=2):: chr4
logical:: there,ANS
real(kind=8):: dx,dy,dz,bx,by,bz
data inunit,punit,ounit /10,11,12/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of output file'
read(*,'(a)')outfile
write(*,*)outfile
open(ounit,file=outfile,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the no of total atoms'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter DNA first and last atom nno u want to analyse'
read(*,*)ifirst,ilast
write(*,*)'ifirst,ilast=====>',ifirst,ilast
write(*,'(1x,a)')'Enter the total no of bases present in DNA'
read(*,*)num_dna_base
write(*,*)'num_dna_base======>',num_dna_base
write(*,'(1x,a)')'Choose from following options'
write(*,'(1x,a)')'Total DNA system: 0'
write(*,'(1x,a)')'Tail-1: 1'
write(*,'(1x,a)')'Tail-2: 2'
write(*,'(1x,a)')'Total DNA except Tail-1 & Tail-2: 3'
read(*,*)ind
write(*,*)'You have selected the option',ind
write(*,'(1x,a)')'Enter the timesteps in ps'
read(*,*)dt
write(*,*)'dt======>',dt
write(*,'(1x,a)')'Enter the skip value'
read(*,*)nskip
write(*,*)'nskip======>',nskip

num_dna_base_pair=num_dna_base/2

write(*,*)'# of DNA base pair====>',num_dna_base_pair

ifirst_tail1_strand1=1
ilast_tail1_strand1=10

ilast_tail1_strand2=num_dna_base
ifirst_tail1_strand2=num_dna_base-9

write(*,*)'Base pair number range for nucleosomal DNA tail-1'
write(*,*)'Strand-1:',ifirst_tail1_strand1,'to',ilast_tail1_strand1
write(*,*)'Strand-2:',ifirst_tail1_strand2,'to',ilast_tail1_strand2

ifirst_tail2_strand1=num_dna_base_pair-9
ilast_tail2_strand1=num_dna_base_pair

ilast_tail2_strand2=num_dna_base-num_dna_base_pair+10
ifirst_tail2_strand2=num_dna_base-num_dna_base_pair+1

write(*,*)'Base pair number range for nucleosomal DNA tail-2'
write(*,*)'Strand-1:',ifirst_tail2_strand1,'to',ilast_tail2_strand1
write(*,*)'Strand-2:',ifirst_tail2_strand2,'to',ilast_tail2_strand2

!Memory allocation

allocate(x(natom),y(natom),z(natom),xr(natom),yr(natom),zr(natom))
allocate(atp(natom),xm(natom))
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of pdbfile'
read(*,'(a)')pdbfile
write(*,*)pdbfile
open(punit,file=pdbfile,status='old',form='formatted',iostat=ios)

np=0
do i=1,natom
read(punit,111)chr2,chr3,nres,xr(i),yr(i),zr(i),chr4
111   format(12x,a4,1x,a4,1x,i4,4x,3f8.3,22x,a2)
if(chr4 == ' H')xm(i)=1.0080
if(chr4 == ' C')xm(i)=12.0110
if(chr4 == ' N')xm(i)=14.0070
if(chr4 == ' O')xm(i)=15.9990
if(chr4 == ' S')xm(i)=32.0600
if(chr4 == ' P')xm(i)=30.9740

if(chr4 == ' H')cycle

if(i < ifirst.or.i > ilast)cycle

if(IND == 0)then
        np=np+1
        atp(np)=i
endif

if(IND == 1)then
        if((nres >= ifirst_tail1_strand1.and.nres <= ilast_tail1_strand1).or.&
        &(nres >= ifirst_tail1_strand2.and.nres <= ilast_tail1_strand2))then
        np=np+1
        atp(np)=i
        endif
endif

if(IND == 2)then
        if((nres >= ifirst_tail2_strand1.and.nres <= ilast_tail2_strand1).or.&
        &(nres >= ifirst_tail2_strand2.and.nres <= ilast_tail2_strand2))then
        np=np+1
        atp(np)=i
        endif
endif

if(IND == 3)then
        if((nres >= ifirst_tail1_strand1.and.nres <= ilast_tail1_strand1).or.&
        &(nres >= ifirst_tail1_strand2.and.nres <= ilast_tail1_strand2).or.&
        &(nres >= ifirst_tail2_strand1.and.nres <= ilast_tail2_strand1).or.&
        &(nres >= ifirst_tail2_strand2.and.nres <= ilast_tail2_strand2))cycle

np=np+1
atp(np)=i
endif

enddo

write(*,*)'# OF SELECTED HEAVY ATOM====>',np
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
nu=0
nr=0
time=0.0
dt=dt*nskip
do inp=1,ninput
write(*,*)'Enter name of trajectory file',inp
write(*,*)dcdfile(inp)
open(inunit,file=dcdfile(inp),status='old',form='unformatted')
read(inunit)dummyc,nframes,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
read(inunit)dummyi,dummyr
write(*,*)"dummyi,dummyr====>",dummyi,dummyr
read(inunit)natom
write(*,*)"natom====>",natom
write(*,*)"nframes===>",nframes
do ii=1,nframes
read(inunit)dx,bx,dy,by,bz,dz
read(inunit)(x(j),j=1,natom)
read(inunit)(y(j),j=1,natom)
read(inunit)(z(j),j=1,natom)
nr=nr+1
if(mod(nr,nskip) /= 0)cycle
nu=nu+1
!------------------------------------------------------------------------------
time=time+dt
call CENTER(natom,np,atp,x,y,z,xm)
call ALIGN(natom,np,atp,x,y,z,xr,yr,zr)

rm=0.0
do i=1,np
j=atp(i)
rx=x(j)-xr(j)
ry=y(j)-yr(j)
rz=z(j)-zr(j)
rm=rm+rx*rx+ry*ry+rz*rz
enddo

rm=rm/real(np)
rm=sqrt(rm)
write(ounit,'(2f12.6)')time/1000.0,rm
        
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
real,dimension(*),intent(inout):: xa,ya,za
real,dimension(*),intent(in):: q
real,dimension(40000):: a,b,c
real:: q0,q2,qd,qv,q1
integer,intent(in):: n

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
