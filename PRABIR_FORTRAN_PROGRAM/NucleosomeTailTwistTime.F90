!Program for calculating twist angle of the nucleosomal tail as 
!a function of time

!The code is written by Prabir Khatua

!Date: 12th November, 2021

!Usage: f95 NucleosomeTailTwistTime.F90 or gfortran NucleosomeTailTwistTime.F90

!       ./a.out  (for interactive run and then enter the inputs interactivly 
!                 as the program will ask)
!       
!       ./a.out<input.inp >output.log &  
!                            (input.inp is the input file containing all the 
!                            required inputs that has to be prepared by the 
!                            user. For this, one needs to open the code or 
!                            better compile the code once in interactive mode 
!                            to know what set of inputs will be required and 
!                            prepare the input file. output.log will print 
!                            general information or error) 

!The code is written to analyse the dcd trajectories
!------------------------------------------------------------------------------
program OrderParameter
real,dimension(:),allocatable:: x,y,z
real(kind=8),dimension(:),allocatable:: mass
integer,dimension(:),allocatable:: BaseIDTail1,BaseIDTail2
integer,dimension(:),allocatable:: SenseDnaResNatom,AntiSenseDnaResNatom,ResID
integer,dimension(:,:),allocatable:: SenseDnaAtomId,AntiSenseDnaAtomId
integer:: dummyi,ounit,punit,NumDnaBase,NumDnaBasePair
character(len=100),dimension(:),allocatable:: dcdfile
character(len=100):: pdbfile,outfile
character(len=4),dimension(:),allocatable:: AtomName,ResName
character(len=1):: ElementName
logical:: there,ANS
real(kind=8):: dx,dy,dz,bx,by,bz
real(kind=8):: time,dt,TwistTail1,TwistTail2
data inunit,punit,ounit /10,11,12/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of output file for contact map'
read(*,'(a)')outfile
write(*,*)outfile
open(ounit,file=outfile,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the no of total atoms, natom'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter the first and last atom no for DNA, dna_first, dna_last'
read(*,*)ifirst_dna,ilast_dna
write(*,*)'ifirst_dna,ilast_dna=====>',ifirst_dna,ilast_dna
write(*,'(1x,a)')'Enter the total no of bases present in DNA'
read(*,*)NumDnaBase
write(*,*)'NumDnaBase====>',NumDnaBase
write(*,'(1x,a)')'Choose from the following options'
write(*,'(1x,a)')'Type 0 if the trajectory is saved as wrapped coordinates'
write(*,'(1x,a)')'Type 1 if the trajectory is not saved as wrapped coordinates'
read(*,*)iwrap
write(*,*)'You have selected the option',iwrap
write(*,'(1x,a)')'Enter the timestep in ps, dt'
read(*,*)dt
write(*,*)'dt=====>',dt
write(*,'(1x,a)')'Enter the skip value, nskip'
read(*,*)nskip
write(*,*)'nskip====>',nskip

NumDnaBasePair = NumDnaBase/2
write(*,*)'# of Nucleosomal DNA base pair===>',NumDnaBasePair

nHISTONE = ilast_histone - ifirst_histone + 1
nDNA = ilast_DNA - ifirst_DNA + 1

write(*,*)'# of DNA atoms========>',nDNA
write(*,*)'# of Histone atoms====>',nHISTONE

allocate(x(natom),y(natom),z(natom),mass(natom))
allocate(SenseDnaResNatom(NumDnaBasePair),AntiSenseDnaResNatom(NumDnaBasePair))
allocate(SenseDnaAtomId(NumDnaBasePair,100))
allocate(AntiSenseDnaAtomId(NumDnaBasePair,100))
allocate(AtomName(natom),ResName(natom),ResID(natom))
allocate(BaseIDTail1(4),BaseIDTail2(4))
!------------------------------------------------------------------------------
write(*,'(1x,a)')'type 4 base index for tail1 twisting on the sense strand'
read(*,*)(BaseIDTail1(i),i=1,4)
write(*,*)(BaseIDTail1(i),i=1,4)
write(*,*)

write(*,'(1x,a)')'type 4 base index for tail2 twisting on the sense strand'
read(*,*)(BaseIDTail2(i),i=1,4)
write(*,*)(BaseIDTail2(i),i=1,4)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of pdbfile'
read(*,'(a)')pdbfile
write(*,*)pdbfile

inquire(file=pdbfile,exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: pdbfile NOT FOUND'
stop
endif

open(punit,file=pdbfile,status='old',form='formatted',iostat=ios)

SenseDnaResNatom(:) = 0
AntiSenseDnaResNatom(:) = 0

do i=1,natom
read(punit,111)AtomName(i),ResName(i),ResID(i),x(i),y(i),z(i)
111   format(12x,a4,1x,a4,1x,i4,4x,3f8.3)

ElementName=adjustl(AtomName(i))


if(ElementName == 'H')mass(i) = 1.0080
if(ElementName == 'C')mass(i) = 12.0110
if(ElementName == 'N')mass(i) = 14.0070
if(ElementName == 'O')mass(i) = 15.9990
if(ElementName == 'S')mass(i) = 32.0600
if(ElementName == 'P')mass(i) = 30.9740

if(ElementName == 'H')cycle

if(i >= ifirst_dna.and.i <= ilast_dna)then
        if(ResID(i) <= NumDnaBasePair)then
                SenseDnaResNatom(ResID(i)) = SenseDnaResNatom(ResID(i)) + 1
                SenseDnaAtomId(ResID(i),SenseDnaResNatom(ResID(i))) = i
        else
                IDBasePair = NumDnaBase + 1 - ResID(i)
                AntiSenseDnaResNatom(IDBasePair) = AntiSenseDnaResNatom(IDBasePair) + 1
                AntiSenseDnaAtomId(IDBasePair,AntiSenseDnaResNatom(IDBasePair)) = i
        endif
endif
enddo
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
write(*,'(a)')'ERROR: dcdfile NOT FOUND'
stop
endif
enddo
!------------------------------------------------------------------------------
nread=0
nused=0

time=0.d0
dt=dt*nskip

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
Call TWIST_ANGLE(BaseIDTail1,NumDnaBasePair,SenseDnaResNatom,SenseDnaAtomId,&
        &AntiSenseDnaResNatom,AntiSenseDnaAtomId,mass,x,y,z,&
        &iwrap,dx,dy,dz,TwistTail1)

Call TWIST_ANGLE(BaseIDTail2,NumDnaBasePair,SenseDnaResNatom,SenseDnaAtomId,&
        &AntiSenseDnaResNatom,AntiSenseDnaAtomId,mass,x,y,z,&
        &iwrap,dx,dy,dz,TwistTail2)

time = time + dt
write(ounit,'(3f12.6)')time/1000000.0,abs(TwistTail1),abs(TwistTail2)
enddo
enddo

end program OrderParameter
!------------------------------------------------------------------------------
subroutine TWIST_ANGLE(BaseID,NumDnaBasePair,SenseDnaResNatom,SenseDnaAtomId,&
        &AntiSenseDnaResNatom,AntiSenseDnaAtomId,mass,x,y,z,&
        &iwrap,BoxX,BoxY,BoxZ,Angle)
integer,intent(in):: NumDnaBasePair,iwrap
integer,dimension(*),intent(in):: BaseID,SenseDnaResNatom
integer,dimension(*),intent(in):: AntiSenseDnaResNatom
integer,dimension(NumDnaBasePair,*),intent(in):: SenseDnaAtomId
integer,dimension(NumDnaBasePair,*),intent(in):: AntiSenseDnaAtomId
real,dimension(*),intent(in):: x,y,z
real(kind=8),dimension(*),intent(in):: mass
real(kind=8),intent(in):: BoxX,BoxY,BoxZ
real(kind=8):: SenseBaseComX,SenseBaseComY,SenseBaseComZ
real(kind=8):: AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ
real(kind=8),dimension(4):: XX,YY,ZZ
real(kind=8),intent(out):: Angle

do ii=1,4

i = BaseID(ii)

call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,x,y,z,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &x,y,z,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

XX(ii) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
YY(ii) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
ZZ(ii) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)
enddo

call DIHEDRAL(XX,YY,ZZ,iwrap,BoxX,BoxY,BoxZ,Angle)

!if(Angle < 0.0)then
!        Angle = 360 + Angle
!else
!        Angle = Angle
!endif

end subroutine TWIST_ANGLE
!------------------------------------------------------------------------------
!subroutine for calculating dihedral angle
!------------------------------------------------------------------------------
subroutine DIHEDRAL(x,y,z,iwrap,bx,by,bz,phi)
real(kind=8),dimension(*),intent(in):: x,y,z
real(kind=8),dimension(3):: rx,ry,rz
real(kind=8),dimension(2):: rxc,ryc,rzc
real(kind=8),intent(out):: phi
real(kind=8):: pi,s,a,b,c,ss
real(kind=8):: vx,vy,vz,dx,dy,dz
real(kind=8),intent(in):: bx,by,bz
integer,intent(in):: iwrap

pi=4.d0*datan(1.d0)

do i=1,3
j = i + 1

rx(i) = x(j) - x(i)
ry(i) = y(j) - y(i)
rz(i) = z(j) - z(i)

if(iwrap == 1)then
rx(i) = rx(i) - anint(rx(i)/bx)*bx
ry(i) = ry(i) - anint(ry(i)/by)*by
rz(i) = rz(i) - anint(rz(i)/bz)*bz
endif

enddo

do i=1,2

j = i + 1

vx = rx(i)
vy = ry(i)
vz = rz(i)
dx = rx(j)
dy = ry(j)
dz = rz(j)

call CROSS_PROD(vx,vy,vz,dx,dy,dz,rxc(i),ryc(i),rzc(i))
enddo

call DOT_PROD(rxc(1),ryc(1),rzc(1),rxc(2),ryc(2),rzc(2),s)
call CROSS_PROD(rxc(1),ryc(1),rzc(1),rxc(2),ryc(2),rzc(2),a,b,c)
call DOT_PROD(rx(2),ry(2),rz(2),a,b,c,ss)   !the value of this dot product
                                            ! will decide whether the
                                            ! didheral angle will be (+)/(-)ve

phi = dacos(s)
phi = 180*phi/pi
if(ss < 0.0)phi = -phi

end subroutine DIHEDRAL
!------------------------------------------------------------------------------
subroutine CROSS_PROD(ax,ay,az,bx,by,bz,cx,cy,cz)   ! C = A X B
real(kind=8),intent(in):: ax,ay,az,bx,by,bz
real(kind=8),intent(out):: cx,cy,cz

cx = ay*bz - by*az                        
cy = bx*az - ax*bz
cz = ax*by - bx*ay

end subroutine CROSS_PROD
!------------------------------------------------------------------------------
subroutine DOT_PROD(aX,ay,aZ,bX,bY,bZ,Angle)   ! A.B = |A|*|B| cos(Angle)
real(kind=8),intent(in):: aX,ay,aZ,bX,bY,bZ    !     = aX*bX + aY*bY + aZ*BZ
real(kind=8),intent(out):: Angle
real(kind=8):: DotX,DotY,DotZ,ModA,ModB

DotX = aX*bX
DotY = aY*bY
DotZ = aZ*bZ

ModA = dsqrt(aX*aX + aY*aY + aZ*aZ)
ModB = dsqrt(bX*bX + bY*bY + bZ*bZ)

Angle = DotX + DotY + DotZ
Angle = Angle/(ModA*ModB)

end subroutine DOT_PROD
!------------------------------------------------------------------------------
subroutine UNIT_VECTOR(aX,aY,aZ,bX,bY,bZ,abX,abY,abZ)  ! AB = B - A
real(kind=8),intent(in):: aX,aY,aZ,bX,bY,bZ
real(kind=8),intent(out):: abX,abY,abZ
real(kind=8):: Modab

abX = bX - aX
abY = bY - aY
abZ = bZ - aZ

Modab = sqrt(abX*abX + abY*abY + abZ*abZ)

abX=abX/Modab
abY=abY/Modab
abZ=abZ/Modab

end subroutine UNIT_VECTOR
!------------------------------------------------------------------------------
subroutine DNA_BASE_COM(i,Nbase,ResNatom,AtomId,mass,x,y,z,ComX,ComY,ComZ)
integer,intent(in):: i,Nbase,ResNatom
integer,dimension(Nbase,*),intent(in):: AtomId
real,dimension(*),intent(in):: x,y,z
real(kind=8),dimension(*),intent(in):: mass
real(kind=8),intent(out):: ComX,ComY,ComZ
real(kind=8):: TotalMass

TotalMass = 0.d0
ComX = 0.d0
ComY = 0.d0 
ComZ = 0.d0

do j=1,ResNatom

ii = AtomId(i,j)

ComX = ComX + x(ii)*mass(ii)
ComY = ComY + y(ii)*mass(ii)
ComZ = ComZ + z(ii)*mass(ii)

TotalMass = TotalMass + mass(ii)

enddo

ComX = ComX/TotalMass
ComY = ComY/TotalMass
ComZ = ComZ/TotalMass

end subroutine DNA_BASE_COM
!------------------------------------------------------------------------------
subroutine COM(Npair,AtomId,mass,x,y,z,ComX,ComY,ComZ)
real,dimension(*),intent(in):: x,y,z
real(kind=8),dimension(*),intent(in):: mass
integer,dimension(*),intent(in):: AtomId
integer,intent(in):: Npair
real(kind=8),intent(out):: ComX,ComY,ComZ
real(kind=8):: TotalMass

TotalMass = 0.d0
ComX = 0.d0
ComY = 0.d0
ComZ = 0.d0

do i=1,Npair

ii = AtomId(i)

ComX = ComX + x(ii)*mass(ii)
ComY = ComY + y(ii)*mass(ii)
ComZ = ComZ + z(ii)*mass(ii)

TotalMass = TotalMass + mass(ii)

enddo

ComX = ComX/TotalMass
ComY = ComY/TotalMass
ComZ = ComZ/TotalMass

end subroutine COM
!-----------------------end of the program-------------------------------------
