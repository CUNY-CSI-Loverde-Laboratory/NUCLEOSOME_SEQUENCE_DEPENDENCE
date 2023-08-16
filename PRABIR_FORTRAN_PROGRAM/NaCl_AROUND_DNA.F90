!Program for calculating # of DNA-Histone contact as a function of time 
!(total, native, and non-native) for core region and non-core regions

!The code is written by Prabir Khatua

!Date: 1st July, 2022

!Usage: f95 DNA_HISTONE_CONTACT_ANALYSIS.F90 or gfortran DNA_HISTONE_CONTACT_ANALYSIS.F90

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
program DNA_HISTONE_CONTACT
real,dimension(:),allocatable:: x,y,z,mass
integer,dimension(:),allocatable:: ResID,DnaAtomId,NaId,ClId
integer:: dummyi,punit,Total,NumNa,NumCl
character(len=200),dimension(:),allocatable:: dcdfile
character(len=200):: pdbfile
character(len=300):: OutFile
character(len=4),dimension(:),allocatable:: AtomName,ResName
character(len=1):: ElementName
logical:: there,ANS
real(kind=8):: dx,dy,dz,bx,by,bz
data inunit,punit /10,11/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of output file'
read(*,'(a)')OutFile
write(*,*)OutFile

open(12,file=OutFile,status='unknown')
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the no of total atoms, natom'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter 1st and last bp no for sense strand'
read(*,*)ifirst_sense,ilast_sense
write(*,*)'ifirst_sense,ilast_sense=====>',ifirst_sense,ilast_sense
write(*,'(1x,a)')'Enter the total no of bases present in DNA'
read(*,*)NumDnaBase
write(*,*)'NumDnaBase====>',NumDnaBase
write(*,'(1x,a)')'Enter the first and last Na+ ion atom id'
read(*,*)ifirst_na,ilast_na
write(*,*)'ifirst_na,ilast_na====>',ifirst_na,ilast_na
write(*,'(1x,a)')'Enter the first and last Cl- ion atom id'
read(*,*)ifirst_cl,ilast_cl
write(*,*)'ifirst_cl,ilast_cl====>',ifirst_cl,ilast_cl
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
write(*,'(1x,a)')'Enter the cut-off distance'
read(*,*)rc
write(*,*)'rc====>',rc

rc = rc * rc

ifirst_antisense = NumDnaBase - ilast_sense + 1
ilast_antisense = NumDnaBase - ifirst_sense + 1

N_na = ilast_na - ifirst_na + 1
N_cl = ilast_cl - ifirst_cl + 1

write(*,*)'N_na,N_cl====>',N_na,N_cl

write(*,*)'Antisense bp no',ifirst_antisense,ilast_antisense

allocate(x(natom),y(natom),z(natom),mass(natom))
allocate(DnaAtomId(natom),NaId(N_na),ClId(N_cl))
allocate(AtomName(natom),ResName(natom),ResID(natom))
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

Npair = 0
NpairNa = 0
NpairCl = 0

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

if((ResID(i) >= ifirst_sense.and.ResID(i) <= ilast_sense).or.&
        &(ResID(i) >= ifirst_antisense.and.ResID(i) <= ilast_antisense))then

        Npair = Npair + 1
        DnaAtomID(Npair) = i
endif

if(i >= ifirst_na.and.i <= ilast_na)then
        if(AtomName(i) == 'Na+ ')then
                NpairNa = NpairNa + 1
                NaId(NpairNa) = i
        endif
endif

if(i >= ifirst_cl.and.i <= ilast_cl)then
        if(AtomName(i) == 'Cl- ')then
                NpairCl = NpairCl + 1
                ClId(NpairCl) = i
        endif
endif

enddo

write(*,*)'# of selected DNA atom===>',Npair
write(*,*)'# of Na+ ions===>',NpairNa
write(*,*)'# of Cl- ions===>',NpairCl
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
nread = 0
nused = 0

time = 0.d0
dt = dt*nskip

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
time = time + dt
call IONS(Npair,NpairNa,DnaAtomID,NaId,x,y,z,iwrap,dx,dy,dz,rc,NumNa)
call IONS(Npair,NpairCl,DnaAtomID,ClId,x,y,z,iwrap,dx,dy,dz,rc,NumCl)

write(12,'(f12.6,1x,2i10)')time/1000000.0,NumNa,NumCl

enddo
enddo

end program DNA_HISTONE_CONTACT
!------------------------------------------------------------------------------
subroutine IONS(NP1,NP2,ID1,ID2,x,y,z,iwrap,dx,dy,dz,rc,Num)
integer,intent(in):: NP1,NP2,iwrap
integer,dimension(*),intent(in):: ID1,ID2
real,dimension(*),intent(in):: x,y,z
real,intent(in):: rc
integer,intent(out):: Num
real(kind=8),intent(in):: dx,dy,dz

Num = 0
do ii = 1,NP1
i = ID1(ii)
do jj = 1,NP2
j = ID2(jj)

rx = x(i) - x(j)
ry = y(i) - y(j)
rz = z(i) - z(j)

if(iwrap == 1)then
rx = rx - dx * anint(rx/dx)
ry = ry - dy * anint(ry/dy)
rz = rz - dz * anint(rz/dz)
endif

r = rx * rx + ry * ry + rz * rz

if(r > rc)cycle

Num = Num + 1

enddo
enddo

end subroutine IONS
!------------------------------------------------------------------------------
