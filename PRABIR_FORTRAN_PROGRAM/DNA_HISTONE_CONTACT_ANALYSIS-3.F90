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
integer,dimension(:),allocatable:: HistoneCoreID,NTailID,CTailID
integer,dimension(:),allocatable:: ResID,DnaAtomId
!integer,dimension(:),allocatable:: NativeCore,NativeNTail,NativeCTail
integer:: dummyi,punit,Total,ArgContact,LysContact,CTAIL
character(len=200),dimension(:),allocatable:: dcdfile
character(len=200):: pdbfile,Prefix
character(len=300):: CoreFile,NTailFile,CTailFile
character(len=4),dimension(:),allocatable:: AtomName,ResName
character(len=1):: ElementName
logical:: there,ANS
real(kind=8):: dx,dy,dz,bx,by,bz
data inunit,punit /10,11/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the prefix for all output files'
read(*,'(a)')Prefix
write(*,*)Prefix

CoreFile = trim(Prefix)//"_core_contact_all_type_time.out"
NTailFile = trim(Prefix)//"_Ntail_contact_all_type_time.out"
CTailFile = trim(Prefix)//"_Ctail_contact_all_type_time.out"

open(12,file=CoreFile,status='unknown')
open(13,file=NTailFile,status='unknown')
open(14,file=CTailFile,status='unknown')
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
write(*,'(1x,a)')'Enter first and last residue no for the Histone subunit core region'
read(*,*)ifirst_histone_core,ilast_histone_core
write(*,*)'ifirst_histone_core,ilast_histone_core====>',ifirst_histone_core,ilast_histone_core
write(*,'(1x,a)')'Enter first and last residue no for the Histone subunit Ntail region'
read(*,*)ifirst_NTail,ilast_NTail
write(*,*)'ifirst_NTail,ilast_NTail====>',ifirst_NTail,ilast_NTail
write(*,'(1x,a)')'Type 1 if there is C-terminal tail or else 0'
read(*,*)CTAIL
write(*,*)CTAIL

if(CTAIL == 1)then
        write(*,'(1x,a)')'Enter first and last residue no for the Histone subunit Ctail region'
        read(*,*)ifirst_CTail,ilast_CTail
        write(*,*)'ifirst_CTail,ilast_CTail====>',ifirst_CTail,ilast_CTail
endif

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
write(*,'(1x,a)')'Enter the cut-off distance for contact'
read(*,*)rc
write(*,*)'rc====>',rc

rc = rc * rc

ifirst_antisense = NumDnaBase - ilast_sense + 1
ilast_antisense = NumDnaBase - ifirst_sense + 1

write(*,*)'Antisense bp no',ifirst_antisense,ilast_antisense

allocate(x(natom),y(natom),z(natom),mass(natom))
allocate(DnaAtomId(natom),HistoneCoreID(natom))
allocate(CTailID(natom),NTailID(natom))
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
NpairCore = 0
NpairNTail = 0
NpairCTail = 0

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

if(ResID(i) >= ifirst_histone_core.and. ResID(i) <= ilast_histone_core)then
        NpairCore = NpairCore + 1
        HistoneCoreID(NpairCore) = i
endif

if(ResID(i) >= ifirst_NTail.and. ResID(i) <= ilast_NTail)then
        NpairNTail = NpairNTail + 1
        NTailID(NpairNTail) = i
endif

if(CTAIL == 1)then
        if(ResID(i) >= ifirst_CTail.and. ResID(i) <= ilast_CTail)then
                NpairCTail = NpairCTail + 1
                CTailID(NpairCTail) = i
        endif
endif

enddo

write(*,*)'# of selected DNA atom===>',Npair
write(*,*)'# of selected core-Histone atom===>',NpairCore
write(*,*)'# of selected N-terminal Histone atom===>',NpairNTail
write(*,*)'# of selected C-terminal Histone atom===>',NpairCTail

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
call CONTACT(Npair,NpairCore,ResName,DnaAtomID,HistoneCoreID,x,y,z,iwrap,dx,dy,dz,rc,Total,ArgContact,LysContact)

write(12,'(f12.6,1x,4i10)')time/1000000.0,LysContact,ArgContact,LysContact+ArgContact,Total

call CONTACT(Npair,NpairNTail,ResName,DnaAtomID,NTailID,x,y,z,iwrap,dx,dy,dz,rc,Total,ArgContact,LysContact)

write(13,'(f12.6,1x,4i10)')time/1000000.0,LysContact,ArgContact,LysContact+ArgContact,Total

if(CTAIL == 1)then
        call CONTACT(Npair,NpairCTail,ResName,DnaAtomID,CTailID,x,y,z,iwrap,dx,dy,dz,rc,Total,ArgContact,LysContact)
        write(14,'(f12.6,1x,4i10)')time/1000000.0,LysContact,ArgContact,LysContact+ArgContact,Total
endif

enddo
enddo

end program DNA_HISTONE_CONTACT
!------------------------------------------------------------------------------
subroutine CONTACT(NP1,NP2,ResName,ID1,ID2,x,y,z,iwrap,dx,dy,dz,rc,Total,ArgContact,LysContact)
integer,intent(in):: NP1,NP2,iwrap
integer,dimension(*),intent(in):: ID1,ID2
character(len=4),dimension(*),intent(in):: ResName
real,dimension(*),intent(in):: x,y,z
real,intent(in):: rc
integer,intent(out):: Total,ArgContact,LysContact
real(kind=8),intent(in):: dx,dy,dz

Total = 0
ArgContact = 0 
LysContact = 0

do ii = 1,NP1
i = ID1(ii)
do jj = 1,NP2
j = ID2(jj)
nCounter = nCounter + 1

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

Total = Total + 1
if(ResName(j) == 'ARG ')ArgContact = ArgContact + 1
if(ResName(j) == 'LYS ')LysContact = LysContact + 1

enddo
enddo

end subroutine CONTACT
!------------------------------------------------------------------------------
