!Program for calculating native contact

!The code is written by Prabir Khatua

!Original code was written for analyzing SAA at OU on: 23rd May, 2019
!The code is modified for nucleosome on: 12th October, 2021


!Usage: f95 DNA_HISTONE_NATIVE_CONTACT.F90 or gfortran DNA_HISTONE_NATIVE_CONTACT.F90
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

!Output file will give the time evolution of fraction of native contacts

!The code requires a reference pdb file that will be used to identify the 
!native contact pairs. The no of atoms in the reference pdb file and trajectory
!dcd file should be same. 

!Definition of native conatcts has been adapted from PNAS, 2013 110 (44)
!17874-17879 

!The value of native contacts depends on the values of smoothing parameter 
!beta and lambda. Generally, one uses the value of beta as 5 Angs^-1. 
!lambda is generally considered as 1.2 and 1.8 for Go and all-atom model, 
!respectively. 

!The unit of beta is in Angs^-1. Therefore, the distance has to be measured in
!Angs unit. As this is gromacs generated trajectories, which writes coordinates
!in nm unit, appropriate unit conversion has been done inside the code. 

!The code is written for analyzing dcd formatted trajectory
!------------------------------------------------------------------------------
!variable declaration

program CONTACT
real,dimension(:),allocatable:: x,y,z
real(kind=8):: dx,dy,dz,nat_con,beta,lambda
real(kind=8),dimension(:),allocatable:: nat_dist
integer,dimension(:),allocatable:: dna_atom_id,histone_atom_id
integer,dimension(:,:),allocatable:: pair
integer:: ounit,runit
character(len=100):: outfile,reffile
character(len=100),dimension(:),allocatable:: dcdfile
character(len=4):: chr2,chr3
character(len=2):: chr4
logical:: there
data ounit,inunit,runit /12,14,16/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Type the name of output file'
read(*,'(a)')outfile
write(*,*)outfile
open(ounit,file=outfile,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------------------------------------
!input information about protein to be analysed

write(*,'(1x,a)')'provide the following protein related inputs as would be asked'

write(*,'(1x,a)')'Enter the no of total atoms, natom'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter the first and last atom no for DNA, dna_first, dna_last'
read(*,*)ifirst_dna,ilast_dna
write(*,*)'ifirst_dna,ilast_dna=====>',ifirst_dna,ilast_dna
write(*,'(1x,a)')'Enter the first and last atom no for histone, histone_first, histone_last'
read(*,*)ifirst_histone,ilast_histone
write(*,'(1x,a)')'Enter the total no of bases present in DNA'
read(*,*)num_dna_base
write(*,*)'num_dna_base======>',num_dna_base

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

!input information about the trajectory to be analysed'

write(*,'(1x,a)')'provide the following trajectory information as would be asked'
write(*,'(1x,a)')'Choose from the following options'
write(*,'(1x,a)')'Type 0 if the trajectory is saved as wrapped coordinates'
write(*,'(1x,a)')'Type 1 if the trajectory is not saved as wrapped coordinates'
read(*,*)iwrap
write(*,*)'You have selected the option',iwrap
write(*,'(1x,a)')'Type the timestep between two frames of the trajectory in ps'
read(*,*)dt
write(*,*)'dt===>',dt
write(*,'(1x,a)')'Enter the skip value'
read(*,*)nskip
write(*,*)'nskip======>',nskip

!input information related to contact map calculation

write(*,'(1x,a)')&
&'Note that generally all heavy atoms are considered for native contact &
&calculation. However, the present code allows users to calculate native &
& contact for different options as would be choosen. '

write(*,'(1x,a)')'Type the cut-off distance to define the contact in nm, rc'
read(*,*)rc
write(*,*)'rc===>',rc
write(*,'(1x,a)')'Choose from the following options'
write(*,'(1x,a)')'Total DNA system: 0'
write(*,'(1x,a)')'Tail-1: 1'
write(*,'(1x,a)')'Tail-2: 2'
write(*,'(1x,a)')'Total DNA except Tail-1 & Tail-2: 3'
read(*,*)ind
write(*,*)'You have selected the option', ind

rc=rc*rc

!input information related to native contact calculation

write(*,'(1x,a)')'Enter the value of beta for smoothing (normally 5 A^-1)'
read(*,*)beta
write(*,*)'beta===>',beta
write(*,'(1x,a)')'Enter the value of lambda for smoothing (normally 1.8)'
read(*,*)lambda
write(*,*)'lambda====>',lambda
!------------------------------------------------------------------------------
!memory allocation

allocate(x(natom),y(natom),z(natom))
allocate(dna_atom_id(natom),histone_atom_id(natom))
!------------------------------------------------------------------------------
!identifying the native contacts from the reference file

write(*,'(1x,a)')'Type the name of reference pdb file that would be &
&required to define the native contact pair'
read(*,'(a)')reffile
write(*,*)reffile
inquire(file=reffile,exist=there)

if(.not.there)then
write(*,'(a)')'ERROR: REFERENCE PDB FILE IS NOT FOUND!!'
stop
endif

open(runit,file=reffile,status='old',form='formatted',iostat=ios)

np_dna=0                   ! np_dna = total # of selected dna atoms
np_histone=0               ! np_histone = total # of selected histone atoms

do i=1,natom
read(runit,111)chr2,chr3,nres,x(i),y(i),z(i),chr4
111   format(12x,a4,1x,a4,1x,i4,4x,3f8.3,22x,a2)

if(chr4 == ' H')cycle

if(i >= ifirst_histone.and.i <= ilast_histone)then
        np_histone=np_histone+1
        histone_atom_id(np_histone)=i
endif


if(i >= ifirst_dna.and.i <= ilast_dna)then
    if(IND == 0)then
    np_dna=np_dna+1
    dna_atom_id(np_dna)=i
    endif 

    if(IND == 1)then
            if((nres >= ifirst_tail1_strand1.and.nres <= ilast_tail1_strand1)&
            &.or.(nres >= ifirst_tail1_strand2.and.&
            &nres <= ilast_tail1_strand2))then
            np_dna=np_dna+1
            dna_atom_id(np_dna)=i
            endif
    endif

    if(IND == 2)then
            if((nres >= ifirst_tail2_strand1.and.nres <= ilast_tail2_strand1)&
            &.or.(nres >= ifirst_tail2_strand2.and.&
            &nres <= ilast_tail2_strand2))then
            np_dna=np_dna+1
            dna_atom_id(np_dna)=i
            endif
    endif

    if(IND == 3)then
            if((nres >= ifirst_tail1_strand1.and.nres <= ilast_tail1_strand1)&
            &.or.(nres >= ifirst_tail1_strand2.and.nres <= ilast_tail1_strand2)&
            &.or.(nres >= ifirst_tail2_strand1.and.nres <= ilast_tail2_strand1)&
            &.or.(nres >= ifirst_tail2_strand2.and.nres &
            &<= ilast_tail2_strand2))cycle

            np_dna=np_dna+1
            dna_atom_id(np_dna)=i
    endif
endif
enddo

close(runit)

write(*,*)'# OF SELECTED DNA HEAVY ATOM====>',np_dna
write(*,*)'# OF SELECTED HISTONE HEAVY ATOM====>',np_histone
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
write(*,'(a)')'ERROR: TRAJECTORY FILE NOT FOUND!!'
stop
endif
enddo
!------------------------------------------------------------------------------
nread=0
nused=0
open(inunit,file=dcdfile(1),status='old',form='unformatted')
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
write(*,*)'nused=====>',nused

if(nused == 1)then
call CONTACT_PAIR(np_dna,np_histone,dna_atom_id,histone_atom_id,x,y,z,dx,dy,dz,iwrap,rc,npair)

write(*,*)'no of atom pairs forming native contact===>',npair

nd=npair
allocate(pair(nd,2),nat_dist(nd))

call PAIR_IDENTITY(np_dna,np_histone,dna_atom_id,histone_atom_id,x,y,z,dx,dy,dz,iwrap,rc,nat_dist,nd,pair,npair)

write(*,*)'# of total native contact atom pair====>',npair
exit
endif
enddo
close(inunit)
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
call NATIVE(pair,nd,npair,x,y,z,nat_dist,dx,dy,dz,iwrap,beta,lambda,nat_con)

time=time+dt
write(ounit,'(2f12.6)')time/1000.0,nat_con
enddo
close(inunit)
enddo

end program CONTACT
!------------------------------------------------------------------------------
subroutine NATIVE(pair,nd,npair,x,y,z,rn,dx,dy,dz,iwrap,beta,lambda,q_all)
real,dimension(*),intent(in):: x,y,z
integer,intent(in):: npair,nd,iwrap
integer,dimension(nd,2),intent(in):: pair
real(kind=8),intent(out):: q_all
real(kind=8),dimension(*),intent(in):: rn
real(kind=8),intent(in):: dx,dy,dz,beta,lambda
real(kind=8):: rr

q_all=0.d0
do ii=1,npair
i=pair(ii,1)
j=pair(ii,2)

rx=x(i)-x(j)
ry=y(i)-y(j)
rz=z(i)-z(j)

if(iwrap == 1)then
rx=rx-dx*anint(rx/dx)
ry=ry-dy*anint(ry/dy)
rz=rz-dz*anint(rz/dz)
endif

r=rx*rx+ry*ry+rz*rz
r=sqrt(r)
rr=1.0+exp(beta*10.0*(r-lambda*rn(ii)))


q_all=q_all+1.0/rr
enddo

q_all=q_all/dble(npair)

end subroutine NATIVE
!------------------------------------------------------------------------------
subroutine PAIR_IDENTITY(np_dna,np_histone,dna_atom_id,histone_atom_id,x,y,z,dx,dy,dz,iwrap,rc,rn,nd,pair,npair)
real,dimension(*),intent(in):: x,y,z
integer,intent(in):: np_dna,np_histone,iwrap
integer,dimension(*),intent(in):: dna_atom_id,histone_atom_id
real,intent(in):: rc
integer,intent(in):: nd
integer,dimension(nd,2),intent(out):: pair
integer,intent(out):: npair
real(kind=8),dimension(*),intent(out):: rn
real(kind=8),intent(in):: dx,dy,dz

npair=0
do ii=1,np_dna
i=dna_atom_id(ii)

do jj=1,np_histone
j=histone_atom_id(jj)

rx=x(i)-x(j)
ry=y(i)-y(j)
rz=z(i)-z(j)

if(iwrap == 1)then
rx=rx-dx*anint(rx/dx)
ry=ry-dy*anint(ry/dy)
rz=rz-dz*anint(rz/dz)
endif

r=rx*rx+ry*ry+rz*rz
if(r > rc)cycle

npair=npair+1
rn(npair)=sqrt(r)
pair(npair,1)=i
pair(npair,2)=j

enddo
enddo

end subroutine PAIR_IDENTITY
!------------------------------------------------------------------------------
subroutine CONTACT_PAIR(np_dna,np_histone,dna_atom_id,histone_atom_id,x,y,z,dx,dy,dz,iwrap,rc,npair)
real,dimension(*),intent(in):: x,y,z
integer,intent(in):: np_dna,np_histone,iwrap
integer,dimension(*),intent(in):: dna_atom_id,histone_atom_id
real,intent(in):: rc
integer,intent(inout):: npair
real(kind=8),intent(in):: dx,dy,dz

npair=0
do ii=1,np_dna
i=dna_atom_id(ii)

do jj=1,np_histone
j=histone_atom_id(jj)

rx=x(i)-x(j)
ry=y(i)-y(j)
rz=z(i)-z(j)

if(iwrap == 1)then
rx=rx-dx*anint(rx/dx)
ry=ry-dy*anint(ry/dy)
rz=rz-dz*anint(rz/dz)
endif

r=rx*rx+ry*ry+rz*rz
if(r <=rc)npair=npair+1
enddo
enddo

end subroutine CONTACT_PAIR
!------------------------------------------------------------------------------

