!Program for calculating block-wise deformation energy score to six DNA helical 
!parameter (twist, roll, tilt, rise, shift, slide) as a function of residue no

!Reference: PNAS, 1998, 95, 1163-1168
!           JPCB, 2018, 122, 11827-11840

!The code is written by Prabir Khatua

!Date: 25th January, 2022

!Usage: f95 NucleosomeElasticForceConstBlock.F90 or gfortran NucleosomeElasticForceConstBlock.F90

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


!The results are presented in units of these variables i.e., Angstrom or 
!degree. The users can multiply or divide by KT to convert the units 
!in energy unit


!Before running this code, one has to generate DNA deformation variables for 
!each of the bases corresponding to each time frame using Curves+ and Canal
!that will generate six output files for twist, roll, tilt, rise, shift, slide

!Before running this code, one has to generate the stiffness matrix from 
!the program named NucleosomeElasticForceConstBlock.F90 

!These output files obtained from Curves+ and Canal and 
!NucleosomeElasticForceConstBlock.F90  will be used as input for this code
!------------------------------------------------------------------------------
program ElasticProperty
real(kind=8),dimension(:,:),allocatable:: x,xav
integer,dimension(:),allocatable:: typ
integer:: SHL0
real(kind=8),dimension(:),allocatable:: ETwist,ERoll,ETilt,ERise 
real(kind=8),dimension(:),allocatable:: EShift,ESlide,EDeform
real(kind=8),dimension(:,:),allocatable:: ETwistAv,ERollAv,ETiltAv
real(kind=8),dimension(:,:),allocatable:: ERiseAv,EShiftAv,ESlideAv,EDeformAv
real(kind=8),dimension(:,:,:),allocatable:: Fconst
character(len=200)::Seq,Prefix 
character(len=200):: InTwist,InRoll,InTilt,InRise,InShift,InSlide
character(len=200):: OutTwist,OutRoll,OutTilt,OutRise,OutShift,OutSlide
character(len=200):: ForceConstFile,OutDeform,IndexFile 
character(len=2),dimension(16):: SeqType
real(kind=8):: KT,En,Average,Bar
real(kind=8),dimension(0:3):: nCluster
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the input as would be asked'

write(*,'(1x,a)')'Type the prefix of all input/output files'
write(*,'(1x,a)')'This program assumes that all the required six input files for &
        & series as generated from curves+/canal has the naming style &
        & prefix_twist.ser, prefix_roll.ser, prefix_tilt.ser and so on for &
        & other parameters. If the file names are not like this, please convert the &
        & file name in this format to be compatible with this program'
read(*,'(a)')Prefix
write(*,*)Prefix

write(*,'(1x,a)')'type DNA sequence in captital letter format for A/T/G/C'
read(*,'(a)')Seq
write(*,*)Seq
write(*,'(1x,a)')'Enter the no of basepair for which the calculation to be done'
read(*,*)Nbase
write(*,*)'Nbase====>',Nbase
write(*,'(1x,a)')'Enter the SHL0 basepair no, SHL0'
read(*,*)SHL0
write(*,*)'SHL0=====>',SHL0

write(*,'(1x,a)')'Enter the no of lines in each of the input files, Nline'
read(*,*)Nline
write(*,*)'Nline====>',Nline
write(*,'(1x,a)')'Enter the temperature in Kelvin'
read(*,*)Temp
write(*,*)'Temp====>',Temp

NBlock = 3
!------------------------------------------------------------------------------
!setting the name of output files with user provided prefix

ForceConstFile = trim(prefix)//"_reswise_stiffness_matrix.dat"
IndexFile = trim(prefix)//"_cluster_index.dat"

OutTwist = trim(prefix)//"_reswise_E-twist_score_cluster.dat"
OutRoll = trim(prefix)//"_reswise_E-roll_score_cluster.dat"
OutTilt = trim(prefix)//"_reswise_E-tilt_score_cluster.dat"
OutRise = trim(prefix)//"_reswise_E-rise_score_cluster.dat"
OutShift = trim(prefix)//"_reswise_E-shift_score_cluster.dat"
OutSlide = trim(prefix)//"_reswise_E-slide_score_cluster.dat"
OutDeform = trim(prefix)//"_reswise_E-deform_score_cluster.dat"

InTwist = trim(prefix)//"_twist.ser"
InRoll = trim(prefix)//"_roll.ser"
InTilt = trim(prefix)//"_tilt.ser"
InRise = trim(prefix)//"_rise.ser"
InShift = trim(prefix)//"_shift.ser"
InSlide = trim(prefix)//"_slide.ser"

open(37,file=ForceConstFile,status='old')
open(38,file=IndexFile,status='old')

open(11,file=OutTwist,status='unknown')
open(12,file=OutRoll,status='unknown')
open(13,file=OutTilt,status='unknown')
open(14,file=OutRise,status='unknown')
open(15,file=OutShift,status='unknown')
open(16,file=OutSlide,status='unknown')
open(17,file=OutDeform,status='unknown')
!------------------------------------------------------------------------------
!memory allocation

allocate(x(Nbase,6),xav(Nbase,6),typ(Nbase),Fconst(Nbase,6,6))
allocate(ETwist(Nbase),ERoll(Nbase),ETilt(Nbase))
allocate(ERise(Nbase),EShift(Nbase),ESlide(Nbase))
allocate(EDeform(Nbase),ETwistAv(Nbase,NBlock),ERollAv(Nbase,NBlock))
allocate(ETiltAv(Nbase,NBlock),ERiseAv(Nbase,NBlock),EShiftAv(Nbase,NBlock))
allocate(ESlideAv(Nbase,NBlock),EDeformAv(Nbase,NBlock)) 
!------------------------------------------------------------------------------
call AVERAGE_VALUE(Nbase,Seq,xav,typ,SeqType)

!reading residue-wise stiffness matrix

do i=1,Nbase-1
read(37,*)
read(37,*)
do mm=1,6
read(37,'(8x,6f12.6)')(Fconst(i,mm,nn),nn=1,6)
enddo
read(37,*)
enddo

KT = 1.987*Temp/1000.d0

ETwistAv(:,:) = 0.d0; ERollAv(:,:) = 0.d0; ETiltAv(:,:) = 0.d0
ERiseAv(:,:) = 0.d0; EShiftAv(:,:) = 0.d0; ESlideAv(:,:) = 0.d0; EDeformAv(:,:) = 0.d0

open(1,file=InTwist,status='old')
open(2,file=InRoll,status='old')
open(3,file=InTilt,status='old')
open(4,file=InRise,status='old')
open(7,file=InShift,status='old')
open(8,file=InSlide,status='old')

nCluster(:) = 0.d0

do i=1,Nline

!reading from Curves+/Canal generated time series data file

read(1,*)j,(x(k,1),k=1,Nbase-1)     ! Twist
read(2,*)j,(x(k,2),k=1,Nbase-1)     ! Roll
read(3,*)j,(x(k,3),k=1,Nbase-1)     ! Tilt
read(4,*)j,(x(k,4),k=1,Nbase-1)     ! Rise
read(7,*)j,(x(k,5),k=1,Nbase-1)     ! Shift
read(8,*)j,(x(k,6),k=1,Nbase-1)     ! Slide

read(38,*)j,iBlock

if(iBlock == 0 )cycle

nCluster(iBlock) = nCluster(iBlock) + 1

EDeform(:) = 0.d0 

do m=1,Nbase-1

do mm=1,6        ! calculation of deformation energy score
do nn=1,6
En = 0.5*Fconst(m,mm,nn)*(x(m,mm) - xav(m,mm))*(x(m,nn) - xav(m,nn))
EDeform(m) = EDeform(m) + En

if(mm == 1.and.nn == 1)ETwist(m) = En
if(mm == 2.and.nn == 2)ERoll(m) = En
if(mm == 3.and.nn == 3)ETilt(m) = En
if(mm == 4.and.nn == 4)ERise(m) = En
if(mm == 5.and.nn == 5)EShift(m) = En
if(mm == 6.and.nn == 6)ESlide(m) = En

enddo
enddo

ETwistAv(m,iBlock)  = ETwistAv(m,iBlock)  + ETwist(m)
ERollAv(m,iBlock)   = ERollAv(m,iBlock)   + ERoll(m)
ETiltAv(m,iBlock)   = ETiltAv(m,iBlock)   + ETilt(m)
ERiseAv(m,iBlock)   = ERiseAv(m,iBlock)   + ERise(m)
EShiftAv(m,iBlock)  = EShiftAv(m,iBlock)  + EShift(m)
ESlideAv(m,iBlock)  = ESlideAv(m,iBlock)  + ESlide(m)
EDeformAv(m,iBlock) = EDeformAv(m,iBlock) + EDeform(m)

enddo
enddo

close(1)
close(2)
close(3)
close(4)
close(5)
close(7)
close(8)
!------------------------------------------------------------------------------
!writing residue-wise force constat data

do i=1,3
write(*,*)i,nCluster(i)
enddo

do m = 1,Nbase - 1

write(11,'(i4,3f8.3)')m,(ETwistAv(m,j)/nCluster(j),j=1,3)
write(12,'(i4,3f8.3)')m,(ERollAv(m,j)/nCluster(j),j=1,3)
write(13,'(i4,3f8.3)')m,(ETiltAv(m,j)/nCluster(j),j=1,3)
write(14,'(i4,3f8.3)')m,(ERiseAv(m,j)/nCluster(j),j=1,3)
write(15,'(i4,3f8.3)')m,(EShiftAv(m,j)/nCluster(j),j=1,3)
write(16,'(i4,3f8.3)')m,(ESlideAv(m,j)/nCluster(j),j=1,3)
write(17,'(i4,3f8.3)')m,(EDeformAv(m,j)/nCluster(j),j=1,3)

enddo

end program ElasticProperty
!------------------------------------------------------------------------------
subroutine AVERAGE_VALUE(Nb,Seq,xav,typ,SeqType)
integer,intent(in):: Nb
character(len=200),intent(in):: Seq
character(len=2),dimension(16),intent(out):: SeqType
real(kind=8),dimension(Nb,6),intent(out):: xav
real(kind=8),dimension(16,6):: RestVal
integer,dimension(*),intent(out):: typ

call REST_VALUE(RestVal,SeqType)

do i=1,Nb-1
if(Seq(i:i+1) == 'AA')then; xav(i,:) = RestVal(1,:); typ(i) = 1; endif
if(Seq(i:i+1) == 'AG')then; xav(i,:) = RestVal(2,:); typ(i) = 2; endif
if(Seq(i:i+1) == 'AT')then; xav(i,:) = RestVal(3,:); typ(i) = 3; endif
if(Seq(i:i+1) == 'AC')then; xav(i,:) = RestVal(4,:); typ(i) = 4; endif
if(Seq(i:i+1) == 'GA')then; xav(i,:) = RestVal(5,:); typ(i) = 5; endif
if(Seq(i:i+1) == 'GG')then; xav(i,:) = RestVal(6,:); typ(i) = 6; endif
if(Seq(i:i+1) == 'GT')then; xav(i,:) = RestVal(7,:); typ(i) = 7; endif
if(Seq(i:i+1) == 'GC')then; xav(i,:) = RestVal(8,:); typ(i) = 8; endif
if(Seq(i:i+1) == 'TA')then; xav(i,:) = RestVal(9,:); typ(i) = 9; endif
if(Seq(i:i+1) == 'TG')then; xav(i,:) = RestVal(10,:); typ(i) = 10; endif
if(Seq(i:i+1) == 'TT')then; xav(i,:) = RestVal(11,:); typ(i) = 11; endif
if(Seq(i:i+1) == 'TC')then; xav(i,:) = RestVal(12,:); typ(i) = 12; endif
if(Seq(i:i+1) == 'CA')then; xav(i,:) = RestVal(13,:); typ(i) = 13; endif
if(Seq(i:i+1) == 'CG')then; xav(i,:) = RestVal(14,:); typ(i) = 14; endif
if(Seq(i:i+1) == 'CT')then; xav(i,:) = RestVal(15,:); typ(i) = 15; endif
if(Seq(i:i+1) == 'CC')then; xav(i,:) = RestVal(16,:); typ(i) = 16; endif

enddo

end subroutine AVERAGE_VALUE
!------------------------------------------------------------------------------
subroutine REST_VALUE(x,Seq)
real(kind=8),dimension(16,6),intent(out):: x
character(len=2),dimension(16),intent(out):: Seq

!These are rest-state or equlilibrium value for DNA step helical parameter 
!as obtained from averaging the data over crystal structures corresponding to 
!many protein-DNA complex. 

!These values are taken from Ref. PNAS, 1998, 95, 1163-1168

!   Twist          Roll            Tilt            Rise             Shift            Slide

x(1,1)  = 35.1; x(1,2)  = 0.7; x(1,3)  = -1.4; x(1,4)  =  3.27; x(1,5)  = -0.03; x(1,6)  =  -0.08;  Seq(1)  = 'AA'
x(2,1)  = 31.9; x(2,2)  = 4.5; x(2,3)  = -1.7; x(2,4)  =  3.34; x(2,5)  =  0.09; x(2,6)  =  -0.25;  Seq(2)  = 'AG'
x(3,1)  = 29.3; x(3,2)  = 1.1; x(3,3)  =  0.0; x(3,4)  =  3.31; x(3,5)  =  0.00; x(3,6)  =  -0.59;  Seq(3)  = 'AT'
x(4,1)  = 31.5; x(4,2)  = 0.7; x(4,3)  = -0.1; x(4,4)  =  3.36; x(4,5)  =  0.13; x(4,6)  =  -0.58;  Seq(4)  = 'AC'
x(5,1)  = 36.3; x(5,2)  = 1.9; x(5,3)  = -1.5; x(5,4)  =  3.37; x(5,5)  = -0.28; x(5,6)  =   0.09;  Seq(5)  = 'GA'
x(6,1)  = 32.9; x(6,2)  = 3.6; x(6,3)  = -0.1; x(6,4)  =  3.42; x(6,5)  =  0.05; x(6,6)  =  -0.22;  Seq(6)  = 'GG'
x(7,1)  = 31.5; x(7,2)  = 0.7; x(7,3)  =  0.1; x(7,4)  =  3.36; x(7,5)  = -0.13; x(7,6)  =  -0.58;  Seq(7)  = 'GT'
x(8,1)  = 33.6; x(8,2)  = 0.3; x(8,3)  =  0.0; x(8,4)  =  3.40; x(8,5)  =  0.00; x(8,6)  =  -0.38;  Seq(8)  = 'GC'
x(9,1)  = 37.8; x(9,2)  = 3.3; x(9,3)  =  0.0; x(9,4)  =  3.42; x(9,5)  =  0.00; x(9,6)  =   0.05;  Seq(9)  = 'TA'
x(10,1) = 37.3; x(10,2) = 4.7; x(10,3) = -0.5; x(10,4) =  3.33; x(10,5) = -0.09; x(10,6) =   0.53;  Seq(10) = 'TG'
x(11,1) = 35.1; x(11,2) = 0.7; x(11,3) =  1.4; x(11,4) =  3.27; x(11,5) =  0.03; x(11,6) =  -0.08;  Seq(11) = 'TT'
x(12,1) = 36.3; x(12,2) = 1.9; x(12,3) =  1.5; x(12,4) =  3.37; x(12,5) =  0.28; x(12,6) =   0.09;  Seq(12) = 'TC'
x(13,1) = 37.3; x(13,2) = 4.7; x(13,3) =  0.5; x(13,4) =  3.33; x(13,5) =  0.09; x(13,6) =   0.53;  Seq(13) = 'CA'
x(14,1) = 36.1; x(14,2) = 5.4; x(14,3) =  0.0; x(14,4) =  3.39; x(14,5) =  0.00; x(14,6) =   0.41;  Seq(14) = 'CG'
x(15,1) = 31.9; x(15,2) = 4.5; x(15,3) =  1.7; x(15,4) =  3.34; x(15,5) = -0.09; x(15,6) =  -0.25;  Seq(15) = 'CT'
x(16,1) = 32.9; x(16,2) = 3.6; x(16,3) =  0.1; x(16,4) =  3.42; x(16,5) = -0.05; x(16,6) =  -0.22;  Seq(16) = 'CC'

end subroutine REST_VALUE
!----------------------end of the program--------------------------------------
