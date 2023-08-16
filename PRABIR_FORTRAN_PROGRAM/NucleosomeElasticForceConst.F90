!Program for calculating force constant corresponding to six DNA helical 
!parameter (twist, roll, tilt, rise, shift, slide) as a function of time

!Reference: PNAS, 1998, 95, 1163-1168
!           JPCB, 2018, 122, 11827-11840

!The code is written by Prabir Khatua

!Date: 20th December, 2021
!Finalized on: 18th January, 2022
!Modified on: 28th January, 2022

!Usage: f95 NucleosomeElasticForceConst.F90 or gfortran NucleosomeElasticForceConst.F90

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

!Stiffness Matrix F is calculated from the inversion of covariance matrix V. 
!V is 6 X 6 covariance matrix calculated from the DNA
!deformation variable (twist, roll, tilt, rise, shift, and slide). 
!The element of V(i,j) = <del_theta(i)del_theta(j)>; where 
!del_theta(i) = theta(i) - theta(i)_0; theta(i) are these six deformation 
!variable. The reference rest-state theta(i)_0 values are taken from 
!above mentioned PNAS paper

!The results are presented in units of these variables i.e., Angstrom or 
!degree. The users can multiply or divide by KT to convert the units 
!in energy unit

!In this code the inverse matrix is calculated using Cholesky Decomposition

!Before running this code, one has to generate DNA deformation variables for 
!each of the bases corresponding to each time frame using Curves+ and Canal
!that will generate six output files for twist, roll, tilt, rise, shift, slide

!These output files obtained from Curves+ and Canal will be used as input 
!for this code

!There will be a series of output files such as residue-wise force constant; 
!time evolution of these force constats for each of the dinucleotide steps; 
!stiffness matrix corresponding to each of the base-steps; scattered data 
!of these six variables for 16 different possible dinucleotide steps that 
!users can use for further analysis
!------------------------------------------------------------------------------
program ElasticProperty
real(kind=8),dimension(:,:),allocatable:: x,xav
integer,dimension(:),allocatable:: typ
integer:: SHL0
real(kind=8),dimension(:),allocatable:: Twist,Roll,Tilt,Rise,Shift,Slide
real(kind=8),dimension(6,6):: Vinv
real(kind=8),dimension(:,:,:),allocatable:: V,Fconst
character(len=200)::Seq,Prefix 
character(len=200):: InTwist,InRoll,InTilt,InRise,InShift,InSlide
character(len=200):: OutRight,OutLeft
character(len=200):: ResAvgOut,ForceConstFile 
character(len=200),dimension(16):: ScatterFile
character(len=2),dimension(16):: SeqType
real(kind=8):: Percentage,KT,En
real(kind=8):: TwistAv,RollAv,TiltAv,RiseAv,ShiftAv,SlideAv
real(kind=8):: TwistBar,RollBar,TiltBar,RiseBar,ShiftBar,SlideBar
logical:: SUCCESS
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
!------------------------------------------------------------------------------
!setting the name of output files with user provided prefix

ForceConstFile = trim(prefix)//"_reswise_stiffness_matrix.dat"

ResAvgOut = trim(prefix)//"_reswise_force-const.dat"

OutRight = trim(prefix)//"_reswise_force_const_right_time.dat"
OutLeft = trim(prefix)//"_reswise_force_const_left_time.dat"

InTwist = trim(prefix)//"_twist.ser"
InRoll = trim(prefix)//"_roll.ser"
InTilt = trim(prefix)//"_tilt.ser"
InRise = trim(prefix)//"_rise.ser"
InShift = trim(prefix)//"_shift.ser"
InSlide = trim(prefix)//"_slide.ser"

ScatterFile(1) = trim(prefix)//"_AA_helical_para_scatter_plot.dat"
ScatterFile(2) = trim(prefix)//"_AG_helical_para_scatter_plot.dat"
ScatterFile(3) = trim(prefix)//"_AT_helical_para_scatter_plot.dat"
ScatterFile(4) = trim(prefix)//"_AC_helical_para_scatter_plot.dat"
ScatterFile(5) = trim(prefix)//"_GA_helical_para_scatter_plot.dat"
ScatterFile(6) = trim(prefix)//"_GG_helical_para_scatter_plot.dat"
ScatterFile(7) = trim(prefix)//"_GT_helical_para_scatter_plot.dat"
ScatterFile(8) = trim(prefix)//"_GC_helical_para_scatter_plot.dat"
ScatterFile(9) = trim(prefix)//"_TA_helical_para_scatter_plot.dat"
ScatterFile(10) = trim(prefix)//"_TG_helical_para_scatter_plot.dat"
ScatterFile(11) = trim(prefix)//"_TT_helical_para_scatter_plot.dat"
ScatterFile(12) = trim(prefix)//"_TC_helical_para_scatter_plot.dat"
ScatterFile(13) = trim(prefix)//"_CA_helical_para_scatter_plot.dat"
ScatterFile(14) = trim(prefix)//"_CG_helical_para_scatter_plot.dat"
ScatterFile(15) = trim(prefix)//"_CT_helical_para_scatter_plot.dat"
ScatterFile(16) = trim(prefix)//"_CC_helical_para_scatter_plot.dat"

open(10,file=ResAvgOut,status='unknown')
open(37,file=ForceConstFile,status='unknown')

open(11,file=OutRight,status='unknown')
open(12,file=OutLeft,status='unknown')

open(1,file=InTwist,status='old')
open(2,file=InRoll,status='old')
open(3,file=InTilt,status='old')
open(4,file=InRise,status='old')
open(7,file=InShift,status='old')
open(8,file=InSlide,status='old')

iunit = 19
do i=1,16
iunit = iunit + 1
open(iunit,file=ScatterFile(i),status='unknown')
enddo

write(11,100)
write(12,100)
100 format('   frame  bp      Twist        Roll        Tilt       Rise        Shift       Slide')
!------------------------------------------------------------------------------
!memory allocation

allocate(x(Nbase,6),xav(Nbase,6),typ(Nbase))
allocate(V(Nbase,6,6),Fconst(Nbase,6,6))
allocate(Twist(Nbase),Roll(Nbase),Tilt(Nbase))
allocate(Rise(Nbase),Shift(Nbase),Slide(Nbase))
!------------------------------------------------------------------------------
call AVERAGE_VALUE(Nbase,Seq,xav,typ,SeqType)

V(:,:,:) = 0.d0
KT = 1.987*Temp/1000.d0

do i=1,Nline

!reading from Curves+/Canal generated time series data file

read(1,*)j,(x(k,1),k=1,Nbase-1)     ! Twist
read(2,*)j,(x(k,2),k=1,Nbase-1)     ! Roll
read(3,*)j,(x(k,3),k=1,Nbase-1)     ! Tilt
read(4,*)j,(x(k,4),k=1,Nbase-1)     ! Rise
read(7,*)j,(x(k,5),k=1,Nbase-1)     ! Shift
read(8,*)j,(x(k,6),k=1,Nbase-1)     ! Slide

do m=1,Nbase-1

!writing scattered data for 16 different dinucleotide steps

if(typ(m) == 1)write(20,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 2)write(21,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 3)write(22,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 4)write(23,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 5)write(24,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 6)write(25,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 7)write(26,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 8)write(27,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 9)write(28,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 10)write(29,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 11)write(30,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 12)write(31,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 13)write(32,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 14)write(33,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 15)write(34,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)
if(typ(m) == 16)write(35,'(i12,1x,i4,1x,6f8.3)')i,m,(x(m,n),n=1,6)

do mm=1,6        !calculation of cavariance matrix
do nn=1,6
En = (x(m,mm) - xav(m,mm))*(x(m,nn) - xav(m,nn))
V(m,mm,nn)=V(m,mm,nn) + En 
enddo
enddo

if(i <= 50)cycle

do mm=1,6
do nn=1,6
Vinv(mm,nn) = V(m,mm,nn)/dble(i)
enddo
enddo

call CHOLDC_INVERSE(Vinv,6,SUCCESS) ! Matrix inversion
!call LU_DECOMP(Vinv,6,SUCCESS)    ! Matrix inversion

if(.not.SUCCESS)cycle

do mm=1,6 ! storing the stiffness matrix Fconst for 
do nn=1,6 ! each of the residues
Fconst(m,mm,nn) = Vinv(mm,nn) 
enddo
enddo

Twist(m) = Vinv(1,1)
Roll(m)  = Vinv(2,2)
Tilt(m)  = Vinv(3,3)
Rise(m)  = Vinv(4,4)
Shift(m) = Vinv(5,5)
Slide(m) = Vinv(6,6)

enddo

if(i <= 50)cycle

!writing the time-series of residue-wise force constat 
!as obtained from matrix inversion

do m=2,SHL0
mm = m - SHL0
write(11,98)i,mm,Twist(m),Roll(m),Tilt(m),Rise(m),Shift(m),Slide(m)
98 format(i7,1x,i4,1x,6f12.6)
enddo
write(11,*)

do m=SHL0,Nbase-1
mm = m - SHL0
write(12,98)i,mm,Twist(m),Roll(m),Tilt(m),Rise(m),Shift(m),Slide(m)
enddo
write(12,*)

enddo

!writing the residue-wise stiffness matrix 

do i=1,Nbase-1
write(37,'(a13,i4)')'Base Step = ',i
write(37,*)
do mm=1,6
write(37,'(8x,6f12.6)')(Fconst(i,mm,nn),nn=1,6)
enddo
write(37,*)
enddo

iunit = 19
do i=1,16
iunit = iunit + 1
write(iunit,'(a7,1x,a9,1x,6a8)')'# frame','BaseStep','Twist','Roll','Tilt','Rise','Shift','Slide'
enddo

!writing residue-wise force constat data

do i = 1,NBase - 1
write(10,'(a2,2x,i4,1x,6f12.6)')SeqType(typ(i)),i,Twist(i),Roll(i),Tilt(i),Rise(i),Shift(i),Slide(i)
enddo

write(10,'(a3,1x,a9,1x,6a12)')'Seq','BaseStep','F-Twist','F-Roll','F-Tilt','F-Rise','F-Shift','F-Slide'

end program ElasticProperty
!------------------------------------------------------------------------------
subroutine LU_DECOMP(a,n,SUCCESS)
!==============================================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!------------------------------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!==============================================================================
integer n
real(kind=8),dimension(n,*),intent(inout):: a
real(kind=8),dimension(n,n):: c
real(kind=8),dimension(n,n):: L,U
real(kind=8),dimension(n):: b,d,x
real(kind=8):: coeff
logical,intent(out):: SUCCESS

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

SUCCESS = .TRUE.
! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      if(a(k,k) == 0.d0)SUCCESS=.FALSE.
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do

!storing the elements of c, inverse of the original matrix
! a into a. So a will be the inverse of the matrix a itself now

do i=1,n
do j=1,n
a(i,j) = c(i,j)
enddo
enddo

end subroutine LU_DECOMP
!------------------------------------------------------------------------------
subroutine CHOLDC_INVERSE(a,n,SUCCESS)
integer,intent(in):: n
real(kind=8),dimension(n,*),intent(inout):: a
!real(kind=8),dimension(n,*),intent(out):: c
!real(kind=8),dimension(n,n):: L
real(kind=8),dimension(n):: p
real(kind=8):: summ
logical,intent(out):: SUCCESS

!Step1: Cholesky decomposition to get lower triangular matrix

SUCCESS = .TRUE.

do i=1,n
summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
if (summ <= 0.0)then
SUCCESS = .FALSE.
endif
p(i)=sqrt(summ)
a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
end do

!Step2: Calculation of inverse of lower triangular matrix obtained in step 1

do i=1,n
a(i,i) = 1.d0 / p(i)
do j = i + 1, n
summ = 0.d0
do k = i, j-1
summ = summ - a(j,k) * a(k,i)
enddo
a(j,i) = summ / p(j)
enddo
enddo


!Step3: Multiply the inverse of lower triangular matrix and 
!its transpose appropriately to get the inverse matrix of 
!the original matrix a

do i = 1, n
do j = i + 1, n
a(i,j) = 0.d0
enddo
enddo

do i = 1, n
a(i,i) = a(i,i) * a(i,i)
do k = i + 1, n
a(i,i) = a(i,i) + a(k,i) * a(k,i)
enddo

do j = i + 1, n
do k = j, n
a(i,j) = a(i,j) + a(k,i) * a(k,j)
enddo
enddo
enddo

do i = 1,  n
do j = 1, i-1
a(i,j) = a(j,i)
enddo
enddo

end subroutine CHOLDC_INVERSE
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
