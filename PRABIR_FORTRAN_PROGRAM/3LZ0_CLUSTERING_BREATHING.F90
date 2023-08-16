program APS
character(len=200):: Loop1,Loop2,OutFile
integer,dimension(0:3):: nCount
real:: nc,n1,n2

write(*,'(1x,a)')'Enter the two input files for Tail1 and Tail2'
read(*,'(a)')Loop1
write(*,*)Loop1
read(*,'(a)')Loop2
write(*,*)Loop2

write(*,'(1x,a)')'Enter the name of output file'
read(*,'(a)')OutFile
write(*,*)OutFile

open(10,file=Loop1,status='old')
open(12,file=Loop2,status='old')
open(14,file=OutFile,status='unknown')

write(*,'(1x,a)')'Enter the cut-off brathing distance'
read(*,*)nc
write(*,*)'nc===>',nc

write(*,'(1x,a)')'Enter the no of data points'
read(*,*)ndata
write(*,*)'ndata====>',ndata

nCount(:) = 0

do i=1,ndata
read(10,*)x,n1
read(12,*)x,n2

if(n1 >= nc.and.n2 >= nc)then
        ID = 3
else
        if(n1 >= nc.and.n2 < nc)then
                ID = 1
        else
                if(n1 < nc.and.n2 >= nc)then
                        ID = 2
                else
                        ID = 0
                endif
        endif
endif

nCount(ID) = nCount(ID) + 1
write(14,'(i8,1x,i2)')i,ID
enddo

do i = 0,3
write(14,'(a13,1x,i2,1x,i8)')'ID,Count====>',i,nCount(i)
enddo

end program APS
