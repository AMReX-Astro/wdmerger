program logreader

implicit none

character(len=30):: filename,varname

integer :: status,i,nvals = 0
real(kind = 8) :: temp

call GETARG(1,filename)

open (unit = 9, file = filename, status = "OLD", action = "READ", iostat = status)

!Skip the first line, it has the variable names in it
read (9,20)
20 format(/)
backspace(9)

!Find out how many lines of data there are
do
   read (9,*,iostat = status)

   if (status < 0) exit
   
   nvals = nvals + 1

enddo
!Rewind and skip the first line again
rewind(9)
read (9,20)
backspace(9)

write(*,*) "Enter variable name. Options are: time, mass, pmass, smass, xmom,"
write(*,*) "ymom, zmom, kineticenergy, potentialenergy, internalenergy, totalenergy,"
write(*,*) "xangmom, yangmom, zangmom, xcom, ycom, zcom, pxcom, pycom, pzcom, sxcom,"
write(*,*) "sycom, szcom, xvel, yvel, zvel, pxvel, pyvel, pzvel, sxvel, syvel, szvel."
read(*,*) varname

!Choose which column to read based on user-inputted variable
if (varname == "time") then
      do i = 1, nvals
         read (9,30) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "mass") then
      do i = 1, nvals
         read (9,40) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "xmom") then
      do i = 1, nvals
         read (9,50) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "ymom") then
      do i = 1, nvals
         read (9,60) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "zmom") then
      do i = 1, nvals
         read (9,70) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "kineticenergy") then
      do i = 1, nvals
         read (9,80) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "potentialenergy") then
      do i = 1, nvals
         read (9,90) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "internalenergy") then
      do i = 1, nvals
         read (9,100) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "totalenergy") then
      do i = 1, nvals
         read (9,110) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "xangmom") then
      do i = 1, nvals
         read (9,120) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "yangmomy") then
      do i = 1, nvals
         read (9,130) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "zangmom") then
      do i = 1, nvals
         read (9,140) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "xcom") then
      do i = 1, nvals
         read (9,150) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "ycom") then
      do i = 1, nvals
         read (9,160) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "zcom") then
      do i = 1, nvals
         read (9,170) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pmass") then
      do i = 1, nvals
         read (9,180) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "smass") then
      do i = 1, nvals
         read (9,190) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pxcom") then
      do i = 1, nvals
         read (9,200) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "sxcom") then
      do i = 1, nvals
         read (9,210) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pycom") then
      do i = 1, nvals
         read (9,220) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "sycom") then
      do i = 1, nvals
         read (9,230) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pzcom") then
      do i = 1, nvals
         read (9,240) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "szcom") then
      do i = 1, nvals
         read (9,250) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "xvel") then
      do i = 1, nvals
         read (9,260) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "yvel") then
      do i = 1, nvals
         read (9,270) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "zvel") then
      do i = 1, nvals
         read (9,280) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pxvel") then
      do i = 1, nvals
         read (9,290) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "sxvel") then
      do i = 1, nvals
         read (9,300) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pyvel") then
      do i = 1, nvals
         read (9,310) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "syvel") then
      do i = 1, nvals
         read (9,320) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "pzvel") then
      do i = 1, nvals
         read (9,330) temp
         write(*,'(ES22.15)') temp
      enddo

   elseif (varname == "szvel") then
      do i = 1, nvals
         read (9,340) temp
         write(*,'(ES22.15)') temp
      enddo

   else
      write(*,*) "Invalid variable name"
      stop

endif

!Where to look for each variable
30 format(T7,E22.0)
40 format(T25,E22.0)
50 format(T48,E22.0)
60 format(T71,E22.0)
70 format(T94,E22.0)
80 format(T117,E22.0)
90 format(T140,E22.0)
100 format(T163,E22.0)
110 format(T186,E22.0)
120 format(T209,E22.0)
130 format(T232,E22.0)
140 format(T255,E22.0)
150 format(T278,E22.0)
160 format(T301,E22.0)
170 format(T324,E22.0)
180 format(T347,E22.0)
190 format(T370,E22.0)
200 format(T393,E22.0)
210 format(T416,E22.0)
220 format(T439,E22.0)
230 format(T462,E22.0)
240 format(T485,E22.0)
250 format(T508,E22.0)
260 format(T531,E22.0)
270 format(T554,E22.0)
280 format(T577,E22.0)
290 format(T600,E22.0)
300 format(T623,E22.0)
310 format(T646,E22.0)
320 format(T669,E22.0)
330 format(T692,E22.0)
340 format(T715,E22.0)

end program logreader
