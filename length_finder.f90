subroutine length(natm)
implicit none
character :: junk1*20
integer :: i,k,natm

integer:: j
Do i=1,4
    open(unit=1, file= "C:\Users\User\Desktop\Fortran Codings\generating PDB file\Trajectory_input.txt")
 READ(1,'(a)') junk1
 k=INDEX(TRIM(junk1),'NUMBER',.FALSE.)
!write(*,*) i
IF (k .GT. 0) then
cycle
end if
 end do
 !write(*,*)  junk1
read(junk1,*) natm
!write(*,*) natm
100 format (a)
end subroutine length
