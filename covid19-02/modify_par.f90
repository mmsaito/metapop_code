function modify_par_version() result (s)
character(len=19):: s
!s = '2018-06-28T17:41:22'  ! edit it to include a part of a path name
s = '2020-04-21-nseg#2'  ! Windows not allows to use colon as a path name.
end function


subroutine init_par(p)
type(param) :: p

p%nGroup   = 1 ! ìññ ÅCì˙ñ{ëSçëìØÇ∂Ç∆Ç∑ÇÈ
p%nBetaJmp = 2

p%betaJmp(1,1)%tAt   = -1d0
p%betaJmp(1,1)%value = 1d0

p%betaJmp(2,1)%tAt   = 30d0   !the sentinel time point
p%betaJmp(2,1)%value = 2d0

p%betaJmp(3,1)%tAt   = 1d300  !the sentinel time point
p%betaJmp(3,1)%value = 1d0
end subroutine

subroutine modify_par(p,i)
type(param) :: p
integer :: i
999 continue
select case (i)
case (0)
p%betaJmp(1,1)%value = p%betaJmp(1,1)%value + rgauss()*5d-2
if (p%betaJmp(1,1)%value < 0d0 .or. p%betaJmp(1,1)%value > 5d0) goto 999

case (1)
p%betaJmp(2,1)%value = p%betaJmp(2,1)%value + rgauss()*5d-2
if (p%betaJmp(2,1)%value < 0d0 .or. p%betaJmp(2,1)%value > 5d0) goto 999


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The below is dead code!
case (2)
p%betaJmp(3,1)%value = p%betaJmp(3,1)%value + rgauss()*5d-2
if (p%betaJmp(3,1)%value < 0d0 .or. p%betaJmp(3,1)%value > 5d0) goto 999
case (3)
p%betaJmp(4,1)%value = p%betaJmp(4,1)%value + rgauss()*5d-2
if (p%betaJmp(4,1)%value < 0d0 .or. p%betaJmp(4,1)%value > 5d0) goto 999
case (4)
p%betaJmp(5,1)%value = p%betaJmp(5,1)%value + rgauss()*5d-2
if (p%betaJmp(5,1)%value < 0d0 .or. p%betaJmp(5,1)%value > 5d0) goto 999
case (5)
p%betaJmp(6,1)%value = p%betaJmp(6,1)%value + rgauss()*5d-2
if (p%betaJmp(6,1)%value < 0d0 .or. p%betaJmp(6,1)%value > 5d0) goto 999
case (6)
p%betaJmp(1,2)%value = p%betaJmp(1,2)%value + rgauss()*5d-2
if (p%betaJmp(1,2)%value < 0d0 .or. p%betaJmp(1,2)%value > 5d0) goto 999
case (7)
p%betaJmp(2,2)%value = p%betaJmp(2,2)%value + rgauss()*5d-2
if (p%betaJmp(2,2)%value < 0d0 .or. p%betaJmp(2,2)%value > 5d0) goto 999
case (8)
p%betaJmp(3,2)%value = p%betaJmp(3,2)%value + rgauss()*5d-2
if (p%betaJmp(3,2)%value < 0d0 .or. p%betaJmp(3,2)%value > 5d0) goto 999
case (9)
p%betaJmp(4,2)%value = p%betaJmp(4,2)%value + rgauss()*5d-2
if (p%betaJmp(4,2)%value < 0d0 .or. p%betaJmp(4,2)%value > 5d0) goto 999
case (10)
p%betaJmp(5,2)%value = p%betaJmp(5,2)%value + rgauss()*5d-2
if (p%betaJmp(5,2)%value < 0d0 .or. p%betaJmp(5,2)%value > 5d0) goto 999
end select
i = modulo(i + 1,2)
end subroutine


