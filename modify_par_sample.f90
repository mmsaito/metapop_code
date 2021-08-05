function modify_par_version() result (s)
character(len=19):: s
s = '2016-10-05T19:18:43'
end function


subroutine init_par(p)
type(param) :: p
p%betaJmp(1,1)%value = 1.0d0
p%betaJmp(2,1)%value = 1.0d0
p%betaJmp(3,1)%value = 1.0d0
p%betaJmp(4,1)%value = 1.0d0
p%betaJmp(5,1)%value = 1.0d0
p%betaJmp(6,1)%value = 1.0d0
p%betaJmp(7,1)%value = 1.0d0
p%betaJmp(8,1)%value = 1.0d0
p%betaJmp(1,2)%value = 1.0d0
p%betaJmp(2,2)%value = 1.0d0
p%betaJmp(3,2)%value = 1.0d0
p%betaJmp(4,2)%value = 1.0d0
p%betaJmp(5,2)%value = 1.0d0
p%betaJmp(6,2)%value = 1.0d0
p%betaJmp(7,2)%value = 1.0d0
p%betaJmp(8,2)%value = 1.0d0
end subroutine


subroutine modify_par(p,i)
type(param) :: p
integer :: i
999 continue
select case (i)
case (0)
p%betaJmp(1,1)%value = p%betaJmp(1,1)%value + rgauss()*5d-2
if (p%betaJmp(1,1)%value < 0d0 .or. p%betaJmp(1,1)%value > 800d0) goto 999
case (1)
p%betaJmp(2,1)%value = p%betaJmp(2,1)%value + rgauss()*5d-2
if (p%betaJmp(2,1)%value < 0d0 .or. p%betaJmp(2,1)%value > 800d0) goto 999
case (2)
p%betaJmp(3,1)%value = p%betaJmp(3,1)%value + rgauss()*5d-2
if (p%betaJmp(3,1)%value < 0d0 .or. p%betaJmp(3,1)%value > 800d0) goto 999
case (3)
p%betaJmp(4,1)%value = p%betaJmp(4,1)%value + rgauss()*5d-2
if (p%betaJmp(4,1)%value < 0d0 .or. p%betaJmp(4,1)%value > 800d0) goto 999
case (4)
p%betaJmp(5,1)%value = p%betaJmp(5,1)%value + rgauss()*5d-2
if (p%betaJmp(5,1)%value < 0d0 .or. p%betaJmp(5,1)%value > 800d0) goto 999
case (5)
p%betaJmp(6,1)%value = p%betaJmp(6,1)%value + rgauss()*5d-2
if (p%betaJmp(6,1)%value < 0d0 .or. p%betaJmp(6,1)%value > 800d0) goto 999
case (6)
p%betaJmp(7,1)%value = p%betaJmp(7,1)%value + rgauss()*5d-2
if (p%betaJmp(7,1)%value < 0d0 .or. p%betaJmp(7,1)%value > 800d0) goto 999
case (7)
p%betaJmp(8,1)%value = p%betaJmp(8,1)%value + rgauss()*5d-2
if (p%betaJmp(8,1)%value < 0d0 .or. p%betaJmp(8,1)%value > 800d0) goto 999
case (8)
p%betaJmp(1,2)%value = p%betaJmp(1,2)%value + rgauss()*5d-2
if (p%betaJmp(1,2)%value < 0d0 .or. p%betaJmp(1,2)%value > 800d0) goto 999
case (9)
p%betaJmp(2,2)%value = p%betaJmp(2,2)%value + rgauss()*5d-2
if (p%betaJmp(2,2)%value < 0d0 .or. p%betaJmp(2,2)%value > 800d0) goto 999
case (10)
p%betaJmp(3,2)%value = p%betaJmp(3,2)%value + rgauss()*5d-2
if (p%betaJmp(3,2)%value < 0d0 .or. p%betaJmp(3,2)%value > 800d0) goto 999
case (11)
p%betaJmp(4,2)%value = p%betaJmp(4,2)%value + rgauss()*5d-2
if (p%betaJmp(4,2)%value < 0d0 .or. p%betaJmp(4,2)%value > 800d0) goto 999
case (12)
p%betaJmp(5,2)%value = p%betaJmp(5,2)%value + rgauss()*5d-2
if (p%betaJmp(5,2)%value < 0d0 .or. p%betaJmp(5,2)%value > 800d0) goto 999
case (13)
p%betaJmp(6,2)%value = p%betaJmp(6,2)%value + rgauss()*5d-2
if (p%betaJmp(6,2)%value < 0d0 .or. p%betaJmp(6,2)%value > 800d0) goto 999
case (14)
p%betaJmp(7,2)%value = p%betaJmp(7,2)%value + rgauss()*5d-2
if (p%betaJmp(7,2)%value < 0d0 .or. p%betaJmp(7,2)%value > 800d0) goto 999
case (15)
p%betaJmp(8,2)%value = p%betaJmp(8,2)%value + rgauss()*5d-2
if (p%betaJmp(8,2)%value < 0d0 .or. p%betaJmp(8,2)%value > 800d0) goto 999
end select
i = modulo(i + 1,16)
end subroutine


