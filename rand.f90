real(8) function rand(dummy)
    real(8), optional :: dummy
    call random_number(rand)
end function

subroutine mywrite(s,t)
  character(len=*) :: s,t
  write(*,*) s
  call flush(6)
end subroutine

subroutine mystop()
  stop
end subroutine
