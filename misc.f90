! a vector among processes are reduced by op into a scalar.
!
! at rank 0: s(0) = op(op(op(v(1),v(2)),v(3), ...))
! at rank 1: s(1) = op(op(op(v(1),v(2)),v(3), ...))
! ...
! s = op(op(op(s(1),s(2)), ...))
subroutine mpi_reduce_lg(vec, ans, sizeof, n, op, root, comm, ierr)
!    use MPI
    implicit none
    include 'mpif.h'
    integer :: sizeof, n, root, comm, ierr
    integer(1) :: vec(sizeof, n)
    integer(1), dimension(sizeof) :: ans, tmp
    integer :: i, hdl_typ, hdl_op
    interface
        subroutine op(x, y)
            integer(1) :: x(*), y(*)
        end subroutine
    end interface
    procedure(), pointer :: ptr_op
    call mpi_type_contiguous(sizeof, MPI_BYTE, hdl_typ, ierr)
    call mpi_type_commit(hdl_typ, ierr)
    ptr_op => inter_op
    call mpi_op_create(ptr_op, .true., hdl_op, ierr)
    tmp = vec(:,1)
    do i = 2, n
        call op(vec(:,i), tmp)
    end do
    call mpi_reduce(tmp, ans, 1, hdl_typ, hdl_op, root, comm, ierr)
    call mpi_type_free(hdl_typ, ierr)
    call mpi_op_free(hdl_op, ierr)
contains
    subroutine inter_op(invec, inoutvec, length, typ)
        integer(1), dimension(sizeof, length) :: invec, inoutvec
        integer :: length, typ
        integer :: i
        do i = 1, length
            call op(invec(:,i), inoutvec(:,i))            
        end do
    end subroutine
end subroutine

! allreduce version
subroutine mpi_allreduce_lg(vec, ans, sizeof, n, op, comm, ierr)
!    use MPI
    implicit none
    include 'mpif.h'
    integer :: sizeof, n, comm, ierr
    integer(1) :: vec(sizeof, n)
    integer(1), dimension(sizeof) :: ans, tmp
    integer :: i, hdl_typ, hdl_op
    interface
        subroutine op(x, y)
            integer(1) :: x(*), y(*)
        end subroutine
    end interface
    procedure(), pointer :: ptr_op
    call mpi_type_contiguous(sizeof, MPI_BYTE, hdl_typ, ierr)
    call mpi_type_commit(hdl_typ, ierr)
    ptr_op => inter_op
    call mpi_op_create(ptr_op, .true., hdl_op, ierr)
    tmp = vec(:,1)
    do i = 2, n
        call op(vec(:,i), tmp)
    end do
    call mpi_allreduce(tmp, ans, 1, hdl_typ, hdl_op, comm, ierr)
    call mpi_type_free(hdl_typ, ierr)
    call mpi_op_free(hdl_op, ierr)
contains
    subroutine inter_op(invec, inoutvec, length, typ)
        integer(1), dimension(sizeof, length) :: invec, inoutvec
        integer :: length, typ
        integer :: i
        do i = 1, length
            call op(invec(:,i), inoutvec(:,i))            
        end do
    end subroutine
end subroutine


module PF_Misc
!    use mpi
! implicit none !ToDo: activate this
    include 'mpif.h'
    real(8), parameter :: pi = 3.141592653589793D0
    type meansd_t
        real(8) :: mean, sd
    end type
    interface
        elemental real(8) function rand(dummy)
            real(8), intent(in), optional :: dummy
            !integer, intent(in) :: i
        end function
        
        pure subroutine mywrite(s,t)
          character(len=*), intent(in) :: s
           character(len=*), intent(out) :: t
        end subroutine

        pure subroutine mystop()
        end subroutine
    end interface

    type result_idx32_qsort
      integer, allocatable :: loc(:,:)
      real(8), allocatable :: values_all(:)
    end type
contains

  integer function realComp(x,y)
    real(8), intent(in) :: x, y
    !write(*,*) x, y
    if (x.gt.y) then
      realComp = 1
    else if (x.lt.y) then
      realComp = -1
    else
      realComp = 0
    end if
  end function

    ! fortran's (a)nint doesn't round to the nearest even.
    ! Note that this code eventually owes to rounding-mode being "to the nearest even", 
    !   which is the default rounding mode.
    elemental function round_to_nearest_even(x) result (z)
      real(8), intent(in) :: x
      real(8), volatile :: z
      real(8), parameter :: maxInt = 4503599627370496d0
      if (x.ge.0) then 
        z = x + maxInt
        z = z - maxInt
      else
        z = x - maxInt
        z = z + maxInt
      end if
    end function

    pure real(8) function rgauss (dummy)
        real(8) :: u1, u2
        real(8), optional, intent(in) :: dummy
        u1 = rand(u1)
        u2 = rand(u2)
        rgauss = sqrt(-2D0*log(u1))*cos(2D0*pi*u2)
    end function

    pure real(8) function lnGamma(z)
        real(8), intent(in) :: z
        lnGamma =&
          5D-1*( log (pi+pi)&
               - log(z)&
               + z*( 2.0*log(z)&
                   + log (z*sinh(1.0/z) + 1.0/810.0/z**6)&
                   - 2.0 &
                   )&
               )
    end function

    pure real(8) function rbinom(n,p)
      real(8), intent(in) :: n, p
      real(8) :: spq, b, a, c, vr, u, us, v
      real(8) :: alpha, m, h, lpq
      ! 2015-11-29: a quick hack to avoid a stuck at p = 0 or 1.
      rbinom = 0d0
      if (p.eq.0d0) then
        return
      else if (p.eq.1d0) then
        rbinom = n
        return
      end if

      !integer, optional :: cnt
      spq = sqrt(n*p*(1D0-p))
      b = 1.15D0 + 2.53D0*spq
      a = -0.0873D0 + 0.0248D0*b + 0.01D0*p
      c = n*p + 0.5D0
      vr = 0.92D0 - 4.2D0/b
!      cnt = 1
   10 continue
        u = rand() - 0.5D0
        us = 0.5 - abs(u)
        rbinom = floor((2D0*a/us + b)*u + c)
        v = rand()
        if (us.ge.0.07D0.and.v.le.vr) return
        if (rbinom.lt.0.or.rbinom.gt.n) goto 10
        alpha = (2.83D0 + 5.1D0/b)*spq
        lpq = log(p/(1D0 - p))
        m = floor((n+1D0)*p)
        h = lnGamma(m + 1D0) + lnGamma(n - m + 1D0)
        v = v*alpha/(a/us/us + b)
        if (v.le.h - lnGamma(rbinom + 1.0D0)&
                   - lnGamma(n - rbinom + 1.0D0)&
                   + (rbinom - m)*lpq) return
!        cnt = cnt + 1
        goto 10
    end function

    pure real(8) function rpoisson (lambda) result (n)
      real(8), intent(in) :: lambda
      real(8) :: s, r
      character(255) :: ss
      s = exp(lambda)
      n = 0D0

      do
        r = rand(r) ! r is passed to suppress the optimization
        s = s * r
!        if (n > 50*lambda) then 
!          write(ss,*)  lambda, s, n
!          call mywrite(ss,ss)
!        end if
        if (s < 1.0) return
        n = n + 1
      end do
    end function

    pure real(8) function lnPo (m, x) 
        real(8), intent(in) :: m, x
        lnPo = -m + x * log(m) - lnGamma (x + 1D0)
    end function

    pure real(8) function lnBinom(n,p,x)
      real(8), intent(in) :: n, p, x
      integer :: k, n_, x_
      real(8) :: hoge(5)

      n_ = anint(n)
      x_ = anint(x)

      lnBinom = sum((/ (log(dble(k)), k = n_ - x_ + 1, n_) /) ) &
                  - sum((/(log(dble(k)), k = 1, x_)/)) &
                    + x*log(p) + (n - x)*log(1d0 - p)
    end function

    type (meansd_t) function meansd(xs)
        real(8), dimension(:) :: xs
        real(8) :: n
        real(8) :: s, mu, var
        integer :: i
        n = size(xs)
        s = 0d0
        do i = 1, size(xs)
            s = s + xs(i)
        end do
        mu = s / n
        s = 0d0
        do i = 1, size(xs)
            s = s + xs(i)**2
        end do
        var = s / n
        meansd = meansd_t (mu, dsqrt(var))
    end function

    elemental logical function isnan(x)
        real(8), intent(in) :: x
        isnan = .not. (x .eq. x)
    end function

    ! inRange: check if i in [b,e]
    pure logical function inRange(i,b,e)
      integer, intent(in) :: i,b,e
      inRangeInt = b.le.i .and. i.le.e
    end function

    pure logical function inIdxRange(i,arr)
      type(*), intent(in) :: arr(:)
      integer, intent(in) :: i
      inIdxRange = inRange(i,lbound(arr,1),ubound(arr,1))
    end function

    pure real(8) function rbinom_approx(n,p) result (v)
      real(8), intent(in) :: n, p
      character(len=32) :: ss

      if (isnan(n)) then; 
        write(ss,*) 'rbinom: n is nan' 
        call mywrite(ss,ss)
        v = n
        return
      else if (isnan(p)) then; 
        write(ss,*) 'rbinom: p is nan' 
        call mywrite(ss,ss)
        v = p
        return
      end if

      if (n*p .gt. 50) then
        v = min(max(aint(n*p + sqrt(n*p)*rgauss(v)),0d0),n)
      else
        v = min(rpoisson (lambda = n*p),n)  !quick hack
      end if
    end function

    !sR, sI, printing numbers
    pure character(len=24) function sR(x)
      real(8), intent(in) :: x
      write(sR,*) x
      sR = adjustl(sR)
    end function

    pure character(len=12) function sI(x)
      integer, intent(in) :: x
      write(sI,*) x
      sI = adjustl(sI)
    end function

    ! short cut of mpi functions
    integer function me()
        integer :: ierr
        call mpi_comm_rank(mpi_comm_world, me, ierr)
    end function

    integer function nproc()
        integer :: ierr
        call mpi_comm_size(mpi_comm_world, nproc, ierr)
    end function

    type (result_idx32_qsort) function mpi_seq_idx32_qsort(x, iroot, comp) result (z)
      real(8), intent(in) :: x(:)
      integer, intent(in) :: iroot
      interface
        integer function comp(x,y)
          real(8), intent(in) :: x,y
        end function
      end interface

      integer :: ierr, n_all
      integer, allocatable :: idxbuf(:)

      n_all = size(x)*nproc()
      allocate(idxbuf(n_all))

      allocate(z%loc(2,n_all))
      allocate(z%values_all(n_all))

      call mpi_gather(x, size(x), MPI_REAL8, z%values_all, size(x), MPI_REAL8, iroot, MPI_COMM_WORLD, ierr)

      idxG = 1
      do iproc = 0, nproc() - 1
        do idx = 1, size(x)
          idxbuf(idxG) = idxG
          z%loc(1,idxG) = iproc
          z%loc(2,idxG) = idx
          idxG = idxG + 1
        end do
      end do

      call idx32_qsort(z%values_all, n_all, 8, idxbuf, realComp, 1)

      z%loc(1,:) = z%loc(1,idxbuf)
      z%loc(2,:) = z%loc(2,idxbuf)
      z%values_all(:) = z%values_all(idxbuf)
    end function

    subroutine mpi_real_quantile(xqs, xs, qs)
      real(8), intent(in) :: xs(:), qs(:)
      real(8), intent(out) :: xqs(:)

      type (result_idx32_qsort) :: ord
      integer :: irank(size(qs)), n_all

      ord = mpi_seq_idx32_qsort(xs, 0, realComp)
      n_all = size(xs)*nproc()
      irank = 1 + qs*(n_all - 1) ! irank = (1 - q)*1 + q*N
      xqs(1:size(irank)) = ord%values_all(irank)
    end subroutine

    real(8) function mpi_real_mean(xs) result (z)
      real(8), intent(in) :: xs(:)
      real(8) :: tmp
      tmp = sum(xs)
      call mpi_reduce(tmp, z, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (me().eq.0) then
        z = z/dble(size(xs)*nproc())
      end if
    end function  

    real(8) function mpi_real_rmse(y,xs) result (z)
      real(8), intent(in) :: y, xs(:)
      real(8) :: tmp, ans
      integer :: ierr
      tmp = sum((xs(:) - y)**2)
      call mpi_reduce(tmp, z, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (me().eq.0) then
        z = sqrt(z/dble(size(xs)*nproc()))
      end if
    end function  
end module
