module PF
    ! to which the PF template library is applied.
    use Model
    implicit none
    real(8), parameter :: epsTime = 1d-2
contains
! The following is useless (simply apply xs in elemntwise)
! And, the function that returns an array may rise
!  SEGMENTATION FALT, if its size is very huge and 
!  complied with Intel Fortran.
!    function prediction (p, xs) result (zs)
!        type (param) :: p
!        type (state) :: xs(:)
!        type (state) :: zs(size(xs))
!        integer :: k
!        integer, parameter :: is(3) = (/1,3,6/)
!
!        k = 2
!        zs = advance(p,xs)
!        zs(:)%logL = xs(:)%logL
!    end function
!
    function getObs(eps, t, o_seq, idx)
        integer :: getObs, idx
        real(8) :: eps, t, t_obs, dt, dtMin
        type (obs) :: o_seq(:)
        dtMin = 1D300
        getObs = 0 
        do while (idx .le. size(o_seq))
            t_obs = obsTimeOf(o_seq(idx))
            dt = abs (t - t_obs)
            if (t_obs < t - eps) then
                continue
            else if (t_obs <= t + eps) then
                if (dt.lt.dtMin) then
                    getObs = idx
                    dtMin  = dt
                    idx = idx + 1
                    exit
                end if
            else
                exit
            end if
            idx = idx + 1
        end do
    end function

    function seq_filter(p, y, xs) result (zs)
        type (param) :: p
        type (obs)   :: y
        type (state) :: xs(:), zs(size(xs))
        real(8), dimension(size(xs)) :: logWs, ws
        real(8) :: np, logWMin, s, c
        integer :: i, j, jj
        np      = size(xs)
        logWs   = xlogOC(xs) 
        logWMin = minval(logWs)
        ws      = xexp(logWs)
        s       = sum(ws) / np
        ws      = xdivs(ws)
        j = 1
        c = 0d0
        do i = 1, size(xs)
            if (isnan(logWs(i))) then 
                write(*,*) 'NaN arisen!!', i
                write(*,*) xs(i)%x(1:6,1)
                write(*,*) xs(i)%x(1:6,2)
                stop
            end if
        end do

        do i = 1, size(xs)
            jj = j + int(round_to_nearest_even(ws(i) + c))
            if (jj - 1 .le. size(zs))&
                zs(j:jj - 1) = xs(i)
            c = ws(i) + c - round_to_nearest_even(ws(i) + c)
            j = jj
        end do
    contains
        elemental real(8) function xlogOC(x)
            type (state), intent(in) :: x
            xlogOC = logOC(p,y,x)
        end function
        elemental real(8) function xexp(logW)
            real(8), intent(in) :: logW
            xexp = exp(logWMin - logW)
        end function
        elemental real(8) function xdivs(w)
            real(8), intent(in) :: w
            xdivs = w / s
        end function
    end function

    subroutine realMin(x,y)
        real(8) :: x, y
        y = min(x,y)
    end subroutine

    subroutine realSum(x,y)
        real(8) :: x, y
        y = x + y
    end subroutine

    function filter(p, y, xs, meanLogL1) result (zs)
        type (param) :: p
        type (obs)   :: y
        type (state) :: xs(:), zs(size(xs))
        type (state), allocatable, dimension(:) :: xs_buf
        ! note size(ws_buf), size(xs_buf) is enough if size(xs) + nproc() - 1
        real(8), dimension(size(xs)) :: logWs, ws
        !real(8), dimension(size(xs)*2) :: ws_buf
        real(8), dimension(:), allocatable :: ws_buf
        real(8) :: np, logWMin, s, c, np_g, w_avg, w_sum_me, np_me
        real(8), allocatable :: np_rank(:)
        integer :: i, j, jj, ierr, stat(MPI_STATUS_SIZE)
        integer :: src, dst, cnt_recv, cnt_send, cnt_new_p, i_buf, n_req, n_sk
        integer :: i_buf0, i0, l_msg, sizeof_state, n_recv_p
        ! to skip the problem of an unlimited number of requests
        integer, dimension(40000) :: req, ierr_req! enough if size(req) > 4*nproc()
        integer, dimension(MPI_STATUS_SIZE, 400000) :: stat_req
        logical :: flag
        type sked
            sequence
            integer :: from, to, n
        end type
        type (sked) :: sk_me(10000), sk ! enough if size(sked) > nproc()
        character(80) :: errmsg

        real(8) :: meanL1, meanLogL1

        allocate(ws_buf(size(xs) + nproc()))
        allocate(xs_buf(size(xs) + nproc()))

        np      = size(xs)
        np_g    = np * nproc()
        logWs   = xlogOC(xs) 
        xs(:)%logL = xs(:)%logL + logWs(:)
   
        call mpi_allreduce_lg(logWs, logWMin, 8, size(logWs), realMin, MPI_COMM_WORLD, ierr)
        ws      = xexp(logWs)
        w_sum_me = sum(ws)

        call mpi_allreduce_lg(ws, meanL1, 8, size(ws), realSum, MPI_COMM_WORLD, ierr)
        meanLogL1 = - (log (meanL1 / dble(np_g)) - logWMin)

        call mpi_barrier(MPI_COMM_WORLD, ierr)
        ! 1. determine the number of particles after filtering in each process
        if (me().ne.0) then
            call mpi_send(w_sum_me,1,MPI_REAL8,0,0,MPI_COMM_WORLD,ierr)
            call mpi_recv(np_me,1,MPI_REAL8,0,989,MPI_COMM_WORLD,stat,ierr)
        else
            allocate(np_rank(0:nproc()-1))
            np_rank(0) = w_sum_me
            do i = 1, nproc() - 1
                call mpi_recv(s,1,MPI_REAL8,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,stat,ierr)
                np_rank(stat(MPI_SOURCE)) = s
            end do
            w_avg = sum(np_rank) / np_g
            np_rank = np_rank / w_avg

            ! adjusting boundaries to make np_rank(:) be integers.
            c = 0d0
            do i = 0, size(np_rank) - 1
                s = np_rank(i) + c
                np_rank(i) = round_to_nearest_even(s)
                c = s - np_rank(i)
            end do

            ! notify the number of particles to each process
            do i = 1, nproc() - 1
                call mpi_send(np_rank(i),1,MPI_REAL8,i,989,MPI_COMM_WORLD,ierr)
            end do 
            np_me = np_rank(0)  ! for myself 
        end if

        ! 2. normilise weigth to be the unit of the number (and adjust to an integer)
        ws = ws / w_sum_me * np_me
        c = 0
        do i = 1, size(ws)
            s = c + ws(i)
            ws(i) = round_to_nearest_even(s)
            c = s - round_to_nearest_even(s)
        end do

        call mpi_barrier(MPI_COMM_WORLD, ierr)
        ! 3. make a schedule of data transfer (from who to whom)
        sizeof_state = loc(xs(2)) - loc(xs(1))
        n_sk = 0
        if (me().eq.0) then
            src = 0
            do dst = 0, nproc() - 1
                cnt_recv = np
                do while (cnt_recv.gt.0)
                    do while (np_rank(src).lt.1); src = src + 1; end do
                    sk%from = src
                    sk%to   = dst
                    sk%n    = min(int(np_rank(src)), cnt_recv)
                    if (sk%from.ne.0) then
                        call mpi_send(sk,3,MPI_INTEGER4,sk%from,1,MPI_COMM_WORLD,ierr)
                    else
                        n_sk = n_sk + 1    
                        sk_me(n_sk) = sk
                    end if
                    cnt_recv = cnt_recv - sk%n
                    np_rank(src) = np_rank(src) - sk%n
                end do
            end do
        else
            cnt_send = np_me
            do while (cnt_send.gt.0) 
                call mpi_recv(sk,3,MPI_INTEGER4,0,1,MPI_COMM_WORLD,stat,ierr)
                n_sk = n_sk + 1
                sk_me(n_sk) = sk
                cnt_send = cnt_send - sk%n
            end do
        end if

        call mpi_barrier(MPI_COMM_WORLD, ierr)
!        write(*,*) me(), ': n_sk = ', n_sk
!        do i = 1, n_sk
!            write(*,*) me(), ': sk(', i, ') = ', sk_me(i)
!        end do
        call mpi_barrier(MPI_COMM_WORLD, ierr)

        ! 3. issue send requests according to given schedules
        i     = 1
        i_buf = 1
        n_req = 0 
        do j = 1, n_sk
            i_buf0 = i_buf
            i0     = i 
            l_msg = 0

            if (sk_me(j)%n .lt. 1) then
                write(*,*) 'pf.90: process ', me(), 'has no probabilistic weight! terminated abnormally'
                stop
            end if

            do while (sk_me(j)%n .gt. 0)
                do while (int(ws(i)).eq.0)
                    if (i.eq.size(ws)) then
                        write(*,*) me(), ': lack of particle: ', sk_me(j)%n, np_me
                        stop
                    end if
                    i = i + 1
                end do
                ws_buf(i_buf) = min(ws(i), dble(sk_me(j)%n))
                xs_buf(i_buf) = xs(i)
                ws(i)         = ws(i)      - ws_buf(i_buf)
                sk_me(j)%n    = sk_me(j)%n - ws_buf(i_buf)
                l_msg = l_msg + 1
                i_buf = i_buf + 1
            end do
            n_req = n_req + 1
            call mpi_isend(ws_buf(i_buf0), l_msg, MPI_REAL8,&
                sk_me(j)%to, 3, MPI_COMM_WORLD, req(n_req), ierr) 
            n_req = n_req + 1
            call mpi_isend(xs_buf(i_buf0), l_msg*sizeof_state, MPI_BYTE,&
                sk_me(j)%to, 4, MPI_COMM_WORLD, req(n_req), ierr)

            !write(*,*) sk_me(j)%to, ',',&
            !    int(sum(ws_buf(i_buf0:i_buf0+l_msg-1))), l_msg, ', from ', me(), ', send'

            if (ierr .ne. 0) then
                call mpi_error_string(ierr, errmsg, ierr)
                write(*,*) me(), ': send error ', errmsg
            end if
        end do

        ! 4. complete requests, with issuing recv requests on demand 
        i = 1
        cnt_recv = 0
        cnt_new_p = 0
        do
            call mpi_iprobe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, flag, stat, ierr)
            if (flag) then
                n_req = n_req + 1
                if (i .gt. size(ws)) then
                    write(*,*) me(), 'subscript out of range arisen!'
                    stop
                end if
                call mpi_irecv(ws(i), size(ws) - i + 1, MPI_REAL8,&
                    stat(MPI_SOURCE), 3, MPI_COMM_WORLD, req(n_req), ierr)
                call mpi_get_count(stat, MPI_REAL8, l_msg, ierr)
                n_req = n_req + 1
                call mpi_irecv(zs(i), (size(zs) - i + 1)*sizeof_state, MPI_BYTE,&
                    stat(MPI_SOURCE), 4, MPI_COMM_WORLD, req(n_req), ierr) 
                cnt_recv = cnt_recv + l_msg
                i = i + l_msg
            end if
            call mpi_testall(n_req, req, flag, stat_req, ierr_req)
            if (flag) then
                cnt_new_p = sum(ws(1:cnt_recv))
                if (cnt_new_p.ge.size(xs)) exit
            end if
        end do
        do i = 1, n_req
            if (req(i).ne.MPI_REQUEST_NULL) then
                write(*,*) me(), ': req isnot MPI_REQUEST_NULL'
            end if
        end do

        ! 5. expand new particles
        j = cnt_recv + 1
        
        ! 2015-12-02: a quick hack to fix a bug
        ! In a very small prorability, max(jj) get greater than size(zs). This is a memory violation.
        ! np_rank and ws may not be consistent strictrically: the code seems to do differnt rounding.
        ! For a mean while, I add a guard to avoid the memory viiolation
        do i = 1, cnt_recv
            jj = j + int(ws(i)) - 2
            j = min(j,size(ws))    ! -- guard --
            jj = min(jj,size(ws))  ! -- guard --
            zs(j:jj) = zs(i)
            j = jj + 1
        end do
        
        if (me().eq.0) deallocate(np_rank)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
    contains
        elemental real(8) function xlogOC(x)
            type (state), intent(in) :: x
            xlogOC = logOC(p,y,x)
        end function
        elemental real(8) function xexp(logW)
            real(8), intent(in) :: logW
            xexp = exp(logWMin - logW)
        end function
        elemental real(8) function xdivs(w)
            real(8), intent(in) :: w
            xdivs = w / s
        end function
    end function

    subroutine filterApp(p, cond, recorder, o_seq, xs, negLogL)
        type (param) :: p
        type (obs)   :: o_seq(:)
        type (state) :: xs(:)!, zs(size(xs))
        integer      :: idx, i, ipar
        real(8)      :: meanLogL1(size(o_seq)), tmp, negLogL
        interface
            logical function cond(o_seq,idx,xs)
                import :: state, obs, param
                type (obs) :: o_seq(:)
                integer :: idx
                type (state) :: xs(:)
            end function
            subroutine recorder(xs,y)
               import :: state, obs, param
               type (state) :: xs(:)
               type (obs), optional :: y
            end subroutine
        end interface
        idx = 1
        xs(:)%logL = 0D0
        meanLogL1 = 0D0
        !write(*,*) size(meanLogL1)
        do while (cond(o_seq, idx, xs))
            xs = advance(p,xs)
            i = getObs(epsTime, timeOf(xs(1)), o_seq, idx)
            if (i > 0) then
                xs = filter(p, o_seq(i), xs, meanLogL1(i))!tmp)
                !meanLogL1(i) = tmp
                call recorder(xs, o_seq(i))
                do ipar = 1,size(xs)
                  xs(ipar) = postObs(p, o_seq(i), xs(ipar),ipar)
                end do
            else
                call recorder(xs)
            end if
            !    write(*,*) i, xs(1)%x(:,1)
        end do
        negLogL = sum(meanLogL1)
    end subroutine

    ! meanLogL returns negative log-likelihood
    real(8) function meanLogL (xs)
        type (state) :: xs(:)
        real(8) :: logL (size(xs)), p(size(xs))
        real(8) :: l_min
        integer :: ierr
        logL = xs%logL
        call mpi_allreduce_lg(logL, l_min, 8, size(logL), realMin, MPI_COMM_WORLD, ierr)
        p = xexp(logL)
        call mpi_reduce_lg(p, meanLogL, 8, size(p), realSum, 0, MPI_COMM_WORLD, ierr)
        meanLogL = - (log (meanLogL / dble(size(xs)*nproc())) - l_min)
    contains
        elemental real(8) function xexp (l)
            real(8), intent(in) :: l
            xexp = exp (l_min - l)
        end function      
    end function

    function mapState(xs) result (x)
        type (state) :: xs(:), x, xlocal
        integer :: n, szState, ierr
        n = size(xs)
        szState = loc(xs(2)) - loc(xs(1))
        ! short hand of mpi_reduce. see misc.f90.
        call mpi_reduce_lg(xs, x, szState, size(xs), xopLHMin, 0, MPI_COMM_WORLD, ierr)
    end function

    subroutine xopLHMin(invec, inoutvec)
        type (state) :: invec, inoutvec
        if (invec%logL .lt. inoutvec%logL) inoutvec = invec
    end subroutine    
end module 
