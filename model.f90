!   - summariser: mean, sd, ...
!   - parameter maker
module Model
  use PF_Misc
  implicit none

  integer, parameter :: idxS = 1
  integer, parameter :: idxI = 2
  integer, parameter :: idxR = 3
  integer, parameter :: idxS2I = 4

  integer, parameter :: nSite = 47

  !integer, parameter :: nWeek = 52
  integer, parameter :: nWeek = 82   
    ! very dirty quick-hack. asap, switch to an allocatable array
    ! plus, time step is no longer week but day.

  !integer, parameter :: nBetaJmp = 2

  logical :: calc_CI = .true.

  type jump_point
    real(8) :: tAt, value 
  end type

  type state
      real(8) :: x(nSite,4), x0(nSite,4)
      real(8) :: logL  ! state & state' are identified in Fortran ver.
      real(8) :: t

      ! Augumented stateを使って betaを推定するときのため
      ! param::beta を単位とする: beta = rbeta * param::beta
      real(8) :: rbeta(nSite,5), sd_rbeta
  end type

  type table_pop
      real(8), dimension(nSite) :: total, man, woman
      character(len=16) :: name(nSite)
  end type

  type table_od_time
    real(8), dimension(nSite,nSite) :: od
    integer    :: t
    integer(1) :: mm,dd
    integer(2) :: yy
  end type

  type param
      real(8) :: rho(nSite,nSite), beta, gamma
      real(8) :: prob_mig_pre(nSite,nSite)
      real(8) :: prob_mig(nSite,nSite)
      real(8) :: dt
      real(8) :: muS
      integer :: nBetaJmp, nGroup
      type(jump_point) :: betaJmp(20,2)
      integer :: groupid(nSite)
      logical :: useAggregatedModel = .false.
      real(8) :: conn_redu = 1d0, contact_pow = 0d0
      logical :: useBinomObsModel = .false.
      real(8) :: probBinomObs = 0.5d0

      ! OD行列が時変なシミュレーション用
      type(table_od_time), pointer :: od_seq(:)

      ! 人口テーブル
      type(table_pop), pointer :: tbl_pop
      ! OD行列が時間依存のときのテーブルへのインデックス
      integer :: i_od_seq
  end type

  type obs
      real(8) :: t, j(nSite), rbeta(nSite)
  end type

contains
    pure real(8) function pr_rate(rate_t) 
      real(8), intent(in) :: rate_t
      character(len=16) :: ss
      pr_rate = 1d0 - exp(-rate_t)
      if (isnan(pr_rate)) call mywrite('*',ss)
    end function

    elemental function advance(p,x) result (z)
        type (state), intent(in) :: x
        type (state) :: z
        type (param), intent(in) :: p

        real(8), dimension(nSite,nSite) :: d_ss, d_ii, d_rr
        real(8) :: d_si, d_ir

        real(8) :: n_present(nSite), tmp, beta(nSite), cr

        integer :: i, j, igr, lag
        character(len=256) ss
        logical :: mask(nSite)

        ! if useAggregatedModel is active, another competitor model is used instead.
        ! this model is used only for a performance test.
        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if (p%useAggregatedModel) then
          z = advance_aggregated(p,x)
          return
        end if
        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        forall (i=1:nSite) n_present(i) = sum(x%x(i,1:3))
        do i = 1, nSite
          do j = 1, nSite
            if (i.ne.j) then
              d_ss(i,j) = rbinom_approx(n = x%x(i,idxS), p = p%prob_mig(i,j))
              d_ii(i,j) = rbinom_approx(n = x%x(i,idxI), p = p%prob_mig(i,j))
              d_rr(i,j) = rbinom_approx(n = x%x(i,idxR), p = p%prob_mig(i,j))
              ! nan-check
              if (isnan(d_ss(i,j))) then
                write(ss,*) i, j, x%x(i,idxS), p%prob_mig(i,j)
                call mywrite(ss,ss)
                call mystop()
              else if (isnan(d_ii(i,j))) then
                write(ss,*) i, j, x%x(i,idxI), p%prob_mig(i,j)
                call mywrite(ss,ss)
              else if (isnan(d_rr(i,j))) then
                write(ss,*) i, j, x%x(i,idxR), p%prob_mig(i,j)
                call mywrite(ss,ss)
              end if
            else
              d_ss(i,j) = 0d0
              d_ii(i,j) = 0d0
              d_rr(i,j) = 0d0
            end if
          end do
        end do

        ! パラメータによる betaのseasonalityの表現
        !   注意: 番兵くんが必要なので，p%nBetaJmp >= 2 でなければならない．
if (.false.) then
        beta(:) = p%beta
        do igr = 1, p%nGroup
          do i = 1, p%nBetaJmp
            if (p%betaJmp(i,igr)%tAt .gt. x%t) exit
          end do
          i = i - 1
          if (p%nBetaJmp.eq.1 .or. p%betaJmp(i,igr)%tAt.le. x%t .and. x%t .lt. p%betaJmp(i+1,igr)%tAt) then
          else
            write(ss,*)'Assertion failure on betaJmpAt', i, igr, p%betaJmp(i,igr)%tAt, x%t, p%betaJmp(i+1,igr)%tAt
            call mywrite(ss, ss)
          end if
  
          beta(igr) = p%beta * p%betaJmp(i,igr)%value
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 502-05-28: quick hack: Augumented state で beta を推定する場合.
else
        beta(:) = p%beta * x%rbeta(:,1)
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!ばぐっとる　beta
        do i = 1, nSite
          igr = p%groupid(i)


! 29-05-25: an adhoc introduction of scale factor on FoI
          cr = (n_present(i)/1.28D7)**p%contact_pow

          d_si = rbinom_approx(n = x%x(i,idxS)&
                       ,p = pr_rate(beta(igr)*x%x(i,idxI)/n_present(i)*cr*p%dt))
          d_ir = rbinom_approx(n = x%x(i,idxI)&
                       ,p = pr_rate(p%gamma*p%dt))

          ! nan-check of d(s->s) 
            if (isnan(d_si)) then
              tmp = beta(igr)*x%x(i,idxI)/n_present(i)*cr*p%dt
              write(ss,*) 'd_si NaN:',  x%x(i,idxI), n_present(i)
              call mywrite(ss,ss)
            end if

 
          z%x(i,idxS) = max(x%x(i,idxS) - d_si - sum(d_ss(i,:)) + sum(d_ss(:,i)),0d0)
          z%x(i,idxI) = max(x%x(i,idxI) + d_si - d_ir - sum(d_ii(i,:)) + sum(d_ii(:,i)),0d0)


          z%x(i,idxR) = max(x%x(i,idxR) + d_ir - sum(d_rr(i,:)) + sum(d_rr(:,i)), 1d0)
          z%x(i,idxS2I) = d_si
        end do

        z%x0 = x%x0 ! for smoothing of initial condition

        z%t = x%t + p%dt

        ! for augmented state to est parames
        lag = size(x%rbeta,2)
        z%rbeta(:,2:lag) = x%rbeta(:,1:(lag-1))  ! for smoothing 
        if (x%sd_rbeta.gt.0) then
          do i = 1, nSite
            z%rbeta(i,1) = x%rbeta(i,1) + rgauss()*x%sd_rbeta
            if (z%rbeta(i,1).lt.0) z%rbeta(i,1) = 0.1
          end do
        else
          z%rbeta(:,1) = x%rbeta(:,1)
        end if
        z%sd_rbeta = x%sd_rbeta
    end function

    ! "flat sir" model as a competitor of the meta-population model
    elemental function advance_aggregated(p,x) result (z)
        type (state), intent(in) :: x
        type (state) :: z
        type (param), intent(in) :: p

        real(8), dimension(nSite,nSite) :: d_ss, d_ii, d_rr
        real(8) :: d_si, d_ir

        real(8) :: n_present(nSite), tmp, beta(2)
        real(8) :: x_agg(4), i_agg_gr(2), beta_times_i, n_agg, scal, d_si_scal, d_ir_scal

        integer :: i, j, igr
        character(len=256) ss

        forall (i=1:nSite) n_present(i) = sum(x%x(i,1:3))
        forall (j=1:4) x_agg(j) = sum(x%x(:,j))
        n_agg = sum(n_present)

        beta(:) = p%beta
        do igr = 1,2 
          do i = 1, p%nBetaJmp
            if (p%betaJmp(i,igr)%tAt .gt. x%t) exit
          end do
          i = i - 1

          if (p%betaJmp(i,igr)%tAt.le. x%t .and. x%t .lt. p%betaJmp(i+1,igr)%tAt) then
          else
            call mywrite('Assertion failure on betaJmpAt', ss)
          end if
  
          beta(igr) = p%beta * p%betaJmp(i,igr)%value
        end do

        i_agg_gr(1) = sum(x%x(:,idxI), mask = p%groupid .eq. 1)
        i_agg_gr(2) = sum(x%x(:,idxI), mask = p%groupid .eq. 2)
        beta_times_i = beta(1)*i_agg_gr(1) + beta(2)*i_agg_gr(2)

        d_si = rbinom_approx(n = x_agg(idxS)&
                            ,p = pr_rate(beta_times_i/n_agg*p%dt))
        d_ir = rbinom_approx(n = x_agg(idxI)&
                            ,p = pr_rate(p%gamma*p%dt))
 
        ! distribute transitters d_si, d_ir to each prefecture,
        ! proportionally to its population.
        do i = 1, nSite
          scal = n_present(i)/n_agg
          d_si_scal = d_si*scal
          d_ir_scal = d_ir*scal
          z%x(i,idxS) = x%x(i,idxS) - d_si_scal
          z%x(i,idxI) = max(x%x(i,idxI) + d_si_scal - d_ir_scal, 0d0)
          z%x(i,idxR) = x%x(i,idxR)
          z%x(i,idxS2I) = d_si_scal
        end do

        z%t = x%t + p%dt
    end function

    !negative log-likelihood
    elemental real(8) function logOC(p, y, x)
        type (param), intent(in) :: p
        type (state), intent(in) :: x
        type (obs), intent(in)   :: y
        integer :: i

        logOC = 0d0
        if (p%useBinomObsModel) then
          do i = 1,nSite
            logOC = logOC &
              - max(-7d2, lnBinom(x = y%j(i), n = x%x(i,idxS2I)/p%probBinomObs, p = p%probBinomObs))
            !use log(e**(-700)) as log(0) = -Infinity for numerical convenience.
          end do
        else
          do i = 1,nSite
            logOC = logOC &
              - lnPo(x = y%j(i), m = max(1D0,x%x(i,idxS2I)))
          end do
        end if
    end function

    function postObs(p,y,x,ipar) result (z)
        type (param), intent(inout) :: p
        type (state), intent(in)    :: x
        type (obs), intent(in)      :: y
        integer, intent(in)         :: ipar
        type (state) :: z

        z = x

        ! 時変ODを使う場合 (粒子番号1で確率表を更新する)
        if (associated(p%od_seq).and.ipar.eq.1) then
          p%i_od_seq = p%i_od_seq + 1
          if (inIdxRange(p%i_od_seq, p%od_seq)) then
            call mk_migration_prob(p, p%od_seq(p%i_od_seq)%od, p%tbl_pop, p%dt)
          else
            p%i_od_seq = p%i_od_seq - 1
          end if
        end if

        ! for a meanwhile, Rt is overwritten by that of the data!
        !   * That y%rbeta(i) is NaN means N/A.
        where (.not.isnan(y%rbeta(:))) z%rbeta(:,1) = y%rbeta(:)
        z%x(:,idxS2I) = 0D0
    end function

    real(8) function timeOf (x)
        type (state) :: x
        timeOf = x%t
    end function

    real(8) function obsTimeOf (y)
        type (obs) :: y
        obsTimeOf = y%t
    end function

    subroutine calc_quantile(xqs, xs, qs)
      type(state), intent(out) :: xqs(:)
      type(state), intent(in) :: xs(:)
      real(8), intent(in) :: qs(:)

      real(8) :: buf(size(xs)), qbuf(size(qs))
      integer :: i, j

      if (size(xqs).lt.size(qs)) then
        write(*,*) 'calc_quantile: buffer size of xqs is not enough.'
        stop
      end if

      do i = 1, nSite
        do j = 1, 4
          buf(:) = xs(:)%x(i,j)
          call mpi_real_quantile(qbuf, buf, qs)
          xqs(:)%x(i,j) = qbuf(:)
        end do
      end do
      do j = 1, size(xs(1)%rbeta,2)   ! = lag
        do i = 1, size(xs(1)%rbeta,1) ! = nSite
          buf(:) = xs(:)%rbeta(i,j)
          call mpi_real_quantile(qbuf, buf, qs)
          xqs(:)%rbeta(i,j) = qbuf(:)
        end do
      end do
      xqs(:)%t = xs(1)%t
    end subroutine

    subroutine calc_mean(z,xs)
      type(state), intent(out) :: z
      type(state), intent(in) :: xs(:)

      real(8) :: buf(size(xs))
      integer :: i, j

      do i = 1, nSite
        do j = 1, 4
          buf(:) = xs(:)%x(i,j)
          z%x(i,j) = mpi_real_mean(buf)
        end do
      end do
      do j = 1, size(xs(1)%rbeta,2)   ! = lag
        do i = 1, size(xs(1)%rbeta,1) ! = nSite
          buf(:) = xs(:)%rbeta(i,j)
          z%rbeta(i,j) = mpi_real_mean(buf)
        end do
      end do
      z%t = xs(1)%t
    end subroutine

    subroutine calc_rmse(rmse, y, xs)
      type(obs), intent(out) :: rmse
      type(obs), intent(in) :: y
      type(state), intent(in) :: xs(:)

      real(8) :: buf(size(xs))
      integer :: i
      do i = 1, nSite
        buf(:) = xs(:)%x(i,idxS2I)
        rmse%j(i) = mpi_real_rmse(y%j(i), buf)
      end do
      rmse%t = y%t
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! additional functions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Building the initial condition & parameters.
    ! You should properly give:
    !   x% x(e), x(i)
    !   p% rr0, rrmax, n, alpha, gamma, eps(1,2), eps(2,1), sdsys
    logical function cond(o_seq,idx,xs)
      type (obs) :: o_seq(:)
      integer :: idx
      type (state) :: xs(:)
      cond = idx .le. size(o_seq)
      !cond = xs(1)%t .lt. 1000
      !write(*,*) idx, size(o_seq)
    end function

    subroutine open_record_files(tbl_pop, outdir, nSite)
      type (table_pop) :: tbl_pop
      integer :: i, nSite
      character(len=*) :: outdir

      if (me().eq.0) then
        open(9600,file=trim(outdir)//'/filter_result.csv')
        write(9600,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9601,file=trim(outdir)//'/filter_result_S.csv')
        write(9601,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9602,file=trim(outdir)//'/filter_result_I.csv')
        write(9602,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9603,file=trim(outdir)//'/filter_result_R.csv')
        write(9603,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9604,file=trim(outdir)//'/filter_result_rmse.csv')
        write(9604,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9605,file=trim(outdir)//'/filter_result_q5.csv')
        write(9605,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9606,file=trim(outdir)//'/filter_result_S_q5.csv')
        write(9606,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9607,file=trim(outdir)//'/filter_result_I_q5.csv')
        write(9607,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9608,file=trim(outdir)//'/filter_result_R_q5.csv')
        write(9608,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9609,file=trim(outdir)//'/filter_result_q95.csv')
        write(9609,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9610,file=trim(outdir)//'/filter_result_S_q95.csv')
        write(9610,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9611,file=trim(outdir)//'/filter_result_I_q95.csv')
        write(9611,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

        open(9612,file=trim(outdir)//'/filter_result_R_q95.csv')
        write(9612,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)
      end if
    end subroutine

    subroutine recorder(xs,y)
       type (state) :: xs(:), x, xg, xqs(2)

       type (obs), optional :: y
       real(8) :: n, ng
       real(8) :: beta(2), beta_g(2)
       integer :: k, m, ierr
       integer, parameter :: is(2) = (/1,6/)

       if (.not.present(y).and.xs(1)%t.lt.728) return

       n = size(xs)
       ng = n*nproc()

      !配列のレイアウト: real(8) :: x(nSite,4)
       do m = 1, 4
         do k = 1, nSite
           x%x(k,m) = sum(xs(:)%x(k,m))
         end do
       end do

       do m = 1, size(xs(1)%rbeta,2)
         do k = 1, nSite
           x%rbeta(k,m) = sum(xs(:)%rbeta(k,m))
         end do
       end do
       call mpi_reduce(x%x    , xg%x    , nSite*4          , MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call mpi_reduce(x%rbeta, xg%rbeta, size(xs(1)%rbeta), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       ! mpi_reduceはループの外に出せるはず．

       if (me().eq.0) then
           x%t = xs(1)%t
           x%x = xg%x/ng
           x%rbeta = xg%rbeta/ng
           write(9600,'(48(G16.8,","))') x%t, x%x(:,idxS2I)
           write(9601,'(48(G16.8,","))') x%t, x%x(:,idxS)
           write(9602,'(48(G16.8,","))') x%t, x%x(:,idxI)
           write(9603,'(48(G16.8,","))') x%t, x%x(:,idxR)
       end if

if (calc_CI) then
       call calc_quantile(xqs, xs, (/0.025d0,0.975d0/))

        ! q5%
       write(9605,'(48(G16.8,","))') x%t, xqs(1)%x(:,idxS2I)
       write(9606,'(48(G16.8,","))') x%t, xqs(1)%x(:,idxS)
       write(9607,'(48(G16.8,","))') x%t, xqs(1)%x(:,idxI)
       write(9608,'(48(G16.8,","))') x%t, xqs(1)%x(:,idxR)
       ! q95%
       write(9609,'(48(G16.8,","))') x%t, xqs(2)%x(:,idxS2I)
       write(9610,'(48(G16.8,","))') x%t, xqs(2)%x(:,idxS)
       write(9611,'(48(G16.8,","))') x%t, xqs(2)%x(:,idxI)
       write(9612,'(48(G16.8,","))') x%t, xqs(2)%x(:,idxR)

       k = size(x%rbeta)
       !write(9700,'(4(G16.8,","))',advance='no' )  x%t, x%rbeta(1), xqs(1)%rbeta(1), xqs(2)%rbeta(1)
       !write(9700,'(4(G16.8,","))',advance='yes')  x%t-k+1, x%rbeta(k), xqs(1)%rbeta(k), xqs(2)%rbeta(k)

       ! on a meanwhile, the uncertainty is ignored!
       write(9700,'(48(G16.8,","))')  x%t, x%rbeta(:,1)
else
       write(9700,'(48(G16.8,","))')  x%t, x%rbeta(:,1)
end if

    end subroutine

    subroutine recorder_no_out(xs,y)
       type (state) :: xs(:)
       type (obs), optional :: y
    end subroutine
 
     subroutine read_pop(table,filename)
       type(table_pop) :: table
       integer :: id, i
       character(len=*) :: filename
       character(len=16) :: nameJP
       open(51,file=filename)
       read(51,*)
       do i = 1, nSite
       read(51,*) id, nameJP, table%name(i), table%total(i), table%man(i), table%woman(i)
       end do
       close(51)
     end subroutine

     subroutine read_matrix_od(table,filename)
       real(8), dimension(nSite,nSite) :: table
       character(len=*) :: filename
       character(len=16) :: nameJP
       integer :: i
       open(51,file=filename)
       read(51,*)
       do i = 1, nSite
         read(51,*) nameJP, table(i,:)
       end do
       close(51)
     end subroutine
  
     integer function seek_matrix_od_set(set,yy,mm,dd)
       type(table_od_time) :: set(:)
       integer :: yy, mm, dd
       integer :: iday, jo, kd, nday
       nday = size(set)
       do iday=1,nday
         seek_matrix_od_set = iday
         if (set(iday)%yy.eq.yy.and.&
            set(iday)%mm.eq.mm.and.&
            set(iday)%dd.eq.dd) return
       end do
       seek_matrix_od_set = -1
     end function

     subroutine read_matrix_od_set(set,filename)
       type(table_od_time), allocatable :: set(:)
       character(len=*) :: filename
       character(len=16) :: nameJP
       character(len=256) :: line
       integer :: idx, iday, jsite,nday
       open(51,file=filename)
       nday = 915  ! tentatively hard-coded
       read(51,*) ! skip the header
       allocate(set(nday))
       do iday=1,nday
         do jsite=1,nSite
           read(51,*,err=9900) idx, set(iday)%yy, set(iday)%mm, set(iday)%dd, nameJP, set(iday)%od(jsite,:)
           if (jsite.eq.1 .and. nameJP.ne.'北海道') then
             write(*,*) 'read_matrix_od_set: alignment is broken:', jsite, trim(nameJP)
           end if
         end do
       end do
       return
9900 continue
       read(51,'(a256)') line
       write(*,*) 'read_matrix_od_set: format error:', trim(line)
       stop 9900
     end subroutine

     ! mitigate a huge number of N/As in od set
     subroutine repair_matrix_od_set(set)
       type(table_od_time), allocatable :: set(:)
       integer :: iday, jo, kd, nday
       nday = size(set)
       do iday=1,nday
         do jo=1,nSite
           do kd=1,nSite
             if (isnan(set(iday)%od(jo,kd))) &
               set(iday)%od(jo,kd)= set(iday)%od(kd,jo)
             if (isnan(set(iday)%od(jo,kd))) &
               set(iday)%od(jo,kd)= 0D0
           end do
         end do
       end do
     end subroutine

     pure subroutine mk_migration_prob(p, matrix_od, tbl_pop, dt)
       ! matrix_od(src, dst): OD-matrix
       type(param), intent(inout) :: p
       real(8), intent(in), dimension(nSite,nSite) :: matrix_od
       type(table_pop), intent(in) :: tbl_pop
       real(8), intent(in) :: dt

       real(8) :: p_esc0, p_esc, redu
       integer :: i, j

       forall(i=1:nSite,j=1:nSite) 
         p%rho(i,j) = matrix_od(i,j) / tbl_pop%total(i)
       end forall
       forall(i=1:nSite) p%rho(i,i) = 0d0

       ! add following fileds to p
       do i = 1, nSite
         p_esc0 = 0d0
         p_esc  = 1d0
         do j = 1, nSite
           p%prob_mig_pre(i,j) = pr_rate(p%conn_redu*p%rho(i,j)*dt)
           p_esc0 = p_esc0 + p%prob_mig_pre(i,j) 
           p_esc  = p_esc * (1d0 - p%prob_mig_pre(i,j)) 
         end do
         p_esc = 1d0 - p_esc
         redu = (p_esc + 1d-300)/(p_esc0 + 1d-300)
         p%prob_mig(i,:) = p%prob_mig_pre(i,:)*redu
       end do

!       if (me().eq.0) then
!         open(unit=9999, file='rho.csv')
!         do i = 1, nSite
!           write(9999, '(46(G14.6,","),G14.6)') p%rho(i,:)
!         end do
!         close(unit=9999)
!
!         open(unit=9999, file='prob_mig_pre.csv')
!         do i = 1, nSite
!           write(9999, '(46(G14.6,","),G14.6)') p%prob_mig_pre(i,:)
!         end do
!         close(unit=9999)
!
!         open(unit=9999, file='prob_mig.csv')
!         do i = 1, nSite
!           write(9999, '(46(G14.6,","),G14.6)') p%prob_mig(i,:)
!         end do
!         close(unit=9999)
!       end if
     end subroutine

     ! 観察データ(新規感染者数時系列)の読み込み
     subroutine read_obs(table,filename)
       type(obs), allocatable :: table(:)
       character(len=*) :: filename
       real(8) :: gomi
       integer :: i, nLine

       open(51,file=filename)
       read(51,*)
       read(51,*)
       nLine = 0
       do
         read(51,*,END=1000)
         nLine = nLine + 1
       end do
1000   allocate(table(nLine)) 
       rewind(51)
       read(51,*)
       read(51,*)
       do i = 1, nLine
         read(51,*) table(i)%t, gomi, table(i)%j 
       end do
       close(51)
     end subroutine

     ! 観察データ + 事前推定値(Rt) の読み込み
     subroutine read_obs_and_preestimate(table,file_obs,file_est)
       type(obs), allocatable :: table(:), table2(:)
       character(len=*) :: file_obs, file_est
       integer :: nLine, nLine2, i

       call read_obs(table,file_obs)
       call read_obs(table2,file_est)

       nLine  = size(table)
       nLine2 = size(table2)

       if (nLine .ne. nLine2) then
         write(0,*) 'obs file (',nLine,'lines) and est file (',nLine2,'lines) are not consistent'  
         stop
       end if

       do i = 1, nLine
         table(i)%rbeta(:) = table2(i)%j(:)
       end do
       deallocate(table2)
     end subroutine

     subroutine read_parameters(p, parfile)
       type (param) :: p
       character(len=*) :: parfile
       integer :: ierr
       character(len=32) :: str(2)

        open(103, file = parfile, status = "old")

        nullify(p%od_seq)

        do
          read(103,*,iostat=ierr) str(1), str(2)
          if (ierr < 0) then 
            exit
          else if (ierr > 0) then
            read(103,*) str(1)
            write(0,*) 'model.f90:read_parameters: failed to read a parameter table in line:' , str(1)
      !      stop
          else
            select case (str(1))
              case ("muS"); read(str(2),*) p%muS
              ! ToDo: Index shifting by one is ugly! Correct later.
              case ("beta1Jmp1"); read(str(2),*) p%betaJmp(1,1)%value
              case ("beta1JmpAt1"); read(str(2),*) p%betaJmp(1,1)%tAt
              case ("beta1Jmp2"); read(str(2),*) p%betaJmp(2,1)%value
              case ("beta1JmpAt2"); read(str(2),*) p%betaJmp(2,1)%tAt
              case ("beta1Jmp3"); read(str(2),*) p%betaJmp(3,1)%value
              case ("beta1JmpAt3"); read(str(2),*) p%betaJmp(3,1)%tAt
              case ("beta1Jmp4"); read(str(2),*) p%betaJmp(4,1)%value
              case ("beta1JmpAt4"); read(str(2),*) p%betaJmp(4,1)%tAt
              case ("beta1Jmp5"); read(str(2),*) p%betaJmp(5,1)%value
              case ("beta1JmpAt5"); read(str(2),*) p%betaJmp(5,1)%tAt
              case ("beta1Jmp6"); read(str(2),*) p%betaJmp(6,1)%value
              case ("beta1JmpAt6"); read(str(2),*) p%betaJmp(6,1)%tAt
              case ("beta1Jmp7"); read(str(2),*) p%betaJmp(7,1)%value
              case ("beta1JmpAt7"); read(str(2),*) p%betaJmp(7,1)%tAt
              case ("beta1Jmp8"); read(str(2),*) p%betaJmp(8,1)%value
              case ("beta1JmpAt8"); read(str(2),*) p%betaJmp(8,1)%tAt

              case ("beta2Jmp1"); read(str(2),*) p%betaJmp(1,2)%value
              case ("beta2JmpAt1"); read(str(2),*) p%betaJmp(1,2)%tAt
              case ("beta2Jmp2"); read(str(2),*) p%betaJmp(2,2)%value
              case ("beta2JmpAt2"); read(str(2),*) p%betaJmp(2,2)%tAt
              case ("beta2Jmp3"); read(str(2),*) p%betaJmp(3,2)%value
              case ("beta2JmpAt3"); read(str(2),*) p%betaJmp(3,2)%tAt
              case ("beta2Jmp4"); read(str(2),*) p%betaJmp(4,2)%value
              case ("beta2JmpAt4"); read(str(2),*) p%betaJmp(4,2)%tAt
              case ("beta2Jmp5"); read(str(2),*) p%betaJmp(5,2)%value
              case ("beta2JmpAt5"); read(str(2),*) p%betaJmp(5,2)%tAt
              case ("beta2Jmp6"); read(str(2),*) p%betaJmp(6,2)%value
              case ("beta2JmpAt6"); read(str(2),*) p%betaJmp(6,2)%tAt
              case ("beta2Jmp7"); read(str(2),*) p%betaJmp(7,2)%value
              case ("beta2JmpAt7"); read(str(2),*) p%betaJmp(7,2)%tAt
              case ("beta2Jmp8"); read(str(2),*) p%betaJmp(8,2)%value
              case ("beta2JmpAt8"); read(str(2),*) p%betaJmp(8,2)%tAt

              case ("useAggregatedModel"); 
                read(str(2),'(I1)') p%useAggregatedModel
                write(*,*) 'WARNING: aggregated model is activated.'
              case default; write(*,*) 'WARNING: unknown parameter name '//trim(str(2))//' is ignored.'
            end select
          end if
        end do




     end subroutine

!    function readObs(fn1, fn2, n, t0)
!        type (obs) :: ys(200)
!        character(*) :: fn1, fn2
!        type (obs), pointer :: readObs(:)
!        real(8) :: n(2), t0
!        integer :: m, mm
!        open(51, file=fn1)
!        open(52, file=fn2)        
!        read(51,*) n(1)
!        read(52,*) n(2)
!        m = 1
!        do
!            read(51,*, END=999) ys(m)%tb, ys(m)%te, ys(m)%j(1)
!            read(52,*, END=999) ys(m)%tb, ys(m)%te, ys(m)%j(2)
!            m = m + 1
!        end do
!999     continue
!        close(51)
!        close(52)
!        mm = m - 1
!        t0 = ys(1)%tb
!
!        m = 1
!        ys(m+1:mm+1) = ys(m:mm) 
!        ys(m)%kind = OBS_START
!        m = m + 1
!        mm = mm + 1
!        do 
!            if (m > mm - 1) exit
!            if (ys(m)%te .eq. ys(m+1)%tb) then   
!                ys(m)%kind = OBS_ENDSTART
!                m = m + 1
!            else
!                ys(m)%kind = OBS_END
!                m = m + 1
!                ys(m+1:mm+1) = ys(m:mm) 
!                ys(m)%kind = OBS_START
!                m = m + 1
!                mm = mm + 1
!            end if
!        end do
!        ys(m)%kind = OBS_END
!        allocate(readObs(mm))
!        readObs(1:mm) = ys(1:mm)
!    end function
end module
