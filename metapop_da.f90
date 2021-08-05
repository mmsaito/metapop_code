! rubella_03b: PFing to a configuration of the parameters
!   - the configuration is read from a file.
!   - the parameters are identical to those of rubella_06
!     - there are four jumps in beta(t).
! ------------------------------------------------------
program rubella_03_lik2_step7
  use Model
  use PF
  use PF_Misc

  type(table_pop), target :: tbl_pop, tbl_vac
  real(8), dimension(nSite,nSite) :: matrix_od
  type(table_od_time), target, allocatable :: od_seq(:)
  type(obs), dimension(:), allocatable :: ys
  type(param) :: p, ps_dummy(2)
! integer, parameter :: npar1 = 64  ! short @ 128procs
!  integer, parameter :: npar1 = 1024  ! standard @ 128procs

!  integer, parameter :: npar1 = 2048   ! R02-04時点現在，通常
  integer, parameter :: npar1 = 2730   ! R02-04時点現在，通常 (3CORES)
!  integer, parameter :: npar1 = 4096 
! integer, parameter :: npar1 = 16384  ! R02-04時点現在，長時間計算

!  integer, parameter :: npar1 = 8192 ! big @ 128procs, standard @ 8proc
!  integer, parameter :: npar1 = 65536 ! super big @ 128@procs
  integer, parameter :: max_n_inf0 = 5
  type (state) :: x, xs(npar1)
  real(8) :: R0, gamma_reduced
  integer :: ierr, sz
  integer, allocatable :: gt(:)
  real(8) :: muS, nlogL
  real(8), dimension(nSite) :: n_out
  character(len=32) :: str(2)
  character(len=256) :: line, parfile
  real(8) :: npopall

  character(len=256) :: outdir, vacfile, inffile
  integer :: irow_inffile = 1


  ! 学習するデータの上限(時点による)と予測する時刻(物理時刻)
  ! 実装上の都合で別の単位を用いる．
  integer :: n_tp_learn
  real(8) :: t_upto

  call mpi_init(ierr)   


  ! フラグ: Confidence Intervalを記録するか: Module Model で定義
  calc_CI = .true.  

  ! ランク毎に異なる乱数シードの設定
  call random_seed(size=sz)
  allocate(gt(sz))
  call random_seed(get=gt)
  gt = gt + me()
  call random_seed(put=gt)
  deallocate(gt)

!  open(1234,file='pois'//trim(sI(me())) //'.csv')
!  do i = 1,100000
!    write(1234,*) rpoisson(lambda = 1d-3)
!  end do
  ! byte-size of parameter structure
  isizeof_param = loc(ps_dummy(2)) - loc(ps_dummy(1))


  ! 感染者数時系列 [+ Rtの事前推定値]
  ! call read_obs(ys)
  call read_obs_and_preestimate(ys, "COVID-19-pref.csv", "Rt_pref.csv" )

  n_tp_upto  = size(ys) ! getopt(p)で上書き可能

  tbl_vac%total(:) = 0d0
  inffile = ""
  irow_infile = 1

  call getopt(p)

  ! 都道府県人口データ
  call read_pop(tbl_pop, "total_pop_japan_pref_2010_nowhole.csv")

  ! OD-matrix
  call read_matrix_od(matrix_od, "transport_all_nowhole_per_day.csv")
!!!!!!!!!!!!!!!!quick hack ここて゛ODを切り替える !!!!!!!!!!!!!!!!!!!!
  !call read_matrix_od(matrix_od, "transport_keikai_reduce#0.5.csv")
  !call read_matrix_od(matrix_od, "transport_keikai_reduce#0.05.csv")
  !call read_matrix_od(matrix_od, "transport_keikai_reduce#0.001.csv")
  !call read_matrix_od(matrix_od, "transport_keikai_reduce#0.csv")

  ! timeseries of OD-matrix
  call read_matrix_od_set(od_seq,'ODflow_pref_daily_nan_2018-2021.csv')
  call repair_matrix_od_set(od_seq)

  R0      = 1d0
  p%muS   = 1d0
  p%dt    = 1d0
  T_inf   = 4d0
  p%gamma = -log(1 - p%dt/T_inf)/p%dt
  ! Note: #[recovery] ~ Bin(#[infective], p = dt/T_inf) will be realized.

  p%beta  = R0/T_inf

  p%nGroup   = nSite
  p%nBetaJmp = 1

  p%betaJmp(1,:)%tAt    = -1d300
  p%betaJmp(1,:)%value  = 0.0d0

  ! xx all prefs belong to the single group 
  ! all groups are now a singleton!
  p%groupid(:) = 1  
  forall (i=1:nSite) p%groupid(i) = i

  ! parameter set: read from a file >>>>>>>>>
  !   * muS, betaJmp(,)を設定する．
  call read_parameters(p, parfile)

  call mpi_barrier(MPI_COMM_WORLD, ierr)


  if (me().eq.0) then
    open(9599,file=trim(outdir)//'/data_referred.csv')
    write(9599,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)
    do i=1,size(ys)
      write(9599,'(48(G16.8,","))') ys(i)%t,ys(i)%j(:)
    end do
    close(9599)
  end if

  if (me().eq.0) then
  open(9999,file=trim(outdir)//"/config_remark_predict.txt") 
   write(9999,*) "#npar*nproc = ", npar1,"*", nproc(), "=", npar1*nproc()
   write(9999,*) "T_inf = ",  T_inf
   write(9999,*) "1/gamma = ", 1/p%gamma 
   do j = 1, 2
   do i = 1, p%nBetaJmp
     write(9999,*) "t = ", p%betaJmp(i,j)%tAt, "- beta = ", p%betaJmp(i,j)%value
   end do
   end do
   write(9999,*) "conn_redu = ", p%conn_redu
   write(9999,*) "vac = ", tbl_vac%total(:)
   if (p%useBinomObsModel) write(9999,*) 'using binomial observation model with p = ', p%probBinomObs
   if (p%useAggregatedModel) write(9999,*) 'usign aggregated model' 
  close(9999)
  end if


  ! Waste a R.N. My env. gives me a super outlier (5.36 from N(0,1)) at the first time.
  do i = 1,5
    gomi = rgauss()
    ! Oh! this write sentense is necessary since you cheat purity of random num. gen!
    write(str(1),*) gomi
  end do

  call mk_migration_prob(p,matrix_od, tbl_pop, p%dt)

! Uncomment out below, if the root modifies the parameters.
! call mpi_bcast(p, isizeof_param, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Way of choosing initial infected cases
  ! (1) estimated from data
  !call make_initial_set_by_smoothing(xs, tbl_pop)


  ! (2) chosen from a file (one setting by one line)
  call read_init_inf(x, irow_inffile, inffile)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 502-05-28: quick hack: Augumented state で beta を推定する場合.
  x%rbeta    = 1.0D0
  x%sd_rbeta = 0.1D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! (2.5) futher initialization for time-varying OD-matrix
  p%i_od_seq = seek_matrix_od_set(od_seq,2020,01,14)
  p%od_seq   => od_seq
  p%tbl_pop  => tbl_pop
  call mk_migration_prob(p, p%od_seq(p%i_od_seq)%od, p%tbl_pop, p%dt)

  ! (3) setting by hand writing
!  do i = 1, 47
!    x%x(i,idxS) = aint(p%muS*tbl_pop%total(i) - tbl_vac%total(i))
!    x%x(i,idxR) = tbl_pop%total(i) - x%x(i,idxS)
!    x%x(i,idxI) = 0d0
!    x%x(i,idxS2I) = 0d0
!  end do
!
!  ! index 13 == Tokyo
!  i = 13 
!
!  if (.not. p%useAggregatedModel) then
!    i0 = 5d0
!    i = 13
!    x%x(i,idxR) = x%x(i,idxR) - i0
!    x%x(i,idxI) = i0
!    x%t = 0d0
!
!    i = 27
!    x%x(i,idxR) = x%x(i,idxR) - i0
!    x%x(i,idxI) = i0
!    x%t = 0d0
!  else
!    i0 = 10d0
!    npopall = sum(tbl_pop%total(1:47))
!    do i = 1,47
!      scal = tbl_pop%total(i)/npopall
!      x%x(i,idxR) = x%x(i,idxR) - i0*scal
!      x%x(i,idxI) = i0*scal
!      x%t = 0d0
!    end do
!  end if


  !write(*,*) '*'
  if (me().eq.0) write(*,*) 'n_tp_learn = ', n_tp_learn

  if (me().eq.0) then
    open(9600,file=trim(outdir)//'/filter_result.csv')
    write(9600,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

    open(9601,file=trim(outdir)//'/filter_result_S.csv')
    write(9601,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

    open(9602,file=trim(outdir)//'/filter_result_I.csv')
    write(9602,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)

    open(9603,file=trim(outdir)//'/filter_result_R.csv')
    write(9603,'(48(a,","))') 'Time', (trim(tbl_pop%name(i)),i=1,nSite)
if (calc_CI) then
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

    open(9700,file=trim(outdir)//'/filter_result_Rt.csv')
    write(9700,*) 'Time,Rt,RtL,RtU,Time2,RtS,RtSL,RtSU'
end if
  end if


  xs(:) = x
  call filterApp(p, cond, recorder, ys(1:n_tp_learn), xs, nLogL)

  if (me().eq.0) write(*,*) 'prediction starts at t = ', xs(1)%t


  write(*,*) 'WARNING: OD has been replaced before prediction.'
  call mk_migration_prob(p, matrix_od, tbl_pop, p%dt)

  call evaluate_prediction_performance(p, ys(n_tp_learn+1:), xs, calc_CI, t_upto)

  if (me().eq.0) then
    do i = 9600,9612
      close(UNIT=i)
    end do
  end if

  if (me().eq.0) then
    open(9999,file=trim(outdir)//'/logL.txt',access='append')
    write(9999,*) p%conn_redu, -nlogL
    close(9999)
  end if

  deallocate(ys)
  deallocate(od_seq)
9999 CONTINUE
  call mpi_finalize(ierr)   
contains
  subroutine make_initial_set_by_smoothing(xs, tbl_pop)
    type (state) :: xs(:)
    type (table_pop) :: tbl_pop
    integer :: i, j, k, n_inf0, ipar
    integer :: ns_dup(nSite)
    integer, allocatable :: indexes(:), initinf(:,:,:), initinf1(:,:)


    if (.not. p%useAggregatedModel) then

      ns_dup = nint(10.0 * tbl_pop%total / minval(tbl_pop%total))
      allocate(indexes(sum(ns_dup)))
      indexes = (/ ((i, j=1, ns_dup(i)), i = 1, nSite) /)

      do ipar = 1, size(xs)
        !n_inf0 = 1 + int(rand()*max_n_inf0)
        n_inf0 = max_n_inf0

        x%t = 0d0
        do i = 1, nSite
          x%x(i,idxS) = aint(p%muS*tbl_pop%total(i) - tbl_vac%total(i))
          x%x(i,idxR) = tbl_pop%total(i) - x%x(i,idxS)
          x%x(i,idxI) = 0d0
          x%x(i,idxS2I) = 0d0
        end do

        do k = 1, n_inf0
          j = 1 + int(rand()*size(indexes))
          i = indexes(j)

          x%x(i,idxR) = x%x(i,idxR) - 1
          x%x(i,idxI) = x%x(i,idxI) + 1
          x%t = 0d0
        end do
        x%x0 = x%x    ! embed init cond to sys vars
        xs(ipar) = x
      end do

      call open_record_files(tbl_pop, outdir, nSite)

      call write_initinf_set(xs, tbl_pop, trim(outdir)//'/init_inf_prior.csv')
      call filterApp(p, cond, recorder_initinf_mid, ys, xs, nLogL)
      call write_initinf_set(xs, tbl_pop, trim(outdir)//'/init_inf_post.csv')
    else
      write(0, *) 'make_initial_set_by_smoothing: unimplemented for aggr model!'
      stop
    end if
  end subroutine

  subroutine recorder_initinf_mid(xs,y)
    type(state) :: xs(:)
    type(obs), optional :: y
    character(len=128) :: filepath

    write(filepath,'(a,"/init_inf_summary_mid_t#",I4.4,".csv")') trim(outdir), int(y%t)
    call write_initinf_summary(xs, tbl_pop, max_n_inf0, filepath)

    !call recorder(xs,y)
    !write(filepath,'(a,"/init_inf_mid_t#",I4.4,".csv")') trim(outdir), int(y%t)
    !call write_initinf_set(xs, tbl_pop, filepath)
  end subroutine

  subroutine write_initinf_summary(xs, tbl_pop, max_ninf, filename)
    type (state) :: xs(:)
    type (table_pop) :: tbl_pop
    character(len=*) :: filename
    integer :: max_ninf
    integer :: counts(max_ninf+1, nSite)
    integer :: i, j, ipar, ninf0

    counts = 0
    do ipar = 1, size(xs)
      do i = 1, nSite
        ninf0 = xs(ipar)%x0(i,idxI) + 1
        counts(ninf0,i) = counts(ninf0,i) + 1
      end do
    end do

    call mpi_reduce(counts, counts, size(counts), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (me().eq.0) then
      open(9650, file=filename)
      write(9650,'(48(",",a))') (trim(tbl_pop%name(i)),i=1,nSite), "Whole"
      do j = 0, max_ninf
        write(9650,'(I3,48(",",I12))') j, (counts(j+1,i), i = 1, nSite)
      end do
      close(9650)
    end if
  end subroutine

  subroutine read_init_inf(x, idx_row, filename)
    type (state) :: x
    integer :: i
    real(8) :: initinf(nSite), initinf_sum, scal, npopall
    character(len=*) :: filename

    open(9650, file=filename) 

    read(9650,*) ! skip the header line
    do i = 1, idx_row - 1
      read(9650,*) 
    end do

    read(9650, *) initinf
    close(9650)

    ! gather once and distribute when using aggregated model
    if (p%useAggregatedModel) then
      initinf_sum = sum(initinf)
      npopall = sum(tbl_pop%total(1:47))
      do i = 1, nSite
        scal = tbl_pop%total(i)/npopall
        initinf(i) = initinf_sum * scal
      end do
    end if

    do i = 1, nSite
      x%x(i,idxS) = aint(p%muS*tbl_pop%total(i) - tbl_vac%total(i))
      x%x(i,idxR) = tbl_pop%total(i) - x%x(i,idxS) - initinf(i)
      x%x(i,idxI) = initinf(i)
      x%x(i,idxS2I) = 0d0
    end do 
  end subroutine

  subroutine write_initinf_set(xs, tbl_pop, filename)
    type (state) :: xs(:)
    type (table_pop) :: tbl_pop
    character(len=*) :: filename
    integer, allocatable :: initinf(:,:,:), initinf1(:,:)

    allocate(initinf1(nSite,size(xs)))
    forall (i=1:nSite, ipar=1:size(xs)) initinf1(i,ipar) = xs(ipar)%x0(i,idxI)

    allocate(initinf(nSite, size(xs), 0:(nproc()-1)))
    !Note: This allocation should be necessary only in the root note. 
    !      But necessary to works well with boundary checking ...

      ! THIS IS INVALID. DATA isnot continuous
    call mpi_gather(initinf1(1,1), nSite*size(xs), MPI_INTEGER&
                   ,initinf(1,1,me()), nSite*size(xs), MPI_INTEGER&
                   , 0, MPI_COMM_WORLD, ierr)

    if (me().eq.0) then
      if (ierr.ne.0) then
        write(0, *) 'fail to gather init cond in each proc.'
        stop
      end if
     
      open(9650, file=filename)
      write(9650,'(48(a,","))') (trim(tbl_pop%name(i)),i=1,nSite), "Whole"
      do j = 0, nproc()-1
        do ipar = 1, size(xs)
          write(9650, '(48(I3,","))') initinf(:,ipar,j), sum(initinf(:,ipar,j))
        end do
      end do
      close(9650)
    end if
  end subroutine

  subroutine evaluate_prediction_performance(p, o_seq, xs, calc_CI, t_upto)
    implicit none
    type (param) :: p
    type (obs)   :: o_seq(:)
    type (obs)   :: rmse
    type (state) :: x, xs(:), xqs(2)!, zs(size(xs))
    integer      :: idx, i, k, n_tp_obs_last
    real(8)      :: t_obs_last, t_upto, t
    logical      :: calc_CI

    integer(2), allocatable :: newInfSeq(:,:,:) ! 2byte int. newinf should be < 30,000.

    ! what should i do?
    ! - at least the prediction up to the given time point.
    ! - the discrepancy (O/C) is evaluated e.g. in terms of RMS
    !  - posterier mean.
    ! - You can use a similar loop to filterApp, but without filtering.

    ! file unit - 9600:
    !      s i r j     rbeta(Rt)
    !-------------------------
    !mean  0 1 2 3     100
    !rmse        4
    !q5    5 6 7 8     100
    !q95   9 10 11 12  100
    idx = 1
    n_tp_obs_last = size(ys)
    t_obs_last    = obsTimeOf(ys(n_tp_obs_last))

    ! 予測フェーズでは Rtの遊走は止める.
    xs(:)%sd_rbeta = 0D0

    open(1234,file=trim(outdir)//'/newinf-'//trim(sI(me()))//'.bin',form='unformatted')

    ! 予測開始時刻 (t_obs_lastと一致するはずだが、システム変数から直接取る)
    t = xs(1)%t
    write(1234) t, p%dt

    ! 「特定の県の軌道」がすぐにとれるほうが便利なので、
    !   * (時刻,粒子番号,県番号)のレイアウトにする。
    !   * 注意: dt = 1の仮定のもと、時点(int)と時刻(real)を近藤している。
    !           シミュレーションステップが変わるとトラブルが起こる。
    allocate(newInfSeq(t:t_upto, 1:npar1, 1:nSite))
    write(1234) size(newInfSeq,1), size(newInfSeq,2), size(newInfSeq,3)
    !write(1234) (size(newInfSeq,i),i=1,3)

    forall (i=1:npar1) newInfSeq(t,i,:) = xs(i)%x(:,idxS2I)

    !do while (cond(o_seq, idx, xs)) !データがなくなるまで
    do while (timeOf(xs(1)).le. t_upto) !指定時刻まで
      xs = advance(p,xs)

      t = xs(1)%t
      forall (i=1:npar1) newInfSeq(t,i,:) = xs(i)%x(:,idxS2I)

      if (me().eq.0) write(*,*) me(), xs(1)%t
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      i = getObs(epsTime, timeOf(xs(1)), o_seq, idx)
      if (i > 0 .or. timeOf(xs(1)) > t_obs_last) then
        call calc_mean(x,xs)
if (calc_CI) then
        if (i > 0) call calc_rmse(rmse, o_seq(i), xs)
        call calc_quantile(xqs, xs, (/0.025d0,0.975d0/))
        if (me().eq.0) write(*,*) xqs(2)%t, xqs(2)%x(13,idxI)
        call flush(6)
else
        if (me().eq.0) write(*,*) xs(1)%t
        call flush(6)
end if
        !ToDo: add here proc. to eval. RMS with x and o_seq(i).
        ! and keep x in e.g. x_mean_seq(:) for plotting. 
        ! Fig(a): mean, Q5&Q95, obs, Fig(b): RMSE
        ! mean
        if (me().eq.0) then
        write(9600,'(48(G16.8,","))') x%t, x%x(:,idxS2I)
        write(9601,'(48(G16.8,","))') x%t, x%x(:,idxS)
        write(9602,'(48(G16.8,","))') x%t, x%x(:,idxI)
        write(9603,'(48(G16.8,","))') x%t, x%x(:,idxR)
if (calc_CI) then
        ! rmse
        if (i > 0) then
          write(9604,'(48(G16.8,","))') x%t, rmse%j(:)
        else
          write(9604,'(48("NA,"))')  !R's Not Available
        endif

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
! for meanwhile keeping older impl as comments
!        write(9700,'(4(G16.8,","))',advance='no' )  x%t, x%rbeta(1), xqs(1)%rbeta(1), xqs(2)%rbeta(1)
!        write(9700,'(4(G16.8,","))',advance='yes')  x%t-k+1, x%rbeta(k), xqs(1)%rbeta(k), xqs(2)%rbeta(k)
        write(9700,'(48(G16.8,","))')  x%t, x%rbeta(:,1)
else
        write(9700,'(48(G16.8,","))')  x%t, x%rbeta(:,1)
end if
        end if
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        ! Here data should NOT be reflected to the state vector!!
        !if (i > 0) xs = postObs(p, o_seq(i), xs)
      end if
    end do

    write(1234) newInfSeq
    close(1234)
    deallocate(newInfSeq)
  end subroutine

  ! currently copied & pasted from rubella_05 (MCMC code). But it is not good manner...
  subroutine getopt(p)
    type(param) :: p
    character(len=64) :: arg
    integer :: n, i

    outdir = '.'  !current directory is default.
  
    call getArg(1,str(1)) 
    read(str(1),*) n_tp_learn
    n_tp_learn = min(n_tp_learn,size(ys))

    call getArg(2,parfile) 
    write(*,*) 'using parameter set in ', trim(parfile)

    n = iargc()
    i = 3 ! a dirty hack to keep a backward compatibility (arg1 & 2)
    do while (i.le.n)
      call getarg(i, arg); i=i+1
      select case (arg)
        case ('-model')
          call getarg(i,arg); i=i+1
          select case (arg)
            case('meta'); p%useAggregatedModel = .false.
            case('aggr'); p%useAggregatedModel = .true.
            case default; write(0,*) 'unknown opt for -model:', trim(arg); stop
          end select
        case ('-conn')
          call getarg(i,arg); i=i+1
          read(arg,*) p%conn_redu
        case ('-outdir')
          call getarg(i,outdir); i=i+1
          call system("mkdir "//trim(outdir))
          write(0,*) 'results will be outputted to '//trim(outdir)
        case ('-vacfile')
          call getarg(i,vacfile); i=i+1
          call read_pop(tbl_vac, vacfile)
        case ('-inffile')
          call getarg(i,inffile); i=i+1
        case ('-at-inffile')
          call getarg(i,arg); i=i+1
          read(arg,*) irow_inffile
        case ('-noci')
          calc_CI = .false.
          if (me().eq.0) write(*,*) 'CIs are not calculated.'
        case ('-contact-pow')
          call getarg(i,arg); i=i+1
          read(arg,*) p%contact_pow 
        case ('-obs-binom')
          p%useBinomObsModel = .true.
          call getarg(i,arg); i=i+1
          read(arg,*) p%probBinomObs
          if (me().eq.0) write(0,*) 'using binomial observation model with p = ', p%probBinomObs
        case ('-t-upto')
          call getarg(i,arg); i=i+1
          read(arg,*) t_upto
        case default
          if (me().eq.0) write(0,*) 'unknown option:', trim(arg); stop
      end select
    end do 
  end subroutine
end program
