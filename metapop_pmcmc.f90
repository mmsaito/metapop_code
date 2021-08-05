!まずきちんと平均を取るコードを書こう．
! ... これはできたっぽいので，つぎにデータとのオーバープロットを
! MCMC parameter searching
!   502-04-07: covid19用に改変
program rubella_03_lik2_step7
  use Model
  use PF

  integer, parameter :: n_retry = 1
  ! to use 16384 particle as a whole:
  !integer, parameter :: n_particle1 = 16   ! nproc = 1024 
  !integer, parameter :: n_particle1 = 64  ! nproc = 256
  !integer, parameter :: n_particle1 = 128  ! nproc = 128
  !integer, parameter :: n_particle1 = 256  ! nproc = 64

  integer, parameter :: n_particle1 = 1024  ! nproc = 16

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2018-06-28  m9(^o^) QUICK HACK TO CALC ONE OF TRACKS, TRACK#6
  !integer, parameter :: n_track = 9
  !real(8) :: temps(n_track) = (/1d0, 2d0, 4d0, 8d0, 16d0, 32d0, 64d0, 128d0, 256d0/)
  integer, parameter :: n_track = 1
  real(8) :: temps(n_track) = (/64d0/)
  !
  ! C A U T I O N :  H A R M F U L C O D I N G  m9(^o^) m9(^o^)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(8) :: recip_temps(n_track)

  type(table_pop) :: tbl_pop
  real(8), dimension(nSite,nSite) :: matrix_od
  type(obs), dimension(:), allocatable :: ys
  type(param) :: p, ps(n_track), psPrev(n_track), pX, ps_dummy(2)
  type (state) :: x, xs(n_particle1)
  real(8) :: R0, gamma_reduced, crit
  real(8), volatile :: gomi
  integer :: ierr, sz
  integer, allocatable :: gt(:)
  real(8) :: nlogLs(n_track), nLogLsPrev(n_track), tmp, nLogL_retry(n_retry)
  real(8), dimension(nSite) :: n_out
  character(len=256) :: ss
  character(len=4)   :: ssint

  integer :: i_mcmc_sched, i, j, k, iswp(3)


  recip_temps = 1d0/temps

  call mpi_init(ierr)   

  ! ランク毎に異なる乱数シードの設定
  call random_seed(size=sz)
  allocate(gt(sz))
  call random_seed(get=gt)
  gt = gt + me()
  call random_seed(put=gt)
  deallocate(gt)

  ! byte-size of parameter structure
  isizeof_param = loc(ps_dummy(2)) - loc(ps_dummy(1))

  ! 都道府県人口データ
  call read_pop(tbl_pop, "total_pop_japan_pref_2010_nowhole.csv")

  ! OD-matrix
  !call read_matrix_od(matrix_od, "transport_all_nowhole.csv") <- unit 1000 per year
  call read_matrix_od(matrix_od, "transport_all_nowhole_per_day.csv")

  ! 都道府県別感染者数時系列
  call read_obs(ys,"COVID-19-pref.csv")

  if (me().eq.0) write(*,*) 'modify_par_version:', modify_par_version()

  call getopt(p)

  !R0 = 0.5d0 # 04-07実施の設定．覚えておくこと．
  R0            = 1d0
  p%muS         = 1d0
  p%dt          = 1d0
  p%gamma       = 1d0/4d0
  gamma_reduced = (1.0 - exp(-p%gamma*p%dt))/p%dt
  p%beta = R0*gamma_reduced

  ! FOI-推定区分数のデフォルト値．
  !   * modify_par.f90::init_par(p)で上書きされる．
  p%nBetaJmp = 1
  p%nGroup   = 1

  ! todo: model switching should be setted via some interface.
  if (p%useAggregatedModel) then
    write(*,*) 'WARNING: aggregated model is activated.'
  end if

  call mk_migration_prob(p,matrix_od, tbl_pop, p%dt)

  call init_par(p)
  
  ! overwrite init par by "good parameters" previously found
  !call read_parameters(p, "parameter_init_290525.csv")

  ! TENTAIVELY HARD CORDED !!!
  p%groupid(:) = 1
  p%groupid((/18,21,24,25,26,27,28,29,30,35,40,41,47/)) = 2

  ! Waste a R.N. My env. gives me a super outlier (5.36 from N(0,1)) at the first time.
  do i = 1,5
    gomi = rgauss()
    ! Oh! this write sentense is necessary since you cheat purity of random num. gen!
    !wrte(*,*) gomi
  end do

  if (me().eq.0) then
  open(9999,file="config_remark.txt")
   write(9999,*) "gamma = ", p%gamma
   write(9999,*) "gamma_reduced = ",  gamma_reduced
   write(9999,*) "R0 = ", R0, ", in other words ", "beta0 = ", p%beta

   write(9999,*) "seed parameter setting:"
   do j = 1, p%nGroup
   do i = 1, p%nBetaJmp
     write(9999,*) "t = ", p%betaJmp(i,j)%tAt, "- beta = ", p%betaJmp(i,j)%value
   end do
   end do

   write(9999,*) "conn_redu = ", p%conn_redu
   write(9999,*) "contact_pow = ", p%contact_pow
  close(9999)



!  open(9633,file='nlogL_modify_beta_281005.csv')
  do k = 1, n_track
    write(ss, "('nlog_track#',I2.2,'_',a'.csv')") k, trim(modify_par_version())
    open(9800 + k, file=ss)
  end do

  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Note you should garantee 16 = 2*p%nBetaJmp
  600 format("muS, ", 16("beta",I1,"Jmp",I1,", ","beta",I1,"JmpAt",I1,", "))
  !write(9633, 600) ((i,j,i,j, j=1, p%nBetaJmp),i=1,2)

  ! duplicate parameter for parallel tempering
  ps(:) = p 

  nLogLs = 1d300
  i_mcmc_sched = 0
  do l = 1, 16384 ! 1024 !65536 * 32

    do k = 1, n_track
      ! sample a parameter set is sampled and shared among nodes.
      if (me().eq.0) then
        psPrev(k) = ps(k)
        nLogLsPrev(k) = nLogLs(k)
        call modify_par(ps(k),i_mcmc_sched)
      end if
      call mpi_bcast(ps(k), isizeof_param, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)

      do i = 1, 47
        x%x(i,idxS) = aint(ps(k)%muS*tbl_pop%total(i))
        x%x(i,idxR) = tbl_pop%total(i) - x%x(i,idxS)
        x%x(i,idxI) = 0d0
        x%x(i,idxS2I) = 0d0
      end do

      if (.not. p%useAggregatedModel) then
        i0 = 1
        i = 14
        x%x(i,idxR) = x%x(i,idxR) - i0
        x%x(i,idxI) = i0
        x%t = 0d0
      else
        i0 = 10d0
        npopall = sum(tbl_pop%total(1:47))
        do i = 1,47
          scal = tbl_pop%total(i)/npopall
          x%x(i,idxR) = x%x(i,idxR) - i0*scal
          x%x(i,idxI) = i0*scal
          x%t = 0d0
        end do
      end if

      do i_retry = 1, n_retry
        xs(:) = x
        call filterApp(ps(k), cond, recorder_no_out, ys, xs, nLogL_retry(i_retry))
        write(*,*) me(),nLogL_retry(i_retry)
      end do
    
      nLogLs(k) = sum(nLogL_retry)/dble(n_retry)

      if (me().eq.0) then
        if (nLogLs(k).lt.nlogLsPrev(k).or.log(rand(gomi)).lt.(nLogLsPrev(k) - nlogLs(k))*recip_temps(k)) then
        else
          ps(k) = psPrev(k)
          nLogLs(k) = nLogLsPrev(k)
        end if
      end if
    end do ! end loop for track number k

    ! swap between tracks
    if (me().eq.0) then
      i = 1 + int(rand(gomi)*n_track)
      j = 1 + int(rand(gomi)*n_track) 
      crit = recip_temps(i)*(nlogLs(i) - nLogLs(j)) + recip_temps(j)*(nLogLs(j) - nLogLs(i))
      iswp(1) = i
      iswp(2) = j
      iswp(3) = rand(gomi).lt.crit
    end if
    call mpi_bcast(iswp, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (iswp(3).ne.0) then
      i = iswp(1)
      j = iswp(2)
      if (i .ne. j) then
        p    = ps(i)    ; ps(i)     = ps(j)    ; ps(j)     = p
        gomi = nlogLs(i); nlogLs(i) = nlogLs(j); nlogLs(j) = gomi
      end if
    end if
   
    if (me().eq.0) then
      do k = 1, n_track
        ! 34 + 1 = 35 = 1(muS) + 8seg*2gr*2(time&value) + 2(nlogL,nlogL/1000)
        ! 4 + 1 = 5 = 1(muS) + 1seg*1gr*2(time&value) + 2(nlogL,nlogL/1000)
        !
        ! * 502-04-08: 出力形式の変更
        !    * 旧: p%betaJmp (p%betaに対する相対値), 新: 再生産数に相当する値
        write(ssint,'(I4)') 1 + p%nBetaJmp*p%nGroup*2 + 2
        write(9800 + k,'('//ssint//'(G14.6,","),G14.6)')&
          p%muS, &
          ((ps(k)%betaJmp(i,igr)%tAt, ps(k)%betaJmp(i,igr)%value*R0, i=1,p%nBetaJmp), igr=1,p%nGroup),&
           nlogLs(k), nlogLs(k)*recip_temps(k)
        flush(9800+k)
      end do
    end if
  end do

  do k = 1, n_track
    close(9800 + k)
  end do

  call mpi_finalize(ierr)   
contains
  subroutine getopt(p)
    type(param) :: p
    character(len=64) :: arg
    integer :: n, i
   
    n = iargc()
    i = 1
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
        case ('-contact-pow')
          call getarg(i,arg); i=i+1
          read(arg,*) p%contact_pow 
        case default
          write(0,*) 'unknown option:', trim(arg); stop
      end select
    end do 
  end subroutine
! The below is a walk around for selecting the MCMC schedular.
! This is very dangerous way. Fix it as fast as possible!
! Which 'modify_par.f90' is included depends on the compile option.
! See Makefile or xxx.mak
  include 'modify_par.f90'
end program
