PROGRAM cdfjz_all
  !!======================================================================                                
  !!                     ***  PROGRAM  cdfjz_all  ***                                                 
  !!======================================================================                                
  !!  ** Purpose : Compute the surface flux of PV into the ocean using both 
  !!			   direct formulae 
  !!               (equations 9 from Maze and Marshall 2011)                                                    
  !!                                                                                                      
  !!  ** Method  : Estimate the potential density tendency, velocity 
  !!               tendency and cross outcrop flow contributions.
  !!----------------------------------------------------------------------
  !! Author: Dhruv Balwada
  !! Date : 26/6/2016 
  !!----------------------------------------------------------------------
  !! ** This code adds additional capabilities for estimating PV fluxes to 
  !! 	CDF tool box for doing analysis on NEMO output.
  !!======================================================================

  USE cdfio
  USE eos
  USE modcdfnames

  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt     ! dummy loop index                          
  INTEGER(KIND=4)                           :: npiglo, npjglo ! size of the domain                        
  INTEGER(KIND=4)                           :: npk, npt       ! size of the domain                        
  INTEGER(KIND=4)                           :: jk             ! vertical index                            
  INTEGER(KIND=4)                           :: narg, iargc    ! browse line                               
  INTEGER(KIND=4)                           :: ncoutu         ! ncid for ugeo file                        
  INTEGER(KIND=4)                           :: ierr           ! error status                              
  INTEGER(KIND=4), DIMENSION(4)             :: ipk            ! levels of output vars                     
  INTEGER(KIND=4), DIMENSION(4)             :: id_varoutu     ! varid for ugeo                            

  REAL(KIND=4)                              :: grav           ! gravity                                   
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim,tim2       ! time counter          
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1u, e2v, ff   ! horiz metrics, coriolis (f-point)         
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e2u       ! horiz metrics                             
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1f, e2f       ! horiz metrics                             
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3             ! vertic metrics                            
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: sigx 		  ! longitude latitude u-point  
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: sigy   		  ! longitude latitude v-point 
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, vn         ! velocity components at time 1                      
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un2, vn2       ! velocity components at time 2                
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fmask   		  ! mask at u and v points            
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask          ! mask at t points                          
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vort           ! vorticity            
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Jz1, Jz2, Jz3  ! surface flux of PV     
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigsurf       ! density at first time 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigsurf2       ! density at second time 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zberno,bx,by   ! bernoulli                                 

  REAL(KIND=4)                              :: zrho0          ! reference density in geos balance         
  REAL(KIND=4)                              :: zohr0          ! reference density in geos balance         
  REAL(KIND=4)                              :: deltat

  CHARACTER(LEN=256)                        :: cf_dfil, cf_dfil2     ! input density files at time 1 and 2             
  CHARACTER(LEN=256)                        :: cf_bfil        ! input bernoulli file
  CHARACTER(LEN=256)                        :: cf_ufil,cf_vfil       ! input velocity files at time 1               
  CHARACTER(LEN=256)                        :: cf_ufil2,cf_vfil2     ! input velocity files at time 2          
  CHARACTER(LEN=256)                        :: cf_out='Jz.nc'       ! outpute file name
  
  TYPE(variable), DIMENSION(4)              :: stypvaru       ! attributes for ugeo        
  LOGICAL                                   :: lchk           ! file existence flag                       
  
  !!----------------------------------------------------------------------                                
  
  CALL ReadCdfNames()

  grav  = 9.81    ! gravity                                                 
  zrho0 = 1025.d0 ! reference density                                                                        
  zohr0 = 1. / zrho0

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfjz_all'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the surface flux of PV into the ocean using both '
     PRINT *,'       direct and bulk formulae. (refer to Maze and Marshall 2011) '
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'      Rho_file 1, Rho_file 2, Bernoulli_file, U_file1, V_file1, U_file2, V_file 2 '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fmsk),' ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - netcdf file : ', TRIM(cf_out)
     STOP
  ENDIF


  CALL getarg(1, cf_dfil)
  CALL getarg(2, cf_dfil2)
  CALL getarg(3, cf_bfil)
  CALL getarg(4, cf_ufil)
  CALL getarg(5, cf_vfil)
  CALL getarg(6, cf_ufil2)
  CALL getarg(7, cf_vfil2)

  ! Check if all input files exist 
  IF (chkfile(cf_dfil) .OR. chkfile(cf_dfil2)) THEN
  	PRINT *, 'Missing density files'
  	STOP
  ELSE IF (chkfile(cf_ufil) .OR. chkfile(cf_vfil) .OR. chkfile(cf_vfil2) .OR. chkfile(cf_ufil2)) THEN
  	PRINT *, 'Missing velocity files'
  	STOP
  ELSE IF (chkfile(cf_bfil)) THEN
  	PRINT *, 'Missing bernoulli file'
  	STOP
  ENDIF

  ! Check for existence of grid definitions
  lchk = chkfile(cn_fhgr) 
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cf_dfil) .OR. lchk
  IF ( lchk ) THEN 
  	PRINT *,'Missing grid files ', TRIM(cn_fmsk),' ',TRIM(cn_fhgr),' or ',TRIM(cn_fzgr)
  	STOP  
  ENDIF                                                                         

  ! Read dimensions of the domain
  npiglo = getdim(cf_dfil, cn_x)
  npjglo = getdim(cf_dfil, cn_y)
  npk    = getdim(cf_dfil, cn_z)
  npt    = getdim(cf_dfil, cn_t)

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( e1u(npiglo,npjglo), e2v(npiglo,npjglo) )
  ALLOCATE ( ff(npiglo,npjglo) )
  ALLOCATE ( un2(npiglo,npjglo), sigx(npiglo,npjglo)  )
  ALLOCATE ( vn2(npiglo,npjglo), sigy(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo) , vn(npiglo,npjglo)  )
  ALLOCATE ( fmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( vort(npiglo,npjglo))
  ALLOCATE ( Jz1(npiglo,npjglo) , Jz2(npiglo,npjglo), Jz3(npiglo,npjglo) )
  ALLOCATE ( zberno(npiglo,npjglo), bx(npiglo,npjglo), by(npiglo,npjglo))
  ALLOCATE ( zsigsurf(npiglo,npjglo) , zsigsurf2(npiglo,npjglo) )
  ALLOCATE ( e3(npiglo,npjglo) )
  ALLOCATE ( tim(npt), tim2(npt) )
  ALLOCATE ( e1f(npiglo,npjglo), e2f(npiglo,npjglo) )

  ! Read the metrics from the mesh_hgr file
  e2u   = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)
  e1v   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
  e1u   = getvar(cn_fhgr, cn_ve1u,  1, npiglo, npjglo)
  e2v   = getvar(cn_fhgr, cn_ve2v,  1, npiglo, npjglo)
  ff    = getvar(cn_fhgr, cn_vff,   1, npiglo, npjglo)
  e1f   = getvar(cn_fhgr, cn_ve1f,  1, npiglo, npjglo)
  e2f   = getvar(cn_fhgr, cn_ve2f,  1, npiglo, npjglo)

  tim    = getvar1d(cf_dfil, cn_vtimec, npt     )
  tim2   = getvar1d(cf_dfil2, cn_vtimec, npt    )
  deltat = tim2(1) - tim(1)
  
  
  ipk(1)                        = 1
  stypvaru(1)%cname             = 'Jzbernoulli1'
  stypvaru(1)%cunits            = 'kg/m3/s2'
  stypvaru(1)%rmissing_value    = 999999.
  stypvaru(1)%valid_min         = -1e-4
  stypvaru(1)%valid_max         = 1e-4
  stypvaru(1)%clong_name        = 'Jz_rho_tendency'
  stypvaru(1)%cshort_name       = 'Jzbernoulli1'
  stypvaru(1)%conline_operation = 'N/A'
  stypvaru(1)%caxis             = 'TZYX'

  ipk(2)                        = 1
  stypvaru(2)%cname             = 'Jzbernoulli2'
  stypvaru(2)%cunits            = 'kg/m3/s2'
  stypvaru(2)%rmissing_value    = 999999.
  stypvaru(2)%valid_min         = -1e-4
  stypvaru(2)%valid_max         = 1e-4
  stypvaru(2)%clong_name        = 'Jz_vel_tendency'
  stypvaru(2)%cshort_name       = 'Jzbernoulli2'
  stypvaru(2)%conline_operation = 'N/A'
  stypvaru(2)%caxis             = 'TZYX'

  ipk(3)                        = 1
  stypvaru(3)%cname             = 'Jzbernoulli3'
  stypvaru(3)%cunits            = 'kg/m3/s2'
  stypvaru(3)%rmissing_value    = 999999.
  stypvaru(3)%valid_min         = -1e-4
  stypvaru(3)%valid_max         = 1e-4
  stypvaru(3)%clong_name        = 'Jz_bernoulli_cross'
  stypvaru(3)%cshort_name       = 'Jzbernoulli3'
  stypvaru(3)%conline_operation = 'N/A'
  stypvaru(3)%caxis             = 'TZYX'

  ipk(4)                        =  1
  stypvaru(4)%cname             = 'surf_rho'
  stypvaru(4)%cunits            = 'kg/m3/s2'
  stypvaru(4)%rmissing_value    = 999999.
  stypvaru(4)%valid_min         = 10
  stypvaru(4)%valid_max         = 40
  stypvaru(4)%clong_name        = 'Surface rho'
  stypvaru(4)%cshort_name       = 'surf_rho'
  stypvaru(4)%conline_operation = 'N/A'
  stypvaru(4)%caxis             = 'TZYX'

  ! create output filesets
  ncoutu = create (cf_out , cf_dfil , npiglo , npjglo , npk)
  ierr   = createvar (ncoutu , stypvaru , 4,      ipk , id_varoutu)
  ierr   = putheadervar (ncoutu , cf_dfil , npiglo , npjglo , npk)
  ierr   = putvar1d (ncoutu , tim , npt , 'T')

  ! Read bernoulli
  zberno = getvar(cf_bfil, 'voberno', 1, npiglo, npjglo)
  zberno = zberno*zohr0  
  
  ! Read Surface density at two times 
  zsigsurf(:,:)  = getvar(cf_dfil, cn_vosigma0, 1  , npiglo, npjglo)
  zsigsurf2(:,:) = getvar(cf_dfil2, cn_vosigma0, 1  , npiglo, npjglo)
  
  ! Read Surface velocity fields at two times
  un(:,:)  =  getvar(cf_ufil , cn_vozocrtx , 1 , npiglo , npjglo)
  vn(:,:)  =  getvar(cf_vfil , cn_vomecrty , 1 , npiglo , npjglo)
  un2(:,:) =  getvar(cf_ufil2 , cn_vozocrtx , 1 , npiglo , npjglo)
  vn2(:,:) =  getvar(cf_vfil2 , cn_vomecrty , 1 , npiglo , npjglo)


  tmask(:,:) = 1.d0
  WHERE( zsigsurf == 0 ) tmask = 0

 ! compute the mask 
  fmask(:,:) = 0.d0    
  DO jj = 1, npjglo - 1
    DO ji = 1, npiglo - 1
        fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
      IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
    ENDDO
  ENDDO 
    

  vort(:,:) = 0.d0
  DO jj = 1, npjglo -1 
    DO ji = 1, npiglo -1   ! vector opt. (e1 ~ dx, e2 ~ dy)
       vort(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                &         - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                &         * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
    END DO
  END DO


 DO jj=2,npjglo-1
   DO ji=2,npiglo-1

      bx(ji,jj) = ( zberno(ji+1,jj  ) - zberno(ji,jj) ) / e1u(ji,jj)*tmask(ji+1,jj)*tmask(ji,jj)
      by(ji,jj) = ( zberno(ji  ,jj+1) - zberno(ji,jj) ) / e2v(ji,jj)*tmask(ji,jj+1)*tmask(ji,jj)

      sigx(ji,jj) = ( zsigsurf(ji+1,jj  ) - zsigsurf(ji,jj) ) / e1u(ji,jj)*tmask(ji+1,jj)*tmask(ji,jj)
      sigy(ji,jj) = ( zsigsurf(ji  ,jj+1) - zsigsurf(ji,jj) ) / e2v(ji,jj)*tmask(ji,jj+1)*tmask(ji,jj)

   ENDDO
 ENDDO

 Jz1 = (ff + vort)*(zsigsurf2-zsigsurf)/deltat
 Jz2 = ((un2-un)/deltat)*sigy - ((vn2-vn)/deltat)*sigx
 Jz3 = (bx)*sigy - (by)*sigx

 ! fill in edges with 0 for now.
 Jz1(1     ,:)      = 0.0
 Jz1(:     ,1)      = 0.0
 Jz1(npiglo,:)      = 0.0
 Jz1(:     ,npjglo) = 0.0
 Jz2(1     ,:)      = 0.0
 Jz2(:     ,1)      = 0.0
 Jz2(npiglo,:)      = 0.0
 Jz2(:     ,npjglo) = 0.0
 Jz3(1     ,:)      = 0.0
 Jz3(:     ,1)      = 0.0
 Jz3(npiglo,:)      = 0.0
 Jz3(:     ,npjglo) = 0.0
 
 ! Write output variables 
 ierr = putvar(ncoutu, id_varoutu(1), Jz1(:,:), 1 , npiglo, npjglo)
 ierr = putvar(ncoutu, id_varoutu(2), Jz2(:,:), 1 , npiglo, npjglo)
 ierr = putvar(ncoutu, id_varoutu(3), Jz3(:,:), 1 , npiglo, npjglo)
 ierr = putvar(ncoutu, id_varoutu(4), zsigsurf(:,:),1 , npiglo, npjglo)
 
 ierr = closeout(ncoutu)

END PROGRAM cdfjz_all
