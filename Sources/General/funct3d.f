      SUBROUTINE funct3d_par (lscreen, ier_flag)
      USE vmec_main
      USE vacmod, ONLY: bsqvac, bsqvac0, raxis_nestor, zaxis_nestor,
     &                  nuv, nuv3
      USE vmec_params, ONLY: ntmax, norm_term_flag
      USE realspace
      USE vforces
      USE xstuff
      USE timer_sub
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact, l_edge
      USE vparams, ONLY: twopi
      USE totzsp_mod
      USE tomnsp_mod
      USE timer_sub
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: SAXLASTNTYPE, ZEROLASTNTYPE,
     &                                SAXPBYLASTNTYPE
      USE blocktridiagonalsolver, ONLY: L_COLSCALE
      USE dbgout

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l0pi, l, ivacskip
      INTEGER :: nvskip0 = 0
      REAL(dp), DIMENSION(mnmax) ::
     &   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(dp), DIMENSION(:,:,:), POINTER :: lu, lv
      REAL(dp) :: presf_ns, delr_mse, delt0
      REAL(dp) :: tbroadon, tbroadoff
      REAL(dp), EXTERNAL :: pmass
      INTEGER :: i, j, k, nsmin, nsmax, m
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: bcastbuf
      INTEGER, DIMENSION(4) :: bbuf
      character(len=255) :: filename
C-----------------------------------------------
      CALL second0 (tfunon)

      print *, "funct3d"

!
!     POINTER ALIASES
!

      nfunct3d = nfunct3d + 1
      lu => pczmn;  lv => pcrmn

      ! dump state vector
      write(filename, 998) "dump_parvmec",                                     &
     &                     "funct3d_xc",                                      &
     &                     ns, nfunct3d
998   format(a,'/',a,'_',i5.5,'_',i5.5,'.dat')
      open(unit=73, file=trim(filename), status="unknown")
      write(73, *) pxc
      close(73)

!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
      ACTIVE1: IF (lactive) THEN
         IF (ictrl_prec2d .EQ. 3) THEN
            CALL SAXPBYLASTNTYPE(one, pxc, one, pxcdot, pgc)
            CALL SAXLASTNTYPE(pgc, pscalxc, pgc)
         ELSE IF (ictrl_prec2d.EQ.1 .AND. l_colscale) THEN
            pgc = (pxc-pxsave)*pcol_scale + pxsave
            CALL SAXLASTNTYPE(pgc, pscalxc, pgc)
         ELSE
            CALL SAXLASTNTYPE(pxc, pscalxc, pgc)
         END IF

!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!

         if (open_dbg_context("totzsp_input", num_eqsolve_retries)) then
           call add_real_5d("gc", 3, ntmax, ns, ntor1, mpol, pgc(:neqs),           &
     &                      order=(/ 4, 5, 3, 2, 1 /) )
           call close_dbg_out()
         end if

!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!

         CALL totzsps_par (pgc, pr1, pru, prv, pz1, pzu, pzv, lu, lv,
     &                     prcon, pzcon, ier_flag)

!
!     ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
!
         IF (lasym) THEN
            CALL totzspa_par (pgc, parmn, pbrmn, pextra3, pazmn, pbzmn,
     &                        pextra4, pblmn, pclmn, pextra1, pextra2)

!        SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!        TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
            CALL symrzl_par (pr1, pru, prv, pz1, pzu, pzv, lu, lv,
     &                       prcon, pzcon, parmn, pbrmn, pextra3, pazmn,
     &                       pbzmn, pextra4, pblmn, pclmn, pextra1,
     &                       pextra2)
         END IF

         if (open_dbg_context("funct3d_geometry", num_eqsolve_retries))        &
     &   then
             call add_real_4d("r1",   ns, 2, nzeta, ntheta3,   pr1,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("ru",   ns, 2, nzeta, ntheta3,   pru,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("rv",   ns, 2, nzeta, ntheta3,   prv,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("z1",   ns, 2, nzeta, ntheta3,   pz1,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("zu",   ns, 2, nzeta, ntheta3,   pzu,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("zv",   ns, 2, nzeta, ntheta3,   pzv,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("lu",   ns, 2, nzeta, ntheta3,    lu,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("lv",   ns, 2, nzeta, ntheta3,    lv,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("rcon", ns, 2, nzeta, ntheta3, prcon,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call add_real_4d("zcon", ns, 2, nzeta, ntheta3, pzcon,            &
     &                        order=(/ 3, 4, 1, 2 /) )
             call close_dbg_out()
         end if

!      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
!      router = r1(ns,0) + r1(ns,1)
!      rinner = r1(l0pi,0) + r1(l0pi,1)
         r00 = pr1(1,1,0)
         z00 = pz1(1,1,0)


!
!     COMPUTE CONSTRAINT RCON, ZCON
!
         nsmin=tlglob; nsmax=trglob
         DO l = nsmin, nsmax
            prcon(:,l,0) = prcon(:,l,0) + prcon(:,l,1)*psqrts(:,l)
            pzcon(:,l,0) = pzcon(:,l,0) + pzcon(:,l,1)*psqrts(:,l)
            pru0(:,l) = pru(:,l,0) + pru(:,l,1)*psqrts(:,l)
            pzu0(:,l) = pzu(:,l,0) + pzu(:,l,1)*psqrts(:,l)
         END DO

!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc. THIS
!     PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
!
!     NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
!     BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
!     SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
!
         IF (ictrl_prec2d .EQ. 2) THEN
            DO l = nsmin, nsmax
               prcon0(:,l) = prcon(:,l,0)
               pzcon0(:,l) = pzcon(:,l,0)
            END DO
         ELSE IF (iter2        .EQ. iter1 .AND.
     &            ivac         .LE. 0     .AND.
     &            ictrl_prec2d .EQ. 0) THEN
            print *, "extrapolate (r,z)con into volume"

#if defined(MPI_OPT)
            ALLOCATE(bcastbuf(2*nznt))
            bcastbuf(1:nznt)=prcon(:,ns,0)
            bcastbuf(nznt+1:2*nznt)=pzcon(:,ns,0)
            CALL second0(tbroadon)
            CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,nranks-1,
     &                     NS_COMM,MPI_ERR)
            CALL second0(tbroadoff)
            broadcast_time = broadcast_time + (tbroadoff-tbroadon)
            prcon(:,ns,0)=bcastbuf(1:nznt)
            pzcon(:,ns,0)=bcastbuf(nznt+1:2*nznt)
            DEALLOCATE(bcastbuf)
#endif
            DO l = nsmin, nsmax
               prcon0(:,l) = prcon(:,ns,0)*psqrts(:,l)**2
               pzcon0(:,l) = pzcon(:,ns,0)*psqrts(:,l)**2
            END DO
         END IF

!
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
         CALL jacobian_par

         IF (irst.EQ.2 .AND. iequi.EQ.0) THEN
           print *, "bad Jacobian -- early return from funct3d"
         end if

!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID

         CALL second0(tbcovon)
         CALL bcovar_par(lu, lv, pxc, ier_flag)
         CALL second0(tbcovoff)
         bcovar_time=bcovar_time+(tbcovoff - tbcovon)

      END IF ACTIVE1

      CALL MPI_BCast( ier_flag, 1, MPI_INTEGER, 0,
     &                RUNVMEC_COMM_WORLD, MPI_ERR) !SAL 070719

#if defined(MPI_OPT)
      bbuf(1)=irst; bbuf(2)=iequi; bbuf(3)=ivac; bbuf(4)=iter2
      CALL MPI_BCast(bbuf,4,MPI_INTEGER,0,RUNVMEC_COMM_WORLD,MPI_ERR)
      irst=bbuf(1); iequi=bbuf(2); ivac=bbuf(3); iter2=bbuf(4)
      CALL MPI_BCast(lfreeb,1,MPI_LOGICAL,0,RUNVMEC_COMM_WORLD,MPI_ERR)
      IF (ier_flag .ne. norm_term_flag) RETURN !SAL 070719
#endif

      IF (irst.EQ.2 .AND. iequi.EQ.0) THEN
         CALL ZEROLASTNTYPE(pgc)
         GOTO 100
      END IF

      timer(tbcov) = timer(tbcov) + (tbcovoff - tbcovon)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).

      IF (lfreeb       .AND.
     &    iter2 .GT. 1 .AND.
     &    iequi .EQ. 0) THEN

         IF (ictrl_prec2d.LE.1 .AND. (fsqr + fsqz).LE.1.e-3_dp) THEN
            print *, "activate Nestor"
            ivac = ivac+1   !decreased from e-1 to e-3 - sph12/04
         END IF

         IF (nvskip0 .EQ. 0) nvskip0 = MAX(1, nvacskip)

         IVAC0: IF (ivac .GE. 0) THEN
            print *, "compute Nestor contribution"
!SPH OFF: 6.20.17
!           IF INITIALLY ON, TURN OFF rcon0, zcon0 SLOWLY
            print *, "reduce (r,z)con"
            IF (lactive) THEN
               IF (ictrl_prec2d .EQ. 2) THEN
                  prcon0(:,nsmin:nsmax) = 0;  pzcon0(:,nsmin:nsmax) = 0
               ELSE IF (ictrl_prec2d .EQ. 0) THEN
                  prcon0(:,nsmin:nsmax) = 0.9_dp*prcon0(:,nsmin:nsmax)
                  pzcon0(:,nsmin:nsmax) = 0.9_dp*pzcon0(:,nsmin:nsmax)
               END IF
            ENDIF
            CALL second0 (tvacon)
            ivacskip = MOD(iter2-iter1,nvacskip)
            IF (ivac .LE. 2) THEN
               print *, "force full Nestor computation"
               ivacskip = 0
            END IF

!           EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
            IF (ivacskip .EQ. 0) THEN
               nvacskip = one/MAX(1.e-1_dp, 1.e11_dp*(fsqr+fsqz))
               nvacskip = MAX(nvacskip, nvskip0)
            END IF

!
!           NORMALLY, WHEN COMPUTING THE HESSIAN, IT IS SUFFICIENT TO
!           COMPUTE THE VARIATIONS IN THE "EXACT" SOLUTION, NOT THE ENTIRE
!           FIELD PERIOD SUM. THUS, FOR ictrl_prec2d >= 2, SET ivacskip = 1
!           FOR ictrl_prec2d = 1 (RUN WITH PRECONDITIONER APPLIED), MUST
!           COMPUTE EXACT VACUUM RESPONSE NOW.
!
!           THE EXCEPTION TO THIS IS IF WE ARE TESTING THE HESSIAN (lHess_exact=T),
!           THEN MUST USE FULL VACUUM CALCULATION TO COMPUTE IT (ivacskip=0)
!
!           lHess_exact = .FALSE.

            IF (ictrl_prec2d .NE. 0) THEN
               IF (lHess_exact .OR. ictrl_prec2d.EQ.2) THEN        !Accurate Hessian
                  ivacskip = 0
               ELSE
                  ivacskip = 1                                     !Fast vacuum calculation used to compute Hessian
               ENDIF
            ENDIF

!          NOTE: pgc contains correct edge values of r,z,l arrays
!                convert_sym, convert_asym have been applied to m=1 modes

            CALL convert_par(rmnc,zmns,lmns,rmns,zmnc,lmnc,pgc)

!          DO NOT UPDATE THIS WHEN USING PRECONDITIONER: BREAKS TRI-DIAGONAL STRUCTURE
            IF (ictrl_prec2d.EQ.0 .OR. ictrl_prec2d.EQ.2) THEN
               raxis_nestor(1:nzeta) = pr1(1:nzeta,1,0)
               zaxis_nestor(1:nzeta) = pz1(1:nzeta,1,0)

#if defined(MPI_OPT)
               ALLOCATE (bcastbuf(2*nzeta))
               bcastbuf(1:nzeta) = raxis_nestor(1:nzeta)
               bcastbuf(nzeta+1:2*nzeta) = zaxis_nestor(1:nzeta)
               CALL second0(tbroadon)
               CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     &                        RUNVMEC_COMM_WORLD,MPI_ERR)
               CALL second0(tbroadoff)
               broadcast_time = broadcast_time + (tbroadoff - tbroadon)
               raxis_nestor(1:nzeta) = bcastbuf(1:nzeta)
               zaxis_nestor(1:nzeta) = bcastbuf(nzeta+1:2*nzeta)
               DEALLOCATE (bcastbuf)
#endif
            END IF

#if defined(MPI_OPT)
            ALLOCATE (bcastbuf(2))
            bcastbuf(1)=rbtor
            bcastbuf(2)=ctor
            IF (lactive) THEN
               CALL second0(tbroadon)
               CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,
     &                        nranks-1,NS_COMM,MPI_ERR)
               CALL second0(tbroadoff)
               broadcast_time = broadcast_time + (tbroadoff -tbroadon)
            END IF

            CALL second0(tbroadon)
            IF (vlactive) THEN
               CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     &                        VAC_COMM,MPI_ERR)
            END IF
            CALL second0(tbroadoff)
            broadcast_time = broadcast_time + (tbroadoff -tbroadon)
            rbtor=bcastbuf(1)
            ctor=bcastbuf(2)
            DEALLOCATE (bcastbuf)
#endif

            IF (vlactive) THEN
               IF (ictrl_prec2d .NE. 3 .OR.
     &             l_edge) THEN
                  CALL vacuum_par (rmnc, rmns, zmns, zmnc, xm, xn,
     &                             ctor, rbtor, pwint_ns, ns, ivacskip,
     &                             ivac, mnmax, ier_flag, lscreen)
                  IF (ictrl_prec2d .EQ. 2) bsqvac0 = bsqvac
               ELSE
                  bsqvac = bsqvac0; ier_flag = 0
               END IF
            END IF

            IF (vnranks .LT. nranks) THEN
#if defined(MPI_OPT)
               CALL MPI_Bcast(bsqvac,SIZE(bsqvac),MPI_REAL8,0,
     &                        NS_COMM,MPI_ERR)
#endif
            END IF

            IF (ier_flag .NE. 0) THEN
               RETURN
            END IF
!
!          RESET FIRST TIME FOR SOFT START
!
            IF (ivac .EQ. 1) THEN
               print *, "first Nestor call - force restart"//
     &                  "via BAD_JACOBIAN logic"
               irst = 2
               delt0 = delt
               CALL restart_iter(delt0)
               irst = 1
            END IF

!
!          IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!          UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!          IF NON-VARIATIONAL FORCES ARE DESIRED
!
!          presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)
!          MUST NOT BREAK TRI-DIAGONAL RADIAL COUPLING: OFFENDS PRECONDITIONER!
            presf_ns = pmass(hs*(ns-1.5_dp))
            IF (presf_ns .NE. zero) THEN
               presf_ns = (pmass(1._dp)/presf_ns) * pres(ns)
            END IF

            print *, "update rbsq"
            DO l = 1, nznt
               bsqsav(l,3) = 1.5_dp*pbzmn_o(l,ns)
     &                     - 0.5_dp*pbzmn_o(l,ns-1)
               pgcon(l,ns) = bsqvac(l) + presf_ns
               rbsq(l) = pgcon(l,ns)*(pr1(l,ns,0) + pr1(l,ns,1))*ohs
               dbsq(l) = ABS(pgcon(l,ns)-bsqsav(l,3))
            END DO

            if (open_dbg_context("rbsq", num_eqsolve_retries)) then
              call add_real_2d("rbsq", nzeta, ntheta3, rbsq)
              call close_dbg_out()
            end if

            IF (ivac .EQ. 1) THEN
               IF (vlactive) THEN
                  bsqsav(:nznt,1) = pbzmn_o(:,ns)
                  bsqsav(:nznt,2) = bsqvac(:nznt)
#if defined(MPI_OPT)
                  CALL MPI_Bcast(bsqsav(:,1),nznt,MPI_REAL8,
     &                           nranks-1,NS_COMM,MPI_ERR)
#endif
               END IF
            ELSE IF (ictrl_prec2d .NE. 3) THEN
#if defined(MPI_OPT)
               CALL MPI_Bcast(bsqsav(:,1),nznt,MPI_REAL8,
     &                        0,NS_COMM,MPI_ERR)
#endif
            END IF

            CALL second0 (tvacoff)
            timer(tvac) = timer(tvac) + (tvacoff - tvacon)
            IF (ictrl_prec2d .GE. 2) THEN
               timer(tvac_2d) = timer(tvac_2d)+ (tvacoff - tvacon)
            END IF
         END IF IVAC0
      END IF
!
!     COMPUTE CONSTRAINT FORCE
!
      ACTIVE2: IF (lactive) THEN
         IF (iequi .NE. 1) THEN
            DO l = nsmin, nsmax
               pextra1(:,l,0) = (prcon(:,l,0) - prcon0(:,l))*pru0(:,l)
     &                        + (pzcon(:,l,0) - pzcon0(:,l))*pzu0(:,l)
            END DO
            CALL alias_par (pgcon, pextra1(:,:,0), pgc, pgc(1+mns),
     &                      pgc(1+2*mns), pextra1(:,:,1))
         ELSE
            IF (lrecon) THEN
               pxc(:ns) = pxc(:ns) + delr_mse
            END IF
            GOTO 100
         END IF

         if (open_dbg_context("constraint_force", num_eqsolve_retries))        &
     &   then

           call add_real_3d("rcon",  ns, nzeta, ntheta3, prcon(:,:,0),         &
     &                      order = (/ 2, 3, 1 /) )
           call add_real_3d("rcon0", ns, nzeta, ntheta3, prcon0,               &
     &                      order = (/ 2, 3, 1 /) )
           call add_real_3d("ru0",   ns, nzeta, ntheta3, pru0,                 &
     &                      order = (/ 2, 3, 1 /) )

           call add_real_3d("zcon",  ns, nzeta, ntheta3, pzcon(:,:,0),         &
     &                      order = (/ 2, 3, 1 /) )
           call add_real_3d("zcon0", ns, nzeta, ntheta3, pzcon0,               &
     &                      order = (/ 2, 3, 1 /) )
           call add_real_3d("zu0",   ns, nzeta, ntheta3, pzu0,                 &
     &                      order = (/ 2, 3, 1 /) )

           call add_real_3d("extra1", ns, nzeta, ntheta3,                      &
     &                      pextra1(:,:,0), order = (/ 2, 3, 1 /) )
           call add_real_3d("gcon",   ns, nzeta, ntheta3,                      &
     &                      pgcon,          order = (/ 2, 3, 1 /) )

           call add_real_3d("gcs", ns, ntor1, mpol,                            &
     &                      pgc(0*mns+1:1*mns), order = (/ 2, 3, 1 /) )
           call add_real_3d("gsc", ns, ntor1, mpol,                            &
     &                      pgc(1*mns+1:2*mns), order = (/ 2, 3, 1 /) )
           call add_real_3d("gcc", ns, ntor1, mpol,                            &
     &                      pgc(2*mns+1:3*mns), order = (/ 2, 3, 1 /) )
           call add_real_3d("gss", ns, ntor1, mpol,                            &
     &                      pextra1(:,:,1),     order = (/ 2, 3, 1 /) )

           call close_dbg_out()
         end if

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
         CALL forces_par

!     SYMMETRIZE FORCES (in u-v space)
!

         IF (lasym) THEN
            CALL symforce_par (parmn, pbrmn, pcrmn, pazmn, pbzmn,
     &                         pczmn, pblmn, pclmn, prcon, pzcon, pr1,
     &                         pru, prv, pz1, pzu, pzv, pextra3,
     &                         pextra4, pextra1, pextra2)
         END IF

!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!

         CALL tomnsps_par (pgc, parmn, pbrmn, pcrmn, pazmn, pbzmn,
     &                     pczmn, pblmn, pclmn, prcon, pzcon)

         IF (lasym) THEN
            CALL tomnspa_par (pgc, pr1, pru, prv, pz1, pzu, pzv,
     &                        pextra3, pextra4, pextra1, pextra2)
         END IF

!================================================================
!
!     COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
!
!================================================================
         CALL second0 (treson)

         CALL SAXLASTNTYPE(pgc, pscalxc, pgc)

         CALL residue_par(pgc, pgc(1+irzloff), pgc(1+2*irzloff))

      END IF ACTIVE2

!NEED THIS ON ALL PROCESSORS IN GROUP (NOT JUST ACTIVE ONES) FOR STOPPING CRITERION IN EVOLVE
#if defined(MPI_OPT)
      IF (gnranks .GT. nranks) THEN
         ALLOCATE(bcastbuf(6))
         bcastbuf(1) = fsqr; bcastbuf(2) = fsqr1
         bcastbuf(3) = fsqz; bcastbuf(4) = fsqz1
         bcastbuf(5) = fsql; bcastbuf(6) = fsql1
         CALL second0(tbroadon)
         CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     &                  RUNVMEC_COMM_WORLD,MPI_ERR)
         CALL second0(tbroadoff)
         broadcast_time = broadcast_time + (tbroadoff -tbroadon)
         fsqr = bcastbuf(1); fsqr1 = bcastbuf(2)
         fsqz = bcastbuf(3); fsqz1 = bcastbuf(4)
         fsql = bcastbuf(5); fsql1 = bcastbuf(6)
         DEALLOCATE(bcastbuf)
      END IF
#endif

!     Force new initial axis guess IF ALLOWED (l_moveaxis=T)
      IF (lmove_axis                   .and.
     &    iter2                .eq  .1 .and.
     &    (fsqr + fsqz + fsql) .gt. 1.E2_dp) THEN
         irst = 4
      END IF

      CALL second0 (tresoff)
      timer(tres) = timer(tres) + (tresoff - treson)

 100  CONTINUE

      CALL second0 (tfunoff)
      timer(tfun) = timer(tfun) + (tfunoff - tfunon)
      IF (ictrl_prec2d .GE. 2) THEN
         timer(tfun_2d) = timer(tfun_2d) + (tfunoff - tfunon)
      END IF

      END SUBROUTINE funct3d_par

      SUBROUTINE funct3d (lscreen, ier_flag)
      USE vmec_main
#ifdef _VACUUM2
      USE vac2_vacmod, ONLY: mf, nf, bsqvac
#else
      USE vacmod, ONLY: bsqvac, bsqvac0, raxis_nestor, zaxis_nestor
#endif
      USE vmec_params, ONLY: ntmax
      USE realspace
      USE vforces
      USE xstuff
      USE timer_sub
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact, l_edge
      USE vparams, ONLY: twopi
      USE totzsp_mod
      USE tomnsp_mod
      USE timer_sub
#ifdef _HBANGLE
      USE angle_constraints, ONLY: getrz, getfrho, xtempa
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(in) :: lscreen
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l0pi, l, lk, ivacskip
      INTEGER :: nvskip0 = 0
      REAL(dp), DIMENSION(mnmax) ::
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(dp), DIMENSION(:), POINTER :: lu, lv
      REAL(dp) :: presf_ns, delr_mse, delt0
      REAL(dp), EXTERNAL :: pmass
!-----------------------------------------------
!
!     POINTER ALIASES
!
      lu => czmn;  lv => crmn

      CALL second0 (tfunon)

!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
!
#ifdef _HBANGLE
!Overwrites rzl_array, but that is OK since gc = rzl_array in CALL, xc preserved
      xtempa = xc
      CALL getrz(xc)
#endif
      IF (ictrl_prec2d .EQ. 3) THEN
         gc(:neqs) = scalxc(:neqs)*(xc(:neqs)+xcdot(:neqs))
      ELSE
         gc(:neqs) = scalxc(:neqs)*xc(:neqs)
      END IF

!
!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!

      CALL totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)
!
!     ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
!
      IF (lasym) THEN
         CALL totzspa (gc, armn, brmn, extra3, azmn, bzmn, extra4,
     &                 blmn, clmn, extra1, extra2)

!        SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!        TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
         CALL symrzl (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon,
     &                armn, brmn, extra3, azmn, bzmn, extra4, blmn,
     &                clmn, extra1, extra2)
      END IF

      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
      router = r1(ns,0) + r1(ns,1)
      rinner = r1(l0pi,0) + r1(l0pi,1)
      r00 = r1(1,0)
      z00 = z1(1,0)

!
!     COMPUTE CONSTRAINT RCON, ZCON
!
#ifndef _HBANGLE
      rcon(:nrzt,0) = rcon(:nrzt,0) + rcon(:nrzt,1)*sqrts(:nrzt)
      zcon(:nrzt,0) = zcon(:nrzt,0) + zcon(:nrzt,1)*sqrts(:nrzt)
#endif
      ru0(:nrzt) = ru(:nrzt,0) + ru(:nrzt,1)*sqrts(:nrzt)
      zu0(:nrzt) = zu(:nrzt,0) + zu(:nrzt,1)*sqrts(:nrzt)

!
!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc. THIS
!     PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
!
!     NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
!     BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
!     SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
!
#ifndef _HBANGLE
      IF (ictrl_prec2d .EQ. 2) THEN
         rcon0(1:nrzt) = rcon(1:nrzt,0)
         zcon0(1:nrzt) = zcon(1:nrzt,0)
         ALLOCATE(bsqvac0(nznt))
      ELSE IF (iter2        .EQ. iter1 .AND.
     &         ivac         .le. 0     .AND.
     &         ictrl_prec2d .EQ. 0) THEN
         DO l = 1, ns
            rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
            zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
         END DO
      END IF
#endif
!
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
      CALL jacobian
      IF (irst  .EQ. 2 .AND.
     &    iequi .EQ. 0) THEN
         gc=0
         GOTO 100
      END IF

!
!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
!
      CALL second0 (tbcovon)
      CALL bcovar (lu, lv)

      CALL second0 (tbcovoff)
      timer(tbcov) = timer(tbcov) + (tbcovoff - tbcovon)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).

      IF (lfreeb .and. iter2.gt.1 .and. iequi.eq.0) THEN
         IF (ictrl_prec2d  .le. 1 .and.
     &       (fsqr + fsqz) .le. 1.e-3_dp) THEN
            ivac = ivac+1   !decreased from e-1 to e-3 - sph12/04
         END IF
         IF (nvskip0 .eq. 0) THEN
            nvskip0 = MAX(1, nvacskip)
         END IF
         IVAC0: IF (ivac .ge. 0) THEN
!SPH OFF: 6.20.17
!           IF INITIALLY ON, TURN OFF rcon0, zcon0 SLOWLY
            IF (ictrl_prec2d .eq. 2) THEN
               rcon0 = 0;  zcon0 = 0
            ELSE IF (ictrl_prec2d .eq. 0) THEN
               rcon0 = 0.9_dp*rcon0;  zcon0 = 0.9_dp*zcon0
            END IF
            CALL second0 (tvacon)
            ivacskip = MOD(iter2-iter1,nvacskip)
            IF (ivac .LE. 2) THEN
               ivacskip = 0
            END IF

!           EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
            IF (ivacskip .eq. 0) THEN
               nvacskip = one/MAX(1.e-1_dp, 1.e11_dp*(fsqr+fsqz))
               nvacskip = MAX(nvacskip, nvskip0)
            END IF

!
!           NORMALLY, WHEN COMPUTING THE HESSIAN, IT IS SUFFICIENT TO
!           COMPUTE THE VARIATIONS IN THE "EXACT" SOLUTION, NOT THE ENTIRE
!           FIELD PERIOD SUM. THUS, FOR ictrl_prec2d >= 2, SET ivacskip = 1
!           FOR ictrl_prec2d = 1 (RUN WITH PRECONDITIONER APPLIED), MUST
!           COMPUTE EXACT VACUUM RESPONSE NOW.
!
!           THE EXCEPTION TO THIS IS IF WE ARE TESTING THE HESSIAN (lHess_exact=T),
!           THEN MUST USE FULL VACUUM CALCULATION TO COMPUTE IT (ivacskip=0)
!
!           lHess_exact = .FALSE.
           IF (ictrl_prec2d .NE. 0) THEN
              IF (lHess_exact .OR. ictrl_prec2d.EQ.2) THEN        !Accurate Hessian
                 ivacskip = 0
              ELSE
                 ivacskip = 1                                     !Fast vacuum calculation used to compute Hessian
              END IF
           ENDIF

!           NOTE: gc contains correct edge values of r,z,l arrays
!                 convert_sym, convert_asym have been applied to m=1 modes
            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, gc, ns)
            IF (ictrl_prec2d.NE.3 .OR. l_edge) THEN
#ifdef _VACUUM2
               CALL vac2_vacuum(rmnc, rmns, zmns, zmnc, xm, xn, ctor,
     &                          ivacskip, mnmax)

#else
!           DO NOT UPDATE THIS WHEN USING PRECONDITIONER: BREAKS TRI-DIAGONAL STRUCTURE
               IF (ictrl_prec2d .EQ. 0 .OR.
     &             ictrl_prec2d .EQ. 2) THEN
                  raxis_nestor(1:nzeta) = r1(1:ns*nzeta:ns,0)
                  zaxis_nestor(1:nzeta) = z1(1:ns*nzeta:ns,0)
               END IF

               CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn,
     &                      ctor, rbtor, wint, ns, ivacskip, ivac,
     &                      mnmax, ier_flag, lscreen)
               IF (ictrl_prec2d .EQ. 2) THEN
                  bsqvac0 = bsqvac
               END IF
#endif
            ELSE
               bsqvac = bsqvac0
               ier_flag = 0
            END IF

            IF (ier_flag .NE. 0) THEN
               GOTO 100
            END IF
!
!           RESET FIRST TIME FOR SOFT START
!
            IF (ivac .eq. 1) THEN
               irst = 2;  delt0 = delt
#ifdef _HBANGLE
               gc = xc; xc = xtempa
#endif
               CALL restart_iter(delt0)
#ifdef _HBANGLE
               xc = gc
#endif
               irst = 1
            END IF

!
!           IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!           UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!           IF NON-VARIATIONAL FORCES ARE DESIRED
!
!           presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)
!           MUST NOT BREAK TRI-DIAGONAL RADIAL COUPLING: OFFENDS PRECONDITIONER!
            presf_ns = pmass(hs*(ns-1.5_dp))
            IF (presf_ns .ne. zero) THEN
               presf_ns = (pmass(one)/presf_ns) * pres(ns)
            END IF

            lk = 0
!            gcon(:nrzt) = r1(:nrzt,0)+sqrts(:nrzt)*r1(:nrzt,1)
!            gcon(1+nrzt) = 0
            DO l = ns, nrzt, ns
               lk = lk + 1
               bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)
#ifdef _ANIMEC
               gcon(l)     = bsqvac(lk) + pperp_ns(lk)
#else
               gcon(l)     = bsqvac(lk) + presf_ns
#endif

               rbsq(lk) = gcon(l)*(r1(l,0) + r1(l,1))*ohs
               dbsq(lk) = ABS(gcon(l)-bsqsav(lk,3))
            END DO
!
!           COMPUTE m=0,n=0 EDGE "pedestals"
!
!!            alphaR = hs*hs*ard(ns,1)
!!            IF (alphaR .ne. zero) alphaR =
!!     &         hs*SUM(wint(ns:nrzt:ns)*zu0(ns:nrzt:ns)*rbsq)/alphaR

!!            PRINT *,' alphaR/r1(ns) = ', alphaR/gcon(ns)

            IF (ivac .eq. 1) THEN
               bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns)
               bsqsav(:nznt,2) = bsqvac(:nznt)
            END IF
            CALL second0 (tvacoff)
            timer(tvac) = timer(tvac) + (tvacoff - tvacon)
            IF (ictrl_prec2d .GE. 2) THEN
               timer(tvac_2d) = timer(tvac_2d)+ (tvacoff - tvacon)
            END IF

         END IF IVAC0
      END IF

!
!     COMPUTE CONSTRAINT FORCE
!
#ifdef _HBANGLE
      gcon = 0
      IF (iequi .EQ. 1) THEN
         GOTO 100
      END IF
#else
      IF (iequi .NE. 1) THEN
         extra1(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt)
     &                   + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
         CALL alias (gcon, extra1(:,0), gc, gc(1+mns), gc(1+2*mns),
     &               extra1(:,1))
      ELSE
         IF (lrecon) xc(:ns) = xc(:ns) + delr_mse
         GOTO 100
      END IF
#endif
!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      CALL forces

!
!     FFT TEST
!
#ifdef _TEST_FOURIER
!      xc = xc*scalxc
!      CALL tomnsps_t (gc, xc, r1, ru, rv, z1, zu, zv)
!      IF (lasym) CALL tomnspa_t (gc, xc, r1, ru, rv, z1, zu, zv)
!      STOP
#endif
!
!     SYMMETRIZE FORCES (in u-v space): NOTE - gc IS SMALL BY FACTOR 2
!     IF lasym=T
!
      IF (lasym) THEN
         CALL symforce (armn, brmn, crmn, azmn, bzmn,
     &                  czmn, blmn, clmn, rcon, zcon, r1, ru, rv,
     &                  z1, zu, zv, extra3, extra4, extra1, extra2)
!NOT NECESSARY (EVEN THOUGH CORRECT)         gc = 2*gc
      END IF

!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!
      CALL tomnsps (gc, armn, brmn, crmn, azmn, bzmn, czmn,
     &              blmn, clmn, rcon, zcon)

      IF (lasym) THEN
         CALL tomnspa (gc, r1, ru, rv, z1, zu, zv,
     &                 extra3, extra4, extra1, extra2)
      END IF

!================================================================
!
!     COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
!
!================================================================
      CALL second0 (treson)

      gc = gc * scalxc    !!IS THIS CORRECT: SPH010214?
#ifdef _HBANGLE
      CALL getfrho(gc)
#endif
      CALL residue (gc, gc(1+irzloff), gc(1+2*irzloff))
!     Force new initial axis guess IF ALLOWED (l_moveaxis=T)
      IF (lmove_axis                  .and.
     &    iter2                .eq. 1 .and.
     &    (fsqr + fsqz + fsql) .gt. 1.E2_dp) THEN
         irst = 4
      END IF

      CALL second0 (tresoff)
      timer(tres) = timer(tres) + (tresoff - treson)

 100  CONTINUE

#ifdef _HBANGLE
      xc = xtempa
#endif
      CALL second0 (tfunoff)
      timer(tfun) = timer(tfun) + (tfunoff - tfunon)
      IF (ictrl_prec2d .GE. 2) THEN
         timer(tfun_2d) = timer(tfun_2d) + (tfunoff - tfunon)
      END IF

      END SUBROUTINE funct3d
