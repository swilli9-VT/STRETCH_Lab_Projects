c ======================================================================
c User Subroutine UMAT for Abaqus HGO material
c By William Snyder, PhD Student, Virginia Tech STRETCH LAB
c ======================================================================
      subroutine umat(
C Write only -
     1  stress,statev,ddsdde,sse,spd,scd,
C Read only -
     3 rpl,ddsddt,drplde,drpldt,stran,dstran,time,dtime,
     4 temp,dtemp,predef,dpred,cmname,ndi,nshr,ntens,nstatv,
     5 props,nprops,coords,drot,pnewdt,celent,F0,F1,
     6 noEl,npt,layer,kspt,jStep,kinc)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),F0(3,3),F1(3,3),
     4 jStep(4)
C
C Field Variable indices
      parameter( 
     $     i_E11 = 1,
     $     i_E22 = 2,
     $     i_E33 = 3,
     $     i_E12 = 4,     
     $     i_E13 = 5,     
     $     i_E23 = 6,  
     $     i_I1 = 7,
     $     i_I3 = 8,
     $     i_I4 = 9,
     $     i_I6 = 10,
     $     i_beta1 = 11,
     $	 i_beta2 = 12,	
     $     i_W = 13,
     $     i_Ei = 14,
     $     i_dam = 15,
     $     i_delete = 16, 
     $     n_sdv = 16)
C
C
C User defined material properties are stored as  
C     props(1) --> Shear Modulus Term, C10 
C     props(2) --> Bulk Modulus Term, D1 
C     props(3) --> Fiber Stiffness Term 1, k1 
C     props(4) --> Fiber Stiffness Term 2, k2
C     props(5) --> Dispersion Term, kappa
C     props(6) --> Number of Fiber Families, N
C     props(7-9) --> Fiber Family 1 Orientation Vector Components 
C     props(10-12) --> Fiber Family 2 Orientation Vector Components
C     props(13-14) --> E_alpha Damage Threshold Constants 
      parameter( 
     $     i_C10 = 1, 
     $     i_D1 = 2, 
     $     i_k1 = 3, 
     $     i_k2 = 4, 
     $     i_kappa = 5,
     $     i_N = 6, 
     $     i_a1x = 7,
     $     i_a1y = 8,
     $     i_a1z = 9,
     $     i_a2x = 10,
     $     i_a2y = 11,
     $     i_a2z = 12,
     $     i_T1 = 13,
     $     i_T2 = 14)
C
C
C  Declare Relevant Variables
      integer ff, aN, i, j, k, l, m, n, iter, insert
      real*8 C10, D1, xk1, xk2, xkappa, xJ, Ei_d, Ei_f, dam,
     $  xI1, EA, EAE, W, T1, T2, eps, xJp, Ei
C 
      dimension eye(3,3), R(3,3), RT(3,3), F(3,3), Bbar(3,3), Fbar(3,3),
     $  a_ref(3,2),a(3,1),aligned(3,3),dispersed(3,3),xI4(2),aoa(3,3),
     $  S_iso(3,3), S_vol(3,3), S_aniso(3,3), Cauchy(3,3), 
     $  indices(2,ntens), e_m(1,3), e_n(1,3), dF(3,3), F_p(3,3), 
     $  tau(3,3), tau_p(3,3), E(3,3),
     $  a_temp(3,1), b_temp(1,1), beta(2)  
      parameter(ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0, FOUR=4.d0,
     $ FIVE=5.d0, SIX=6.d0, SEVEN=7.d0, EIGHT = 8.d0, NINE=9.d0, 
     $ half=0.5d0)
      
      integer flag_i, flag_c
      real*8 tear_init, el_thickness, tear_len
      dimension tear_ic(2), tear_update(2), ndel_ic(2), ndel_update(2),
     $    ndel_curnt(2) 
      common /sharedf/ tear_init, el_thickness, tear_len, tear_ic,
     $  tear_update, ndel_ic, ndel_update, ndel_curnt, flag_i, flag_c
C      
      if ( time(2) .eq. ZERO ) then  
         if (nstatv .lt. n_sdv) then 
            call stdb_abqerr(-2,'Subroutine UMAT requires the '// 
     *           'specification of %I state variables. Check the '// 
     *           'definition of *DEPVAR in the input file.', 
     *           n_sdv,ZERO,' ') 
            call xit 
         end if 
      end if
C
      if ((noel .EQ. 1) .and. (npt .EQ. 1)) then
          tear_update = ZERO
          ndel_update = 0
      end if
C          
      ! HGO material properties and dispersion properties
      C10 = props(i_C10)
      D1 = props(i_D1)
      xk1 = props(i_k1)
      xk2 = props(i_k2) 
      xkappa = props(i_kappa) 
      aN = props(i_N)
C
      ! Fiber vectors in the reference configuration
      a_ref = reshape((/ props(i_a1x), props(i_a1y), props(i_a1z),
     $                   props(i_a2x), props(i_a2y), props(i_a2z)  /),
     $                   shape(a_ref))
C
      ! Set strain energy failure threshholds
      T1 = props(i_T1)
      T2 = props(i_T2)
C
      ! Initialize Identity Matrix
      eye = reshape((/ one,zero,zero,zero,one,zero,zero,zero,one /),
     $            shape(eye))
C
      ! Initialize elasticity tensor indices
      indices = reshape((/ 1, 1, 2, 2, 3, 3, 
     $                      1, 2, 1, 3, 2, 3 /), shape(indices))
C
      ! Initialize stress storage
      E = ZERO
      Ei = ZERO
      S_iso = ZERO
      S_vol = ZERO
      S_aniso = ZERO
      Cauchy = ZERO
      stress = ZERO
      ddsdde = ZERO
C        
C
      ! check if current integration point has been deleted 
      ! before doing field and state variable computations
      if (statev(i_delete) .NE. ZERO) then 
C
          ! in UMAT with a local coordinate system the deformation
          ! gradient, F, comes in the local form, F_local = R^T*F*R, 
          ! such that F = U*R rather than the typical F = R*U.
          ! We can use Newton iterations for polar decomposition to  
          ! get the rotational tensor, R, even in this atypical form.
C
C         Newton's Method Polar Decomposition
C         -----------------------------------
C
          R = F1 ! Begin with F as given on the first iteration
	    do i=1,20
		    ! Calculate the inverse determinant of the matrix
              detinv = 1/(R(1,1)*R(2,2)*R(3,3) - R(1,1)*R(2,3)*R(3,2)
     $                  - R(1,2)*R(2,1)*R(3,3) + R(1,2)*R(2,3)*R(3,1)
     $                  + R(1,3)*R(2,1)*R(3,2) - R(1,3)*R(2,2)*R(3,1))
C
              ! Calculate the inverse of the matrix
              RT(1,1) = +detinv * (R(2,2)*R(3,3) - R(2,3)*R(3,2))
              RT(2,1) = -detinv * (R(2,1)*R(3,3) - R(2,3)*R(3,1))
              RT(3,1) = +detinv * (R(2,1)*R(3,2) - R(2,2)*R(3,1))
              RT(1,2) = -detinv * (R(1,2)*R(3,3) - R(1,3)*R(3,2))
              RT(2,2) = +detinv * (R(1,1)*R(3,3) - R(1,3)*R(3,1))
              RT(3,2) = -detinv * (R(1,1)*R(3,2) - R(1,2)*R(3,1))
              RT(1,3) = +detinv * (R(1,2)*R(2,3) - R(1,3)*R(2,2))
              RT(2,3) = -detinv * (R(1,1)*R(2,3) - R(1,3)*R(2,1))
              RT(3,3) = +detinv * (R(1,1)*R(2,2) - R(1,2)*R(2,1))
C
              ! update value of R
              R = half*(R+transpose(RT))
          enddo
          !Apply correction to F
          F = matmul(F1, transpose(R))
C         
          !Record fiber orientation in the current configuration
		do ff=1,aN
			a_temp = reshape(a_ref(:,ff), shape(a_temp))
			b_temp = matmul(reshape((/zero, one, zero/), shape(e_m)), 
     $			matmul(F,a_temp))/norm2(matmul(F,a_temp))
			beta(ff) = b_temp(1,1)
		end do
          statev(i_beta1) = acosd(beta(1))
		statev(i_beta2) = acosd(beta(2))
C
          !Populate Green Strain Tensor
          E = half*(matmul(transpose(F), F) - eye)
          do i = 1,ntens 
              statev(i) = E(indices(1,i), indices(2,i))
          end do
C
C         End polar decomposition
C         -----------------------
C
C         
          ! Calculate determinant of the stretch tensor, xJ
          ! NOTE: det(F) = det(U*R) = det(U)*det(R) = det(U)*1 = det(U)
          xJ = F(1,1)*( F(2,2)*F(3,3) - F(2,3)*F(3,2) ) +
     $        F(1,2)*( F(2,3)*F(3,1) - F(2,1)*F(3,3) ) +
     $        F(1,3)*( F(3,2)*F(2,1) - F(3,1)*F(2,2) )
C
          !Calculate and store the third strain invariant, I3 = J^2
          statev(i_I3) = xJ**TWO
C
C       
C         Calculate isotropic strain energy and stress components
C         -------------------------------------------------------
C
          ! Calculate the deviatoric left Cauchy-Green Strain Tensor
          ! NOTE: B = F*F^T = (U*R)*(U*R)^T = U*R*R^T*U^T = U*U^T = U^2
          Bbar = MATMUL(F,TRANSPOSE(F)) / (xJ**(TWO/THREE))
C
          ! Calculate and record the first deviatoric strain invariant, I1_bar = tr(Bbar)
          xI1 = Bbar(1,1)+Bbar(2,2)+Bbar(3,3)
          statev(i_I1) = xI1
C
          ! Compute the isotropic component of the strain energy density
          W = C10*(xI1 - THREE)
C
          ! Compute the isotropic component of the Cauchy Stress
          S_iso = TWO*C10*(Bbar - (xI1/THREE)*eye)/xJ
C
C       
C         Calculate volumetric strain energy and stress components
C         --------------------------------------------------------
C
          ! Compute the volumetric component of the strain energy density
          W = W + (ONE/D1)*( (xJ**TWO - ONE)/TWO - log(xJ) )
C
          ! Compute the volumetric component of the Cauchy Stress
          S_vol = (ONE/D1)*(xJ - (ONE/xJ))*eye
C
C       
C         Calculate anisotropic strain energy and stress components
C         ---------------------------------------------------------
C
          if (xk1 .NE. 0.0) then ! skip these calculations if fibers have 0 stiffness
C
              ! calculate dispersed component of fiber stress
              dispersed = xkappa*(Bbar - (xI1/THREE)*eye)           
C              
              Fbar = F/(xJ**(ONE/THREE))
C
              do ff = 1, aN
                  aoa = ZERO
C                  
                  ! Put fiber orientations into the corotational current configuration
                  ! NOTE: a = U*a_ref = F*R^T*a_ref = U*R*R^T*a_ref
                  a = matmul(Fbar,reshape(a_ref(:,ff), shape(a)))
C                  
                  ! Calculate deviatoric fiber structural tensor: 
                  ! a \otimes a = a*a^T -> (a*a^T)*xJ^(-2/3) = abar*abar^T
                  aoa = MATMUL(a,transpose(a))
C          
                  ! Calculate the fiber pseudo-invariant, I4_bar = tr(abar*abar^T)
                  xI4(ff) = aoa(1,1) + aoa(2,2) + aoa(3,3) !- ONE
C
                  ! calculate aligned components of fiber stress
                  aligned = (ONE - THREE*xkappa)*
     $                      (aoa - (xI4(ff)/THREE)*eye)
C
                  !Calculate strain-like quantity, EA
                  EA = xkappa*(xI1 - THREE) + 
     $                (ONE - THREE*xkappa)*(xI4(ff) - ONE)
                  ! Apply Macauley Bracket to EA to ensure fibers do not contribute in compression
                  EA = HALF*(abs(EA) + EA)
                  ! Calculate exponential function of EA
                  EAE = EA*exp(xk2*EA**TWO)
C
                  ! Compute anisotropic strain energy density
                  W = W + (xk1/TWO*xk2)*(exp(xk2*EA**TWO) - ONE)
C
                  Ei = Ei + EA
                  ! Compute anisotropic stress component
                  S_aniso = S_aniso + TWO*xk1*EAE*
     $                     (dispersed + aligned)/xJ
              end do
C
	        statev(i_W) = W
              statev(i_Ei) = Ei
              statev(i_I4) = xI4(1)
              statev(i_I6) = xI4(2)
C
          end if
C     
          Ei_d = T1*(tear_len/tear_init)**(TWO/THREE)
          Ei_f = T2*(tear_len/tear_init)**(TWO/THREE)
C
	  if ((coords(3) .GT. 7.d0) .and. (flag_i .EQ. 1)) then
		Ei_d = 1.3d0*Ei_d
		Ei_f = 1.3d0*Ei_f
	  else if ((coords(3) .LE. 7.d0) .and. (flag_c .EQ. 1)) then
		Ei_d = 1.3d0*Ei_d
		Ei_f = 1.3d0*Ei_f
	  end if
C
          if (Ei .LT. Ei_d) then
                dam = ZERO
          else if (Ei .LE. Ei_f) then
                dam = ((Ei - Ei_d)/(Ei_f - Ei_d))**TWO
          else 
                dam = ONE
          end if
          
          if (dam .GT. statev(i_dam)) then
	        statev(i_dam) = dam
	    else
		    dam = statev(i_dam)
          end if
C		  
	  if (dam .LE. 0.99d0) then
            statev(i_delete) = 1
        else
            statev(i_delete) = 0
		  if (coords(3) .GT. 6.d0) then
                tear_update(1) = tear_update(1) + celent
                ndel_update(1) = ndel_update(1) + 1
		  else
                tear_update(2) = tear_update(2) + celent
                ndel_update(2) = ndel_update(2) + 1
		  end if
	  end if
C		
        Cauchy = (ONE - dam)*(S_iso + S_vol + S_aniso)
C
        do i = 1,ntens 
            stress(i) = Cauchy(indices(1,i), indices(2,i))
        end do
C
C
C         Calculate the elasticity tensor, DDSDDE
C         ---------------------------------------
C
        eps = 2.d-8 !perturbation
        tau = xJ*Cauchy !kirchhoff Stress
          
        do k = 1, ntens
              ! set up for perturbation on deformation gradient
	      m = indices(1,k)
	      n = indices(2,k)
	      e_m = reshape((/eye(m,1), eye(m,2), eye(m,3)/),
     $                        shape(e_m))
	      e_n = reshape((/eye(n,1), eye(n,2), eye(n,3)/),
     $                        shape(e_n))
C
              !Construct perturbation matrix
	      dF = half*eps*(matmul(transpose(e_m),matmul(e_n,F))
     $             + matmul(transpose(e_n),matmul(e_m,F)))
C
            ! Perturb F
            F_p = F + dF
            
            E = half*(matmul(transpose(F_p), F_p) - eye)
C
C
              !Calculate perturbed Jacobian
            xJp = F_p(1,1)*( F_p(2,2)*F_p(3,3) - F_p(2,3)*F_p(3,2) ) +
     $        F_p(1,2)*( F_p(2,3)*F_p(3,1) - F_p(2,1)*F_p(3,3) ) +
     $        F_p(1,3)*( F_p(3,2)*F_p(2,1) - F_p(3,1)*F_p(2,2) )
C              
            Bbar = MATMUL(F_p,TRANSPOSE(F_p)) / (xJp**(TWO/THREE))
C             
            xI1 = Bbar(1,1)+Bbar(2,2)+Bbar(3,3)
C             
            tau_p = TWO*C10*(Bbar - (xI1/THREE)*eye) 
     $                + (ONE/D1)*(xJp**TWO - ONE)*eye
	 
	      ! Reset strain-like quantity E_alpha for failure model update
	      Ei = ZERO
C          
            if (xk1.NE.0.0) then ! skip these calculations if fibers have 0 stiffness
C
                  ! calculate dispersed component of fiber stress
                  dispersed = xkappa*(Bbar - (xI1/THREE)*eye)         
C              
                  Fbar = F_p/(xJp**(ONE/THREE))
                  do ff = 1, aN
                      aoa = ZERO
                      
                      a = matmul(Fbar,reshape(a_ref(:,ff), shape(a)))
                      ! Calculate deviatoric fiber structural tensor: 
                      ! a \otimes a = a*a^T -> (a*a^T)*xJ^(-2/3) = abar*abar^T
                      aoa = MATMUL(a,transpose(a))
C          
                      ! Calculate the fiber pseudo-invariant, I4_bar = tr(abar*abar^T)
                      ! NOTE: I4 >= 1, enforced by I4 = 0.5*[|I4-1| + (I4-1)] + 1
                      xI4(ff) = aoa(1,1) + aoa(2,2) + aoa(3,3) !- ONE
                      !xI4(ff) = HALF*(abs(xI4(ff)) + xI4(ff)) + ONE
C
                      ! calculate aligned components of fiber stress
                      aligned = (ONE - THREE*xkappa)*
     $                      (aoa - (xI4(ff)/THREE)*eye)
C
                      !Calculate strain-like quantity, EA
                      EA = xkappa*(xI1 - THREE) + 
     $                (ONE - THREE*xkappa)*(xI4(ff) - ONE)
                      ! Apply Macauley Bracket to EA to ensure fibers do not contribute in compression
                      EA = HALF*(abs(EA) + EA)
                      ! Calculate exponential function of EA
                      EAE = EA*exp(xk2*EA**TWO)
C
                      ! Compute anisotropic stress component
                      tau_p = tau_p + TWO*xk1*EAE*(dispersed + aligned)
					  
		            ! Compute anisotropic strain energy density
		            Ei = Ei + EA
                  end do
              end if
C  
              if (Ei .LT. Ei_d) then
              	    dam = ZERO
              else if (Ei .LE. Ei_f) then
                    dam = ((Ei - Ei_d)/(Ei_f - Ei_d))**TWO
              else 
              	    dam = ONE
              end if
C              
              if (dam .LT. statev(i_dam)) then
                  dam = statev(i_dam)
              end if
C              
	        tau_p = (ONE - dam)*tau_p
C
              do l = 1,ntens
                  i = indices(1,l)
                  j = indices(2,l)
                  ddsdde(l,k) = one/(xJ*eps)*(tau_p(i,j)-tau(i,j))
              end do
        end do
C
      end if                      
C
      end subroutine umat
C      
C     
C      
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
C
      include 'aba_param.inc'
C
      parameter(ZERO=0.d0)
      integer i
      dimension time(2)
      
      integer flag_i, flag_c
      real*8 tear_init, el_thickness, tear_len
      dimension tear_ic(2), tear_update(2), ndel_ic(2), ndel_update(2),
     $    ndel_curnt(2) 
      common /sharedf/ tear_init, el_thickness, tear_len, tear_ic,
     $  tear_update, ndel_ic, ndel_update, ndel_curnt, flag_i, flag_c
C
      if ( lop .EQ. 0 ) then
C
          tear_init = 3.1d0
          el_thickness = 2.d0
          tear_len = tear_init
          tear_ic = ZERO
          ndel_ic = 0
	  flag_i = 0
	  flag_c = 0
          open(10, 
     $      file='D:\Users\Will\Tear_Propagation_Project\'//
     $      'tear_log.txt', status = 'replace', action='write')
          close(10)
C          
      else if (lop .EQ. 1) then
C
          tear_update = ZERO
          ndel_update = 0
C
      else if (lop .EQ. 2) then
C
          ndel_ic = ndel_ic + ndel_update
          ndel_curnt = ndel_curnt + ndel_update
          tear_ic = tear_ic + (tear_update/el_thickness)
C          
          open(3, 
     $      file='D:\Users\Will\Tear_Propagation_Project\'//
     $      'tear_log.txt', position='append', status = 'old', 
     $      action='write')
          write(3,*) 'Increment:', kinc, 'Step Time:', time(1)
          write(3,*) ndel_ic
          write(3,*) tear_ic
	  write(3,*) flag_i, flag_c
          write(3,*) ''
C          
C
          if (ndel_curnt(1) .GE. 2) then
                  tear_len = tear_init + tear_ic(1) + tear_ic(2)
                  ndel_curnt = 0
                  write(3,*) 'Updated tear length:'
                  write(3,fmt='(g15.8)') tear_len
                  write(3,*) ''
		  flag_i = 1
		  flag_c = 0
          else if (ndel_curnt(2) .GE. 2) then
                  tear_len = tear_init + tear_ic(1) + tear_ic(2)
                  ndel_curnt = 0
                  write(3,*) 'Updated tear length:'
                  write(3,fmt='(g15.8)') tear_len
                  write(3,*) ''
		  flag_i = 0
		  flag_c = 1
          end if
C
C          
          close(3)
C          
      else if (lop .EQ. 6) then
          open(6, 
     $      file='D:\Users\Will\Tear_Propagation_Project\'//
     $      'tear_log.txt', position='append', status = 'old', 
     $      action='write')
          write(6,*) 'Final tear length:'
          write(6,fmt='(g15.8)') tear_len
          close(6)
      end if
C      
      return
      end subroutine uexternaldb
          
