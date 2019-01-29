!###################################################################### 
! Module to wrapp the computation of matrix H, G, D to python 
! Reference:
!  Numerical solution of acoustic scattering by finite perforated elastic plates
!  A. V. G. Cavalieri, W. R. Wolf and J. W. Jaworski - PRSA 2016
! More details about the formulation can be found in the reference
!###################################################################### 
module mntmat
  
  implicit none
  
contains
  
  subroutine hgmatrix(k0,x1,y1,x2,y2,xc,yc,n1,n2,ds,H,G, n)
    real(kind=8), intent(in) :: k0		
    real(kind=8),dimension(n), intent(in) :: x1 ,y1, x2 ,y2, xc, yc, n1, n2, ds
    integer(kind=4), intent(in) :: n
    complex(kind=8), dimension(n,n), intent(out) :: H, G
    !f2py depend(n)  x1 ,y1, x2 ,y2, xc, yc, n1, n2, ds, H, G
    integer(kind=4) i, j ,k, ng
    real(kind=8) xgauss(3), ygauss(3), wgauss(3)
    real(kind=8) dx, dy, arg,  auxx, auxy,  pi
    complex(kind=8) Ha0, Ha1, ii, axG, dGdn1, dGdn2 , argdx, argdy
    
    
    pi = dacos(-1.0d0)
    
    ! Define complex number
    ii = (0.0d0,1.0d0)
    
    ng = 3
    xgauss(1) = - dsqrt(3.0d0/5.0d0)
    xgauss(2) = 0.0d0
    xgauss(3) = dsqrt(3.0d0/5.0d0)
    ygauss(1) = - dsqrt(3.0d0/5.0d0)
    ygauss(2) = 0.0d0
    ygauss(3) =  dsqrt(3.0d0/5.0d0)
    wgauss(1) = 5.0d0/9.0d0
    wgauss(2) = 8.0d0/9.0d0
    wgauss(3) = 5.0d0/9.0d0
    
    H(:,:) = (0.0d0,0.0d0)
    G(:,:) = (0.0d0,0.0d0)
    do i = 1,n
       do j = 1,n
          if (i .ne. j) then
             axG = 0.0d0
             dGdn1 = 0.0d0
             dGdn2 = 0.0d0
             do k = 1,ng
                auxx = ((x2(j) - x1(j))/2.0d0) * xgauss(k) + ( x2(j) + x1(j) )/2.0d0
                auxy = ((y2(j) - y1(j))/2.0d0) * ygauss(k) + ( y2(j) + y1(j) )/2.0d0
                dx = xc(i) - auxx
                dy = yc(i) - auxy
                arg = k0 * dsqrt(dx**2.0d0 + dy**2.0d0)
                Ha0 = dcmplx( BESSEL_JN(0,arg), BESSEL_YN(0,arg))
                Ha1 = dcmplx( BESSEL_JN(1,arg), BESSEL_YN(1,arg))
                
                argdx = (ii / 4.0d0 ) * -Ha1 * ( - k0 * dx )/( dsqrt(dx**2.0d0 + dy**2.0d0))
                argdy = (ii / 4.0d0 ) * -Ha1 * ( - k0 * dy )/( dsqrt(dx**2.0d0 + dy**2.0d0))
                
                axG = axG + wgauss(k) * ( (ii/4.0d0) * Ha0 )
                dGdn1 = dGdn1 + wgauss(k) * argdx
                dGdn2 = dGdn2 + wgauss(k) * argdy
             enddo
             axG = axG * ds(j)/2.0d0
             dGdn1 = dGdn1 * ds(j)/2.0d0
             dGdn2 = dGdn2 * ds(j)/2.0d0
             
             H(i,j) = n1(j) * dGdn1 + n2(j) * dGdn2
             G(i,j) = axG
          else
             G(i,j) = ( ((1.0d0 - 0.5772156649d0 - dlog(0.25d0 * k0 * ds(i)))/(2.0 * pi )) + 0.25*ii ) * ds(i)
             H(i,j) = 0.5d0 
          endif
       end do
    end do
  end subroutine hgmatrix
  
  subroutine poroelastic(nplate,k0, alphaH, epsilon, omega, n2, ds, beta, phi,D,n,nm)
    integer(kind=4), intent(in) :: n,nm, nplate
    real(kind=8), intent(in) :: k0,alphaH,epsilon,omega
    real(kind=8),dimension(n), intent(in) :: n2,ds
    real(kind=8),dimension(nm), intent(in) :: beta
    real(kind=8),dimension(n,nm), intent(in) :: phi
    !f2py	integer(kind=4), intent(in), depend(phi) :: n=shape(phi,0), nm=shape(phi,1)	
    real(kind=8), dimension(n,n), intent(out) :: D
    integer(kind=4) i,j,k
    real(kind=8) pi,R,kr,aux1, aux2, auxE
    
    pi = dacos(-1.0d0)
    R = 1.0d-3
    kr = 4.0d0/pi
    
    aux1 = 0.0d0
    aux2 = 0.0d0
    D(:,:) = 0.0d0
    
    aux1 = ( 1.0d0 + alphaH*kr ) * epsilon * k0**5.0d0/omega**6.0d0 
    aux2 =  (k0**4.0d0)/(( 1.0d0 - alphaH )*omega**4.0d0 )
    
    do i = 1,nplate
       do k = 1,nplate 
          auxE = 0.0d0
          do j = 1,nm
             auxE = auxE + ( phi(k,j) * phi(i,j) ) / ( beta(j)**4.0d0 - aux2 )
          enddo
          D(i,k) =  aux1 * n2(i) * n2(k) * ds(k) * auxE
          
          if (k == i .or. k == nplate-i) then
             D(i,k) = D(i,k) - 0.5d0 * (alphaH / R) * kr * n2(i) * ( n2(k) )
          endif
       end do
    end do
    close(1)
  end subroutine poroelastic
  
  subroutine hgobs(k0,xc,yc,x1,y1,x2,y2,n1,n2,ds,Hobs,Gobs,nobs,n)
    integer(kind=4), intent(in) :: nobs, n
    real(kind=8), intent(in) :: k0
    real(kind=8),dimension(nobs), intent(in) :: xc,yc
    real(kind=8),dimension(n), intent(in) :: x1 ,y1, x2 ,y2, n1, n2, ds
    complex(kind=8), dimension(nobs,n), intent(out) :: Hobs, Gobs
    integer(kind=4) i, j ,k, ng
    real(kind=8) xgauss(3), ygauss(3), wgauss(3)
    real(kind=8) dx, dy, arg,  auxx, auxy,  pi
    complex(kind=8) Ha0, Ha1, ii, axG, dGdn1, dGdn2 , argdx, argdy
    
    
    pi = dacos(-1.0d0)
    
    ! Define complex number
    ii = (0.0d0,1.0d0)
    
    ng = 3
    xgauss(1) = - dsqrt(3.0d0/5.0d0)
    xgauss(2) = 0.0d0
    xgauss(3) = dsqrt(3.0d0/5.0d0)
    ygauss(1) = - dsqrt(3.0d0/5.0d0)
    ygauss(2) = 0.0d0
    ygauss(3) =  dsqrt(3.0d0/5.0d0)
    wgauss(1) = 5.0d0/9.0d0
    wgauss(2) = 8.0d0/9.0d0
    wgauss(3) = 5.0d0/9.0d0
    
    Hobs(:,:) = (0.0d0,0.0d0)
    Gobs(:,:) = (0.0d0,0.0d0)
    do i = 1,nobs
       do j = 1,n
          axG = 0.0d0
          dGdn1 = 0.0d0
          dGdn2 = 0.0d0
          do k = 1,ng
             auxx = ((x2(j) - x1(j))/2.0d0) * xgauss(k) + ( x2(j) + x1(j) )/2.0d0
             auxy = ((y2(j) - y1(j))/2.0d0) * ygauss(k) + ( y2(j) + y1(j) )/2.0d0
             dx = xc(i) - auxx
             dy = yc(i) - auxy
             arg = k0 * dsqrt(dx**2.0d0 + dy**2.0d0)
             Ha0 = dcmplx( BESSEL_JN(0,arg), BESSEL_YN(0,arg))
             Ha1 = dcmplx( BESSEL_JN(1,arg), BESSEL_YN(1,arg))
             
             argdx = (ii / 4.0d0 ) * -Ha1 * ( - k0 * dx )/( dsqrt(dx**2.0d0 + dy**2.0d0))
             argdy = (ii / 4.0d0 ) * -Ha1 * ( - k0 * dy )/( dsqrt(dx**2.0d0 + dy**2.0d0))
             
             axG = axG + wgauss(k) * ( (ii/4.0d0) * Ha0 )
             dGdn1 = dGdn1 + wgauss(k) * argdx
             dGdn2 = dGdn2 + wgauss(k) * argdy
          enddo
          axG = axG * ds(j)/2.0d0
          dGdn1 = dGdn1 * ds(j)/2.0d0
          dGdn2 = dGdn2 * ds(j)/2.0d0
          
          Hobs(i,j) = n1(j) * dGdn1 + n2(j) * dGdn2
          Gobs(i,j) = axG
          
       end do
    end do
  end subroutine hgobs
  
  subroutine field(k0,z1,z2,x1,y1,x2,y2,n1,n2,ds,dpdn,p,x,y,pfield,nx,ny,n)
    integer(kind=4), intent(in) :: nx,ny, n
    real(kind=8), intent(in) :: k0,z1,z2
    real(kind=8),dimension(n), intent(in) :: x1 ,y1, x2 ,y2, n1, n2, ds, dpdn, p 
    real(kind=8),dimension(nx,ny), intent(in) :: x,y
    !f2py	integer(kind=4), intent(in), depend(x) :: nx=shape(x,0), ny=shape(x,1)
    !f2py	integer(kind=4), intent(in), depend(y) :: nx=shape(y,0), ny=shape(y,1)	
    complex(kind=8), dimension(nx,ny), intent(out) :: pfield
    
    integer(kind=4) i, j, jj ,k, ng
    real(kind=8) xgauss(3), ygauss(3), wgauss(3)
    real(kind=8) dx, dy, arg,  auxx, auxy,  pi, xo, yo
    complex(kind=8) Ha0, Ha1, ii, argdx, argdy, Gf
    complex(kind=8), dimension(nx) :: Gobs
    complex(kind=8), dimension(nx,ny) :: Hobs
    
    
    pi = dacos(-1.0d0)
    
    ! Define complex number
    ii = (0.0d0,1.0d0)
    
    ng = 3
    xgauss(1) = - dsqrt(3.0d0/5.0d0)
    xgauss(2) = 0.0d0
    xgauss(3) = dsqrt(3.0d0/5.0d0)
    ygauss(1) = - dsqrt(3.0d0/5.0d0)
    ygauss(2) = 0.0d0
    ygauss(3) =  dsqrt(3.0d0/5.0d0)
    wgauss(1) = 5.0d0/9.0d0
    wgauss(2) = 8.0d0/9.0d0
    wgauss(3) = 5.0d0/9.0d0
    
    Hobs(:,:) = (0.0d0,0.0d0)
    Gobs(:) = (0.0d0,0.0d0)
    pfield(:,:) = (0.0d0,0.0d0)
    
    do jj = 1,ny
       do j = 1,nx
          xo = x(j,jj)
          yo = y(j,jj)
          ! Loop over the panels
          do i = 1,n
             
             Gobs(j) = (0.0d0,0.0d0)
             do k = 1,ng
                auxx = ((x2(i) - x1(i))/2.0d0) * xgauss(k) + ( x2(i) + x1(i) )/2.0d0
                auxy = ((y2(i) - y1(i))/2.0d0) * ygauss(k) + ( y2(i) + y1(i) )/2.0d0
                dx = xo - auxx
                dy = yo - auxy
                arg = k0 * dsqrt(dx**2.0d0 + dy**2.0d0)
                Ha0 = dcmplx( BESSEL_JN(0,arg), BESSEL_YN(0,arg))
                Ha1 = dcmplx( BESSEL_JN(1,arg), BESSEL_YN(1,arg))
                
                argdx = (ii / 4.0d0 ) * -Ha1 * ( - k0 * dx )/( dsqrt(dx**2.0d0 + dy**2.0d0))
                argdy = (ii / 4.0d0 ) * -Ha1 * ( - k0 * dy )/( dsqrt(dx**2.0d0 + dy**2.0d0))
                
                Gobs(j) = Gobs(j) + ds(i)*0.5d0*(wgauss(k) * ( (ii/4.0d0) * Ha0 ))*dpdn(i)
                
                Hobs(j,jj) = Hobs(j,jj) +  ds(i)*0.5d0*p(i)*wgauss(k)*(n1(i)*argdx + n2(i)*argdy)
                
             enddo
          end do
          dx = xo - z1
          dy = yo - z2
          arg = k0 * dsqrt(dx**2.0d0 + dy**2.0d0)
          Ha0 = dcmplx( BESSEL_JN(0,arg), BESSEL_YN(0,arg))
          Gf = (ii/4.0d0) * Ha0 
          
          pfield(j,jj) = Gobs(j) - Hobs(j,jj) - Gf	        
       end do
    enddo
  end subroutine field
  
  
end module mntmat
