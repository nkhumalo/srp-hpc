*
* $Id: lattice.f,v 1.7 2003-07-12 22:12:31 bylaska Exp $
*

      real*8 function lattice_wcut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_wcut = wcut
      return
      end

      real*8 function lattice_ecut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_ecut = ecut
      return
      end

      real*8 function lattice_ggcut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_ggcut = 2.0d0*ecut
      return
      end

      real*8 function lattice_wggcut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_wggcut = 2.0d0*wcut
      return
      end



      real*8 function lattice_omega()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_omega = omega
      return
      end

      real*8 function lattice_unita(i,j)
      implicit none
      integer i,j

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_unita = unita(i,j)
      return
      end


      real*8 function lattice_unitg(i,j)
      implicit none
      integer i,j

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_unitg = unitg(i,j)
      return
      end



      subroutine lattice_init()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

*     **** local variables ****
      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8  gx,gy,gz,gg
      real*8  gg1,gg2,gg3
      real*8  ecut0,wcut0

*     **** external functions ****
      integer  control_ngrid
      real*8   control_unita,control_ecut,control_wcut
      external control_ngrid
      external control_unita,control_ecut,control_wcut
        
      ecut0 = control_ecut()
      wcut0 = control_wcut()
      unita(1,1) = control_unita(1,1)
      unita(2,1) = control_unita(2,1)
      unita(3,1) = control_unita(3,1)
      unita(1,2) = control_unita(1,2)
      unita(2,2) = control_unita(2,2)
      unita(3,2) = control_unita(3,2)
      unita(1,3) = control_unita(1,3)
      unita(2,3) = control_unita(2,3)
      unita(3,3) = control_unita(3,3)
      call get_cube(unita,unitg,omega)


*     *** set the ecut variable ***
c     call D3dB_nx(1,nx)
c     call D3dB_ny(1,ny)
c     call D3dB_nz(1,nz)
      nx = control_ngrid(1)
      ny = control_ngrid(2)
      nz = control_ngrid(3)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      gx = unitg(1,1)*dble(nxh)
      gy = unitg(2,1)*dble(nxh)
      gz = unitg(3,1)*dble(nxh)
      gg1 = gx*gx + gy*gy + gz*gz

      gx = unitg(1,2)*dble(nyh)
      gy = unitg(2,2)*dble(nyh)
      gz = unitg(3,2)*dble(nyh)
      gg2 = gx*gx + gy*gy + gz*gz

      gx = unitg(1,3)*dble(nzh)
      gy = unitg(2,3)*dble(nzh)
      gz = unitg(3,3)*dble(nzh)
      gg3 = gx*gx + gy*gy + gz*gz

      gg = gg1
      if (gg2.lt.gg) gg=gg2
      if (gg3.lt.gg) gg=gg3

      ecut = 0.5d0*gg
      if (ecut0.lt.ecut) then
         ecut = ecut0
      end if

      wcut = ecut
      if (wcut0.lt.wcut) then
         wcut = wcut0
      end if

      return
      end



*     *******************************
*     *                             *
*     *         lattice_r_grid      *
*     *                             *
*     *******************************
*
*     This routine computes coordinates of grid points in
*     the unit cell
*
*     Uses -
*          Parallel_taskid --- processor number
*          D3dB_nx --- number of grid points in direction 1
*          D3dB_ny --- number of grid points in direction 2
*          D3dB_nz --- number of grid points in direction 2
*          lattice_unita -- primitive lattice vectors in real space
*
*     Exit -
*          r  --- coordinates of grid points (Rx,Ry,Rz)
*
*
      subroutine lattice_r_grid(r)
      implicit none
      real*8 r(3,*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,q,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3)

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita


*     **** constants ****
      call Parallel_taskid(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/np1
         a(i,2) = lattice_unita(i,2)/np2
         a(i,3) = lattice_unita(i,3)/np3
      end do

      call dcopy(3*n2ft3d,0.0d0,0,r,1)

*     **** grid points in coordination space ****
      do k3 = -nph3, nph3-1
        do k2 = -nph2, nph2-1
          do k1 = -nph1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               call D3dB_ktoqp(1,k+1,q,p)
               if (p .eq. taskid) then
                  index = (q-1)*(np1+2)*np2
     >                  + j    *(np1+2)
     >                  + i+1
c               r(1,index) = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
c               r(2,index) = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
c               r(3,index) = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3

*               **** reverse y and z ****
                r(1,index) = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
                r(2,index) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
                r(3,index) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

               end if
          end do
        end do
      end do

      return
      end


      subroutine lattice_r_grid_sym(r)
      implicit none
      real*8 r(3,*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,q,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3)

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita


*     **** constants ****
      call Parallel_taskid(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/np1
         a(i,2) = lattice_unita(i,2)/np2
         a(i,3) = lattice_unita(i,3)/np3
      end do

      call dcopy(3*n2ft3d,0.0d0,0,r,1)

*     **** grid points in coordination space ****
      do k3 = -nph3+1, nph3-1
        do k2 = -nph2+1, nph2-1
          do k1 = -nph1+1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               call D3dB_ktoqp(1,k+1,q,p)
               if (p .eq. taskid) then
                  index = (q-1)*(np1+2)*np2
     >                  + j    *(np1+2)
     >                  + i+1
c               r(1,index) = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
c               r(2,index) = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
c               r(3,index) = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3

*               **** reverse y and z ****
                r(1,index) = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
                r(2,index) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
                r(3,index) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

               end if
          end do
        end do
      end do

      return
      end

      subroutine lattice_mask_sym(r)
      implicit none
      real*8 r(*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,q,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3


*     **** constants ****
      call Parallel_taskid(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2


      call dcopy(n2ft3d,0.0d0,0,r,1)

*     **** grid points in coordination space ****
      do k3 = -nph3+1, nph3-1
        do k2 = -nph2+1, nph2-1
          do k1 = -nph1+1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               call D3dB_ktoqp(1,k+1,q,p)
               if (p .eq. taskid) then
                  index = (q-1)*(np1+2)*np2
     >                  + j    *(np1+2)
     >                  + i+1

                r(index) =  1.0d0
               end if
          end do
        end do
      end do

      return
      end




      subroutine get_cube(unita,unitg,volume)

******************************************************************************
*                                                                            *
*     This routine computes primitive vectors both in coordination           *
*     space and in reciporocal space and the volume of primitive cell.       *
*                                                                            *
*     Inputs:                                                                *
*             type --- type of cube (1=SC, 2=FCC, 3=BCC, 4=linear)           *
*             unit --- lattice constants                                     *
*                                                                            *
*     Outputs:                                                               *
*             volume --- volume of primitive cell                            *
*             unita  --- primitive vectors in coordination space             *
*             unitg  --- primitive vectors in reciprocal space               *
*                                                                            *
*     Library:  DSCAL from BLAS                                              *
*                                                                            *
*     Last modification:  7/03/93  by R. Kawai                               *
*                                                                            *
******************************************************************************

      implicit none

*     ------------------
*     argument variables
*     ------------------
      double precision unita(3,3), unitg(3,3)
      double precision volume

*     ---------------
*     local variables
*     ---------------
      double precision twopi

      twopi = 8.0d0*datan(1.0d0)


*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      volume = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/volume,unitg,1)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      volume=dabs(volume)

      return
      end

*     *******************************
*     *                             *
*     *     lattice_abc_abg         *
*     *                             *
*     *******************************
*
*     This routine computes a,b,c,alpha,beta,gamma.
*
      subroutine lattice_abc_abg(a,b,c,alpha,beta,gamma)
      implicit none
      real*8 a,b,c
      real*8 alpha,beta,gamma

*     *** local variables ****
      real*8 d2,pi

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita

*     **** determine a,b,c,alpha,beta,gmma ***
      pi = 4.0d0*datan(1.0d0)
      a = dsqrt(lattice_unita(1,1)**2 
     >        + lattice_unita(2,1)**2 
     >        + lattice_unita(3,1)**2)
      b = dsqrt(lattice_unita(1,2)**2 
     >        + lattice_unita(2,2)**2 
     >        + lattice_unita(3,2)**2)
      c = dsqrt(lattice_unita(1,3)**2 
     >        + lattice_unita(2,3)**2 
     >        + lattice_unita(3,3)**2)
 
      d2 = (lattice_unita(1,2)-lattice_unita(1,3))**2 
     >   + (lattice_unita(2,2)-lattice_unita(2,3))**2 
     >   + (lattice_unita(3,2)-lattice_unita(3,3))**2
      alpha = (b*b + c*c - d2)/(2.0d0*b*c)
      alpha = dacos(alpha)*180.0d0/pi
 
      d2 = (lattice_unita(1,3)-lattice_unita(1,1))**2 
     >   + (lattice_unita(2,3)-lattice_unita(2,1))**2 
     >   + (lattice_unita(3,3)-lattice_unita(3,1))**2
      beta = (c*c + a*a - d2)/(2.0d0*c*a)
      beta = dacos(beta)*180.0d0/pi
 
      d2 = (lattice_unita(1,1)-lattice_unita(1,2))**2 
     >   + (lattice_unita(2,1)-lattice_unita(2,2))**2 
     >   + (lattice_unita(3,1)-lattice_unita(3,2))**2
      gamma = (a*a + b*b - d2)/(2.0d0*a*b)
      gamma = dacos(gamma)*180.0d0/pi

      return
      end
