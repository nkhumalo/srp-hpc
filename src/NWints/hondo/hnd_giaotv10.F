      subroutine hnd_giaotv10(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,coord,zan,nat,nint,tv10,scr,lscr)
c
c $Id$
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr)
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension tv10(nint,3),coord(3,nat),zan(nat)
c
c Use scratch for temps in integral calculation
c Scratch needs are
c 6*(Li+2)*(Lj+3) + 6*(Li+2)*(Lj+1)*((Li+Lj+1)/2+1)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_giaotv101(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,coord,zan,nat,nint,tv10,scr(gh01(1)),
     2 scr(gh01(12)),scr(gh01(2)),scr(gh01(3)),scr(gh01(4)),
     3 scr(gh01(5)),scr(gh01(6)),scr(gh01(7)),scr(gh01(8)),scr(gh01(9)),
     4 scr(gh01(10)),scr(gh01(11)))
c
      return
      end
c
      subroutine hnd_giaotv101(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1                 expj,coefj,j_nprim,j_ngen,Lj,coord,zan,nat,nint,
     2                 tv10,xs,ys,zs,xxs,yys,zzs,xv,yv,zv,xxv,yyv,zzv)
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
      integer NoKinetic             ! FA-11-08-10
      common /skipKinetic/NoKinetic ! FA-11-08-10

      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,ccx,ccy,ccz
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension tv10(nint,3),coord(3,nat),zan(nat)
      dimension  xs(Li+2,Lj+3),ys(Li+2,Lj+3),zs(Li+2,Lj+3)
      dimension  xv(Li+2,Lj+1,*),yv(Li+2,Lj+1,*),zv(Li+2,Lj+1,*)
      dimension xxs(Li+2,Lj+3),yys(Li+2,Lj+3),zzs(Li+2,Lj+3)
      dimension xxv(Li+2,Lj+1,*),yyv(Li+2,Lj+1,*),zzv(Li+2,Lj+1,*)
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*3,0.0d0,0,tv10,1)
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      lit = Li + 1
      maxi = lit*(lit+1)/2
      litmod=lit+1
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      ljt = Lj + 1
      maxj = ljt*(ljt+1)/2
      ljtmod=ljt+2
c
      qijx=yi*zj-zi*yj
      qijy=zi*xj-xi*zj
      qijz=xi*yj-yi*xj
      tijx=xi-xj
      tijy=yi-yj
      tijz=zi-zj
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      nroots=(1+lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_giaotv10: need higher Rys rt',nroots,INT_ERR)
      endif
c
c     ----- i primitive -----
c
      do ig=1, i_nprim
      ai=expi(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=coefi(ig,i_ngen)
c
c     ----- j primitive -----
c
      do jg=1,j_nprim
      aj=expj(jg)
      aa=ai+aj
      aa1=1.0d0/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) goto 1000
      fac= exp(-dum)
      csj=coefj(jg,j_ngen)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      cij=csi*csj*fac
c
c     ----- kinetic energy -----
c
      if (NoKinetic .eq. 1) goto 20 ! SKIP KINETIC ENERGY-FA-11-08-10
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do j=1,ljtmod
         nj=j
         do i=1,litmod
            ni=i
            call hnd_sxyz
            xs(i,j)=xint*t
            ys(i,j)=yint*t
            zs(i,j)=zint*t
         enddo
         do i=1,lit
            xxs(i,j)=xs(i+1,j)
            yys(i,j)=ys(i+1,j)
            zzs(i,j)=zs(i+1,j)
         enddo
      enddo
c
      do i=1,lit
         xv(i,1,1)=(xs(i,1)-xs(i,3)*(aj+aj))*aj
         yv(i,1,1)=(ys(i,1)-ys(i,3)*(aj+aj))*aj
         zv(i,1,1)=(zs(i,1)-zs(i,3)*(aj+aj))*aj
         xxv(i,1,1)=(xxs(i,1)-xxs(i,3)*(aj+aj))*aj
         yyv(i,1,1)=(yys(i,1)-yys(i,3)*(aj+aj))*aj
         zzv(i,1,1)=(zzs(i,1)-zzs(i,3)*(aj+aj))*aj
         if (ljt.gt.1) then
            xv(i,2,1)=(xs(i,2)*dble(2+2-1)-xs(i,4)*(aj+aj))*aj
            yv(i,2,1)=(ys(i,2)*dble(2+2-1)-ys(i,4)*(aj+aj))*aj
            zv(i,2,1)=(zs(i,2)*dble(2+2-1)-zs(i,4)*(aj+aj))*aj
            xxv(i,2,1)=(xxs(i,2)*dble(2+2-1)-xxs(i,4)*(aj+aj))*aj
            yyv(i,2,1)=(yys(i,2)*dble(2+2-1)-yys(i,4)*(aj+aj))*aj
            zzv(i,2,1)=(zzs(i,2)*dble(2+2-1)-zzs(i,4)*(aj+aj))*aj
         endif
         do j=3,ljt
            xv(i,j,1)=(xs(i,j)*dble(j+j-1)-xs(i,j+2)*(aj+aj))*aj
     1                 -xs(i,j-2)*dble(((j-1)*(j-2))/2)
            yv(i,j,1)=(ys(i,j)*dble(j+j-1)-ys(i,j+2)*(aj+aj))*aj
     1                 -ys(i,j-2)*dble(((j-1)*(j-2))/2)
            zv(i,j,1)=(zs(i,j)*dble(j+j-1)-zs(i,j+2)*(aj+aj))*aj
     1                 -zs(i,j-2)*dble(((j-1)*(j-2))/2)
            xxv(i,j,1)=(xxs(i,j)*dble(j+j-1)-xxs(i,j+2)*(aj+aj))*aj
     1                 -xxs(i,j-2)*dble(((j-1)*(j-2))/2)
            yyv(i,j,1)=(yys(i,j)*dble(j+j-1)-yys(i,j+2)*(aj+aj))*aj
     1                 -yys(i,j-2)*dble(((j-1)*(j-2))/2)
            zzv(i,j,1)=(zzs(i,j)*dble(j+j-1)-zzs(i,j+2)*(aj+aj))*aj
     1                 -zzs(i,j-2)*dble(((j-1)*(j-2))/2)
         enddo
      enddo
c
      ij=0
      do j=1,maxj
         call getNxyz(Lj,j,Nxyz)
         jx = Nxyz(1) + 1
         jy = Nxyz(2) + 1
         jz = Nxyz(3) + 1
         do i=1,maxi
            call getNxyz(Li,i,Nxyz)
            ix = Nxyz(1) + 1
            iy = Nxyz(2) + 1
            iz = Nxyz(3) + 1
            dum = xv(ix,jx,1)* ys(iy,jy)  * zs(iz,jz)
     1          + xs(ix,jx)  * yv(iy,jy,1)* zs(iz,jz)
     2          + xs(ix,jx)  * ys(iy,jy)  * zv(iz,jz,1)
            dumx=xxv(ix,jx,1)* ys(iy,jy)  * zs(iz,jz)
     1          +xxs(ix,jx)  * yv(iy,jy,1)* zs(iz,jz)
     2          +xxs(ix,jx)  * ys(iy,jy)  * zv(iz,jz,1)
            dumy= xv(ix,jx,1)*yys(iy,jy)  * zs(iz,jz)
     1          + xs(ix,jx)  *yyv(iy,jy,1)* zs(iz,jz)
     2          + xs(ix,jx)  *yys(iy,jy)  * zv(iz,jz,1)
            dumz= xv(ix,jx,1)* ys(iy,jy)  *zzs(iz,jz)
     1          + xs(ix,jx)  * yv(iy,jy,1)*zzs(iz,jz)
     2          + xs(ix,jx)  * ys(iy,jy)  *zzv(iz,jz,1)
            ij=ij+1
            tv10(ij,1)=tv10(ij,1)+cij*(qijx*dum+tijy*dumz-tijz*dumy)
            tv10(ij,2)=tv10(ij,2)+cij*(qijy*dum+tijz*dumx-tijx*dumz)
            tv10(ij,3)=tv10(ij,3)+cij*(qijz*dum+tijx*dumy-tijy*dumx)
         enddo
      enddo
 20   continue ! SKIP KINETIC ENERGY-FA-11-08-10
c
c     ----- nuclear attraction -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do ic=1,nat
      znuc=-zan(ic)
      cx=coord(1,ic)
      cy=coord(2,ic)
      cz=coord(3,ic)
      yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      call hnd_droot
      do iroot=1,nroots
         uu=u9(iroot)*aa
         ww=w9(iroot)*znuc
         tt=1.0d0/(aa+uu)
         t = sqrt(tt)
         x0=(aax+uu*cx)*tt
         y0=(aay+uu*cy)*tt
         z0=(aaz+uu*cz)*tt
         do j=1,ljt
            nj=j
            do i=1,litmod
               ni=i
               call hnd_sxyz
               xv(i,j,iroot)=xint
               yv(i,j,iroot)=yint
               zv(i,j,iroot)=zint*ww
            enddo
            do i=1,lit
               xxv(i,j,iroot)=xv(i+1,j,iroot)
               yyv(i,j,iroot)=yv(i+1,j,iroot)
               zzv(i,j,iroot)=zv(i+1,j,iroot)
            enddo
         enddo
      enddo
c
      ij=0
      do j=1,maxj
         call getNxyz(Lj,j,Nxyz)
         jx = Nxyz(1) + 1
         jy = Nxyz(2) + 1
         jz = Nxyz(3) + 1
         do i=1,maxi
            call getNxyz(Li,i,Nxyz)
            ix = Nxyz(1) + 1
            iy = Nxyz(2) + 1
            iz = Nxyz(3) + 1
            dum =0.0d0
            dumx=0.0d0
            dumy=0.0d0
            dumz=0.0d0
            do iroot=1,nroots
            dum =dum + xv(ix,jx,iroot)* yv(iy,jy,iroot)* zv(iz,jz,iroot)
            dumx=dumx+xxv(ix,jx,iroot)* yv(iy,jy,iroot)* zv(iz,jz,iroot)
            dumy=dumy+ xv(ix,jx,iroot)*yyv(iy,jy,iroot)* zv(iz,jz,iroot)
            dumz=dumz+ xv(ix,jx,iroot)* yv(iy,jy,iroot)*zzv(iz,jz,iroot)
            enddo
            dum =dum *(aa1*pi212)
            dumx=dumx*(aa1*pi212)
            dumy=dumy*(aa1*pi212)
            dumz=dumz*(aa1*pi212)
            ij=ij+1
            tv10(ij,1)=tv10(ij,1)+cij*(qijx*dum+tijy*dumz-tijz*dumy)
            tv10(ij,2)=tv10(ij,2)+cij*(qijy*dum+tijz*dumx-tijx*dumz)
            tv10(ij,3)=tv10(ij,3)+cij*(qijz*dum+tijx*dumy-tijy*dumx)
         enddo
      enddo
c
c     ----- end loop over atoms -----
c
      enddo
c
c     ----- end loop over primitives -----
c
 1000 continue
      enddo
      enddo
c
      return
 9997 format(' in -gitv10- the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i2,/,' lit,ljt,nroots= ',3i3)
      end
