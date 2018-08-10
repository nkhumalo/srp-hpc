      logical function ana_rdfram(x,ix,w)
c
c $Id: ana_rdfram.F 19708 2010-10-29 18:04:21Z d3y133 $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3)
      integer ix(msa)
c
      character*108 card
      integer i,j,k,lq
      character*255 fname
c
      logical succes
      real*8 timp
      integer inwm,inwa
      integer flags(7)
      integer ipos, isize, ipos2
cDEBUG
      character*255 message
      character*63  string1,string2
      double precision util_clock_gettime
      double precision tautracegettime
      integer, save :: call_count
      data call_count/0/
cDEBUG
c
      lxw=.false.
      lvw=.false.
      lfw=.false.
      lsx=.false.
      lsv=.false.
      lfs=.false.
c
      card = " "
c
      if(me.eq.0) then
  100 continue
      if(fmttrj.eq.'trj') then
      inquire(unit=lfntrj,pos=ipos)
c     ipos2 = 0
c     inquire(unit=lfntrj,pos=ipos,size=isize)
c     do while (ipos+214+25*nsa.gt.isize)
c       read(unit=lfntrj,pos=ipos+ipos2,err=6001,end=6001,fmt='(a)')card
c6001   inquire(unit=lfntrj,size=isize)
c       read(unit=lfntrj,pos=ipos,err=6002,end=6002,fmt='(a)')card
c6002   continue
c       ipos2=ipos2+1
c       if (mod(ipos2,1000).eq.0) write(*,*)'HvD: size = ',isize
c     enddo
    1 continue
 6000 read(lfntrj,1000,pos=ipos,err=6000,end=6000) card
cDEBUG
c     write(*,*)'HvD: ana_rdfram: ipos,card=',ipos,trim(card)
cDEBUG
 1000 format(a5)
      if(card(1:5).ne.'frame') then
c       close(lfntrj)
c       call ana_rdhdr(byte_mb(i_snam))
        ipos = ipos + 1
        card = " "
        goto 1
      endif
      inquire(unit=lfntrj,pos=ipos)
 6010 read(lfntrj,1001,pos=ipos,err=6010,end=6010) timp,temp,pres,
     +     datum,tijd
 1001 format(2f12.6,e12.5,2a10)
      if(lclean.and.timp.lt.timr) then
      write(*,2221) timp,temp,pres,datum,tijd
 2221 format(' skip  ',2f12.6,e12.5,2a10)
      goto 1
      endif
      if(util_print('frames_read',print_debug))
     + write(*,2201) timp,temp,pres,datum,tijd
 2201 format(' frame ',2f12.6,e12.5,2a10)
      if(timp.lt.timr) timoff=time
      time=timp+timoff
      timr=timp
      inquire(unit=lfntrj,pos=ipos)
 6020 read(lfntrj,1100,pos=ipos,err=6020,end=6020) card
 1100 format(a36)
      if(card(54:54).eq.'.') then
      read(card,1012) box
 1012 format(f12.6,36x,f12.6,36x,f12.6)
      else
      read(card,1002) box(1)
 1002 format(f12.6)
      inquire(unit=lfntrj,pos=ipos)
 6030 read(lfntrj,1022,pos=ipos,err=6030,end=6030) box(2),box(3)
 1022 format(12x,f12.6,/,24x,f12.6)
      endif
      inquire(unit=lfntrj,pos=ipos)
 6040 read(lfntrj,1006,pos=ipos,err=6040,end=6040)
     +     lxw,lvw,lfw,lpw,lsx,lsv,lfs,lps,inwm,inwa,nsa
c6040 read(lfntrj,1110,pos=ipos,err=6040,end=6040) card
c1110 format(a38)
c     if(card(8:8).eq.'F'.or.card(8:8).eq.'T') then
c     read(card,1006) lxw,lvw,lfw,lpw,lsx,lsv,lfs,lps,inwm,inwa,nsa
c     elseif(card(6:6).eq.'F'.or.card(6:6).eq.'T') then
c     read(card,1007) lxw,lvw,lfw,lsx,lsv,lfs,inwm,inwa,nsa
c     lpw=.false.
c     lps=.false.
c     else
c     read(card,1003) lxw,lvw,lsx,lsv,inwm,inwa,nsa
c     lfw=.false.
c     lfs=.false.
c     endif
      if(inwm.gt.0) then
      nwm=inwm
      nwa=inwa
      endif
      if(nwm.gt.mwm)
     + call md_abort('Restart file not compatible with trajectory',nwm)
      if(util_print('frames_read',print_debug)) then
      write(*,'(a,f12.6,3i10,1x,4l1)')
     + ' Frame read:',time,inwm,inwa,nsa,lxw,lvw,lsx,lsv
      endif
 1003 format(4l1,3i10,2l1)
 1006 format(8l1,3i10)
 1007 format(6l1,3i10)
      if(nsa.gt.msa) call md_abort('Error in trj 1',0)
      if((lxw.or.lvw).and.nwm.gt.0) then
      do 2 i=1,nwm
      do 3 j=1,nwa
      inquire(unit=lfntrj,pos=ipos)
 6050 read(lfntrj,1004,pos=ipos,err=6050,end=6050) (w(i,j,k),k=1,3)
 1004 format(3f8.3)
    3 continue
    2 continue
      endif
      if(lsx.and.nsa.gt.0) then
      do 4 i=1,nsa
      inquire(unit=lfntrj,pos=ipos)
 6060 read(lfntrj,1005,pos=ipos,err=6060,end=6060) (x(i,j),j=1,3)
 1005 format(3f8.3)
      ix(i)=i
    4 continue
      endif
      succes=.true.
      goto 9000
      elseif(fmttrj.eq.'sco'.or.fmttrj.eq.'coo') then
cDEBUG
      write(*,*)'HvD: ana_rdfram: trouble 1'
      call md_abort('Trouble 1',0)
cDEBUG
      datum='          '
      tijd='          '
   11 continue
      inquire(unit=lfntrj,pos=ipos)
 6070 read(lfntrj,1000,pos=ipos,err=6070,end=6070) card
      if(card(1:5).ne.'time:') goto 11
      inquire(unit=lfntrj,pos=ipos)
 6080 read(lfntrj,2001,pos=ipos,err=6080,end=6080) timp
 2001 format(f12.6)
      if(timp.lt.timr) timoff=time
      time=timp+timoff
      timr=timp
   12 continue
      inquire(unit=lfntrj,pos=ipos)
 6090 read(lfntrj,1000,pos=ipos,err=6090,end=6090) card
      if(card(1:15).ne.'box_dimensions:') goto 12
      inquire(unit=lfntrj,pos=ipos)
 6100 read(lfntrj,2002,pos=ipos,err=6100,end=6100) box
 2002 format(3f12.6)
      lxw=fmttrj.eq.'coo'
      lvw=.false.
      lsx=.true.
      lsv=.false.
   13 continue
      inquire(unit=lfntrj,pos=ipos)
 6110 read(lfntrj,1000,pos=ipos,err=6110,end=6110) card
      if(card(1:10).ne.'atom_list:') goto 13
      do 14 i=1,nsa
      inquire(unit=lfntrj,pos=ipos)
 6120 read(lfntrj,1005,pos=ipos,err=6120,end=6120) (x(i,j),j=1,3)
      x(i,1)=x(i,1)*1.0d-1
      x(i,2)=x(i,2)*1.0d-1
      x(i,3)=x(i,3)*1.0d-1
      ix(i)=i
   14 continue
      if(lxw) then
      do 15 i=1,nwm
      do 16 j=1,nwa
      inquire(unit=lfntrj,pos=ipos)
 6130 read(lfntrj,1005,pos=ipos,err=6130,end=6130) (w(i,j,k),k=1,3)
   16 continue
   15 continue
      endif
      succes=.true.
      goto 9000
      else
      call md_abort('Trajectory file format error',0)
      endif
c
 9999 continue 
c
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      if(lrdf) call ana_rdfwrt()
c
      fname=filtrj
c
      lq=index(filtrj,'?')
      if(lq.eq.0) goto 9998
c
      iscof=iscof+1
      if(iscof.gt.ilast) goto 9998
      write(cnum,'(i3.3)') iscof
c
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9998,access="stream",form="formatted")
      write(*,3333) fname(1:index(fname,' ')-1)
 3333 format(/,' Opening trj file ',a)
c
      if(lrdf) call ana_rdfhdr(int_mb(i_rdf))
      goto 100
c
 9998 continue
      succes=.false.
c
      endif
 9000 continue
cDEBUG Simulate inter workflow communication for Wei
      if (succes) then
c     write(string1,'(e14.6)')util_clock_gettime()
      write(string1,'(e14.6)')tautracegettime()
      write(string2,'(f20.6)')timp
      if (call_count.eq.0) then
      message = '{"event-type":"input","name":"ANA_RDFRAM","time":"'//
     &trim(adjustl(string1))//'","node-id": "4", "thread-id": "0", '//
     &'"filename":"trajectory-'//trim(adjustl(string2))//'"}'
      else
      message = ',{"event-type":"input","name":"ANA_RDFRAM","time":"'//
     &trim(adjustl(string1))//'","node-id": "4", "thread-id": "0", '//
     &'"filename":"trajectory-'//trim(adjustl(string2))//'"}'
      endif
      write(80,*)trim(message)
      call_count = call_count + 1
      endif
cDEBUG
c
c     center and fold and rotate
c     --------------------------
c
      if(me.eq.0) then
      if(icent.gt.0) then
      if(.not.ltop) call md_abort('No topology for center',0)
      call ana_center_fold(x,w,int_mb(i_imol),nwm)
      endif
      if(nrot.gt.0) then
      call ana_rotate(x,w,nwm)
      endif
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d03,succes,ma_sizeof(mt_log,1,mt_byte),0)
      call ga_brdcst(mag_d04,ix,msa*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d05,x,3*msa*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d06,w,3*mwm*mwa*ma_sizeof(mt_dbl,1,mt_byte),0)
      do 234 i=1,6
      flags(i)=0
  234 continue
      if(lxw) flags(1)=1
      if(lvw) flags(2)=1
      if(lfw) flags(3)=1
      if(lsx) flags(4)=1
      if(lsv) flags(5)=1
      if(lfs) flags(6)=1
      flags(7)=nwm
      call ga_brdcst(mag_d07,flags,7*ma_sizeof(mt_int,1,mt_byte),0)
      lxw=flags(1).ne.0
      lvw=flags(2).ne.0
      lfw=flags(3).ne.0
      lsx=flags(4).ne.0
      lsv=flags(5).ne.0
      lfs=flags(6).ne.0
      nwm=flags(7)
      endif
c
      ana_rdfram=succes
c
      return
      end
      subroutine ana_center_fold(xs,xw,imol,nw)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer imol(msa),nw
c
      integer i,j,k,nmol,ncount
      real*8 xc(3),dx,factor
c
      real*8 boxh(3)
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      ncount=0
c
      nmol=0
      do 1 i=1,nsa
      nmol=max(nmol,imol(i))
      if(imol(i).ge.icent.and.imol(i).le.jcent) then
      ncount=ncount+1
      xc(1)=xc(1)+xs(i,1)
      xc(2)=xc(2)+xs(i,2)
      xc(3)=xc(3)+xs(i,3)
      endif
    1 continue
      if(ncount.eq.0) call md_abort('No atom in molecule',icent)
      xc(1)=xc(1)/dble(ncount)
      xc(2)=xc(2)/dble(ncount)
      xc(3)=xc(3)/dble(ncount)
c
      do 2 i=1,nsa
      xs(i,1)=xs(i,1)-xc(1)
      xs(i,2)=xs(i,2)-xc(2)
      xs(i,3)=xs(i,3)-xc(3)
    2 continue
      do 3 j=1,nwa
      do 4 i=1,nw
      xw(i,j,1)=xw(i,j,1)-xc(1)
      xw(i,j,2)=xw(i,j,2)-xc(2)
      xw(i,j,3)=xw(i,j,3)-xc(3)
    4 continue
    3 continue
c
      do 5 j=1,nmol
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      ncount=0
      do 6 i=1,nsa
      nmol=max(nmol,imol(i))
      if(imol(i).eq.j) then
      ncount=ncount+1
      xc(1)=xc(1)+xs(i,1)
      xc(2)=xc(2)+xs(i,2)
      xc(3)=xc(3)+xs(i,3)
      endif
    6 continue
      if(ncount.eq.0) goto 5
      xc(1)=xc(1)/dble(ncount)
      xc(2)=xc(2)/dble(ncount)
      xc(3)=xc(3)/dble(ncount)
      do 7 k=1,3
      dx=zero
      if(xc(k).lt.-boxh(k)) dx=box(k)
      if(xc(k).gt.boxh(k)) dx=-box(k)
      do 8 i=1,nsa
      if(imol(i).eq.j) xs(i,k)=xs(i,k)+dx
    8 continue
    7 continue
    5 continue
      if(nwa.gt.0.and.nw.gt.0) then
      factor=one/dble(nwa)
      do 9 i=1,nw
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      do 10 j=1,nwa
      xc(1)=xc(1)+xw(i,j,1)
      xc(2)=xc(2)+xw(i,j,2)
      xc(3)=xc(3)+xw(i,j,3)
   10 continue
      xc(1)=xc(1)*factor
      xc(2)=xc(2)*factor
      xc(3)=xc(3)*factor
      do 11 k=1,3
      dx=zero
      if(xc(k).lt.-boxh(k)) dx=box(k)
      if(xc(k).gt.boxh(k)) dx=-box(k)
      do 12 j=1,nwa
      xw(i,j,k)=xw(i,j,k)+dx
   12 continue
   11 continue
    9 continue
      endif
c
      return
      end
      subroutine ana_rotate(xs,xw,nw)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer nw
c
      real*8 rv(3),rw(3),rotang,rx(3),ry(3),pifac
      integer i,j,k
c
      pifac=tan(1.0d0)/90.0d0
c
      do 1 k=1,nrot
      rv(1)=0.0d0
      rv(2)=0.0d0
      rv(3)=0.0d0
      rw(1)=0.0d0
      rw(2)=0.0d0
      rw(3)=0.0d0
      rw(irot(k))=1.0d0
      rotang=arot(k)*pifac
      do 2 i=1,nsa
      rx(1)=xs(i,1)
      rx(2)=xs(i,2)
      rx(3)=xs(i,3)
      ry(1)=xs(i,1)
      ry(2)=xs(i,2)
      ry(3)=xs(i,3)
      call rotate(rv,rw,rotang,ry,rx)
      xs(i,1)=rx(1)
      xs(i,2)=rx(2)
      xs(i,3)=rx(3)
    2 continue
      do 3 i=1,nw
      do 4 j=1,nwa
      rx(1)=xw(i,j,1)
      rx(2)=xw(i,j,2)
      rx(3)=xw(i,j,3)
      ry(1)=xw(i,j,1)
      ry(2)=xw(i,j,2)
      ry(3)=xw(i,j,3)
      call rotate(rv,rw,rotang,ry,rx)
      xw(i,j,1)=rx(1)
      xw(i,j,2)=rx(2)
      xw(i,j,3)=rx(3)
    4 continue
    3 continue
    1 continue
c
      return
      end

