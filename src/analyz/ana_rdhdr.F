      subroutine ana_rdhdr(sgmnam)
c
c $Id: ana_rdhdr.F 19708 2010-10-29 18:04:21Z d3y133 $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(msa)
      character*80 card
      character*255 fname
      integer i,lq
      integer ipos
c
      integer ib,jb,nsb,nwb
c
      if(me.eq.0) then
c
      iscof=ifrst
      timoff=0.0d0
      time=0.0d0
      timr=0.0d0
      write(cnum,'(i3.3)') iscof
c
      fname=filtrj
      lq=index(filtrj,'?')
      if(lq.gt.0) then
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
      endif
c
      if(fmttrj.eq.'trj') then
      write(*,3333) fname(1:index(fname,' ')-1)
 3333 format(/,' Opening trj file ',a)
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9999,access="stream",form="formatted")
      rewind(lfntrj)
      if(lrdf) call ana_rdfhdr(int_mb(i_rdf))
    1 continue
      inquire(unit=lfntrj,pos=ipos)
 6000 read(lfntrj,1000,pos=ipos,err=6000,end=6000) card
 1000 format(a6)
      if(card(1:6).ne.'header') goto 1
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      inquire(unit=lfntrj,pos=ipos)
 6010 read(lfntrj,1003,pos=ipos,err=6010,end=6010) nwa,nsa,nsb,nwb
      if (nsa.eq.0) goto 6010
cDEBUG
      write(*,*)'HvD: ana_rdhdr: nwa,nsa,nsb,nwb=',nwa,nsa,nsb,nwb
cDEBUG
      write(44,1003) nwa,nsa,nsb,nwb
 1003 format(4i10)
      if(nsa.gt.msa) call md_abort('Error in trj',0)
      if(nwa.gt.0) then
      do 2 i=1,nwa
      inquire(unit=lfntrj,pos=ipos)
 6020 read(lfntrj,1004,pos=ipos,err=6020,end=6020) wnam(i)
      write(44,1004) wnam(i)
 1004 format(a16)
    2 continue
      endif
      if(nsa.gt.0) then
      do 3 i=1,nsa
      inquire(unit=lfntrj,pos=ipos)
 6030 read(lfntrj,1005,pos=ipos,err=6030,end=6030) sgmnam(i)
      write(44,1005) sgmnam(i)
 1005 format(a16)
    3 continue
      endif
      if(nsb+nwb.gt.0) then
      do 4 i=1,nsb+nwb
      inquire(unit=lfntrj,pos=ipos)
 6040 read(lfntrj,1001,pos=ipos,err=6040,end=6040) ib,jb
      write(44,1001) ib,jb
 1001 format(2i8)
    4 continue
      endif
      close(unit=44)
      elseif(fmttrj.eq.'sco'.or.fmttrj.eq.'coo') then
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9999,access="stream",form="formatted")
      write(*,3333) fname(1:index(fname,' ')-1)
   11 continue
      inquire(unit=lfntrj,pos=ipos)
 6050 read(lfntrj,1000,pos=ipos,err=6050,end=6050) card
      if(card(1:18).eq.'num_solvent_atoms:')then
      read(card(19:26),'(i8)') nwa
      goto 11
      endif
      if(card(1:17).eq.'num_solute_atoms:') then
      read(card(18:25),'(i8)') nsa
      goto 11
      endif
      if(card(1:10).eq.'num_bonds:') then
      read(card(11:19),'(i9)') nsb
      else
      if(card(1:5).ne.'time:') goto 11
      if(.not.ltop) call md_abort('Format error trajectory file',0)
      endif
      rewind(lfntrj)
      nwa=0
      if(.not.ltop) then
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      write(44,1003) nwa,nsa,nsb,0
      endif
      if(nsa.gt.msa) call md_abort('Error in trj',0)
      if(fmttrj.ne.'coo') then
   12 continue
      inquire(unit=lfntrj,pos=ipos)
 6060 read(lfntrj,1000,pos=ipos,err=6060,end=6060) card
      if(card(1:11).ne.'atom_names:') goto 12
      do 13 i=1,nsa
      if(ltop) then
      inquire(unit=lfntrj,pos=ipos)
 6070 read(lfntrj,1000,pos=ipos,err=6070,end=6070) card
      else
      inquire(unit=lfntrj,pos=ipos)
 6080 read(lfntrj,1005,pos=ipos,err=6080,end=6080) sgmnam(i)
      write(44,1005) sgmnam(i)
      endif
   13 continue
      if(nsb.gt.0) then
   14 continue
      inquire(unit=lfntrj,pos=ipos)
 6090 read(lfntrj,1000,pos=ipos,err=6090,end=6090) card
      if(card(1:10).ne.'bond_list:') goto 14
      do 15 i=1,nsb
      inquire(unit=lfntrj,pos=ipos)
 6100 read(lfntrj,1001,pos=ipos,err=6100,end=6100) ib,jb
      if(.not.ltop) write(44,1001) ib,jb
   15 continue
      endif
      endif
      if(.not.ltop) close(unit=44)
      else
      call md_abort('Trajectory file format error',0)
      endif
c
      endif
c
      return
c
 9997 continue
      call md_abort('End of file reading trajectory file',0)
 9998 continue
      call md_abort('Error reading trajectory file',0)
 9999 continue
      call md_abort('Failed to open trajectory file',0)
      return
      end
