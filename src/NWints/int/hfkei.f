      Subroutine hfkei(alpha,E,Tab,Ti,Nint,NPP,La,Lb,Li,canAB)
c $Id: hfkei.f,v 1.2 1994-04-04 20:31:07 d3e129 Exp $

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Exponents

      Dimension alpha(2,NPP)

c--> Kinetic Energy Integrals

      Dimension Tab(Nint)

c--> Scratch Space

      Dimension Nxyz(3),Ti(NPP)
c
c Compute the kinetic energy integrals.
c
c     Formula:                                                       
c
c            1  /  Ia,Ib   Ja,Jb   Ka,Kb                       \
c     Tab =  - | T      Ey      Ez      + "Y-term" + "Z-term"   |
c            2  \  X       0       0                           /
c                                                                     
c      i,j         i-1,j-1      i-1,j+1       i+1,j-1       i+1,j+1   
c     T      = ijEx       - 2ibEx      - 2ajEx       + 4abEx          
c      X           0             0            0             0
c                                                                     
c******************************************************************************
  
c Initialize the block of KEIs.

      do 10 nn = 1,Nint
       Tab(nn) = 0.D0
   10 continue

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

c Loop over shell components.

      nn = 0

      do 420 ma = 1,La2

c Define the angular momentum indices for shell "A".

       call getNxyz(La,ma,Nxyz)

       Ia = Nxyz(1)
       Ja = Nxyz(2)
       Ka = Nxyz(3)

       if( canAB )then
        mb_limit = ma
       else
        mb_limit = Lb2
       end if

       do 410 mb = 1,mb_limit

c Define the angular momentum indices for shell "B".

        call getNxyz(Lb,mb,Nxyz)

        Ib = Nxyz(1)
        Jb = Nxyz(2)
        Kb = Nxyz(3)

        nn = nn + 1
  
c Build Tx.
  
        if( Ia.gt.0 .and. Ib.gt.0 )then
         do 100 m = 1,NPP
          Ti(m) =   0.5D0*(        Ia*Ib        )*E(1,m,0,Ia-1,Ib-1)
     &            -       (        Ia*alpha(2,m))*E(1,m,0,Ia-1,Ib+1)
     &            -       (alpha(1,m)*Ib        )*E(1,m,0,Ia+1,Ib-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
  100    continue
        else if( Ia.gt.0 )then
         do 110 m = 1,NPP
          Ti(m) = -       (        Ia*alpha(2,m))*E(1,m,0,Ia-1,Ib+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
  110    continue
        else if( Ib.gt.0 )then
         do 120 m = 1,NPP
          Ti(m) = -       (alpha(1,m)*Ib        )*E(1,m,0,Ia+1,Ib-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
  120    continue
        else
         do 130 m = 1,NPP
          Ti(m) =   2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
  130    continue
        end if
  
c Add Tx*Ey*Ez to Tab
  
        do 140 m = 1,NPP
         Tab(nn) = Tab(nn) + Ti(m)*E(2,m,0,Ja,Jb)*E(3,m,0,Ka,Kb)
  140   continue
  
c Build Ty.
  
        if( Ja.gt.0 .and. Jb.gt.0 )then
         do 200 m = 1,NPP
          Ti(m) =   0.5D0*(        Ja*Jb        )*E(2,m,0,Ja-1,Jb-1)
     &            -       (        Ja*alpha(2,m))*E(2,m,0,Ja-1,Jb+1)
     &            -       (alpha(1,m)*Jb        )*E(2,m,0,Ja+1,Jb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
  200    continue
        else if( Ja.gt.0 )then
         do 210 m = 1,NPP
          Ti(m) = -       (        Ja*alpha(2,m))*E(2,m,0,Ja-1,Jb+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
  210    continue
        else if( Jb.gt.0 )then
         do 220 m = 1,NPP
          Ti(m) = -       (alpha(1,m)*Jb        )*E(2,m,0,Ja+1,Jb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
  220    continue
        else
         do 230 m = 1,NPP
          Ti(m) =   2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
  230    continue
        end if
  
c Add Ex*Ty*Ez to Tab.
  
        do 240 m = 1,NPP
         Tab(nn) = Tab(nn) + E(1,m,0,Ia,Ib)*Ti(m)*E(3,m,0,Ka,Kb)
  240   continue
  
c Build Tz.
  
        if( Ka.gt.0 .and. Kb.gt.0 )then
         do 300 m = 1,NPP
          Ti(m) =   0.5D0*(        Ka*Kb        )*E(3,m,0,Ka-1,Kb-1)
     &            -       (        Ka*alpha(2,m))*E(3,m,0,Ka-1,Kb+1)
     &            -       (alpha(1,m)*Kb        )*E(3,m,0,Ka+1,Kb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
  300    continue
        else if( Ka.gt.0 )then
         do 310 m = 1,NPP
          Ti(m) = -       (        Ka*alpha(2,m))*E(3,m,0,Ka-1,Kb+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
  310    continue
        else if( Kb.gt.0 )then
         do 320 m = 1,NPP
          Ti(m) = -       (alpha(1,m)*Kb        )*E(3,m,0,Ka+1,Kb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
  320    continue
        else
         do 330 m = 1,NPP
          Ti(m) =   2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
  330    continue
        end if
  
c Add Ex*Ey*Tz to Tab.
  
        do 340 m = 1,NPP
         Tab(nn) = Tab(nn) + E(1,m,0,Ia,Ib)*E(2,m,0,Ja,Jb)*Ti(m)
  340   continue

  410  continue

  420 continue
  
      end
