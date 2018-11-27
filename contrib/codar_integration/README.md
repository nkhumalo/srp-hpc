# ADIOS integration with NWChem and adaptive sampling

## A. Workflow of NWChem -> ADIOS -> adaptive sampling

* *test_staging.sh*: script for running the whole workflow with ADIOS staging (DATASPACES)
* *adios_read_ws_checker.py*: comparing the coordinates of ADIOS's output and the original NWChem MD output
* *md_compress_adaptive_adios_dataspaces.py*: adaptive sampling code with ADIOS reader
* *dataspaces.conf*: DATASPACES configuration file


## B. ADIOS integration with adaptive sampling

Code:
```
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank()
      size = comm.Get_size()

      ad.read_init("DATASPACES", comm, "verbose=3;")

      f = ad.file("nwchem_xyz.bp", "DATASPACES", comm, True, timeout_sec = 20.0)

      print("RANK:",rank)

      streaming_next = 0
      current_step = 0
```

Comments: MPI and ADIOS initialization

Code:
```
      v = f.vars['TWATOMS']
      twatoms = v.read(nsteps=1)
      v = f.vars['TSATOMS']
      tsatoms = v.read(nsteps=1)

      natoms=twatoms+tsatoms

      if rank == 0:
          mds2 = MDTrSampler(natoms, 2, conv_size = 10, n_samples=32, batch_size=50, manifold_size=25)
```

Comments: compute the total number of atoms and initalize the adaptive sampling algorithm

Code:
```
      v = f.vars['WID']
      Wmolecules = v.dims[0]
      v = f.vars['WLX']
      Watoms = v.dims[0]
      v = f.vars['SLX']
      Satoms = v.dims[0]

      localWatoms=int(Watoms/size)
      localWmolecules=int(Wmolecules/size)
      offsetWatoms=localWatoms*rank
      offsetWmolecules=localWmolecules*rank
      localSatoms=int(Satoms/size)
      offsetSatoms=localSatoms*rank
      if rank==(size-1):
          localWatoms=Watoms-localWatoms*(size-1)
          localWmolecules=Wmolecules-localWmolecules*(size-1)
          localSatoms=Satoms-localSatoms*(size-1)
```

Comments: compute the local number of water and solute atoms and offset

Code:
```
      comm.Barrier()
      if rank == 0:
          print(">>STEP:", current_step)

      if rank == 0:
          Wx = np.empty(Watoms, dtype=np.float64)
          Wy = np.empty(Watoms, dtype=np.float64)
          Wz = np.empty(Watoms, dtype=np.float64)
          Wi = np.empty(Wmolecules, dtype=np.int32)
          Sx = np.empty(Satoms, dtype=np.float64)
          Sy = np.empty(Satoms, dtype=np.float64)
          Sz = np.empty(Satoms, dtype=np.float64)
          Si = np.empty(Satoms, dtype=np.int32)
          ordered_x = np.empty(Watoms+Satoms, dtype=np.float64)
          ordered_y = np.empty(Watoms+Satoms, dtype=np.float64)
          ordered_z = np.empty(Watoms+Satoms, dtype=np.float64)
      else:
          Wx = None
          Wy = None
          Wz = None
          Wi = None
          Sx = None
          Sy = None
          Sz = None
          Si = None
          ordered_x = None
          ordered_y = None
          ordered_z = None

      comm.Barrier()
      vx = f.var['WLX']
      vy = f.var['WLY']
      vz = f.var['WLZ']
      vi = f.var['WID']
      x_val = vx.read(offset=(offsetWatoms,),count=(localWatoms,),nsteps=1)
      y_val = vy.read(offset=(offsetWatoms,),count=(localWatoms,),nsteps=1)
      z_val = vz.read(offset=(offsetWatoms,),count=(localWatoms,),nsteps=1)
      i_val = vi.read(offset=(offsetWmolecules,),count=(localWmolecules,),nsteps=1)
      sendWatoms = np.array(comm.gather(localWatoms, 0))
      comm.Gatherv(x_val, (Wx, sendWatoms), 0)
      comm.Gatherv(y_val, (Wy, sendWatoms), 0)
      comm.Gatherv(z_val, (Wz, sendWatoms), 0)
      sendWmolecules = np.array(comm.gather(localWmolecules, 0))
      comm.Gatherv(i_val, (Wi, sendWmolecules), 0)

      comm.Barrier()
      vx = f.var['SLX']
      vy = f.var['SLY']
      vz = f.var['SLZ']
      vi = f.var['SID']
      x_val = vx.read(offset=(offsetSatoms,),count=(localSatoms,),nsteps=1)
      y_val = vy.read(offset=(offsetSatoms,),count=(localSatoms,),nsteps=1)
      z_val = vz.read(offset=(offsetSatoms,),count=(localSatoms,),nsteps=1)
      i_val = vi.read(offset=(offsetSatoms,),count=(localSatoms,),nsteps=1)

      sendSatoms = np.array(comm.gather(localSatoms, 0))
      comm.Gatherv(x_val, (Sx, sendSatoms), 0)
      comm.Gatherv(y_val, (Sy, sendSatoms), 0)
      comm.Gatherv(z_val, (Sz, sendSatoms), 0)
      comm.Gatherv(i_val, (Si, sendSatoms), 0)
```

Comments: use ADIOS to read the water and solute atoms locally and then gather them into MPI rank 0

Code:
```
      if rank == 0:
        nwa = int(Watoms/Wmolecules)
        k = 0
        for i in range(Wmolecules):
            des = (Wi[i]-1)*nwa
            for j in range(nwa):
                ordered_x[des] = Wx[k]
                ordered_y[des] = Wy[k]
                ordered_z[des] = Wz[k]
                des = des + 1
                k = k + 1
       k = 0
       for i in range(Satoms):
            des = Watoms - 1 + Si[i]
            ordered_x[des] = Sx[k]
            ordered_y[des] = Sy[k]
            ordered_z[des] = Sz[k]
            k = k + 1

       master_xyz = np.array([ordered_x,ordered_y,ordered_z]).T
       mds2.adaptive_sampling_step(master_xyz[:,:]) 
```

Comments: order the water and solute atoms and send the ordered coordinates to adaptive sampling

## C. ADIOS integration with NWChem

The integration code can be found in [https://github.com/hjjvandam/nwchem-1](https://github.com/hjjvandam/nwchem-1) on the branch "pretauadio".
 
This ADIOS integration with NWChem can work with either MPI mode or DATASPACES mode.
The MPI mode outputs the trajectory file with the ADIOS BP format. It can be queried
by using bpls command (e.g., `bpls -l xx.bp`, `bpls` is installed from ADIOS). The
DATASPACES mode outputs the data into a dataspaces server and then connect to the
analysis code. There will be no file generated during the execution. To set the
ADIOS mode in NWChem, you can use
```
  set "sp:adios_mode" "DATASPACES"
```
in the `.nw` file.

The following sections present the details of how ADIOS integrates with NWChem.

### 1. nwchem-1/nwchem_make_env.sh

Code: export BLAS_LIB="-lsci_gnu_71_mpi -lsci_gnu_71 `adios_config -l -f`"

Comments: ADIOS compilation flags

### 2. nwchem-1/src/config/makefile.h

Code: INCLUDES = -I. $(LIB_INCLUDES) -I$(INCDIR) $(INCPATH) $(shell adios_config -c -f)

Comments: ADIOS compilation flags

### 3. nwchem-1/src/nwmd/md_main.F

Code:
```
      call sp_wrttrj(lfntrj,lxw,.false.,.false.,lxs,.false.,.false.,
      + stime,pres,temp,tempw,temps,
      + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
      + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
      + dbl_mb(i_fs),mdstep,mdacq)
```

Comments: add the code of passing the timestep (mdstep) to sp_init which calls the main ADIOS integration code

### 4. nwchem-1/src/space/sp_init.F

Code:
```
      allocate(nwmnlist(np))
      allocate(nsanlist(np))
      allocate(templist(np))
c
      call ga_brdcst(10,numw,ma_sizeof(MT_INT,1,MT_BYTE),0)
      call ga_brdcst(10,nums,ma_sizeof(MT_INT,1,MT_BYTE),0)
c
      call ga_sync
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_access(ga_ip, ilp, ihp, jlp, jhp, localindex, localld)
      localipl(1:mbox,1:mip2)=>int_mb(localindex:localindex+mbox*mip2)
      nwmn=localipl(1,2)
      nsan=ipl(2,2)
```

Comments: allocate memory, broadcast the number of total water and solute atoms to different MPI ranks and give the local number of water and solute atoms to nwmn and nsan

Code:
```
      do i = 1, np
        templist(i)=0
      enddo
      templist(me+1)=nwmn
      call ga_sync
      call ga_igop(10,templist,np,'+')
      do i = 1, np
        nwmnlist(i)=templist(i)
      enddo
c
      do i = 1, np
        templist(i)=0
      enddo
      templist(me+1)=nsan
      call ga_sync
      call ga_igop(10,templist,np,'+')
      do i = 1, np
        nsanlist(i)=templist(i)
      enddo
c
      do i = 2, np
        nwmnlist(i)=nwmnlist(i)+nwmnlist(i-1)
        nsanlist(i)=nsanlist(i)+nsanlist(i-1)
      enddo
```

Comments: based on the local number of water and solute atoms, all MPI ranks compute its offset in the global view. For example, 4 MPI ranks have water atoms of 1, 2, 3, 4, then their offset would be 0, 1, 3, 6

Code:
```
      if(nwmn.gt.0) then
c
      call ga_distribution(ga_iw,me,ili,ihi,jli,jhi)
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      call sp_unpackw(nwmn,iwl,iwlp)
      allocate(localiwl(nwmn))
      do i = 1, nwmn
c      localiwl(i)=iwl(i,lwgmn)-number
      localiwl(i)=iwl(i,lwgmn)
      enddo
c
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_access(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,localindex,
     +localld)
      localxw(1:mwm,1:3,1:mwa)=>dbl_mb(localindex:localindex+mwm*3*mwa)
      allocate(localWx(nwa*nwmn))
      allocate(localWy(nwa*nwmn))
      allocate(localWz(nwa*nwmn))
      j = 1
      do i = 1, nwmn
        do k = 1, nwa
          localWx(j) = localxw(i,1,k)
          localWy(j) = localxw(i,2,k)
          localWz(j) = localxw(i,3,k)
          j = j + 1
        enddo
      enddo
c
      endif
```

Comments: load the actual coordinates of water atoms to local allocated memory from the GA space

code:
```
      if(nsan.gt.0) then
c 
      call ga_distribution(ga_is,me,ili,ihi,jli,jhi)
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call sp_unpack(nsan,isl,islp)
      allocate(localisl(nsan))
      do i = 1, nsan
c      localisl(i)=isl(i,lsgan)-number
      localisl(i)=isl(i,lsgan)
      enddo
c
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
      call ga_access(ga_s,ils,ils+nsan-1,jls,jls+2,localindex,localld)
      localxs(1:msa,1:3)=>dbl_mb(localindex:localindex+3*msa)
      allocate(localSx(nsan))
      allocate(localSy(nsan))
      allocate(localSz(nsan))
      do i = 1, nsan
        localSx(i) = localxs(i,1)
        localSy(i) = localxs(i,2)
        localSz(i) = localxs(i,3)
      enddo
c      print*,'-----me',me,localSx(1:nsan)
c
      endif
```
 
Comments: load the actual coordinates of solute atoms to local allocated memory from the GA space

Code:
```
      localWmolecules=nwmn
      localWatoms=localWmolecules*nwa
      localSatoms=nsan
c
      totalWmolecules=numw
      totalWatoms=totalWmolecules*nwa
      totalSatoms=nums
c
      offsetWmolecules=nwmnlist(me)
      offsetWatoms=nwmnlist(me)*nwa
      offsetSatoms=nsanlist(me)
```

Comments: compute loca, total and offset in terms of number of atoms rather than moleculars, water is indexed by molecular and nwa is usually 3 (each molecular has 3 atoms), solute is indexed by atom 

Code:
```
      adios_comm_world=mpi_comm_world
      call adios_init_noxml(adios_comm_world, adios_err)
      call adios_set_max_buffer_size(i10)
      call adios_declare_group(m_adios_group, "coordinates", "iter",
     +adios_stat_default, adios_err)
      call adios_select_method(m_adios_group, "MPI", "verbose=2",
     +"",adios_err)
```

Comments: ADIOS initialization

Code:
```
      call adios_define_var(m_adios_group, "LWATOMS", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "TWATOMS", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "OWATOMS", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "LWMOLECULES", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "TWMOLECULES", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "OWMOLECULES", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "LSATOMS", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "TSATOMS", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "OSATOMS", "", i2,
     +"", "", "", varid);
      call adios_define_var(m_adios_group, "WLX", "", adios_double,
     +"LWATOMS", "TWATOMS", "OWATOMS", varid);
      call adios_define_var(m_adios_group, "WLY", "", adios_double,
     +"LWATOMS", "TWATOMS", "OWATOMS", varid);
      call adios_define_var(m_adios_group, "WLZ", "", adios_double,
     +"LWATOMS", "TWATOMS", "OWATOMS", varid);
      call adios_define_var(m_adios_group, "WID", "", adios_integer,
     +"LWMOLECULES", "TWMOLECULES", "OWMOLECULES", varid);
      call adios_define_var(m_adios_group, "SLX", "", adios_double,
     +"LSATOMS", "TSATOMS", "OSATOMS", varid);
      call adios_define_var(m_adios_group, "SLY", "", adios_double,
     +"LSATOMS", "TSATOMS", "OSATOMS", varid);
      call adios_define_var(m_adios_group, "SLZ", "", adios_double,
     +"LSATOMS", "TSATOMS", "OSATOMS", varid);
      call adios_define_var(m_adios_group, "SID", "", adios_integer,
     +"LSATOMS", "TSATOMS", "OSATOMS", varid);
```

Comments: ADIOS definition of varials, for example, WLX, WLY, WLZ represent x, y, z coordinates of each water atom, WID represents the index number of the water atoms

Code:
```
      call adios_open(adios_handle, "coordinates", "nwchem_xyz.bp",
     +"a", adios_comm_world, adios_err)
      call adios_write (adios_handle, "LWATOMS", localWatoms,
     +adios_err)
      call adios_write (adios_handle, "TWATOMS", totalWatoms,
     +adios_err)
      call adios_write (adios_handle, "OWATOMS", offsetWatoms,
     +adios_err)
      call adios_write (adios_handle, "LWMOLECULES", localWmolecules,
     +adios_err)
      call adios_write (adios_handle, "TWMOLECULES", totalWmolecules,
     +adios_err)
      call adios_write (adios_handle, "OWMOLECULES", offsetWmolecules,
     +adios_err)
      call adios_write (adios_handle, "LSATOMS", localSatoms,
     +adios_err)
      call adios_write (adios_handle, "TSATOMS", totalSatoms,
     +adios_err)
      call adios_write (adios_handle, "OSATOMS", offsetSatoms,
     +adios_err)
      call adios_write (adios_handle, "WLX", localWx(:), adios_err)
      call adios_write (adios_handle, "WLY", localWy(:), adios_err)
      call adios_write (adios_handle, "WLZ", localWz(:), adios_err)
      call adios_write (adios_handle, "WID", localiwl(:), adios_err)
      call adios_write (adios_handle, "SLX", localSx(:), adios_err)
      call adios_write (adios_handle, "SLY", localSy(:), adios_err)
      call adios_write (adios_handle, "SLZ", localSz(:), adios_err)
      call adios_write (adios_handle, "SID", localisl(:), adios_err)
      call adios_close (adios_handle, adios_err)
```

Comments: actual write operation of ADIOS for each output timestep

Code:
```
      if(mdstep.eq.maxstep) then
      adiosrank = me
      call adios_finalize(adiosrank, adios_err)
      endif
c
      call util_flush(lfntrj)

      if(nwmn.gt.0) then
      deallocate(localiwl)
      deallocate(localWx)
      deallocate(localWy)
      deallocate(localWz)
      endif

      if(nsan.gt.0) then
      deallocate(localisl)
      deallocate(localSx)
      deallocate(localSy)
      deallocate(localSz)
      endif
c
      deallocate(nwmnlist)
      deallocate(nsanlist)
```

Comments: release allocated memory


