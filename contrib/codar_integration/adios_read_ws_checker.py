#!/usr/bin/env python
import adios_mpi as ad
import numpy as np
from mpi4py import MPI
import sys
import math

np.set_printoptions(threshold=np.inf)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    fn = open('ethanol_md.trj')

ad.read_init("DATASPACES", comm, "verbose=3;")

f = ad.file("nwchem_xyz.bp", "DATASPACES", comm, True, timeout_sec = 20.0)

print("RANK:",rank)

streaming_next = 0
current_step = 0

while streaming_next==0:
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

    comm.Barrier()
    if rank == 0:
        print(">>STEP:", current_step)
    comm.Barrier()
#    print("Rank - Tmolecules, Lmolecules: ",rank,nmolecules,localmolecules)
#    print("Rank - Tatoms, Latoms: ",rank,natoms,localatoms)



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
#        if current_step == 0:
#            print(ordered_x)

        
        line = fn.readline()
        words = line.split()
        while words[0] != 'TFFFTFFF' and line:
            line = fn.readline()
            words = line.split()
            print ("LINE: ",line)
        for i in range(Watoms+Satoms):
            line = fn.readline()
            words = line.split()
            print ("Data: ",words[0],round(ordered_x[i],3))
            if float(words[0]) != round(ordered_x[i],3):
                print('Wrong LX: step - position',current_step,i)
            if float(words[1]) != round(ordered_y[i],3):
                print('Wrong LY: step - position',current_step,i)
            if float(words[2]) != round(ordered_z[i],3):
                print('Wrong LZ: step - position',current_step,i)
#            print(words[0],words[1],words[2])

 

    comm.Barrier()
    streaming_next = f.advance()
    current_step += 1



f.close()

if rank == 0:
    fn.close()

#print("\n>>> Done.\n")





