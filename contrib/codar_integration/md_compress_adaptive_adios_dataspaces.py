import numpy as np
import math

from sklearn.decomposition import TruncatedSVD
from sklearn.utils.extmath import svd_flip, randomized_svd

from scipy.sparse.linalg import svds
from scipy.misc import imfilter
import scipy.ndimage as ndi

import scipy.io

from heapq import heappush, heappop

from random import uniform

import matplotlib
import matplotlib.pyplot as plt


from pprint import pprint


# ADIOS and MPI modules
import adios_mpi as ad
from mpi4py import MPI
import sys


def  svd_wrapper (Y, k, method='svds'):
    if method is 'svds':
        Ut, St, Vt = svds(Y, k)
        idx = np.argsort(St)[::-1]
        St = St[idx] # have issue with sorting zero singular values
        Ut, Vt = svd_flip(Ut[:, idx], Vt[idx])
    elif method is 'random':
        Ut, St, Vt = randomized_svd(Y, k)
    else:
        Ut, St, Vt = np.linalg.svd(Y, full_matrices=False)
        # now truncate it to k
        Ut = Ut[:, :k]
        St = np.diag(St[:k])
        Vt = Vt[:k, :]

    return Ut, St, Vt


class WeightedReservoirSampler:
    def __init__(self, reservoir_size):
        self.m = reservoir_size

        self.reservoir = []

    def add(self, weight, data):
        rw = uniform(0.0,1.0)**(1.0/weight)
        if(len(self.reservoir)< self.m):
            heappush(self.reservoir, (rw,data))
        elif rw > self.reservoir[0][0]:
            heappop(self.reservoir)
            heappush(self.reservoir, (rw,data))

class MDTrSampler:
    def __init__(self, n_atoms, n_dim = 2, conv_size = 50, n_samples=1000, batch_size=100, manifold_size=64):
        self.n_dim = n_dim
        self.conv_size = conv_size
        self.n_atoms = n_atoms
        self.batch_size = batch_size
        self.l = manifold_size
        self.n_samples = n_samples
        self.bq = np.zeros((batch_size,n_atoms,3)) # batch queue
        self.total_samples = 0
        self.bq_index = 0
        self.Btp = np.zeros((n_atoms*3,self.l))
        self.strm_smplr = WeightedReservoirSampler(n_samples)
        self.last_Vt = None

    def traj_char(self, c, s):
        x = c.shape[0]
        # scaled data <-- not so much meaningful as we're doing l2 normalization later
        ps = np.mat(c[:,:self.n_dim]) * np.mat(np.diag(s[:self.n_dim]))
        # gradient or changes of each data point
        psd = ps[0:(x-1),:self.n_dim] - ps[1:x, :self.n_dim];

        # convoluted gradient --> smoothed gradient
        psdm = ndi.convolve(psd, np.ones((self.conv_size,2))/self.conv_size*2);
        # L2 normalized smoothed gradients --> now we focus on the angle of gradient only as atoms may move different speed
        npsdm = np.divide(psdm, np.mat(np.sum(np.abs(psdm)**2,axis=-1)**(1./2)).T  * np.mat(np.ones((1,self.n_dim))))
        # Angle smoothing with the assumption that the overall angle can not be radically changed
        nmpsdm = ndi.convolve(npsdm, np.ones((self.conv_size,2))/self.conv_size*2);
        # Smoothed angle changes (accelerations)
        psdd = nmpsdm[0:(x-2)] - nmpsdm[1:(x-1)]
        # results for gradient of changes
        psdu = np.sum(abs(psd),axis=1)
        # results of smoothed normalized acceleration
        psddu = np.sum(abs(psdd),axis=1)
        prob_dist = (abs(psddu) / np.sum(psddu))
        return nmpsdm, psddu, prob_dist

    def strmML(self, t2):
        n_t = t2.shape[1]
        Ct = np.concatenate( (self.Btp, t2), axis = 1)
        Ct = Ct[:, ~(Ct==0).all(0)]

        Ut, St, Vt = svd_wrapper(Ct, self.l) #SVD_l(matrix)
        Ut_l = Ut[:, :self.l]
        St_l = St[:self.l] - St[self.l-1] # to be adaptive, added singular substraction
        Vt_l = Vt[:self.l, -n_t:]
        self.Btp= np.dot(Ut_l, np.diag(St_l))

        return Ut_l, St_l, Vt_l

    def batch_sampling(self, trace):
        (x,y,z) = trace.shape
        self.sampling_rate = float(self.n_samples) / float(x)
        t2 = trace.reshape((x, y*z),order='C')
        c,s,v = svd_wrapper(t2, 2, method='random')
        nmpsdm, psddu, prob_dist = self.traj_char(c,s)
        total = 10**-10 + prob_dist[0]
        sampling_entries = int(x * self.sampling_rate)
        target = np.zeros((sampling_entries, y, z))
        time_stamps = np.zeros(sampling_entries)
        #output initialization: the first one should be added always
        target[0,:,:] = trace[0,:,:]
        time_stamps[0] = 0;
        t_idx = 1;
        for i in range(2,x):
            total = total + prob_dist[i-2];
            if(total > 1 / float(sampling_entries - 1)):
                target[t_idx,:,:] = trace[i,:,:]
                time_stamps[t_idx] = i
                total = total - 1 / float(sampling_entries-1);
                t_idx = t_idx + 1;
        return target, time_stamps


    def adaptive_sampling_update(self):
        t2 = self.bq[0:self.bq_index,:,:].reshape((self.bq_index, self.n_atoms*3),order='C').T # 3 is for three dimensional coordinates
        Ut, s, VTt = self.strmML(t2)
        if self.last_Vt is None: # for the first time
            c = VTt.T
            self.strm_smplr.add(10**10,[0, self.bq[0,:,:]])
            nmpsdm, psddu, prob_dist = self.traj_char(c, s)
            sidx = self.total_samples - self.bq_index
            for i in range(2,self.bq_index):
                self.strm_smplr.add(psddu[i-2], [sidx+i, self.bq[i,:,:]])
        else:
            c =  np.concatenate( (self.last_Vt, VTt.T), axis = 0)
            nmpsdm, psddu, prob_dist = self.traj_char(c, s)
            sidx = self.total_samples - self.bq_index
            for i in range(0,self.bq_index):
                self.strm_smplr.add(psddu[i], [sidx+i, self.bq[i,:,:]])

        self.last_Vt = VTt.T[-2:,:] # for last two to start compare against


    def adaptive_sampling_step(self, dataframe_t):
        #pprint(dataframe_t)
        self.total_samples = self.total_samples + 1
        self.bq[self.bq_index] = dataframe_t
        self.bq_index = self.bq_index+1
        if self.bq_index == (self.batch_size):
            self.adaptive_sampling_update()
            self.bq_index = 0



np.set_printoptions(threshold=np.inf)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


ad.read_init("DATASPACES", comm, "verbose=3;")

f = ad.file("nwchem_xyz.bp", "DATASPACES", comm, True, timeout_sec = 20.0)

print("RANK:",rank)

streaming_next = 0
current_step = 0

v = f.vars['TWATOMS']
twatoms = v.read(nsteps=1)
v = f.vars['TSATOMS']
tsatoms = v.read(nsteps=1)

natoms=twatoms+tsatoms

if rank == 0:
    mds2 = MDTrSampler(natoms, 2, conv_size = 10, n_samples=32, batch_size=50, manifold_size=25)

while streaming_next == 0:
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

        master_xyz = np.array([ordered_x,ordered_y,ordered_z]).T
        mds2.adaptive_sampling_step(master_xyz[:,:])

    comm.Barrier()
    streaming_next = f.advance()
    current_step += 1


f.close()
ad.read_finalize("DATASPACES")
comm.Barrier()



if rank == 0:
    adaptive_samples = sorted([ mds2.strm_smplr.reservoir[i][1][0] for i in range(mds2.n_samples)])
    pprint(adaptive_samples)

