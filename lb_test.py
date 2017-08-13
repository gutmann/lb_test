#!/usr/bin/env python
#
# 2D Lattice Boltzmann (BGK) model of a fluid.
#  c3  c2   c1  D2Q9 model. At each timestep, particle densities propagate
#    \  |  /    outwards in the directions indicated in the figure. An
#  c4 -c9 - c0  equivalent 'equilibrium' density is found, and the densities
#    /  |  \    relax towards that state, in a proportion governed by omega.
#  c5  c6   c7      Iain Haslam, March 2006.
#
import numpy as np
import matplotlib.pyplot as plt

# set up constants for calculations
omega   = 0.5         # relaxation term
density = 1.0
t1      = 4.0/9       # non-movement coefficient
t2      = 1.0/9       # moving perpendicular
t3      = 1.0/36      # moving diagonal
c_squ   = 1/3.0
nx      = 110          # number of x grid cells
ny      = 110          # number of y grid cells


def latb2d(F=None, BOUND=None, slow=False, ux=None, uy=None, maxsteps=4000, minsteps=100):
    # F is the array of probability of particle movement in each direction
    # nx x ny x 9 possible directions of movement (or non-movement)
    if F is None:
        F = np.zeros((9,nx,ny)) + density / 9.0
    FEQ=F # F equilibrium initially equals F

    # matrix size
    msize = nx * ny

    # CI is an array of offsets into the matrix.  each element moves forward by one nx x ny grid
    #  note sure if this should be np.arange(8) or (1,9), or???
    CI = np.arange(9) * msize

    # what is bound?
    if BOUND is None:
        # BOUND = np.random.rand(nx,ny) > 0.8
        BOUND = np.zeros((nx,ny))
        height = int(nx/2)
        position = int(ny/5)
        BOUND[height/2:height,position:position+2] = 1
        BOUND[0,:] = 1
        BOUND[-1,:] = 1
        BOUND = BOUND > 0.5

    # index into the data array "ON" = grid cells on a boundary(?)
    ON = np.where(BOUND)
    ONx = ON[0]
    ONy = ON[1]
    ON = ONx*ny + ONy

    print(ON.shape)
    TO_REFLECT=np.array([ON+CI[0], ON+CI[1], ON+CI[2], ON+CI[3],
                ON+CI[4], ON+CI[5], ON+CI[6], ON+CI[7]])
    REFLECTED= np.array([ON+CI[4], ON+CI[5], ON+CI[6], ON+CI[7],
                ON+CI[0], ON+CI[1], ON+CI[2], ON+CI[3]])
    TO_REFLECT = TO_REFLECT.flat[:]
    REFLECTED = REFLECTED.flat[:]

    # F[[0,1,7],:,:]*=1.1
    # F[[0,1,7],:,:]*=1.1
    # F[[3,4,5],:,:]/=1.1
    # F[[3,4,5],:,:]/=1.1
    F[[0,1,7],0, 1:-1]*=1.1
    F[[0,1,7],-1,1:-1]*=1.1
    F[[3,4,5],0, 1:-1]/=1.1
    F[[3,4,5],-1,1:-1]/=1.1
    F[[0,1,7],:, 0]*=1.1
    F[[0,1,7],:,-1]*=1.1
    F[[3,4,5],:, 0]/=1.1
    F[[3,4,5],:,-1]/=1.1
    # F[4]/=1.1
    # F[4]/=1.1
    # F[0]*=1.1
    # F[0]*=1.1
    # F[:,0,:]=0
    # F[:,-1,:]=0
    # F[[5,6,7],0,:]=0
    # F[[5,6,7],-1,:]=0
    # F[[3,2,1],0,:]=0
    # F[[3,2,1],-1,:]=0

    avu = 1
    prevavu = 1
    ts = 0
    # increase inlet pressure more if there are fewer inlet nodes active
    deltaU = 1e-7 #* np.sqrt((1 / max(1e-10,1- np.mean(BOUND[0,:]))))
    numactivenodes = np.sum( ~BOUND )

    while ((ts < maxsteps) and (1e-10 < abs((prevavu-avu)/avu))) or (ts < minsteps):

        # Propagate
        F[0,1:-1,1:-1] = F[0,2:  ,1:-1]
        F[1,1:-1,1:-1] = F[1,2:  ,2:  ]
        F[2,1:-1,1:-1] = F[2,1:-1,2:  ]
        F[3,1:-1,1:-1] = F[3,0:-2,2:  ]
        F[4,1:-1,1:-1] = F[4,0:-2,1:-1]
        F[5,1:-1,1:-1] = F[5,0:-2,0:-2]
        F[6,1:-1,1:-1] = F[6,1:-1,0:-2]
        F[7,1:-1,1:-1] = F[7,2:  ,0:-2]

        BOUNCEDBACK = F.flat[TO_REFLECT] # Densities bouncing back at next timestep
        DENSITY = np.sum(F,axis=0)

        UX = (np.sum(F[[0, 1, 7], :,:],axis=0) - np.sum(F[[3, 4, 5], :,:], axis=0)) / DENSITY
        UY = (np.sum(F[[1, 2, 3], :,:],axis=0) - np.sum(F[[5, 6, 7], :,:], axis=0)) / DENSITY

        UX[:,0] += deltaU # Increase inlet pressure

        UX.flat[ON] = 0
        UY.flat[ON] = 0
        DENSITY.flat[ON] = 0

        U_SQU = UX**2 + UY**2
        U_C2  = UX + UY
        U_C4  = -1*UX + UY
        U_C6  = -1*U_C2
        U_C8  = -1*U_C4

        # Calculate equilibrium distribution: stationary
        FEQ[8,:,:]=t1 * DENSITY*(1-U_SQU/(2*c_squ))
        # nearest-neighbours
        FEQ[0,:,:]=t2 * DENSITY * (1 + UX/c_squ + 0.5*(UX/c_squ)**2 - U_SQU/(2*c_squ))
        FEQ[2,:,:]=t2 * DENSITY * (1 + UY/c_squ + 0.5*(UY/c_squ)**2 - U_SQU/(2*c_squ))
        FEQ[4,:,:]=t2 * DENSITY * (1 - UX/c_squ + 0.5*(UX/c_squ)**2 - U_SQU/(2*c_squ))
        FEQ[6,:,:]=t2 * DENSITY * (1 - UY/c_squ + 0.5*(UY/c_squ)**2 - U_SQU/(2*c_squ))

        # corner neighbours
        FEQ[1,:,:] = t3 * DENSITY * (1 + U_C2/c_squ + 0.5*(U_C2/c_squ)**2 - U_SQU/(2*c_squ))
        FEQ[3,:,:] = t3 * DENSITY * (1 + U_C4/c_squ + 0.5*(U_C4/c_squ)**2 - U_SQU/(2*c_squ))
        FEQ[5,:,:] = t3 * DENSITY * (1 + U_C6/c_squ + 0.5*(U_C6/c_squ)**2 - U_SQU/(2*c_squ))
        FEQ[7,:,:] = t3 * DENSITY * (1 + U_C8/c_squ + 0.5*(U_C8/c_squ)**2 - U_SQU/(2*c_squ))

        F = omega*FEQ+(1-omega)*F

        F.flat[REFLECTED] = BOUNCEDBACK

        prevavu = avu
        avu = max(np.sum(UX) / numactivenodes, 1e-10)
        ts += 1

        if ((ts % 10) == 0) and slow:
            colormax = np.max(np.abs(UX))
            plt.clf()
            plt.subplot(2,2,1)
            plt.imshow(UX, cmap=plt.cm.RdBu)
            plt.clim(-colormax,colormax)
            plt.colorbar()

            plt.subplot(2,2,2)
            plt.imshow(UY, cmap=plt.cm.RdBu)
            plt.clim(-colormax,colormax)
            plt.colorbar()

            plt.subplot(2,2,3)
            plt.quiver(UX[::5,::5], UY[::5,::5])

            plt.subplot(2,2,4)
            plt.imshow(DENSITY)
            plt.clim(0.8,1.0)
            plt.colorbar()

            plt.show()
            plt.pause(0.1)

    return F


def main():
    print("This file contains python code translated from older IDL and Matlab Lattice Boltzmann code. ")
    print("Running...")

    z = np.random.rand(ny)*nx/2
    BOUND=np.zeros((nx,ny))
    for i in range(ny):
        z[i] = np.mean(z[max(i-2,0):min(i+2,ny-1)])
        BOUND[:int(z[i]),i]=1
    BOUND = BOUND>0.5


    BOUND[:]=False
    BOUND[ nx/3:-nx/3,ny/2] = True
    F = latb2d(slow=True, BOUND=BOUND)

    print("Completed!")

if __name__ == '__main__':
    main()
