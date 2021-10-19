import numpy
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

eV=27.2114
"function that does everything"
def GC_QMC(data,Vdata,mbeg,mend,mstep,Theta,nelec):
    mu = numpy.arange(mbeg,mend,mstep)
    nn = numpy.arange(-nelec+2,nelec-1,2)
    N = numpy.zeros((mu.shape[0],Theta))
    E = numpy.zeros((mu.shape[0],Theta))
    VE = numpy.zeros((mu.shape[0],Theta))

    for i,m in enumerate(mu):
        for theta in range(Theta):
            a=True
            for n,ni in enumerate(nn):
                if (data[n+1,theta]-data[n+1-1,theta])/2.<=m and (data[n+1+1,theta]-data[n+1,theta])/2.>m and a:
#                    a = False
                    N[i,theta]=96+ni
                    E[i,theta]=data[n+1,theta]
                    VE[i,theta]=Vdata[n+1,theta]
    N[N == 0.] = numpy.nan
    E[E == 0.] = numpy.nan
    VN = numpy.nanmean(N**2, axis=1)
    N = numpy.nanmean(N, axis=1)
    E = numpy.nanmean(E, axis=1)
    VE = numpy.sqrt(numpy.sum(VE**2, axis=1))/Theta #to verify
    VN = numpy.sqrt((VN-N*N)/(Theta-1))
    a=0.
    for i,m in enumerate(mu):
        if i<mu.shape[0]-2:
            if N[i]==96.: #E[i] == E[i+1] == E[i+2] == E[i+3]:
                #print m*27.2114
                a+=mstep
    print(a*27.2114)
    plt.plot(mu,N,'o')
    plt.show()
    return numpy.stack((mu,N,E,VE,VN)) #,(a)*27.2114

"load data for different P and T"
M=32
nconf=1
name='rmc'
nelec=4
corr_name='rmc_corr'
data=[]
Vdata=[]
for i in range(nconf):
    data.append(numpy.stack((numpy.average(numpy.loadtxt(('c%d/%s%d.np.1.excite' % (i+1,name,k)),
                                                               usecols=(range(0,M*2,2)), skiprows=100),axis=0)
                              for k in range(-nelec,nelec+1,2))
                            ))
    Vdata.append(numpy.stack((numpy.std(numpy.loadtxt(('c%d/%s%d.np.1.excite' % (i+1,name, k)),
                                                          usecols=(range(0,M*2,2)), skiprows=100),axis=0)
                               for k in range(-nelec,nelec+1,2))
                             ))
corr=numpy.loadtxt('%s.out' % corr_name)
#saving size corrected energies for each twist, number of electrons, and nuclear configuation
#corr * 2 factor of 2 comes from kinetic energy
data_corr=(numpy.asarray(data)+numpy.transpose(numpy.tile(corr*2,(nconf,M,1)),axes=(0,2,1)))
data=numpy.asarray(data)
numpy.save('%s_corr' % name,data_corr)
for i in range(nelec+1):
    plt.plot(data[nconf-1,i],'o')

plt.show()
GC_data=[]
GC_data_ns=[]
for i in range(nconf):
    GC_data.append(GC_QMC(data_corr[i],Vdata[i],mbeg=0.0,mend=1.3,mstep=0.007,Theta=M,nelec=nelec))

    GC_data_ns.append(GC_QMC(data[i],Vdata[i],mbeg=0.0,mend=1.3,mstep=0.007,Theta=M,nelec=nelec))

GC_data=numpy.asarray(GC_data)
numpy.savetxt('%s_gcdat' % name ,numpy.transpose(GC_data[0]))

