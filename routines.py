
import numpy
import math
import random

def AverageAndVar(data):
    N = data.shape[0]
    average = numpy.average(data,axis = 0)
    Var = numpy.sum(data*data, axis=0)/N - average*average
    return numpy.concatenate(([average], [numpy.sqrt(Var/(N/2.))]), axis=0)


def BlockingAver(data,k):
    results = [AverageAndVar(data)]
    for i in range(k):
        N = data.shape[0]
        if N%2 != 0:
            data = numpy.delete(data,N-1,0)
            N = data.shape[0]
        data = numpy.mean(data.reshape(N/2,2,5),1)
        results = numpy.concatenate((results, [AverageAndVar(data)]), axis=0)
    return results
        
def CorrFunc(data,k):
    N = data.shape[0]
    tmp = AverageAndVar(data)
    average = tmp[0]
    Var = tmp[1]**2*N
    Corr = [Var]
    Tau = [numpy.repeat(1./2., data.shape[1])]# + numpy.sum(Corr/Var,axis=0)]
    for i in range(1,k):
        Corr = numpy.concatenate((Corr, [numpy.sum((data[0:N-i]-average)*(data[i:N]-average),axis=0)/(N-i)]))
        Tau = numpy.concatenate((Tau,[1./2. + numpy.sum(Corr/Var,axis=0)]),axis=0)
        Error = numpy.sqrt(Var*2*Tau/N)
    return Corr, Tau, Error

def JackKnife(bdata):
    JJK = [numpy.zeros(5)]
    N = bdata.shape[0]
    for i in range(N):
        JiJK = numpy.average(numpy.delete(bdata,i,0),axis=0)
        JJK = numpy.concatenate((JJK,[JiJK/JiJK[0]]),axis=0)
    JJK = numpy.delete(JJK,0,0)
    JJK = numpy.delete(JJK,0,1)
    Ji = numpy.average(bdata,axis=0)
    JestJKcap = N*numpy.delete(Ji/Ji[0],0) - (N-1)*numpy.average(JJK,axis=0)
    VarJJK = numpy.average((-JJK+JestJKcap)*(-JJK+JestJKcap),axis=0)
    ErrorJK = numpy.sqrt((N-1)*VarJJK)
    return ErrorJK

def EnergyTotal(data,N,L,rc):
    E = 0
    for i in range(N-1):
        for j in range(i+1,N):
            r = data[i] - data[j]
            r = r - L*numpy.rint(r/L)
            r = numpy.sqrt(sum(r*r))
            if r < L*rc:
                E = E + numpy.exp(-r)/r
    return E

def EnergyPP(i,data,N,L,rc):
    E = 0
    for j in range(N):
        if j!=i:
            r = data[i] - data[j]
            r = r - L*numpy.rint(r/L)
            r = numpy.sqrt(sum(r*r))
            if r < L*rc:
                E = E + numpy.exp(-r)/r
    return E

def VirialCalc(data,N,L,rc):
    P = 0    
    for i in range(N-1):
        for j in range(i+1,N):
            r = data[i] - data[j]
            r = r - L*numpy.rint(r/L)
            r = numpy.sqrt(sum(r*r))
            if r < L*rc:
                P = P + numpy.exp(-r)*(1./r+1)
    return P

def MC_simulations(beta,N_iter,N,L,d,rc):
    A = 0
    P = L*numpy.random.random_sample((N,3))
    E = EnergyTotal(P,N,L,rc)
    E_array = numpy.zeros(N_iter)
    Virial = numpy.zeros(N_iter)
    delta = L/d
    hist = numpy.zeros(100)
    for i in range(N_iter):
        Virial[i] = VirialCalc(P,N,L,rc)
        for j in range(N):
            d = delta*(numpy.random.random_sample(3)-0.5)
            E_old = EnergyPP(j,P,N,L,rc)
            P[j] = P[j] + d
            E_new = EnergyPP(j,P,N,L,rc)
            if random.random() < math.exp(-beta*(E_new-E_old)):
                if i > 100 and i <= 300:
                    A += 1.
                E = E + E_new - E_old
            else:
                P[j] = P[j] - d
        E_array[i] = E * 1.5 / N
        if i%500==0:
            print(i)
        if i>(N_iter-100):
            hist += g_of_r(P,N,L,rc)
    return E_array, Virial, A/(200*N), hist/100

#calculating g(r)
def g_of_r(data,N,L,rc): # pass coordinates
    Nhist = 100
    dr = L/(2.*Nhist)
    hist = numpy.zeros(Nhist)    # creating histogram array
    for i in range(N-1):
        for j in range(i+1,N):
            r = data[i] - data[j]
            r = r - L*numpy.rint(r/L)
            r = numpy.sqrt(sum(r*r))
            if r < L*rc:
                ig = math.floor(r/dr)
                hist[ig] = hist[ig] + 2  # add two particles for each interaction
 
    for i in range(Nhist):
        hist[i] = hist[i]/(N*(N/L**3)*(4.0/3.0)*math.pi*dr**3*((i+1)**3-i**3))  # normalazing 
    return hist