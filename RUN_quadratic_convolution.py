# -*- coding: utf-8 -*-
"""
Created on Fri Aug 04 17:26:28 2017

@author: mhl14
"""

#import training data
import scipy.io as sio
data= sio.loadmat('8_3_pixel.mat');
stim = data['pixel']; 
#stim = float(stim);
spikes = data['response_compressed'];
del data;

#run algorithm
from quadratic_convolution_updated import gradDescent
gradDescent('7_4_pixel',spikes,stim,jack=2,fsize=(16,16),
        pixelNorm=True,extrapSteps=10,filepath=None,model='softplus',
        maxIts=10,maxHours=None,perm=True,overwrite=True,Njack=4,
        start='stim_rand',nlags=1,splits=None,LRType='DecayRate',LRParams = {})

#extract result from valid set
from numpy import load
filename="8_3_pixel_QuadraticSoftPlus_valid_1.dat.npy";
data1 = load(filename);

filename="8_3_pixel_QuadraticSoftPlus_valid_2.dat.npy";
data2 = load(filename);

filename="8_3_pixel_QuadraticSoftPlus_valid_3.dat.npy";
data3 = load(filename);

filename="8_3_pixel_QuadraticSoftPlus_valid_4.dat.npy";
data4 = load(filename);

#average result
a1 = (data1[0] +data2[0]+data3[0]+data4[0])/4
v1 = (data1[1] +data2[1]+data3[1]+data4[1])/4
J = (data1[2]+data2[2]+data3[2]+data4[2])/4
a2 = (data1[3] +data2[3]+data3[3]+data4[3])/4
v2 = (data1[4]+data2[4]+data3[4]+data4[4])/4
d = (data1[5] +data2[5]+data3[5]+data4[5])/4

#plot result
import matplotlib.pyplot as plt
plt.imshow(v1[:,:,0].reshape([16,16]), cmap='jet')
sio.savemat('h.mat',{"h_mean":v1.reshape([16,16])})
plt.imshow(J[:,:,0,:,:,0].reshape ([256,256]), cmap='jet')
sio.savemat('J.mat',{"J_mean":J.reshape([1,256*256])})
#plt.imshow(v2[:,:,:].reshape ([17,65]), cmap='jet')


#import testing data
import scipy.io as sio
data= sio.loadmat('testing.mat');
stim = data['pixel']; 
#stim = float(stim);
spikes = data['tresponse'];
del data;

#calculate repsonse
from quadratic_convolution_updated import Resp,respSP,gridStim,normStim
import matplotlib.pyplot as plt
P= a1,v1,J,a2,v2,d; 
fsize=(16,16,1); nlags=1;pixelNorm=True;
stim = normStim(stim,pixelNorm)[0]
S = gridStim(stim,fsize,nlags)
R = Resp(S,P,respSP)
plt.plot(spikes); plt.plot(R)
plt.legend(['real spikes','predicted spikes'],loc=1)

#calculate correlation coefficient
from numpy import corrcoef
corr=corrcoef(spikes[1:3994,0], R[1:3994])[0,1]






#old pesudo code



#run predicting algorithm
#from numpy import tensordot as tdot
#from numpy import zeros,transpose
#from quadratic_convolution_updated import normStim,gridStim,logistic,softPlus,respSP
#f1,f2 = logistic,softPlus
#ndim = v1.ndim
#r1=zeros((6143,1,1,1))
#r2=zeros((6143))
#prob=zeros((6143))
#for j in range(0,2000):
#    r1[j,...] = f1(a1+tdot(v1,S[j,...],2*(range(ndim),))+(tdot(J,S[j,...],
#      2*(range(ndim),))*S).sum(tuple(range(ndim)))) 
#    r2[j] = f2(a2+(r1[j,...]*transpose(v2)).sum())
#    prob[j]=d*r2[j]

    
#from numpy.lib.stride_tricks import as_strided
#from numpy.random import RandomState
#from Params import Params
#from numpy import ndarray
#P=Params(a1,v1,J,a2,v2,d)
#Ntrials = stim.shape[-1]-nlags+1
#p = ndarray(stim.shape[-1])
#Nvalid = Ntrials / Njack
#RS = RandomState(0)
#p = ndarray(stim.shape[-1])
#R = array([respSP(S[j,...],P) for j in p])

#with open(filename) as f:
#    for line in f:
#        data= line.append([float(n) for n in line.strip().split('          ')])

#n=31;
#data2= [line[i:i+n] for i in range(0, len(line), n)]

#data=[len(line)]    
#for i in range(1000):
#        data=line.strip().split("          ")
#datalist= asarray(dataset_array)
            
#dataset_array = []
#for item in line.split("          "):
#    dataset_array.append(item)
            
#data=line.split('           ')
    
#data=append(line, length(line))
#data=float(line)

#train1 = fromfile("inputV2_QuadraticSoftPlus_train_1.dat", sep=" ",count=-1)
#valid4 = fromfile("inputV2_QuadraticSoftPlus_train_4.dat", dtype=float, count=-1)

#train1 = str(train1);
#train1s = list(train1[0]);

#for line in train1:
#        data = get_value(line)
        
#sio.savemat('train4.mat', {"train4":train4})
#sio.savemat('valid4.mat', {"valid4":valid4})

#jack=1;fsize=(16,16);pixelNorm=True;extrapSteps=10;filepath=None;
#model='softplus';maxIts=None;
#maxHours=None;perm=True;overwrite=True;Njack=4;start='stim_uniform';
#nlags=10;splits=None;LRType='DecayRate';LRParams = {};

#from numpy import fromfile,prod,float64, split, asarray
#from matplotlib import pyplot as py 
#import string 
#from differential_evolution import DETrainer
#from Params import Params
#import math

#data = sio.loadmat('input64x1x16875.mat');
#stim = data['pixel'];
#spikes = data['response_compressed'];
#stim = data['mov2']; 
#stim = float64(stim);
#spikes = data['psth2'];

#fsize=(16,16);
#gsize=(4,5);
#shapes = [(1,),fsize,2*fsize,(1,),gsize,(1,)]

#assert shapes is not None
#dtype = tuple
#params=dtype(file);
##with open(params,'r') as f:
#    train = [p.copy() for p in params]
##    [fromfile(f,count=prod(s),dtype=dtype).reshape(s) for s in shapes]
#del file

#file = "inputV2_QuadraticSoftPlus_train_1.dat"
#with open(file,'rb') as f1:
#    train = fromfile(f1, dtype=float, count=-1)
##    train = [fromfile(file,count=prod(s),dtype=float).reshape(s) for s in shapes]
##    train = [fromfile(file,count=prod(s),dtype=float64).reshape(s) for s in shapes]
##    sio.savemat('train1.mat', {"train4":train})
#del file

#file = "inputV2_QuadraticSoftPlus_valid_1.dat"
#with open(file,'rb') as f2: 
##    valid = [fromfile(file,count=prod(s),dtype=float).reshape(s) for s in shapes]
#    valid = fromfile(f2, dtype=float, count=-1)
##    sio.savemat('valid1.mat', {"valid4":valid})
#del file

#with open(file,'r') as f1:
#    size = [prod(s) for s in shapes]
#    size.insert(0,0)
#    assert train.size == sum(size)
#    ind = [sum(size[:j+1]) for j in range(len(size))]
#    train = [train[ind[j]:ind[j+1]].reshape(shapes[j]).astype(float) for j in range(len(size)-1)]

    