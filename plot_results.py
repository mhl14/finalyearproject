# -*- coding: utf-8 -*-
"""
Created on Tue May 29 16:14:22 2018

@author: edwardlaw
"""
#extract result from train set
from numpy import load
filename="8_3_pixel_QuadraticSoftPlus_train_1.dat.npy";
data1 = load(filename);

filename="8_3_pixel_QuadraticSoftPlus_train_2.dat.npy";
data2 = load(filename);

filename="8_3_pixel_QuadraticSoftPlus_train_3.dat.npy";
data3 = load(filename);

filename="8_3_pixel_QuadraticSoftPlus_train_4.dat.npy";
data4 = load(filename);

#avrage result
a1 = (data1[0] +data2[0]+data3[0]+data4[0])/4
v1 = (data1[1] +data2[1]+data3[1]+data4[1])/4
J = (data1[2]+data2[2]+data3[2]+data4[2])/4
a2 = (data1[3] +data2[3]+data3[3]+data4[3])/4
v2 = (data1[4]+data2[4]+data3[4]+data4[4])/4
d = (data1[5] +data2[5]+data3[5]+data4[5])/4

#plot results
import matplotlib.pyplot as plt
import scipy.io as sio
plt.imshow(v1[:,:,0].reshape([16,16]), cmap='jet')
sio.savemat('h.mat',{"h_mean":v1})
plt.imshow(J[:,:,0,:,:,0].reshape ([256,256]), cmap='jet')
sio.savemat('J2.mat',{"J_mean":J.reshape([256,256])})
plt.imshow(v2[:,:,:].reshape ([17,13]), cmap='jet')


#plot individual results
a1=data1[0]
v1=data1[1]
J=data1[2]
a2=data1[3]
v2=data1[4]
d=data1[5]
import matplotlib.pyplot as plt
plt.imshow(v1[:,:,0].reshape ([16,20]), cmap='jet')
plt.imshow(J[:,:,0,:,:,0].reshape ([256, 400]), cmap='jet')
plt.imshow(v2[:,:,:].reshape ([17,13]), cmap='jet')

a1=data2[0]
v1=data2[1]
J=data2[2]
a2=data2[3]
v2=data2[4]
d=data2[5]
plt.imshow(v1[:,:,0].reshape ([16,20]), cmap='jet')
plt.imshow(J[:,:,0,:,:,0].reshape ([256,400]), cmap='jet')
plt.imshow(v2[:,:,:].reshape ([17,65]), cmap='jet')

a1=data3[0]
v1=data3[1]
J=data3[2]
a2=data3[3]
v2=data3[4]
d=data3[5]
plt.imshow(v1[:,:,0].reshape ([16,20]), cmap='jet')
plt.imshow(J[:,:,0,:,:,0].reshape ([256,400]), cmap='jet')
plt.imshow(v2[:,:,:].reshape ([17,65]), cmap='jet')

a1=data4[0]
v1=data4[1]
J=data4[2]
a2=data4[3]
v2=data4[4]
d=data4[5]
plt.imshow(v1[:,:,0].reshape ([16,20]), cmap='jet')
plt.imshow(J[:,:,0,:,:,0].reshape ([256,400]), cmap='jet')
plt.imshow(v2[:,:,:].reshape ([17,65]), cmap='jet')



