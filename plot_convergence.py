#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 12:15:49 2021

@author: relogu
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.mplot3d import Axes3D

path = '/home/relogu/OneDrive/UNIBO/Magistrale/Physics of Complex Systems/project/convergence/'

onlyfiles = glob.glob(path+'*.out')

conv = pd.read_csv(onlyfiles[0])
i=0
for file in onlyfiles[1:]:
    i+=1
    df = pd.read_csv(file)
    conv = conv.append(df, ignore_index = True)
    
conv['modX'] = np.sqrt(conv['Xprec_i']*conv['Xprec_i'] + \
                       conv['Xprec_j']*conv['Xprec_j'] + \
                       conv['Xprec_k']*conv['Xprec_k']) + 1.19e-7

conv1 = conv[conv['N']==100]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon_v$")
sns.scatterplot(x='beta', y='Vprec', hue='sigma', data=conv1)
plt.axvline(0.5, 0, 1,color='Red')
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\log\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='beta', y='Vprec', hue='sigma', data=conv1)
plt.axvline(0.5, 0, 1,color='Red')
plt.show()
conv1['modX'] = conv1['modX']/conv1['sigma']
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon_x$")
sns.scatterplot(x='beta', y='modX', hue='sigma', data=conv1)
plt.axvline(0.5, 0, 1,color='Red')
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\log\epsilon_x$")
ax.set(yscale="log")
sns.scatterplot(x='beta', y='modX', hue='sigma', data=conv1)
plt.axvline(0.5, 0, 1,color='Red')
plt.show()

conv2=conv[conv['sigma']==1]
#plt.ion()
sns.set(style = "darkgrid")
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x = conv1['beta']
y1 = conv1['Vprec']
y2 = conv1['modX']
z = conv1['sigma']

plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon$")
ax.set_zlabel(r"$\sigma$")

M = int(np.max(y1,0))+2
N = int(np.max(z,0))+1

yy, zz = np.meshgrid(range(M), range(N))
g2 = ax.plot_surface(0.5, yy, zz, alpha=0.2, color='red')

g0 = ax.scatter(x, y1, z)
g1 = ax.scatter(x, y2, z)
plt.legend((g0, g1), (r'$\epsilon_v$', r'$\epsilon_x$'))


plt.show()

conv2=conv[conv['sigma']==1]
#plt.ion()
sns.set(style = "darkgrid")
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x = conv2['beta']
y1 = conv2['Vprec']
y2 = conv2['modX']
z = conv2['N']

ax.set_xlabel(r"$\beta$")
ax.set_ylabel(r"$\epsilon$")
ax.set_zlabel("k")

M = int(np.max(y1,0))+2

yy, zz = np.meshgrid(range(M), range(100))
g2 = ax.plot_surface(0.5, yy, zz, alpha=0.2, color='red')

g0 = ax.scatter(x, y1, z)
g1 = ax.scatter(x, y2, z)
plt.legend((g0, g1), (r'$\epsilon_v$', r'$\epsilon_x$'))

plt.show()

conv3 = conv[conv['N']==100]
conv3 = conv3[conv3['sigma']==1]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon$")
g0 = plt.scatter(x=conv3['beta'], y=conv3['Vprec'])
g1 = plt.scatter(x=conv3['beta'], y=conv3['modX'])
plt.legend((g0, g1), (r'$\epsilon_v$', r'$\epsilon_x$'))
plt.axvline(0.5, 0, 1,color='Red')
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\log\epsilon$")
g0 = plt.scatter(x=conv3['beta'], y=conv3['Vprec'])
g1 = plt.scatter(x=conv3['beta'], y=conv3['modX'])
plt.legend((g0, g1), (r'$\epsilon_v$', r'$\epsilon_x$'))
ax.set(yscale="log")
plt.axvline(0.5, 0, 1,color='Red')
plt.show()