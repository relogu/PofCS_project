#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 12:15:49 2021

@author: relogu
"""
#%% importations
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import time
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
#%% functions
proj_path = '/home/relogu/OneDrive/UNIBO/Magistrale/Physics of Complex Systems/project/'

def make_array(df, colname=None, N=None):
    """Util to convert array entries of the results."""
    if N is None and colname is None:
        return np.fromstring(df,sep=' ')
    else :
        return np.fromstring(df[colname],sep=' ').reshape(df[N], 3)
    
def build_dataframe(path, last_only=False):
    """Retrieve the dataframe containing the results of the simulations."""
    onlyfiles = glob.glob(path+'*.dat')
    print('File retrieved')
    print('Reading files')
    df = pd.read_csv(onlyfiles[0], index_col=False)
    if last_only: df = df[-1:]
    i=0
    for file in onlyfiles[1:]:
        i+=1
        if last_only: df = df.append(pd.read_csv(file)[-1:], ignore_index = True)
        else: df = df.append(pd.read_csv(file), ignore_index = True)
        
    print('Dataframe built up, operating...')
    df['[space_conv]'] = df['[space_conv]'].apply(make_array)
    df['[vel_conv]'] = df['[vel_conv]'].apply(make_array)
    df['space_norm'] = df['[space_conv]'].apply(np.linalg.norm)
    df['vel_norm'] = df['[vel_conv]'].apply(np.linalg.norm)
    df.loc[df['space_norm']==0, 'space_norm']+=1.0e-16
    df.loc[df['vel_norm']==0, 'vel_norm']+=1.0e-16
    df['[[positions]]'] = df.apply(make_array, axis=1, colname='[[positions]]', N='N')
    df['[[velocities]]'] = df.apply(make_array, axis=1, colname='[[velocities]]', N='N')
    df['tau'] = df['iter']/np.power(df['sigma'], 2*df['beta'])
    df['[[Positions]]'] = df['[[positions]]']/df['sigma']
    df['[[Velocities]]'] = df['[[velocities]]']/df['sigma']
    df['[Space_conv]'] = df['[space_conv]']/df['sigma']
    df['[Vel_conv]'] = df['[vel_conv]']/df['sigma']
    df['Space_norm'] = df['[Space_conv]'].apply(np.linalg.norm)
    df['Vel_norm'] = df['[Vel_conv]'].apply(np.linalg.norm)
    print('Done.')
    return df

def traslate_to_baricenter(list_of_iterations):
    """Traslate the given list of position to its baricenter."""
    for i in range(len(list_of_iterations)):
        bx=0
        by=0
        bz=0
        for j in range(len(list_of_iterations[i])):
            bx+=list_of_iterations[i][j][0]
            by+=list_of_iterations[i][j][1]
            bz+=list_of_iterations[i][j][2]
        bx = bx/list_of_iterations[i].shape[0]
        by = by/list_of_iterations[i].shape[0]
        bz = bz/list_of_iterations[i].shape[0]
        for j in range(len(df[i])):
            list_of_iterations[i][j][0] -= bx
            list_of_iterations[i][j][1] -= by
            list_of_iterations[i][j][2] -= bz
    return list_of_iterations
    

def animate_scatters(iteration, data, scatters):
    """
    Update the data held by the scatter plot and therefore animates it.
    
    Args
    ----
        iteration (int): Current iteration of the animation
        data (list): List of the data positions at each iteration.
        scatters (list): List of all the scatters (One per element)
        
    Returns
    -------
        list: List of scatters (One per element) with new coordinates
    """
    for i in range(data[0].shape[0]):
        scatters[i]._offsets3d = (data[iteration][i,0:1], data[iteration][i,1:2], data[iteration][i,2:])
    return scatters

def positions_evolutions(filename, data, save=False):
    """
    Create the 3D figure and animates it with the input data.
    
    Args
    ----
        data (list): List of the data positions at each iteration.
        save (bool): Whether to save the recording of the animation. (Default to False).
    """
    # Attaching 3D axis to the figure
    fig = plt.figure()
    ax = Axes3D(fig)

    # Initialize scatters    
    scatters = [ ax.scatter(data[0][i,0:1], data[0][i,1:2], data[0][i,2:]) for i in range(data[0].shape[0]) ]

    # Number of iterations
    iterations = len(data)
    for iteration in range(iterations):
        for i in range(data[0].shape[0]):
            scatters[i]._offsets3d = (data[iteration][i,0:1], data[iteration][i,1:2], data[iteration][i,2:])
        
    
    for d in data:
        Mx=0
        My=0
        Mz=0
        mx=1000
        my=1000
        mz=1000
        for dd in d:
            Mx = max(dd[0], Mx)
            My = max(dd[1], My)
            Mz = max(dd[2], Mz)
            mx = min(dd[0], mx)
            my = min(dd[1], my)
            mz = min(dd[2], mz)

    # Setting the axes properties
    ax.set_xlim3d([mx-0.1, Mx+0.1])
    ax.set_xlabel('X')

    ax.set_ylim3d([mx-0.1, Mx+0.1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([my-0.1, My+0.1])
    ax.set_zlabel('Z')

    ax.set_title('3D Flock Evolution')

    # Provide starting angle for the view.
    ax.view_init(25, 10)
    
    ani = animation.FuncAnimation(fig, animate_scatters, iterations, fargs=(data, scatters),
                                  interval=100, blit=True, repeat=True)
    
    if save:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800, extra_args=['-vcodec', 'libx264'])
        ani.save(filename+'.mp4', writer=writer)
        
    #plt.show()
    
def flock_diameter(pos):
    """Compute the flock diameter given the array with position."""
    ret = 0
    for p1 in pos:
        for p2 in pos:
            tmp = p1-p2
            ret = max(np.linalg.norm(tmp), ret)
    return ret
    
    


#%% study on beta

path=proj_path+'results1/'
results = build_dataframe(path, True)

tmp = results[results['N']==10]
beta = tmp[tmp['sigma']==1].reset_index()

space_norms = []
vel_norms = []
betas = []
for b in beta['beta'].unique():
    betas.append(b)
    tmp = beta[beta['beta']==b].reset_index()
    space_norms.append(float(tmp['space_norm']))
    vel_norms.append(float(tmp['vel_norm']))
df = pd.DataFrame({'beta': betas,
                   'space_norm': space_norms,
                   'vel_norm': vel_norms})

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set(yscale="log")
ax.set_ylabel(r'$\left\Vert d^{max}\right\Vert$')
ax.relim()
g0 = plt.scatter(x=df['beta'], y=df['space_norm'], marker='.')
#sns.scatterplot(x='beta', y='space_norm', data=df)
g1 = plt.axvline(0.5, 0, 1,color='Red')
g2 = plt.axhline(float(df.loc[df['beta']==0.5, 'space_norm']), 0, 1,color='Green')
plt.legend((g0, g1, g2),
           (r'$\left\Vert d^{max}\right\Vert(\beta)$',
            r'$\beta = 0.5$',
            r'$\left\Vert d^{max}\right\Vert$'
            +'={:.2e}'.format(float(df.loc[df['beta']==0.5, 'space_norm']))))
           
plt.show()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r'$\left\Vert v^{max}\right\Vert$')
ax.set(yscale="log")
ax.relim()
g0 = plt.scatter(x=df['beta'], y=df['vel_norm'], marker='.')
#sns.scatterplot(x='beta', y='vel_norm', data=df)
g1 = plt.axvline(0.5, 0, 1,color='Red')
g2 = plt.axhline(float(df.loc[df['beta']==0.5, 'vel_norm']),
            0, 1,color='Green')
#g3 = plt.axhline(1, 0, 1,color='darkorange')
plt.legend((g0, g1, g2),#, g3),
           (r'$\left\Vert v^{max}\right\Vert(\beta)$',
            r'$\beta = 0.5$',
            r'$\left\Vert v^{max}\right\Vert$'
            +'={:.2e}'.format(float(df.loc[df['beta']==0.5, 'vel_norm']))))#,
            #r'$\log\left\Vert d^{max}\right\Vert=1.0$'))
plt.show()

del beta
del space_norms
del vel_norms
del betas

#%% study on sigma

path=proj_path+'results3/'
results = build_dataframe(path, True)
tmp = results[results['N']==10]
sigma = tmp[tmp['beta']==0.4].reset_index()

diameters = []
velocities = []
vel_norms = []
sigmas = []
inv_sigmas = []
sigmas_beta = []
for s in sigma['sigma'].unique():
    sigmas.append(s)
    inv_sigmas.append(1/s)
    tmp = sigma[sigma['sigma']==s].reset_index()
    sigmas_beta.append((pow(s,2*float(tmp['beta']-1))))
    pos = np.array(tmp['[[positions]]'])[0]
    velocities.append(np.array(tmp['[[velocities]]'])[0][0])
    diameters.append(flock_diameter(pos))

df = pd.DataFrame({'sigma': sigmas,
                   'inv_sigma': inv_sigmas,
                   'sigma_beta': sigmas_beta,
                   'velocity': velocities,
                   'diameter': diameters})
df['vel_norm'] = df['velocity'].apply(np.linalg.norm)

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\sigma^{-1}$')
ax.set_ylabel('D')
ax.set(xscale="log")
g0 = plt.scatter(x=df['inv_sigma'], y=df['diameter'], marker='.')
g1 = plt.axvline(1, 0, 1,color='Red')
plt.legend([g0, g1], [r'$D(\sigma)$', r'$\sigma = 1.0$'])
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\sigma^{-1}$')
ax.set_ylabel(r'$\left\Vert\widehat{v}\right\Vert$')
ax.set(xscale="log")
g0 = plt.scatter(x=df['inv_sigma'], y=df['vel_norm'], marker='.')
g1 = plt.axvline(1, 0, 1,color='Red')
plt.legend([g0, g1], [r'$\left\Vert\widehat{v}\right\Vert(\sigma)$', r'$\sigma = 1.0$'])
plt.show()

del sigma
del diameters
del velocities
del vel_norms
del sigmas


#%% study on birds

path=proj_path+'results2/'
results = build_dataframe(path, True)
tmp = results[results['sigma']==1]
birds = tmp[tmp['beta']==0.4].reset_index()

diameters = []
velocities = []
vel_norms = []
Ns = []
for b in birds['N'].unique():
    Ns.append(b)
    tmp = birds[birds['N']==b].reset_index()
    pos = np.array(tmp['[[positions]]'])[0]
    velocities.append(np.array(tmp['[[velocities]]'])[0][0])
    diameters.append(flock_diameter(pos))

df = pd.DataFrame({'N': Ns,
                   'velocity': velocities,
                   'diameter': diameters})
df['vel_norm'] = df['velocity'].apply(np.linalg.norm)

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$k$')
ax.set_ylabel(r'$D$')
ax.relim()
g0 = plt.scatter(x=df['N'], y=df['diameter'], marker='.')
#sns.scatterplot(x='N', y='diameter', data=df)
plt.legend([g0], [r'$D(k)$'])
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$k$')
ax.set_ylabel(r'$\left\Vert\widehat{v}\right\Vert$')
ax.relim()
g0 = plt.scatter(x=df['N'], y=df['vel_norm'], marker='.')
#sns.scatterplot(x='N', y='vel_norm', data=df)
plt.legend([g0], [r'$\left\Vert\widehat{v}\right\Vert(k)$'])
plt.show()

#%% study on acceptance

path=proj_path+'results4/'
results = build_dataframe(path, True)
tmp = results[results['beta']==0.5]
acceptance = tmp[tmp['sigma']==1].reset_index()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$t_{last}$')
ax.set_ylabel(r'$\left\Vert d^{max}\right\Vert$')
ax.relim()
ax.set(yscale="log")
g0 = plt.scatter(x=acceptance['iter'], y=acceptance['space_norm'], marker='.')
#sns.scatterplot(x='iter', y='space_norm', data=acceptance)
plt.legend([g0], [r'$\left\Vert d^{max}\right\Vert(t_{last})$'])
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$t_{last}$')
ax.set_ylabel(r'$\left\Vert v^{max}\right\Vert$')
ax.relim()
ax.set(yscale="log")
g0 = plt.scatter(x=acceptance['iter'], y=acceptance['vel_norm'], marker='.')
#sns.scatterplot(x='iter', y='space_norm', data=acceptance)
plt.legend([g0], [r'$\left\Vert v^{max}\right\Vert(t_{last})$'])
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()

path=proj_path+'results5/'
results = build_dataframe(path, True)
tmp = results[results['beta']==0.5]
acceptance = tmp[tmp['sigma']==1].reset_index()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$k$')
ax.set_ylabel(r'$\left\Vert d^{max}\right\Vert$')
ax.relim()
#ax.set(yscale="log")
g0 = plt.scatter(x=acceptance['N'], y=acceptance['space_norm'], marker='.')
#sns.scatterplot(x='iter', y='space_norm', data=acceptance)
plt.legend([g0], [r'$\left\Vert d^{max}\right\Vert(k)$'])
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()

sns.set(style = "whitegrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$k$')
ax.set_ylabel(r'$\left\Vert v^{max}\right\Vert$')
ax.relim()
#ax.set(yscale="log")
g0 = plt.scatter(x=acceptance['N'], y=acceptance['vel_norm'], marker='.')
#sns.scatterplot(x='iter', y='space_norm', data=acceptance)
plt.legend([g0], [r'$\left\Vert v^{max}\right\Vert(k)$'])
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()


#%% moving birds animation
path=proj_path+'results2/'
results = build_dataframe(path)
tmp = results[results['sigma']==1]
birds = tmp[tmp['beta']==0.4].reset_index()

for k in birds['N'].unique():
    if k%10==0:
        print(k)
        tmp = birds[birds['N']==k]
        filename = proj_path+'/images/animation_N'+str(k)
        df = list(tmp['[[positions]]'])
        df = traslate_to_baricenter(df)
        positions_evolutions(filename, df, True)

#%%

# velocity convergence criterion
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20,10))
plt.xlabel('t')
ax.set_ylabel(r'$\epsilon_v$')
ax.set(yscale="log")
ax.relim()
sns.scatterplot(x='tau', y='vel_norm', hue='N',
                data=birds, legend='full')
plt.show()

# spatial convergence criterion
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20,10))
plt.xlabel('t')
ax.set_ylabel(r'$\epsilon_v$')
ax.set(yscale="log")
ax.relim()
sns.scatterplot(x='tau', y='space_norm', hue='N',
                data=birds, legend='full')
plt.show()

#%%
conv_1 = results[results['N']==10]

sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='space_norm', hue='beta', data=conv_1)
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
ax.relim()
sns.scatterplot(x='iter', y='vel_norm', hue='beta', data=conv_1)
plt.show()

#%% scale law
conv_1 = results[results['N']==10]
conv_1 = conv_1[conv_1['beta']==0.4]
conv_1 = conv_1[conv_1['sigma']!=1]


sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel('t')
ax.set_ylabel(r"$\epsilon_x$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='space_norm', hue='sigma', data=conv_1)
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel('t')
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='Vel_norm', hue='sigma', data=conv_1)

sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\tau$")
ax.set_ylabel('t')
ax.set(yscale="log")
sns.scatterplot(x='iter', y='Space_norm', hue='sigma', data=conv_1)
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\tau$")
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='Vel_norm', hue='sigma', data=conv_1)
plt.show()

sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\epsilon_x$")
ax.set_ylabel(r"$\epsilon_X$")
ax.set(yscale="log")
sns.scatterplot(x='space_norm', y='Space_norm', hue='sigma', data=conv_1)
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\epsilon_V$")
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='vel_norm', y='Vel_norm', hue='sigma', data=conv_1)
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel("t")
ax.set_ylabel(r"$\tau$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='tau', hue='sigma', data=conv_1)
plt.show()

sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\tau$")
ax.set_ylabel('t')
ax.set(yscale="log")
sns.scatterplot(x='tau', y='Space_norm', hue='sigma', data=conv_1)
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\tau$")
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='tau', y='Vel_norm', hue='sigma', data=conv_1)
plt.show()

sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\tau$")
ax.set_ylabel('t')
ax.set(yscale="log")
sns.scatterplot(x='iter', y='Space_norm', hue='sigma', data=conv_1)
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r"$\tau$")
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='Vel_norm', hue='sigma', data=conv_1)
plt.show()

#%%

conv_1 = conv_1[conv_1['sigma']==1]

sns.set(style = "darkgrid")
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='space_norm', hue='N', data=conv_1)
plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
plt.xlabel(r'$\beta$')
ax.set_ylabel(r"$\epsilon_v$")
ax.set(yscale="log")
sns.scatterplot(x='iter', y='vel_norm', hue='N', data=conv_1)
plt.show()


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
conv1['modX'] = conv1['modX']/np.sqrt(conv1['sigma'])
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