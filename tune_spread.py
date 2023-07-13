import pylab as plt
import numpy as np
from PySCRDT.PySCRDT import PySCRDT
from matplotlib import colors
import resonance_lines
import input_parameters 
import matplotlib


s = PySCRDT()

keys=dir(input_parameters)

if 'bunchLength_ns' in keys:
    bunchLength=(input_parameters.bunchLength_ns/4.)*0.3
elif 'bunchLength_m' in keys:
    bunchLength = input_parameters.bunchLength_m
else:
    bunchLength = None

if 'bF' in keys:
    bF=input_parameters.bF
else:
    bF = None


s.setParameters(
    intensity=input_parameters.intensity,
    bunchLength=bunchLength,
    emittance_x=input_parameters.emittance_x,
    emittance_y=input_parameters.emittance_y, 
    dpp_rms=input_parameters.dpp_rms,  
    bF=bF
)

s.prepareData(twissFile=input_parameters.twiss_file)

if ('b' in keys) and ('g' in keys):
    s.updateParameters(b=input_parameters.b)
    s.updateParameters(g=input_parameters.g)
elif 'g' in keys:
    s.updateParameters(g=input_parameters.g)
    s.updateParameters(b=np.sqrt(1-1/input_parameters.g**2))
elif 'b' in keys:
    s.updateParameters(b=input_parameters.b)
    s.updateParameters(g=1/np.sqrt(1-input_parameters.b**2))

if 'bunchLength_ns' in keys:
    s.updateParameters(bunchLength=bunchLength*s.parameters['b'])

if 'ro' in keys:
    s.updateParameters(ro=input_parameters.ro)
if 'harmonic' in keys:
    s.updateParameters(harmonic=input_parameters.harmonic)
# caclulate detuning coefficients using potentials up to 20th order (needed for up to 3 sigma particles)

detuning=[]
# order in x
for i in range(0,21,2):
    # order in y
    for j in range(0,21,2): 
        if (i==0) and (j==0):
            pass
        elif i+j<21:
            s.setOrder([int(i),int(j),'any'])
            s.potential()
            s.detuning()
            detuning.append([i,j,s.getDetuning()])
    
detuning=np.array(detuning)

#  initialize grid for calculation
 
s_N=6
s_max=3
theta_N=5

def initial_xy_polar(s_max, s_N, theta_N):
    return np.array(
        [
            [(s*np.cos(theta), s*np.sin(theta)) for s in np.linspace(0, s_max, s_N+1)]
            for theta in np.linspace(0, np.pi/2., theta_N)
        ])

S = initial_xy_polar(s_max=s_max, s_N=s_N, theta_N=theta_N)

#  estimate tunes from the detuning coefficients

en_x=s.parameters['emittance_x']
en_y=s.parameters['emittance_y']
beta=s.parameters['b']
gamma=s.parameters['g']
J_x=S[:,:,0]**2*en_x/2./beta/gamma
J_y=S[:,:,1]**2*en_y/2./beta/gamma

Qx,Qy=input_parameters.Qh,input_parameters.Qv 

for x_q,y_q,detuning_coef in detuning:
    if x_q:
        Qx+=x_q/2.*detuning_coef*(J_x**(x_q/2.-1))*(J_y**(y_q/2.))
    if y_q:
        Qy+=y_q/2.*detuning_coef*(J_y**(y_q/2.-1))*(J_x**(x_q/2.))

Q = np.dstack(
    (
        [qx.tolist() + [input_parameters.Qh] for qx in Qx],
        [qy.tolist() + [input_parameters.Qv] for qy in Qy],
    )
)

Q[:,:,0] += 0.00
Q[:,:,1] += 0.00

sigma_x = Q.shape[0]-1
sigma_y = Q.shape[1]-1
p1 = Q[:-1, :-1, :].reshape(sigma_x*sigma_y, 2)[:, :]
p2 = Q[1:, :-1, :].reshape(sigma_x*sigma_y, 2)[:]
p3 = Q[1:, 1:, :].reshape(sigma_x*sigma_y, 2)[:]
p4 = Q[:-1, 1:, :].reshape(sigma_x*sigma_y, 2)[:]


# do the plotting

cmap_base = plt.cm.hot
c_indcs = np.int_(np.linspace(0.1,0.6,s_N+1)*cmap_base.N)
cmap = colors.ListedColormap([cmap_base(c_indx) for c_indx in c_indcs])

# Stack endpoints to form polygons
Polygons = np.transpose(np.stack((p1, p2, p3, p4)), (1, 0, 2))
patches = list(map(matplotlib.patches.Polygon, Polygons))
p_collection = matplotlib.collections.PatchCollection(
#     patches, edgecolor='grey', linewidth=1,
    patches, edgecolor='k', linewidth=0.5,
#     facecolors=[],
    facecolors=cmap.colors,
#     facecolors=['SkyBlue'],
    alpha=0.7
)

detuning_y=[i for i in np.where(detuning[:,1]==2)[0] if i in np.where(detuning[:,0]==0)[0]][0]
detuning_x=[i for i in np.where(detuning[:,0]==2)[0] if i in np.where(detuning[:,1]==0)[0]][0]
print('dQx = ', str(detuning[detuning_x,2]))
print('dQy = ', str(detuning[detuning_y,2]))

fig,ax = plt.subplots(1,figsize=(4.5,4.5))
tune_diagram = resonance_lines.resonance_lines(input_parameters.plot_range[0],
            input_parameters.plot_range[1], np.arange(1,input_parameters.plot_order+1), input_parameters.periodicity)
tune_diagram.plot_resonance(figure_object=fig)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.add_collection(p_collection)
ax.set_aspect('equal')
plt.savefig(input_parameters.figure)
plt.show()