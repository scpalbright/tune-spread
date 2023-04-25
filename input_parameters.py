#input parameters
intensity       =      25e10
#bunchLength_ns =       150 #4 sigma in ns
#bunchLength_m	=	    10.67 #1 sigma in m
emittance_x     =       2.1e-6
emittance_y     =       1.7e-6
dpp_rms         =       1.85e-3
#ro             =       # classical radius : default proton change for ions accordingly
bF              =       0.24 # None for Gaussian
Qh		=	4.17    # to estimate detuning and plot
Qv		=	4.23    # to estimate detuning and plot
plot_range  =   [[3.9,4.3],[3.9,4.3]]   # range in Qh & Qv for the plot
plot_order  =   4   # order of resonances to plot
periodicity =   16  # periodicity of ring for the colorcode of the plot
figure		=	'figure.png'
twiss_file  =   'twiss_files/twiss_PSB'
#b  =   # relativistic beta  :  if desired energy different than in twiss file
#g  =   # relativistic gamma :  if desired energy different than in twiss file
