#input parameters
intensity       =       50e10
bunchLength_ns  =       271 #4 sigma in ns
#bunchLength_m	=	10.67 #1 sigma in m
emittance_x     =       5.4e-6
emittance_y     =       3.5e-6
dpp_rms         =       1.31e-3
bF              =       0.174 # None for Gaussian
Qh		=	4.08    # to estimate detuning and plot
Qv		=	4.13    # to estimate detuning and plot
plot_range  =   [[3.9,4.2],[3.9,4.2]]   # range in Qh & Qv for the plot
plot_order  =   4   # order of resonances to plot
periodicity =   16  # periodicity of ring for the colorcode of the plot
figure		=	'figure.png'
twiss_file  =   'twiss_files/twiss_PSB'
