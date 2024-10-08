import pylab as plt
import numpy as np
from PySCRDT.PySCRDT import PySCRDT
from matplotlib import colors
import resonance_lines
import input_parameters
import matplotlib

s = PySCRDT()

keys = dir(input_parameters)

if 'bunch_length_ns' in keys:
    bunch_length = (input_parameters.bunch_length_ns/4.) * 0.3
elif 'bunch_length_m' in keys:
    bunch_length = input_parameters.bunch_length_m
else:
    bunch_length = None

if 'bF' in keys:
    bF = input_parameters.bF
else:
    bF = None

s.intensity = input_parameters.intensity
s.bunch_length = bunch_length
s.emittance_x = input_parameters.emittance_x
s.emittance_y = input_parameters.emittance_y
s.dpp_rms = input_parameters.dpp_rms
s.bunching_factor = bF

s.prepare_data(twiss_file = input_parameters.twiss_file)

if ('b' in keys) and ('g' in keys):
    s.b = input_parameters.b
    s.g = input_parameters.g

elif 'g' in keys:
    s.g = input_parameters.g
    s.b = np.sqrt(1 - 1/input_parameters.g**2)

elif 'b' in keys:
    s.b = input_parameters.b
    s.g = 1/np.sqrt(1 - input_parameters.b**2)

if 'bunchLength_ns' in keys:
    s.bunch_length = bunch_length*s.b

if 'ro' in keys:
    s.ro = input_parameters.ro

if 'harmonic' in keys:
    s.harmonic = input_parameters.harmonic
# caclulate detuning coefficients using potentials up to 20th order (needed for up to 3 sigma particles)

dQx = s.calc_detuning(2, 0, 'any')
dQy = s.calc_detuning(0, 2, 'any')

print(f"Direct dQx = {dQx}")
print(f"Direct dQy = {dQy}")
