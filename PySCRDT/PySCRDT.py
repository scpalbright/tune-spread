"""
#
#   PySCRDT
#   A module to calculate the resonance driving terms from the space
#   charge potential
#
#   Version :   1.1
#               pre-calculated potentials for faster evaluation
#   Author  : F. Asvesta
#   Contact : fasvesta .at. cern .dot. ch
#
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals, annotations)

import dataclasses as dc
import numpy as np
import sympy as sy
import dill

from icecream import ic
import inspect
import sys
import time

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal, Optional, List, Tuple, Dict, Callable
    from sympy import Expr


__version   = 1.1
__PyVersion = ["2.7"]
__author    = ["Foteini Asvesta"]
__contact   = ["fasvesta .at. cern .dot. ch"]

@dc.dataclass
class SCParameters:

    intensity: float = 41E10
    bunch_length: float = 5.96
    ro: float = 1.5347e-18
    emittance_x: float = 2E-6
    emittance_y: float = 1.1E-6
    dpp_rms: float = 0.5E-3
    dpp: float = 0
    bF: Optional[float] = None
    harmonic: int = 1
    b: Optional[float] = None
    g: Optional[float] = None
    C: Optional[float] = None


class PySCRDT(object):
    """
    class for the calculation of the rdts
    Returns: PySCRDT instance
    """

    def __init__(self, parameters: bool | str = False,
                 mode: Optional[Literal[3] | Literal[5]] = None,
                 twiss_file: Optional[str] = None, order: Optional[List] = None,
                 twiss_table_xsuite: Optional[str] = None):
        """
        Initialization function
        Input :  parameters : [bool|str]  Parameters needed for the
                                          calculations (default=False)
                                          if True the default values of
                                          [set_parameters] are used if
                                          str parameters are read in a
                                          file
                 mode       : [3|5]       Resonance description mode
                                          (default=None)
                                          If None, len(order) is used.
                 order      : [list]      Resonance order and harmonic
                                          (default=None)
#TODO:                                    If None, TO BE COMFIRMED
                 twiss_file  : [str]      MAD-X Twiss file
                                          (default=None)
#TODO:                                    If None, TO BE CONFIRMED
        """

        self.x, self.y, self.t= sy.symbols('x y t')

        self.a = sy.Symbol('a', positive=True, real=True)
        self.b = sy.Symbol('b', positive=True, real=True)
        self.D = sy.Symbol('D', positive=True, real=True)

        self.fx = sy.Symbol('fx', positive=True, real=True)
        self.fy = sy.Symbol('fy', positive=True, real=True)

        self.V = None
        self.K = None
        self.data = None
        self.factor = None
        self.factor_d = None
        self.rdt = None
        self.rdt_d = None
        self.feed = False
        self._mode = None
        self._order = None

        self._potential_functions: Dict[Tuple[int, int], Callable] = {}
        try:
            self.load_potential_functions()
        except RuntimeError:
            print("Pre-calculated potential functions not available")

        if type(parameters) is str:
            self.parameters=None
            self.read_parameters(parameters)

        else:

            if parameters:
                self.set_parameters()
            else:
                self.parameters=None
                print("# PySCRDT : Set parameters with function [set_parameters]"
                      +" or read parameters with [read_parameters]")

        if twiss_file is None:
            if twiss_table_xsuite is None:
                print("# PySCRDT : Import Madx twiss file using the function "
                      + "[prepare_data] or X-suite Twiss from method "
                      + "[twiss_table_xsuite]")
            else:
                self.load_twiss_from_xsuite(twiss_table_xsuite)

        else:
            self.prepare_data(twiss_file)

        if (order is None) and (mode is None):
            print("# PySCRDT : Set order in [set_order]")
        elif order is None:
            print("# PySCRDT : Set order in [set_order]")
        else:
            if mode is None:
                self.set_mode(len(order))
            else:
                self.set_mode(mode)
            self.set_order(order)

    @property
    def mode(self) -> int:
        return self._mode

    @mode.setter
    def mode(self, value: Optional[Literal[3] | Literal[5]]):
        if value not in (3, 5):
            raise ValueError("# PySCRDT::set_mode: Mode needs to be 3 or 5 "
                          + "depending on the desired resonance description")
        else:
            self._mode = value


    def set_mode(self, mode: Optional[Literal[3] | Literal[5]] = None):
        """
        Sets the mode for the characterization of resonances
        Input :  mode  : [3|5]
        Returns: void
        """

        if mode is None:
            print("# PySCRDT : Set mode in [set_mode]")
        else:
            self.mode = mode


    @property
    def order(self) -> Tuple[int, int, int] | Tuple[int, int, int, int, int]:

        if self._mode == 3:
            return (self.m, self.n, self.l)
        elif self._mode == 5:
            return (self.h, self.i, self.j, self.k, self.l)
        elif self._mode is None:
            return self._mode
        else:
            raise RuntimeError("Order not available, requires mode == 3 or"
                               + f" mode==5, but mode=={self._mode}")

    @order.setter
    def order(self, value: Tuple[int | str, ...]):

        if self._mode is None:
            self.mode = len(value)

        elif self._mode != len(value):
            raise ValueError("len(order) must match mode.  Mode is "
                             + f"{self._mode}, len(order) is {len(value)}")


        for v in value:
            if not isinstance(v, (int, str)):
                raise TypeError("All values of order must be of type int")

        self.h = None
        self.i = None
        self.j = None
        self.k = None
        self.l = None
        self.m = None
        self.n = None

        if len(value) == 3:
            self.m, self.n, self.l = value
        elif len(value) == 5:
            self.h, self.i, self.j, self.k, self.l = value


    def set_order(self, args):
        """
        Sets the Resonance Orders
        In "3 mode"
        Input :  m  : [int] order in H
                 n  : [int] order in V
                 l  : [int|'any'] harmonic of the resonance
        In "5 mode"
        Input :  h  : [int] characteristic of H
                 i  : [int] characteristic of H
                 j  : [int] characteristic of V
                 k  : [int] characteristic of V
                 l  : [int|'any'] harmonic of the resonance
        Returns: void
        """

        self.order = args
        self.factor = None
        self.factor_d = None


    @property
    def parameters(self) -> Dict[str, float | int]:
        return dc.asdict(self._parameters)

    @parameters.setter
    def parameters(self, value: Dict[str, float | int] | SCParameters):

        if value is None:
            self._parameters = value
        elif isinstance(value, SCParameters):
            self._parameters = value
        else:
            self._parameters = SCParameters(**value)


    def set_parameters(self, intensity: float = 41e10, bunch_length: float = 5.96,
                      ro: float = 1.5347e-18, emittance_x: float = 2e-6,
                      emittance_y: float = 1.1e-6, dpp_rms: float = 0.5e-3,
                      dpp: float = 0.0, bF: Optional[float] = None,
                      harmonic: int = 1):
        """
        Sets the parameters for the calculation:
        Input :  intensity  : [float] bunch intensity in ppb
                              (Default=41e10)
                 bunch_length: [float] RMS bunch length in m
                              (Default=5.96)
                 ro         : [float] classical particle radious in m
                              (Default=1.5347e-18 {proton})
                 emittance_x: [float] normalized horizontal emittance in
                                      m*rad
                              (Default=2e-6)
                 emittance_y: [float] normalized vertical emittance in
                                      m*rad
                              (Default=1.1e-6)
                 dpp_rms    : [float] RMS Dp/p
                              (Default=0.5e-3)
                 dpp        : [float] Single particle Dp/p
                              (Default=0)
                 bF         : [float] Bunching factor
                              (Default=None)
#TODO:                        If none, TO BE CONFIRMED
                 harmonic   : [int]   Harmonic number, # of buckets
                              (Default=1)
        """

        self.parameters={'intensity':intensity, 'bunch_length':bunch_length,
                         'ro':ro, 'emittance_x':emittance_x,
                         'emittance_y':emittance_y, 'dpp_rms':dpp_rms,
                         'dpp':dpp, 'bF':bF, 'harmonic':harmonic}


    def read_parameters(self, input_file: str):
        # TODO: Understand input
        """
        Reads the parameters from an input file:
        Input Format: intensity   = [float] # bunch intensity in ppb (Default=41e10)
                      bunch_length = [float] # RMS bunch length in m  (Default=5.96)
                      ro          = [float] # classical particle radious in m (Default=1.5347e-18 {proton})
                      emittance_x = [float] # normalized horizontal emittance in m*rad (Default=2e-6)
                      emittance_y = [float] # normalized vertical emittance in m*rad (Default=1.1e-6)
                      dpp_rms     = [float] # RMS Dp/p (Default=0.5e-3)
                      dpp         = [float] # Single particle Dp/p (Default=0.5e-3)
                      bF          = [float] # Bunching factor (Default=None)
                      harmonic    = [int]   # Harmonic number,# of buckets (Default=1)
        Returns: void
        """

        params = np.genfromtxt(input_file, dtype=str)

        if self.parameters is None:
            self.set_parameters()

        if len(np.shape(params)) == 1:
            if params[0] not in self.parameters.keys():
                raise IOError(f"# PySCRDT::read_parameters: {params[0]}"
                              + " not recognized [check_writing]")
            else:
                self.parameters[params[0]]=float(params[2])

        else:
            for i in enumerate(params):
                if params[i[0]][0] not in self.parameters.keys():
                    raise IOError("# PySCRDT::read_parameters: "
                                   + f"{params[i[0]][0]} not recognized "
                                   + "[check_writing]")
                else:
                    self.parameters[params[i[0]][0]]=float(params[i[0]][2])

        if self.data is not None:
            self.beam_size()
            self.ksc()


    def load_potential_functions(self):
        """
        Attempt to load the precalculated potential functions

        Raises:
            RuntimeError: If the functions cannot be loaded, a
                          RuntimeError is raised.

        Returns:

        """
        ic()
        try:
            with open(__file__[:__file__.find('PySCRDT.py')]
                        +'potentialsPy3','rb') as f:
                pots = dill.load(f)

        except:
            try:
                with open(__file__[:__file__.find('PySCRDT.py')]
                          + 'potentialsPy2','rb') as f:
                    pots = dill.load(f)
            except:
                raise RuntimeError("Cannot load potentials")

        for p in pots:
            self._potential_functions[(p[0], p[1])] = p[2]


    def potential(self, feed_down: bool = False, look_up=True):
        #TODO: Check look_up
        """
        Calculates the space charge potential for the given resonance
        order
        Inputs : feedDown : [bool] needed when single particle Dp/p non 0
                            (default=False)
        """
        if self.mode == 5:
            self.m = self.h+self.i
            self.n = self.j+self.k

        if (self.m is None) or (self.n is None):
            raise IOError("# PySCRDT::potential: You need to define "
                          + "resonance order in [set_order]")

        if self.m%2 != 0 and (feed_down == False):
            raise IOError("# PySCRDT::potential: Space charge potential"
                          + " contains only even orders without Dp/p (order "
                          + f"given {self.m}), change the order in [set_order]")

        if self.n%2 != 0:
            raise IOError("# PySCRDT::potential: Space charge potential "
                          + f"contains only even orders (order given {self.m})"
                          + ", change the order in [set_order]")

        if (self.m+self.n < 21) and (feed_down == False) and (look_up == True):
            if (self.m == 32 and self.n == 8) or (self.m == 8 and self.n == 32):
                pass
            else:
                try:
                    self.f = self._potential_functions[(self.m, self.n)]
                except KeyError:
                    look_up = False

        # TODO: Make new pre-calculator
        if (self.m+self.n > 21) or (feed_down == True) or (look_up == False):
            V = ((-1 + sy.exp(-self.x**2 / (self.t + 2*self.a**2)
                              -self.y**2 / (self.t + 2*self.b**2)))
                  / sy.sqrt((self.t + 2*self.a**2)
                           *(self.t + 2*self.b**2)))

            if self.m>self.n:
                if feed_down:
                    p1 = sy.series(V, self.x, 0, abs(self.m)+2).removeO()
                else:
                    p1 = sy.series(V, self.x, 0, abs(self.m)+1).removeO()

                p2 = sy.series(p1, self.y, 0, abs(self.n)+1).removeO()
                termy = sy.collect(p2, self.y, evaluate=False)
                termpowy = termy[self.y**abs(self.n)]

                if feed_down:
                    termpowy = sy.expand(termpowy.subs(self.x, self.x+self.D))

                termx = sy.collect(termpowy, self.x, evaluate = False)
                termpowx = termx[self.x**abs(self.m)]
                sterm = sy.simplify(termpowx)

            else:
                p1 = sy.series(V, self.y, 0, abs(self.n)+1).removeO()

                if feed_down:
                    p2 = sy.series(p1, self.x, 0, abs(self.m)+2).removeO()
                else:
                    p2 = sy.series(p1, self.x, 0, abs(self.m)+1).removeO()

                termx = sy.collect(p2, self.x, evaluate=False)

                if feed_down:
                    termx = sy.expand(termx.subs(self.x, self.x+self.D))

                termpowx = termx[self.x**abs(self.m)]
                termy = sy.collect(termpowx, self.y, evaluate=False)
                termpowy = termy[self.y**abs(self.n)]
                sterm = sy.simplify(termpowy)

            res = sy.integrate(sterm, (self.t, 0, sy.oo)).doit()
            result = res.doit()
            self.V = sy.simplify(result)
            self.f = ic(sy.lambdify((self.a, self.b, self.D), self.V))


    def ksc(self):
        """
        Calculates the space charge perveance Ksc from the parameters
        dictionary
        """

        if self.parameters is None:
            raise IOError("# PySCRDT::ksc: You need to define "
                          + "parameters in [set_parameters]")

        if self.data is None:
            raise IOError("# PySCRDT::ksc: You need to define Madx "
                          + "twiss file in [prepare_data]")

        params = self._parameters

        if params.bF:
            self.K = (2*params.intensity*params.ro
                      *(params.harmonic/params.bF)
                      /(params.C * params.b**2 * params.g**3))

        else:
            self.K = (2*params.intensity*params.ro
                      /(np.sqrt(2*np.pi)*params.bunch_length
                        * params.b**2 * params.g**3))


    def beam_size(self):
        """
        Calculates the transverse beam sizes from the parameters
        dictionary and the twiss file
        """

        if self.parameters is None:
            raise IOError("# PySCRDT::beamSize: You need to define "
                          + "parameters in [set_parameters]")

        if self.data is None:
            raise IOError("# PySCRDT::ksc: You need to define Madx "
                          + "twiss file in [prepare_data]")

        params = self._parameters

        self.sx = np.sqrt(params.emittance_x*self.data[:,1]
                          /(params.b*params.g)
                          +(params.dpp_rms*self.data[:,3])**2)
        self.sy = np.sqrt(params.emittance_y*self.data[:,2]
                          /(params.b*params.g)
                          +(params.dpp_rms*self.data[:,4])**2)


    def prepare_data(self, twiss_file: Optional[str] = None):
        """
        Prepares the data from a MADX Twiss file including at least
        {s, betx, bety, dx, dy, mux, muy, l}
        Inputs : twiss_file : [str] twiss file (default=None)
        """

        if twiss_file is None:
            raise IOError("# PySCRDT::prepare_data: You need to define "
                          + "Madx twiss file in [prepare_data]")
        if self.parameters is None:
            raise IOError("# PySCRDT::prepare_data: You need to define "
                          + "parameters in [set_parameters]")

        with open(twiss_file, 'r') as f:
            for line in enumerate(f.readlines()):
                if line[1][0]=='*':
                    skip_header_nr = line[0]
                elif line[1][0]=='$':
                    skip_rows_nr = line[0]+1
                    break

        params = np.genfromtxt(twiss_file, max_rows = 40, dtype = str)

        for i in enumerate(params):
            if params[i[0]][1] == 'GAMMA':
                self._parameters.g = float(params[i[0]][3])
                self._parameters.b = np.sqrt(1 - 1/self._parameters.g**2)

            elif params[i[0]][1] == 'LENGTH':
                self._parameters.C = float(params[i[0]][3])

            elif params[i[0]][1] == 'Q1':
                self.actualQx = float(params[i[0]][3])

            elif params[i[0]][1] == 'Q2':
                self.actualQy = float(params[i[0]][3])

        # 45 originally, below 47
        header = np.genfromtxt(twiss_file, skip_header = skip_header_nr,
                               max_rows = 1, dtype = str)

        data = np.loadtxt(twiss_file, skiprows = skip_rows_nr,
                          usecols = (np.where(header == 'S')[0][0] - 1,
                                     np.where(header == 'BETX')[0][0] - 1,
                                     np.where(header == 'BETY')[0][0] - 1,
                                     np.where(header == 'DX')[0][0] - 1,
                                     np.where(header == 'DY')[0][0] - 1,
                                     np.where(header == 'MUX')[0][0] - 1,
                                     np.where(header == 'MUY')[0][0] - 1,
                                     np.where(header == 'L')[0][0] - 1))

        s = np.linspace(0, self._parameters.C, 100000)

        data2 = np.zeros((100000, 8))

        data2[:,1] = np.square(np.interp(s, data[:,0], np.sqrt(data[:,1])))
        data2[:,2] = np.square(np.interp(s, data[:,0], np.sqrt(data[:,2])))
        data2[:,3] = np.interp(s,data[:,0], self._parameters.b*data[:,3])
        data2[:,4] = np.interp(s,data[:,0], self._parameters.b*data[:,4])
        data2[:,5] = np.interp(s,data[:,0], data[:,5])
        data2[:,6] = np.interp(s,data[:,0], data[:,6])
        data2[:,7] += self._parameters.C/len(s)
        data2[:,0] = s
        self.data = data2

        self.beam_size()
        self.ksc()


    def load_twiss_from_xsuite(self, twiss_table_xsuite: Optional[dict] = None):
        """
        Instead of using a MADX Twiss file, load Twiss data generated
        from X-suite tracker including at least
        {s, betx, bety, dx, dy, mux, muy, l}
        Inputs : twiss_table_xsuite : [dict] twiss table (default=None)

        Just like prepare_data method, increases resolution of Twiss table by interpolation
        """

        if twiss_table_xsuite is None:
            raise IOError("# PySCRDT::load_twiss_from_xsuite: You need to "
                          + "define Xsuite twiss table in "
                          + "[load_twiss_from_xsuite]")

        if self.parameters is None:
            raise IOError("# PySCRDT::load_twiss_from_xsuite: You need to "
                          + "define parameters in [set_parameters]")

        # Define parameters from Twiss table
        self._parameters.g = twiss_table_xsuite['particle_on_co'].gamma0[0]
        self._parameters.b = twiss_table_xsuite['particle_on_co'].beta0[0]
        self._parameters.C = twiss_table_xsuite['circumference']
        self.actualQx = twiss_table_xsuite['qx']
        self.actualQy = twiss_table_xsuite['qy']

        # Set up data for increased resolution by interpolation
        #columns = ['s', 'betx', 'bety', 'dx', 'dy', 'mux', 'muy', 'l']
        s = np.linspace(0, self._parameters.C, 100000)

        data2 = np.zeros((100000, 8))

        data2[:,1] = np.square(np.interp(s, twiss_table_xsuite['s'],
                                         np.sqrt(twiss_table_xsuite['betx'])))
        data2[:,2] = np.square(np.interp(s, twiss_table_xsuite['s'],
                                         np.sqrt(twiss_table_xsuite['bety'])))
        data2[:,3] = np.interp(s, twiss_table_xsuite['s'],
                               self._parameters.b*twiss_table_xsuite['dx'])
        data2[:,4] = np.interp(s, twiss_table_xsuite['s'],
                               self._parameters.b*twiss_table_xsuite['dy'])
        data2[:,5] = np.interp(s, twiss_table_xsuite['s'],
                               twiss_table_xsuite['mux'])
        data2[:,6] = np.interp(s, twiss_table_xsuite['s'],
                               twiss_table_xsuite['muy'])
        data2[:,7] += self._parameters.C/len(s)
        data2[:,0] = s

        self.data=data2
        self.beam_size()
        self.ksc()


    def re_indexing(self, factor: dict) -> dict:
        """
        Auxiliary Function
        Returns: [dict]
        """

        self.dictionary={}

        for i in factor.keys():
            #if len(i.args)==0:
            self.dictionary[i]=factor[i]
            #else:
            #    self.dictionary[sy.exp(i.args[0]/1.)]=factor[i]

        return self.dictionary


    def calculate_factor(self, detuning: bool = False) -> float:
        """
        Auxiliary Function
        Returns: [float]
        """

        if detuning == True:
            if self.m == 0:
                det1 = sy.cos(self.fy)**abs(self.n)
                det3 = det1.rewrite(sy.exp)
                det2 = sy.expand(det3)
                self.factor_d = float(sy.collect(det2, sy.exp(1j*self.fy),
                                                 evaluate=False)[1])

            elif self.n == 0:
                det1 = sy.cos(self.fx)**abs(self.m)
                det3 = det1.rewrite(sy.exp)
                det2 = sy.expand(det3)
                self.factor_d = float(sy.collect(det2, sy.exp(1j*self.fx),
                                                 evaluate=False)[1])

            else:
                det1 = (sy.cos(self.fx)**abs(self.m)
                        *sy.cos(self.fy)**abs(self.n))

                det3 = det1.rewrite(sy.exp)
                det2 = sy.expand(det3)
                factor1 = sy.collect(det2, sy.exp(1j*self.fx),
                                     evaluate=False)[1]
                self.factor_d = float(sy.collect(factor1, sy.exp(1j*self.fy),
                                               evaluate=False)[1])

            return self.factor_d

        else:
            if self.mode == 5:
                self.m = self.h-self.i
                self.n = self.j-self.k

            if self.m == 0:
                det1 = sy.cos(self.fy)**abs(self.n)
                det3 = det1.rewrite(sy.exp)
                det2 = sy.expand(det3)
                factor = sy.collect(det2, sy.exp(1j*self.fy), evaluate=False)
                dictionary = self.re_indexing(factor)
                self.factor = float(2. * dictionary[sy.exp(1j*self.fy)
                                                    **float(abs(self.n))])

            elif self.n == 0:
                det1 = sy.cos(self.fx)**abs(self.m)
                det3 = det1.rewrite(sy.exp)
                det2 = sy.expand(det3)
                factor = sy.collect(det2, sy.exp(1j*self.fx), evaluate=False)
                dictionary = self.re_indexing(factor)
                self.factor = float(2. * dictionary[sy.exp(1j*self.fx)
                                                    **float(abs(self.m))])

            else:
                det1 = (sy.cos(self.fx)**abs(self.m)
                        *sy.cos(self.fy)**abs(self.n))
                det3 = det1.rewrite(sy.exp)
                det2 = sy.expand(det3)
                factor1 = sy.collect(det2, sy.exp(1j*self.fx), evaluate=False)
                dictionary = self.re_indexing(factor1)
                factor1 = dictionary[sy.exp(1j*self.fx)**float(abs(self.m))]
                factor2 = sy.collect(factor1, sy.exp(1j*self.fy),
                                     evaluate=False)
                dictionary = self.re_indexing(factor2)
                self.factor = float(2. * dictionary[sy.exp(1j*self.fy)
                                                    **float(abs(self.n))])

            return self.factor


    def resonance_driving_terms(self, feed_down: bool = False):
        """
        Calculates the resonance driving terms for the resonance requested
        Returns: Void
        """

        self.feed = feed_down

        if self.V is None:
            self.potential(feed_down = self.feed)

        if self.data is None:
            raise IOError("# PySCRDT::resonance_driving_terms: You need to run"
                          +" [prepare_data] first")

        if self.K is None:
            self.ksc()

        if self.factor is None:
            self.calculate_factor()

        if self.mode == 3:
            if self.l == 'any':
                self.rdt_s = (self.factor*self.data[:,7]*self.K/2./(2*np.pi)
                              * (np.sqrt(2*self.data[:,1]) ** abs(self.m))
                              * (np.sqrt(2*self.data[:,2]) ** abs(self.n))
                              * self.f(self.sx, self.sy,
                                       self.parameters['dpp']*self.data[:,3])
                                       * np.exp(1j*2*np.pi
                                                *(self.m*self.data[:,5]
                                                  + self.n*self.data[:,6])))
            else:
                self.rdt_s = (self.data[:,7]*self.factor*self.K/2./(2*np.pi)
                              * (np.sqrt(2*self.data[:,1]) ** abs(self.m))
                              * (np.sqrt(2*self.data[:,2]) ** abs(self.n))
                              * self.f(self.sx, self.sy,
                                       self.parameters['dpp']*self.data[:,3])
                                       * np.exp(
                                           1j*(self.m*2*np.pi*self.data[:,5]
                                               + self.n*2*np.pi*self.data[:,6]
                                               + (self.l-self.m * self.actualQx
                                                  - self.n * self.actualQy)
                                                  * 2*np.pi * self.data[:,0]
                                                  / self.parameters['C'])))

        else:
            if self.l == 'any':
                self.rdt_s = (self.factor*self.data[:,7]*self.K/2./(2*np.pi)
                              * (np.sqrt(2*self.data[:,1]) ** abs(self.h
                                                                  +self.i))
                              * (np.sqrt(2*self.data[:,2]) ** abs(self.j
                                                                  +self.k))
                                 * self.f(self.sx, self.sy,
                                          self.parameters['dpp']*self.data[:,3])
                                          *np.exp(1j*2*np.pi * ((self.h-self.i)
                                                  * self.data[:,5]
                                                  + (self.j-self.k)
                                                  * self.data[:,6])))
            else:
                self.rdt_s = (self.data[:,7]*self.factor*self.K/2./(2*np.pi)
                              * (np.sqrt(2*self.data[:,1]) ** abs(self.h
                                                                  +self.i))
                              * (np.sqrt(2*self.data[:,2]) ** abs(self.j
                                                                  +self.k))
                              * self.f(self.sx, self.sy,
                                       self.parameters['dpp']*self.data[:,3])
                                       * np.exp(1j*((self.h-self.i)*2*np.pi
                                                    * self.data[:,5]
                                                    + (self.j-self.k)*2*np.pi
                                                    * self.data[:,6]
                                                    + (self.l - (self.h-self.i)
                                                       * self.actualQx
                                                       - (self.j - self.k)
                                                       * self.actualQy)
                                                    * 2*np.pi*self.data[:,0]
                                                      /self.parameters['C'])))

        self.rdt = sum(self.rdt_s)


    def detuning(self):
        """
        Calculates the non linear detuning terms
        Returns: Void
        """

        if self.V is None:
            self.potential()

        if self.data is None:
            raise IOError("# PySCRDT::resonance_driving_terms: You need "
                          + "to run [prepare_data] first")

        if self.K is None:
            self.ksc()

        if self.factor_d is None:
            self.calculate_factor(detuning = True)

        if self.mode == 3:
            self.rdt_s_d = (self.factor_d*self.data[:,7]*self.K/2./(2*np.pi)
                            * (np.sqrt(2*self.data[:,1]) ** abs(self.m))
                            * (np.sqrt(2*self.data[:,2]) ** abs(self.n))
                            * self.f(self.sx, self.sy,
                                     self.parameters['dpp']*self.data[:,3]))

        else:
            self.rdt_s_d = (self.factor_d*self.data[:,7]*self.K/2./(2*np.pi)
                            * (np.sqrt(2*self.data[:,1]) ** abs(self.h+self.i))
                            * (np.sqrt(2*self.data[:,2]) ** abs(self.j+self.k))
                            * self.f(self.sx, self.sy,
                                     self.parameters['dpp']*self.data[:,3]))

        self.rdt_d = sum(self.rdt_s_d)


    def update_parameters(self, **kwargs):
        """
        Updates the parameter dictionary
        Input :  any of  'intensity'
                         'bunchLength'
                         'ro'
                         'emittance_x'
                         'emittance_y'
                         'dpp_rms'
                         'dpp'
                         'b'
                         'g'
                         'bF'
                         'harmonic'
        """

        if kwargs is not None:
            for key, value in kwargs.items():
                if key not in self.parameters.keys():
                    raise IOError(f"# PySCRDT::update_parameters: {key} "
                                  + "not recognized [check_writing]")
                else:
                    setattr(self._parameters, key, value)

            if self.data is not None:
                self.beam_size()
                self.ksc()


    def get_parameters(self) -> dict:
        """
        Returns the parameters dictionary
        Returns: [dict] the parameters dictionary
        """

        if self.parameters is None:
            raise IOError("# PySCRDT::get_parameters: You need to define"
                          + "parameters in [set_parameters]|[read_parameters]")
        else:
            return self.parameters


    def get_working_point(self) -> Tuple[float, float]:
        """
        Returns the tunes (Qx, Qy)
        Returns: [tuple]
        """

        if self.data is None:
            raise IOError("# PySCRDT::get_working_point: You need to define Madx"
                          + "twiss file in [prepare_data]")
        else:
            return self.actualQx, self.actualQy


    def get_order(self) -> Tuple[int]:
        """
        Returns the orders of the resonances
        Returns: [tuple]
        """

        if self.mode is None:
            raise IOError("# PySCRDT::get_order: You need to define resonance "
                          + "mode description in [set_mode]")

        elif self.mode == 3:
            if self.m is None:
                raise IOError("# PySCRDT::get_order: You need to define the "
                              + "order in [set_order]")
            else:
                return self.m, self.n, self.l

        elif self.mode == 5:
            if self.h is None:
                raise IOError("# PySCRDT::get_order: You need to define the "
                              + "order in [set_order]")
            else:
                return self.h, self.i, self.j, self.k, self.l


    def get_mode(self) -> int:
        """
        Returns the resonance mode description
        Returns: [int]
        """

        if self.mode is None:
            raise IOError("# PySCRDT::get_mode: You need to define resonance "
                          + "mode description in [set_mode]")
        else:
            return self.mode


    def get_ksc(self) -> float:
        """
        Returns the space charge perveance Ksc
        Returns: [float]
        """

        if self.K is None:
            self.ksc()

        return self.K


    def get_potential(self) -> Expr:
        """
        Returns the potential V
        Returns: [sympy expression]
        """
        if self.V is None:
            self.potential()

        return self.V


    def get_resonance_driving_terms(self, feed_down: Optional[bool] = False)\
                                                                       -> dict:
        """
        Returns the RDTs
        Inputs : feedDown : [bool] needed only if [resonance_driving_terms]
                            has not been already used
                            (default=False)
        Returns: [dict] the resonance driving terms
        """

        self.feed = feed_down

        if self.rdt is None:
            self.resonance_driving_terms(feed_down = self.feed)

        self.RDT={'RDT':self.rdt, 'Amplitude': abs(self.rdt),
                  'Phase': np.angle(self.rdt)}

        return self.RDT


    def get_detuning(self) -> dict:
        """
        Returns the RDTs
        Returns: [dict] the resonance driving terms
        """

        if self.rdt_d is None:
            self.detuning()

        return self.rdt_d


    def check_writing(self) -> dict:
        """
        Returns the correct writing format for setting or updating parameters
        Returns: [dict]
        """

        return {'Set & Update':['intensity', 'bunch_length', 'emittance_x',
                                'emittance_y', 'dpp_rms', 'dpp', 'ro'],
                                 'Update only': ['b', 'g']}
