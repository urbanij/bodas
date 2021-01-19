# -*- coding: utf-8 -*-

#  _               _
# | |__   ___   __| | __ _ ___
# | '_ \ / _ \ / _` |/ _` / __|
# | |_) | (_) | (_| | (_| \__ \
# |_.__/ \___/ \__,_|\__,_|___/


"""
Bodas library
~~~~~~~~~~~~~~~~~~~~~

Bodas is a library written in Python, sitting on top of Sympy, for asymptotic Bode plots.
Basic usage:

   >>> import bodas
   >>> bodas.plot('1 / (1+s/5000) ')

:copyright: (c) 2021 by Francesco Urbani <https://urbanij.github.io>
:license: see LICENSE for more details.
"""



import sympy
from sympy.abc import s
from sympy.physics.control.lti import TransferFunction

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (9,10)

from typing import List

NUM_POINTS = 3000
w = np.logspace(-3, 8, NUM_POINTS)


class Tf:
    """docstring for TF"""
    def __init__(self, tf):
        self.tf: TransferFunction = tf
        self.H: str = tf.num / tf.den
        self.H_tf = sympy.lambdify(s, self.H, "numpy")
        # print(f"{self.tf.zeros()=}") # f-string debug variant is a Python3.8 feature 
        # print(f"{self.tf.poles()=}")
        
        # for visualization purposes
        self._abs_roots: List[sympy.core.numbers.Float] = list(map( abs, self.tf.zeros() + self.tf.poles() ))
        self._min_omega: sympy.core.numbers.Float       = max(w[0],  min(self._abs_roots)/100)
        self._max_omega: sympy.core.numbers.Float       = min(w[-1], max(self._abs_roots)*100)

    def _addSingularityContributionMagPlot(self, omega, singularity_type) -> np.ndarray:
        """
        """
        d_sign = {'pole': -1, 'zero': +1} 
        if omega == 0:
            f = d_sign[singularity_type] * 20 * np.log10(w)
        else:
            f = d_sign[singularity_type] * 20 * np.log10(w/float(omega)) * np.heaviside(w - float(omega), 1)
        return f

    def _addSingularityContributionPhasePlot(self, omega: complex, singularity_type: str, style: str = 'vertical') -> np.ndarray:
        
        d_sign = {'pole': -1, 'zero': +1} 

        sign = d_sign[singularity_type] * (1 if omega.real <= 0 else -1)
        abs_omega = abs(omega)

        if abs_omega == 0:
            f = sign * 90
        else:
            if style == 'vertical':
                f = sign * 90 * np.heaviside(w - float(abs_omega), 0)
            elif style == 'sloped':
                f = sign * \
                    (
                        45 * np.log10(10*w/abs_omega) * 
                        # ^-- 45 means 45 degrees per decade i.e. π/4 per decade.
                            ( np.heaviside(w * 10 - abs_omega, 0) - np.heaviside(0.1 * w - abs_omega, 0) ) + \
                            #  ^---- this line is basically a rect function spanning a decade before and after the actual singularity.
                        90 * np.heaviside(0.1 * w - abs_omega, 0) \
                    )
        return f

    def _calcLastNonZeroCoefficient(self, poly) -> float:
        """ https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly.EC """
        if sympy.degree(poly) > 0:
            print(f"{sympy.poly(poly)=}")
            print(f"{sympy.poly(poly).EC()=}")
            lnzc = sympy.poly(poly).EC()
        else:
            lnzc = poly
        return lnzc

    def _buildMagnitudeAsymptotes(self):
        f = 0 * w

        # constant
        f += 20 * np.log10( float( self._calcLastNonZeroCoefficient(self.tf.num) / self._calcLastNonZeroCoefficient(self.tf.den) ) )

        # zeros and poles
        for root in self.tf.zeros(): f += self._addSingularityContributionMagPlot(float(abs(root)), 'zero')
        for root in self.tf.poles(): f += self._addSingularityContributionMagPlot(float(abs(root)), 'pole')
        return f


    def _buildPhaseAsymptotes(self, style: str):
        f = 0 * w
        
        # add shitf correction here
        init_angle = np.angle(self.H_tf(1j * w), deg=True)[0] # starting angle of the actual bode plot function.
        f += 0 #init_angle
        
        for root in self.tf.zeros(): f += self._addSingularityContributionPhasePlot(complex(root), 'zero', style)
        for root in self.tf.poles(): f += self._addSingularityContributionPhasePlot(complex(root), 'pole', style)
        
        # f += init_angle
        return f
    
    def plot(self):

        LINEWIDTH_ACTUAL_PLOT, LINEWIDTH_ASYMP_PLOT = 0.6, 0.8
        
        H = self.H_tf(1j * w)
        H_db = 20 * np.log10(np.abs(H))
        H_phase = np.angle(H, deg=True)
        
        mag_asymptotes = self._buildMagnitudeAsymptotes()
        phase_asymptotes_sloped = self._buildPhaseAsymptotes('sloped')
        phase_asymptotes_vertical = self._buildPhaseAsymptotes('vertical')

        plt.figure(1)
        plt.suptitle(f"Bode plot: {self.H}")

        plt.subplot(211)
        plt.semilogx(w, H_db, 
            color="blue", 
            linestyle="dashed",
            linewidth=LINEWIDTH_ACTUAL_PLOT, 
            label=f'${self.H}$')
        plt.semilogx(w, mag_asymptotes, 
            color="red", 
            linestyle="solid",
            linewidth=LINEWIDTH_ASYMP_PLOT)
        
        plt.xlim(float(self._min_omega), float(self._max_omega))
        plt.ylabel("Magnitude (dB)")
        plt.grid(True, which='both', color='#9c9b97', linestyle='-.', linewidth=0.2)
        

        plt.subplot(212)
        plt.semilogx(w, H_phase, 
            color="blue", 
            linestyle="dashed",
            linewidth=LINEWIDTH_ACTUAL_PLOT, 
            label=f'${self.H}$')
        plt.semilogx(w, phase_asymptotes_sloped, 
            color="red", 
            linestyle="solid",
            linewidth=LINEWIDTH_ASYMP_PLOT, 
            label=f'${self.H}$')
        plt.semilogx(w, phase_asymptotes_vertical, 
            color="black", 
            linestyle=":",
            linewidth=0.4*LINEWIDTH_ASYMP_PLOT)
        
        plt.xlim(float(self._min_omega), float(self._max_omega))
        plt.ylabel("Phase (deg)")
        plt.xlabel("$\omega$ (rad/s)")
        plt.grid(True, which='both', color='#9c9b97', linestyle='-.', linewidth=0.2)
        
        plt.show()


if __name__ == '__main__':
    H = '((1+s/892) * (1-s/12.4))/((231+s))'
    tf = Tf(TransferFunction(*sympy.fraction(H), s))
    tf.plot()


