# -*- coding: utf-8 -*-

#  _               _
# | |__   ___   __| | __ _ ___
# | '_ \ / _ \ / _` |/ _` / __|
# | |_) | (_) | (_| | (_| \__ \
# |_.__/ \___/ \__,_|\__,_|___/


"""
Bodas library
~~~~~~~~~~~~~~~~~~~~~

Bodas is a library written in Python, sitting on top of SymPy, for asymptotic Bode plots.
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

from typing import List

__title__ = 'bodas'
__version__ = '0.0.2'
__author__ = u'Francesco Urbani'


DEBUG = False

NUM_POINTS = 3000
w = np.logspace(-3, 8, NUM_POINTS)


class Tf:
    """docstring for TF"""
    def __init__(self, tf):
        self.tf: TransferFunction = tf
        self.H: str = tf.num / tf.den
        self.H_tf = sympy.lambdify(s, self.H, "numpy")
        if DEBUG:
            print(f"{self.tf.zeros()=}") # f-string debug version is a Python>=3.8 feature 
            print(f"{self.tf.poles()=}")
        
        # for visualization purposes
        self._abs_roots: List[sympy.core.numbers.Float] = list(map( abs, self.tf.zeros() + self.tf.poles() ))
        self._min_omega: sympy.core.numbers.Float       = max(w[0],  min(self._abs_roots)/100)
        self._max_omega: sympy.core.numbers.Float       = min(w[-1], max(self._abs_roots)*100) if \
                                                          min(w[-1], max(self._abs_roots)*100) > self._min_omega else w[-1] 
                                                            # otherwise in 1/s self._max_omega would be 0 (greater than self._min_omega)

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
        """
        """
        """ https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly.EC """
        if sympy.degree(poly) > 0:
            # print(f"{sympy.poly(poly)=}")
            # print(f"{sympy.poly(poly).EC()=}")
            lnzc = sympy.poly(poly).EC()
        else:
            lnzc = poly
        return lnzc

    def _buildMagnitudeAsymptotes(self):
        """
        """
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
        # print(f"{init_angle=}")

        f += 0 #init_angle
        
        for root in self.tf.zeros(): f += self._addSingularityContributionPhasePlot(complex(root), 'zero', style)
        for root in self.tf.poles(): f += self._addSingularityContributionPhasePlot(complex(root), 'pole', style)
        return f
    
    def plot(self):
        """
        """
        LINEWIDTH_ACTUAL_PLOT, LINEWIDTH_ASYMP_PLOT = 0.6, 0.8
        
        H = self.H_tf(1j * w)
        H_db = 20 * np.log10(np.abs(H))
        H_phase = np.angle(H, deg=True)
        
        mag_asymptotes = self._buildMagnitudeAsymptotes()
        phase_asymptotes_sloped = self._buildPhaseAsymptotes('sloped')
        phase_asymptotes_vertical = self._buildPhaseAsymptotes('vertical')

        # breakpoint()
    
        fig = plt.figure(figsize=(9,8))
        fig.canvas.set_window_title("github.com/urbanij/bodas")
        plt.suptitle("Bode plot of\n" + \
                     "$\\frac{{{0}}}{{{1}}}$".format(
                                                str(self.tf.num).replace('**', '^').replace('*','⋅'),
                                                str(self.tf.den).replace('**', '^').replace('*','⋅')))
        
        MAJOR_PLOT_ROW_SPAN = 10
        ax1 = plt.subplot2grid(shape=(MAJOR_PLOT_ROW_SPAN*2 + 1, 1), 
            loc=(0, 0), 
            rowspan=MAJOR_PLOT_ROW_SPAN)
        plt.semilogx(w, H_db, 
            color="blue", 
            linestyle="dashed",
            linewidth=LINEWIDTH_ACTUAL_PLOT, 
            label='actual')
        plt.semilogx(w, mag_asymptotes, 
            color="red", 
            linestyle="solid",
            linewidth=LINEWIDTH_ASYMP_PLOT,
            label='asymptotic')
        ax1.set_xticks([])
        # plt.ylim()
        plt.ylabel("Magnitude (dB)")
        plt.grid(True, which='both', color='#786E74', linestyle='-.', linewidth=0.18)
        plt.legend()

        ax2 = plt.subplot2grid(shape=(MAJOR_PLOT_ROW_SPAN*2 + 1, 1), 
            loc=(MAJOR_PLOT_ROW_SPAN, 0), 
            rowspan=1, 
            sharex=ax1)
        plt.scatter(
            list(map( abs, self.tf.zeros())),
            [0]*len(self.tf.zeros()),
            marker='o',
            facecolors="none", 
            color="red")
        plt.scatter(
            list(map( abs, self.tf.poles())),
            [0]*len(self.tf.poles()),
            marker='x',
            color="red")
        plt.setp(ax1.get_xticklabels(), visible=True)
        plt.ylim(-5,5)
        # plt.axes.get_xaxis().set_visible(False)  # remove the x-axis and its ticks
        ax2.yaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        # ax2.set_aspect(.2)

        ax3 = plt.subplot2grid(shape=(MAJOR_PLOT_ROW_SPAN*2 + 1, 1), 
            loc=(MAJOR_PLOT_ROW_SPAN+1, 0), 
            rowspan=MAJOR_PLOT_ROW_SPAN, 
            sharex=ax1)
        plt.semilogx(w, H_phase, 
            color="blue", 
            linestyle="dashed",
            linewidth=LINEWIDTH_ACTUAL_PLOT, 
            label='actual')
        plt.semilogx(w, phase_asymptotes_sloped, 
            color="red", 
            linestyle="solid",
            linewidth=LINEWIDTH_ASYMP_PLOT, 
            label='asymptotic')
        plt.semilogx(w, phase_asymptotes_vertical, 
            color="red", 
            linestyle=":",
            linewidth=0.88*LINEWIDTH_ASYMP_PLOT)
        
        # plt.xlim(float(self._min_omega), float(self._max_omega))
        plt.xlim(w[0], w[-1])

        plt.ylabel("Phase (deg)")
        plt.xlabel("$\omega$ (rad/s)")
        plt.grid(True, which='both', color='#786E74', linestyle='-.', linewidth=0.18)
        plt.legend()

        plt.show()
        # plt.savefig(f"_bodas.png")


def plot(H):
    """
    """
    if type(H) == str:
        H_str = H
    else: # e.g. sympy.core.mul.Mul
        H_str = str(H)

    tf = Tf(TransferFunction(*sympy.fraction( H_str ), s))
    tf.plot()

