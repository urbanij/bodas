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

__title__ = 'bodas'
__version__ = '0.0.3'
__author__ = u'Francesco Urbani'

from .main import plot