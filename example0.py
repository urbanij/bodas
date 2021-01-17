import bodas
import sympy
from sympy.abc import s
from sympy.physics.control.lti import TransferFunction

H = '((1+s/892) * (1-s/12.4))/((231+s))'
bodas.Tf(TransferFunction(*sympy.fraction(H), s)).plot()
