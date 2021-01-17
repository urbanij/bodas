import bodas
import sympy
from sympy.abc import s
from sympy.physics.control.lti import TransferFunction

H = ' 100 * (1+s/100)*(1+s/553) / (s * (1+s) * (1+s/823) * (1+s/9822) ) '
#H = '100/s**5' # s**3 the phase is shifted
H = '140*(1+s/6043)/s'
#H = '140*(1+s/6043)/(s*(s+544))'
#H = '1/(1+s*0.001/1000+s**2/1000**2)'
#H = '1/(1+s*0.5/1000+s**2/1000**2)'
H = '1/((s**2) * (1+s/9000))' # phase doesnt workd
# H = '1/(s * (1+s/9000))'
# H = '140/(1+s/100)'
H = '140/(1+s/100)**2'
H = '((1+s/892) * (1-s/12.4))/((231+s))'
#H = '1/(1+s*0.5/1000+s**2/1000**2)'

#H = '-36249.85081752*(-3000000.0 - 300*(32312500.0 + 3.3e+15/(s*(33000.0 + 21276595.7446809/s)))/(11575.0 + 702127659574.468/(s*(33000.0 + 21276595.7446809/s))))/((1 + 1.12109783699823*(3000000.0 + 300*(32312500.0 + 3.3e+15/(s*(33000.0 + 21276595.7446809/s)))/(11575.0 + 702127659574.468/(s*(33000.0 + 21276595.7446809/s))))/((10072.2789115646 + 702127659574.468/(s*(33000.0 + 21276595.7446809/s)))*(10323.7484982265 + 1000000.0/s)))*(10072.2789115646 + 702127659574.468/(s*(33000.0 + 21276595.7446809/s)))*(87256.1403508772 + 301*(10099009.9009901 + 1000000000.0/s)/(11099.0099009901 + 1000000.0/s)))'
bodas.Tf(TransferFunction(*sympy.fraction(H), s)).plot()
