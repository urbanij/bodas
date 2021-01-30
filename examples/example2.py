"""
Wed Jan 27 2021 10:45:23 am CET

"""
import bodas
import sympy

"""matlab
s = tf('s');
w0 = 10000;

% just bode
for Q = [0.01, 0.1, 0.25, 0.5, 0.8, 1, 2, 10]
    bode(1 / (1 + s * Q/(2*w0) + (s/w0)^2));
    hold on;
end

% bode with asymptotes
for Q = [0.01, 0.1, 0.25, 0.5, 0.8, 1, 2, 10]
    % requires [this](https://nl.mathworks.com/matlabcentral/fileexchange/10183-bode-plot-with-asymptotes)
    asymptotic_bode(1 / (1 + s * Q/(2*w0) + (s/w0)^2));
    hold on;
end
"""


w0 = 10_000
s = sympy.Symbol('s')

"""
for Q in (0.01, 0.1, 0.25, 0.5, 0.8, 1, 2, 10):
    bodas.plot( 1 / (1 + s * Q/(2*w0) + (s/w0)**2) )
"""

Q = 0.2
bodas.plot( 1 / (1 + s * Q/(2*w0) + (s/w0)**2) )