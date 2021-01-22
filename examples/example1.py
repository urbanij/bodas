import sympy
from sympy.abc import s
from sympy.physics.control.lti import TransferFunction
import bodas

import multiprocessing

def proc_func(H):
    # print(H)
    tf = bodas.Tf(TransferFunction(*sympy.fraction(H), s))
    tf.plot()


def main():

    some_transfer_functions = [
        ' 100 * (1+s/100)*(1+s/5533) / (s * (1+s) * (1+s/823) * (1+s/9822) ) ',
        '1/(1+s*0.5/1000+s**2/1000**2)',
        '1/(s**2 * (1+s/9000))', 
        '0.14*s/(1+s/100)**2',
        '((1+s/8592) * (1-s/12.4))/((231+s)*(42-s))',
        '((1+s/8592)*(1+s/12.5)*(1+s/938) * (1-s/121))/((1+s/563)*(42-s))',
        '140*(1+s/3)/(1+s/100)**2',
        '1/((1+s/12)*(1+s/1233)*(1+s/91)*(1+s/9123)*(1-s/89212)*(1+s/2))'
    ]
    
    p = []
    for i in range(len(some_transfer_functions)):
        p.append( multiprocessing.Process(target=proc_func, args=(some_transfer_functions[i],)) )

    for proc in p:
        proc.start()

    for proc in p:
        proc.join()


if __name__ == '__main__':
    main()
