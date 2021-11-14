#Equation of State Solver by Newton-Raphson Method

from sympy import *
import numpy as np
e1 = str('Peng-Robinson')
e2 = str('Redlich-Kwong')
e3 = str('van der Waals')



def mathSolver(sel, r, t, p, a, b):
#in this function, we have to define some parameters:
# -sel corresponds to which equation will de solved(1 for Peng-Robinson, 2 for Redlich-Kwong, and 3
# for van der Waals
# -num is the initial value of x requested by Newton-Raphson method
# -MAXITER is the superior limit of interactions
# -r is the constant of gases
# -t is the temperature
# -p is the pressure
# -a is the constant for atraction in real gases(depends on the parameterization of selected equation)
# -b is the constant for repulsion in real gases(depends on the parameterization of selected equation)
    num = 20
    MAXITER = 1000
    if sel == 3:
        #eq = x**3 - (r * t/p + b)*x**2 + (a/p)*x - a*b/p
        x = var('x')
        f = Function('f')
        A = 1
        B = -(r*t/p)-b
        C = a/p
        D = - a*b/p
        f = A*x**3 + B*x**2 + C*x + D
        print(f'f(x): {f}')
        d = diff(f, x)
        print(f'df/dx: {d}')
        L = range(1, MAXITER+1)
        iteraction = 0
        xnew = num
        for i in L:
            xnewNew = xnew
            if d.subs(x, xnewNew) != 0:
                xnewNew = xnewNew - (f.subs(x, xnewNew)/d.subs(x, xnewNew))
                error = xnewNew - xnew
                xnew = xnewNew
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75

        for i in L:
            xnewNewb = -xnew
            if d.subs(x, xnewNewb) != 0:
                xnewNewb = xnewNewb - (f.subs(x, xnewNewb) / d.subs(x, xnewNewb))
                error = xnewNewb - xnew
                xnew = -xnewNewb
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75



    if sel == 2:
        #eq = x**3 - (r * t/p + b)*x**2 + (a/p)*x - a*b/p
        x = var('x')
        f = Function('f')
        A = 1
        B = - r*t/p
        C = - (b*r*t/p) - b**2 + (a/(p*(t**0.5)))
        D = - a*b/((t**0.5)*p)
        f = A*x**3 + B*x**2 + C*x + D
        print(f'f(x): {f}')
        d = diff(f, x)
        print(f'df/dx: {d}')
        L = range(1, MAXITER+1)
        iteraction = 0
        xnew = num
        for i in L:
            xnewNew = xnew
            if d.subs(x, xnewNew) != 0:
                xnewNew = xnewNew - (f.subs(x, xnewNew)/d.subs(x, xnewNew))
                error = xnewNew - xnew
                xnew = xnewNew
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75

        for i in L:
            xnewNewb = -xnew
            if d.subs(x, xnewNewb) != 0:
                xnewNewb = xnewNewb - (f.subs(x, xnewNewb) / d.subs(x, xnewNewb))
                error = xnewNewb - xnew
                xnew = -xnewNewb
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75

        for i in L:
            xnewNewc = -xnew/100
            if d.subs(x, xnewNewc) != 0:
                xnewNewc = xnewNewc - (f.subs(x, xnewNewc) / d.subs(x, xnewNewc))
                error = xnewNewc - xnew
                xnew = -xnewNewc/100
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75


    if sel == 1:
        #eq = x**3 - (r * t/p + b)*x**2 + (a/p)*x - a*b/p
        x = var('x')
        f = Function('f')
        A = 1
        B = b - (r*t/p)
        C = (a-(a*r*t*b)-(3*(b**2)*p))/p
        D = (((b**3)*p)+ (r*t*(b**2)) - (a*b))/p
        f = A*x**3 + B*x**2 + C*x + D
        print(f'f(x): {f}')
        d = diff(f, x)
        print(f'df/dx: {d}')
        L = range(1, MAXITER+1)
        iteraction = 0
        xnew = num
        for i in L:
            xnewNew = xnew
            if d.subs(x, xnewNew) != 0:
                xnewNew = xnewNew - (f.subs(x, xnewNew)/d.subs(x, xnewNew))
                error = xnewNew - xnew
                xnew = xnewNew
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75

        for i in L:
            xnewNewb = -xnew
            if d.subs(x, xnewNewb) != 0:
                xnewNewb = xnewNewb - (f.subs(x, xnewNewb) / d.subs(x, xnewNewb))
                error = xnewNewb - xnew
                xnew = -xnewNewb
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75

        for i in L:
            xnewNewc = -xnew/100
            if d.subs(x, xnewNewc) != 0:
                xnewNewc = xnewNewc - (f.subs(x, xnewNewc) / d.subs(x, xnewNewc))
                error = xnewNewc - xnew
                xnew = -xnewNewc/100
                iteraction = i
            else:
                iteraction = MAXITER + 1
                break
            if abs(error) <= 1e-6:
                break
        if iteraction > MAXITER:
            iteraction = 0.25
        elif iteraction == MAXITER:
            iteraction = 0.75

    print('-*'*30)
    title = 'RESULTS'
    print(title.center(60))
    print('-*'*30)
    print('Roots: %.4f and %.4f' %(xnewNew, xnewNewb))
