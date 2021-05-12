import numpy as np
import matplotlib.pyplot as plt
import math


def Ux0(t):
    return 0


def Uxl(t):
    return 0


def U(x):
    return 0

def ksi2(x):
    return 2*np.exp(-x)*np.sin(x)


def Solution(x, t):
    return np.exp(-t-x) * np.sin(x)*np.sin(2*t)


def progonka(a, b, c, d, s):
    P = np.zeros(s)
    Q = np.zeros(s)

    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]

    k = s - 1

    for i in range(1, s):
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1])
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])
    P[k] = 0
    Q[k] = (d[k] - a[k] * Q[k - 1]) / (b[k] + a[k] * P[k - 1])

    x = np.zeros(s)
    x[k] = Q[k]

    for i in range(s - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]
    return x

x0 = 0
xl = math.pi
# t = 2
param_a = 1
param_c = -3
param_b = 2


def autofill(x0, space_step, m, n, param_a, time_step, aprox_f):
    Uarray = np.zeros([n, m])

    tmp_x = x0
    for j in range(m):
        Uarray[0][j] = U(tmp_x)
		
        if aprox_f == 1:
            Uarray[1][j] = U(tmp_x) + ksi2(tmp_x)*time_step
			
        #if aprox_f == 2:
            #Uarray[1][j] = U(tmp_x) + ksi2(tmp_x)*time_step
				
        print(tmp_x)
        tmp_x += space_step
    
    return Uarray


def explicit(t, m, n, aprox, aprox_f):
    x0 = 0
    xl = math.pi

    space_step = (xl - x0) / (m - 1)
    time_step = t / (n - 1)


    Uarray = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)

    sigma = param_a**2 * time_step**2 / space_step**2
    print("time_step = ",time_step)
    print("sigma = ",sigma)
	
    for k in range(1, n - 1):
        for j in range(1, m - 1):
		
            Uarray[k + 1][j] = \
                Uarray[k][j + 1] *(sigma + time_step**2*param_b/(2*space_step))/(time_step+1) +\
                Uarray[k][j] * (2 - 2*sigma + param_c*time_step**2) /(time_step+1) + \
                Uarray[k][j - 1] *(sigma - time_step**2*param_b/(2*space_step))/(time_step+1) + \
                (Uarray[k - 1][j] * (-1)*(1-time_step))/(time_step+1)

        if aprox == 1:
            Uarray[k + 1][0] = 0
            Uarray[k + 1][m - 1] = 0
    for j in range(0,m):
        print(Uarray[49][j])
    return Uarray


def implicit(t, m, n, aprox, aprox_f):
    x0 = 0
    xl = math.pi

    space_step = (xl - x0) / (m - 1)
    time_step = t / (n - 1)


    Uarray = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)

    sigma = param_a**2 * time_step**2 / space_step**2

    alpha = 0
    betta = 1
    gamma = 0
    delta = 1
    tmp = time_step**2*param_b/(2*space_step)
    print("time_step =",time_step)
    print("sigma =",sigma) 

    for k in range(1, n - 1):
        a = np.zeros(m)
        b = np.zeros(m)
        c = np.zeros(m)
        d = np.zeros(m)

        for j in range(1, m - 1):
            a[j] = sigma - tmp
            b[j] = -(time_step + 1 + 2*sigma-param_c*time_step**2)
            c[j] = sigma + tmp
            d[j] = Uarray[k - 1][j]*(1-time_step) - 2*Uarray[k][j]
        if aprox == 1:
            b[0] = betta - alpha / space_step
            c[0] = alpha / space_step
            d[0] = Ux0((k + 1) * time_step)

            a[m - 1] = - gamma / space_step
            b[m - 1] = delta + gamma / space_step
            d[m - 1] = Uxl((k + 1) * time_step)


        Y = progonka(a, b, c, d, m)
        Uarray[k + 1] = Y
    for i in range(0, m):
        print('{:.30f}'.format(Uarray[20][i]))
    return Uarray

def Error(x0,x1,t,m,n,a,c,space_step,time_step,aprox,aprox_f):
    X = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)
    U = explicit(t, m, n, aprox, aprox_f)
    #E = np.zeros(m)
    for i in range(1,25):
         for j in range(0,1):
             #E = abs(U[i][j]- Solution(X[j],time_step*i,a,c))
             print(time_step*i)
             #print(abs(Ur[i][0]-np.sin(x0)))
             print((list(map('{:.30f}'.format,abs(U[i][1] - Solution(X[1],time_step*i))))),sep=" ")
             #print('{:.30f}'.format(abs(Ur[i][0]-np.sin(x0))))
             #print(abs(Ur[i][5]- Solution(X[5],time_step*i,a,c)))

aprox = 1
aprox_f = 1
m = 50
n = 50
t = 1
space_step = (xl - x0) / (m - 1)
time_step = t / (n - 1)
a = 1
c = -3
b = 2
x0 = 0
xl = math.pi
#explicit(t, m, n, aprox, aprox_f)
#implicit(t, m, n, aprox, aprox_f)
Error(x0,xl,t,m,n,a,c,space_step,time_step,aprox,aprox_f)