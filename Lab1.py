
import numpy as np
import math as math


def Ux1(t, a,c):
    return np.exp((c-a) * t)


def Ux2(t, a,c):
    return np.exp((c-a) * t)


def U(x):
    return np.sin(x)


def Solution(x, t, a,c):
    return np.exp((c-a) * t) * np.sin(x)


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

# данные
# param_a = 1
# x0 = 0
# xl = mth.pi
# t = 2
# интервал с 21 по х и 2001 по t
# aprox = 2
# space_step = 0.157
# time_step = 0.001
#


def autofill(x0, space_step, m, n):
    Uarray = np.zeros([n, m])
    tmp_x = x0
    for j in range(0,m):
        Uarray[0][j] = U(tmp_x)
        tmp_x += space_step
        #print(tmp_x)
        print(tmp_x)

    return Uarray

def explicit(param_a, param_c, space_step, time_step, m, n):
    x0 = 0
    xl = math.pi/2


    Uarray = autofill(x0, space_step, m, n)

    sigma = (param_a**2 * time_step) / (space_step**2)

    tmp_time = 0
    for k in range(1, n):
        for j in range(1,m-1):
             Uarray[k][j] = sigma * (Uarray[k - 1][j + 1] + Uarray[k - 1][j - 1]) + (1 - 2 * sigma + param_c*time_step) * Uarray[k - 1][j]
        Uarray[k][m - 1] = Ux2(tmp_time, param_a, param_c)
        Uarray[k][0] = Uarray[k][1] - space_step * Ux1(tmp_time, param_a, param_c)
        

        tmp_time += time_step
        print(tmp_time)
    print(sigma)
    #print(time_step)
    for i in range(0,m):
        #print(Uarray[n-1][i])
        #print(format(Uarray[n-1][i],"f"))
        print("{0:.20f}".format(Uarray[1][i]))
    return Uarray

   


def implicit(param_a, param_c, space_step, time_step, m, n):
    x0 = 0
    xl = math.pi/2


    Uarray = autofill(x0, space_step, m, n)

    sigma = param_a**2 * time_step / space_step**2

    for k in range(0,1):
        a = np.zeros(m)
        b = np.zeros(m)
        c = np.zeros(m)
        d = np.zeros(m)

        # alpha = 1
        # beta = 0
        # gamma = 1
        # tetta = 0
        for j in range(1, m-1):
            a[j] = sigma
            b[j] = -(1 + 2 * sigma - param_c*time_step)
            c[j] = sigma
            d[j] = -Uarray[k][j]
       
        b[0] = -1 / space_step
        c[0] = 1 / space_step
        d[0] = Ux1((k + 1) * time_step,param_a,param_c)

        a[m - 1] = 0
        b[m - 1] = -1
        d[m - 1] = -Ux2((k+1)*time_step,param_a,param_c) 
        Y = progonka(a, b, c, d, m)
        Uarray[k + 1] = Y
        #print(time_step*(k+1))
        #print(Y)
        for i in range(0, m):
           print('{:.20f}'.format(Y[i]))
    return Uarray
    
        
        



def KN(param_a,param_c, space_step, time_step, m, n):
    x0 = 0
    xl = math.pi/2


    Uarray = autofill(x0, space_step, m, n)

    sigma = param_a**2 * time_step / space_step**2
    #print((xl-x0)/space_step)

    for k in range(n-1):
        a = np.zeros(m)
        b = np.zeros(m)
        c = np.zeros(m)
        d = np.zeros(m)

        # alpha = 1
        # beta = 0
        # gamma = 1
        # tetta = 0

        for j in range(1, m - 1):
            a[j] = -sigma / 2
            b[j] = (1 + sigma - (param_c*time_step)/2)
            c[j] = -sigma / 2
            d[j] = (sigma / 2) * Uarray[k][j + 1] + (1 - sigma + (param_c*time_step)/2) * Uarray[k][j] +  (Uarray[k][j - 1] )*(sigma / 2) 
        a[m - 1] = 0
        b[m - 1] = -1
        d[m - 1] = -Ux2((k+1)*time_step,param_a,param_c) 
        b[0] = -1 / space_step
        c[0] = 1 / space_step
        d[0] = Ux1((k + 1) * time_step,param_a,param_c)
        Y = progonka(a, b, c, d, m)
        Uarray[k + 1] = Y
        #print(time_step*(k+1))
        for i in range(0,m):
            Uarray[k+1][i] = Y[i]
            #print(Y[i])
            print('{:.20f}'.format(Y[i]))
        return Uarray

def Error(x0,x1,m,n,a,c,space_step,time):
    X = autofill(x0, space_step, m, n)
    Ur = implicit(a,c, space_step, time_step, m, n)
    #E = np.zeros(m)
    for i in range(1,n):
         for j in range(0,1):
             #E = abs(U[i][j]- Solution(X[j],time_step*i,a,c))
             print(time_step*i)
             #print(abs(Ur[i][0]-np.sin(x0)))
             print((list(map('{:.30f}'.format,abs(Ur[i][1] - Solution(X[1],time_step*i,a,c))))))
             #print('{:.30f}'.format(abs(Ur[i][0]-np.sin(x0))))
             #print(abs(Ur[i][5]- Solution(X[5],time_step*i,a,c)))
x0 = 0
xl = math.pi/2
m = 11
n = 11
a = 1
c = -2
space_step = (xl - x0) / (m - 1)
time_step = 10/(n-1)
#autofill(x0, space_step, m, n)
#implicit(a,c, space_step, time_step, m, n)
#explicit(a, c, space_step, time_step, m, n)
#print(time_step)
KN(a,c, space_step, time_step, m, n)
#Error(x0,xl,m,n,a,c,space_step,time_step)



