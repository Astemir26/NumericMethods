import numpy as np

def Analitic(x,y):
    return y*np.sin(x)


def Fi1(y):
    return 0



def Fi2(y):
    return y



def Fi3(x):
    return np.sin(x)


def Fi4(x):
    return 0


def Get_Norma(Cur,Prev,len_x,len_y):
    my_max = 0
    for i in range (0,len_x):
        for j in range(0,len_y):
            if (abs(Cur[i][j] - Prev[i][j]) > my_max):
                my_max = abs(Cur[i][j] - Prev[i][j])
    return my_max


iterat = 0


def Simple_Iteration_Method(hx, hy, epsilon):
    len_x = int(np.pi/(2*hx)+1)
    len_y = int(1/(hy)+1)
    print(len_x,len_y)
    x = np.zeros((len_x,len_y))
    y = np.zeros((len_x,len_y))
    Prev_T = np.zeros((len_x,len_y))
    Cur_T = np.zeros((len_x,len_y))
    
    for i in range(0,len_x):
        x[i] = hx*i
    for j in range(0,len_y):
        y[j] = hy*j
    for j in range(0,len_y):
        for i in range(0,len_x):
            Cur_T[i][j] = 0
		
	#интерполяция
    for j in range(0,len_y):
        coeff = (Fi2(y[j]) - Fi1(y[j]))/(len_x)
        addition = Fi1(y[j])
        print(coeff)	
        for i in range(0,len_x):
            Cur_T[i][j] = coeff*i

    for j in range(0,len_y):
        Cur_T[0][j] = Fi1(y[j])
        Cur_T[len_x - 1][j] = Fi2(y[j])
	

    while(1):
        for i in range(0,len_x):
            for j in range(0,len_y):
                Prev_T[j]=Cur_T[i][j]
        for i in range(1,len_x-1):
            Cur_T[i][0] = Prev_T[i][1] - Fi3(x[i])*hy
            for j in range(1,len_y-1):
                Cur_T[i][j] = (hy**2)*(Prev_T[i-1][j] + Prev_T[i+1][j] + \
                                       (hx**2) * (Prev_T[i][j+1] + Prev_T[i][j-1]))/\
                    (2*(hx**2)+2*(hy**2)-(hx*hy)**2)
            
            Cur_T[i][len_y-1] = (Fi4(x[i])+Prev_T[i][len_y-2])/(1-hy)		

        if (Get_Norma(Cur_T, Prev_T, len_x, len_y)<=epsilon):
            break				
        for i in range(0,len_x):
            for j in range(0,len_y):
                print(Cur_T[i][j]-Analitic(x[i],y[j]))
        
		
def Zeidel_Method(hx, hy, epsilon):
    len_x = np.pi/(2*hx)+1
    len_y = 1/(hy)+1
    len_x = int(len_x)
    len_y = int(len_y)
    print(len_x,len_y)
    
    x = np.zeros((len_x,len_y))
    y = np.zeros((len_x,len_y))
    Prev_T = np.zeros((len_x,len_y))
    Cur_T = np.zeros((len_x,len_y))
    for i in range(0,len_x):
        x[i] = hx*i
	
	
    for j in range(0,len_y):
        y[j] = hy*j
    for j in range(0,len_y):
        for i in range(0,len_x):
            Cur_T[i][j] = 0

#//интерполяция	
    for j in range(0,len_y):
        coeff = (Fi2(y[j]) - Fi1(y[j])) / len_x
        addition = Fi1(y[j])
        for i in range(0,len_x):
            Cur_T[i][j] = coeff[i]*i

    for j in range(0,len_y):
        Cur_T[0][j] = Fi1(y[j])
        Cur_T[len_x - 1][j] = Fi2(y[j])

    while(1):
        for i in range(0,len_x):
            for j in range(0,len_y):
                Prev_T[j]=Cur_T[i][j]

        for i in range(1 ,len_x-1):
            Cur_T[i][0] = Prev_T[i][1] - Fi3(x[i])*hy
            for j in range(1,len_y-1):
                Cur_T[i][j] = (hy**2)*(Cur_T[i-1][j] + Prev_T[i+1][j] + \
                                       (hx**2) * (Prev_T[i][j+1]+Cur_T[i][j-1]))/\
                    (2*(hx**2)+2*(hy**2)-(hx*hy)**2)
            Cur_T[i][len_y-1] = Fi4(x[i]) + Prev_T[i][len_y-2]/(1-hy)
		
        if(Get_Norma(Cur_T, Prev_T,len_x,len_y) <= epsilon):
            break
        for i in range(0,len_x):
            for j in range(1,len_y-1):
                print(Cur_T[i][j]-Analitic(x[i],y[j]))
def Relaxation_Method(hx, hy, epsilon):
    len_x = np.pi/(2*hx)+1
    len_y = 1/(hy)+1
    len_x = int(len_x)
    len_y = int(len_y)
    print(len_x,len_y)
    x = np.zeros((len_x,len_y))
    y = np.zeros((len_x,len_y))
    Prev_T = np.zeros((len_x,len_y))
    Cur_T = np.zeros((len_x,len_y))
    for i in range(0,len_x):
        x[i] = hx*i
	
    for j in range(0,len_y):
        y[j] = hy*j
	
    for j in range(0,len_y):

        for i in range(0,len_x):
            Cur_T[i][j] = 0
	#интерполяция
    for j in range(0,len_y):
        coeff = (Fi2(y[j]) - Fi1(y[j])) / len_x	
        addition = Fi1(y[j])		

        for i in range(0,len_x):
            Cur_T[i][j] = coeff*i

    for j in range(0,len_y):
        Cur_T[0][j] = Fi1(y[j])
        Cur_T[len_x - 1][j] = Fi2(y[j])
        relax_param = 1.5
    while(1):
        for i in range(0,len_x):
            for j in range(0,len_y):
                Prev_T[j]=Cur_T[i][j]
        for i in range(1,len_x-1):
            Cur_T[i][0] = Prev_T[i][1] - Fi3(x[i])*hy
            for j in range(1,len_y-1):
                M = (hy**2)*(Cur_T[i-1][j] + Prev_T[i+1][j] + \
                                       (hx**2) * (Prev_T[i][j+1]+Cur_T[i][j-1]))/\
                    (2*(hx**2)+2*(hy**2)-(hx*hy)**2)
                Cur_T[i][j] = (1 - relax_param) * Prev_T[i][j] + relax_param * M
            Cur_T[i][len_y-1] = (Fi4(x[i])+Prev_T[i][len_y-2])/(1-hy)
		
        if(Get_Norma(Cur_T, Prev_T,len_x,len_y) <= epsilon):
            break
        for i in range(0,len_x):
            for j in range(0,len_y):
                print(Cur_T[i][j]-Analitic(x[i],y[j]))

Simple_Iteration_Method(0.1, 0.1,0.1)
#Zeidel_Method(0.1, 0.1,0.1)
#Relaxation_Method(0.01, 0.1,0.1)