


import numpy as np
import matplotlib.pyplot as plt
import math



#Setando os parâmetros do modelo de acordo com o enunciado

N=5 #numero de elementos
V0=1 #voltagem da placa superior
D=0.001 #distancia entre placas
L=0.1 #tamanho da placa
E0 = 8.854187817*10**-12 #constante de permissividade

DELTA = L/N



v = []

#preenchendo o vetor das voltagens
#v=v0 apenas para a placa de cima, ou seja, para os N^2 ultimos elementos
v=[0 if i<N**2 else V0 for i in range(2*N**2)]



#configurando o vetor de posições
#cada entrada do vetor possui as coordenadas x,y e z do centro de um elemento da placa
r=[]
for k in range(2):
    for i in range(N):
        for j in range(N):
            r.append([j*DELTA+0.5*DELTA,i*DELTA+0.5*DELTA,k*D])

print(r[:10])



#funcao que calcula os elementos da matriz de impedancia a partir dos indices dos elementos, do vetor de posições e do delta
def impedance(r,delta,i,j):
    if i!=j:
        d = ((r[i][0] - r[j][0])**2 + (r[i][1] - r[j][1])**2 + (r[i][2] - r[j][2])**2)**(0.5)
        return 1/(4 * math.pi * E0) * (delta ** 2)/d

    return delta/(math.pi * E0) * math.log(1 + 2**(0.5))



z=[]

#preenchendo a matriz de impedancia
for i in range(2*N**2):
    z.append([])
    for j in range(2*N**2):
        z[i].append(impedance(r,DELTA,i,j))



#resolvendo o sistema linear para acharmos os coeficientes an
A = np.linalg.solve(z,v)

print(A[:10])



#funcao para plotar o grafico da densidade de carga em uma placa
def plot_plate(Xs,Ys,Zs,n,title):
    #trocando o formato de array 1d de formato n^2x1 para uma matriz 2d de formato nxn para podermos colocar no grafico de wireframe
    shaped_Xs = np.reshape(Xs,(n,n))
    shaped_Ys = np.reshape(Ys,(n,n))
    shaped_Zs = np.reshape(Zs,(n,n))

    ax = plt.axes(projection='3d')
    ax.plot_wireframe(shaped_Xs, shaped_Ys, shaped_Zs, color='black')
    ax.set_title(title)
    plt.show()



#placa inferior

new_r=np.array(r)

#vetores dos pontos x e y
#apenas a primeira metade dos pontos que correspondem a primeira placa
Xs = new_r[:len(new_r)//2,0]
Ys = new_r[:len(new_r)//2,1]
Zs = np.array(A[:len(A)//2])

plot_plate(Xs,Ys,Zs,N,"Bottom plate")



#placa superior

#vetores dos pontos x e y
#apenas a segunda metade dos pontos que correspondem a segunda placa
Xs = new_r[len(new_r)//2:,0]
Ys = new_r[len(new_r)//2:,1]
Zs = np.array(A[len(A)//2:]) 

plot_plate(Xs,Ys,Zs,N,"Top plate")


#Criando uma funcao genérica para fazer todos os calculos a partir de um N

def calc_density(v0,d,l,n):
    v = []

    delta = l/n

    #preenchendo o vetor das voltagens
    #v=v0 apenas para a placa de cima, ou seja, para os N^2 ultimos elementos
    v=[0 if i<n**2 else v0 for i in range(2*n**2)]

    #configurando o vetor de posições
    #cada entrada do vetor possui as coordenadas x,y e z do centro de um elemento da placa
    r=[]
    for k in range(2):
        for i in range(n):
            for j in range(n):
                r.append([j*delta+0.5*delta,i*delta+0.5*delta,k*d])

    z=[]

    #preenchendo a matriz de impedancia
    for i in range(2*n**2):
        z.append([])
        for j in range(2*n**2):
            z[i].append(impedance(r,delta,i,j))

    a = np.linalg.solve(z,v)

    #retornando os coeficientes an e o vetor de posições
    return (a,np.array(r))

 
# ### Agora podemos analisar vários valores de N diferentes


n=15
(an,r) = calc_density(V0,D,L,n)

Xs = r[:len(r)//2,0]
Ys = r[:len(r)//2,1]
Zs = np.array(an[:len(an)//2])

plot_plate(Xs,Ys,Zs,n,"Bottom plate ("+str(n)+")")

Xs = r[len(r)//2:,0]
Ys = r[len(r)//2:,1]
Zs = np.array(an[len(an)//2:]) 

plot_plate(Xs,Ys,Zs,n,"Top plate ("+str(n)+")")



n=30
(an,r) = calc_density(V0,D,L,n)

Xs = r[:len(r)//2,0]
Ys = r[:len(r)//2,1]
Zs = np.array(an[:len(an)//2])

plot_plate(Xs,Ys,Zs,n,"Bottom plate ("+str(n)+")")

Xs = r[len(r)//2:,0]
Ys = r[len(r)//2:,1]
Zs = np.array(an[len(an)//2:]) 

plot_plate(Xs,Ys,Zs,n,"Top plate ("+str(n)+")")



n=45
(an,r) = calc_density(V0,D,L,n)

Xs = r[:len(r)//2,0]
Ys = r[:len(r)//2,1]
Zs = np.array(an[:len(an)//2])

plot_plate(Xs,Ys,Zs,n,"Bottom plate ("+str(n)+")")

Xs = r[len(r)//2:,0]
Ys = r[len(r)//2:,1]
Zs = np.array(an[len(an)//2:]) 

plot_plate(Xs,Ys,Zs,n,"Top plate ("+str(n)+")")

 
# ### Podemos ver que a partir de N=30 o modelo converge bem para o resultado esperado, placas com cargas contrárias e grande concentração nas bordas.


#Capacitancia teorica
E0*L**2/D



n=10
(an,r) = calc_density(V0,D,L,n)

q = ((L/n)**2)*an[:len(an)//2].sum()
capacitance = q/V0
print(capacitance)



n=20
(an,r) = calc_density(V0,D,L,n)

q = ((L/n)**2)*an[:len(an)//2].sum()
capacitance = q/V0
print(capacitance)



n=30
(an,r) = calc_density(V0,D,L,n)

q = ((L/n)**2)*an[len(an)//2:].sum()
capacitance = q/V0
print(capacitance)



n=60
(an,r) = calc_density(V0,D,L,n)

q = ((L/n)**2)*an[len(an)//2:].sum()
capacitance = q/V0
print(capacitance)

 
# ### Aqui calculamos a capacitancia para diversos N diferentes e plotamos o gráfico do erro (diferença entre valor teórico e calculado).


#calculo dos coeficientes a_n
ans = []
ns = np.array(range(5,70,5))
for n in ns:
    (an,r) = calc_density(V0,D,L,n)

    ans.append((an,n))



#plot dos erros comparados com o modelo teórico
error = []
for an_tuple in ans:
    (an,n) = an_tuple
    q = ((L/n)**2)*an[len(an)//2:].sum()
    capacitance = abs(q)/V0
    error.append(capacitance-E0*L**2/D)

plt.plot(ns,error, label="erro")
plt.plot(ns,np.zeros(len(ns)), label="0")
plt.legend(loc="upper left")

 
# ### Como é possivel observar no gráfico, com o aumento de N o resultado converge para o valor analítico.

