# Solução do problema de sedimentação de uma esfera em regime de baixo Reynolds 
# via Método de Runge-Kutta de quarta ordem

# Importando a biblioteca NumPy
import numpy as np
import matplotlib.pyplot as plt

# Variáveis de entrada
p_s = 2500  # densidade da esfera (vidro) em [kg/mˆ3]
p_f = 1000  # densidade do fluido (fluido leve, ex. óleo/glicerina diluída) em [kg/mˆ3]
eta = 0.05  # viscosidade do fluido em [Pa·s]
g = 9.81  # gravidade em [m/sˆ2]
a = 0.002  # raio da esfera em [m]
# variar o a de 0.0005 a 0.003 m
v = np.zeros((100,2)) # matriz de tempo e velocidade - primeira coluna = tempo, segunda coluna = velocidade
v[0,0] = 0  # condição inicial t = 0
v[0,1] = 0  # condição inicial v = 0
h = 0.1  # passo escolhido
tolerancia = 0.01  # tolerância escolhida para o resultado
i = 0  # variável de contagem das iterações

    #U_s = (2*(a**2)*(p_s-p_f)*g)/(9*eta)  # velocidade de Stokes
    #Re_s = (p_f*U_s*a)/eta  # número de Reynolds de partícula
    #St = (2*a*p_s*U_s)/(9*eta)  # número de Stokes

# Definição da função a ser resolvida (equação adimensionalizada do movimento)
def f(v,Re_s,St):
    dvdt = (1 - v - (3/8)*Re_s*(v**2))/St
    return dvdt

# Criação da função que aplica o método Runge-Kutta de quarta ordem
def rk(t,v,i,Re_s,St):
    # Cálculo da Relação de Recorrência
    k1 = f(v,Re_s,St)
    x_2 = t + (1/2)*h
    y_2 = v + (1/2)*k1*h
    k2 = f(y_2,Re_s,St)
    y_3 = v + (1/2)*k2*h
    k3 = f(y_3,Re_s,St)
    x_4 = t + h
    y_4 = v + k3*h
    k4 = f(y_4,Re_s,St)
    v_novo = v + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h
    t_novo = (i+1)*h
    return t_novo,v_novo


# 1. Para o caso de Re → 0

Re_s = 0
erro = 1
t_final = 8
St = [0.2, 0.5, 1.0, 2.0]
# Loop de aplicação do método numérico
for k in St:
    n_passos = int(np.ceil(t_final / h))  # número de passos
    v = np.zeros((n_passos+1,2)) # matriz de tempo e velocidade - primeira coluna = tempo, segunda coluna = velocidade
    v[0,0] = 0  # condição inicial t = 0
    v[0,1] = 0  # condição inicial v = 0
    solucao = np.zeros((n_passos+1,1)) # matriz da solução analítica
    # Loop do Runge-Kutta para cada número de Stokes
    for i in range(n_passos):
        t,v = rk(v[i,0],v[i,1],i,Re_s,St)
        v[(i+1), 0] = t
        v[(i+1),1] = v
        # Cálculo da soluação analítica
        solucao[i+1,1] = 1 - np.exp(-(v[i+1,0])/St)
