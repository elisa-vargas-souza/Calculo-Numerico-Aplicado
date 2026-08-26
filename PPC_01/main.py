# Solução do problema de sedimentação de uma esfera em regime de baixo Reynolds 
# via Método de Runge-Kutta de quarta ordem

# Importando a biblioteca NumPy
import numpy as np

# dvdt = (1 - v - (3/8)*Re_s*(v**2))/St

# Variáveis de entrada
v = np.zeros((100,2)) # matriz de velocidades - primeira coluna = v_x, segunda coluna = v_y
v[0,0] = 0  # condição inicial v_x = 0
v[0,1] = 0  # condição inicial v_y = 0
h = 0.1  #
tolerancia = 0.01  # tolerância escolhida para o resultado
i = 0  # variável de contagem das iterações

# Loop de aplicação do método Runge-Kutta
while erro > tolerancia:
  # Cálculo da Relação de Recorrência
  k1 = f(v[i,0],v[i,1])
  x_2 = v[i,0] + (1/2)*h
  y_2 = v[i,1] + (1/2)*k1*h
  k2 = f(x_2,y_2)
  y_3 = v[i,1] + (1/2)*k2*h
  k3 = f(x_2,y_3)
  x_4 = v[i,0] + h
  y_4 = v[i,1] + k3*h
  k4 = f(x_4,y_4)
  # Loop de cálculo das novas componentes da velocidade
  for k in range(2):
    v[i+1,k] = v[i,k] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h
  # Cálculo do erro
