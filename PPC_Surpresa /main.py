# PPC Surpresa - Gauss-Siedel

import numpy as np
import matplotlib.pyplot as plt

## Inputs
n = int(input()) # n° de nós internos em cada direção
T_0 = 20 # temperatura inicial em °C
tol = 10e-6 # tolerância adotada


## Criando a Matriz T
T = np.full((n+2,n+2), T_0)

for j in range (1,n+1):
    T[n+1,j] = 100 # temperatura inicial na borda inferior em °C

## Atualização Iterativa
iter = 0 # contador de iterações
erro = 1
iter_max =  100000 # limite de iterações

while erro > tol:
    erro = 0

    for j in range (1,n+1):
        for i in range (1,n+1):
            antigo = T[i,j]

            T[i,j] = (T[i-1,j] + T[i+1,j] + T[i,j-1] + T[i,j+1])/4

            diferenca = abs(T[i,j]-antigo)
            erro = max(erro,diferenca)
    
    iter = iter+1

    if iter > iter_max:
        print('Não houve convergência')
        break


## Gráfico

plt.imshow(T, cmap='viridis')
plt.title('Distribuição de Temperaturas')
plt.colorbar()
plt.show()

## Informações do processo
print('Número de iterações: ', iter)
print('Erro final: ', erro)
print('temperatura no ponto médio: ', T[3,3])

## 5) Questões rápidas:

# a) as temperaturas nas bordas não são atualizadas, pois a iteração de Gauss-Siedel
# calcula os novos pontos a partir de uma média de seus pontos vizinhos e as bordas não
# possuem 4 pontos vizinhos
