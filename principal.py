# -*- coding: utf-8 -*-

from funcoes_cheb2 import *

# Boas Vindas
print("\nBem Vindo(a)!\n")
print("--------------------------------------------------------------------")
print("Tipos de Filtros Chebyshev 2:")
print("(PB) - Passa-Baixa\n(PA) - Passa-Alta")
print("(PF) - Passa-Faixa\n(RF) - Rejeita-Faixa")
print("--------------------------------------------------------------------")

# Pega o tipo de filtro e os pontos de projeto
tipo = input("Digite a SIGLA do tipo que deseja: ")
print("\nAgora vamos aos Pontos de Projeto ... ")
if tipo == "PB" or tipo == "PA":
    Wp = float(input("Digite a Frequência de Passagem (Wp): "))
    Ws = float(input("Digite a Frequência de Rejeição (Ws): "))
    Ap = float(input("Digite a Atenuação de Passagem (Ap): "))
    As = float(input("Digite a Atenuação de Rejeição (As): "))
elif tipo == "PF" or tipo == "RF":
    Wp1 = float(input("Digite a Frequência de Passagem (Wp1): "))
    Wp2 = float(input("Digite a Frequência de Passagem (Wp2): "))
    Ws1 = float(input("Digite a Frequência de Rejeição (Ws1): "))
    Ws2 = float(input("Digite a Frequência de Rejeição (Ws2): "))
    Ap = float(input("Digite a Atenuação de Passagem (Ap): "))
    As = float(input("Digite a Atenuação de Rejeição (As): "))

# Inicializa um Objeto da Classe Butterworth
if tipo == "PB" or tipo == "PA":
    filtro = chebyshev(tipo, Wp, Ws, Ap, As)
elif tipo == "PF" or tipo == "RF":
    filtro = chebyshev(tipo, Wp1, Ws1, Ap, As, Wp2, Ws2)

print("\n--------------------------------------------------------------------")
print("As definições do seu filtro são as seguintes:")
print("--------------------------------------------------------------------")

# Define e exibe a Constante de Proporcionalidade
e = filtro.constProp()
print(f'(e)  - Constante de Proporcionalidade: {e:.5f}')

# Define e exibe frequência de ressonancia e bandas de passagem
if tipo == "PF" or tipo == "RF":
    Wo = filtro.freq_ress()
    print(f'(Wo) - Frequência de Ressonância: {Wo:.5f}')
    
    Bp, Bs = filtro.bandas()
    print(f'(Bp) - Banda de Passagem: {Bp:.5f}')
    print(f'(Bs) - Banda de Rejeição: {Bs:.5f}')

# Define e exibe ordem do filtro
n, N = filtro.ordem()
print(f'(N)  - Ordem: {N} ({n:.5f})')

# Define e exibe frequência de corte
if tipo == "PB" or tipo == "PA":
    Wc = filtro.freq_corte()
    print(f'(Wc) - Frequência de Corte: {Wc:.5f}')
elif tipo == "PF" or tipo == "RF":
    Wc1, Wc2 = filtro.freq_corte()
    print(f'(Wc1) - Frequência de Corte 1: {Wc1:.5f}')
    print(f'(Wc2) - Frequência de Corte 2: {Wc2:.5f}')
    
# Define e exibe Função de Transferência
H = filtro.func_tranf()
print("--------------------------------------------------------------------")
print("Sua Função de Transferência H(s) é:")
print(H)
print("--------------------------------------------------------------------")

# Plotar Gráficos
print("Deseja Plotar os gráficos de Bode?")
resposta = input("'s' ou 'n' (sem aspas): ")
if resposta == 's':
    filtro.plotar()
elif resposta == 'n':
    print("Gráfico não iniciado!")
else:
    print("Resposta não identificada.")
print("--------------------------------------------------------------------\n")
