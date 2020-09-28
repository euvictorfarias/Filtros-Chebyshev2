# -*- coding: utf-8 -*-

#Importa Bibliotecas Necessárias
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np


class chebyshev:
    
    # Iniciação do Objeto e seus parâmetros
    def __init__(self, tipo, Wp1, Ws1, Ap, As, Wp2 = 0, Ws2 = 0):
        self.tipo = tipo
        self.Wp = Wp1
        self.Ws = Ws1
        self.Ap = Ap
        self.As = As
        if tipo == "PF" or tipo == "RF":
            self.Wp1 = Wp1
            self.Wp2 = Wp2
            self.Ws1 = Ws1
            self.Ws2 = Ws2
        
    
    # Constante de Proporcionalidade
    def constProp(self):
        e = 1 / np.sqrt(pow(10, (-0.1*self.As)) - 1)
        self.e = e
        return e
    
    
    # Bandas de Passagem
    def bandas(self):
        Bp = self.Wp2 - self.Wp1
        Bs = self.Ws2 - self.Ws1
        self.Bs = Bs
        self.Bp = Bp
        return Bp, Bs
    
    
    # Essa função define e retorna a ordem do filtro
    def ordem(self):
        if self.tipo == "PB":
            n = np.arccosh(np.sqrt(pow(10, (-0.1*self.As)) - 1) / 
                           np.sqrt(pow(10, (-0.1*self.Ap)) - 1)
                           ) / np.arccosh(self.Ws/self.Wp)
        elif self.tipo == "PA":
            n = np.arccosh(np.sqrt(pow(10, (-0.1*self.As)) - 1) / 
                           np.sqrt(pow(10, (-0.1*self.Ap)) - 1)
                           ) / np.arccosh(self.Wp/self.Ws)
        elif self.tipo == "PF":
            n = np.arccosh(np.sqrt(pow(10, (-0.1*self.As)) - 1) / 
                           np.sqrt(pow(10, (-0.1*self.Ap)) - 1)
                           ) / np.arccosh(self.Bs/self.Bp)
        elif self.tipo == "RF":
            n = np.arccosh(np.sqrt(pow(10, (-0.1*self.As)) - 1) / 
                           np.sqrt(pow(10, (-0.1*self.Ap)) - 1)
                           ) / np.arccosh(self.Bp/self.Bs)
        N = int(np.ceil(n))
        
        if N % 2 == 0:
            Nv = N / 2
        else:
            Nv = (N - 1) / 2
        
        self.N = N
        self.Nv = int(Nv)
        return n, N
    
    
    # Essa função define e retorna a frequência de corte do filtro
    def freq_corte(self):
        if self.tipo == "PB" or self.tipo == "PA":
            Wc = self.Ws / np.cosh( (1/self.N) * np.arccosh(1/self.e * 
                np.sqrt(pow(10, -0.1*self.Ap)-1)) )
            self.Wc = Wc
            return Wc
        elif self.tipo == "PF" or self.tipo == "RF":
            Wc1 = self.Ws1 / np.cosh( (1/self.N) * np.arccosh(1/self.e* 
                np.sqrt(pow(10, -0.1*self.Ap)-1)) )
            Wc2 = self.Ws2 / np.cosh( (1/self.N) * np.arccosh(1/self.e* 
                np.sqrt(pow(10, -0.1*self.Ap)-1)) )
            self.Wc1 = Wc1
            self.Wc2 = Wc2
            return Wc1, Wc2
    
    
    # Frequência de Ressonância
    def freq_ress(self):
        Wo = np.sqrt(self.Ws1*self.Ws2)
        self.Wo = Wo
        return Wo
    
    
    # Essa função define e retorna as raízes do denominador da FT
    def raizes_unit(self):
        Sk = list()
        Ok = list()
        Wk = list()
        D = list()
        for k in range(1, self.N+1):
            D.append(pow(np.sinh((1/self.N) * np.arcsinh(1/self.e)), 2) * 
                  pow(np.sin(np.pi/(2*self.N) * (2*k - 1)), 2) + 
                  pow(np.cosh((1/self.N) * np.arcsinh(1/self.e)), 2) * 
                  pow(np.cos(np.pi/(2*self.N) * (2*k - 1)), 2)
                  )
            Ok.append(-np.sinh((1/self.N) * np.arcsinh(1/self.e)) * 
                      np.sin(np.pi/(2*self.N)*(2*k - 1)) / D[k-1]
            )
            Wk.append(-np.cosh((1/self.N) * np.arcsinh(1/self.e)) * 
                      np.cos(np.pi/(2*self.N)*(2*k - 1)) / D[k-1]
            )
            Sk.append(complex(Ok[k-1], Wk[k-1]))
        self.Sk = Sk
        return Sk
    
    
    # Essa função define e retorna a FT do filtro
    def func_tranf(self):
        chebyshev.raizes_unit(self)
        poli_polos = np.poly(self.Sk)
        polos = poli_polos.real
        poli_zeros = list()
        zeros = list()
        
        for i in range(1, self.Nv+1):
            Wz = 1 / np.cos((2*i - 1) * np.pi / (2*self.N))
            poli_zeros.append(complex(0, -Wz))
            poli_zeros.append(complex(0, Wz))
        
        poli_zeros = np.poly(poli_zeros)
        zeros = poli_zeros.real
        aux = polos[-1]
        
        for i in range(0, len(polos)):
            polos[i] = polos[i] * zeros[-1]
        for i in range(0, len(zeros)):
            zeros[i] = zeros[i] * aux
        
        if self.tipo == "PB":
            num, den = signal.lp2lp(zeros, polos, self.Ws)
        elif self.tipo == "PA":
            num, den = signal.lp2hp(zeros, polos, self.Ws)
        elif self.tipo == "PF":
            num, den = signal.lp2bp(zeros, polos, self.Wo, self.Bs)
        elif self.tipo == "RF":
            num, den = signal.lp2bs(zeros, polos, self.Wo, self.Bs)
        H = signal.TransferFunction(num, den)
        self.H = H
        return H


    # Essa função plota o Diagrama de Bode
    def plotar(self):
        # Plotagem do Módulo
        w, y, phase = self.H.bode(w = np.arange(0, 15000, step = 1))
        plt.figure(1)
        plt.grid(True)
        plt.xlim(0, 10000)
        plt.ylim(-100, 0)
        plt.plot(w, y)
        
        # Plotagem da Fase
        plt.figure(2)
        plt.grid(True)
        plt.xlim(0, 6000)
        plt.ylim(-360, 0)
        plt.plot(w, phase, 'r')
    
    
    
    
    
    
    
    
    
    
    
    