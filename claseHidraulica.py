# -*- coding: utf-8 -*-

import numpy as np
import math

pi = math.pi
 
""" Tuberías a presión """
class redPresion:
    
    def __init__(self, Q, I, D, k, nu, g=9.80665, tolerancia=1e-6):
        self.Q = Q # Caudal
        self.I = I # Pendiente motriz
        self.D = D # Diámetro
        self.k = k # Aspereza
        self.nu = nu # Viscosidad cinemática
        self.g = g # Aceleración gravedad
        self.tolerancia = 1e-6 # Diferencia máxima admitida entre iteracciones en los cálculos

    """
    Colebrook-White:
        - QWC cálculo del caudal
        - IWC cálculo de la pendiente motriz
        - DWC cálculo del diámetro
    """    
    def QWC(self):
        self.Q = np.where(np.logical_or(self.I==0, self.D==0), 0, pi*np.square(self.D)/4*np.sqrt(2*self.g*self.I*self.D)*(-2*np.log10(self.k/3.7/self.D + 2.51*self.nu/self.D/np.sqrt(2*self.g*self.I*self.D))))       

    def IWC(self):        
        IAnt = np.ones_like(self.Q)*1
        self.I = np.ones_like(self.Q)*1e-6

        while (np.abs(self.I-IAnt) > self.tolerancia).any():
            IAnt = self.I
            self.I = np.where(np.logical_or(self.Q==0, self.D==0), 0, np.where(self.Q/pi/self.D*4/self.nu < 2400, 128*self.nu*self.Q/pi/self.g/np.power(self.D, 4), np.square(4*self.Q/pi/np.square(self.D)/(-2*np.log10(self.k/3.7/self.D + 2.51*self.nu/self.D/np.sqrt(2*self.g*IAnt*self.D))))/2/self.g/self.D))
        
    def DWC(self):
        DAnt = np.ones_like(self.Q)*0
        self.D = np.ones_like(self.Q)*1
    
        while (np.abs(self.D-DAnt) > self.tolerancia).any():
            DAnt = self.D
            self.D = np.where(np.logical_or(self.Q==0, self.I==0), 0, np.power(4*self.Q/pi/np.sqrt(2*self.g*self.I)/(-2*np.log10(self.k/3.7/self.D + 2.51*self.nu/DAnt/np.sqrt(2*self.g*self.I*DAnt))),2/5))

""" Sistemas de bombeo """
class bombeo:
    
    def __init__(self, a, b, N):
        self.a = a # Vector con los coeficientes a0, a1 y a2 de la curva Q-H
        self.b = b # Vector con los coeficientes b1 y b2 de la curva Q-eta
        self.N = N # Vector con las relaciones de velocidades rotación de cada bomba, de la primera a la última. 0 indica bomba parada, 1 indica bomba a velocidad de rotación nominal
        
""" Canales de sección trapecial """
class canalTrapecial:
    
    def __init__(self, Q=None, b=None, z=None, n=None, I=None, I0=None, y=None, QvPG=None, y1PG=None, QvPD=None, y1PD=None, pPG=None, lPG=0, pPD=None, yC=None, g=9.80665, tolerancia=1e-6):
        self.Q = Q # Caudal
        self.I = I # Pendiente motriz
        self.b = b # Ancho de la solera
        self.z = z # Inclinación quijeros 1V:zH
        self.n = n # Coeficiente Ec. Manning
        self.I0 = I0 # Pendiente de la solera
        self.y = y # Calado
        self.QvPG = QvPG # Caudal vertedero pared gruesa
        self.y1PG = y1PG # Calado en la aproximación al vertedero de pared gruesa
        self.QvPD = QvPD # Caudal vertedero pared delgada
        self.y1PD = y1PD # Caudal en la aproximación al vertedero de pared delgada
        self.pPG = pPG # Altura de la coronación del vertedero de pared gruesa
        self.lPG = lPG # Longitud de la coronación del vertedero de pared gruesa
        self.pPD = pPD # Altura de la coronación del vertedero de pared delgada
        self.yC = yC # Calado crítico
        self.g = g # Aceleración gravedad
        self.tolerancia = tolerancia # Diferencia máxima admitida entre iteracciones en los cálculos
    
    """
    Manning para sección trapecial:
        - QManning cálculo del caudal
        - IManning cálculo de la pendiente motriz
        - yManning cálculo del calado
    """
    def QManning(self):
        self.Q = np.where(np.logical_or(self.I==0, self.y==0), 0, 1/self.n*np.power(self.b*self.y+self.z*np.square(self.y), 5/3)/np.power(self.b+self.y*np.sqrt(1+np.square(self.z)), 2/3)*np.sqrt(self.I))
        
    def IManning(self):
        try:
           self.I = np.where(np.logical_or(self.Q==0, self.y==0), 0., np.square(self.Q*self.n/np.power(self.b*self.y+self.z*np.square(self.y), 5/3)*np.power(self.b+2*self.y*np.sqrt(1+np.square(self.z)), 2/3)))
        except:
            self.I = np.zeros_like(self.Q)
            
    def yManning(self):
        self.y = self.y # Pendiente
    
    # Cálculo del calado crítico
    def yCrit(self):
        self.yC = np.ones_like(self.Q)
        yCAnt = np.zeros_like(self.Q)
        while (abs(yCAnt-self.yC) > self.tolerancia).any():
            yCAnt = self.yC
            self.yC = np.where(self.Q==0, 0, np.power(np.square(self.Q)*(self.b+2*self.z*yCAnt)/self.g/np.power(self.b+self.z*yCAnt, 3), 1/3))
    
    """
    Vertedero de pared gruesa:
        - y3Lim Cálculo del calado aguas abajo de aforador de pared gruesa correspondiente al límite de aforo modular
        - QvertPG Cálculo del caudal (no sumergido)
    """                    
    def y3Lim(self):
        self.yCrit()
        
        # Sánchez y Juana (2022)  DOI: https://doi.org/10.4995/ia.2022.16788
        ImpUPeso = np.square(self.pPG+self.yC)*(3*self.b+2*self.z*(self.pPG+self.yC))/6+np.square(self.QvPG)/self.g/((self.b+2*self.z*self.pPG)*self.yC+self.z*np.square(self.yC))
        
        y3LimAnt = np.zeros_like(self.QvPG)
        self.y3Lim = self.pPG + self.yC
        while (abs(y3LimAnt-self.y3Lim) > self.tolerancia).any():
            y3LimAnt = self.y3Lim
            ImpUPeso2 = np.square(self.QvPG)/self.g/(self.b*y3LimAnt+self.z*np.square(self.y3Lim))
            ImpUPeso1 = (3*self.b+2*self.z*y3LimAnt)/6
            self.y3Lim = np.power((ImpUPeso-ImpUPeso2)/ImpUPeso1, 1/2)
        
    
    def QvertPG(self):
        QAnt = np.ones_like(self.y1PG)*10*self.tolerancia
        self.QvPG = np.zeros_like(self.y1PG)
        while (abs(QAnt-self.QvPG) > self.tolerancia).any():
            self.Q = self.QvPG
            self.yCrit()            
            self.y=self.yC
            self.IManning()
            hf = self.lPG*self.I
            QAnt = self.QvPG
            self.QvPG = self.b*np.sqrt(self.g)*np.power(2/3*(self.y1PG-self.pPG+np.square(QAnt/(self.b*self.y1PG))/2/self.g-hf),3/2)
            
            
    """
    Vertedero de  pared delgada:
        - QvertPD Cálculo del caudal (no sumergido)
    """
    def QvertPD(self):
        QAnt = np.ones_like(self.y1PD)
        self.QvPD = np.zeros_like(self.y1PD)
        while (abs(QAnt-self.QvPD) > self.tolerancia).any():
            QAnt = self.QvPD
            
            # Consideración u omisión de la altura cinética en la aproximación
            
            # SumCinet = np.square(QAnt/(self.b*self.y1PD))/2/self.g
            SumCinet = 0*np.square(QAnt/(self.b*self.y1PD))/2/self.g #Ver Rouse - Elementary Mechanics of Fluids p.93
            
            # Coeficiente de contracción de Rouse para vertederos rectangulares
            Cc = 0.611+0.075*(self.y1PD-self.pPD)/self.pPD
            
            self.QvPD = self.b*Cc*np.sqrt(2*self.g)*2/3*(np.power(self.y1PD-self.pPD+SumCinet, 3/2)-np.power(SumCinet, 3/2))
            
""" Canales de sección circular """
class canalCircular:
    
    def __init__(self, Q, I, y, D, n, g=9.80665):
        self.Q = Q
        
""" Riego por goteo """
class unidadGoteo:
    
    def __init__(self, Dr, Dpr, Lr, Lpr, sg, sr, leg, ler, g=9.80665):
        self.Dr = Dr

