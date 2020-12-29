# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 15:50:18 2020

@author: MOHAMMAD SHARIQUE AHMAD
sharique.ahmad@ieee.org

"""

# BER RAYLEIGH FADING CHANNEL
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nr
import matplotlib
matplotlib.use('Agg')
fig, ax = plt.subplots()
plt.style.use('dark_background')

lBlocks = 10000;
nBlocks = 10000;
SNRdB = np.arange(1.0,45.0,3.0); #1 step size 3 and 45 not included
BERSimulation = np.zeros(len(SNRdB));  # for simualtion
BERTheoritical = np.zeros(len(SNRdB)); # for theorotical
SNR = 10**(SNRdB/10);


for blk in range (nBlocks):
    h = (nr.normal(0.0, 1.0, lBlocks)+1j*nr.normal(0.0, 1.0, lBlocks))/np.sqrt(2);
    noise = nr.normal(0.0, 1.0, lBlocks)+1j*nr.normal(0.0, 1.0, lBlocks);
    Sym = 2*nr.randint(2,size=lBlocks)-1;
    for K in range(len(SNRdB)):
        Tx = np.sqrt(SNR[K])*Sym;
        Rx = h*Tx + noise;
        DecBits = 2*(np.real(np.conj(h)*Rx)>0)-1;
        BERSimulation[K]=BERSimulation[K] + np.sum(DecBits != Sym);
    
BERSimulation = BERSimulation/lBlocks/nBlocks;
BERTheoritical = 1/2*(1-np.sqrt(SNR/(2+SNR)));
BERa = 1/2/SNR;
plt.yscale('log')
plt.plot(SNRdB, BERSimulation,'m-');
plt.plot(SNRdB, BERTheoritical,'ro');
plt.plot(SNRdB, BERa,ms=40, lw=4, alpha=0.9, mfc='orange''bs');
#plt.grid(1,which='both')
plt.suptitle('BER for Rayleigh Fading Channel')
plt.legend(["Simulation", "Theory","Approx"], loc ="best");
plt.xlabel('SNR(dB)---->')
plt.ylabel('Bit Error Rate ---->') 
fig.text(0.95, 0.05, 'Copyright Sharique',
         fontsize=20, color='white',
         ha='right', va='bottom', alpha=0.5)

