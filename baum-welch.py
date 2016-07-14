from mpmath import *
import numpy as np
import matplotlib.pyplot as plt

slownik={"S":0,"N":1,"L":2, "R":3,"E":4}
#s="SNNLNNNNNNNLNNNNNRNNNNRNNNNLNNNNNNE"
#end_state="E"

emisja=matrix([[1.0,0.0,0.0,0.0,0.0],[0.0,0.96,0.036,0.004,0.0],[0.0,0.96,0.004,0.036,0.0],[0.0,0.0,0.0,0.0,1.0]])
przejscie=[[mpf(0.0),mpf(1.0),mpf(0.0),mpf(0.0)],[mpf(0.0),mpf(0.9998),mpf(0.0002),0.0],[0.0,0.0,mpf(0.9998),mpf(0.0002)],[0.0,0.0,0.0,1.0]]
#opens file
def robiliste(plik):
    o=open(plik)
    seqs=[]
    for line in o.readlines():
        seqs.append(line.strip())
    return seqs
#backward-forward algorithm    
def alfa(emisja,przejscie,seq,slownik):
    macierzalfa=matrix(4,len(seq))
    macierzbeta=matrix(4,len(seq))
    gdziezaczynamy=[1,0,0,0]
    for t in range(len(seq)):
        for i in range(4):
            if t == 0:
                macierzalfa[i, 0] = gdziezaczynamy[i] * emisja[i, 0]
            else:
                macierzalfa[i, t] = sum([macierzalfa[j, t-1]*przejscie[j][i] for j in range(4)]) * emisja[i, slownik[seq[t]]]
                
    return macierzalfa
            
def beta(emisja,przejscie,seq,slownik):
    macierzbeta=matrix(4,len(seq))
    gdziezaczynamy=[1,0,0,0]
    for t in range(len(seq)):
        for i in range(4):
            if t == 0:
                macierzbeta[i, len(seq)-1] = 1 
            else:
                tt = len(seq)-t-1
                macierzbeta[i, tt] = sum([emisja[j, slownik[seq[tt+1]]] * macierzbeta[j, tt+1] * przejscie[i][j] for j in range(4)])
    return macierzbeta
    
#Baum-welch algorithm    
def bw(emisja,przejscie,slownik):
    
    listaphi=[]
    oldest=mpf(0.0)
    est=mpf(1.0)
    eps=mpf("1e-6")
    seqs=robiliste("plik")
    for seq in seqs:
        #print seqs
        while True:
            a=alfa(emisja, przejscie,seq, slownik)
            b=beta(emisja, przejscie,seq, slownik)
            dziele= fsum([a[1, t]*b[1, t] for t in range(1, len(seq)-1)])
            oldest=mpf(est)
            est=fsum([a[i, 0] * b[i, 0] for i in range(4)])
            if fabs(est-oldest) > eps:
                phi=[]
                for t in range(len(seq)-1):
                    #print seq
                    phit=(a[1, t]*przejscie[1][2]*emisja[1, slownik[seq[t+1]]]*b[2, t+1]) / est
                    phi.append(phit)
                break
            przejscie[1][1] = (przejscie[1][1]*fsum([a[1, t]*emisja[1, slownik[seq[t+1]]] for t in range(1, len(seq)-1)])) / dziele
            przejscie[1][2] = 1 - przejscie[1][1]
            przejscie[2][3] = mpf(przejscie[1][2])
            przejscie[2][2] = mpf(przejscie[1][1])
            k1=[k for k,e in enumerate(list(seq)) if e == "N"]
            k2=[k for k,e in enumerate(list(seq)) if e == "L"]
            k3=[k for k,e in enumerate(list(seq)) if e == "R"]
            emisja[1, 1]= fsum([a[1, k] * b[1, k] for k in k1]) / dziele
            emisja[1, 2]= fsum([a[1, k] * b[1, k] for k in k2]) / dziele
            emisja[1, 3]= fsum([a[1, k] * b[1, k] for k in k3]) / dziele
            emisja[2, 1]=mpf(emisja[1, 1])
            emisja[2, 2]=mpf(emisja[1, 3])
            emisja[2, 3]=mpf(emisja[1, 2])
        listaphi.append(phi)
        #print phi
    print listaphi
    #print len(listaphi)
#makes a plot    
        
def robiwykresy():
    lista = bw(emisja,przejscie,slownik)
    for i in lista:
        plt.plot(i)
        zmienna=str(lista.index(i))
        plt.savefig(zmienna)
        plt.clf()

        
        
robiwykresy()
#robiwykresy() 
#robiliste("plik")

#bw(emisja,przejscie,slownik)
