from bacteria import bacteria
from chemiotaxis import chemiotaxis
import matplotlib.pyplot as plt

import pandas as pd
import time
import numpy 
import numpy as np

fitness_prograssion = []
nfe_progression = []


poblacion = []
path = r"multifasta.fasta"
# numeroDeBacterias = np.random.randint(5, 10) #5
# numRandomBacteria = np.random.randint(1, 10) #1
# iteraciones = 150
# tumbo = np.random.randint(1, 3) #1                                              #numero de gaps a insertar 
# nado = np.random.randint(3, 10) #3
# chemio = chemiotaxis()
# veryBest = bacteria(path)                #mejor bacteria   
# tempBacteria = bacteria(path)            #bacteria temporal para validaciones
# original = bacteria(path)                #bacteria original sin gaps
# globalNFE = 0      #numero de evaluaciones de la funcion objetivo

# dAttr= round(np.random.uniform(0.1, 6), 1) #0.1
# wAttr= round(np.random.uniform(0.1, 6), 1) #0.2
# hRep=dAttr
# wRep= np.random.randint(10, 40)    #10



numeroDeBacterias = 7 #5
numRandomBacteria = 4 #1
iteraciones = 50
tumbo = 1 #1                                              #numero de gaps a insertar 
nado = 3 #3
chemio = chemiotaxis()
veryBest = bacteria(path)                #mejor bacteria   
tempBacteria = bacteria(path)            #bacteria temporal para validaciones
original = bacteria(path)                #bacteria original sin gaps
globalNFE = 0      #numero de evaluaciones de la funcion objetivo

dAttr= 0.1 #0.1
wAttr= 0.2 #0.2
hRep=dAttr
wRep= 10   #10


columns = ["time_stamp" , "epoca" ,  "numero_random_bacterias" , "numero_de_bacterias" , "numero_iteracion" , "tumbo" , "nado" , "dAttr" , "wAttr" , "wRep" , "fitness", "interaccion", "NFE"]
body = [] 
time_stamp = pd.Timestamp(time.time(), unit='s', tz='US/Pacific')
epoca = 1
itecion = 0


df_historical_algorithm = pd.read_csv("algorithm_historical_progrssion.csv")

def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    
def validaSecuencias(path, veryBest):
    #clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    #descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
    #tempBacteria.tumboNado(1)    

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return
      


for i in range(numeroDeBacterias):                                            #poblacion inicial
    poblacion.append(bacteria(path))


for _ in range(iteraciones):                                                  #numero de iteraciones  
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        #bacteria.tumboNado(nado)
        bacteria.autoEvalua()  
        #print("blosumScore: ",bacteria.blosumScore)
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)                 #d_attr, w_attr, h_rep, w_rep):
    globalNFE += chemio.parcialNFE 
    best = max(poblacion, key=lambda x: x.fitness)
    if (veryBest == None) or (best.fitness > veryBest.fitness):
         clonaBest(veryBest, best)
    print("interaccion: ",veryBest.interaction,"fitness: ",veryBest.fitness, " NFE:",globalNFE )
    fitness_prograssion.append(veryBest.fitness)
    nfe_progression.append(globalNFE)
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)                #inserta  bacterias aleatorias
    print("poblacion: ",len(poblacion))

    body.append([ time_stamp,  epoca , numRandomBacteria , numeroDeBacterias, itecion , tumbo , nado , dAttr , wAttr , wRep , veryBest.fitness, veryBest.interaction, globalNFE])
    itecion += 1

    if   fitness_prograssion.count(veryBest.fitness) >= 30:
        wRep += 1
        fitness_prograssion = []
        print("-------------------------------NEW wAttr ASSIGNED-----------------------------------------------------")

    if   fitness_prograssion.count(veryBest.fitness) >= 20:
        numRandomBacteria = np.random.randint(7, 15) #1
        tumbo = np.random.randint(10, 20) #1                                              #numero de gaps a insertar 
        nado = np.random.randint(3, 10) #3
        dAttr= round(np.random.uniform(3.24, 6), 1) #0.1
        wAttr= round(np.random.uniform(3.58, 6), 1) #0.2
        wRep= np.random.randint(15, 40)    #10
        fitness_prograssion = []
        print("-------------------------------NEW VALUES ASSIGNED-----------------------------------------------------")



veryBest.showGenome()
validaSecuencias(path, veryBest)


df_algorithm = pd.DataFrame(body , columns= columns)

df_historical_algorithm = pd.concat([df_algorithm , df_historical_algorithm])

df_historical_algorithm.to_csv("algorithm_historical_progrssion.csv" , index =  False)


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

ax1.set_title("Fitness Progression")
ax1.plot([i for i in range(len(fitness_prograssion))], fitness_prograssion, color="green", label='Fitness')
ax1.set_ylabel("Progression")
ax1.set_xlabel("Iterations")
ax1.legend()

ax2.set_title("NFE Progression")
ax2.plot([i for i in range(len(nfe_progression))], nfe_progression, color="red", label='NFE')
ax2.set_ylabel("Progression")
ax2.set_xlabel("Iterations")
ax2.legend()

plt.tight_layout()

plt.show()
