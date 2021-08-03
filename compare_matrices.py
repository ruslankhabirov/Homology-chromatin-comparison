#!/usr/bin/env python

from collections import defaultdict
from math import sqrt, log
import numpy as np
import matplotlib.pyplot as plt

chunk = 10000

def mass_center_function(signal):
    mass_x = 0
    mass_y = 0
    for x,y in signal.keys():
        mass_x += x*singnal[(x,y)]
        mass_y += y*singnal[(x,y)]                     
        
    total_mass = sum(signal.values())
    
    return (mass_x // total_mass, mass_y // total_mass)

def tanimoto_coeff(signal1, signal2):
    length_signal1 = len(signal1)
    length_signal2 = len(signal2)
    signal1 = set(signal1.keys())
    signal2 = set(signal2.keys())
    matches = len(signal1 & signal2)
    
    tanimoto = matches // (length_signal1 + length_signal2 - matches)
    
    return tanimoto

def readcoords(f, skipline=False):
    if skipline:
        f.readline()
    signal = {}
    for line in f.readlines():
        items = [item for item in line.split() if item != ""]
        x,y = int(items[0])//chunk, int(items[1])//chunk
        signal[(x,y)] = float(items[2])
    return signal

def writecoords(f, signals):
    for (x,y) in signals.keys():
        f.write("{}\t{}\t{}\n".format(x,y,signals[(x,y)] ))
 

def generate_neib(x,y):
    return set ( [ (x,y),
                    (x-1,y),(x+1,y),(x,y-1),(x,y+1),
                    (x-1,y-1),(x-1,y+1),(x+1,y-1),(x+1,y+1), 
                    (x-2,y),(x+2,y),(x,y-2),(x,y+2),   ] )
    
    
def plotheatmap(signals, filename):
    print("Populating matrix for ", filename)
    coords = set(signals.keys())
    sizex = max( (x for x,y in coords) )
    sizey = max( (y for x,y in coords) )
    size = max(sizex, sizey)
    # Converting to square matrix, logarithming, setting 0 for diagonal elements
    matrix = np.zeros((size, size))
    for x,y in coords:
        matrix[x-1,y-1] = log(signals[(x,y)] + 1) if x != y else 0
        matrix[y-1,x-1] = matrix[x-1,y-1]
    #print("Normalizing ", filename)
    #matrix = np.log(matrix + 1)
    
    print("Plotting heatmap for ", filename)
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.savefig(filename)
    
f1 = open("Human.chr21.fa.gz.21.seqcorr.txt")
f2 = open("chr21_10kb.RAWobserved")


print("Loading dataset1")
signals1 = readcoords(f1, skipline=True) 
plotheatmap(signals1,'data1.png')

print("Loading dataset2")
signals2 = readcoords(f2)
plotheatmap(signals2,'data2.png')

coords1 = set(signals1.keys())
coords2 = set(signals2.keys())

results = defaultdict(float)

D=2
i=1

L = len(coords1)
print("Start neibors processing")
for x1,y1 in coords1:
    # generating coordinates of all neibors
    neibors1 = generate_neib(x1,y1)
    matches = neibors1 & coords2 # now we know all points in signal2 which are in signal1 neibors
    
    for x2,y2 in matches:
        results[(x1,y1)] += signals1[(x1,y1)]*signals2[(x2,y2)]/(abs(x2-x1)+abs(y2-y1)+1) 
        
        
    #weightsum = sum( ( signals1[(x1,y1)]*signals2[(x2,y2)]/(abs(x2-x1)+abs(y2-y1)+1) for x2,y2 in matches ) )
    #results[(x1,y1)] += weightsum
    
    i+=1
    if i % 10000 == 0 :
        print("Processed {}/{} records".format(i,L))
    
print("Generated {} records".format(len(results)))
        
f3 = open("results_match.txt","w")
writecoords(f3,results)
plotheatmap(results,'data_results.png')
    
