import sys
import pandas as pd
from scipy.spatial import distance
import numpy as np
import mySpkmeansModule as km

def main(args):

    #check if arguments are valid
    if not (len(sys.argv) == 4):
        print("Error few args")
        return

    k = args[1]
    try: #put in c doc
        k = int(k)
        if k < 0:
            print("invalid k")
            return
    except:
        print("invalid k")
        return

    goal=args[2]
    fileName=args[3]
    maxIter = 300


    pointsNewDim=km.runCFlow(k,goal,fileName) #if not spk, then equals an empty list. if spk equals Tmat as list of lists
    if len(pointsNewDim)==0:
        return

    if(k==0):
        k=len(pointsNewDim[0])

    d = k
    
    points = np.array(pointsNewDim)
    numOfPoints = len(points)
    ###from here the same
    #computes the correct centroids
    centroids = calcCentroids(points, k) ##centroids list of tuples[index,data]  ; data is ndarray too

    centroidsIndexes= [t[0] for t in centroids]

    dListCen=[t[1].tolist() for t in centroids]
    dListPoints=[]

    for pnt in points:
        pList=[]
        for num in pnt:
            pList.append(num)
        dListPoints.append(pList)
    

    indexesToPrint = ""

    for i in centroidsIndexes:
       indexesToPrint = indexesToPrint + str(i) + ","
    indexesToPrint = indexesToPrint[:-1]

    print(indexesToPrint)

    finalMeans = km.fit(k, maxIter, numOfPoints, d, dListPoints, dListCen)  


def calcCentroids(points, k): # points is numpy ; #centroids list of tuples[index,data]  ; data is ndarray too
    np.random.seed(0)
    firstI = np.random.choice(len(points))
    centroids = [(firstI,points[firstI])]
    dists=[float("inf")]*len(points)
    for i in range(k-1):
        dists = getDistances(points, centroids,dists)
        probs = getProbabilities(dists)
        chosenCentroid = np.random.choice(len(points), p=probs)
        centroids.append((chosenCentroid, points[chosenCentroid]))

    return centroids

def getDistances(points, cents, dists):

    for i,p in enumerate(points):
        dst = distance.euclidean(p, cents[-1][1]) ** 2
        dists[i]=min(dists[i],dst)
    return dists


def getProbabilities(dists):
    sumD = sum(dists)
    probs = [d/sumD for d in dists]
    return probs

main(sys.argv)
