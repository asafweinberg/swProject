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
    # filesIndex = 2

    # if len(args) == 5:
    #     filesIndex = 3
    #     maxIter = args[2]

    #     try:
    #         maxIter = int(maxIter)
    #         if maxIter <= 0:
    #             print("invalid max iterations")
    #             return
    #     except:
    #         print("invalid max iterations")
    #         return
    
    # fileName1 = args[filesIndex]
    # fileName2 = args[filesIndex + 1]
    
    #reads the points from files

    pointsNewDim=km.runCFlow(k,goal,fileName) #if not spk, then equals an empty list. if spk equals Tmat as list of lists
    if len(pointsNewDim)==0:
        return
    # data, success = readDataSp(fileName) #returns panda

    # if(success == False):
    #     return
    
    # if k >= len(data):
    #     print("error: k >= n") SOMEWHERE TO CHECK
    #     return


    # for i in range(len(pointsNewDim)):
    #     for j in range(len(pointsNewDim[0])):
    #         print(str(round(pointsNewDim[i][j], 4)) + ',', end="")
    #     print("\n")

    # print("END")

    if(k==0):
        k=len(pointsNewDim[0])

    d = k
    
    points = np.array(pointsNewDim)
    numOfPoints = len(points)
    ###from here the same
    #computes the correct centroids
    centroids = calcCentroids(points, k) ##centroids list of tuples[index,data]  ; data is ndarray too

    # for tup in centroids:
    #     for n in tup[1]:
    #         print(str(round(n, 4)) + ',', end="")
    #     print("\n")

    # print("END")

    centroidsIndexes= [t[0] for t in centroids]

    # print(centroidsIndexes)
    dListCen=[t[1].tolist() for t in centroids]
    dListPoints=[]

    for pnt in points:
        pList=[]
        for num in pnt:
            pList.append(num)
        dListPoints.append(pList)

    # for p in dListCen:
    #     for n in p:
    #         print(str(round(n, 4)) + ',', end="")
    #     print("\n")
    
    finalMeans = km.fit(k, maxIter, numOfPoints, d, dListPoints, dListCen)
    # print("CHECK")
    indexesToPrint = ""

    for i in centroidsIndexes:
       indexesToPrint = indexesToPrint + str(i) + ","
    indexesToPrint = indexesToPrint[:-1]

    print(indexesToPrint)

    for fm in finalMeans[:-1]:
        for num in fm[:-1]:
            print(str(round(num, 4)) + ',', end="")
        print(str(round(fm[-1], 4)))
    
    fm=finalMeans[-1]
    for num in fm[:-1]:
        print(str(round(num, 4)) + ',', end="")
    print(str(round(fm[-1], 4)), end="")


    


def readData(fileName1, fileName2):
    data = None
    try:
        data1 = pd.read_csv(fileName1, header=None)
        data1.rename(columns={0 :'key'}, inplace=True)
        data2 = pd.read_csv(fileName2, header=None)
        data2.rename(columns={0 :'key'}, inplace=True)

    except:
        print("Error reading file")
        return None, False
    
    data = pd.merge(data1, data2, on = 'key')
    sortedPoints = data.sort_values(by=['key'])
    sortedPoints = sortedPoints.set_index('key')
    return sortedPoints, True

def readDataSp(fileName):
    data = None
    try:
        data = pd.read_csv(fileName, header=None)
        data.rename(columns={0 :'key'}, inplace=True)
        # data2 = pd.read_csv(fileName2, header=None)
        # data2.rename(columns={0 :'key'}, inplace=True)

    except:
        print("Error reading file")
        return None, False
    
    # data = pd.merge(data1, data2, on = 'key')
    sortedPoints = data.sort_values(by=['key'])
    sortedPoints = sortedPoints.set_index('key')
    return sortedPoints, True    


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
