from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import time
from sklearn.manifold import MDS, TSNE
import json
import math
import pickle

def computeScore(row):
    '''compute a single edit score wrt reference'''
    type1, ref1, alt1 = row[3:5+1]
    if type1 == 'snp':
        return 1
    else:
        return abs(len(ref1) - len(alt1))


def distRef(data):
    #chr1, start1, stop1, type1, ref1, alt1 = data[0]
    score = sum(computeScore(i) for i in data)
    return score


def distPair(data1, data2):
    '''data1 - list of tuples/rows,
    compute dist assuming they are in the same chromosome'''
    i, j = 0, 0
    score = 0
    while i<len(data1) and j<len(data2):
        #i,jth row
        chr1, start1, stop1, type1, ref1, alt1 = data1[i]
        chr2, start2, stop2, type2, ref2, alt2 = data2[j]
        
        if start1 == start2:
            if type1 == type2 == 'snp':
                if alt1 != alt2:
                    score += 1
            elif type1 == type2 == 'deletion':
                score += abs( len(ref1) - len(ref2) )
            elif type1 == type2 == 'insertion':
                score += abs( len(alt1) - len(alt2) )
            else: # different type
                #TODO more detailed diff?
                score += max(len(alt1), len(alt2), 
                            len(ref1), len(ref2))
            i+=1
            j+=1
        else:
            #count the diff of the smaller start
            if start1 < start2:
                score += computeScore(data1[i])
                i+=1

            elif start1 > start2:#TODO
                score += computeScore(data2[j])
                j+=1

    if i == len(data1) and j < len(data2):
        for j in range(j, len(data2)):
            score += computeScore(data2[j])
    if j == len(data2) and i < len(data1):
        for i in range(i, len(data1)):
            score += computeScore(data1[i])
        #TODO if 2 overlap partially?
    return score


def pklDump(obj, fn):
    with open(fn, 'w') as f:
        pickle.dump(obj, f)


def pklLoad(fn):
    with open(fn) as f:
        return pickle.load(f)
        

dataDict = {}
def getData(fn):
    if fn in dataDict:
        return dataDict[fn]
    t0 = time()
    with open(fn) as f:
        data = [ l[:-1].split(',') for l in f.readlines()]
    data = [(r[0], int(r[1]), int(r[2]), r[3], r[4], r[5]) for r in data]
    dataDict[fn] = data
    print 'file'+str(len(dataDict)), fn, time()-t0
    return data
 
        
def test1(chrNames='12345CM'):
    t0 = time()
    try:
        dist = pklLoad('pkl/dist_'+ ''.join(chrNames) +'.pkl')
    except IOError:
        dist = {}
    fileList = sorted(glob('data/*.txt'))
    count = 0
    for fn1 in fileList:
        eid1 = fn1.split('/')[-1].split('.')[0]
        for fn2 in fileList:
            eid2 = fn2.split('/')[-1].split('.')[0]
            if eid1 >= eid2: continue
            
            if (eid1,eid2) not in dist:
                data1 = getData(fn1)
                data2 = getData(fn2)
            
                #assuming chrNames same in data{1,2}
                '''reshpae data by chromosomes'''
                data1 = [[ l for l in data1 if l[0]==chrName ] for chrName in chrNames]
                data2 = [[ l for l in data2 if l[0]==chrName ] for chrName in chrNames]
            
            if (eid1,'ref') not in dist:
                #sum across all chromosomes
                dist[eid1,'ref'] = sum(distRef(d1) for d1 in data1)
            if (eid2,'ref') not in dist:
                dist[eid2,'ref'] = sum(distRef(d2) for d2 in data2)
            if (eid1, eid2) not in dist:
                dist[eid1, eid2] = sum(distPair(d1, d2) for d1,d2 in zip(data1, data2))
            count += 1
            print '%s: %s, %s, %d' % (chrNames, eid1, eid2, dist[eid1, eid2])

            if count % 10 == 0:
                pklDump(dist, 'pkl/dist_'+ ''.join(chrNames) +'.pkl')
                print 'saved dist', len(dist)
    pklDump(dist, 'pkl/dist_'+ ''.join(chrNames) +'.pkl')
    print 'saved dist', len(dist)
    print 'runtime', time()-t0
    return matrixOf(dist)


def matrixOf(dist):
    '''store dist dictionary into matrix form
    args    dist - distionary of pairwise distance
    return  m - 2d numpy array,
            eids - list of experiment id strings
    '''
    eids = sorted(set(  k[0] for k in dist.keys()  ))
    eids += ['ref']
    m = np.zeros([len(eids), len(eids)] )
    for i in range(len(eids)):
        for j in range(len(eids)):
            if i<j:
                m[i,j] = dist[eids[i], eids[j]]
            elif i>j:
                m[i,j] = dist[eids[j], eids[i]]

    return m, np.array(eids)


def plotMatrix(m, eids):
    plt.imshow(m)
    #plt.xticks(range(0,len(eids),5),eids[::5], rotation='vertical')
    #plt.yticks(range(0,len(eids),5),eids[::5])
    #plt.xticks(range(len(eids)),eids, rotation='vertical')
    #plt.yticks(range(len(eids)),eids)
    plt.colorbar()
    

def printGeoLocation():
    loc = {}
    for fn in glob('experiment/experiment_*/annotations.csv'):
        eid = fn.split('/')[-2].split('_')[-1]
        with open(fn) as f:
            l = f.readlines()
            lon, lat = l[3], l[2]
            lon = lon.split(',')[2].replace('"', '')
            lat = lat.split(',')[2].replace('"', '')
            print "%s, %s" % (lat, lon)
            loc[eid] = (lon, lat)
    print loc


def test2(m, eids, side=3, subplot=1, title=''):
    for i in range(len(m[0])-1):
        order = sorted( zip(m[i:, i:][0], range(len(m[0])-i)) )
        order =  [o[1] for o in order]
        m[i:,i:] = m[i:,i:][order]
        m[i:,i:] = m[i:,i:][:,order]
        eids[i:] = eids[i:][order]

    plt.figure(fig1.number)
    #plt.subplot(side,side,subplot)
    #ax = plt.gca()
    ##plt.subplot(121)
    plotMatrix(m, eids)
    plt.title(title)
    
    plt.figure(fig2.number)
    #plt.subplot(side,side,subplot)
    ##plt.subplot(122)
    plotScatter(m, eids)
    plt.title(title)


def plotScatter(m, eids):
    model = MDS(n_components=2,
                dissimilarity='precomputed',
                max_iter=15000)
    x = model.fit_transform(m)
    '''
    model = TSNE()
    x = model.fit_transform(x)'''
    for i in range(len(eids)):
        plt.scatter(x[i,0], x[i,1], c='C0')
        #plt.text(x[i,0], x[i,1], eids[i])#use xx in eid='13xx' for tag
    plt.yticks([])
    plt.axis('equal')
    
    

if __name__ == '__main__':
    
    chrNames = 'C'
    '''    
    for chrName in chrNames:
        m, eids= test1([chrName, ])
        print m, eids
    m, eids= test1()
    print m, eids
    '''

    fig1 = plt.figure()
    fig2 = plt.figure()
    side = math.ceil((len(chrNames)) **0.5)
    '''   
    for i,chrName in enumerate(chrNames):
        dist = pklLoad('pkl/dist_'+chrName+'.pkl')
        # take a small subset
        #dist = dict(i for i in dist.items() if i[0][0]<'1390' and (i[0][1]<'1390' or i[0][1]=='ref'))
        print len(dist)
        m, eids = matrixOf(dist)
        test2(m, eids, side=side, subplot=i+1, title='chr'+chrName)
    '''
    dist = pklLoad('pkl/dist_12345CM.pkl')
    m, eids = matrixOf(dist)
    test2(m, eids, subplot=side**2, title='all')
    plt.show()
