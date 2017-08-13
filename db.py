import subprocess as sp
import pandas as pd
import numpy as np
from natsort import natsorted

EXP_PATH_PREFIX = 'experiment/experiment_'

def getRawData(dirName= EXP_PATH_PREFIX +'10838', limit=9e15):
    x = sp.Popen('ibis -d %s -q "select ' % dirName+\
                'chr,start,stop,type,ref,alt ' +\
                'where 1=1 order by chr start limit %d"' % limit,
                shell=True, 
                stderr=sp.PIPE,
                )
    gen = (i.replace('\n','')
            .replace('"','')
            .split(', ') for i in x.stderr.readlines()[3:-1])
    df = pd.DataFrame(gen, dtype=np.int64,
            columns= ['chr', 'start', 'stop', 'type','ref', 'alt'])
    return df


def getCdf(series):
    '''input - a pandas series'''
    l = series.as_matrix()
    cdf = []
    for i,s in enumerate(l):
        if i==0 or s != l[i-1]:
            cdf.append([s, i+1])
        else:
            cdf[-1][1] = i+1
    return cdf


def getCdfs(eid):
    cdfs = {}
    df = getRawData(EXP_PATH_PREFIX + str(eid))
    chrs = set(df['chr'])
    for c in chrs:
        start = df[df['chr']==c]['start']
        cdf = getCdf(start)
        cdfs[c] = cdf
    return cdfs


def test_snp_density():
    df = getRawData()
    #one chromosome
    df = df[df['chr']==df['chr'][0]]
    start = df['start']
    density = len(start) / float(start.max())
    print density


def test_cdf():
    df = getRawData()
    df = df[df['chr']==df['chr'][0]]
    start = df['start']
    print start
    cdf = getCdf(start)

    cdf = np.array(cdf)
    plt.step(cdf[:,0], cdf[:,1], where = 'mid')
    plt.show()

    
def test_cdfs():
    cdfs = getCdfs('10838')
    print cdfs.keys()

def test_scatterPlot():
    df = getRawData(EXP_PATH_PREFIX + '10838')
    df = df[df['chr']==df['chr'][0]]
    start1 = df['start']
    
    df = getRawData(EXP_PATH_PREFIX + '10813')
    df = df[df['chr']==df['chr'][0]]
    start2 = df['start']
    
    df = getRawData(EXP_PATH_PREFIX + '10841')
    df = df[df['chr']==df['chr'][0]]
    start3 = df['start']
    

    h1 = np.histogram(start1, bins=np.arange(0, np.max(start1), 20000))
    h2 = np.histogram(start2, bins=np.arange(0, np.max(start2), 20000))
    h3 = np.histogram(start3, bins=np.arange(0, np.max(start3), 20000))

    cm = plt.cm.get_cmap('viridis')


    plt.subplot(332)
    plt.plot(h2[0], h1[0], c='gray', alpha=0.1)
    plt.scatter(h2[0], h1[0], c=h1[1][:-1], cmap=cm)
    
    plt.subplot(333)
    plt.plot(h3[0], h1[0], c='gray', alpha=0.1)
    plt.scatter(h3[0], h1[0], c=h1[1][:-1], cmap=cm)
    
    plt.subplot(336)
    plt.plot(h3[0], h2[0], c='gray', alpha=0.1)
    plt.scatter(h3[0], h2[0], c=h1[1][:-1], cmap=cm)
    
    
    plt.show()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_scatterPlot()
