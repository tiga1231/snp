import subprocess as sp
import pandas as pd
import numpy as np
from natsort import natsorted

def getAllDataFromExp(dirName = 'experiment_10838', limit=9e15):
    x = sp.Popen('ibis -d %s -q "select ' % dirName+\
                ' chr,start,stop,type,ref,alt ' +\
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



