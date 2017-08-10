from flask import Flask
from flask import render_template, url_for

from db import getAllDataFromExp as getData
import json
import numpy as np
app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/exp/<int:eid>_<int:binsize>.json')
def exp(eid, binsize):
    df = getData('experiment/experiment_%d' % eid)
    chrs = df['chr'].unique()
    res = {}
    for c in chrs:
        start = df[df['chr']==c]['start']
        hist, bin_edges = np.histogram(start, bins=range(0,start.max(), binsize))
        res[c] = {  'hist':list(hist), 
                    'bin_edges':list(bin_edges) }
    return json.dumps(res)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    res = exp(10838, 2000)

