from flask import Flask
from flask import render_template, url_for

from db import getCdfs
import json
import numpy as np
app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/exp/<int:eid>.json')
def exp(eid):
    cdfs = getCdfs(eid)
    return json.dumps(cdfs)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    res = exp(10838, 2000)

