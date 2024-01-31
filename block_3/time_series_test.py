import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob


plots = os.listdir('time_series/')

dfs = list()

for plot in plots:
    #if 'p_' in plot or 'Q_' in plot:
        
        
    if '.csv' in plot:
        dfs.append(pd.read_csv(f'time_series/{plot}', encoding='windows-1251'))
            

rawData = dfs
parametersNames = [data.columns.tolist()[2] for data in rawData]
coordLabel = [data.columns.tolist()[1] for data in rawData]
plotsCount = len(parametersNames) 

fig = plt.figure(figsize=(8, 5))

axes = [plt.subplot(plotsCount, 1, _ + 1) for _ in range(plotsCount)]

def init_func():
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel("Время")
        axes[i].set_ylabel(parametersNames[i])
        coordData = rawData[i][coordLabel[i]]
        paramData = rawData[i][parametersNames[i]]
        xLim = [min(coordData) - 0.1 * (max(coordData) - min(coordData)), max(coordData) + 0.1 * (max(coordData) - min(coordData))]
        yLim = [min(paramData) - 0.1 * (max(paramData) - min(paramData)), max(paramData) + 0.1 * (max(paramData) - min(paramData))]
        axes[i].set_xlim(xLim)
        axes[i].set_ylim(yLim)
    ini_draw()
        
plots = list()
def ini_draw(step=0):
    for i in range(len(parametersNames)):
        coordData = rawData[i][coordLabel[i]]
        paramData = rawData[i][parametersNames[i]]
        plots.append(axes[i].plot(coordData, paramData, color='b'))

init_func()

plt.subplots_adjust(left=0.06, bottom=0.06, right=0.971, top=0.973, 
         wspace=0.2, hspace=0.24)

plt.show()


