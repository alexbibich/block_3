import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

rawData = pd.read_csv('output/profiles.csv', encoding='windows-1251')
parametersNames = rawData.columns.tolist()[2:]
[timeLabel, coordLabel] = rawData.columns.tolist()[:2]
plotsCount = len(parametersNames) 

fig = plt.figure()
axes = [plt.subplot(plotsCount, 1, _ + 1) for _ in range(plotsCount)]

data_skip = len(set(rawData[coordLabel]))
time_moments = list(set(rawData[timeLabel]))

def init_func():
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel(coordLabel)
        axes[i].set_ylabel(parametersNames[i])
        coordData = rawData[coordLabel]
        paramData = rawData[parametersNames[i]]
        xLim = [min(coordData) - 0.1 * (max(coordData) - min(coordData)), max(coordData) + 0.1 * (max(coordData) - min(coordData))]
        yLim = [min(paramData) - 0.1 * (max(paramData) - min(paramData)), max(paramData) + 0.1 * (max(paramData) - min(paramData))]
        axes[i].set_xlim(xLim)
        axes[i].set_ylim(yLim)

def update(step):
    init_func()
    for i in range(len(parametersNames)):
        coordData = rawData[coordLabel]
        paramData = rawData[parametersNames[i]]
        axes[i].plot(coordData[step * data_skip: (step + 1) * data_skip], paramData[step * data_skip: (step + 1) * data_skip], color='b')

init_func()
update(0)

ax_time = plt.axes([0.15, 0.001, 0.65, 0.03])
time_slider = Slider(ax_time, 'Steps', 0, len(time_moments), valstep=1)

time_slider.on_changed(update)

#plt.subplots_adjust(bottom=0.12)
plt.show()


