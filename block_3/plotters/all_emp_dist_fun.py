import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats


def find_prob_x(quant, prob, p, q):
    min_flag = True
    for i in range(len(prob)):
        if prob[i] >= q and min_flag:
            i_min = i
            min_flag = False
        if prob[i] > p + q:
            i_max = i - 1
            return quant[i_min], quant[i_max]
        

folders_path = '../research_out/Quasistationary/'
folders = [folder + '/' for folder in os.listdir(folders_path)]
filename = 'diff_press_pout.csv'

dfs = []
for folder in folders:
    dfs.append(pd.read_csv(folders_path + folder + filename, encoding='windows-1251'))

[time_moments, name_OX] = dfs[0].columns.tolist()[:2]

parameters_names = [df.columns.tolist()[2] for df in dfs]

x_right = dfs[0][parameters_names[0]].max()
x_left = dfs[0][parameters_names[0]].min()
y_bot = 0
y_top = 0
interval = 30000

axes = [plt.subplot(len(parameters_names)//2, len(parameters_names)//2, _ + 1) for _ in range(len(parameters_names))]

def find_y_top(df, name, n):
    df['bins'] = pd.cut(df[name], bins=n)
    freqs = df.groupby(['bins'])[name].count().reset_index()
    return freqs[name].max()

for name, df in zip(parameters_names, dfs):
    if x_right < df[name].max():
        x_right = df[name].max()
    if x_left > df[name].min():
        x_left = df[name].min()
    
def init_func():
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel('Погрешность давления, Па')
        axes[i].set_ylabel(parameters_names[i])
        axes[i].set_xlim(-200, 200)

def draw_fun():
    for i in range(len(parameters_names)):
        res = stats.ecdf(dfs[i][parameters_names[i]])
        axes[i].ecdf(dfs[i][parameters_names[i]])

        p = 0.90
        q = (1 - p) / 2

        [imin, imax] = find_prob_x(res.cdf.quantiles, res.cdf.probabilities, p, q)

        fsize = 11
        msize = 5
        axes[i].text(imin - 50, 0.08, f'{imin}', fontsize=fsize)
        axes[i].text(imin + 150, 0.90, f'{imax}', fontsize=fsize)
        axes[i].text(-180, 0.5, f'delta = {imax - imin}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
        axes[i].text(-190, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
        
        axes[i].axhline(y = q, color = 'r', linestyle = '--') 
        axes[i].axhline(y = p + q, color = 'r', linestyle = '--') 
        axes[i].plot(imin, q, 'bo', markersize=msize)
        axes[i].plot(imax, p + q, 'bo', markersize=msize)
        
       

init_func()

draw_fun()

plt.subplots_adjust(left=0.05, bottom=0.06, right=0.976, top=0.967, 
         wspace=0.2, hspace=0.136)

plt.show()

