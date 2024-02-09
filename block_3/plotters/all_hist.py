import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


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
axes = [plt.subplot(len(parameters_names), 1, _ + 1) for _ in range(len(parameters_names))]

def find_y_top(df, name, n):
    df['bins'] = pd.cut(df[name], bins=n)
    freqs = df.groupby(['bins'])[name].count().reset_index()
    return freqs[name].max()

for name, df in zip(parameters_names, dfs):
    if x_right < df[name].max():
        x_right = df[name].max()
    if x_left > df[name].min():
        x_left = df[name].min()

    bins_count = int((df[name].max() - df[name].min()) / interval)
    freq = find_y_top(df, name, bins_count)
    if freq > y_top:
        y_top = freq
    
def init_func():
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel(name_OX)
        axes[i].set_ylabel(parameters_names[i])
        axes[i].set_xlim(x_left, x_right)
        axes[i].set_ylim(y_bot, y_top + 100)

def draw_fun():
    for i in range(len(parameters_names)):
        bins_count = int((dfs[i][parameters_names[i]].max() - dfs[i][parameters_names[i]].min()) / interval)
        top_pos = y_top - 400
        axes[i].hist(dfs[i][parameters_names[i]], bins=bins_count, color='skyblue', edgecolor='black')
        axes[i].text(x_right - 200000, top_pos, f'СКО: {dfs[i][parameters_names[i]].std():.4f}', fontsize=12, bbox={'facecolor': 'white', 'alpha': 1})

init_func()

draw_fun()

plt.subplots_adjust(left=0.05, bottom=0.06, right=0.976, top=0.967, 
         wspace=0.2, hspace=0.136)

plt.show()

#data_clean = data[name].tolist()
#data_clean.sort()

#median = int((len(data_clean) - 1) / 2)

#Q1 = int(median / 2)

#Q3 = median + Q1

#IQR = data_clean[Q3] - data_clean[Q1]

#data_itog = []

#for da in data_clean:
    #if (da >= Q1 - 1.5 * IQR) and (da <= Q3 + 1.5 * IQR):
        #data_itog.append(da)