import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


diff_pout = 'output/diff_press_pout.csv'

data = pd.read_csv(diff_pout, encoding='windows-1251')
name = data.columns.tolist()[2]

bins_count = int((data[name].max() - data[name].min()) / 10000)

plt.hist(data[name], bins=bins_count, color='skyblue', edgecolor='black')
plt.grid(visible=True)
print(f'Смещение: {data[name].skew()}\nМедиана: {data[name].median()}\nСреднее: {data[name].mean()}')
print(f'СКО: {data[name].std()}')

plt.show()
