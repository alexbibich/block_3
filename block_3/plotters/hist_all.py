import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


folders_path = '../research_out/Quasistationary/'
folders = [folder for folder in os.listdir(folders_path)]

data = pd.read_csv(diff_pout, encoding='windows-1251')
name = data.columns.tolist()[2]

bins_count = int((data[name].max() - data[name].min()) / 30000)

plt.hist(data[name], bins=bins_count, color='skyblue', edgecolor='black')
plt.grid(visible=True)
#Смещение: {data[name].skew()}\n
print(f'Медиана: {data[name].median()}\nСреднее: {data[name].mean()}')
print(f'СКО: {data[name].std()}')

plt.show()

