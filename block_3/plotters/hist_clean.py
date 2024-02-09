import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


path = '../research_out/Quasistationary/'

plots = os.listdir(path)
while True:
    for i in range(len(plots)):
        print(f'{i+1}. {plots[i]}')
    try:
        choose = input("Выберите файл: ")
        if len(choose) != 1:
            print("неверный ввод")
            continue
        choice = int(choose) - 1
        file = plots[choice]
        break
    except:
        print("неверный ввод")

data = pd.read_csv(path + file + '/diff_press_pout.csv', encoding='windows-1251')
name = data.columns.tolist()[2]

data_clean = data[name].tolist()
data_clean.sort()

median = int((len(data_clean) - 1) / 2)

Q1 = int(median / 2)

Q3 = median + Q1

IQR = data_clean[Q3] - data_clean[Q1]

data_itog = []

for da in data_clean:
    if (da >= Q1 - 1.5 * IQR) and (da <= Q3 + 1.5 * IQR):
        data_itog.append(da)

bins_count = int((max(data_itog) - min(data_itog)) / 30000)

plt.hist(data_itog, bins=bins_count, color='skyblue', edgecolor='black')
plt.grid(visible=True)
#Смещение: {data[name].skew()}\n
print(f'Медиана: {data[name].median()}\nСреднее: {data[name].mean()}')
print(f'СКО: {data[name].std()}')

plt.show()

