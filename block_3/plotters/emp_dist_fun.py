import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from statsmodels.distributions.empirical_distribution import ECDF


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

df = pd.read_csv(path + file + '/diff_press_pout.csv', encoding='windows-1251')
name = df.columns.tolist()[2]

ecdf = ECDF(df[name])
plt.plot(ecdf.x, ecdf.y)
plt.grid(visible=True)

plt.xlabel('x')
plt.ylabel('F(x)')

plt.axhline(y = 0.05, color = 'r', linestyle = '--') 
plt.axhline(y = 0.95, color = 'r', linestyle = '--') 

plt.show()
