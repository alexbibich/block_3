import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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

test = [0, 2, 2, 6]

sns.kdeplot(test, cumulative=True)

