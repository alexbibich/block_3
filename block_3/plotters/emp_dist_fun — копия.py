import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats

def find_prob_x(quant, prob):
    min_flag = True
    for i in range(len(prob)):
        if prob[i] >= 0.05 and min_flag:
            i_min = i
            min_flag = False
        if prob[i] > 0.95:
            i_max = i - 1
            return quant[i_min], quant[i_max]
        
        

path = '../research_out/Quasistationary/'

plots = os.listdir(path)
while True:
    for i in range(len(plots)):
        print(f'{i+1}. {plots[i]}')
    try:
        #choose = input("Выберите файл: ")
        choose = '1'
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

res = stats.ecdf(df[name])
res.cdf.plot()
plt.grid(visible=True)

[imin, imax] = find_prob_x(res.cdf.quantiles, res.cdf.probabilities)
print(f'{imin} --- {imax}')

plt.xlabel('x')
plt.ylabel('F(x)')

plt.text(imin - 100, 0.08, f'{imin}', fontsize=9)
plt.text(imin, 0.05, '.', fontsize=60)
plt.text(imin + 150, 0.90, f'{imax}', fontsize=9)

plt.axhline(y = 0.05, color = 'r', linestyle = '--') 
plt.axhline(y = 0.95, color = 'r', linestyle = '--') 

plt.show()
