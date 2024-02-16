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

res = stats.ecdf(df[name])
#res.cdf.plot()
plt.ecdf(df[name])
plt.grid(visible=True)

p = 0.90
q = (1 - p) / 2

[imin, imax] = find_prob_x(res.cdf.quantiles, res.cdf.probabilities, p, q)


plt.xlabel('x')
plt.ylabel('F(x)')
fsize = 11
msize = 5
plt.text(imin - 50, 0.08, f'{imin}', fontsize=fsize)
plt.text(imin + 150, 0.90, f'{imax}', fontsize=fsize)
plt.text(-180, 0.5, f'delta = {imax - imin}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
plt.text(-190, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})

plt.xlim(-200, 200)

plt.axhline(y = q, color = 'r', linestyle = '--') 
plt.axhline(y = p + q, color = 'r', linestyle = '--') 
plt.plot(imin, q, 'bo', markersize=msize)
plt.plot(imax, p + q, 'bo', markersize=msize)

plt.show()
