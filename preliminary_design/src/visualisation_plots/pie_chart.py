import matplotlib.pyplot as plt
import numpy as np

facecolor = '#eaeaf2'
font_color = '#525252'
hfont = {'fontname': 'Calibri'}
outer_labels = ['Cost', 'Power', 'Reliability', 'Market value']
soham_labels = ['Component', 'Maintenance', 'Operational']
soham_values = [5.39, 19.12, 75.49]
inner_labels = ['1', '2', '3', '1', '2', '3',
              '4', '5', '1', '2', '3', '4', '5',
              '1', '2', '1']


colors = ['#ffc09f', '#ffee93', '#adf7b6', '#a0ced9', '#fcf5c7']
#inner_colors = cmap([1, 2, 3, 5, 6, 7, 8, 9 , 10, 11, 12, 13, 14, 15, 16])

size = 0.3
outer_vals = [30, 30, 25, 15]
inner_vals = np.array([10.4, 10.4, 5.2, 2.6, 5.2, 5.2, 5.2, 7.8, 4, 6, 2, 2, 6, 7.2, 10.8, 10])
plt.pie(soham_values, radius=1, colors=['#ffc09f', '#adf7b6', '#a0ced9'], pctdistance=0.6, autopct='%1.1f%%', labeldistance=1.05, labels=soham_labels, wedgeprops=dict(edgecolor='grey'))
#        wedgeprops=dict(width=size, edgecolor='w')
#plt.pie(inner_vals, radius=1-size, colors=inner_colors,
#       wedgeprops=dict(width=size, edgecolor='w'), labels=inner_labels, labeldistance=0.7)
#plt.show()
plt.savefig(f'Costs_pie.png', dpi=360)
