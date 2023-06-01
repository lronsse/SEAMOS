import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

facecolor = '#eaeaf2'
font_color = '#525252'
hfont = {'fontname': 'Calibri'}
colors = ['tab:red', 'tab:red',
          'tab:blue', 'tab:blue', 'tab:blue', 'tab:blue',
          'tab:green', 'tab:green', 'tab:green',
          'tab:orange', 'tab:orange']
custom_lines = [Line2D([0], [0], color='tab:orange', lw=4),
                Line2D([0], [0], color='tab:green', lw=4),
                Line2D([0], [0], color='tab:blue', lw=4),
                Line2D([0], [0], color='tab:red', lw=4)]


outer = ['Power', 'Cost', 'Reliability', 'Market']
inner = ['Navigation', 'Monitoring', 'System', 'Msintenance', 'Operation', 'R&D', 'Wave stability', 'Mechanisms',
         'Reliability', 'Time to market', 'Value']
values = [0.87, 0.21, 0.17, 1.23, 0.00, 2.33, 1.33, 0.33, 0.33, 1.00, 1.00]
# df = pd.DataFrame(values, index=outer, columns=inner)
plt.barh(inner, values, color=colors)
plt.legend(custom_lines, ['Market', 'Reliability', 'Cost', 'Power'])
plt.title('Grade variance')
#plt.axvline(1.5, color='black', linestyle='--')
plt.tight_layout()
'''
print(df)
df.plot(kind='bar', figsize=(10, 4))

ax = plt.gca()
pos = []
for bar in ax.patches:
    pos.append(bar.get_x() + bar.get_width() / 2.)

ax.set_xticks(pos, minor=True)
lab = []
for i in range(len(pos)):
    l = df.columns.values[i // len(df.index.values)]
    lab.append(l)

ax.set_xticklabels(lab, minor=True)
ax.tick_params(axis='x', which='major', pad=15, size=0)
plt.setp(ax.get_xticklabels(), rotation=0)'''

#plt.show()
plt.savefig(f'subcriteria_bar.png', dpi=360)
