import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


facecolor = '#eaeaf2'
font_color = '#525252'
hfont = {'fontname': 'Calibri'}
colors = ['tab:red', 'tab:red',
          'tab:blue', 'tab:blue', 'tab:blue', 'tab:blue',
          'tab:green', 'tab:green', 'tab:green',
          'tab:orange', 'tab:orange']

species = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10 - End of life profits')
penguin_means = {
    "Cost subscription": (-65008, -41407, -41407, -41407, -41407, -41407, -41407, -41407, -41407, -41407),
    "Potential subscription revenue": (54709, 54709, 54709, 54709, 54709, 54709, 54709, 54709, 54709, 54709),
    "Potential renting revenue": (48000, 48000, 48000, 48000, 48000, 48000, 48000, 48000, 48000, 48000,),
}
x = np.arange(len(species))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0
fig, ax = plt.subplots(layout='constrained')
for attribute, measurement in penguin_means.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

plt.title('Cash flows of subscription and renting models')
plt.xlabel('Year')
plt.ylabel('Euro [€]')
plt.axhline(color='black')
plt.ylim(-80000, 80000)

plt.show()
#plt.savefig(f'subcriteria_bar.png', dpi=360)
