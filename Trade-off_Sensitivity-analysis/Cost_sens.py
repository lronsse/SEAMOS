import matplotlib.pyplot as plt
import numpy as np
import copy
import math


def trade(x, w_t):
    return np.average(x, weights=w_t)


def generate_latex_table(row_titles, col_titles, values, weights):
    row_titles = [str(i) for i in row_titles]
    col_titles = [str(i) for i in col_titles]
    formatted_grades = []
    total = []
    for i in range(len(values)):
        total.append(trade(values[i], weights))
        vals = []
        for j in range(len(values[i])):
            if int(math.floor(values[i][j])) == 1:
                vals.append('\cellone{'+str(values[i][j])+'}')
            elif int(math.floor(values[i][j])) == 2:
                vals.append('\celltwo{'+str(values[i][j])+'}')
            elif int(math.floor(values[i][j])) == 3:
                vals.append('\cellthree{'+str(values[i][j])+'}')
            elif int(math.floor(values[i][j])) == 4:
                vals.append('\cellfour{'+str(values[i][j])+'}')
            elif int(math.floor(values[i][j])) == 5:
                vals.append('\cellfive{'+str(values[i][j])+'}')
        formatted_grades.append(vals)
    num_rows = len(row_titles)
    tot_width = 12.25
    width = [12.25 * weights[i] / 100 for i in range(len(weights))]
    width = ["{" + str(width[i]) + "cm}" for i in range(len(width))]
    width = [f"p{width[i]}|" for i in range(len(width))]
    # Begin LaTeX table
    latex_table = "\\begin{table}[H]" + "\n" + "\centering" + "\n" + "\caption{Fill in caption}" + "\n" + \
                  "\label{tab:label}" + "\n" + "\\resizebox{\\textwidth}{!}{\n"
    latex_table += "\\begin{tabular}{|l|" + ''.join(width) + "p{1.5cm}|}\\\ \hline\n"
    # Column titles
    latex_table += "Concepts & " + " & ".join(col_titles) + " & \\textbf{Total}" + "\\\ \hline\n"

    # Rows with values
    for i in range(num_rows):
        latex_table += row_titles[i] + " & " + " & ".join(formatted_grades[i]) + f"{total[i]}" + "\\\ \hline\n"

    # End LaTeX table
    latex_table += "\end{tabular}\n" + "}\n" + "\end{table}"
    return latex_table





def sens(v_puff, v_hyb, v_ms, w):
    w = [w[i] / 100 for i in range(len(w))]
    res_1, res_2, res_3 = [], [], []
    for i in range(len(w)):
        weights = copy.deepcopy(w)
        w_init = 0.4
        sub_1, sub_2, sub_3 = [], [], []

        while w_init >= 0:  # Iteration over weights & calculating final trade value for each weight
            for j in range(len(w)):
                weights[j] = w[j] / (1 - w[i]) * (1 - w_init)
            weights[i] = w_init
            res_hyb = trade(v_hyb, weights)
            res_puff = trade(v_puff, weights)
            res_ms = trade(v_ms, weights)
            sub_1.append(res_hyb)
            sub_2.append(res_puff)
            sub_3.append(res_ms)
            w_init -= 0.05
            w_init = round(w_init, 4)  # remove float innacuracies
        res_1.append(sub_1)
        res_2.append(sub_2)
        res_3.append(sub_3)
        results = [res_1, res_2, res_3]
    return results


criteria = ['Power', 'Cost', 'Reliability', 'Market']
weights = [30, 30, 25, 15]
systems = ['Puffin', 'Hybrid', 'Multi-system']
grades = [[3.06, 1.97, 2.97, 3], [2.18, 2.76, 2.98, 3.5], [2.23, 2.31, 2.75, 3.5]]
tab = generate_latex_table(systems, criteria, grades, weights)  # generate latex code for table
print(tab)
v_hyb = grades[1]
v_puff = grades[0]
v_ms = grades[2]

results = sens(v_hyb, v_puff, v_ms, weights)  # generate results of sensitivity analysis

colors = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue', 'tab:purple']
for i in range(len(criteria)):
    fig = plt.figure()
    for j in range(3):
        values = results[j][i]
        values.reverse()
        slope = (values[-1] - values[0]) / 40
        plt.plot(np.arange(0, 45, 5), values, color=colors[j], label=systems[j]+f' | {(slope * 100):2.2f}')
        # plt.vlines(0, values[-1] - (slope * 10), values[-1], color='black', linestyles='-', linewidth=1)
        # plt.hlines(values[-1] - (slope * 10), 0, 10, color='black', linestyles='-', linewidth=1)  # slope lines
    plt.title(label=criteria[i])
    plt.axvline(weights[i], label='Original value')
    plt.axvline((weights[i] + 10), color='black', linestyle='--', linewidth=1)
    plt.axvline((weights[i] - 10), color='black', linestyle='--', linewidth=1)
    plt.axvspan(weights[i] - 10, weights[i] + 10, facecolor='b', alpha=0.1)
    plt.xlabel(xlabel='Weight [%]')
    plt.ylabel(ylabel='Final value [-]')
    plt.legend()
    plt.grid()
    plt.show()
    # plt.savefig(f'{criteria[i]}_diagram.png', dpi=360)

