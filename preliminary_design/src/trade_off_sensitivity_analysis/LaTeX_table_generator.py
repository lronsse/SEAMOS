def generate_latex_table(row_titles, col_titles, values, weights):
    row_titles = [str(i) for i in row_titles]
    col_titles = [str(i) for i in col_titles]
    for i in range(len(values)):
        for j in range(len(values[i])):
            print(type(values[i][j]))
            if round(values[i][j]) == 1:
                values[i][j] = '\cellone{'+str(values[i][j])+'}'
            if round(values[i][j]) == 2:
                values[i][j] = '\celltwo{'+str(values[i][j])+'}'
            if round(values[i][j]) == 3:
                values[i][j] = '\cellthree{'+str(values[i][j])+'}'
            if round(values[i][j]) == 4:
                values[i][j] = '\cellfour{'+str(values[i][j])+'}'
            if round(values[i][j]) == 5:
                values[i][j] = '\cellfive{'+str(values[i][j])+'}'
    num_rows = len(row_titles)
    tot_width = 12.25
    width = [12.25 * weights[i] / 100 for i in range(len(weights))]
    width = ["{" + str(width[i]) + "cm}" for i in range(len(width))]
    width = [f"p{width[i]}|" for i in range(len(width))]
    # Begin LaTeX table
    latex_table = "\\begin{table}[H]" + "\n" + "\centering" + "\n" + "\caption{Fill in caption}" + "\n" + \
                  "\label{tab:label}" + "\n" + "\\resizebox{\\textwidth}{!}{\n"
    latex_table += "\\begin{tabular}{|l|" + ''.join(width) + "}p{1.5cm}\\hline\n"
    # Column titles
    latex_table += "Concepts & " + " & ".join(col_titles) + " \\\\\n"
    latex_table += "\\hline\n"

    # Rows with values
    for i in range(num_rows):
        latex_table += row_titles[i] + " & " + " & ".join(values[i]) + "\\hline\n"

    # End LaTeX table
    latex_table += "\end{tabular}\n" + "}\n" + "\end{table}"



    return latex_table


criteria = ['Power', 'Cost', 'Reliability', 'Market']
weights = [30, 30, 25, 15]
systems = ['Puffin', 'Hybrid', 'Multi-system']
grades = [[5.2, 5.3, 5, 5], [5, 5, 5, 5], [5, 5, 5, 5], [5, 5, 5, 5]]

tab = generate_latex_table(systems, criteria, grades, weights)

print(tab)
