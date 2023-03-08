# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import csv

x = []
y = []
u = []
v = []

with open('data/plot.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            (pos_x, pos_y, arr_x, arr_y) = (
                int(row[0]), int(row[1]), int(row[2]), int(row[3]))
            x.append(pos_x)
            y.append(pos_y)
            u.append(arr_x)
            v.append(arr_y)
            line_count += 1
    print(f'Processed {line_count} lines.')

colors = []
for index in range(len(x)):
    colors.append(u[index] * u[index] + v[index] * v[index])


norm = Normalize()
norm.autoscale(colors)

colormap = cm.plasma

# Shifting the values
xmax = max(x)
ymax = max(y)

for index in range(len(x)):
    x[index] = (x[index] + (xmax / 2)) % (xmax + 1)
    y[index] = (y[index] + (xmax / 3)) % (ymax + 1)

# Plotting Vector Field with quiver() function
plt.quiver(x, y, u, v, color=colormap(norm(colors)), scale_units='xy', scale=max((max(u), max(v))),
           width=0.01, headwidth=1, headlength=1)
plt.title('Vector Field')

# Setting boundary limits
plt.xlim(min(x), max(x))
plt.ylim(min(y), max(y))

# Show plot with grid
plt.grid()
plt.show()
