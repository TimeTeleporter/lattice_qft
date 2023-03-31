import numpy as np
import matplotlib.pyplot as plt
import csv

index: int = 1

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

x_max = max(x)
