import typing
import numpy as np
import matplotlib.pyplot as plt
import csv

data_path: str = 'data/plot_data/'
images_path: str = 'images/'
results_path: str = 'data/results.csv'

indices = [5]

plot_difference: bool = True
plot_energy: bool = True
plot_correlation: bool = True

for index in indices:
    # We specify the index that we want to plot
    d: int
    size: int
    max_x: int
    max_y: int
    max_t: int
    temp: float
    comptype: str
    action: float
    # action_error: float
    energy_data: bool
    difference_data: bool
    correlation_data: bool

    # Then we read the result data of the given index
    with open(results_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                # print(f'Column names are {", ".join(row)}')
                line_count = line_count
            elif int(row[0]) == index:
                d = int(row[1])
                size = int(row[2])
                max_x = int(row[3])
                max_y = int(row[4])
                max_t = int(row[5])
                temp = float(row[6])
                comptype = str(row[7])
                action = float(row[8])
                # action_error: float = float(row[9])
                energy_data = (row[10] == 'true')
                difference_data = (row[11] == 'true')
                difference_data = (row[12] == 'true')
            line_count += 1

    plots = []
    if difference_data & plot_difference:
        plots.append('difference')
    if energy_data & plot_energy:
        plots.append('energy')
    if correlation_data & plot_correlation:
        plots.append('correlation')

    for plot in plots:
        if plot == 'difference' | plot == 'energy':
            # Next we plot the fields for each direction
            for direction in range(d):
                data: str = data_path + plot + '_' + \
                    str(index) + '_' + str(direction) + '.csv'
                image: str = images_path + plot + '_' + \
                    str(index) + '_' + str(direction) + '.png'
                print('Plotting ' + data + ' to ' + image)

                x = []
                y = []
                t = []
                v = []

                with open(data) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    line_count = 0
                    for row in csv_reader:
                        if line_count == 0:
                            # print(f'Column names are {", ".join(row)}')
                            line_count = line_count
                        else:
                            (pos_x, pos_y, pos_t, val) = (
                                int(row[0]), int(row[1]), int(row[2]), float(row[3]))
                            x.append(pos_x)
                            y.append(pos_y)
                            t.append(pos_t)
                            v.append(val)
                        line_count += 1

                # Check coherence
                assert max_x == max(x) + 1
                assert max_y == max(y) + 1
                assert max_t == max(t) + 1

                Z = np.zeros([max_x, max_y])

                for i in range(len(x)):
                    x_plot = int(y[i] + int(max_y / 2)) % max_y
                    y_plot = int(x[i] + int(max_x / 4)) % max_x
                    Z[x_plot, y_plot] += v[i]

                """
                for entry in Z:
                    entry = entry
                """

                fig, ax = plt.subplots()

                pos = ax.imshow(Z)
                ax.set_xlabel('x')
                ax.set_ylabel('y')

                fig.colorbar(pos, ax=ax)

                fig.savefig(image)
        if plot == 'correlation':
