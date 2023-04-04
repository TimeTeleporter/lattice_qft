import typing
import numpy as np
import matplotlib.pyplot as plt
import csv

difference_path = 'data/plot_data/difference_'
energy_path = 'data/plot_data/energy_'

results_path: str = 'data/results.csv'

indices = [0, 1]

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
    bonds_data: bool

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
                bonds_data = (row[11] == 'true')
            line_count += 1

    paths = []
    if bonds_data:
        paths.append(difference_path)
    if energy_data:
        paths.append(energy_path)

    for obs in paths:
        # Next we plot the fields for each direction
        for direction in range(d):
            path: str = obs + str(index) + '_' + str(direction) + '.csv'
            print('Printing to ' + path)

            x = []
            y = []
            t = []
            v = []

            with open(path) as csv_file:
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
                #Z[x[i], y[i]] += x[i]
                Z[x[i], y[i]] += v[i]

            for entry in Z:
                entry = entry

            fig, ax = plt.subplots()

            ax.imshow(Z)
            ax.set_xlabel('x')
            ax.set_ylabel('y')

            fig.savefig('images/difference_' + str(index) +
                        '_' + str(direction) + '.png')
