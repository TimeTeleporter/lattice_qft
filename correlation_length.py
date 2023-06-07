import typing
from typing import List
import numpy as np
import csv
import math

# Folder where all the correlation function files are stored
CORRELATION_FUNCTION_PATH: str = 'data/plot_data/'
RESULTS_PATH: str = 'data/results.csv'

# Correlation lenghts to be calculated
indices = [0]


def calculate_correlation_length(correlation_function: List[float], max_t: int) -> float:
    p1: float = 2.0 * math.pi / float(max_t)
    g1: float = discrete_fourier_transform(correlation_function, p1)
    g2: float = discrete_fourier_transform(correlation_function, p1 * 2.0)
    cos1: float = math.cos(p1)
    cos2: float = math.cos(p1 * 2.0)
    cosh: float = (g1 * cos1 - g2 * cos2) / (g1 - g2)
    print(cosh)
    return math.acosh(cosh)


def discrete_fourier_transform(correlation_function: List[float], momentum: float) -> float:
    correlation_function = enumerate(correlation_function)
    correlation_function = map(
        lambda tup: tup[1] * math.cos(momentum * tup[0]), correlation_function)
    return sum(correlation_function)


for index in indices:
    # We read the corresponding files
    path: str = CORRELATION_FUNCTION_PATH + \
        'correlation_' + str(index) + '.csv'

    correlation_function = []

    # Reading the correlation functioin
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            correlation_function.append(float(row[0]))
            line_count += 1

    # Reading the simulation data
    max_t: int

    with open(RESULTS_PATH) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count = line_count
            elif int(row[0]) == index:
                # d = int(row[1])
                # size = int(row[2])
                # max_x = int(row[3])
                # max_y = int(row[4])
                max_t = int(row[5])
                # temp = float(row[6])
                # comptype = str(row[7])
                # action = float(row[8])
                # action_error: float = float(row[9])
                # energy_data = (row[10] == 'true')
                # difference_data = (row[11] == 'true')
                # difference_data = (row[12] == 'true')
            line_count += 1

    print(correlation_function)
    print(max_t)

    correlation_length = calculate_correlation_length(
        correlation_function, max_t)

    print(correlation_length)
