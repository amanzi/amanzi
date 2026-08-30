#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def read_samples(filename):
    data = np.loadtxt(filename)

    x_data = data[:, 0]
    y_data = data[:, 1]

    return x_data, y_data 


def plot_samples(filename):
    X, Y = read_samples(filename)

    fig, ax = plt.subplots(figsize=(6, 7))

    ax.scatter(X, Y, color='black', s=0.05)
    ax.set_xlim(-50.0, 1100.0)
    ax.set_ylim(280.0, 950.0)

    # Read and plot phase boundary line
    line = np.loadtxt("phase_boundary.dat")

    x_line = line[:, 0]
    y_line = line[:, 1]

    ax.plot(x_line, y_line,
            color="black", linestyle="-", linewidth=2.0, label="special line")

    # Read and plot spinodal boundary line
    line = np.loadtxt("spinodal_boundary.dat")

    x_line = line[:, 0]
    y_line = line[:, 1]

    ax.plot(x_line, y_line,
            color="black", linestyle="-.", linewidth=2.0, label="special line")

    ax.set_xlabel("density")
    ax.set_ylabel("temperature")
    # ax.set_title("Isolines of f(x, y)")
    ax.grid(alpha=0.2)


    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_samples("samples.dat")

