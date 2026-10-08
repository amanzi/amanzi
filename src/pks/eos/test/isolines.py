#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def read_grid(filename):
    data = np.loadtxt(filename)

    x_data = data[:, 0]
    y_data = data[:, 1]
    f_data = data[:, 2] # 3 - spline

    x = np.unique(x_data)
    y = np.unique(y_data)

    nx = len(x)
    ny = len(y)

    if nx * ny != len(data):
        raise ValueError(
            "Input data do not form a complete rectangular grid."
        )

    z = np.full((ny, nx), np.nan)

    # This tracks whether a grid point was actually supplied
    present = np.zeros((ny, nx), dtype=bool)

    x_index = {value: i for i, value in enumerate(x)}
    y_index = {value: j for j, value in enumerate(y)}

    for xi, yi, fi in zip(x_data, y_data, f_data):
        i = x_index[xi]
        j = y_index[yi]

        z[j, i] = fi
        present[j, i] = True

    # Check only whether coordinates are missing,
    # not whether f(x,y) happens to be NaN.
    if not np.all(present):
        raise ValueError("Some grid points are missing.")

    X, Y = np.meshgrid(x, y)

    return X, Y, z


def filtered_levels(zmin, zmax, n_candidates=50, min_fraction=0.06, power=2):
    # Work in transformed coordinates where polynomial growth is flattened
    a = zmin**(1.0 / power)
    b = zmax**(1.0 / power)

    candidates_t = np.linspace(a, b, n_candidates)

    kept_t = [candidates_t[0]]

    min_spacing = min_fraction * (b - a)

    for v in candidates_t[1:]:
        if v - kept_t[-1] >= min_spacing:
            kept_t.append(v)

    return np.asarray(kept_t)**power


def filtered_levels_exp(zmin, zmax, n_candidates=50, min_fraction=0.06):
    a = np.log(zmin)
    b = np.log(zmax)

    candidates_t = np.linspace(a, b, n_candidates)

    kept_t = [candidates_t[0]]

    min_spacing = min_fraction * (b - a)

    for v in candidates_t[1:]:
        if v - kept_t[-1] >= min_spacing:
            kept_t.append(v)

    return np.exp(np.asarray(kept_t))


def plot_isolines(filename):
    X, Y, Z = read_grid(filename)
    zmin = Z.min()
    zmax = Z.max()
    print(zmin, zmax)

    Z = abs(Z);
    print("Protting absolute value of the field...")

    Z_masked = np.ma.masked_invalid(Z)

    zmin = Z_masked.min()
    zmax = Z_masked.max()

    # levels = np.geomspace(zmin, zmax, 40)
    levels = filtered_levels(zmin, zmax, n_candidates=100, min_fraction=0.02, power=4)
    # levels = filtered_levels_exp(zmin, zmax, n_candidates=100, min_fraction=0.02)

    fig, ax = plt.subplots(figsize=(6, 7))

    # Choose contour levels.
    # Passing an integer asks Matplotlib to choose approximately
    # this many useful levels.
    contours = ax.contour(X, Y, Z_masked,
                          levels=levels,
                          linewidths=0.8, linestyles="solid", colors="black")

    # Print values directly on isolines.
    ax.clabel(contours, inline=True, fontsize=9, fmt="%.2g")

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
            color="black", linestyle="dotted", linewidth=2.0, label="special line")

    ax.set_xlabel("pressure")
    ax.set_ylabel("enthalpy")
    #ax.set_xlabel("density")
    #ax.set_ylabel("temperature")
    # ax.set_title("Isolines of f(x, y)")
    ax.grid(alpha=0.2)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_isolines("field.dat")

