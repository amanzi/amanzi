import numpy as np
import plotly.graph_objects as go

data = np.loadtxt("eos_iapws95_spline.dat")

rho = data[:, 0]
temperature = data[:, 1]
spline_phi = data[:, 2]

rho_unique = np.unique(rho)
temperature_unique = np.unique(temperature)

nrho = len(rho_unique)
nT = len(temperature_unique)

RHO = rho.reshape(nT, nrho)
TEMP = temperature.reshape(nT, nrho)
PHI = spline_phi.reshape(nT, nrho)

fig = go.Figure(
    data=[
        go.Surface(
            x=RHO,
            y=TEMP,
            z=PHI,
            connectgaps=False,
            contours={
                "z": {
                    "show": True,
                    "usecolormap": True,
                    "project_z": True,
                }
            },
            hovertemplate=(
                "rho = %{x:.6g} kg/m^3<br>"
                "T = %{y:.6g} K<br>"
                "phi_r = %{z:.8e}<extra></extra>"
            ),
        )
    ]
)

fig.update_layout(
    title="Splined residual Helmholtz energy",
    scene={
        "xaxis_title": "Density rho [kg/m^3]",
        "yaxis_title": "Temperature T [K]",
        "zaxis_title": "phi_r",
        "aspectmode": "auto",
        "camera": {
            "eye": {
                "x": 1.5,
                "y": 1.5,
                "z": 1.0,
            }
        },
    },
    width=1100,
    height=800,
)

fig.write_html(
    "helmholtz_interactive.html",
    include_plotlyjs=True,
)

fig.show()

