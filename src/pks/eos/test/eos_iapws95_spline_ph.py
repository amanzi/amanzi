import numpy as np
import plotly.graph_objects as go

data = np.loadtxt("eos_iapws95_spline_ph.dat")

pressure = data[:, 0]
enthalpy = data[:, 1]
spline_s = data[:, 2]

p_unique = np.unique(pressure)
h_unique = np.unique(enthalpy)

np = len(p_unique)
nh = len(h_unique)

P = pressure.reshape(nh, np)
H = enthalpy.reshape(nh, np)
S = spline_s.reshape(nh, np)

fig = go.Figure(
    data=[
        go.Surface(
            x=P,
            y=H,
            z=S,
            connectgaps=False,
            contours={
                "z": {
                    "show": True,
                    "usecolormap": True,
                    "project_z": True,
                }
            },
            hovertemplate=(
                "p = %{x:.6g} MPa<br>"
                "h = %{y:.6g} kJ<br>"
                "entropy = %{z:.8g}<extra></extra>"
            ),
        )
    ]
)

fig.update_layout(
    title="Splined entropy",
    scene={
        "xaxis_title": "Pressure p [MPa]",
        "yaxis_title": "Enthalpy h [kJ]",
        "zaxis_title": "entropy",
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
    "entropy_interactive.html",
    include_plotlyjs=True,
)

fig.show()

