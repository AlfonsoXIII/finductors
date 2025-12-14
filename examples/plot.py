import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

data = np.genfromtxt("bfield.csv", delimiter=",")

x  = data[:, 0]
z  = data[:, 2]
Bx = data[:, 3]
By = data[:, 4]
Bz = data[:, 5]

Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
Bmag[Bmag <= 0.0] = np.min(Bmag[Bmag > 0.0])

xmin, xmax = x.min(), x.max()
zmin, zmax = z.min(), z.max()

plt.figure(figsize=(7, 6))

sc = plt.scatter(
    x, z,
    c=Bmag,
    cmap="viridis",
    norm=LogNorm(vmin=Bmag.min(), vmax=Bmag.max()),
    s=40
)

plt.colorbar(sc, label="|B| (log)")
plt.xlabel("X")
plt.ylabel("Z")
plt.title("Campo magnético en el plano XZ (escala log)")

plt.xlim(xmin, xmax)
plt.ylim(zmin, zmax)

plt.margins(0.0)
plt.tight_layout()
plt.show()

