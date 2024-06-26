import numpy as np
for i,D in enumerate([10,30]):
    for f in range(1,30):
        filename = f"L-SRTDE_F{f}_D{D}.txt"
        X = np.loadtxt(filename)
        np.savetxt(filename,X.T)