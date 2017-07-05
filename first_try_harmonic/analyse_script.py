import matplotlib.pyplot as plt
import numpy as np

PATH_TO_DATA = "/home/ilyas/space/data/tmp/harmonic/"

data = []
data.append(np.loadtxt(PATH_TO_DATA + "out_0"))
data.append(np.loadtxt(PATH_TO_DATA + "out_1"))
data.append(np.loadtxt(PATH_TO_DATA + "out_2"))
data.append(np.loadtxt(PATH_TO_DATA + "out_3"))


f, ax = plt.subplots(4, 1, sharex=True)
for i in range(4):
    data[i][data[i][:, 1] > 1e12, 1] = 1e12
    ax[i].plot(data[i][:, 0], data[i][:, 1])

f.show()
