import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(211, facecolor='#FFFFCC')

# fixme - load from the file saved
data = np.loadtxt(os.path.splitext(__file__)[0] + '.dat').T
time = data[1] / 1000 # to ns
rmsd = data[2]

ax.plot(time, rmsd, 'o--', label="all")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"RMSD ($\AA$)")
ax.legend(loc="best").set_draggable(True)


ax2 = fig.add_subplot(212)
ax2.hist(rmsd)


plt.show()

# PYMOL INTERAL: The following lines ensure that the interactive features in PyMOL work
pymol_x_axis = time
pymol_y_axis = rmsd
pymol_plot_ax = ax
pymol_hist_ax = ax2
# -------------------------------------------------------------------------------------