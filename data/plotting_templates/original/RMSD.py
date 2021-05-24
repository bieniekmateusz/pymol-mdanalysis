import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def plot(use_pylustrator):
    # use the pylustrator only if it is asked for
    if use_pylustrator:
        # import pylustrator
        import pylustrator
        # activate pylustrator
        pylustrator.start()

    fig = plt.figure(figsize=(8, 6))
    rmsd_ax = fig.add_subplot(211, facecolor='#FFFFCC')

    # fixme - load from the file saved
    data = np.loadtxt(os.path.splitext(__file__)[0] + '.dat').T
    time = data[1] / 1000  # to ns
    rmsd = data[2]

    rmsd_ax.plot(time, rmsd, 'o--', label="all")
    rmsd_ax.set_xlabel("time (ps)")
    rmsd_ax.set_ylabel(r"RMSD ($\AA$)")
    rmsd_ax.legend(loc="best").set_draggable(True)

    rmsd_hist_ax = fig.add_subplot(212)
    rmsd_hist_ax.hist(rmsd)

    plt.show()

    return time, rmsd, rmsd_ax, rmsd_hist_ax, fig


if __name__ == "__main__":
    # use the pylustrator only if it is asked for
    use_pylustrator = False
    if len(sys.argv) > 1 and sys.argv[1] == 'pylustrator':
        # now import pylustrator
        use_pylustrator = True
    plot(use_pylustrator)