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
    ax = fig.add_subplot(111, facecolor='#FFFFCC')

    # fixme - load from the file saved
    data = np.loadtxt(os.path.splitext(__file__)[0] + '.dat').T
    resids = data[0]
    rmsf = data[1]

    ax.plot(resids, rmsf, 'o--', label="all")
    ax.set_xlabel("resids ")
    ax.set_ylabel(r"RMSF ($\AA$)")

    plt.show()

    return resids, rmsf, ax, fig


if __name__ == "__main__":
    # use the pylustrator only if it is asked for
    use_pylustrator = False
    if len(sys.argv) > 1 and sys.argv[1] == 'pylustrator':
        # now import pylustrator
        use_pylustrator = True
    plot(use_pylustrator)