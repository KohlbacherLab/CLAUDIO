import matplotlib.pyplot as plt
import numpy as np


def create_homo_signal_histograms(data, filename, output_directory):
    # create histograms for homology signal, e.g. for adjacency of interaction sites and peptide overlaps
    # between/including interacting residues
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    # adjacency histogram
    adj_data = data["homo_adjacency"].copy()
    bins = [round(x * 0.1, 1) for x in range(11)]
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    freq, _, patches = plt.hist(adj_data, bins=bins, alpha=.3)
    plt.xlabel(f"relative adjacency")
    bin_centers = np.diff(bins) * .5 + bins[:-1]
    i = 0
    for fr, x, patch in zip(freq, bin_centers, patches):
        height = int(freq[i])
        plt.annotate(str(height) if height > 0 else '', xy=(x, height), xytext=(0, .2),
                     textcoords="offset points", ha="center", va="bottom")
        i += 1
    plt.ylabel("frequency")
    x_labels = [str(x) for x in bins]
    plt.xticks(bins, x_labels)
    plt.title("Frequency of relative interaction site adjacencies " 
              "$\left[adj = 1-\\frac{|pos_a - pos_b|}{len_{seq}}\\right]$")
    plt.savefig(f"{output_directory}{filename}_adj_hist.png")

    # Clear figure
    plt.clf()

    # histogram of peptide overlaps between/including interacting residues
    overl_data = data["homo_int_overl"].copy()
    bins = [round(x * 0.1, 1) for x in range(11)]
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    freq, _, patches = plt.hist(overl_data, bins=bins, alpha=.3)
    bin_centers = np.diff(bins) * .5 + bins[:-1]
    i = 0
    for fr, x, patch in zip(freq, bin_centers, patches):
        height = int(freq[i])
        plt.annotate(str(height) if height > 0 else '', xy=(x, height), xytext=(0, .2),
                     textcoords="offset points", ha="center", va="bottom")
        i += 1
    plt.xlabel(f"relative overlap")
    plt.ylabel("frequency")
    x_labels = [str(x) for x in bins]
    plt.xticks(bins, x_labels)
    plt.title("Relative peptide overlaps (" + '$n_{\\neg 0}$' + f" = {len(overl_data[overl_data > 0].index)})")
    plt.savefig(f"{output_directory}{filename}_int_ovl.png")

    # Clear figure
    plt.clf()

    # histogram of peptide overlaps between/including interacting residues
    overl_data = data["homo_pep_overl"].copy()
    bins = [0, .5,  1]
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    freq, _, patches = plt.hist(overl_data.astype(int), bins=bins, alpha=.3)
    bin_centers = np.diff(bins) * .5 + bins[:-1]
    i = 0
    for fr, x, patch in zip(freq, bin_centers, patches):
        height = int(freq[i])
        plt.annotate(str(height) if height > 0 else '', xy=(x, height), xytext=(0, .2),
                     textcoords="offset points", ha="center", va="bottom")
        i += 1
    plt.ylabel("frequency")
    x_labels = [False, True]
    plt.xticks(bin_centers, x_labels)
    plt.title("Boolean for whether peptides overlap")
    plt.savefig(f"{output_directory}{filename}_pep_ovl.png")