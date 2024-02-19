import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
import shutil
from matplotlib.colors import LogNorm
import numpy as np

rcParams['text.usetex']= True if shutil.which('latex') else False


def plot_fig_6_trap_density_variation():
    k_B = 8.617e-5  # eV/K
    traps_properties = {
        "D1": {
            "A_0": 6.18e-3,
            "E_A": 0.24,
            "n_max": 6.8e25,
            "K": 1.9e28,
        },
        "D2": {
            "A_0": 6.18e-3,
            "E_A": 0.24,
            "n_max": 6.2e25,
            "K": 8.0e27,
        },
        "D3": {
            "A_0": 6.18e-3,
            "E_A": 0.30,
            "n_max": 5.3e25,
            "K": 3.5e27,
        },
        "D4": {
            "A_0": 6.18e-3,
            "E_A": 0.30,
            "n_max": 4.9e25,
            "K": 7.0e27,
        },
        "D5": {
            "A_0": 0,
            "E_A": 0,
            "n_max": 2.0e25,
            "K": 1.2e26,
        },
    }
    def A(T, A_0, E_A):
        return A_0 * np.exp(-E_A / k_B / T)
    def n(t, K, damage_rate, n_max, A):
        tau = (
            damage_rate * K / n_max + A
        ) ** -1  # characteristic time for trap creation/annealing
        n_infinity = (A / (damage_rate * K) + 1 / n_max) ** -1  # equilibrium density
        return n_infinity * (1 - np.exp(-t / tau))
    n_max = 6.8e25
    K = 4e27
    one_fpy = 1 * 365 * 24 * 3600
    colours = ["tab:blue", "tab:red"]
    for T, colour in zip([295, 800], colours):
        damage_rate_range = np.logspace(-3, 3) / one_fpy
        densities = 0
        for trap in traps_properties.values():
            A_0 = trap["A_0"]
            E_A = trap["E_A"]
            n_max = trap["n_max"]
            K = trap["K"]
            A_trap = A(T, A_0, E_A)
            densities_trap = n(one_fpy, K, damage_rate_range, n_max, A_trap)
            densities = densities + densities_trap
        plt.loglog(
            damage_rate_range * one_fpy, densities, label="{} K".format(T), color=colour
        )
    plt.xlabel("Damage rate (dpa/fpy)")
    plt.ylabel("Trap density (m$^{-3}$)")
    plt.legend(loc="lower right")
    plt.xlim(1e-03, 1e03)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # # add twin x axis for damage rate
    # ax1 = plt.gca()
    # ax2 = ax1.twiny()
    # ax2.set_xscale("log")
    # ax2.set_xlim([lim for lim in ax1.get_xlim()])
    # ax2.set_xlabel("Damage (dpa)")


def plot_fig_7_inventory_transient_and_distribution():
    dpa_values = np.geomspace(1e-05, 1e02, num=8)

    # ##### plotting ##### #
    # define colourbar
    norm = LogNorm(vmin=min(dpa_values), vmax=max(dpa_values))
    colorbar = cm.viridis
    sm = plt.cm.ScalarMappable(cmap=colorbar, norm=norm)
    colours = [colorbar(norm(dpa)) for dpa in dpa_values]

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=([6.4, 9.6]))

    # transient plot
    results_folder = "data/festim_model_results/transient_figure/"
    for dpa, colour in zip(dpa_values, colours):
        data_file = results_folder + "dpa={:.1e}/T=700/derived_quantities.csv".format(dpa)
        data = np.genfromtxt(data_file, delimiter=",", names=True)
        inv_per_dpa = data["Total_retention_volume_1"]
        times = data["ts"]
        axs[0].plot(times, inv_per_dpa, color=colour)

    # standard case
    standard_file = results_folder + "dpa=0.0e+00/T=700/derived_quantities.csv"
    data = np.genfromtxt(standard_file, delimiter=",", names=True)
    inv_per_dpa = data["Total_retention_volume_1"]
    times = data["ts"]
    axs[0].plot(times, inv_per_dpa, color="black", label=r"undamaged")

    day = 3600 * 24
    axs[0].vlines(day, 1e15, 1e23, color="black", alpha=0.5, linestyle="dashed")
    axs[0].annotate("24h", [day * 1.2, 1e16], color="black", alpha=0.5)

    axs[0].set_ylabel(r"T inventory (m$^{-2}$)")
    axs[0].set_xlabel(r"Time (s)")
    axs[0].set_ylim(1e15, 1e23)
    axs[0].set_xlim(left=1e-01)
    axs[0].set_yscale("log")
    axs[0].set_xscale("log")
    axs[0].legend()
    axs[0].spines["top"].set_visible(False)
    axs[0].spines["right"].set_visible(False)

    # ##### retention profile ##### #
    results_folder = "data/festim_model_results/profiles/"
    data_file = results_folder + "retention_profile_dpa=1.0e+02.csv"
    x_data = np.genfromtxt(data_file, delimiter=",", names=True)
    x_values = x_data["arc_length"]
    x_plot = x_values * 1000

    inventories = []
    for dpa in dpa_values:
        data_file = results_folder + "retention_profile_dpa={:.1e}.csv".format(dpa)
        data = np.genfromtxt(data_file, delimiter=",", names=True)
        inv = data["retention"]
        inventories.append(inv)

    for inv, dpa, colour in zip(inventories, dpa_values, colours):
        axs[1].plot(x_plot, inv, color=colour)

    # standard case
    standard_file = results_folder + "retention_profile_dpa=0.0e+00.csv"
    data = np.genfromtxt(standard_file, delimiter=",", names=True)
    inv = data["retention"]
    axs[1].plot(x_plot, inv, color="black")

    axs[1].set_xlabel(r"x (mm)")
    axs[1].set_ylabel(r"T retention (T m$^{-3}$)")
    axs[1].set_xlim(0, 2)
    axs[1].set_yscale("log")
    axs[1].set_ylim(inv[-1], 1e26)
    axs[1].spines["right"].set_visible(False)
    axs[1].spines["top"].set_visible(False)

    plt.tight_layout()
    fig.colorbar(sm, label=r"Damage rate (dpa/fpy)", ax=axs, shrink=0.5)


def plot_fig_8_inventory_variation():
    temperature_values = np.linspace(600, 1300, 50)
    dpa_values = np.geomspace(1e-05, 1e02, 8)

    inventories = []
    inventories_no_damage = []

    for dpa in dpa_values:
        inventory_per_dpa = []
        for T in temperature_values:
            data_file = "data/festim_model_results/dpa={:.1e}/T={:.0f}/derived_quantities.csv".format(dpa, T)
            data = np.genfromtxt(data_file, delimiter=",", names=True)
            
            inventory = data["Total_retention_volume_1"][-1]
            inventory_per_dpa.append(inventory)
        inventories.append(inventory_per_dpa)

    # undamaged values
    for T in temperature_values:
        data_file = "data/festim_model_results/dpa=0.0e+00/T={:.0f}/derived_quantities.csv".format(T)
        data = np.genfromtxt(data_file, delimiter=",", names=True)
        
        inventory = data["Total_retention_volume_1"][-1]
        inventories_no_damage.append(inventory)
        
    inventories = np.array(inventories)
    inventories_no_damage = np.array(inventories_no_damage)
    
    # normalise inventories
    normalised_inventories = []
    for dpa_case in inventories:
        norm = dpa_case/inventories_no_damage    
        normalised_inventories.append(norm)

    # colourbar
    norm = LogNorm(vmin=min(dpa_values), vmax=max(dpa_values))
    colorbar = cm.viridis
    sm = plt.cm.ScalarMappable(cmap=colorbar, norm=norm)
    colours = [colorbar(norm(dpa)) for dpa in dpa_values]


    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=([6.4, 9.6]))

    for case, colour in zip(inventories, colours):
        axs[0].plot(temperature_values, case, color=colour)
    axs[0].plot(temperature_values, inventories_no_damage, color="black")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(r"T Inventory (m-3)")
    axs[0].set_xlim(600, 1300)
    axs[0].spines["top"].set_visible(False)
    axs[0].spines["right"].set_visible(False)
    
    for case, colour in zip(normalised_inventories, colours):
        axs[1].plot(temperature_values, case, color=colour)
    axs[1].set_yscale("log")
    axs[1].set_ylabel(r"Relative T Inventory")
    axs[1].set_xlabel(r"Temperature (K)")
    axs[1].set_xlim(600, 1300)
    axs[1].spines["top"].set_visible(False)
    axs[1].spines["right"].set_visible(False)
    
    plt.tight_layout()
    fig.colorbar(sm, label=r"Damage rate (dpa/fpy)", ax=axs, shrink=0.5)


if __name__ == "__main__":
    
    plot_fig_6_trap_density_variation()
    plot_fig_7_inventory_transient_and_distribution()
    plot_fig_8_inventory_variation()
    
    plt.show()
