import numpy as np
import sympy as sp
from scipy.integrate import odeint
import festim as F

# common variables
defect_type_1_densities = [0.230, 0.230, 0.225, 0.153, 0.107]
defect_type_2_densities = [0.290, 0.290, 0.280, 0.280, 0.189]
defect_type_3_densities = [0.05, 0.05, 0.05, 0.05, 0.06]

trap_D1_densities = [0, 1.0e24, 3.5e24, 2.2e25, 4.8e25, 5.4e25, 5.5e25, 6.8e25]
trap_D2_densities = [0, 2.5e24, 5.0e24, 1.5e25, 3.8e25, 4.4e25, 4.6e25, 6.1e25]
trap_D3_densities = [0, 1.0e24, 2.5e24, 6.5e24, 2.6e25, 3.6e25, 4.0e25, 5.0e25]
trap_D4_densities = [0, 1.0e24, 1.9e24, 2.1e25, 3.6e25, 3.9e25, 4.5e25, 5.0e25]
trap_D5_densities = [0, 2.0e23, 1.6e24, 6.0e24, 1.1e25, 1.4e25, 1.7e25, 2.0e25]

# defect type I
A_0_1 = 6.1838e-03
E_A_1 = 0.24

trap_D1_K = 9.0e26
trap_D1_n_max = 6.9e25
trap_D2_K = 4.2e26
trap_D2_n_max = 7.0e25

# defect type II
A_0_2 = 6.1838e-03
E_A_2 = 0.30

trap_D3_K = 2.5e26
trap_D3_n_max = 6.0e25
trap_D4_K = 5.0e26
trap_D4_n_max = 4.7e25

# defect type III
A_0_3 = 0
E_A_3 = 1

trap_D5_K = 1.0e26
trap_D5_n_max = 2.0e25

k_B = F.k_B


def generate_fig_2_annealed_trap_fitting_data():
    """
    orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
    subsequently fitted by E.Hodille et al, available at https://doi.org/10.1088/1741-4326/aa5aa5
    """

    annealing_time = 7200

    # ##### standard variables ##### #
    t_annealing = np.linspace(0, annealing_time, num=1000)
    T_values = np.linspace(1, 800, num=1000)

    A_0_1 = 6.1834e-03
    E_A_1 = 0.24
    E_A_2 = 0.30

    # dummy values for no annealing
    A_0_2 = 0
    E_A_3 = 1

    # dummy values for damage parameters
    phi = 0
    K = 1
    n_max = 1

    n_0_defect_1, n_0_defect_2, n_0_defect_3 = (
        defect_type_1_densities[0],
        defect_type_2_densities[0],
        defect_type_3_densities[0],
    )

    (
        annealed_defect_1_densities,
        annealed_defect_2_densities,
        annealed_defect_3_densities,
    ) = ([], [], [])

    for T in T_values:
        # defect type 1
        defect_1_extra_args = (phi, K, n_max, A_0_1, E_A_1, T)
        n_traps_annleaing_defect_1 = odeint(
            neutron_trap_creation_numerical,
            n_0_defect_1,
            t_annealing,
            args=defect_1_extra_args,
        )
        annealed_defect_1_densities.append(n_traps_annleaing_defect_1[-1])

        # defect type 2
        defect_2_extra_args = (phi, K, n_max, A_0_1, E_A_2, T)
        n_traps_annleaing_defect_2 = odeint(
            neutron_trap_creation_numerical,
            n_0_defect_2,
            t_annealing,
            args=defect_2_extra_args,
        )
        annealed_defect_2_densities.append(n_traps_annleaing_defect_2[-1])

        # defect type 3
        defect_3_extra_args = (phi, K, n_max, A_0_2, E_A_3, T)
        n_traps_annleaing_defect_3 = odeint(
            neutron_trap_creation_numerical,
            n_0_defect_3,
            t_annealing,
            args=defect_3_extra_args,
        )
        annealed_defect_3_densities.append(n_traps_annleaing_defect_3[-1])

    # exporting
    np.savetxt("data/annealed_defect_1_densities.txt", annealed_defect_1_densities)
    np.savetxt("data/annealed_defect_2_densities.txt", annealed_defect_2_densities)
    np.savetxt("data/annealed_defect_3_densities.txt", annealed_defect_3_densities)


def generate_fig_4_TDS_fitting_data():
    # 0 dpa values
    print("running 0 dpa values")
    festim_sim(
        n1=0,
        n2=0,
        n3=0,
        n4=0,
        n5=0,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0/",
    )

    # 0.001 dpa values
    print("running 0.001 dpa values")
    festim_sim(
        n1=trap_D1_densities[1],
        n2=trap_D2_densities[1],
        n3=trap_D3_densities[1],
        n4=trap_D4_densities[1],
        n5=trap_D5_densities[1],
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.001/",
    )

    # 0.005 dpa values
    print("running 0.005 dpa values")
    festim_sim(
        n1=trap_D1_densities[2],
        n2=trap_D2_densities[2],
        n3=trap_D3_densities[2],
        n4=trap_D4_densities[2],
        n5=trap_D5_densities[2],
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.005/",
    )

    # 0.023 dpa values
    print("running 0.023 dpa values")
    festim_sim(
        n1=trap_D1_densities[3],
        n2=trap_D2_densities[3],
        n3=trap_D3_densities[3],
        n4=trap_D4_densities[3],
        n5=trap_D5_densities[3],
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.023/",
    )

    # 0.1 dpa values
    print("running case 0.1 dpa")
    festim_sim(
        n1=trap_D1_densities[4],
        n2=trap_D2_densities[4],
        n3=trap_D3_densities[4],
        n4=trap_D4_densities[4],
        n5=trap_D5_densities[4],
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.1/",
    )

    # 0.23 dpa values
    print("running 0.23 dpa values")
    festim_sim(
        n1=trap_D1_densities[5],
        n2=trap_D2_densities[5],
        n3=trap_D3_densities[5],
        n4=trap_D4_densities[5],
        n5=trap_D5_densities[5],
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.23/",
    )

    # 0.5 dpa values
    print("running 0.5 dpa values")
    festim_sim(
        n1=trap_D1_densities[6],
        n2=trap_D2_densities[6],
        n3=trap_D3_densities[6],
        n4=trap_D4_densities[6],
        n5=trap_D5_densities[6],
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.5/",
    )

    # 2.5 dpa values
    print("running 2.5 dpa values")
    festim_sim(
        n1=trap_D1_densities[7],
        n2=trap_D2_densities[7],
        n3=trap_D3_densities[7],
        n4=trap_D4_densities[7],
        n5=trap_D5_densities[7],
        results_foldername="testing/dpa_2.5/",
    )


def generate_fig_5_damaged_trap_fitting_data():
    """
    TDS data from T.Swartz-Selinger, currently unpublished
    """
    phi = 8.9e-05
    t_damage = int(3 / phi)
    t = np.linspace(0, t_damage, 1000)
    n_0 = 0
    T_damage = 800  # K

    # trap 1
    trap1_extra_args = (phi, trap_D1_K, trap_D1_n_max, A_0_1, E_A_1, T_damage)
    n_trap1_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap1_extra_args
    )
    # trap 2
    trap2_extra_args = (phi, trap_D2_K, trap_D2_n_max, A_0_1, E_A_1, T_damage)
    n_trap2_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap2_extra_args
    )
    # trap 3
    trap3_extra_args = (phi, trap_D3_K, trap_D3_n_max, A_0_2, E_A_2, T_damage)
    n_trap3_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap3_extra_args
    )
    # trap 4
    trap4_extra_args = (phi, trap_D4_K, trap_D4_n_max, A_0_2, E_A_2, T_damage)
    n_trap4_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap4_extra_args
    )
    # trap 5
    trap5_extra_args = (phi, trap_D5_K, trap_D5_n_max, A_0_3, E_A_3, T_damage)
    n_trap5_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap5_extra_args
    )

    # exporting
    np.savetxt("data/damage_trap_D1_fitting.txt", n_trap1_damaged)
    np.savetxt("data/damage_trap_D2_fitting.txt", n_trap2_damaged)
    np.savetxt("data/damage_trap_D3_fitting.txt", n_trap3_damaged)
    np.savetxt("data/damage_trap_D4_fitting.txt", n_trap4_damaged)
    np.savetxt("data/damage_trap_D5_fitting.txt", n_trap5_damaged)


def neutron_trap_creation_numerical(
    n, t, phi=9.64e-7, K=3.5e28, n_max=1e40, A_0=6.1838e-03, E_A=0.2792, T=298
):
    """
    Temporal evolution of n with resepct to time, for a given value of n and t

    Args:
        n (float): number of traps (m-3)
        t (float): time of simulation (s)
        phi (float): damage per second (dpa s-1). Defaults to 9.64e-07
        K (float): trap creation factor (traps dpa-1). Defaults to 1e28 traps
            s-1
        n_max (float): maximum traps per unit damage (m-3). Defaults to 1e40
            m-3
        A_0 (float): trap annealing factor (s-1). Defaults to 1e-02 s-1
        E_A (float): Annealing activation energy (eV). Defaults to 0.1034 eV
        T (float): the annealing temperature (K). Defaults to 298 K

    Returns:
        float: dn/dt
    """
    dndt = phi * K * (1 - (n / n_max)) - A_0 * np.exp(-E_A / (k_B * T)) * n
    return dndt


def festim_sim(
    n1=1,
    n2=1,
    n3=1,
    n4=1,
    n5=1,
    results_foldername="Results/",
):
    """Runs a FESTIM simulation with a custom mesh generator created with the
    automatic_vertices function.

    Args:
        n1 (float): trap_D1_density in m-3
        n2 (float): trap_D2_density in m-3
        n3 (float): trap_D3_density in m-3
        n4 (float): trap_D4_density in m-3
        n5 (float): trap_D5_density in m-3
        results_foldername (str) results folder location
    """
    fluence = 1.5e25
    implantation_time = 72 * 3600
    flux = fluence / implantation_time
    resting_time = 0.5 * 24 * 3600
    exposure_temp = 370
    resting_temp = 295
    ramp = 3 / 60
    tds_time = (1000 - 300) / ramp
    size = 8e-04
    atom_density_W = 6.3222e28

    center = 0.7e-9
    width = 0.5e-9
    distribution = (
        1 / (width * (2 * 3.14) ** 0.5) * sp.exp(-0.5 * ((F.x - center) / width) ** 2)
    )

    my_model = F.Simulation(log_level=40)

    # define materials
    tungsten = F.Material(
        id=1,
        borders=[0, size],
        D_0=1.6e-07,
        E_D=0.28,
    )
    my_model.materials = F.Materials([tungsten])

    # define traps
    damage_dist = 1 / (1 + sp.exp((F.x - 2.5e-06) / 5e-07))

    my_traps = [
        F.Trap(
            k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
            E_k=tungsten.E_D,
            p_0=1e13,
            E_p=1.04,
            density=2.4e22,
            materials=tungsten,
        ),
    ]

    if n1 == 0:
        my_model.traps = F.Traps(my_traps)
    else:
        my_traps.extend(
            [
                F.Trap(
                    k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                    E_k=tungsten.E_D,
                    p_0=1e13,
                    E_p=1.15,
                    density=n1 * damage_dist,
                    materials=tungsten,
                ),
                F.Trap(
                    k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                    E_k=tungsten.E_D,
                    p_0=1e13,
                    E_p=1.35,
                    density=n2 * damage_dist,
                    materials=tungsten,
                ),
                F.Trap(
                    k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                    E_k=tungsten.E_D,
                    p_0=1e13,
                    E_p=1.65,
                    density=n3 * damage_dist,
                    materials=tungsten,
                ),
                F.Trap(
                    k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                    E_k=tungsten.E_D,
                    p_0=1e13,
                    E_p=1.85,
                    density=n4 * damage_dist,
                    materials=tungsten,
                ),
                F.Trap(
                    k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                    E_k=tungsten.E_D,
                    p_0=1e13,
                    E_p=2.05,
                    density=n5 * damage_dist,
                    materials=tungsten,
                ),
            ]
        )
        my_model.traps = F.Traps(my_traps)

    # define mesh
    vertices = np.concatenate(
        [
            np.linspace(0, 3e-9, num=100),
            np.linspace(3e-9, 8e-6, num=100),
            np.linspace(8e-6, 8e-5, num=100),
            np.linspace(8e-5, 8e-4, num=100),
        ]
    )
    my_model.mesh = F.MeshFromVertices(vertices)

    # define temperature
    my_model.T = F.Temperature(
        value=sp.Piecewise(
            (exposure_temp, F.t < implantation_time),
            (
                resting_temp,
                (F.t >= implantation_time) & (F.t < implantation_time + resting_time),
            ),
            (
                300 + ramp * (F.t - (implantation_time + resting_time)),
                F.t >= implantation_time + resting_time,
            ),
            (0, True),
        )
    )

    # define boundary conditions
    my_model.boundary_conditions = [
        F.DirichletBC(surfaces=[1, 2], field="solute", value=0)
    ]

    # define sources
    my_model.sources = [
        F.Source(
            volume=1,
            field=0,
            value=sp.Piecewise(
                (flux * distribution, F.t < implantation_time),
                (0, True),
            ),
        )
    ]

    # define exports
    folder_results = results_foldername
    my_derived_quantities = F.DerivedQuantities(
        filename=folder_results + "last.csv",
        nb_iterations_between_exports=1,
    )

    average_T = F.AverageVolume("T", volume=1)
    H_flux_left = F.HydrogenFlux(surface=1)
    H_flux_right = F.HydrogenFlux(surface=2)
    solute = F.TotalVolume("solute", volume=1)
    retention = F.TotalVolume("retention", volume=1)
    trap_1 = F.TotalVolume("1", volume=1)
    my_derived_quantities.derived_quantities = [
        average_T,
        H_flux_left,
        H_flux_right,
        solute,
        retention,
        trap_1,
    ]
    if n1 != 0:
        my_derived_quantities.derived_quantities.extend(
            [
                F.TotalVolume("2", volume=1),
                F.TotalVolume("3", volume=1),
                F.TotalVolume("4", volume=1),
                F.TotalVolume("5", volume=1),
                F.TotalVolume("6", volume=1),
            ]
        )

    my_model.exports = F.Exports(
        [
            my_derived_quantities,
        ]
    )

    # define settings
    my_model.dt = F.Stepsize(
        1,
        stepsize_change_ratio=1.1,
        t_stop=implantation_time + resting_time * 0.5,
        dt_min=1e-1,
        stepsize_stop_max=50,
    )

    my_model.settings = F.Settings(
        absolute_tolerance=1e10,
        relative_tolerance=1e-10,
        final_time=implantation_time + resting_time + tds_time,
        transient=True,
        maximum_iterations=30,
    )

    my_model.initialise()
    my_model.run()


if __name__ == "__main__":
    # generate_fig_2_annealed_trap_fitting_data()
    generate_fig_4_TDS_fitting_data()
    # generate_fig_5_damaged_trap_fitting_data()
