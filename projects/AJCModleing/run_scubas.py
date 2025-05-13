import datetime as dt
import os

from cable import SCUBASModel
from fetch_cable_routes import (
    calculate_bathymetry_byLITHO1,
    calculate_conductive_profiles_with_distance,
    get_cable_route,
    plot_bathymatry,
    plot_routes,
)
from fetch_data import clean_B_fields, load_speadas
from loguru import logger
from parse_supermag import get_nearest_station
from utils import create_from_lat_lon


def run_May2024_storm():
    os.makedirs("simulation/May2024/", exist_ok=True)
    dates = [dt.datetime(2024, 5, 10, 12), dt.datetime(2024, 5, 12)]
    themis_fgm, themis_mom = load_speadas(dates)

    # Create dSegments and Profiles
    # 1. From cable routes run following methods
    #   a. get_cable_route
    #   b. calculate_bathymetry_byLITHO1
    #   c. calculate_conductive_profiles
    # This will help to check how d, profiles will help computing segment files
    o = get_cable_route()
    dpths_pts = calculate_bathymetry_byLITHO1(o)
    # 2. Manually check and assign the segment files to each profile
    # Probably remove some points in 'd' and then send for profilings and check via plots
    # Ids = [0,2,4,6,8,10,14,16,20,31]
    Id_st, Id_nt = (
        [0, 2, 4, 6, 8, 10, 14, 16, 20, 31],
        [0, 0, 2, 4, 6, 8, 10, 14, 16, 20],
    )
    plot_routes(o, dpths_pts.iloc[Id_st])
    profiles = []
    for ids, idn in zip(Id_st, Id_nt):
        ids = [ids, ids]
        dp = dpths_pts.iloc[ids]
        dpn = dpths_pts.iloc[idn]
        profiles.extend(calculate_conductive_profiles_with_distance(dp, dpn))
    plot_bathymatry(profiles)

    if True:
        # 2.b assign nearest mag stations and segment files
        segment_files, dSegments, to_profiles = (
            [],
            [],
            [p["profile"] for p in profiles],
        )
        for _, d in dpths_pts.iterrows():
            _, file = get_nearest_station(d.geolats, d.geolongs)
            segment_files.append([file])
            dSegments.append([d.geolats, d.geolongs])

        # 3. Conduuct a mag data cleaining by invoking  clean_B_fields
        #       (will make sure files are readable)
        clean_B_fields(
            ["CNB", "CTA", "GUA", "KAK"],
            [
                ["datasets/IM.cnb.min.txt"],
                ["datasets/IM.cta.min.txt"],
                ["datasets/IM.gua.min.txt"],
                ["datasets/IM.kak.min.txt"],
            ],
        )
        # There will be -1 in profiles which is correct number.
        logger.info(
            f" Len(files): {len(segment_files)} Len(segments): {len(dSegments)} Len(profiles): {len(to_profiles)}"
        )

        # 4. Run SCUBASModel for specific AJC cable
        #   Note: Make sure you have modified get_cable_informations
        cable = create_from_lat_lon(
            dSegments[:-1],
            to_profiles,
            left_active_termination=None,
            right_active_termination=None,
        )
        model = SCUBASModel(cable_structure=cable, segment_files=segment_files)
        model.initialize_TL()
        # 5. Run the folloding lines to finally run the code and save (plot).
        model.run_cable_segment()
        data = model.cable.tot_params.reset_index().copy()
        data.to_csv(
            "datasets/SCUBAS-AJC-4-Stations.csv", float_format="%g", index=False
        )

        # 6. Plotting the datasets
        import matplotlib.dates as mdates
        from plots import TimeSeriesPlot

        ts = TimeSeriesPlot(
            dates,
            major_locator=mdates.HourLocator(byhour=range(0, 24, 6)),
            minor_locator=mdates.HourLocator(byhour=range(0, 24, 3)),
            fig_title="AJC Simulated Voltage",
            text_size=15,
            num_subplots=1,
        )
        data.set_index("Time", inplace=True)
        ax = ts.add_voltage(
            data, ylim=[-200, 200], xlabel="Hours Since 10 May 2024, 12 UT"
        )
        ts.save("simulation/May2024/SCUBAS-Simulation-Compare.png")
        ts.close()
    return


def plot_magnetometer_datasets():
    import matplotlib.dates as mdates
    from plots import TimeSeriesPlot

    stns = ["CNB", "CTA", "GUA", "KAK"]

    datasets = [
        clean_B_fields([s], [[f"datasets/IM.{s.lower()}.min.txt"]])[s] for s in stns
    ]

    dates = [dt.datetime(2024, 5, 10, 12), dt.datetime(2024, 5, 12)]
    ts = TimeSeriesPlot(
        dates,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 6)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 3)),
        fig_title="Observations from magnetometers",
        text_size=15,
        num_subplots=3,
    )
    ax = ts._add_axis()
    for d, col, stn in zip(datasets, ["r", "g", "k", "b"], stns):
        ax.plot(d.index, d.X, label=f"{stn}", color=col, ls="-", lw=0.8)
    ax.legend(loc=2)
    ax.set_ylim(-500, 500)
    ax.set_ylabel("N-S components, nT")

    ax = ts._add_axis()
    for d, col, stn in zip(datasets, ["r", "g", "k", "b"], stns):
        ax.plot(d.index, d.Y, label=f"{stn}", color=col, ls="-", lw=0.8)
    ax.set_ylim(-500, 500)
    ax.set_ylabel("E-W components, nT")

    ax = ts._add_axis()
    for d, col, stn in zip(datasets, ["r", "g", "k", "b"], stns):
        ax.plot(d.index, d.Z, label=f"{stn}", color=col, ls="-", lw=0.8)
    ax.set_ylim(-500, 500)
    ax.set_ylabel("Z components, nT")
    ax.set_xlabel("Time since 12 UT on 10 May 2024, UT")

    ts.save("simulation/May2024/magnetometers.png")
    ts.close()
    return


if __name__ == "__main__":
    # plot_magnetometer_datasets()
    run_May2024_storm()
