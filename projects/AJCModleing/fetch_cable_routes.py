import datetime as dt
from types import SimpleNamespace

import numpy as np
import pandas as pd
import requests
from geopy.distance import geodesic


def calculate_bathymetry_byLITHO1(o, distance_interval=600):
    # 1. Compute 10 points between 2 points
    # 2. Compute water depth for each locations
    # repete steps 1/2 for all points
    geolats, geolongs = [], []
    for xy in o.geometry["coordinates"]:
        lons, lats = np.array(xy)[:, 0], np.array(xy)[:, 1]
        total_distance = geodesic((lats[0], lons[0]), (lats[-1], lons[-1])).km
        if total_distance > distance_interval:
            for i in range(len(lats) - 1):
                td_km = geodesic((lats[i], lons[i]), (lats[i + 1], lons[i + 1])).km
            geolats.extend(lats)
            geolongs.extend(lons)
    d = pd.DataFrame()
    d["geolats"], d["geolongs"] = geolats, geolongs
    d = d.sort_values(by="geolats")
    return d


def calculate_conductive_profiles_with_distance(dp, dpn, base_name="AJC"):
    from scubas.conductivity import ConductivityProfile
    from scubas.datasets import Site

    cp = ConductivityProfile()
    profiles = []
    bin_n = (dpn.geolats, dpn.geolongs)
    for i in range(len(dp) - 1):
        bin_i, bin_j = (
            (dp.geolats.iloc[i], dp.geolongs.iloc[i]),
            (dp.geolats.iloc[i + 1], dp.geolongs.iloc[i + 1]),
        )
        ipts = cp.get_interpolation_points(bin_i, bin_j)
        profile = cp._compile_profile_(ipts)
        profile = Site.init(
            1.0 / profile["resistivity"].to_numpy(dtype=float),
            profile["thickness"].to_numpy(dtype=float) * 1e3,  # Convert to m
            profile["name"],
            "",
            base_name + f"_{i}",
        )
        td_km = geodesic(bin_i, bin_n).km
        profiles.append(
            dict(
                profile=profile,
                bin_i=bin_i,
                bin_j=bin_n,
                td_km=td_km,
                depth=profile.get_thicknesses(0),
            )
        )
    return profiles


def calculate_conductive_profiles(d, base_name="AJC"):
    from scubas.conductivity import ConductivityProfile
    from scubas.datasets import Site

    cp = ConductivityProfile()
    profiles = []
    for i in range(len(d) - 1):
        bin_i, bin_j = (
            (d.geolats.iloc[i], d.geolongs.iloc[i]),
            (d.geolats.iloc[i + 1], d.geolongs.iloc[i + 1]),
        )
        ipts = cp.get_interpolation_points(bin_i, bin_j)
        profile = cp._compile_profile_(ipts)
        profile = Site.init(
            1.0 / profile["resistivity"].to_numpy(dtype=float),
            profile["thickness"].to_numpy(dtype=float) * 1e3,  # Convert to m
            profile["name"],
            "",
            base_name + f"_{i}",
        )
        td_km = geodesic(bin_i, bin_j).km
        profiles.append(
            dict(
                profile=profile,
                bin_i=bin_i,
                bin_j=bin_j,
                td_km=td_km,
                depth=profile.get_thicknesses(0),
            )
        )
    return profiles


def plot_routes(o, geo, fname="figures/ajc_routes.png", d=dt.datetime(1958, 2, 11)):
    from fan import CartoDataOverlay

    cb = CartoDataOverlay(
        date=d,
        central_longitude=130,
        central_latitude=20,
        extent=[110, 170, -50, 50],
        plt_lats=np.arange(-90, 80, 10),
    )
    ax = cb.add_axes()
    for xy in o.geometry["coordinates"]:
        xy = np.array(xy)
        lon, lat = xy[:, 0], xy[:, 1]
        xyz = cb.proj.transform_points(cb.geo, lon, lat)
        ax.plot(xyz[:, 0], xyz[:, 1], ls="-", lw=0.8, color="r", transform=cb.proj)
    xyz = cb.proj.transform_points(
        cb.geo, np.array(geo.geolongs), np.array(geo.geolats)
    )
    ax.plot(xyz[:, 0], xyz[:, 1], ".", ms=0.8, color="b", transform=cb.proj)
    o = pd.read_csv("datasets/20250211-20-39-supermag.csv", parse_dates=["Date_UTC"])
    iagas = o.IAGA.unique()
    for iaga in iagas:
        x = o[o.IAGA == iaga]
        # print([x["GEOLON"].tolist()[0]], [x["GEOLAT"].tolist()[0]])
        Lon, Lat = [x["GEOLON"].tolist()[0]], [x["GEOLAT"].tolist()[0]]
        xyz = cb.proj.transform_points(cb.geo, np.array(Lon), np.array(Lat))
        ax.scatter(xyz[:, 0], xyz[:, 1], s=4, color="m", marker="D", transform=cb.proj)
    cb.save(fname)
    cb.close()
    return


def get_cable_route(
    url="https://www.submarinecablemap.com/api/v3/cable/cable-geo.json",
    name_key="australia-japan-cable-ajc",
):
    o = None
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        for d in data["features"]:
            if d["properties"]["id"] == name_key:
                o = SimpleNamespace(**d)
    return o


def plot_bathymatry(profiles):
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "Tahoma",
        "DejaVu Sans",
        "Lucida Grande",
        "Verdana",
    ]

    fig = plt.figure(figsize=(5, 2), dpi=300)
    ax = fig.add_subplot(111)
    distance, depths = [], []
    for profile in profiles:
        distance.append(profile["td_km"])
        depths.append(profile["depth"])
    ax.plot(
        np.cumsum(distance),
        np.array(depths) / 1e3,
        ls="-",
        lw=0.8,
        color="r",
    )
    ax.invert_yaxis()
    ax.set_ylim(8, 0)
    ax.set_xlabel("Distance, km")
    ax.set_ylabel("Depths, km")
    ax.set_xlim(0, np.cumsum(distance)[-1])
    fig.savefig("figures/ajc_route_bathymetry.png", bbox_inches="tight")
    print(np.cumsum(distance), depths)
    return


if __name__ == "__main__":
    o = get_cable_route()
    d = calculate_bathymetry_byLITHO1(o)
    profiles = calculate_conductive_profiles(d)
    plot_routes(o, d)
    plot_bathymatry(profiles)
