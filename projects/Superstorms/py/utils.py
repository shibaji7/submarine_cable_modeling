from types import SimpleNamespace

from scubas.datasets import PROFILES


def create_from_lat_lon(
    dSegments,
    profiles,
    width=1.0,
    flim=[1e-6, 1e0],
    left_active_termination=None,
    right_active_termination=None,
):
    cable_seg = []
    for i in range(len(dSegments) - 1):
        initial = dSegments[i]
        final = dSegments[i + 1]
        cable_seg.append(
            dict(
                initial=dict(lat=initial[0], lon=initial[1]),
                final=dict(lat=final[0], lon=final[1]),
                sec_id=f"Sec-{i}",
                site=profiles[0],
                active_termination=dict(
                    right=None,
                    left=None,
                ),
            )
        )
    if left_active_termination:
        cable_seg[0]["active_termination"] = dict(
            right=None,
            left=left_active_termination,
        )
    if right_active_termination:
        cable_seg[-1]["active_termination"] = dict(
            right=right_active_termination,
            left=None,
        )
    cable = SimpleNamespace(**dict(cable_seg=cable_seg))
    for seg in cable.cable_seg:
        seg["center"] = dict(
            lat=0.5 * (seg["initial"]["lat"] + seg["final"]["lat"]),
            lon=0.5 * (seg["initial"]["lon"] + seg["final"]["lon"]),
        )
        seg["width"], seg["flim"] = width, flim
    return cable


def get_cable_informations(kind="TAT-8", width=1.0, flim=[1e-6, 1e0]):
    if kind == "TAT-8":
        land50 = PROFILES.CS_E
        land50.layers[0].thickness = 50
        cable = SimpleNamespace(
            **dict(
                cable_seg=[
                    dict(
                        initial=dict(lat=39.6, lon=-74.33),
                        final=dict(lat=38.79, lon=-72.62),
                        sec_id="CS-W",
                        site=PROFILES.CS_W,
                        active_termination=dict(
                            right=None,
                            left=PROFILES.LD,
                        ),
                    ),
                    dict(
                        initial=dict(lat=38.79, lon=-72.62),
                        final=dict(lat=37.11, lon=-68.94),
                        sec_id="DO-1",
                        site=PROFILES.DO_1,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=37.11, lon=-68.94),
                        final=dict(lat=39.80, lon=-48.20),
                        sec_id="DO-2",
                        site=PROFILES.DO_2,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=39.80, lon=-48.20),
                        final=dict(lat=40.81, lon=-45.19),
                        sec_id="DO-3",
                        site=PROFILES.DO_3,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=40.81, lon=-45.19),
                        final=dict(lat=43.15, lon=-39.16),
                        sec_id="DO-4",
                        site=PROFILES.DO_4,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=43.15, lon=-39.16),
                        final=dict(lat=44.83, lon=-34.48),
                        sec_id="DO-5",
                        site=PROFILES.DO_5,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=44.83, lon=-34.48),
                        final=dict(lat=46.51, lon=-22.43),
                        sec_id="MAR",
                        site=PROFILES.MAR,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=46.51, lon=-22.43),
                        final=dict(lat=47.85, lon=-9.05),
                        sec_id="DO-6",
                        site=PROFILES.DO_6,
                        active_termination=dict(
                            right=None,
                            left=None,
                        ),
                    ),
                    dict(
                        initial=dict(lat=47.85, lon=-9.05),
                        final=dict(lat=50.79, lon=-4.55),
                        sec_id="CS-E",
                        site=PROFILES.CS_E,
                        active_termination=dict(
                            right=land50,
                            left=None,
                        ),
                    ),
                ]
            )
        )
    for seg in cable.cable_seg:
        seg["center"] = dict(
            lat=0.5 * (seg["initial"]["lat"] + seg["final"]["lat"]),
            lon=0.5 * (seg["initial"]["lon"] + seg["final"]["lon"]),
        )
        seg["width"], seg["flim"] = width, flim
    return cable
