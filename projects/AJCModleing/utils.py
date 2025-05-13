from types import SimpleNamespace


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
