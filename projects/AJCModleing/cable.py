import datetime as dt
from loguru import logger
import matplotlib.dates as mdates

from scubas.cables import TransmissionLine, Cable

# from plots import TimeSeriesPlot
from utils import get_cable_informations
from fetch_data import clean_B_fields

class SCUBASModel(object):
    
    def __init__(
            self, 
            cable_name="TAT-1",
            cable_structure=get_cable_informations(),
            segment_files=[],
        ):
        self.cable_name = cable_name
        self.cable_structure = cable_structure
        self.segment_files = segment_files
        logger.info(f"Initialize {cable_name}")
        return
    
    def read_stations(
            self, 
            stns=["FRD", "STJ", "HAD"], 
            stn_files=[
                ["dataset/May2024/frd20240510psec.sec.txt"],
                ["dataset/May2024/had20240510psec.sec.txt"], 
                ["dataset/May2024/stj20240510psec.sec.txt"]
            ]
        ):
        self.frames = clean_B_fields(stns, stn_files)
        return
    
    def initialize_TL(self):
        self.tlines = []
        for i, seg in enumerate(self.cable_structure.cable_seg):
            self.tlines.append(
                TransmissionLine(
                    sec_id=seg["sec_id"],
                    directed_length=dict(
                        edge_locations=dict(
                            initial=seg["initial"], 
                            final=seg["final"]
                        )
                    ),
                    elec_params=dict(
                        site=seg["site"],
                        width=seg["width"],
                        flim=seg["flim"],
                    ),
                    active_termination=seg["active_termination"],
                ).compile_oml(self.segment_files[i]),
            )
        return
    
    def run_cable_segment(self):
        # Running the cable operation
        logger.info(f"Components: {self.tlines[0].components}")
        self.cable = Cable(self.tlines, self.tlines[0].components)
        return
    
    def plot_TS_with_others(
        self, date_lim, fname, fig_title, vlines,
        themis_fgm, themis_mom, stations,
        major_locator=mdates.MinuteLocator(byminute=range(0, 60, 5)),
        minor_locator=mdates.MinuteLocator(byminute=range(0, 60, 1)), 
        text_size=15
    ):
        ts = TimeSeriesPlot(
            date_lim,
            major_locator=major_locator,
            minor_locator=minor_locator,
            fig_title=fig_title, 
            text_size=text_size,
            num_subplots=3,
        )
        if len(themis_fgm) and len(themis_mom):
            ts.add_vlines(
                ts.add_themis(themis_fgm, themis_mom, ["thc_fgs_gsm", "pdyn"]),
                vlines=vlines, colors=["r"]
            )
        ts.add_vlines(
            ts.add_mag(self.frames, stations, ylim=[-200, 500]), 
            vlines=vlines, colors=["r"]
        )
        ts.add_vlines(
            ts.add_voltage(self.cable.tot_params, xlabel="Minutes since 17 UT", ylim=[-200, 100]),
            vlines=vlines, colors=["r"]
        )
        ts.save(fname)
        ts.close()
        return
    
def get_cable_segments():
    return