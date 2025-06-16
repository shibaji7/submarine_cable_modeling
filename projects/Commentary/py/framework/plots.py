from typing import Optional, Sequence

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import scienceplots
plt.style.use(["science", "ieee"])
plt.rcParams.update(
    {
        "text.usetex": False,
    }
)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import matplotlib as mpl

mpl.rc("font", size=15)


class StackPlots:
    def __init__(
        self,
        nrows: int,
        ncols: int,
        dpi: int = 1000,
        figsize=(4, 3),
        sharex=True,
        sharey=True,
    ):
        self.nrows = nrows
        self.ncols = ncols
        self.sharex = (sharex,)
        self.sharey = (sharey,)
        self.fig = plt.figure(
            figsize=(self.ncols * figsize[0], self.nrows * figsize[1]),
            dpi=dpi,
        )
        self.axes = []
        self.fig.subplots_adjust(hspace=0.5, wspace=0.3)

    def get_axes(self, grid, polar=False):
        ax = plt.subplot2grid(
            (self.nrows, self.ncols),
            grid["start"],
            colspan=grid["colspan"],
            rowspan=grid["rowspan"],
            polar=polar,
        )
        self.axes.append(ax)
        return ax

    def plot_stack_plots(
        self,
        time: Sequence,
        value: Sequence,
        grid: dict,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        label: Optional[str] = None,
        text: Optional[str] = None,
        tag: Optional[str] = None,
        ylim: Optional[Sequence] = None,
        xlim: Optional[Sequence] = None,
        color: str = "blue",
        lw: float = 0.8,
        ls: str = "-",
        ax: Optional[plt.Axes] = None,
        ylabel_color: Optional[str] = "k",
        datetime: Optional[bool] = False,
    ) -> tuple:
        """
        Plot a stack of plots with the given time and value data.
        :param time: Time data for the x-axis
        :param value: Value data for the y-axis
        :param title: Title for the plot
        :param xlabel: Label for the x-axis
        :param ylabel: Label for the y-axis
        """
        if ax is None:
            ax = self.get_axes(grid)
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel, color=ylabel_color)
        if ylim:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim([min(value), max(value)])
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim([time[0], time[-1]])
        if datetime:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.plot(time, value, color=color, linewidth=lw, ls=ls, label=label)
        if text:
            ax.text(0.05, 1.05, text, ha="left", va="center", transform=ax.transAxes)
        if tag:
            ax.text(0.05, 0.95, tag, ha="left", va="center", transform=ax.transAxes)

        plt.tight_layout()
        return self.fig, ax

    def save_fig(self, filename: str):
        """
        Save the figure to a file.
        :param filename: Filename to save the figure
        """
        self.fig.savefig(filename, bbox_inches="tight")
        return

    def close(self):
        """
        Close the figure.
        """
        plt.close(self.fig)
        return

    def plot_dirctional_plots(
        self,
        theta: Sequence,
        r: Sequence,
        grid: dict,
        title: Optional[str] = None,
        text: Optional[str] = None,
        tag: Optional[str] = None,
        color: str = "black",
        ax: Optional[plt.Axes] = None,
        rlims: Optional[Sequence] = [0, 1],
        rticks: Optional[Sequence] = [0, 0.5, 1.0],
        theta_ticks: Optional[Sequence] = [0, np.pi / 2, np.pi, 3 * np.pi / 2],
        cable_angle: Optional[float] = None,
        text_loc: dict = dict(x=-0.1, y=1.05, ha="left", va="center",),
    ):
        """
        Plot directional plots.
        """
        if ax is None:
            ax = self.get_axes(grid, polar=True)
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            if cable_angle is not None:
                ax.plot(
                    np.deg2rad([cable_angle, cable_angle + 180]),
                    [1, 1],
                    lw=1.2,
                    ls="-",
                    color="m",
                )
                ax.plot(
                    np.deg2rad([cable_angle + 90, cable_angle + 270]),
                    [1, 1],
                    lw=0.6,
                    ls="--",
                    color="k",
                )
        if title:
            ax.set_title(title, fontdict=dict(size=12))
        ax.plot(np.deg2rad(theta), r, color=color, lw=0.9, ls="-")
        ax.set_rticks(rticks)
        ax.set_xticks(theta_ticks)
        ax.set_rmax(1)
        ax.set_rmin(rlims[0])
        if text:
            ax.text(
                text_loc["x"],
                text_loc["y"],
                text,
                ha=text_loc["ha"],
                va=text_loc["va"],
                transform=ax.transAxes,
                color=color,
            )
        if tag:
            ax.text(0.05, 0.95, tag, ha="left", va="center", transform=ax.transAxes)

        return self.fig, ax
