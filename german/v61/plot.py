"""Reproduce all V61 plots, tables and results: python plot.py (from any cwd).

Unless explicitly specified, an uncertainty is one unit in the last
*written* decimal place of a measurement. For polarisation, the power uncertainty is
0.01 mW, matching the last written decimal place. Decimal strings are retained before conversion.
Lengths recorded in cm use a standard uncertainty of 1 mm (0.1 cm);
lengths recorded in mm use 0.5 mm.
For stability power measurements, a fixed uncertainty of 0.1 mW applies.
For longitudinal peak frequencies, no measurement uncertainty is specified
or included; comparisons retain only the resonator-length uncertainty.
Fits use absolute uncertainties, including x uncertainties and common
nuisance parameters. Covariances are never rescaled to force chi2/ndf = 1.
"""

from dataclasses import dataclass
from decimal import Decimal
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
LENGTH_SIGMA_MM = 0.5
LENGTH_SIGMA_CM = 0.1
BUILD = ROOT / "build"
BUILD.mkdir(exist_ok=True)
for variable, directory in (("MPLCONFIGDIR", "matplotlib-cache"),
                            ("TEXMFCACHE", "tex-cache"), ("TEXMFVAR", "tex-cache")):
    cache = BUILD / directory
    cache.mkdir(exist_ok=True)
    os.environ.setdefault(variable, str(cache))

import matplotlib
import numpy as np
from PIL import Image, ImageOps
from pillow_heif import register_heif_opener
from scipy.constants import atomic_mass, c, k
from scipy.optimize import curve_fit
from uncertainties import UFloat, correlated_values, ufloat
from uncertainties import unumpy as unp


# Retain the project's LaTeX typography, independent of the caller's cwd.
matplotlib.use("pgf")
matplotlib.rc_file(ROOT.parent / "matplotlibrc")
matplotlib.rcParams["pgf.preamble"] = (
    r"\input{" + str(ROOT.parent / "header-matplotlib.tex") + "}"
)
import matplotlib.pyplot as plt


@dataclass
class Measurements:
    text: list[list[str]]
    values: np.ndarray
    sigma: np.ndarray
    uncertain: np.ndarray
    header: str


def last_digit(token):
    """E.g. '0.020' -> 0.001; '0' and '160' -> 1, not 10 or 100."""
    return float(Decimal(10) ** Decimal(token).as_tuple().exponent)


def read_data(name, *, sigma_overrides=None):
    lines = (ROOT / "data" / name).read_text(encoding="utf-8").splitlines()
    rows = [line.split("#")[0].split() for line in lines
            if line.split("#")[0].strip()]
    values = np.array(rows, dtype=float)
    sigma = np.array([[last_digit(token) for token in row] for row in rows])
    for column, uncertainty in (sigma_overrides or {}).items():
        sigma[:, column] = uncertainty
    return Measurements(rows, values, sigma, unp.uarray(values, sigma),
                        "\n".join(line for line in lines if line.startswith("#")))


@dataclass
class PeakSpectrum:
    length_text: str
    frequency_text: list[str]
    length: UFloat  # Metres; the source file records centimetres.
    frequencies: tuple[float, ...]  # MHz; no frequency uncertainty is specified.


def read_peak_spectra(name):
    """Read a variable number of peaks per resonator length, without padding."""
    spectra = []
    for line in (ROOT / "data" / name).read_text(encoding="utf-8").splitlines():
        tokens = line.split("#")[0].split()
        if not tokens:
            continue
        if len(tokens) < 2:
            raise ValueError("Each resonator length needs at least one peak.")
        length = ufloat(float(tokens[0]), LENGTH_SIGMA_CM) / 100
        spectra.append(PeakSpectrum(tokens[0], tokens[1:], length,
                                    tuple(float(token) for token in tokens[1:])))
    return spectra


def mean_mode_spacing(frequencies):
    """Arithmetic mean of f_n/n, without an assumed frequency uncertainty."""
    return sum(frequency / n for n, frequency in enumerate(frequencies, 1)) / len(frequencies)


@dataclass
class Fit:
    values: np.ndarray
    covariance: np.ndarray
    parameters: tuple
    chi2: float
    dof: int
    full_values: np.ndarray

    @property
    def reduced_chi2(self):
        return self.chi2 / self.dof


def fit_curve(model, xdata, observed, sigma, initial, parameter_count=None):
    """Weighted curve_fit followed by correlated UFloat parameters.

    ``observed`` includes any measured x coordinates or common distances,
    and ``model`` predicts those measurements along with the measured y.
    Their true values are extra fit parameters. The corresponding block of
    curve_fit's full pcov therefore includes their uncertainty contribution.
    """
    popt, pcov = curve_fit(
        model, xdata, observed, p0=initial, sigma=sigma, absolute_sigma=True,
        method="trf", jac="3-point", x_scale="jac", max_nfev=5000,
        ftol=1e-11, xtol=1e-11, gtol=1e-10,
    )
    if not np.isfinite(pcov).all():
        raise RuntimeError("curve_fit could not determine a finite covariance")
    residual = (model(xdata, *popt) - observed) / sigma
    n = len(popt) if parameter_count is None else parameter_count
    covariance = pcov[:n, :n]
    return Fit(popt[:n], covariance, tuple(correlated_values(popt[:n], covariance)),
               float(residual @ residual), len(observed) - len(popt), popt)


def errors_in_variables(model, data, initial):
    x, y = data.values.T
    sx, sy = data.sigma.T
    n = len(initial)

    def joint_model(_index, *parameters):
        true_x = parameters[n:]
        return np.concatenate((true_x, model(np.asarray(true_x), *parameters[:n])))

    observed = np.r_[x, y]
    sigma = np.r_[sx, sy]
    index = np.arange(len(observed))
    fit = fit_curve(joint_model, index, observed, sigma, np.r_[initial, x], n)
    if model is gaussian:
        # Close to the maximum an uncertain x can lie on either flank.
        # Check reflected latent coordinates to avoid a local branch minimum.
        for _ in range(2):
            for i in np.flatnonzero(np.abs(x - fit.values[1]) < 2 * sx):
                start = fit.full_values.copy()
                start[n + i] = 2 * fit.values[1] - start[n + i]
                alternative = fit_curve(joint_model, index, observed, sigma, start, n)
                if alternative.chi2 < fit.chi2:
                    fit = alternative
    return fit


def gaussian(x, amplitude, centre, width, background):
    return background + amplitude * np.exp(-2 * ((x - centre) / width) ** 2)


def hermite(x, amplitude, centre, width, background):
    t = (x - centre) / width
    # amplitude is the height of either lobe above the background.
    return background + amplitude * np.e * 2 * t**2 * np.exp(-2 * t**2)


def malus(angle, amplitude, theta, background):
    return background + amplitude * np.cos(np.deg2rad(angle - theta)) ** 2


def grating_position(order, wavelength_nm, centre, distance, density):
    """Screen position in mm from the grating equation; density in lines/mm."""
    sine = np.asarray(order) * wavelength_nm * density * 1e-6
    return centre + distance * sine / np.sqrt(1 - sine**2)


def fit_wavelength(data, order, density, distance, distance_sigma):
    """Fit one grating series, including its measured screen distance once."""
    observed = np.r_[data.values[:, 0], distance]
    sigma = np.r_[data.sigma[:, 0], distance_sigma]

    def model(_index, wavelength_nm, centre, true_distance):
        return np.r_[grating_position(order, wavelength_nm, centre, true_distance, density),
                     true_distance]

    initial = [650, data.values[np.flatnonzero(order == 0)[0], 0], distance]
    return fit_curve(model, np.arange(len(observed)), observed, sigma, initial)


RESULTS = {}
MACROS = {}


def quantity(name, value, unit=""):
    text = format(value, ".2u").replace("+/-", "+-")
    MACROS[name] = rf"\qty{{{text}}}{{{unit}}}" if unit else rf"\num{{{text}}}"
    RESULTS[name] = {"value": float(value.n), "sigma": float(value.s), "unit": unit}


def nominal_quantity(name, value, unit):
    """Report a value whose measurement uncertainty has not been specified."""
    MACROS[name] = rf"\qty{{{value:.2f}}}{{{unit}}}"
    RESULTS[name] = {"value": float(value), "sigma": None, "unit": unit}


def number(name, value, decimals=2):
    MACROS[name] = rf"\num{{{value:.{decimals}f}}}"
    RESULTS[name] = float(value)


def fit_statistics(prefix, fit):
    number(prefix + "Chi", fit.reduced_chi2)
    number(prefix + "Dof", fit.dof, 0)
    RESULTS[prefix + "Covariance"] = fit.covariance.tolist()
    RESULTS[prefix + "ChiSquare"] = fit.chi2


def save(fig, filename):
    fig.savefig(BUILD / filename)
    plt.close(fig)


def measured(ax, data, **kwargs):
    ax.errorbar(*data.values.T, xerr=data.sigma[:, 0], yerr=data.sigma[:, 1],
                fmt="o", markersize=3, capsize=2, linewidth=0.8, **kwargs)


def fit_plot(name, data, model, fit, kind):
    fig, ax = plt.subplots(figsize=(5.78, 3.6))
    measured(ax, data, color="black", label="Messwerte")
    x = np.linspace(data.values[:, 0].min() - 1, data.values[:, 0].max() + 1, 600)
    ax.plot(x, model(x, *fit.values), color="tab:blue", label="Fit")
    if kind == "polarisation":
        ax.set(xlabel=r"$\varphi \mathbin{/} \unit{\degree}$",
               ylabel=r"$P \mathbin{/} \unit{\milli\watt}$")
    else:
        ax.set(xlabel=r"$x \mathbin{/} \unit{\milli\meter}$",
               ylabel=r"$P \mathbin{/} \unit{\milli\watt}$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    save(fig, name)


def raw(token):
    return rf"\num{{{token}}}"


def raw_uncertainty(token):
    sigma = str(Decimal(10) ** Decimal(token).as_tuple().exponent)
    return rf"\num{{{token}+-{sigma}}}"


def tex_uncertain(value):
    return r"\num{" + format(value, ".2u").replace("+/-", "+-") + "}"


def table(name, caption, label, headings, rows, *, upper_headings=None, alignment="r", colspec=None):
    header_rows = [upper_headings, headings] if upper_headings else [headings]
    header_style = "row{1,2}" if upper_headings else "row{1}"
    text = [r"\begin{table}[H]", r"\centering", r"\small", rf"\caption{{{caption}}}",
            rf"\label{{{label}}}",
            r"\begin{tblr}{colspec={" + (colspec or " ".join(alignment for _ in headings))
            + r"}, rowsep=1pt, " + header_style + r"={guard, mode=math}}", r"\toprule"]
    text += [" & ".join(row) + r" \\" for row in header_rows]
    text += [r"\midrule"]
    text += [" & ".join(row) + r" \\" for row in rows]
    text += [r"\bottomrule", r"\end{tblr}", r"\end{table}", ""]
    (BUILD / name).write_text("\n".join(text), encoding="utf-8")


def two_halves(name, data, caption, label, x_heading, y_heading, *, show_uncertainties=True):
    half = (len(data.text) + 1) // 2
    rows = []
    for i in range(half):
        row = []
        for j in (i, i + half):
            format_value = raw_uncertainty if show_uncertainties else raw
            row += ([raw(data.text[j][0]), format_value(data.text[j][1])]
                    if j < len(data.text) else ["", ""])
        rows.append(row)
    table(name, caption, label, [x_heading, y_heading] * 2, rows)


def main():
    # Convert the original photographs for LaTeX without cropping or retouching.
    register_heif_opener()
    for mode in ("00", "10"):
        with Image.open(ROOT / "data" / f"{mode}.HEIC") as photo:
            ImageOps.exif_transpose(photo).convert("RGB").save(BUILD / f"tem{mode}_foto.png")

    zero = read_data("00_mode.txt", sigma_overrides={0: LENGTH_SIGMA_MM})
    first = read_data("10_mode.txt", sigma_overrides={0: LENGTH_SIGMA_MM})
    polarisation = read_data("polarisation.txt")
    stability = [read_data(f"stabilitätsbedingung_{i}.txt",
                           sigma_overrides={0: LENGTH_SIGMA_CM, 1: 0.1}) for i in (1, 2)]
    longitudinal = read_peak_spectra("longitudinalmoden.txt")
    wavelength = [read_data(f"wellenlänge_{i}.txt", sigma_overrides={0: LENGTH_SIGMA_MM}) for i in (1, 2)]

    for prefix, data, model, initial, kind, filename in (
        ("vZero", zero, gaussian, [1, -7, 10, 0], "gaussian", "tem00.pdf"),
        ("vFirst", first, hermite, [.17, -10, 10, 0], "hermite", "tem10.pdf"),
        ("vPol", polarisation, malus, [4.17, 80, 0], "polarisation", "polarisation.pdf"),
    ):
        fit = errors_in_variables(model, data, initial)
        fit_statistics(prefix, fit)
        names = ("A", "Theta", "B") if kind == "polarisation" else ("A", "Centre", "Width", "B")
        units = ([r"\milli\watt", r"\degree", r"\milli\watt"]
                 if kind == "polarisation" else
                 [r"\milli\watt", r"\milli\meter", r"\milli\meter", r"\milli\watt"])
        for name, parameter, unit in zip(names, fit.parameters, units):
            quantity(prefix + name, parameter, unit)
        if kind == "hermite":
            _, centre, width, _ = fit.parameters
            quantity("vFirstLeft", centre - width / np.sqrt(2), r"\milli\meter")
            quantity("vFirstRight", centre + width / np.sqrt(2), r"\milli\meter")
            quantity("vFirstSeparation", np.sqrt(2) * width, r"\milli\meter")
        if kind == "polarisation":
            quantity("vPolMinimum", fit.parameters[1] + 90, r"\degree")
            quantity("vPolMaximum", fit.parameters[0] + fit.parameters[2], r"\milli\watt")
        fit_plot(filename, data, model, fit, kind)

    for data, name, label, caption in (
        (zero, "table_tem00.tex", "tab:v61-tem00", r"Messwerte des TEM$_{00}$-Profils: gemessene Leistung $P$ als Funktion der transversalen Photodiodenposition $x$. Für die Positionsunsicherheit gilt $\Delta x=\qty{0.5}{\milli\meter}$. Die Leistungsunsicherheiten folgen der letzten notierten Stelle."),
        (first, "table_tem10.tex", "tab:v61-tem10", r"Messwerte des TEM$_{10}$-Profils: gemessene Leistung $P$ als Funktion der transversalen Photodiodenposition $x$. Die Positionsunsicherheit beträgt $\Delta x=\qty{0.5}{\milli\meter}$. Die Leistungsunsicherheiten folgen der letzten notierten Stelle."),
    ):
        two_halves(name, data, caption, label,
                   r"x/\unit{\milli\meter}", r"P/\unit{\milli\watt}")
    two_halves("table_polarisation.tex", polarisation,
               r"Polarisationsmessung: gemessene Leistung $P$ in Abhängigkeit vom Winkel $\varphi$ des Polarisationsfilters. Die Winkelunsicherheit beträgt $\Delta\varphi=\qty{1}{\degree}$, die Leistungsunsicherheit $\Delta P=\qty{0.01}{\milli\watt}$.",
               "tab:v61-polarisation", r"\varphi/\unit{\degree}", r"P/\unit{\milli\watt}",
               show_uncertainties=False)

    rows = []
    for i in range(max(len(d.text) for d in stability)):
        row = []
        for data in stability:
            if i < len(data.text):
                L, P = data.text[i]
                row += [raw(L), raw(P)]
            else:
                row += ["", ""]
        rows.append(row)
    table("table_stability.tex",
          r"Messdaten zur Resonatorstabilität: gemessene Leistung $P$ bei der Resonatorlänge $L$. Die Indizes $\mathrm{kk}$ und $\mathrm{pk}$ bezeichnen die konkav-konkave bzw. plan-konkave Spiegelkonfiguration. Für die in Zentimetern notierten Längenmessungen gilt $\Delta L=\qty{1}{\milli\meter}$ und für alle Leistungsmessungen gilt $\Delta P=\qty{0.1}{\milli\watt}$.",
          "tab:v61-stabilitaet", [r"L_{\mathrm{kk}}/\unit{\centi\meter}",
          r"P_{\mathrm{kk}}/\unit{\milli\watt}", r"L_{\mathrm{pk}}/\unit{\centi\meter}",
          r"P_{\mathrm{pk}}/\unit{\milli\watt}"], rows)

    fig, (ax, ag) = plt.subplots(2, 1, sharex=True, figsize=(5.78, 4.0),
                                gridspec_kw={"height_ratios": [1, 1]})
    L_grid = np.linspace(0, 300, 1000)
    colors = ["tab:blue", "tab:orange"]
    for i, (data, color) in enumerate(zip(stability, colors)):
        label = "konkav--konkav" if i == 0 else "plan--konkav"
        measured(ax, data, color=color, label=label)
        g = (1 - L_grid / 140) ** (2 if i == 0 else 1)
        ag.plot(L_grid, g, color=color, label=label)
    ag.axhspan(0, 1, color="tab:green", alpha=.10)
    for axis in (ax, ag):
        axis.axvline(140, color=colors[1], linestyle=":", linewidth=1)
        axis.axvline(280, color=colors[0], linestyle=":", linewidth=1)
        axis.grid(alpha=.25)
    ax.set(ylabel=r"$P \mathbin{/} \unit{\milli\watt}$")
    ag.set(xlabel=r"$L \mathbin{/} \unit{\centi\meter}$", ylabel=r"$g_1 g_2$",
           ylim=(-1.2, 1.4), xlim=(0, 300))
    ax.legend(fontsize=9)
    ag.legend(fontsize=9)
    save(fig, "stabilitaet.pdf")
    quantity("vStableLast", stability[0].uncertain[-1, 0], r"\centi\meter")
    quantity("vPlanLast", stability[1].uncertain[2, 0], r"\centi\meter")

    # Each f_n/n estimates the same spacing. Average within each spectrum;
    # the measured length enters the theoretical comparison once, with its error.
    L = np.array([spectrum.length.n for spectrum in longitudinal])
    sL = np.array([spectrum.length.s for spectrum in longitudinal])
    spacings = np.array([mean_mode_spacing(spectrum.frequencies)
                         for spectrum in longitudinal])
    expected = np.array([c / (2 * spectrum.length) / 1e6
                         for spectrum in longitudinal], dtype=object)
    differences = spacings - expected
    relative = (spacings / expected - 1) * 100
    nominal_quantity("vLongSpacingMin", spacings.min(), r"\mega\hertz")
    nominal_quantity("vLongSpacingMax", spacings.max(), r"\mega\hertz")
    number("vLongRelativeMin", min(value.n for value in relative))
    number("vLongRelativeMax", max(value.n for value in relative))

    peak_count = max(len(spectrum.frequencies) for spectrum in longitudinal)
    rows = [[raw(spectrum.length_text), *map(raw, spectrum.frequency_text),
             *([r"\textemdash"] * (peak_count - len(spectrum.frequencies)))]
            for spectrum in longitudinal]
    table("table_longitudinal.tex",
          r"Resonatorlängen $L$ und im FFT-Spektrum abgelesene Frequenzpositionen $f_n$ der Schwebungspeaks. $n$ nummeriert die Peaks nach aufsteigender Frequenz. Für die Längenunsicherheit gilt $\Delta L=\qty{0.1}{\centi\meter}$. Ein Strich kennzeichnet eine fehlende Peakangabe.",
          "tab:v61-longitudinal", [r"L/\unit{\centi\meter}",
          *[rf"f_{n}/\unit{{\mega\hertz}}" for n in range(1, peak_count + 1)]], rows)
    rows = [[raw(spectrum.length_text), raw(str(len(spectrum.frequencies))),
             raw(f"{spacing:.2f}"), tex_uncertain(theory), tex_uncertain(difference)]
            for spectrum, spacing, theory, difference
            in zip(longitudinal, spacings, expected, differences)]
    table("table_mode_spacings.tex",
          r"Modenabstände bei der Resonatorlänge $L$: $\overline{\Delta\nu}$ ist der aus den Messwerten ermittelte mittlere Modenabstand, $\Delta\nu_{\mathrm{theo}}=c/(2L)$ der theoretische Modenabstand mit der Lichtgeschwindigkeit $c$. $N$ ist die Anzahl der berücksichtigten Peaks und $d_\nu=\overline{\Delta\nu}-\Delta\nu_{\mathrm{theo}}$ die Abweichung vom theoretischen Wert.",
          "tab:v61-mode-spacings", [r"L/\unit{\centi\meter}", r"N",
          r"\overline{\Delta\nu}/\unit{\mega\hertz}",
          r"\Delta\nu_{\mathrm{theo}}/\unit{\mega\hertz}",
          r"d_\nu/\unit{\mega\hertz}"], rows)

    def value_and_sigma(value):
        if isinstance(value, UFloat):
            return {"value": float(value.n), "sigma": float(value.s)}
        return {"value": float(value), "sigma": None}

    RESULTS["vLongitudinal"] = [
        {"length_cm": value_and_sigma(spectrum.length * 100),
         "peak_count": len(spectrum.frequencies),
         "frequencies_MHz": [value_and_sigma(frequency) for frequency in spectrum.frequencies],
         "normalized_frequencies_MHz": [value_and_sigma(frequency / n)
                                        for n, frequency in enumerate(spectrum.frequencies, 1)],
         "spacing_MHz": value_and_sigma(spacing), "theory_MHz": value_and_sigma(theory),
         "difference_MHz": value_and_sigma(difference),
         "relative_difference_percent": value_and_sigma(deviation),
         "uncertainty_sources": ["resonator_length"]}
        for spectrum, spacing, theory, difference, deviation
        in zip(longitudinal, spacings, expected, differences, relative)
    ]

    fig, (ax, ar) = plt.subplots(2, 1, sharex=True, figsize=(5.78, 5.0),
                                gridspec_kw={"height_ratios": [2, 1]})
    grid = np.linspace(L.min() * .93, L.max() * 1.05, 500)
    ax.errorbar(L, spacings, xerr=sL,
                fmt="o", markersize=3, color="tab:blue", capsize=2,
                label="Messwerte")
    ax.plot(grid, c / (2e6 * grid), color="tab:orange", label=r"Theorie $c/(2L)$")
    ar.errorbar(L, unp.nominal_values(differences), xerr=sL,
                yerr=unp.std_devs(differences), fmt="o", capsize=2,
                color="tab:blue", markersize=3)
    ar.axhline(0, color="black", linewidth=.8)
    ax.set(ylabel=r"$\Delta\nu \mathbin{/} \unit{\mega\hertz}$")
    ar.set(xlabel=r"$L \mathbin{/} \unit{\meter}$",
           ylabel=r"$d_\nu \mathbin{/} \unit{\mega\hertz}$")
    ax.legend(fontsize=9)
    for axis in (ax, ar):
        axis.grid(alpha=.25)
    save(fig, "longitudinalmoden.pdf")
    doppler = np.sqrt(8 * k * 293.15 * np.log(2) / (20.1797 * atomic_mass)) / 632.8e-9
    number("vDoppler", doppler / 1e9, 2)
    number("vDopplerRatioMin", doppler / (spacings.max() * 1e6), 1)
    number("vDopplerRatioMax", doppler / (spacings.min() * 1e6), 1)

    # Grating densities are apparatus nominal values, not measured quantities.
    densities = [float(re.search(r"Gitter:\s*([\d.]+)", d.header)[1]) for d in wavelength]
    distances = []
    for data in wavelength:
        if "Positionen [mm]" not in data.header:
            raise ValueError("Check wavelength position units; mm must be confirmed.")
        token = re.search(r"Schirm\s+([\d.,]+)cm", data.header)[1].replace(",", ".")
        distances.append((float(token) * 10, LENGTH_SIGMA_CM * 10))
    orders = []
    for data in wavelength:
        zero = np.flatnonzero(np.isclose(data.values[:, 0], 6.0))
        if len(zero) != 1:
            raise ValueError("Each grating series must contain one zero-order maximum at 6 mm.")
        orders.append(zero[0] - np.arange(len(data.values)))
    wave_fits = []
    for name, data, order, density, (distance, distance_sigma) in zip(
            ("One", "Two"), wavelength, orders, densities, distances):
        fit = fit_wavelength(data, order, density, distance, distance_sigma)
        wave_fits.append(fit)
        prefix = f"vWave{name}"
        quantity(prefix, fit.parameters[0], r"\nano\meter")
        quantity(prefix + "Centre", fit.parameters[1], r"\milli\meter")
        quantity(prefix + "Distance", fit.parameters[2], r"\milli\meter")
        quantity(prefix + "Deviation", (fit.parameters[0] / 632.8 - 1) * 100, r"\percent")
        fit_statistics(prefix, fit)
        number(prefix + "ReferenceSigma", abs(fit.parameters[0].n - 632.8) / fit.parameters[0].s)
        quantity(f"vScreenDistance{name}", ufloat(distance, distance_sigma) / 10,
                 r"\centi\meter")
    positions_by_order = [dict(zip(order, data.text))
                          for order, data in zip(orders, wavelength)]
    rows = [[raw(str(int(m))),
             *[positions[m][0] + "+-" + str(LENGTH_SIGMA_MM)
               if m in positions else r"\SetCell{guard} ---"
               for positions in positions_by_order]]
            for m in sorted(set(np.concatenate(orders)), reverse=True)]
    distance_caption = (r"Der Abstand $d$ zwischen Gitter und Schirm beträgt jeweils \mbox{$d=\vScreenDistanceOne$}."
                        if distances[0] == distances[1] else
                        r"Die Schirmabstände betragen \mbox{$d_1=\vScreenDistanceOne$} und \mbox{$d_2=\vScreenDistanceTwo$}.")
    table("table_wavelength.tex",
          r"Gemessene Positionen $y_m$ der Gittermaxima der Beugungsordnung $m$ für die beiden Gitter mit der Liniendichte $N$. " + distance_caption,
          "tab:v61-wellenlaenge",
          [r"m", r"y_m/\unit{\milli\meter}", r"y_m/\unit{\milli\meter}"], rows,
          upper_headings=[r"\text{Beugungsordnung}",
                          r"N=\qty[per-mode=reciprocal]{80}{\per\milli\meter}",
                          r"N=\qty[per-mode=reciprocal]{100}{\per\milli\meter}"],
          alignment="c",
          colspec=r"c S[table-format=-2.1(1),table-align-uncertainty=true] S[table-format=-2.1(1),table-align-uncertainty=true]")
    fig, axes = plt.subplots(1, 2, figsize=(5.78, 3.4), sharey=True)
    for ax, data, order, density, fit, color in zip(
            axes, wavelength, orders, densities, wave_fits, colors):
        grid = np.linspace(order.min() - .05, order.max() + .05, 300)
        ax.errorbar(order, data.values[:, 0], yerr=data.sigma[:, 0],
                    fmt="o", color="black", markersize=3, capsize=2, label="Messwerte")
        ax.plot(grid, grating_position(grid, *fit.values, density),
                color=color, label="Fit")
        ax.set(xlabel=r"Beugungsordnung $m$", title=rf"$N={density:g}\,\unit[per-mode=reciprocal]{{\per\milli\meter}}$",
               xticks=sorted(order))
        ax.legend(fontsize=7, loc="upper left")
        ax.grid(alpha=.25)
    axes[0].set_ylabel(r"$y \mathbin{/} \unit{\milli\meter}$")
    save(fig, "wellenlaenge.pdf")

    macro_text = "% Generated by plot.py; do not edit manually.\n"
    macro_text += "\n".join(rf"\newcommand{{\{name}}}{{{value}}}" for name, value in MACROS.items()) + "\n"
    (BUILD / "results.tex").write_text(macro_text, encoding="utf-8")
    (BUILD / "results.json").write_text(json.dumps(RESULTS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({name: value for name, value in RESULTS.items() if "Covariance" not in name},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
