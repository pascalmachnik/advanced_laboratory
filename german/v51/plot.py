import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from measurement import readData, Measurement
from scipy.signal import find_peaks
import pandas as pd
import uncertainties as unc
from uncertainties import unumpy as unp

# example usage of measurement.py
# m = readData(
#    filepath="data/inv_amp_1",
#    unit_freq="kHz",           #options are hz, khz, mhz
#    unit_volt="mV",            #options are mV, V
#    experiment="inv_op_amp",   #options are inv_op_amp or int_diff
#    Vin=5,   #in V             #when chosing inv_op_amp set R1 and R2    for int_diff set R and C like below
#    R1=1000, #in Ohm           #R = 1000
#    R2=10    #in Ohm           #C = 5e-6
# )
#
# print(m)       # returns  "Measurement with 40 points, Vin=5V, Experiment=inv_op_amp"
# Measurement class contains  frequencies in hz, voltages in V, Vin in V, experiment as a string


def makePlot1(measurement):
    frequencies = measurement.frequencies
    voltages_out = measurement.voltages
    Vin = measurement.Vin
    R1 = measurement.R1
    R2 = measurement.R2

    # Compute Gain (linear and dB)
    gain = np.abs(voltages_out / Vin)
    gain = np.clip(gain, 1e-12, None)  # prevent log10(0)
    gain_db = 20 * np.log10(gain)

    # Bode Plot
    fig, ax = plt.subplots(layout="constrained")
    ax.grid(which='both', linestyle='--')
    ax.set_ylabel(r"Verstärkung $\mathbin{/} \unit{\dB}$")
    ax.set_xlabel(r"Frequenz $\mathbin{/} \unit{\hertz}$")
    ax.set_title('Bode-Diagramm - invertierender Verstärker')

    # --- Fit a low-pass function to get cutoff ---
    # A simple first-order low-pass gain model (in dB):
    # See https://alan.ece.gatech.edu/ECE3040/Lectures/Lecture29-OP%20Amp%20Frequency%20Response.pdf  page 12
    def lowpass_gain(f, A0, fc):
        return 20 * np.log10(np.abs(A0 / np.sqrt(1 + (f / fc) ** 2)))

    popt, _ = curve_fit(lowpass_gain, frequencies, gain_db, p0=[100, 10000])
    A0_fit, fc_fit = popt
    bw_fit = fc_fit  # For first-order low-pass, cutoff == bandwidth

    # Plot fitted curve
    f_fit = np.logspace(np.log10(min(frequencies)),
                        np.log10(max(frequencies)), 500)

    # Plot fitted curve and data
    ax.semilogx(frequencies, gain_db, 'o-',
                label=r"Gemessene Verstärkung $\mathbin{/} \unit{\dB}$")
    ax.semilogx(f_fit, lowpass_gain(f_fit, *popt), 'r--',
                label=rf"Fit: $A_0={{{A0_fit:.1f}}}$, $f_c={{{fc_fit:.1f}}} \unit{{\hertz}}$")
    ax.legend()

    # Print fit results
    print(f"Angepasste Verstärkung (A0): {A0_fit:.2f}")
    print(f"Grenzfrequenz (fc): {fc_fit:.2f} Hz")
    print(f"Bandbreite: {bw_fit:.2f} Hz")
    print(f"Erwartete Verstärkung -R2/R1: {R2/R1}")
    print()

    return fig, ax, A0_fit, fc_fit, bw_fit


def makePlot2(measurement, mode="integrator"):
    frequencies = measurement.frequencies
    voltages_out = measurement.voltages
    Vin = measurement.Vin
    R = measurement.R
    C = measurement.C

    # Compute Gain (linear and dB)
    gain = np.abs(voltages_out / Vin)
    gain = np.clip(gain, 1e-12, None)  # prevent log10(0)
    gain_db = 20 * np.log10(gain)

    f_fit = np.logspace(np.log10(min(frequencies)),
                        np.log10(max(frequencies)), 500)

    if mode == "integrator":

        def model(f, A0):
            return 20 * np.log10(np.abs(A0 / f))

        # Divide Vin by two to get the Amplitude (Oscilicope behaves this way)
        A0_theo = 1 / (2 * np.pi * R * C)
        title = "Bode-Diagramm - Integrator"
    elif mode == "differentiator":

        def model(f, A0):
            return 20 * np.log10(np.abs(A0 * f))

        A0_theo = 2 * np.pi * R * C
        title = "Bode-Diagramm - Differenzierer"
    else:
        raise ValueError(
            "mode must be either 'integrator' or 'differentiator'")

    popt, _ = curve_fit(model, frequencies, gain_db, p0=[A0_theo])
    A0_fit = popt[0]

    fig, ax = plt.subplots(layout="constrained")
    ax.grid(which='both', linestyle='--')
    ax.set_ylabel(r"Verstärkung $\mathbin{/} \unit{\dB}$")
    ax.set_xlabel(r"Frequenz $\mathbin{/} \unit{\hertz}$")
    ax.set_title(title)
    ax.semilogx(f_fit, model(f_fit, A0_fit), 'r--',
                label=f"Fit: $A_0 = {{{A0_fit:.2e}}}$")
    ax.semilogx(f_fit, model(f_fit, A0_theo), 'g:',
                label=rf"Theorie: $A_0 = {{{A0_theo:.2e}}}$")
    ax.semilogx(frequencies, gain_db, 'o',
                label=r"Gemessene Verstärkung $\mathbin{/} \unit{\dB}$")
    ax.legend()
    print("Abweichung vom Theoretischen Wert:", (A0_theo - A0_fit) / A0_theo)
    return fig, ax, A0_fit


# Inverting Amp 1
inv_1 = readData(
    filepath="data/inv_amp_1",
    unit_freq="kHz",
    unit_volt="V",
    experiment="inv_op_amp",
    Vin=0.3,
    R1=1e3,
    R2=100e3,
)

fig, ax, A0_fit, fc_fit, bw_fit = makePlot1(inv_1)
fig.savefig("build/inv_1.pdf")


# Inverting Amp 2
inv_2 = readData(
    filepath="data/inv_amp_2",
    unit_freq="kHz",
    unit_volt="V",
    experiment="inv_op_amp",
    Vin=0.1,
    R1=1e3,
    R2=150e3,
)
fig, ax, A0_fit, fc_fit, bw_fit = makePlot1(inv_2)
fig.savefig("build/inv_2.pdf")

# Inverting Amp 3
inv_3 = readData(
    filepath="data/inv_amp_3",
    unit_freq="kHz",
    unit_volt="mV",
    experiment="inv_op_amp",
    Vin=0.1,
    R1=1e3,
    R2=10e3,
)
fig, ax, A0_fit, fc_fit, bw_fit = makePlot1(inv_3)
fig.savefig("build/inv_3.pdf")

# Integrator
R, C = 10e3, 100e-9
integrator = readData(
    filepath="data/integrator",
    unit_freq="Hz",
    unit_volt="V",
    experiment="int_diff",
    Vin=1,
    R=R,
    C=C,
)
tau = R * C
print(f"RC Time Constant τ = {tau:.2e} s")
fig, ax, A0_fit = makePlot2(integrator, "integrator")
fig.savefig("build/integrator.pdf")

# Differentiator
R, C = 100e3, 22e-9
differentiator = readData(
    filepath="data/differentiator",
    unit_freq="kHz",
    unit_volt="V",
    experiment="int_diff",
    Vin=90e-3,
    R=R,
    C=C,
)
tau = R * C
print(f"RC Time Constant τ = {tau:.2e} s")
fig, ax, A0_fit = makePlot2(differentiator, "differentiator")
fig.savefig("build/differentiator.pdf")
print()

# Schmitt Trigger
R_1 = 10e3
R_2 = 100e3
V_in = 2.989
V_out = 28.1
T = R_2 / R_1 * V_in
abweichung = (T - V_out) / T
print("Der Theoretische Schwellenwert in Volt ist", T)
print("Der Experimentelle Schwellenwert in Volt ist", V_out)
print("Abweichung von", abweichung)
print()

# Generator
R_1 = 10e3
R_2 = 100e3
R_3 = 1e3
C = 1e-6
f_theo = R_2 / (4*C*R_1*R_3)
f_exp = 1.67e3
# amplitude_theo = * R_1 / R_2
# amplitude_exp =
print("Die Theoretische Frequenz ist", f_theo)
print("Die Experimentelle Frequenz ist", f_exp)
print("Die Abweichung beträgt", (f_theo - f_exp) / f_theo)
print()

# Generator with variing amplitude
R, C = 10e3, 100e-9
generator = pd.read_csv("data/scope_10.csv")

t = generator["x-axis"]  # Time in seconds
U = generator["2"]  # Voltage in Volts

mask = (t < 0.025) * (t > 0)  # Select the first 25ms

t = t[mask]  # Select the first 25ms
U = U[mask]  # Select the first 25ms

U = U - np.mean(U)  # Remove DC offset
U_abs = np.abs(U)  # Take absolute value

peaks, _ = find_peaks(U_abs, height=0, distance=10)  # Find peaks above 0V

U_ln = np.log(U_abs.iloc[peaks])  # Take the logarithm of the peaks


# Fit the peaks with a linear function
def generator_model(t, a, b):
    return a * t + b  # Exponential decay


popt, pcov = curve_fit(generator_model, t.iloc[peaks], U_ln)
a, b = unc.correlated_values(popt, pcov)
print(f"Fitted parameters for generator: a = {a:.2f}, b = {b:.2f}")
print(f"Tau = {1/a} s")  # Time constant
print(f"Tau_theo = {20 * R * C}")
print(f"T_theo = {2 * np.pi * R * C} s")  # Time constant


fig, ax = plt.subplots(layout="constrained")
ax.plot(t.iloc[peaks] * 1e3, U_ln, "o", label="Generator Ausgangsamplitude")
ax.plot(t.iloc[peaks] * 1e3,
        generator_model(t.iloc[peaks], *popt), label="Fit")
ax.set(
    xlabel=r"$t \mathbin{/} \unit{\milli\second}$",
    ylabel=r"$\ln\left(\frac{A}{\unit{\volt}}\right)$",
)
ax.grid(which="both", linestyle="--")
ax.legend()
fig.savefig("build/generator_tau.pdf")


fig, ax = plt.subplots(layout="constrained")
ax.plot(t * 1e3, U * 1e3, label="Generator Ausgangsamplitude")
fig.savefig("build/generator.pdf")

# write the peaks in a txt file

data = np.array([t.iloc[peaks] * 1e3, U_abs.iloc[peaks] * 1e3, U_ln])
np.savetxt("build/peaks.txt", np.column_stack(data), fmt="%.2f")

# plot csv data

csv_data = pd.read_csv("data/scope_10.csv")
t = csv_data["x-axis"]  # Time in seconds
U1 = csv_data["1"]  # Voltage in Volts
U2 = csv_data["2"]  # Voltage in Volts

fig, ax = plt.subplots(layout="constrained")
ax.plot(t, U1, label="U1")
ax.plot(t, U2, label="U2")
ax.set(
    xlabel=r"$t \mathbin{/} \unit{\milli\second}$",
    ylabel=r"$U \mathbin{/} \unit{\volt}$",
)
ax.grid(which="both", linestyle="--")
ax.legend()
fig.savefig("build/scope_7.pdf")
