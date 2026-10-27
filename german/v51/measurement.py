import numpy as np

class Measurement:
    def __init__(self, frequencies, voltages, unit_freq, unit_volt, experiment, Vin, R1=None, R2=None, R=None, C=None):
        self.frequencies = self.convert_frequencies(frequencies, unit_freq)
        self.voltages = self.convert_voltages(voltages, unit_volt)
        self.Vin = Vin
        self.experiment = experiment

        if experiment == "inv_op_amp":
            if R1 is None or R2 is None:
                raise ValueError("R1 and R2 must be provided for inverting amp.")
            self.R1 = R1
            self.R2 = R2
        elif experiment == "int_diff":
            if R is None or C is None:
                raise ValueError("R and C must be provided Integrator or Differentiator.")
            self.R = R
            self.C = C
        else:
            raise ValueError("Unknown experiment type. Choose 'inv_op_amp' or 'int_diff'.")

    @staticmethod
    def convert_frequencies(freqs, unit):
        unit_factors = {
            'hz': 1,
            'khz': 1e3,
            'mhz': 1e6
        }
        factor = unit_factors.get(unit.lower())
        if factor is None:
            raise ValueError(f"Unknown frequency unit: {unit}")
        return factor * freqs

    @staticmethod
    def convert_voltages(volts, unit):
        unit_factors = {
            'v': 1,
            'mv': 1e-3
        }
        factor = unit_factors.get(unit.lower())
        if factor is None:
            raise ValueError(f"Unknown voltage unit: {unit}")
        return factor * volts

    def __repr__(self):
        return f"<Measurement with {len(self.frequencies)} points, Vin={self.Vin}V, Experiment={self.experiment}>"
    
def readData(filepath, unit_freq, unit_volt, experiment, Vin, R1=None, R2=None, R=None, C=None):
    volts, freqs = np.genfromtxt(filepath, unpack=True)

    return Measurement(
        frequencies=freqs,
        voltages=volts,
        unit_freq=unit_freq,
        unit_volt=unit_volt,
        experiment=experiment,
        Vin=Vin,
        R1=R1,
        R2=R2,
        R=R,
        C=C
    )
