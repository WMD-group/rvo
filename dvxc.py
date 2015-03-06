#! /usr/bin/env python

#dvxc.py

# Copyright 2014 Jonathan M. Skelton and Adam J. Jackson

#COMMENT: EoS data should be a CSV file with each row containing (V, E, p), and no header rows.
#COMMENT: We use VASP, so we assume V is in A^3, E is in eV, and p is in kbar.

#COMMENT: An example :-
# For reference, with the PBEsol EoS ("PbS-EoS-PBEsol.csv"), v0 = 50.86 A^3.
# With HSE 06, the pressure at this volume is 26.24 kbar (2.62 GPa).
# One iteration of dvxc gives a Delta V_xc of 2.28 A^3, and hence a predicted HSE 06 equilibrium volume of 53.14 A^3.
# The correct value (from a HSE 06 EoS fitting) is 52.97 A^3, and the pressure at the corrected volume is -2.03 kbar (-0.2 GPa).


import csv;
import math;

import numpy as np;
import matplotlib as mpl;
import matplotlib.pyplot as plt;

from matplotlib.ticker import FuncFormatter;
from scipy.optimize import curve_fit;

from collections import namedtuple

import interpolate_cell

_EVPerCubicAngstromInGPa = 160.217656;


mpl.rc('font', **{ 'family' : 'sans-serif', 'size' : 10, 'sans-serif' : 'Arial' });
mpl.rc('lines', **{ 'linewidth' : 0.5 });

_FontSize = 10;

_Label1DpFormatter = FuncFormatter(lambda x, pos : "{0:.1f}".format(x));


def murnaghan_fit(eData, vData, initialK0Guess = 50.0, initialK0PrimeGuess = 5.0):
    def _fit_equation(x, e0, v0, k0, kPrime0):
        return e0 + k0 * v0 * ((((x / v0) ** (1.0 - kPrime0)) / (kPrime0 * (kPrime0 - 1)))
                               + (x / (kPrime0 * v0)) - (1 / (kPrime0 - 1)));

    minIndex = np.argmin(eData);

    e0Init = eData[minIndex];
    v0Init = vData[minIndex];

    (e0, v0, k0, kPrime0), pcov = curve_fit(_fit_equation, vData, eData, [e0Init, v0Init, initialK0Guess, initialK0PrimeGuess]);
    eRMS = math.sqrt(np.mean((eData - _fit_equation(vData, e0, v0, k0, kPrime0)) ** 2));

    MurnaghanFit = namedtuple('MurnaghanFit', 'e0 v0 k0 kPrime0 eRMS')
    return MurnaghanFit(e0, v0, k0, kPrime0, eRMS);

def murnaghan_pressure(v, fit):
    """Calculate the pressure given a volume and fitted parameters for a Murnaghan EOS

    Arguments:
        v: unit cell volume in cubic angstroms (AA^3)

        fit: namedtuple containing fit parameters in units of eV and AA^3
            format: namedtuple('MurnaghanFit', 'e0 v0 k0 kPrime0 eRMS')

    Returns:
        p: pressure in kbar (10^8 Pa)

    p(v) = (k0/k'0) * [(v0/v) ^ k'0  -1 ]

    """
    p = (fit.k0/fit.kPrime0) * ((fit.v0 / v) ** fit.kPrime0 - 1)
    p = p * _EVPerCubicAngstromInGPa * 10 # Convert to kbar / A3
    return p
    

def _ApplyDVxc_linear(pressure, vpPolyFit):
    #COMMENT: pressure is the pressure of the "target" functional.
    #COMMENT: vpPolyFit is the coefficients of a polynomial fit to *volume* as a function of *pressure*.
    return np.polyval(vpPolyFit, 0.0) - np.polyval(vpPolyFit, pressure);

def apply_dvxc_murnaghan(pressure, murnaghan_params):
    """
    Estimate volume change to minimise absolute pressure, given reference fit to Murnaghan EoS

    Arguments:
        pressure: External pressure in units of kPa (float)

        murnaghan_params: namedtuple containing fit parameters in units of eV and AA^3
            format: namedtuple('MurnaghanFit', 'e0 v0 k0 kPrime0 eRMS')

    Returns:
        Volume change in AA^3

    From Murnaghan EOS, dP/dv = (-k0/v) * (/v0)**(-kprime0)
    Tangent at P=0 therefore has slope m = DP/Dv = (-k0/v0)
    v0-v ~= (0-P)/m
    
    """
    m = -murnaghan_params.k0 / murnaghan_params.v0
    m = m *_EVPerCubicAngstromInGPa * 10 # Convert to kbar / A3
    return -pressure/m

def _SaveVPFitPlot(pValues, vValues, murnaghan_params, fileName):

    plt.figure(figsize = (8.6 / 2.54, 6 / 2.54))

    plt.plot(pValues / 10.0, vValues, color = 'k', marker = '^', markersize = 5, label = "EoS", fillstyle = 'none')

    xMin, xMax = plt.xlim()

    fitOverlayValues = np.linspace(xMin * 10.0, xMax * 10.0, 100)
    
    e0, v0, k0, kPrime0, eRMS = murnaghan_params
    murnaghan_pValues = (k0/kPrime0)*((vValues/v0)**(-kPrime0) -1) * _EVPerCubicAngstromInGPa
    plt.plot(murnaghan_pValues, vValues, color = 'k', marker = 'x', markersize = 3,
             linestyle = 'none', label = "Murnaghan fit")
    
    plt.xlabel(r"$p$ / GPa", fontweight = 'bold')
    plt.ylabel(r"$V$ / $\AA^{3}$", fontweight = 'bold')

    plt.legend(loc = 'upper right', frameon = False, prop = { 'size' : _FontSize })

    axes = plt.gca()

    for spine in axes.spines.values():
        spine.set_linewidth(0.5)

    plt.tight_layout()

    plt.savefig(fileName, format = 'png', dpi = 300)

    plt.close()


def import_VEP(EoSData):
    vValues, eValues, pValues = [], [], []

    with open(EoSData, 'r') as inputReader:
        inputReaderCSV = csv.reader(inputReader)

        for row in inputReaderCSV:
            if row[0][0] != "#":
                vValue, eValue, pValue = row

                vValues.append(float(vValue))
                eValues.append(float(eValue))
                pValues.append(float(pValue))

    vValues, eValues, pValues = np.array(vValues), np.array(eValues), np.array(pValues)

    return (vValues, eValues, pValues)

def import_VE(EoSData):
    vValues, eValues = [], []

    with open(EoSData, 'r') as inputReader:
        inputReaderCSV = csv.reader(inputReader)

        for row in inputReaderCSV:
            if row[0][0] != "#":
                vValue, eValue = row

                vValues.append(float(vValue))
                eValues.append(float(eValue))

    vValues, eValues = np.array(vValues), np.array(eValues)

    return (vValues, eValues)

def import_lattice(EoSData):
    """"
    Import data from 10-column CSV file in format:
    Energy, ax,ay,az,bx,by,bz,cx,cy,cz

    arguments:
        EoSData: CSV file path as string

    returns:
        vValues: list of volumes in cubic angstroms
        eValues: list of energies in eV
        lattice_vectors: list of lattice vectors as 3x3 numpy arrays

    """

    # Double precision is used here for the all-electron energy TODO
    # Is the program losing precision in places due to use of floats?
    # Floats only have ~ 7.d.p of precision, fine for normalised
    # energies but not for total energies. Perhaps a warning is needed
    # for large energies?
    
    vValues, eValues, lattice_vectors = [], [], []
    with open(EoSData, 'r') as inputReader:
        inputReaderCSV = csv.reader(inputReader)

        for row in inputReaderCSV:
            if row[0][0] != "#":
                eValue = np.double(row[0])
                eValues.append(eValue)
                lattice = np.array(row[1:10]).reshape([3,3]).astype('float')
                lattice_vectors.append(lattice)
                vValue = interpolate_cell.cell_volume(lattice)
                vValues.append(vValue)
    return(vValues, eValues, lattice_vectors)


def estimate_vectors(lattice_vectors, target_volume,
                     vValues, verbose=False):
    if verbose:
        def vprint(*args):
            for arg in args:
                print(arg)
    else:
        def vprint(*args):
            pass

    vprint("  -> Proposed volume of improved cell:" +
           "{0:.3f} A^3 ".format(target_volume))
            
    vprint("  -> Estimating new lattice vectors:")

    # We wrap the sorted volumes with their original indices
    # e.g. [(v1,1),(v3,3),(v2,2)]
    sorted_volumes = zip(vValues,range(len(vValues)))
    sorted_volumes.sort()
    for i in range(len(sorted_volumes)):
        if (sorted_volumes[i][0] <= target_volume and
            sorted_volumes[i+1][0] > target_volume):
            cell1 = lattice_vectors[sorted_volumes[i][1]]
            cell2 = lattice_vectors[sorted_volumes[i+1][1]]
            break
    else: # The target must lie outside the reference range
        print("    (Extrapolating beyond existing volume range.)\n")
        if target_volume < sorted_volumes[0][0]:
            cell1 = lattice_vectors[sorted_volumes[0][1]]
            cell2 = lattice_vectors[sorted_volumes[1][1]]
        elif target_volume >= sorted_volumes[-1][1]:
            cell1 = lattice_vectors[sorted_volumes[-2][1]]
            cell2 = lattice_vectors[sorted_volumes[-1][1]]
        else:
            raise Exception("Trouble finding suitable reference cells")

    vprint("    Lattice vectors used for interpolation:")
    vprint("    v = {0:5f} A^3\n".format(
        interpolate_cell.cell_volume(cell1)))
    vprint(pretty_vectors(cell1, padding=8))
    vprint("    v = {0:5f} A^3\n".format(
        interpolate_cell.cell_volume(cell2)))
    vprint(pretty_vectors(cell2, padding=8))

    # Solve interpolation factor
    Lambda = interpolate_cell.solve_v(cell1,cell2,target_volume)
    vprint("    Interpolation factor Lambda = {0:4f}".format(Lambda))
    new_cell = interpolate_cell.interpolate_cell(cell1,cell2,Lambda)
    vprint("    Suggested lattice vectors for next iteration:\n")
    vprint(pretty_vectors(new_cell, padding=8))
    return new_cell

def pretty_vectors(cell,padding=0):
    """POSCAR-friendly formatting for lattice vector matrices"""
    pretty_string=""
    for row in range(3):
        pretty_string += "{0}{1:4f}  {2:4f}  {3:4f}\n".format(
            " "*padding, cell[row,0], cell[row,1], cell[row,2])
    return pretty_string
                
def main(filename, pressure, current_volume):
    #EoSData = "PbS-EoS-PBEsol.csv";
    #TestFunctionalPressure = 26.24;

    EoSData = filename
    TestFunctionalPressure = pressure

    print("Reading EoS from \"{0}\"...".format(EoSData))
    # Check input file format: note that row[0][0] is used to check
    # for comment character as CSV reader always returns a list.
    with open(EoSData,'r') as f:
        inputReaderCSV = csv.reader(f)
        for row in inputReaderCSV:
            if row[0][0] != '#':
                data_file_width = len(row)
                break
    if data_file_width == 3:        
        print("Importing data from Volume-Energy-Pressure table")
        vValues, eValues, pValues = import_VEP(EoSData)
        lattice_vectors = False

        headerLine = "{0: >8}    {1: >7}    {2: >6}".format(
                      "V / A^3", "E / eV", "p / GPa")
        print(headerLine)
        print("-" * len(headerLine))
        for v, e, p in zip(vValues, eValues, pValues):
            print("{0: >8.2f}    {1: >7.2f}    {2: >6.2f}".format(
                   v, e, p))
        print('');
    elif data_file_width == 2:
        print("Importing data from Volume-Energy table")
        vValues, eValues = import_VE(EoSData)
        pValues, lattice_vectors = False, False

        headerLine = "{0: >8}    {1: >7}".format("V / A^3", "E / eV")
        print(headerLine)
        print("-" * len(headerLine))
        for v, e in zip(vValues, eValues):
            print("{0: >8.2f}    {1: >7.2f}".format(v, e))
        print('')
    elif data_file_width == 10:
        print("Importing data from Energy-lattice vectors table" +
              " and calculating volumes")
        vValues, eValues, lattice_vectors = import_lattice(EoSData)
        pValues = False

        headerLine = "{0: >8}    {1: >7}".format("V / A^3", "E / eV")
        print(headerLine);
        print("-" * len(headerLine));
        for v, e in zip(vValues, eValues):
            print("{0: >8.2f}    {1: >7.2f}".format(v, e));
        print('')
    else:
        print("Data file type not recognised")
        return False
    print('');


    print("Performing EoS fit...");

    murnaghan_params = murnaghan_fit(eValues, vValues);
    print("  -> V0 = {0:.2f} A^3".format(murnaghan_params.v0));
    print("  -> E0 = {0:.2f} eV".format(murnaghan_params.e0));
    print("  -> K0 = {0:.2f} GPa".format(murnaghan_params.k0 * _EVPerCubicAngstromInGPa));
    print("  -> K'0 = {0:.2f}".format(murnaghan_params.kPrime0));
    print("  -> RMS = {0:.2e}".format(murnaghan_params.eRMS));

    print('');

    print("Performing dVxc...");


    DVxc_Murnaghan = apply_dvxc_murnaghan(pressure, murnaghan_params)
    print("  -> Your volume offset (\"DeltaVxc\") by Murnaghan fitting is:")
    print("      ---> ** DVxc = {0:.2f} A^3 ** <---".format(DVxc_Murnaghan))
    print('')

    if current_volume and not lattice_vectors:
        print("  -> Proposed volume of improved cell:" +
              "{0:.3f} A^3 ".format(current_volume + DVxc_Murnaghan))
    elif lattice_vectors and not current_volume:
        print("You have provided lattice vectors as part of the \n"
              "reference data. If you provide the volume of the \n"
              "test calculation with the -v argument, this program \n"
              "can estimate new lattice vectors by interpolating to\n"
              "the target volume. \n")
    elif lattice_vectors and current_volume:
        target_volume = current_volume + DVxc_Murnaghan
        new_cell = estimate_vectors(lattice_vectors, target_volume,
                                    vValues, verbose=True)
    
    if pValues != False:
        print("Saving a plot of the dVxc EoS fit to \"Fit.png\"...");
        _SaveVPFitPlot(pValues, vValues, murnaghan_params, "Fit.png");
    
        print('');        

    print("Thank you for using dVxc!");
    print('');
    
    print("There isn't a dVxc paper yet, but in the meantime you could cite:");
    print("J. Buckeridge et al., Comp. Phys. Commun. 185, 330-338 (2014)");
    print('');

if __name__ == "__main__":
    try:
        import readline
        readline.parse_and_bind("tab: complete")
    except:
        pass

    import argparse
    parser = argparse.ArgumentParser(description="Delta V_xc; affordable volume corrections for high-level DFT calculations")
    parser.add_argument('-f', '--file', action="store", type=str,
                        dest="filename", default=False,
                        help=("csv file containing E-V curve with functional A"+
                              " in units of eV and cubic angstroms"))
    parser.add_argument('-p', '--pressure', action="store", type=float,
                        dest="pressure", default=False,
                        help="Pressure at current volume for functional B in GPa")
    parser.add_argument('-v', '--volume', action="store", type=float,
                        dest="current_volume", default=False,
                        help="Value of current volume; used to estimate new lattice vectors when suitable input file is provided.")

    args = parser.parse_args()

    if args.filename:
        filename = args.filename
    else: 
        filename = raw_input("Please give path/filename of CSV file containing E-V curve: ")

    if args.pressure:
        pressure = args.pressure
    else:
        pressure = float(raw_input("Please input pressure at current volume in GPa: "))

    main(filename,pressure,args.current_volume)
