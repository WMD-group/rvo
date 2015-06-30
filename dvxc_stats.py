#! /usr/bin/env python

# dvxc_stats.py
# Copyright 2014 Adam J. Jackson and Jonathan M. Skelton

from __future__ import print_function

import os
import csv
import math
import numpy as np
from collections import namedtuple

import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib.ticker import FuncFormatter

import dvxc
from dvxc import EVPerCubicAngstromInGPa

### set up data files ###

# src_dir = os.path.dirname(__file__)
# # Append a trailing slash to make coherent directory name - this would select
# # the root directory in the case of no prefix, so we need to check
# if src_dir:
#     src_dir = src_dir + '/'
# data_dir = src_dir + '../data/'

### Plotting configuration ###
functional_colours = {
    'LDA'     : (  0,   0, 204),
    'PW91'    : (204,   0,   0),
    'PBE'     : (  0, 153,   0),
    'PBEsol'  : (112,  48, 160),
    'TPSS'    : (255, 192,   0),
    'revTPSS' : (214,   0, 147),
    'PBE+D2'  : (  0, 176, 240),
    'B3LYP'   : (100, 100,  50),
    'HSE06'   : ( 50, 80,  50)
    }

functional_markers = { 'LDA' : "^", 'PW91' : "o", 'PBE' : "s", 'PBEsol' : "D", 
                       'TPSS' : "x", 'revTPSS' : "p", 'PBE+D2' : "+" ,
                       'B3LYP' : "<", 'HSE06' : ">"}

iter_plot_bar_hatches = [None, "//", "/", "---", 'o', 'oo']

# Global initialisation of Matplotlib settings and variables/objects; the font
# size is kept as a variable, because calls to plt.legend() don't seem to use
# the rc settings.  COMMENT: Note that the line width *is* hard-coded in
# places, so changing the variable here won't update everything on the plots.

mpl.rc('font', **{ 'family' : 'sans-serif', 'size' : 10, 
                   'sans-serif' : 'Arial' })
mpl.rc('lines', **{ 'linewidth' : 0.5 })
font_size = 10

global data_dir

def main(verbosity=False, to_plot=["none"], data_dir="",
         to_write=["none"], compounds=False, compact=False):
    
    if data_dir[-1] != '/':
        data_dir += '/'
    
    all_data_sets = [
        ("PbS", data_dir + u"PbS-EVPData.csv"),
        ("PbTe", data_dir + u"PbTe-EVPData.csv"),
        ("ZnS", data_dir + u"ZnS-EVPData.csv"),
        ("ZnTe", data_dir +u"ZnTe-EVPData.csv")
        ]

    
    ### Set up verbose printing ###
    if verbosity:
        def vprint(*args):
            for arg in args:
                print(arg, end=' ')
                print('')
    else:
        def vprint(*args):
            pass

    ### Let "none" override other output selections
    if "none" in to_plot:
        to_plot = ["none"]
    if "none" in to_write:
        to_write = ["none"]

    ### Trim down data_sets if requested: ###
    # global data_sets # Python is weird about global variables:
    #                  # this line is needed by the "if" statement but not the
    #                  # "for" loop!
    if compounds:
        data_sets = [entry for entry in all_data_sets if entry[0] in compounds]
    else:
        data_sets = all_data_sets

    # Main loop: Simulated DVXC over selected materials
    for output_prefix, filename in data_sets:
        ### Import data and find optimal volumes ###
        print("Processing: {0}\n".format(filename))
        vprint("  -> Reading input data...")
        functionals, data = read_dataset(filename)
        vprint("  -> Determining Veq by Murnaghan fitting...")
        eos_fits = {}
        for functional in functionals:
            murnaghan_params = dvxc.murnaghan_fit(data[functional].e_values,
                                                  data[functional].v_values)
            vprint("    -> {0}: {1:.2f} A^3 (RMS = {2:.2e})".format(
                functional, murnaghan_params.v0, murnaghan_params.eRMS))
            eos_fits[functional] = murnaghan_params
        vprint("\n")

        if "eos" in to_plot or "all" in to_plot:
            plot_filename = "{0}_EoS.png".format(output_prefix)
            print("  -> Plotting {0}...".format(plot_filename))
            plot_EVP_data(functionals, data, eos_fits, 
                          plot_filename=plot_filename, compact=compact)
        if "eos" in to_write or "all" in to_write:
            out_filename = "{0}_EoS.csv".format(output_prefix)
            print("  -> Outputting {0}...".format(out_filename))
            with open(out_filename, 'w') as f:
                csv_writer = csv.writer(f, delimiter = ',',quotechar = '\"',
                                        quoting = csv.QUOTE_ALL)
                csv_writer.writerow(["Functional", "E0 / eV", "V0 / A^3", 
                                     "K0 / GPa", "K'0", "RMS Error"])
                for functional in functionals:
                    fit = eos_fits[functional]
                    csv_writer.writerow([functional, fit.e0, fit.v0,
                                         fit.k0 * EVPerCubicAngstromInGPa,
                                         fit.kPrime0, fit.eRMS])

        ### Carry out delta p xc method for all functionals ###
        print("  -> Applying Delta V xc correction to all " + 
              "combinations of functionals...")
        DvxcResult = namedtuple('DvxcResult', 'functional delta_v residual_p')
        dvxc_results={}
        for test_functional in functionals:
            initial_v_values = data[test_functional].v_values
            initial_p_values = data[test_functional].p_values
            delta_v = {}
            residual_p = {}
            # Compare with all other functionals
            for ref_functional in (f for f in functionals if f != functional):
                dv = [dvxc.apply_dvxc_murnaghan(p, eos_fits[ref_functional]) \
                      for p in initial_p_values]
                p = [dvxc.murnaghan_pressure(v + dv_v,
                                             eos_fits[test_functional]) \
                     for (v, dv_v) in zip(initial_v_values, dv)]
                delta_v.update({ref_functional: dv})
                residual_p.update({ref_functional: p})     
                dvxc_results.update({test_functional: DvxcResult(
                    test_functional, delta_v, residual_p)})

        if "dvxc" in to_write:
            ### write out table of delta p xc corrections ###
            out_filename = "{0}_dvxc.csv".format(output_prefix)
            vprint("  --> Writing to file {0}...".format(out_filename))
            with open(out_filename, 'w') as f:        
                csv_writer = csv.writer(f, delimiter = ',')
                csv_writer.writerow(["Test functional", "Reference functional",
                                     "v / AA^3", "p / kbar", "Delta V / AA^3", 
                                     "residual pressure / kbar"])
                for test_functional in functionals:
                    for ref_functional in dvxc_results[test_functional].delta_v:
                        # Zip together numbers and concatenate to labels
                        for row in zip(
                            data[test_functional].v_values,
                    data[test_functional].p_values,
                    dvxc_results[test_functional].delta_v[ref_functional],
                    dvxc_results[test_functional].residual_p[ref_functional]
                        ):
                        # print(row)
                        # print(list(row))
                            csv_writer.writerow([test_functional, 
                                                 ref_functional] +
                                                list(row))
      
        # Iterative procedure     
        print("  -> Applying correction iteratively to all combinations")
        iter_results = {}
        for functional in functionals:
            iter_results.update({functional:{}})
            for ref_functional in (f for f in functionals if f!= functional):
                iter_result = sim_iterative_dvpx(eos_fits[functional],
                                                 eos_fits[ref_functional],
                                                 max_iter=5)
                iter_results[functional].update({ref_functional:iter_result})

        for option in ("iterative_p", "iterative_v", "all"):
            if option in to_write:
                print("Sorry, data file not yet available for iterations")

        if "iterative_p" in to_plot or "all" in to_plot:
            filename = "{0}_Iterative_P.png".format(output_prefix)
            print("  -> Plotting iterative P results to {0}".format(filename))
            plot_iterative_results(functionals, iter_results, eos_fits,
                                   filename, mode='pressure')
        if "iterative_v" in to_plot or "all" in to_plot:
            filename = "{0}_Iterative_V.png".format(output_prefix)
            print("  -> Plotting iterative V results to {0}".format(filename))
            plot_iterative_results(functionals, iter_results, eos_fits,
                                   filename, mode='volume')

        ### Sensitivity to initial distance ###
        # Calculate volume shift corresponding to a dvxc calculation at each
        # point of one data set, calculate residual pressure associated with
        # this point from fit
        iterations = 4

        # Check if this is required:
        if ("v_sensitivity" in to_plot or
            "all" in to_plot or
            "v_sensitivity" in to_write or
            "all" in to_write):

            for (test_functional, ref_functional) in (
                    ("B3LYP", "PW91"),
                    ("HSE06", "PBE")
            ):
                (v_init, v_corrected) = v_sensitivity(eos_fits[test_functional], 
                                                      eos_fits[ref_functional],
                                                      range_percent=3,
                                                      iterations=iterations)
                if "v_sensitivity" in to_plot or "all" in to_plot:
                    filename = "{0}_{1}_Initial_V_sensitivity.png".format(
                    output_prefix, test_functional)
                    print(("  -> Plotting V sensitivity for {0} " +
                           "(ref. {1}) to {2}").format(test_functional,
                                                       ref_functional,
                                                       filename))
                    plot_v_sensitivity(v_init,v_corrected,
                                       test_functional,ref_functional,
                                       v0=eos_fits[test_functional].v0,
                                       filename=filename)
                if "v_sensitivity" in to_write or "all" in to_write:
                    filename = "{0}_{1}_Initial_V_sensitivity.csv".format(
                    output_prefix, test_functional)
                    print(("  -> Writing V sensitivity data for {0} " +
                           "(ref. {1}) to {2}").format(test_functional,
                                                       ref_functional,
                                                       filename))
                with open(filename,'w') as f:
                    csv_writer = csv.writer(f, delimiter = ',')
                    header = (["v_initial / AA"] + 
                              ["v_corrected (iteration {0})".format(x) 
                               for x in range(1,iterations+1)])
                    csv_writer.writerow(header)
                    for i in range(len(v_init)):
                        csv_writer.writerow(
                            [v_init[i]] + list(v_corrected[i,:])
                        )

        print('\n')

    if "v_bandgap" in to_plot or "all" in to_plot:
        plot_bandgaps([x[0] for x in data_sets], do_plot=True, 
                      plot_filename="v-bandgaps.png", compact=compact)

def read_dataset(filename):
    """Read in data from CSV file

    Arguments: "filename" (text string containing path to file)

    Returns: tuple (functionals, datasets)

        where "functionals" is a list of functional names and "datasets" is a
        set of "DataSet" namedtuple objects containing v_values, e_values and
        p_values (numpy arrays of volumes, energies and pressures in units of
        eV and angstrom)
    
    Required file format:
    Line 1: Functional1,,,Functional2,,, ... ,FunctionalN,,
    Line 2: (unused, contain units for readability)
    Remaining lines (data): volume, energy, p_ext, volume, energy, p_ext ...

    """
    
    with open(filename, 'r') as f:
        input_reader_csv = csv.reader(f)

        # Read functional names from first line
        functionals = next(input_reader_csv)[::3]
        # Skip second line
        next(input_reader_csv)
        #fill datasets
        datasets = {}
        DataSet = namedtuple('DataSet', 'v_values e_values p_values')
        for functional in functionals:
            datasets[functional] = DataSet([],[],[])
            
        for row in input_reader_csv:
            for i, name in enumerate(functionals):
                datasets[name].v_values.append(float(row[i*3]))
                datasets[name].e_values.append(float(row[i*3+1]))
                datasets[name].p_values.append(float(row[i*3+2]))

    # Finally, convert to numpy arrays
    for functional in functionals:
        datasets[functional] = DataSet(
            np.array(datasets[functional].v_values),
            np.array(datasets[functional].e_values),
            np.array(datasets[functional].p_values)
            )
        
    return (functionals, datasets)

def plot_EVP_data(functionals, evp_data, eos_fits, plot_filename=False,
                  compact=False):
    """
    Plot energy-volume and pressure-volume curves

    Arguments:
        functionals: list of named functionals

        evp_data: dict of namedtuple objects containing lists of values:
            {'functional1':DataSet([v1,v2,v3...],[e1,e2,e3...],[p1,p2,p3...],
             'functional2': ...}
            where DataSet = namedtuple('DataSet', 'v_values e_values p_values')

        eos_fits: dict of namedtuple objects containing Murnaghan parameters:
            {'functional1': MurnaghanFit('e0, v0, k0, kPrime0, eRMS'),
             'functional2': ...}  where
            MurnaghanFit = namedtuple('MurnaghanFit', 'e0 v0 k0 kPrime0 eRMS')
    
        plot_filename: File name/path for plotted output. If False, plot is
            displayed onscreen instead.

        compact: [boolean] If True, show 2 plots instead of 4 (for publication)
    
    """
    if compact:
        fig_dimensions = (8 / 2.54, 14 / 2.54)
    else:
        fig_dimensions = (17.2 / 2.54, 12 / 2.54)

    plt.figure(figsize = fig_dimensions)
    subplot_axes = []

    ### Plot 1: E-V curve ###
    if compact:
        plt.subplot(2,1,1)
    else:
        plt.subplot(2,2,1)
    for functional in functionals:
        (r,g,b), marker = (functional_colours[functional], 
                           functional_markers[functional])
        plt.plot(evp_data[functional].v_values,
                 evp_data[functional].e_values - eos_fits[functional].e0,
                 color = (r / 255.0, g / 255.0, b / 255.0),
                 marker = marker, fillstyle = 'none',
                 label=functional)
    plt.xlabel(r"$V$ / $\AA^{3}$", fontweight = 'bold')
    plt.ylabel("$E - E_{0}$ / eV", fontweight = 'bold')

    x_ticks = plt.xticks()[0]
    x_tick_spacing = x_ticks[1] - x_ticks[0]
    y_min, y_max = plt.ylim()
    plt.ylim(0.0, y_max)
    axes=plt.gca()
    subplot_axes.append(axes)

    if compact:
        handles, labels = axes.get_legend_handles_labels()
        # Print legends for first half of list
        # N.B. // is notation for integer floor division
        axes.legend(handles[:len(handles) // 2], labels[:len(handles) // 2],
                loc  = 'upper right', frameon = False, 
                prop = {'size' : font_size})

    
    ### Plot 2: Relative volume ###
    if not compact:
        plt.subplot(2,2,3)
        x_min, x_max = None, None
        for functional in functionals:
            v_values_adjusted = (evp_data[functional].v_values 
                                 - eos_fits[functional].v0)
            e_values_adjusted = (evp_data[functional].e_values 
                                 - eos_fits[functional].e0)
            (r,g,b), marker = (functional_colours[functional], 
                               functional_markers[functional])
            plt.plot(v_values_adjusted, e_values_adjusted,
                     color = (r / 255.0, g /255.0, b / 255.0),
                     marker=marker, fillstyle='none')

            x_min = np.min(v_values_adjusted) if not x_min \
                    else min(x_min, np.min(v_values_adjusted))
            x_max = np.max(v_values_adjusted) if not x_max \
                    else max(x_max, np.max(v_values_adjusted))

        plt.xlabel(r"$V - V_{0}$ / $\AA^{3}$", fontweight = 'bold')
        plt.ylabel("$E - E_{0}$ / eV", fontweight = 'bold')
        plot_x_min = math.floor(x_min / x_tick_spacing) * x_tick_spacing
        plot_x_max = math.ceil(x_max / x_tick_spacing) * x_tick_spacing
        plt.xlim(plot_x_min, plot_x_max)
        y_min, y_max = plt.ylim()
        plt.ylim(0.0, y_max)
        subplot_axes.append(plt.gca())

    ### Plot 3: P-V curve ###
    if compact:
        plt.subplot(2,1,2)
    else:
        plt.subplot(2,2,2)

    def zeroline():
        plt.axhline(y=0,color='k',linestyle='-')
    zeroline()

    for functional in functionals:
        (r,g,b), marker = (functional_colours[functional], 
            functional_markers[functional])
        plt.plot(evp_data[functional].v_values,
            evp_data[functional].p_values / 10,
                     label = functional,
                     color = (r / 255.0, g / 255.0, b / 255.0),
                     marker = marker, fillstyle = 'none' )
       
    plt.xlabel(r"$V$ / $\AA^{3}$", fontweight = 'bold')
    plt.ylabel("$p_{ext}$ / GPa", fontweight = 'bold')
    axes = plt.gca()
    handles, labels = axes.get_legend_handles_labels()

    if compact:
        # Print legends for second half of list
        axes.legend(handles[len(handles) // 2:], 
                    labels[len(handles) // 2:], 
                    loc = 'upper right', frameon = False, 
                    prop = { 'size' : font_size })
    else:
        # Print legends for first half of list
        # N.B. // is notation for integer floor division
        axes.legend(handles[:len(handles) // 2], labels[:len(handles) // 2],
                loc  = 'upper right', frameon = False, 
                prop = {'size' : font_size})

    subplot_axes.append(axes)

    ### Plot 4: P-V for normalised V ###
    if not compact:
        plt.subplot(2,2,4)
        x_min, x_max = None, None
        zeroline()
        
        for functional in functionals:
            v_values_adjusted = (evp_data[functional].v_values 
                                 - eos_fits[functional].v0)
            (r,g,b), marker = (functional_colours[functional], 
                               functional_markers[functional])
            x_min = np.min(v_values_adjusted) if not x_min \
                    else min(x_min, np.min(v_values_adjusted))
            x_max = np.max(v_values_adjusted) if not x_max \
                    else max(x_max, np.max(v_values_adjusted))

            plt.plot(evp_data[functional].v_values - eos_fits[functional].v0,
                     evp_data[functional].p_values / 10.0,
                     label=functional,
                     color = (r / 255.0, g /255.0, b / 255.0),
                     marker=marker, fillstyle='none')
   
        plt.xlabel(r"$V - V_{0}$ / $\AA^{3}$", fontweight = 'bold')
        plt.ylabel("$p_{ext}$ / GPa", fontweight = 'bold')
        plot_x_min = math.floor(x_min / x_tick_spacing) * x_tick_spacing
        plot_x_max = math.ceil(x_max / x_tick_spacing) * x_tick_spacing
        plt.xlim(plot_x_min, plot_x_max)

        axes = plt.gca()
        handles, labels = axes.get_legend_handles_labels()
        axes.legend(handles[len(handles) // 2:], 
                    labels[len(handles) // 2:], 
                    loc = 'upper right', frameon = False, 
                    prop = { 'size' : font_size })
        subplot_axes.append(axes)

    ### Finish up plots ###

    for axes in subplot_axes:
        for spine in axes.spines.values():
            spine.set_linewidth(0.5)
    plt.tight_layout()
    if plot_filename:
        plt.savefig(plot_filename, format = 'png', dpi = 300)
    else:
        plt.show()
    plt.close()

def sim_iterative_dvpx(test_functional_murnaghan_fit,
                       ref_functional_murnaghan_fit,
                       max_iter = 5):
    """Simulate the iterative application of delta V xc procedure

    Rather than carry out new calculations for every combination of
    functionals, we use a Murnaghan fit for both the test functional and the
    reference functional. Calculated pressure values from the reference
    functional fit are a stand-in for what would be a quantum chemical
    calculation at the volume predicted by the delta V xc procedure.

    Arguments:

        test_functional_murnaghan_fit: a MurnaghanFit namedtuple, containing
            the fitted parameters for a Murnaghan EoS describing the P-V curve
            to be minimised. (i.e. results from an expensive DFT functional.)

        ref_functional_murnaghan_fit: a corresponding MurnaghanFit namedtuple,
            parameterising an EoS to be used in the delta V_xc correction.
            (i.e. results from an inexpensive DFT functional.)

        max_iter: Integer number of iterative steps to carry out.

    Returns:

        v_iterations: A list [v0, v1, v2 ... ] of volumes corresponding to the
           iterative steps taken.
    
        p_iterations: A list [p0, p1, p2 ... ] of pressures corresponding to
           iterative steps taken. Beware that these are in the units of the
           EoS, typically eV AA^-3.

    """ 
    try:
        max_iter = int(max_iter)
    except ValueError:
        raise Exception(
            "Couldn't use '{0}' as a number of iterations".format(max_iter))

    # Initialise arrays and variables
    v = ref_functional_murnaghan_fit.v0
    p = dvxc.murnaghan_pressure(v, test_functional_murnaghan_fit)
    v_iterations = [v] + [False]*(max_iter)
    p_iterations = [p] + [False]*(max_iter)

    for index in range(1,max_iter+1):
        delta_v = dvxc.apply_dvxc_murnaghan(p, ref_functional_murnaghan_fit)
        v = v + delta_v
        p = dvxc.murnaghan_pressure(v, test_functional_murnaghan_fit)
        v_iterations[index] = v
        p_iterations[index] = p

    return (v_iterations, p_iterations)

def plot_iterative_results(functionals, iter_results, eos_fits,
                           filename, mode='pressure', n_columns=3):
    """Plot values from iterative application of dvxc

    Arguments:
        functionals: list of functionals to include in plot

        iter_results: dict of dict of tuples containing result lists
            {test_functional: {ref_functional: ([v0,v1,v2...],
                                                [p0,p1,p2...])
                               , ...}
             , ...}

        filename: string containing target plot file path

        mode: string 'volume' or 'pressure', denoting property to plot
              [default value = 'pressure']

        n_columns: width of plot array [default value = 3]

    Raises:
        Exception "Plotting mode [mode] not accepted." if invalid mode

        Exception "Couldn't convert [n_columns] to a meaningful number of
        columns." if n_columns can't be converted to an int.

        Exception "No data for test functional [functional name]"; functional
            appears in list "functionals" but not "iter_results"

        Exception "[test functional] has not been tested with [ref functional]"
            ; [test functional] appears in iter_results, but does not appear
            in the results for [ref functional].

    """
    # Check input
    supported_modes = ('volume', 'pressure')
    if mode not in supported_modes:
        raise Exception("Plotting mode '{0}' not accepted.".format(mode) +
                        "Valid output modes: {0}".format(supported_modes))
    try: 
        n_columns = int(n_columns)
    except ValueError:
        raise Exception("Couldn't convert '{0}' to ".format(n_columns)
                        + "a meaningful number of columns.")

    # Check requested set of functionals are in results
    for test_functional in functionals:
        try:
            test_results = iter_results[test_functional]
        except KeyError:
            raise Exception(
                "No data for test functional {0}".format(test_functional))
        for ref_functional in (x for x in functionals if x != test_functional):
            try:
                ref_result = iter_results[test_functional][ref_functional]
            except ValueError:
                raise Exception("{0} has not been".format(test_functional) +
                                "tested with ref functional " +
                                "{0}".format(ref_functional))

    # Work out number of rows given n_columns
    n_rows, remainder = divmod(len(iter_results), n_columns)
    if remainder: n_rows += 1

    plt.figure(figsize = (25.8 / 2.54, (6.5 * n_rows) / 2.54))
    subplot_axes = []

    # Take number of bars for first combination, assume(!) rest are consistent
    n_bars = len(iter_results[functionals[0]][functionals[1]][0])
    n_groups = len(functionals)

    # Set other plot parameters
    bar_width = 0.8

    # Work out bar positions
    left_edges_base = [(bar_width 
                        + i * n_bars * bar_width + 
                        i * bar_width) for i in range(n_groups)]

    # Get plotting
    for subplot_index, test_functional in enumerate(functionals):
        plot_data = iter_results[test_functional]

    # Modification to (potentially) make the layout look prettier: If the last
    # plot is on its own in a row with an odd number of columns, put it in the
    # middle.

        if (n_columns % 2 != 0 
            and subplot_index == len(functionals) - 1 
            and len(functionals) % n_columns == 1):
            plt.subplot(n_rows, n_columns, 
                        subplot_index + 1 + (n_columns - 1) // 2)
        else:
            plt.subplot(n_rows, n_columns, subplot_index + 1)            
    
        for group_index, ref_functional in enumerate(
            [f for f in functionals if f != test_functional]):
            
            if mode == 'volume':
                volumes = plot_data[ref_functional][0]
                v0 = eos_fits[test_functional].v0
                bar_heights = [abs(v - v0) for v in volumes]
            elif mode == 'pressure':
                bar_heights = [abs(p)/10.0 for p in plot_data[ref_functional][1]]
            
            left_edges = [(left_edges_base[group_index] 
                           + bar_index * bar_width
                       ) for bar_index in range(len(bar_heights))]

            r, g, b = functional_colours[ref_functional]
            bar_group_colour = (r / 255.0, g / 255.0, b / 255.0)

            # Individual bars are plotted separately in order to set the hatch
            # styles.  The "bottom" argument needs to be set to a small (but
            # non-zero) number due to some Matplotlib issues with plt.bar() and
            # log scale y axes. 10^-3 was chosen as 0.001 % / 0.001 GPa is a
            # very high error margin for predicting volumes/pressures!
            lower_lim = 1E-5
              
            for bar_index in range(len(bar_heights)):
                plt.bar(left_edges[bar_index], 
                        bar_heights[bar_index], 
                        bottom = lower_lim, color = bar_group_colour, 
                        edgecolor = bar_group_colour, linewidth = 0.5, 
                        fill = iter_plot_bar_hatches[bar_index] == None, 
                        hatch = iter_plot_bar_hatches[bar_index])
                # Add outline in independent colour
                plt.bar(left_edges[bar_index], bar_heights[bar_index], 
                        bottom = lower_lim, edgecolor = 'k', linewidth = 0.5, fill = False)

            plt.title("Reference: {0}".format(test_functional), 
                      size = font_size, fontweight = 'bold')

        if mode == 'volume':
            plt.ylabel(r"abs($\Delta V$) / %", fontweight = 'bold')
        elif mode == 'pressure':
            plt.ylabel(r"abs($p$) / GPa", fontweight = 'bold')


        plt.xticks([value + (bar_width * n_bars) / 2.0 for value in left_edges_base], 
                   [f for f in functionals if f != test_functional], 
                   rotation = 30, ha='center')

        # plt.xlim(0.0, left_edges_base[-1] + 
        #                 n_bars * bar_width + bar_width)
        plt.xlim(0.0, left_edges_base[-1])

        subplot_axes.append(plt.gca())


    for axes in subplot_axes:
        axes.set_ylim(lower_lim, 10.0)
        axes.set_yscale('log', noposy = 'clip')
        for spine in axes.spines.values():
            spine.set_linewidth(0.5)
        
    plt.tight_layout()
        
    plt.savefig(filename, format = 'png', dpi = 300)

def v_sensitivity(test_functional_fit, ref_functional_fit, 
                  range_percent=3, points=100, iterations=4):
    v0 = test_functional_fit.v0
    v_upper = v0 * (1.+range_percent/100.)
    v_lower = v0 * (1.-range_percent/100.)

    v_init = np.linspace(v_lower,v_upper,points)

    v_corrected = []
    v = v_init.copy()
    for i in range(iterations):
        p = dvxc.murnaghan_pressure(v, test_functional_fit)
        delta_v = dvxc.apply_dvxc_murnaghan(p, ref_functional_fit)
        v = v + delta_v
        v_corrected.append(v)

    return (v_init, np.array(v_corrected).T)

def plot_v_sensitivity(v_init,v_corrected, test_functional,ref_functional,
                       v0=False, filename=False):
    """
    Plot showing convergence over iterations and range of initial volumes


    Arguments:
   
        v_init: 1D numpy array of volume values

        v_corrected: list of 1D numpy arrays, corresponding to iterative
            volume estimates

        test_functional: String containing name of test functional

        ref_functional: String containing name of reference functional

        v0: Optionally provide actual minimum volume. X axis is rescaled to
            show this.

        filename: Target file for plot export. If not provided, graph is drawn
            to screen
    """
    if v0:
        x = v_init/v0
        y = v_corrected/v0
        xlabel = "Initial volume / Optimum volume"
        ylabel = "Corrected volume / Optimum volume"
    else:
        x = v_init
        y = v_corrected
        xlabel = "Initial volume / $\AA^3$"
        ylabel = "Corrected volume / $\AA^3$"

    plt.figure(figsize = (8 / 2.54, 8 / 2.54))
    plt.plot(x, y)

    axes = plt.gca()
    for spine in axes.spines.values():
        spine.set_linewidth(0.5)

    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    axes.yaxis.set_major_formatter(y_formatter)

    plt.title("Corrected functional: {0}\n".format(test_functional) +
              "Reference functional: {0}".format(ref_functional))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(range(1,len(v_corrected)+1), loc='best', title="Iterations")

    plt.tight_layout()
    
    if filename:
        plt.savefig(filename)
    else:
        plt.show()

def plot_bandgaps(materials,do_plot=True,do_write=False,plot_filename=False,
                  compact=False):
    data, kpoints = {}, {}
    BandGapData = namedtuple("BandGapData","factor v eg kpoint")

    markers = {'PbS':'^', 'PbTe':'>', 'ZnS':'o', 'ZnTe':'x'}
    colours  =  {'PbS':'b', 'PbTe':'g', 'ZnS':'r', 'ZnTe':'k'}
    functionals = ["B3LYP","HSE06"]
    for functional in functionals:
        data.update({functional: {}})
        for material in materials:
            with open(data_dir + 
                      "{0} - {1} - Direct Bandgaps.csv".format(material,
                                                               functional),
                      "r") as f:
                csv_reader = csv.reader(f)
                # skip first line
                next(csv_reader)
                # First row
                text_row = next(csv_reader)
                row = [float(x) for x in text_row]
                factor, v, eg = [row[0]], [row[1]], [row[2]]
                kpoint = row[3:]
                for text_row in csv_reader:
                    row = [float(x) for x in text_row]
                    factor.append(row[0])                
                    v.append(row[1])
                    eg.append(row[2])
            data[functional].update({material:BandGapData(factor,v,eg,kpoint)})

    if do_plot:
        if compact:
            fig_dimensions = (8 / 2.54, 14 / 2.54)
        else:
            fig_dimensions = (17.2 / 2.54, 12 / 2.54)
        plt.figure(figsize = fig_dimensions)

        subplot_axes = []
        collect_factors, collect_eg = [], [] # List all values to get ranges

        for i, functional in enumerate(functionals):
            if compact:
                plt.subplot(len(functionals),1,i)
            else:
                plt.subplot(1,len(functionals),i)
            
            for material in materials:
                # Scale factor is in each dimension; cube to obtain volume change
                x = np.power(data[functional][material].factor, 3)
                y = data[functional][material].eg

                plt.plot(x,y, label=material,
                         color=colours[material], marker=markers[material])
                collect_factors += data[functional][material].factor
                collect_eg += data[functional][material].eg


            plt.xlabel(r"$v/v_0$",fontweight = 'bold')
            plt.ylabel("Bandgap / eV",fontweight = 'bold')
            plt.title(functional,fontweight = 'bold')
            plt.legend(loc="best", frameon = False,
                       prop = {'size' : font_size}, ncol=2)

            subplot_axes.append(plt.gca())
        x_min, x_max = float(min(collect_factors))**3, float(max(collect_factors))**3
        y_min, y_max = 0, math.ceil(float(max(collect_eg)))
        
        for axes in subplot_axes:
            axes.set_xlim((x_min, x_max))
            axes.set_ylim((y_min, y_max))
            for spine in axes.spines.values():
                spine.set_linewidth(0.5)

        plt.tight_layout()


        if plot_filename:
            plt.savefig(plot_filename, format='png',dpi=300)
        else:
            plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simulate application of DVXC method to a dataset of binary chalcogenides.")
    parser.add_argument("data_directory", help="path to 'binary_chalcogenides' folder")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("-p", "--plot", 
                        choices=["eos","iterative_p","iterative_v",
                                 "v_sensitivity","v_bandgap", "none", "all"],
                        nargs='*', default=["all"],
                        help="Output plots")
    parser.add_argument("-w", "--write",
                        choices=["eos","dvxc","iterative_p","v_bandgap",
                                 "v_sensitivity","none","all"],
                        nargs='*', default=["all"],
                        help="Write CSV files")
    parser.add_argument("-m", "--material",
                        choices=["PbS","PbTe","ZnS","ZnTe"],
                        nargs='*', default=False)
    parser.add_argument("-c", "--compact", action="store_true",
                        help="Output compact versions of some plots")
    args = parser.parse_args()

    main(data_dir=args.data_directory, verbosity=args.verbose, to_plot=args.plot, to_write=args.write,
         compounds=args.material, compact=args.compact)
    
