# This is a hack for now, so I can run the script on OSX for data generation
# and on Ubuntu for plotting. For some reason, the pik1 processing is much
# faster on OSX. (And I failed to install all the proper python libraries
# for creating figures / running deva on OSX)
try:
    import matplotlib
    matplotlib.use('Qt4Agg')
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
except:
    pass

import numpy as np

import os
WAIS = os.getenv('WAIS')
import sys
sys.path.append('%s/syst/linux/src' % (WAIS))

try:
    import deva
    import deva.basemapUtilities
    import deva.devaUtilities
    import deva.radarUtilities
    import deva.utilities
except:
    pass

import pik1_utils


def plot_all_utig_data():
    '''
    Plots all UTIG transects on a simple basemap, with transects
    color-coded by season. 
    '''
    fig = Figure((24, 20))
    canvas = FigureCanvas(fig)
    ax = fig.add_axes([0,0,1,1])
    bgs = deva.basemapUtilities.make_background_dict(fig, ax)
    bg = bgs['modis_simple']
    bg.set_background()
    gls = deva.devaUtilities.make_grounding_line_dict()
    deva.devaUtilities.set_grounding_line(ax, gls, 'modis')

    ax.axis('equal')
    ax.set_xlim([-3000000, 3000000])
    ax.set_ylim([-2500000, 2500000])

    ax.tick_params(which='both', bottom=False, top=False, left=False, right=False)
    for side in ['bottom', 'top', 'left', 'right']:
        ax.spines[side].set_visible(False)

    transects = deva.devaUtilities.load_transects(antarctic=True)
    season_lookup = deva.utilities.SeasonLookup()
    for pst, data in transects.iteritems():
        season,_ = season_lookup.get_season(pst)
        if season is None:
            print "No season found for %s" % (pst)
            continue
        elif season in ['ASE1', '2001', 'ICP1', 'ICP2', 'ICP3', 'ICP4', 'ICP5']:
            color = 'k'
            zorder = 3
        elif season in ['ICP6']:
            color = 'darkgrey'
            zorder = 2
        else:
            color = 'lightgrey'
            zorder = 1
        ax.plot(data[:,1], data[:,2], color=color, linewidth=1.0, zorder=zorder)
    
    canvas.print_figure('../FinalReport/figures/all_data.png')


def plot_filter_freqs():
    '''
    Shows the various frequencies within the system, and what the filter does.
    '''
    # Each bin is one sample; sampled at 50Mhz
    nbins = 3200

    # Time of each sample, in us
    tt = np.arange(0, nbins/50., 1/50.)
    # 70MHz local oscillator signal
    y70 = np.sin(2*np.pi*tt*70)
    fft70 = np.fft.fft(y70, n=nbins)
    # 140MHz - LO's first harmonic
    y140 = np.sin(2*np.pi*tt*140)
    fft140 = np.fft.fft(y140, n=nbins)

    # Reference chirp used for actual pik1 dechirping
    ref_chirp = pik1_utils.generate_reference_chirp()
    y_ref = np.zeros(nbins)
    y_ref[100:100+len(ref_chirp)] = ref_chirp
    fft_ref = np.fft.fft(y_ref, n=nbins)

    # theoretical reference chirp
    theoretical_chirp = pik1_utils.generate_theoretical_chirp()
    y_theory = np.zeros(nbins)
    y_theory[100:100+len(theoretical_chirp)] = theoretical_chirp
    fft_theory = np.fft.fft(y_theory, n=nbins)    

    # Frequency of fft bins, in MHz
    fftfreq = np.fft.fftfreq(nbins, 1/50.)

    # ideal chirp frequency content
    fft_ideal = [1.0 if 2.5 <= np.abs(elem) <= 17.5 else 0.0 for elem in fftfreq]

    # Hanning filter
    hfilter = pik1_utils.generate_hfilter(nbins)

    fig = Figure((15,3.5))
    canvas = FigureCanvas(fig)
    ax = fig.add_axes([0.05, 0.25, 0.9, 0.7])

    ax.plot(fftfreq, np.abs(fft_ref)/np.max(np.abs(fft_ref)), 
            linewidth=1, color='black', label='pik1 reference chirp')
    ax.plot(fftfreq, np.abs(fft_theory)/np.max(np.abs(fft_theory)), 
            linewidth=2, color='darkgrey', label='theoretical reference chirp')
    ax.plot(fftfreq, np.abs(fft_ideal)/np.max(np.abs(fft_ideal)), 
            linewidth=2, color='lightgrey', label='theoretical reference chirp')
    ax.plot(fftfreq, np.abs(fft70)/np.max(np.abs(fft70)), 
            linewidth=2, color='blue', label='FFT of 70MHz signal')
    ax.plot(fftfreq, np.abs(fft140)/np.max(np.abs(fft140)), 
            linewidth=2, color='blue', label='FFT of 1st harmonic of LO')
    ax.plot(fftfreq, hfilter/np.max(hfilter), linewidth=2, color='red',
            label='Hamming filter')
    ax.set_ylim([0, 1.05])
    ax.set_xlim([-25, 25])
    ax.set_xlabel('Frequency (MHz)', fontsize=24)
    #ax.legend()
    ax.tick_params(which='both', bottom=True, top=False, left=False, right=False,
                   labelbottom=True, labeltop=False, labelleft=False, labelright=False,
                   labelsize=18)
    for side in ['top', 'left', 'right']:
        ax.spines[side].set_visible(False)

    canvas.print_figure('../FinalReport/figures/filter_frequencies.png')


def generate_data_products():
    '''
    I need a variety of intermediate data products in order to create 
    the figures.
    '''
    TOT_bounds, VCD_bounds, THW_bounds = find_quiet_regions()
    VCD_rawfile = WAIS + '/orig/xlob/VCD/JKB2g/DVD01a/RADnh3/bxds'
    TOT_rawfile = WAIS + '/orig/xlob/TOT/JKB2d/X16a/RADnh3/bxds'
    THW_rawfile = WAIS + '/orig/xlob/THW/SJB2/DRP02a/RADjh1/bxds1'
    
    insamples = 3437
    outsamples = 3200
    channel = 2
    blanking = 200
    scale = 1000

    num_sweeps = pik1_utils.get_num_traces(TOT_rawfile, insamples)
    reference_chirp = pik1_utils.generate_reference_chirp()


    
    # # First, I want raw data, just for grins. (For this, blanking should be 0)
    # # plotting range on this is 25000 to 90000
    # input_sweeps = pik1_utils.load_raw_nh3_subset_gen(
    #     TOT_rawfile, channel, np.arange(0, num_sweeps, 50), insamples, outsamples, 0)
    # pik1_utils.save_file(input_sweeps, 'TOT_raw', scale)

    # # Only dechirping, no stacking
    # # plotting range on this is 115000 to 200000
    # input_sweeps = pik1_utils.load_raw_nh3_subset_gen(
    #     TOT_rawfile, channel, np.arange(0, num_sweeps, 50), insamples, outsamples, blanking)
    # dechirped_traces = pik1_utils.dechirp_gen(input_sweeps, reference_chirp,
    #                                           outsamples, use_hamming=False)
    # pik1_utils.save_file(dechirped_traces, 'TOT_dechirped', scale)

    # # Dechirped and filtered
    # # It is important to do this first for better comparisons of w/ and w/o
    # # coherent stacking. It also has the nice property of suppressing the
    # # horizontal LO stripes until I actually want to talk about them.
    # input_sweeps = pik1_utils.load_raw_nh3_subset_gen(
    #     TOT_rawfile, channel, np.arange(0, num_sweeps, 50), insamples, outsamples, blanking)
    # dechirped_traces = pik1_utils.dechirp_gen(input_sweeps, reference_chirp,
    #                                           outsamples, use_hamming=True)
    # pik1_utils.save_file(dechirped_traces, 'TOT_filtered', scale)

    # # Dechirping and incoherent stacking
    # # For a fair comparison with coherent stacking (w/r/t horizontal resolution)
    # # the total depth needs to be the same.
    # input_sweeps = pik1_utils.load_raw_nh3_gen(TOT_rawfile, channel, insamples,
    #                                            outsamples, blanking)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 50)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True)
    # pik1_utils.save_file(dechirped_traces, 'TOT_costacked', scale)

    # # Dechirping, incoherent and coherent stacking
    # input_sweeps = pik1_utils.load_raw_nh3_gen(TOT_rawfile, channel, insamples,
    #                                            outsamples, blanking)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'TOT_stacked', scale)

    
    # trace_start, trace_end = 50*TOT_bounds[0], 50*TOT_bounds[1]
    # correction = pik1_utils.find_LO_params(TOT_rawfile, 2, 3437, trace_start, trace_end)
    # print "LO correction: ", correction
    # print "mag = %r, phs = %r" % (np.abs(correction), np.angle(correction))
    # LO_correction = correction*10 # correct for coherently stacking 10X
    # input_sweeps = pik1_utils.load_raw_nh3_gen(TOT_rawfile, channel, insamples,
    #                                        outsamples, 0)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True, 
    #                                           LO_correction=LO_correction)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'TOT_no_blanking', scale)

    # trace_start, trace_end = 50*VCD_bounds[0], 50*VCD_bounds[1]
    # correction = pik1_utils.find_LO_params(VCD_rawfile, 2, 3437, trace_start, trace_end)
    # print "LO correction: ", correction
    # print "mag = %r, phs = %r" % (np.abs(correction), np.angle(correction))
    # LO_correction = correction*10 # correct for coherently stacking 10X
    # input_sweeps = pik1_utils.load_raw_nh3_gen(VCD_rawfile, channel, insamples,
    #                                        outsamples, 0)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True, 
    #                                           LO_correction=LO_correction)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'VCD_no_blanking', scale)

    # # TODO: CHANNEL1 IS BORKED. TRY AGAIN W/ LO
    # input_sweeps = pik1_utils.load_raw_nh3_gen(TOT_rawfile, 1, insamples,
    #                                        outsamples, 0)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'TOT_ch1_no_blanking', scale)

    # input_sweeps = pik1_utils.load_raw_nh3_gen(VCD_rawfile, 1, insamples,
    #                                        outsamples, 0)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'VCD_ch1_no_blanking', scale)    
    
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'TOT_LO', scale)
    # trace_start, trace_end = 50*TOT_bounds[0], 50*TOT_bounds[1]
    # correction = pik1_utils.find_LO_params(TOT_rawfile, 2, 3437, trace_start, trace_end)
    # print "LO correction: ", correction
    # print "mag = %r, phs = %r" % (np.abs(correction), np.angle(correction))
    # LO_correction = correction*10 # correct for coherently stacking 10X
    # input_sweeps = pik1_utils.load_raw_nh3_gen(TOT_rawfile, channel, insamples,
    #                                        outsamples, blanking)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True, 
    #                                           LO_correction=LO_correction)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'TOT_LO', scale)

    # trace_start, trace_end = 50*VCD_bounds[0], 50*VCD_bounds[1]
    # correction = pik1_utils.find_LO_params(VCD_rawfile, 2, 3437, trace_start, trace_end)
    # print "LO correction: ", correction
    # print "mag = %r, phs = %r" % (np.abs(correction), np.angle(correction))
    # LO_correction = correction*10 # correct for coherently stacking 10X
    # input_sweeps = pik1_utils.load_raw_nh3_gen(VCD_rawfile, channel, insamples,
    #                                        outsamples, blanking)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 10)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True, 
    #                                           LO_correction=LO_correction)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # pik1_utils.save_file(inco_stacks, 'VCD_LO', scale)

    # trace_start, trace_end = 50*THW_bounds[0], 50*THW_bounds[1]
    # correction = pik1_utils.find_LO_params(THW_rawfile, None, 3200, trace_start, trace_end, nh3=False)
    # print "LO correction: ", correction
    # print "mag = %r, phs = %r" % (np.abs(correction), np.angle(correction))
    # LO_correction = correction*10 # correct for coherently stacking 10X
    # input_sweeps = pik1_utils.load_raw_jh1_gen(THW_rawfile, 3200, blanking=150)
    # coh_stacks = pik1_utils.coherent_stack_gen(input_sweeps, 50)
    # dechirped_traces = pik1_utils.dechirp_gen(coh_stacks, reference_chirp,
    #                                           outsamples, use_hamming=True,
    #                                           LO_correction=LO_correction)
    # inco_stacks = pik1_utils.inco_stack_gen(dechirped_traces, 5)
    # #pik1_utils.save_file(coh_stacks, 'THW_LO', scale)
    # pik1_utils.save_file(inco_stacks, 'THW_LO', scale)


def plot_radar(ax, filename, num_samples, bounds=None, clim=None):
    '''
    ax: axis to plot to
    filename: radargram to be plotted (in same file format as pik1)
    bounds: [min_trace, max_trace, min_sample, max_sample] of region to plot.
            If None, plot whole transect.
    clim: [min_counts, max_counts] for the data. If None, use 
          automatically-generated bounds.
    '''
    num_traces = pik1_utils.get_num_traces(filename, num_samples)
    radar_data = np.memmap(filename, '>i4', mode='r', shape=(num_traces, num_samples))
    # TODO: How important is the radar skip for generating images?
    if bounds is None:
        plot_data = radar_data
    else:
        min_trace, max_trace, min_sample, max_sample = bounds
        plot_data = radar_data[min_trace:max_trace, min_sample:max_sample]

    radar_plot = ax.imshow(plot_data.T, aspect='auto', interpolation='nearest')

    if bounds is not None:
        min_trace, max_trace, min_sample, max_sample = bounds
        radar_plot.set_extent([min_trace, max_trace, max_sample, min_sample])
    radar_plot.set_cmap('gray')
    if clim is not None:
        radar_plot.set_clim(clim)

def plot_bounding_box(ax, bounds, text, linewidth=8):
    '''
    * ax: axis to add box to
    * bounds: [xmin, xmax, ymin, ymax] of bounding box to draw
    * text: text to add to bottom left of box
    '''
    xmin, xmax, ymin, ymax = bounds
    ax.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 
            color='red', linewidth=linewidth)
    # The bounds are messed up b/c radar data has neg y axis.
    ax.text(xmin+50, ymax-50, text, fontsize=100, color='red')

def plot_data_products():

    # Define all bounds first s.t. we can draw the context boxes on the 
    # dechirped image.

    # The raw bounds are for the full image
    raw_min_sample = 150
    raw_max_sample = 2300
    raw_bounds = [0, 9075, raw_min_sample, raw_max_sample]

    # The filtering needs to zoom in on the surface
    filter_bounds = [1700, 5900, 425, 650]

    # For the SNR improvements of incoherent stacking, look at the layers that just pop out.
    layer_bounds = [7400, 8580, 200, 1750]
    # For filtered, clim = [110000, 195000] and for stacked, clim=[140000, 225000]

    # Zooming in on the crevasses shows speckle nicely
    incoherent_bounds = [2880, 4900, 850, 1700]
    #The clim for this is best is the coherent is [150000, 234000] and incoherent is [160000, 234000]

    # Appropriate color limits depend on the processing used
    raw_clim = [25000, 90000]
    dechirped_clim = [115000, 200000]
    filtered_clim = [115000, 200000]

    # First, generate the raw figure that requires raw data + dechirped data
    # over the full transect.
    raw_shape = (50, 17)
    #raw_shape = (10, 17./5) # This is ugly ... 
    
    # fig_raw = Figure(raw_shape, dpi=150)
    # canvas_raw = FigureCanvas(fig_raw)
    # ax_raw = fig_raw.add_axes([0, 0, 1, 1])
    # ax_raw.axis('off')
    # plot_radar(ax_raw, 'TOT_raw', 3200, raw_bounds, raw_clim)
    # canvas_raw.print_figure('../FinalReport/figures/TOT_raw_full.jpg')


    multiple_bounds = [3050, 5325, 325, 2675]
    fig_multiples = Figure(raw_shape, dpi=150)
    canvas_multiples = FigureCanvas(fig_multiples)
    ax_multiples = fig_multiples.add_axes([0, 0, 1, 1])
    ax_multiples.axis('off')
    plot_radar(ax_multiples, 'TOT_no_blanking', 3200, multiple_bounds, clim=[140000, 230000])
    ax_multiples.text(3950, 900, 'surface multiple', color='red', fontsize=70,
                      horizontalalignment='left', verticalalignment='bottom')
    ax_multiples.text(3950, 2380, 'basal multiple', color='red', fontsize=70,
            horizontalalignment='left', verticalalignment='top')
    canvas_multiples.print_figure('../FinalReport/figures/TOT_multiples.jpg')

    # fig_labeled = Figure(raw_shape, dpi=150)
    # canvas_labeled = FigureCanvas(fig_labeled)
    # ax_labeled = fig_labeled.add_axes([0, 0, 1, 1])
    # ax_labeled.axis('off')
    # plot_radar(ax_labeled, 'TOT_LO', 3200, raw_bounds, clim=[136000, 233000])
    # xlim = ax_labeled.get_xlim()
    # ylim = ax_labeled.get_ylim()
    # plot_bounding_box(ax_labeled, filter_bounds, '2', linewidth=8)
    # plot_bounding_box(ax_labeled, layer_bounds, '3', linewidth=8)
    # plot_bounding_box(ax_labeled, incoherent_bounds, '4', linewidth=8)
    # # This would require the image to extend further in samples, which in turn
    # # would require a different aspect ratio. 
    # #plot_bounding_box(ax_labeled, multiple_bounds, '5', linewidth=8)
    # dc_bounds = [6000, 9075, 200, 800]    
    # plot_bounding_box(ax_labeled, dc_bounds, '1', linewidth=8)
    # trace_idx = 3900
    # ax_labeled.plot([trace_idx, trace_idx], raw_bounds[2:], 'r--', linewidth=8)
    # ax_labeled.set_xlim(xlim)
    # ax_labeled.set_ylim(ylim)
    # canvas_labeled.print_figure('../FinalReport/figures/TOT_pik1_labeled.jpg')

    # fig_dechirped = Figure(raw_shape, dpi=150)
    # canvas_dechirped = FigureCanvas(fig_dechirped)
    # ax_dechirped = fig_dechirped.add_axes([0, 0, 1, 1])
    # ax_dechirped.axis('off')
    # plot_radar(ax_dechirped, 'TOT_dechirped', 3200, raw_bounds, dechirped_clim)
    # # xlim = ax_dechirped.get_xlim()
    # # ylim = ax_dechirped.get_ylim()
    # # plot_bounding_box(ax_dechirped, filter_bounds, '1', linewidth=8)
    # # plot_bounding_box(ax_dechirped, layer_bounds, '2', linewidth=8)
    # # plot_bounding_box(ax_dechirped, incoherent_bounds, '3', linewidth=8)
    # # dc_bounds = [6000, 9075, 200, 800]    
    # # plot_bounding_box(ax_dechirped, dc_bounds, '4', linewidth=8)
    # # trace_idx = 3900
    # # ax_dechirped.plot([trace_idx, trace_idx], raw_bounds[2:], 'r--', linewidth=8)
    # # ax_dechirped.set_xlim(xlim)
    # # ax_dechirped.set_ylim(ylim)
    # canvas_dechirped.print_figure('../FinalReport/figures/TOT_dechirped_full.jpg')

    # # Next, generate filtered figure - zoom in on the surface
    # filter_shape = (50, 10)
    # fig_dechirped_zoom1 = Figure(filter_shape)
    # canvas_dechirped_zoom1 = FigureCanvas(fig_dechirped_zoom1)
    # ax_dechirped_zoom1 = fig_dechirped_zoom1.add_axes([0, 0, 1, 1])
    # ax_dechirped_zoom1.axis('off')
    # plot_radar(ax_dechirped_zoom1, 'TOT_dechirped', 3200, filter_bounds, dechirped_clim)
    # canvas_dechirped_zoom1.print_figure('../FinalReport/figures/TOT_dechirped_zoom1.jpg')

    # fig_filtered_zoom1 = Figure(filter_shape)
    # canvas_filtered_zoom1 = FigureCanvas(fig_filtered_zoom1)
    # ax_filtered_zoom1 = fig_filtered_zoom1.add_axes([0, 0, 1, 1])
    # ax_filtered_zoom1.axis('off')
    # plot_radar(ax_filtered_zoom1, 'TOT_filtered', 3200, filter_bounds, filtered_clim)
    # canvas_filtered_zoom1.print_figure('../FinalReport/figures/TOT_filtered_zoom1.jpg')


    # # Now to show the horizontal and vertical resolution improvements for coherent stacking
    # # The layers will be good for SNR; what will be good for horizontal resolution?
    # layer_shape = (25, 25)
    # fig_filtered_zoom2 = Figure(layer_shape)
    # canvas_filtered_zoom2 = FigureCanvas(fig_filtered_zoom2)
    # ax_filtered_zoom2 = fig_filtered_zoom2.add_axes([0, 0, 1, 1])
    # ax_filtered_zoom2.axis('off')
    # plot_radar(ax_filtered_zoom2, 'TOT_filtered', 3200, layer_bounds, [110000, 195000])
    # canvas_filtered_zoom2.print_figure('../FinalReport/figures/TOT_filtered_zoom2.jpg')

    # fig_costacked_zoom2 = Figure(layer_shape)
    # canvas_costacked_zoom2 = FigureCanvas(fig_costacked_zoom2)
    # ax_costacked_zoom2 = fig_costacked_zoom2.add_axes([0, 0, 1, 1])
    # ax_costacked_zoom2.axis('off')
    # plot_radar(ax_costacked_zoom2, 'TOT_costacked', 3200, layer_bounds, [140000, 225000])
    # canvas_costacked_zoom2.print_figure('../FinalReport/figures/TOT_costacked_zoom2.jpg')

    # # Finally, show the speckle improvements from incoherent stacking
    # speckle_shape = (25, 25)
    # fig_costacked_zoom3 = Figure(speckle_shape)
    # canvas_costacked_zoom3 = FigureCanvas(fig_costacked_zoom3)
    # ax_costacked_zoom3 = fig_costacked_zoom3.add_axes([0, 0, 1, 1])
    # ax_costacked_zoom3.axis('off')
    # plot_radar(ax_costacked_zoom3, 'TOT_costacked', 3200, incoherent_bounds, [150000, 234000])
    # canvas_costacked_zoom3.print_figure('../FinalReport/figures/TOT_costacked_zoom3.jpg')

    # fig_stacked_zoom3 = Figure(speckle_shape)
    # canvas_stacked_zoom3 = FigureCanvas(fig_stacked_zoom3)
    # ax_stacked_zoom3 = fig_stacked_zoom3.add_axes([0, 0, 1, 1])
    # ax_stacked_zoom3.axis('off')
    # plot_radar(ax_stacked_zoom3, 'TOT_stacked', 3200, incoherent_bounds, [160000, 234000])
    # canvas_stacked_zoom3.print_figure('../FinalReport/figures/TOT_stacked_zoom3.jpg')

def plot_before_after():
    
    shape = (36, 27)
    fig_old = Figure(shape, dpi=150)
    canvas_old = FigureCanvas(fig_old)
    ax_old = fig_old.add_axes([0, 0, 1, 1])
    ax_old.axis('off')
    old_filename = WAIS + '/targ/xtra/ICP3/CMP/pik1.RADnh3/TOT/JKB2d/X16a/MagLoResInco2'
    plot_radar(ax_old, old_filename, 3200, clim=[103000, 200000])
    canvas_old.print_figure('../FinalReport/figures/TOT_X16a_old.jpg')

    fig_new = Figure(shape, dpi=150)
    canvas_new = FigureCanvas(fig_new)
    ax_new = fig_new.add_axes([0, 0, 1, 1])
    ax_new.axis('off')
    plot_radar(ax_new, 'TOT_LO', 3200, clim=[136000, 233000])
    canvas_new.print_figure('../FinalReport/figures/TOT_X16a_new.jpg')


def plot_DC():    

    shape = (36, 18)

    TOT_bounds = [6000, 9075, 200, 800]
    
    fig_TOT_ch1 = Figure(shape, dpi=150)
    canvas_TOT_ch1 = FigureCanvas(fig_TOT_ch1)
    ax_TOT_ch1 = fig_TOT_ch1.add_axes([0, 0, 1, 1])
    ax_TOT_ch1.axis('off')
    old_filename = WAIS + '/targ/xtra/ICP3/CMP/pik1.RADnh3/TOT/JKB2d/X16a/MagLoResInco1'
    plot_radar(ax_TOT_ch1, old_filename, 3200, bounds=TOT_bounds, clim=[80000, 180000])
    canvas_TOT_ch1.print_figure('../FinalReport/figures/TOT_ch1_DC.jpg')

    fig_TOT_ch2 = Figure(shape, dpi=150)
    canvas_TOT_ch2 = FigureCanvas(fig_TOT_ch2)
    ax_TOT_ch2 = fig_TOT_ch2.add_axes([0, 0, 1, 1])
    ax_TOT_ch2.axis('off')
    plot_radar(ax_TOT_ch2, 'TOT_no_blanking', 3200, bounds=TOT_bounds, clim=[130000, 230000])
    canvas_TOT_ch2.print_figure('../FinalReport/figures/TOT_ch2_DC.jpg')

    VCD_bounds = [9600, 10470, 200, 800]
    
    fig_VCD_ch1 = Figure(shape, dpi=150)
    canvas_VCD_ch1 = FigureCanvas(fig_VCD_ch1)
    ax_VCD_ch1 = fig_VCD_ch1.add_axes([0, 0, 1, 1])
    ax_VCD_ch1.axis('off')
    old_filename = WAIS + '/targ/xtra/ICP4/CMP/pik1.RADnh3/VCD/JKB2g/DVD01a/MagLoResInco1'
    plot_radar(ax_VCD_ch1, old_filename, 3200, bounds=VCD_bounds, clim=[80000, 180000])
    canvas_VCD_ch1.print_figure('../FinalReport/figures/VCD_ch1_DC.jpg')

    fig_VCD_ch2 = Figure(shape, dpi=150)
    canvas_VCD_ch2 = FigureCanvas(fig_VCD_ch2)
    ax_VCD_ch2 = fig_VCD_ch2.add_axes([0, 0, 1, 1])
    ax_VCD_ch2.axis('off')
    plot_radar(ax_VCD_ch2, 'VCD_no_blanking', 3200, bounds=VCD_bounds, clim=[140000, 230000])
    canvas_VCD_ch2.print_figure('../FinalReport/figures/VCD_ch2_DC.jpg')    


    
def plot_LO_horiz_stripes():
    '''
    This uses data that has been processed through pik1 but w/ the hanning filter
    disabled s.t. the stripes are more readily apparent.
    '''
    fig = Figure((10, 4))
    canvas = FigureCanvas(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    plot_radar(ax, 'TOT_stacked_nofilter', 3200, None, [135000, 230000])
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    TOT_bounds, VCD_bounds, THW_bounds = find_quiet_regions()
    ax.vlines(TOT_bounds[0:2], 0, 3200, colors='red', linewidth=3, linestyles='dashed')    
    plot_bounding_box(ax, TOT_bounds, '', linewidth=4)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    canvas.print_figure('../FinalReport/figures/TOT_LO_stripes_d.jpg')

    zoom_bounds = [5000, 9000, 3000, 3200]
    zoom_fig = Figure((2.5, 9))
    zoom_canvas = FigureCanvas(zoom_fig)
    zoom_ax = zoom_fig.add_axes([0, 0, 1, 1])
    zoom_ax.axis('off')
    plot_radar(zoom_ax, 'TOT_stacked_nofilter', 3200, zoom_bounds, [135000, 230000])
    zoom_canvas.print_figure('../FinalReport/figures/TOT_LO_stripes_zoom.jpg')



def plot_quiet_regions():
    # Plot the region that the noise was calculated from...
    # For TOT...
    TOT_bounds, VCD_bounds, THW_bounds = find_quiet_regions()

    # TOT/JKB2d/X16a gives: 
    # mag = 32809.224658469, phs = -0.90421798501485484
    # VCD/JKB2g/DVD01a gives:
    # mag = 15720.217174332585, phs = -0.98350090576267946
    # THW/SJB2/DRP02a gives:
    # 26158.900202734963, phs = 1.6808311318828895
    
    TOT_fig = Figure((10, 8))
    TOT_canvas = FigureCanvas(TOT_fig)
    TOT_ax = TOT_fig.add_axes([0, 0, 1, 1])
    TOT_ax.axis('off')
    plot_radar(TOT_ax, 'TOT_LO', 3200, None, [135000, 234000])
    xlim = TOT_ax.get_xlim()
    ylim = TOT_ax.get_ylim()
    TOT_ax.vlines(TOT_bounds[0:2], 0, 3200, colors='red', linewidth=3, linestyles='dashed')
    plot_bounding_box(TOT_ax, TOT_bounds, '', linewidth=4)
    TOT_ax.set_xlim(xlim)
    TOT_ax.set_ylim(ylim)
    TOT_canvas.print_figure('../FinalReport/figures/TOT_quiet_region.jpg')

    VCD_fig = Figure((10, 8))
    VCD_canvas = FigureCanvas(VCD_fig)
    VCD_ax = VCD_fig.add_axes([0, 0, 1, 1])
    VCD_ax.axis('off')
    plot_radar(VCD_ax, 'VCD_LO', 3200, None, [135000, 234000])
    xlim = VCD_ax.get_xlim()
    ylim = VCD_ax.get_ylim()
    VCD_ax.vlines(VCD_bounds[0:2], 0, 3200, colors='red', linewidth=3, linestyles='dashed')
    plot_bounding_box(VCD_ax, VCD_bounds, '', linewidth=4)
    VCD_ax.set_xlim(xlim)
    VCD_ax.set_ylim(ylim)
    VCD_canvas.print_figure('../FinalReport/figures/VCD_quiet_region.jpg')

    THW_fig = Figure((10, 8))
    THW_canvas = FigureCanvas(THW_fig)
    THW_ax = THW_fig.add_axes([0, 0, 1, 1])
    THW_ax.axis('off')
    plot_radar(THW_ax, 'THW_LO', 3200, None, [135000, 234000])
    xlim = THW_ax.get_xlim()
    ylim = THW_ax.get_ylim()
    THW_ax.vlines(THW_bounds[0:2], 0, 3200, colors='red', linewidth=3, linestyles='dashed')
    plot_bounding_box(THW_ax, THW_bounds, '', linewidth=4)
    THW_ax.set_xlim(xlim)
    THW_ax.set_ylim(ylim)
    THW_canvas.print_figure('../FinalReport/figures/THW_quiet_region.jpg')


    
def plot_LO_fft():
    # pik1_utils does all of this ...
    filename = WAIS + '/orig/xlob/TOT/JKB2d/X16a/RADnh3/bxds'
    
    trace_start, trace_end = 405450, 415450 # This is 8109 - 8309 in pik1 traces
    correction = pik1_utils.find_LO_params(filename, 2, 3437, trace_start, trace_end)
    print "LO correction: ", correction
    print "mag = %r, phs = %r" % (np.abs(correction), np.angle(correction))


    # Plot the fft of the full traces
    input_sweeps = pik1_utils.load_raw_nh3_subset_gen(
        filename, 2, np.arange(trace_start, trace_end), 3437, 3200, 200)
    input_stacked = pik1_utils.coherent_stack_gen(input_sweeps, 50)
    radar_data = np.array([trace for trace in input_stacked])
    radar_fft = np.fft.fft(radar_data, n=3200)
    fftfreq3200 = np.fft.fftfreq(3200, 1/50.)
    fig2 = Figure((10, 8))
    canvas2 = FigureCanvas(fig2)
    ax2 = fig2.add_axes([0, 0, 1, 1])
    ax2.imshow(np.abs(radar_fft), cmap='gray', aspect='auto')
    canvas2.print_figure('../FinalReport/figures/data_fft.jpg')

    # Plot the fft of the bottom portion only
    input_sweeps = pik1_utils.load_raw_nh3_subset_gen(
        filename, 2, np.arange(trace_start, trace_end), 3437, 3437, 200)
    input_stacked = pik1_utils.coherent_stack_gen(input_sweeps, 50)
    noise_data = np.array([trace[-800:] for trace in input_stacked])
    noise_fft = np.fft.fft(noise_data, n=800)
    fftfreq800 = np.fft.fftfreq(800, 1/50.)
    fig3 = Figure((10, 8))
    canvas3 = FigureCanvas(fig3)
    ax3 = fig3.add_axes([0, 0, 1, 1])
    ax3.imshow(np.abs(noise_fft), cmap='gray', aspect='auto')
    canvas3.print_figure('../FinalReport/figures/noise_fft.jpg')


    # TODO: Maybe just to nice clean line plots showing full transect and bottom?
    fig4 = Figure((15, 5))
    canvas4 = FigureCanvas(fig4)
    ax4 = fig4.add_axes([0.01, 0.2, 0.98, 0.8])
    abs_noise_fft = np.sum(np.abs(noise_fft), axis=0)
    # cancel out constant term.
    abs_noise_fft[0] = 0 
    plot_noise_fft = abs_noise_fft/np.max(abs_noise_fft)
    ax4.plot(fftfreq800[1:400], plot_noise_fft[1:400], color='black')
    ax4.plot(fftfreq800[401:800], plot_noise_fft[401:800], color='black')
    # Cancel out the two biggest peaks.
    abs_noise_fft[160] = 0
    abs_noise_fft[640] = 0
    plot_noise_fft2 = abs_noise_fft/np.max(abs_noise_fft)
    ax4.plot(fftfreq800[1:400], plot_noise_fft2[1:400], color='darkgrey')
    ax4.plot(fftfreq800[401:800], plot_noise_fft2[401:800], color='darkgrey')
    abs_radar_fft = np.sum(np.abs(radar_fft), axis=0)
    plot_radar_fft = abs_radar_fft/np.max(abs_radar_fft)
    # Plot in two chunks to avoid the awkward line across the figure.
    ax4.plot(fftfreq3200[1:1600], plot_radar_fft[1:1600], color='red')
    ax4.plot(fftfreq3200[1601:3200], plot_radar_fft[1601:3200], color='red')

    ax.set_ylim([0, 1.05])
    ax4.set_xlim([-25, 25])
    ax4.tick_params(which='both', bottom=True, top=False, left=False, right=False,
                    labelbottom=True, labeltop=False, labelleft=False, labelright=False,
                    labelsize=18)
    for side in ['top', 'left', 'right']:
        ax4.spines[side].set_visible(False)
    ax4.set_xlabel('Frequency (MHz)', fontsize=24)

    canvas4.print_figure('../FinalReport/figures/LO_fft.jpg')


def find_quiet_regions():
    '''
    Automaticaly find quiet regions. By inspection, all are good choices =)
    Returns the coordinates for the pik1 product. 
    '''
    
    TOT_filename = WAIS + '/orig/xlob/TOT/JKB2d/X16a/RADnh3/bxds'
    #TOT_trace_start, TOT_trace_end = pik1_utils.find_quiet_region(TOT_filename, 2, 3437, 3200)
    TOT_trace_start, TOT_trace_end = 405450, 415450 # This is 8109 - 8309 in pik1 traces
    TOT_bounds = (TOT_trace_start/50, TOT_trace_end/50, 2400, 3200)

    VCD_filename = WAIS + '/orig/xlob/VCD/JKB2g/DVD01a/RADnh3/bxds'
    #VCD_trace_start, VCD_trace_end = pik1_utils.find_quiet_region(VCD_filename, 2, 3437, 3200)
    VCD_trace_start, VCD_trace_end = 464600, 474600 # (9292 - 9492 in pik1 traces)
    VCD_bounds = (VCD_trace_start/50, VCD_trace_end/50, 2400, 3200)
        
    # TODO: whoops - I only have one of the two channels here :-\
    THW_filename = WAIS + '/orig/xlob/THW/SJB2/DRP02a/RADjh1/bxds1'
    #THW_trace_start, THW_trace_end = pik1_utils.find_quiet_region(THW_filename, 2, 3200, 3200)
    THW_trace_start, THW_trace_end = 42500, 52500 #(850 - 1050 in pik1 traces)
    THW_bounds = (THW_trace_start/50, THW_trace_end/50, 2400, 3200)
    
    return TOT_bounds, VCD_bounds, THW_bounds


def plot_trace():
    # Trace 3900 in pik1 looks pretty nice.
    trace_idx = 3900
    num_in_samples = 3437
    num_out_samples = 3200
    
    raw_y = 2 #(where to plot the raw data ...)
    raw_filename = WAIS + '/orig/xlob/TOT/JKB2d/X16a/RADnh3/bxds'
    # For plotting, I tried raw dB of trace ... I don't like it as much
    # If I try it again, gotta be careful - abs(elem) can cause overflow errors if we don't cast to float first!
    sweep_gen = pik1_utils.load_raw_nh3_subset_gen(
        raw_filename, 2, [50*trace_idx], num_in_samples, num_out_samples, 0)
    raw_trace = [elem for elem in sweep_gen][0]
    raw_min = np.min(raw_trace)
    raw_range = 1.0*np.max(raw_trace) - np.min(raw_trace)
    norm_raw_trace_ch2 = [raw_y + (1.0*elem - raw_min) / raw_range for elem in raw_trace]

    sweep_gen = pik1_utils.load_raw_nh3_subset_gen(
        raw_filename, 1, [50*trace_idx], num_in_samples, num_out_samples, 0)
    raw_trace = [elem for elem in sweep_gen][0]
    raw_min = np.min(raw_trace)
    raw_range = 1.0*np.max(raw_trace) - np.min(raw_trace)
    norm_raw_trace_ch1 = [raw_y + 1 + (1.0*elem - raw_min) / raw_range for elem in raw_trace]


    pik1_filename = 'TOT_LO'
    num_traces = pik1_utils.get_num_traces(pik1_filename, num_in_samples)
    pik1_data = np.memmap(pik1_filename, '>i4', mode='r', shape=(num_traces, num_out_samples))
    processed_trace = pik1_data[trace_idx,:]
    trace_min = np.min(processed_trace)
    trace_range = np.max(processed_trace) - np.min(processed_trace)
    norm_processed_trace = [1.0*(elem - trace_min) / trace_range for elem in processed_trace]

    fig = Figure((30, 15))
    canvas = FigureCanvas(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.minorticks_off()
    ax.tick_params(which='both', 
                    bottom=False, top=False, left=False, right=False,
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for side in ['bottom', 'top', 'left', 'right']:
        ax.spines[side].set_visible(False)

    ax.set_xlim([0, 3200])
    ylim = [0, raw_y + 2]
    ax.set_ylim(ylim)

    bar_y1 = 1.15
    bar_y2 = 1.25
    bar_y3 = 1.35
    text_y1 = 1.0
    text_y2 = 1.3
    text_y3 = 1.4

    fontsize=50
    xshift = 66 # WTF. Why do I have to offset the raw pulse to make it line up with the dechirped?
    
    # reference chirp 50-100
    ax.plot([50, 100], [1.65, 1.65], color='red', linewidth=15)
    ax.text(25, 1.7, 'chirp', fontsize=50, color='red',
            horizontalalignment='left', verticalalignment='bottom')

    # "main bang" 
    ax.plot([100, 540], [1.35, 1.35], '--', color='red', linewidth=15)
    ax.plot([100, 400], [1.35, 1.35], color='red', linewidth=15)
    ax.text(250, 1.4, 'main bang', fontsize=50, color='red',
            horizontalalignment='center', verticalalignment='bottom')
    
    # hardware blanking from 1-134 (in theory ... I don't trust this from the resulting plots
    ax.plot([xshift+0, xshift+134], [1.15, 1.15], color='red', linewidth=15)
    ax.text(xshift, 0.9, 'HW blanking', fontsize=50, color='red',
            horizontalalignment='left', verticalalignment='bottom')
    
    # software blanking from 1-200
    ax.plot([0, 200], [0.85, 0.85], color='red', linewidth=15)
    ax.text(25, 0.6, 'SW blanking', fontsize=50, color='red',
            horizontalalignment='left', verticalalignment='bottom')

    # Surface at 555
    ax.plot([555, 555], ylim, color='lightgray', linewidth=15)
    ax.text(575, 1.3, 'ice surface', fontsize=50, color='black',
            horizontalalignment='left', verticalalignment='bottom')

    # surface multiple at 950
    ax.plot([950, 950], ylim, color='lightgray', linewidth=15)
    ax.text(970, 0.9, 'surface \nmultiple', fontsize=50, color='black',
            horizontalalignment='left', verticalalignment='bottom')
    
    # Crevasses at 1262
    ax.plot([1262, 1262], ylim, color='lightgray', linewidth=15)
    ax.text(1282, 1.7, 'crevasse', fontsize=50, color='black',
            horizontalalignment='left', verticalalignment='bottom')

    # water at 1378
    ax.plot([1378, 1378], ylim, color='lightgray', linewidth=15)
    ax.text(1398, 0.85, 'water', fontsize=50, color='black',
            horizontalalignment='left', verticalalignment='bottom')

    # off-nadir energy at 1400 - 1850
    ax.plot([1400, 2100], [bar_y2, bar_y2], '--', color='red', linewidth=15)
    ax.plot([1500, 1850], [bar_y2, bar_y2], color='red', linewidth=15)
    ax.text(1750, text_y2, 'off-nadir energy', fontsize=50, color='red',
            horizontalalignment='center', verticalalignment='bottom')

    # bed multiple at 2204
    ax.plot([2204, 2204], ylim, color='lightgray', linewidth=15)
    ax.text(2224, 0.5, 'bed multiple', fontsize=50, color='black',
            horizontalalignment='left', verticalalignment='bottom')
    
    # noise at 2400 - 3200 (Can I call this receiver noise?)
    ax.plot([2250, 3200], [bar_y2, bar_y2], '--', color='red', linewidth=15)
    ax.plot([2400, 3200], [bar_y2, bar_y2], color='red', linewidth=15)
    ax.text(2725, text_y2, 'receiver noise', fontsize=50, color='red',
            horizontalalignment='center', verticalalignment='bottom')

    
    # Ideas:
    # * have 3 different y values, for things that apply to one, the other, or both?
    # * shade the background to make it clearer what's going on?
    #   or at least draw light vertical grey bars for the point events?
    # * Be sure to show where this trace comes from on a transect ... 
    
    # ax1 = fig.add_axes([0.0, 0.05, 0.5, 0.9])
    # ax1.set_ylim([3200, 0])
    # ax1.minorticks_off()
    # ax1.tick_params(which='both', 
    #                 bottom=False, top=False, left=False, right=False,
    #                 labelbottom=False, labeltop=False, labelleft=False, labelright=False)
                    
    # for side in ['bottom', 'top', 'left', 'right']:
    #     ax1.spines[side].set_visible(False)    
    
    # ax2 = fig.add_axes([0.5, 0.05, 0.5, 0.9])
    # ax2.plot(processed_trace, range(len(processed_trace)), 'k', markersize=1)

    # ax2.set_ylim([3200, 0])
    # ax2.minorticks_on()
    # ax2.tick_params(which='both', direction='inout', labelsize=18, 
    #                 bottom=False, top=False, left=True, right=False,
    #                 labelbottom=False, labeltop=False, labelleft=True, labelright=False)
                    
    # for side in ['bottom', 'top', 'right']:
    #     ax2.spines[side].set_visible(False)

    ax.plot(range(len(norm_processed_trace)), norm_processed_trace, 'k.', markersize=6)
    #ax.plot(range(len(norm_processed_trace)), norm_processed_trace, color='lightgrey', linewidth=1)
    ax.plot(np.arange(xshift, xshift+len(norm_raw_trace_ch2)), norm_raw_trace_ch2, 'k.', markersize=6)
    ax.plot(np.arange(xshift, xshift+len(norm_raw_trace_ch1)), norm_raw_trace_ch1, 'k.', markersize=6)

    canvas.print_figure('../FinalReport/figures/trace.jpg')    
    
if __name__=="__main__":
    # Plots every line flown by UTIG on an outline of Antarctica.
    #plot_all_utig_data()

    # Generates the figure showing the frequency spectrum of our system. 
    #plot_filter_freqs()

    # Generates the figure showing the raw and pik1 trace with events labeled.
    #plot_trace() 

    # Generates figure showing the FFT of full traces and the selected quiet region.
    #plot_LO_fft()

    # Plots the reference chirp used in deconvolution on top of theoretical chirp.
    #plot_reference_chirp()

    # Where all of the various intermediate data products are generated.
    # Any of them can be opened using radarFigure's --filename option. 
    #generate_data_products()
    
    plot_data_products()

    # Shows the surface artifact in ch2 but not ch1
    #plot_DC()
    
    #plot_LO_horiz_stripes()

    #plot_quiet_regions()
    
    #plot_before_after()



    
