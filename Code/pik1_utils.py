import argparse
import numpy as np
import os

# TODO: Add params for removing the LO noise ...

# TODO: I'd really like to make the code auto-detect and then load the data ...
def get_num_traces(filename, num_samples, nh3=True):
    filesize = os.stat(filename).st_size
    if nh3:
        # RADnh3
        num_traces = filesize/(4*num_samples+8)
    else:
        # RADjh1
        num_traces = filesize/(2*num_samples)
    return num_traces


def load_raw_nh3_subset_gen(filename, channel, traces, num_in_samples=3437, num_out_samples=3200, blanking=200):
    '''
    Load in a subset of the traces from the input file. This is very useful for
    efficiently debugging and generating figures.

    Relies on the fact that HiCARS2 data comes in order, so we can index directly into the file.
    * filename - full path to bxds file
    * channel - 1 or 2
    * traces - zero-index list of traces to pull out of the file.
    * num_in_samples - number of samples recorded. 
    * num_out_samples - how many samples to output. 
    * blanking - zero out this many samples.
    '''
    num_traces = get_num_traces(filename, num_in_samples)
    print "%s has %d traces" % (filename, num_traces)
    with open(filename, 'r') as fp:
        for trace in traces:
            # seek relative to beginning of file. If this is too slow, 
            # try incrementing from previously sought position.
            offset = 8 + trace * (8 + 4*num_in_samples) + (channel-1)*2*num_in_samples
            fp.seek(offset, os.SEEK_SET)
            data = np.fromfile(fp, dtype='>i2', count=num_in_samples)
            data[0:blanking] = 0            
            yield data[0:num_out_samples]

            
def load_raw_nh3_gen(filename, channel, num_in_samples=3437, num_out_samples=3200, blanking=200):
    '''
    Relies on the fact that HiCARS2 data comes in order, so we can index directly into the file.
    * filename - full path to bxds file
    * channel - 1 or 2
    * traces - zero-index list of traces to pull out of the file.
    * num_in_samples - number of samples recorded. 
    * num_out_samples - how many samples to output. 
    * blanking - zero out this many samples.
    '''
    num_traces = get_num_traces(filename, num_in_samples)
    print "file has %d traces" % (num_traces)
    with open(filename, 'r') as fp:
        # while True doesn't work is OSX ... it just keeps going.
        for idx in range(num_traces):
            # Skip the header
            fp.seek(8, os.SEEK_CUR)
            
            if channel == 2:
                fp.seek(2*num_in_samples, os.SEEK_CUR)
            data = np.fromfile(fp, dtype='>i2', count=num_in_samples)
            data[0:blanking] = 0            
            if channel == 1:
                fp.seek(2*num_in_samples, os.SEEK_CUR)
            yield data[0:num_out_samples]
        print fp.tell()


def load_raw_jh1_subset_gen(filename, traces, num_samples=3200, blanking=150):
    '''
    Load in a subset of the traces from the input file. This is very useful for
    efficiently debugging and generating figures.

    Relies on the fact that HiCARS2 data comes in order, so we can index directly into the file.
    * filename - full path to bxds file
    * traces - zero-index list of traces to pull out of the file.
    * num_samples - number of samples recorded. 
    * blanking - zero out this many samples.
    '''
    num_traces = get_num_traces(filename, num_samples, nh3=False)
    print "%s has %d traces" % (filename, num_traces)
    with open(filename, 'r') as fp:
        for trace in traces:
            # seek relative to beginning of file. If this is too slow, 
            # try incrementing from previously sought position.
            offset = trace * 2 * num_samples
            fp.seek(offset, os.SEEK_SET)
            data = np.fromfile(fp, dtype='<i2', count=num_samples)
            data[0:blanking] = 0            
            yield data

        
def load_raw_jh1_gen(filename, num_samples=3200, blanking=50):
    '''
    Relies on the fact that ATRS/HiCARS data comes in order, so we can index directly into the file.
    Note that this has the opposite endianness from the newer data!!!
    * filename - full path to bxds file; no header info, only one channel per file
    * num_samples - number of samples recorded. 
    * blanking - zero out this many samples.
    '''
    num_traces = get_num_traces(filename, num_samples, nh3=False)
    print "file has %d traces" % (num_traces)
    with open(filename, 'r') as fp:
        for idx in range(num_traces):
            data = np.fromfile(fp, dtype='<i2', count=num_samples)
            data[0:blanking] = 0.0
            yield data


def load_raw_gen(filename, channel=2, num_in_samples=3437, num_out_samples=3200, blanking=200, nh3=True):
    '''
    # TODO: I'd really like to make the code auto-detect and then load the data ...
    '''
    if nh3:
        return load_raw_nh3_gen(filename, channel, num_in_samples, num_out_samples, blanking)
    else:
        return load_raw_jh1_gen(filename, num_out_samples, blanking)
        
def coherent_stack_gen(trace_gen, depth):
    '''
    Returns an image where consecutive traces have been added
    '''
    stack = []
    for trace in trace_gen:
        if len(stack) < depth:
            stack.append(trace) # I wonder how much slower this is than just creating the array one, which is totally trivial.
        else:  # only want to run this if we have a full stack
            yield np.sum(stack, 0)
            stack = []
            stack.append(trace)


def inco_stack_gen(trace_gen, depth):
    '''
    Incoherently stack incoming data. Assumes that data is complex.
    '''
    stack = []
    for trace in trace_gen:
        if len(stack) < depth:
            stack.append(np.abs(trace))
        else:
            yield np.sum(stack, 0)
            stack = []
            stack.append(trace)

            
# The reference chirp used for pulse compression. It is sampled at 50MHz.
# We transmit a linear chirp centered at 60MHz with 15MHz bandwidth.
# TODO: Say something about _why_ it gives better results?
def generate_reference_chirp():
    '''
    Simply returns the reference chirp that I found in the original pyk1 code.
    TODO: Figure out why this gives better results than the theoretical chirp
    or the reference chirp that is routed from the amplifier?
    '''
    reference_chirp = np.array([-63, -92, -109, -75, -87, -50, -116, -154, -22, 68, 141, -610, 1461, 3807, -6147, -5375, 10651, -4412, -9810, 15386, -3070, -14499, 15130, 3677, -15935, 3743, 13362, -5884, -13301, 8455, 12542, -8744, -11977, 5105, 13754, -961, -14342, -5184, 10294, 12194, -2709, -14352, -8807, 5965, 15350, 8368, -6605, -14990, -11515, 196, 11276, 15490, 10300, -645, -10730, -15307, -13379, -7342, 377, 7264, 11662, 13435, 11530, 3243, -4865, -4427, -3233, -4000, -2472, -2498, -2361, -1230, -1311, -618, -578, -569, -121, -319, 206, 328, 436, 613, 318, 514, 353, 277, 221, 34, 250, 132, 199, 189, 75, 190, 65, 106, 19, -64, -14, -117])
    return reference_chirp


def generate_theoretical_chirp():
    '''
    Generates a 1us long chirp with center frequency 60MHz and 15MHz bandwidth,
    sampled at 50MHz. This matches the theoretically-transmitted one subsampled
    at actual digitizer frequency.
    (The generated chirp is sampled at 200MHz)
    '''
    # Initial phase
    p0 = np.pi/2
    # Starting frequency (MHz)
    f0 = 52.5
    # Bandwidth (MHz)
    bw = 15
    # Sample times (us)
    times = np.arange(0, 1, 1.0/50)
    chirp = [np.sin(p0 + 2*np.pi*(tt*f0 + 0.5*bw*tt*tt)) for tt in times]
    return chirp


def generate_hfilter(nbins):
    '''
    Generates frequency-space hanning filter from 2.5-17.5 MHz.
    '''
    min_freq = round(2.5 * nbins / 50)
    max_freq = round(17.5 * nbins / 50)
    hamming = np.sin(np.linspace(0, 1, num=max_freq-min_freq+1) * np.pi)
    hfilter = np.hstack((np.zeros(min_freq),
                         hamming*2,
                         np.zeros(nbins-2*max_freq-1),
                         np.zeros(hamming.size), np.zeros(min_freq-1)))
    return hfilter


def dechirp_gen(data, rchirp, nbins, use_hamming=True, LO_correction=None, bandpass=False):
    '''

    '''
    reference_chirp = generate_reference_chirp()
    if bandpass:
        flip_chirp = reference_chirp
    else:
        flip_chirp = np.flipud(reference_chirp)

    fft_chirp = np.fft.fft(flip_chirp, n=nbins)
        
    if use_hamming:
        hfilter = generate_hfilter(nbins)
        fft_filtered = np.multiply(fft_chirp, hfilter)
    else:
        fft_filtered = fft_chirp

    for elem in data:
        fft_elem = np.fft.fft(elem, n=nbins)
        if LO_correction is not None:
            fft_elem[640] -= LO_correction
        # TODO: This is where to add in the LO correction ...
        yield np.fft.ifft(np.multiply(fft_filtered, fft_elem), n=nbins)

        
# TODO: Auto calculate the frequency correction using the quietest
# 25k raw traces of the image? I currently use the bottom 800 samples;
# no reason not to go smaller and use the input samples
# (we currently throw out the bottom 227)
# Alternatively, would it be enough to use the first 50, before the chirp starts?
# The slowest signal is 10MHz, which is at 1/5 the number of samples.
# However, with so few samples, we wouldn't get a good FFT approximation ...


def save_file(traces, filename, mag_scale):
    '''
    * Traces: generator that yields numpy array of samples
    * filename: full path of file to save to
    * mag_scale: how many counts per dB of trace. 
    '''
    with open(filename, 'wb') as fd:
        for trace in traces:
            scaled_mag = np.int32(mag_scale * 20 * np.log10(np.abs(trace)))
            scaled_mag.byteswap(True)
            scaled_mag.tofile(fd)


def find_LO_params(filename, channel, insamples, trace_start, trace_end, nh3=True):
    '''
    Attempts to automatically determine the magnitude and phase of the LO's
    energy leakage. Only the component at 10MHz matters for pik1.
    # TODO: If I want to use the samples from 3200-3437, I have to properly adjust the phase!
    '''
    if nh3:
        input_sweeps = load_raw_nh3_subset_gen(
            filename, channel, np.arange(trace_start, trace_end), insamples, 3200, 200)
    else:
        input_sweeps = load_raw_jh1_subset_gen(filename, np.arange(trace_start, trace_end))
            
    stack_depth = trace_end - trace_start
    noise_data = [trace[-800:] for trace in input_sweeps]
    noise_fft = np.fft.fft(np.sum(noise_data, axis=0), n=800)
    # Extract the component at 10MHz
    # Multiply by 4 to account for switching between a 800 and 3200-bin FFT.
    # Divide by number of traces that were summed to get per-trace value.
    correction = 4*noise_fft[160]/(stack_depth)
    return correction


def find_quiet_region(filename, channel, insamples, outsamples, nh3=True):
    '''
    Searches through entire file for block of traces with the quietest bottom
    part to use as the baseline for local oscillator noise.
    Returns the bounds = (start, end)  in raw traces, using python's indexing convention.
    (use data[start:end] or range(start, end))
    Divide by 50 to convert to indices for the pik1 product.
    '''
    num_traces = get_num_traces(filename, insamples, nh3)
    # Require region to be 10k (raw) traces wide
    width = 10000
    # Require region to be 800 samples high. The 
    nbins_noise = 800
    # We only support conversion for 3200 bins in pik1. If we need others,
    # it is possible to extend it, but should chose the number of noise bins
    # to evenly divide the output samples.
    assert outsamples == 3200
    nbins_pik1 = outsamples
    # TODO: I need to improve find_LO_params to adjust phase before it's OK
    # to calculate quietest region including samples after 3200.
    #trace_gen = load_raw_nh3_gen(filename, channel, insamples, insamples)
    trace_gen = load_raw_gen(filename, channel=channel,
                             num_in_samples=insamples, num_out_samples=outsamples,
                             nh3=nh3)
    stack_depth = 50 # required for efficiency.
    stack_gen = coherent_stack_gen(trace_gen, stack_depth)
    trace_idx = 0
    trace_bottom = np.nan*np.zeros(num_traces/stack_depth)
    Rchirp = np.fft.fft(np.flipud(generate_reference_chirp()), n=outsamples)
    for trace in stack_gen:
        # TODO: This is hacky. Use function. Have to take last 3200 samples of trace ...
        trace_dechirp = np.fft.ifft(np.multiply(np.fft.fft(trace[-outsamples:], n=outsamples), Rchirp))
        # I'm not sure what's theoretically right here...linear or log scale? 
        trace_bottom[trace_idx] = np.sum(np.abs(trace_dechirp[-nbins_noise:]))
        trace_idx += 1
        if trace_idx % 1000 == 0:
            print trace_idx
    print "Cycled through traces"
    print trace_bottom
    trace_sum = np.nan*np.zeros(num_traces/stack_depth - width/stack_depth + 1)
    for trace_idx in np.arange(0, len(trace_sum)):
        trace_sum[trace_idx] = np.sum(trace_bottom[trace_idx:trace_idx+width/stack_depth])
    print "computed sums"
    print trace_sum
    # The quietest region is 
    min_trace_idx = np.argmin(trace_sum)

    trace_start = stack_depth * min_trace_idx
    trace_end = stack_depth*(min_trace_idx+width/stack_depth)
    print "Quietest region is %d - %d (%d - %d)" % (trace_start, trace_end, trace_start/50, trace_end/50)
    return (trace_start, trace_end)


def run_pik1(args):
    # TODO: Add in automatically calculating the LO correction 
    # and passing it to dechirp_gen. (Don't forget to multiply by args.stackdepth)

    # if args.jh1:
    #     input_sweeps = load_raw_jh1_gen(args.infile, args.outsamples, args.blanking)
    # else:
    #     input_sweeps = load_raw_nh3_gen(args.infile, args.channel, args.insamples, args.outsamples, args.blanking)
    input_sweeps = load_raw_gen(args.infile, args.channel, args.insamples, args.outsamples, args.blanking, nh3=(not args.jh1))
    coh_stacks = coherent_stack_gen(input_sweeps, args.stackdepth)
    dechirped_traces = dechirp_gen(coh_stacks, generate_reference_chirp(),
                                   args.outsamples, use_hamming=True)
    inco_stacks = inco_stack_gen(dechirped_traces, args.incodepth)
    save_file(inco_stacks, args.outfile, args.scale)
            

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Pulse compress radar data')
    parser.add_argument('--infile', required=True,
                        help='Full path to 2-byte radar file. (bxds from breakout)')
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--channel', type=int, required=True)
    parser.add_argument('--insamples', type=int, default=3437,
                        help='How many samples in each input trace.')
    parser.add_argument('--outsamples', type=int, default=3200,
                        help='How many samples to write in each output trace.')
    parser.add_argument('--blanking', type=int, default=150,
                        help='How many samples to blank at start of each trace.')
    parser.add_argument('--stackdepth', type=int, required=True)
    parser.add_argument('--incodepth', type=int, required=True)
    parser.add_argument('--scale', type=int, default=1000,
                        help='Output scale. (Counts per dB)')
    parser.add_argument('--jh1', action='store_true')
    args = parser.parse_args()
    run_pik1(args)
