#!/usr/bin/env python3
# shmfall.py
# Real-time time-aligned waterfall plot for WIN-format (Urabe 1994) seismic data.
#
# This script processes multi-channel seismic time series formatted in the
# WIN (Waveform INtegration) format as defined by Urabe (1994). It is designed
# to operate on the plain-text output produced by **shmdump**, which converts
# WIN shared-memory segments into a human-readable stream.
#
# The script reads timestamped 1-second WIN packets from stdin, time-aligns
# all channels using their packet timestamps, and visualizes the continuous
# data stream as a scrolling waterfall plot (time vs. channels) using matplotlib.
#
# Key characteristics:
#   • Input stream is assumed to be WIN-format records decoded by shmdump.
#   • Each packet begins with a timestamp line (YY MM DD hh mm ss msec),
#     followed by any number of channel records:
#         CH_ID   NSAMPLES   v0 v1 v2 ...
#   • All channels are resampled/decimated to a fixed target sampling rate.
#   • Samples are mapped into a circular time-aligned buffer using the
#     absolute timestamp, allowing late packets to be placed correctly
#     within the visible time range.
#
# This implementation is intended for real-time monitoring systems that
# use WIN shared memory (WIN/EX) and the standard WIN utilities, and is
# compatible with shmdump’s deterministic text output.
#
# Reference:
#   Urabe, T. (1994). "A common format for multi-channel earthquake waveform data."
#
# Features:
# - Time alignment using absolute timestamps (late packets are placed at
#   their correct time positions as long as they are within the visible window).
# - Fixed sampling rate per second; packets with too many/few samples are
#   decimated/resampled to match the expected sampling rate.
# - Optional normalization (z-score / robust z-score).
# - Optional envelope (Hilbert) or RMS processing.
# - Optional channel sorting by longitude using a channel table.
# - Circular (ring) buffer in time, with old data dropped outside the visible window.
#
# Expected input format per "packet":
#   1) A timestamp line with 7 fields:
#        YY MM DD hh mm ss msec
#      (YY is year - 2000; i.e., 25 -> 2025)
#
#   2) One or more "channel lines" of the form:
#        CH_ID NSAMPLES v0 v1 v2 ...
#      where CH_ID is a channel identifier (string),
#      NSAMPLES is the number of samples for this 1-second block,
#      and v0.. are integer sample values.
#
#   This pattern (timestamp line + multiple channel lines) repeats over time.
#
# Note:
# - The script assumes a fixed target sampling rate SAMPLING_RATE (e.g., 100 Hz).
# - If NSAMPLES != SAMPLING_RATE, the data is resampled/decimated to SAMPLING_RATE.
# - The "absolute sample index" is derived from the packet timestamp and a
#   reference time t0 (the first timestamp encountered).
# - Data older than the current head index minus MAX_SAMPLES is discarded.
#
# Usage example:
#   shmdump -tq 11 | python shmfall.py -f channel_table.txt -d 60 -n 0 -e 0
#
# Copyright (c) 2025 Hironori Otsuka
# This software is released under the MIT License.

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import sys
import threading
from collections import deque   # kept for potential backward compatibility / future use
from datetime import datetime, timedelta
import argparse
from scipy.signal import decimate, resample, hilbert
import math

'''
import matplotlib
# Try to select an interactive backend if available.
# If it fails (no Tk/Qt installed), it will fall back silently.
for _backend in ("QtAgg", "TkAgg", "GTK3Agg"):
    try:
        matplotlib.use(_backend)
        break
    except Exception:
        continue
'''
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# =====================
#  Basic processing utilities
# =====================

def root_mean_square(data):
    """
    Compute the RMS amplitude along the time axis using a sliding window.

    Parameters
    ----------
    data : np.ndarray
        2D array with shape (n_traces, n_samples).

    Returns
    -------
    np.ndarray
        RMS of the data with the same shape as input, padded with NaNs
        at both ends so that the RMS window is centered.
    """
    window_size = 10

    # If the trace is shorter than the RMS window, just return a float copy
    # (this is a safety fallback).
    if data.shape[1] < window_size:
        return data.astype(float)

    windows = sliding_window_view(data, window_shape=window_size, axis=1)
    rms = np.sqrt(np.mean(windows**2, axis=2))

    pad_before = window_size // 2
    pad_after = window_size - pad_before - 1

    return np.pad(
        rms,
        ((0, 0), (pad_before, pad_after)),
        mode='constant',
        constant_values=np.nan
    )


def haversine(lat1, lon1, lat2, lon2):
    """
    Compute great-circle distance between two points given in degrees,
    using the haversine formula.

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : float
        Latitude and longitude in degrees.

    Returns
    -------
    float
        Distance in kilometers.
    """
    R = 6371.0  # Earth radius [km]

    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (
        math.sin(dphi / 2) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    )
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def calc_envelope(data):
    """
    Compute the signal envelope (instantaneous amplitude) using the Hilbert transform.

    Parameters
    ----------
    data : np.ndarray
        2D array (n_traces, n_samples).

    Returns
    -------
    np.ndarray
        Envelope (absolute value of analytic signal), same shape as input.
    """
    analytic = hilbert(data, axis=1)   # Hilbert transform across time axis
    env = np.abs(analytic)             # Instantaneous amplitude (envelope)
    return env


def zscore(x, axis=None):
    """
    Standard z-score normalization.

    Parameters
    ----------
    x : np.ndarray
        Input array.
    axis : int or tuple of ints, optional
        Axis or axes along which to compute mean and std.

    Returns
    -------
    np.ndarray
        Z-scored array with mean 0 and std 1 along the given axis.
    """
    mean = x.mean(axis=axis, keepdims=True)
    std = np.std(x, axis=axis, keepdims=True)

    # Protect against division by zero.
    std[std == 0] = 1e-8

    return (x - mean) / std


def robust_zscore(x, axis=None):
    """
    Robust z-score normalization using median and MAD.

    Parameters
    ----------
    x : np.ndarray
        Input array.
    axis : int or tuple of ints, optional
        Axis or axes along which to compute median and MAD.

    Returns
    -------
    np.ndarray
        Robust z-scored array where outliers are downweighted.
    """
    median = np.median(x, axis=axis, keepdims=True)
    mad = np.median(np.abs(x - median), axis=axis, keepdims=True)

    # Protect against division by zero.
    mad[mad == 0] = 1e-8

    return (x - median) / mad


def adjust_samples(samples, expected_samples, sampling_rate):
    """
    Adjust the number of samples within a 1-second packet to match
    the target sampling rate.

    If the packet contains more samples than the target sampling_rate,
    it is downsampled via decimation (with an integer factor).
    If the packet contains fewer samples than the target sampling_rate,
    it is resampled using FFT-based resampling.

    Parameters
    ----------
    samples : np.ndarray
        1D array of sample values for a given channel and 1-second block.
    expected_samples : int
        NSAMPLES written in the input line for this block.
    sampling_rate : int
        Target sampling rate (e.g., 100 samples per second).

    Returns
    -------
    np.ndarray
        1D array of length approximately equal to sampling_rate.
    """
    if expected_samples > sampling_rate:
        factor = int(expected_samples / sampling_rate)
        if factor > 1:
            samples = decimate(samples, factor)
        # If factor == 1, leave as is.
    elif expected_samples < sampling_rate:
        samples = resample(samples, sampling_rate)
    return samples


# =====================
#   StreamReader
# =====================

class StreamReader:
    """
    Read timestamped stream data from stdin, time-align it, and write it
    into a 2D ring buffer: buffer[num_channels, MAX_SAMPLES].

    Time-alignment model
    --------------------
    - t0: reference timestamp (the first packet timestamp encountered).
    - packet_time: timestamp of the currently processed packet.
    - head_idx: latest absolute sample index (monotonically increasing).
      For SAMPLING_RATE=100 Hz, 1 second corresponds to +100 in index.

    For each packet:
    - Compute delta_t = (packet_time - t0) in seconds.
    - Compute base_idx = round(delta_t * sampling_rate).
    - Each k-th sample in the packet is mapped to absolute index idx = base_idx + k.

    Circular buffer mapping
    -----------------------
    - The ring buffer stores the last MAX_SAMPLES samples.
    - Each absolute index idx is mapped into the buffer via:
        buf_pos = idx % MAX_SAMPLES
    - Samples older than (head_idx - MAX_SAMPLES) are discarded, i.e.,
      not written into the buffer.

    This allows late-arriving packets to be rendered at the correct time
    positions as long as they still fall within the visible window.
    """
    def __init__(self, buffer, ch_index, sampling_rate, max_samples):
        """
        Parameters
        ----------
        buffer : np.ndarray
            2D array of shape (num_channels, MAX_SAMPLES) used as a ring buffer.
        ch_index : dict
            Mapping from channel ID (string) to row index (int).
        sampling_rate : int
            Target sampling rate in Hz.
        max_samples : int
            Total number of samples stored in the ring buffer
            (DURATION * sampling_rate).
        """
        self.buffer = buffer              # shape: (num_channels, MAX_SAMPLES)
        self.ch_index = ch_index          # ch_id -> row index
        self.sampling_rate = sampling_rate
        self.max_samples = max_samples
        self.current_time = "00:00:00"
        self._stop = False

        self.t0 = None                    # reference timestamp
        self.packet_time = None           # timestamp of the most recent time line
        self.head_idx = -1                # latest absolute sample index

    def stop(self):
        """Request the reading loop to stop."""
        self._stop = True

    def _update_time(self, data):
        """
        Parse a timestamp line (7 elements) and update current_time, t0, and
        packet_time.

        The timestamp line has the form:
            YY MM DD hh mm ss msec
        """
        year = int(data[0]) + 2000
        month, day, hour, minute, second = map(int, data[1:6])
        millisecond = int(data[6])

        base_dt = datetime(year, month, day, hour, minute, second)
        timestamp = base_dt + timedelta(milliseconds=millisecond)

        self.current_time = timestamp.strftime("%Y-%m-%d %H:%M:%S")

        # If we have no reference time yet, use this packet as t0.
        if self.t0 is None:
            self.t0 = timestamp

        # This timestamp applies to subsequent channel lines until
        # the next timestamp line is read.
        self.packet_time = timestamp

    def _write_channel_samples(self, ch_id, expected_samples, samples):
        """
        Write one channel's samples for a 1-second block into the ring buffer.

        Parameters
        ----------
        ch_id : str
            Channel identifier.
        expected_samples : int
            NSAMPLES given in the input line.
        samples : np.ndarray
            1D array of raw sample values.
        """
        # Skip channels that are not listed in the channel table.
        if ch_id not in self.ch_index:
            return

        # Adjust the number of samples to the target sampling rate.
        samples = adjust_samples(samples, expected_samples, self.sampling_rate)

        # If we don't have a valid timestamp yet, skip writing.
        if self.packet_time is None or self.t0 is None:
            return

        # Compute absolute index of the first sample in this packet.
        delta_t = (self.packet_time - self.t0).total_seconds()  # [sec]
        base_idx = int(round(delta_t * self.sampling_rate))

        ch_row = self.ch_index[ch_id]

        for k, value in enumerate(samples):
            idx = base_idx + k

            # Discard samples that are too old to be visible
            # in the current ring buffer window.
            if self.head_idx >= 0 and idx < self.head_idx - self.max_samples:
                continue

            # Update head index if this is the most recent sample so far.
            if idx > self.head_idx:
                self.head_idx = idx

            buf_pos = idx % self.max_samples
            self.buffer[ch_row, buf_pos] = value

    def read_stream(self):
        """
        Main reading loop.

        Reads lines from stdin, distinguishes between timestamp lines
        and channel lines, and writes data into the ring buffer.
        """
        while not self._stop:
            line = sys.stdin.readline()
            if not line:
                continue
            line = line.strip()
            if not line:
                continue

            data = line.split()

            # Timestamp line (7 fields)
            if len(data) == 7:
                self._update_time(data)

            # Channel data line (CH_ID, NSAMPLES, samples...)
            elif len(data) > 2:
                ch_id = data[0].upper()
                try:
                    expected_samples = int(data[1])
                except ValueError:
                    # Invalid NSAMPLES; skip this line.
                    continue
                try:
                    samples = np.array([int(x) for x in data[2:]], dtype=float)
                except ValueError:
                    # Non-integer sample in the list; skip this line.
                    continue

                self._write_channel_samples(ch_id, expected_samples, samples)

    def read_sharemem(self):
        """
        Placeholder for future shared-memory reader implementation.
        """
        raise NotImplementedError("Shared memory reading is not implemented now.")


# =====================
#   Main plotting logic
# =====================

def parse_data_stream(channel_table, vmin, vmax, normalize, envelope,
                      ignore_missing, sort_channels, ref_point, _save_debug):
    """
    Main function that:
    - Reads the channel table and determines channel order.
    - Creates the ring buffer and the StreamReader.
    - Sets up the matplotlib figure and animation.
    - Periodically updates the waterfall image from the ring buffer.
    """
    global matrix  # kept for compatibility / debug (the actual data is in buffer)

    channel_ids = []

    # Example reference coordinates (Yonaguni Island).
#    LAT_0, LON_0 = 24.4545, 122.9325
    LAT_0, LON_0 = ref_point

    # Load channel table and optionally sort by longitude.
    if sort_channels >= 0:
        with open(channel_table, 'r') as f:
            lines = [line.split() for line in f]

        # Expect channel table lines with at least 15 columns where:
        #   [0]   : channel ID
        #   [13]  : latitude
        #   [14]  : longitude
        valid_lines = [
            (line[0], float(line[13]), float(line[14]))
            for line in lines if len(line) > 14
        ]
        
        if sort_channels == 0:
            # Sort by longitude.
            valid_lines.sort(key=lambda x: x[2])
        elif sort_channels == 1:
            # Alternative example: sort by distance from a reference point.
            valid_lines.sort(
             key=lambda x: haversine(LAT_0, LON_0, x[1], x[2])
            )

        channel_ids = [ch[0].upper() for ch in valid_lines]
    else:
        # Use channel IDs in the order listed in the table.
        with open(channel_table, 'r') as f:
            channel_ids = [line.split()[0].upper() for line in f]

    num_channels = len(channel_ids)
    ch_index = {ch: i for i, ch in enumerate(channel_ids)}

    # 2D ring buffer for time-aligned data.
    # Initialize with zeros (interpreted as "no signal" for ignore_missing).
    buffer = np.zeros((num_channels, MAX_SAMPLES), dtype=float)

    # StreamReader shares the ring buffer and channel index.
    reader = StreamReader(
        buffer=buffer,
        ch_index=ch_index,
        sampling_rate=SAMPLING_RATE,
        max_samples=MAX_SAMPLES
    )

    fig, ax = plt.subplots(figsize=(9, 6), dpi=100)
    matrix = np.zeros((num_channels, MAX_SAMPLES), dtype=float)  # debug / legacy only

    fig._saved_debug = _save_debug

    # Configure y-axis ticks in seconds based on total duration.
    if DURATION <= 15:          # <= 15 seconds: tick every 1 second
        tick_interval_sec = 1
    elif DURATION <= 120:       # <= 2 minutes: tick every 10 seconds
        tick_interval_sec = 10
    elif DURATION <= 360:       # <= 6 minutes: tick every 30 seconds
        tick_interval_sec = 30
    else:                       # > 6 minutes: tick every 60 seconds
        tick_interval_sec = 60

    tick_interval_samples = tick_interval_sec * SAMPLING_RATE
    y_ticks = np.arange(0, MAX_SAMPLES + 1, tick_interval_samples)
    y_labels = [f"{int(i / SAMPLING_RATE)}" for i in y_ticks]

    ax.set_ylim(y_ticks[-1] + 1, 0)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    ax.set_xlim(0, num_channels)

    # Select colormap based on envelope option (same as original behavior).
    cmap = 'rainbow' if envelope >= 0 else 'seismic'

    # Initial image; data will be updated by the animation.
    im = ax.imshow(
        matrix.T,
        aspect='auto',
        interpolation='spline16',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )

    # Colorbar label based on normalization / envelope settings.
    labelstring = ""
    if normalize == 0:
        labelstring += "Z-score"
    elif normalize == 1:
        labelstring += "Robust Z-score"
    else:
        labelstring += "Raw value"

    if envelope == 0:
        labelstring += " (Envelope)"
    elif envelope == 1:
        labelstring += " (RMS)"

    plt.colorbar(im, label=labelstring)
    plt.title("Real-time Waterfall")

    print("Matplotlib backend:", plt.get_backend(), flush=True) # debug info

    # -------------------------
    # Animation update function
    # -------------------------
    def update_matrix(_):
        nonlocal buffer
        num_ch, width = buffer.shape

        # Reorder the ring buffer so that columns go from "oldest" to "newest"
        # within the visible window.
        if reader.head_idx < 0:
            # No data yet; just show the initial buffer.
            image_matrix = buffer.copy()
        else:
            # Position in the ring buffer of the latest sample.
            head_pos = reader.head_idx % width

            # Construct an index array that starts just after head_pos
            # and wraps around. The resulting order is:
            #   [head_pos+1, ..., width-1, 0, 1, ..., head_pos]
            idx = np.arange(width)
            rolled_idx = (head_pos + 1 + idx) % width

            image_matrix = buffer[:, rolled_idx]

        # Optionally drop channels that are entirely zero in the visible window.
        if ignore_missing:
            mask = ~np.all(image_matrix == 0, axis=1)
            image_matrix = image_matrix[mask, :]
            if image_matrix.size == 0:
                # If all channels are zero, draw a dummy image to avoid errors.
                image_matrix = np.ones((num_ch, width), dtype=float)
            # Update x-axis range to match the visible number of traces
            visible_ch = image_matrix.shape[0]
            ax.set_xlim(0, visible_ch)
#            print(f"Visible channels: {visible_ch}/{num_ch}")
        else:
            # keep original full range
            ax.set_xlim(0, num_ch)

        # Apply normalization, if requested.
        if normalize == 0:
            image_matrix = zscore(image_matrix, axis=1)
        elif normalize == 1:
            image_matrix = robust_zscore(image_matrix, axis=1)

        # Apply envelope or RMS processing, if requested.
        if envelope == 0:
            image_matrix = calc_envelope(image_matrix)
        elif envelope == 1:
            image_matrix = root_mean_square(image_matrix)

        # Update plot title with current time (from last timestamp line).
        ax.set_title(
            f"Real-time Waterfall - Time: {reader.current_time}"
        )

        # Flip vertically so that the latest time appears at the bottom
        # (keeping backward compatibility with previous visualization).
        im.set_data(np.flipud(image_matrix.T))

        # Save a single debug frame once (optional cross-environment comparison).
        if not getattr(fig, "_saved_debug", False):
            # Make a debug_matrix with channels along vertical axis
            debug_matrix = np.flipud(image_matrix.T)

            fig.canvas.draw()
            fig.savefig("shmfall_latest_debug.png")

            # ---------------------------------------------------------------
            # This block is to create a image for the 
            # original paper of this software.
            # Normally, this is skipped.
            # Let external module create a sliced debug image
            #try:
            #    import shmfall_debugview
            #    shmfall_debugview.plot_matrix_slice(
            #        debug_matrix,
            #        sampling_rate=SAMPLING_RATE,
            #        vmin=vmin,
            #        vmax=vmax,
            #        cbar_label=labelstring,
            #        cmap=cmap,
            #    )
            #except Exception as exc:
            #    print(f"[debugview] failed: {exc}", file=sys.stderr)
            # ---------------------------------------------------------------
            
            fig._saved_debug = True

        return im,

    # Start the reading thread (daemon so it ends with the main program).
    thread = threading.Thread(target=reader.read_stream, daemon=True)
    thread.start()

    ani = animation.FuncAnimation(
        fig, update_matrix, interval=UPDATE_INTERVAL, blit=False
    )
    plt.xlabel("Channels")
    plt.ylabel("Elapsed Time [sec]")

    try:
        plt.show()
    except KeyboardInterrupt:
        print("Ctrl+C pressed. Exiting...")
    finally:
        reader.stop()
        thread.join(timeout=1)
        plt.close(fig)


# =====================
#   Entry point
# =====================

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        '-f', '--channel_table',
        required=True,
        help="Path to channel table file."
    )
    p.add_argument(
        '-v', '--vminmax',
        help="Color scale of plot. Format is [vmin:vmax]"
    )
    p.add_argument(
        '-d', '--duration',
        type=int,
        default=60,
        help="Duration of the plot (seconds)"
    )
    p.add_argument(
        '-n', '--normalize',
        type=int,
        choices=[0, 1],
        default=-1,
        help="Normalize with Z-score [0] or robust Z-score [1]"
    )
    p.add_argument(
        '-e', '--envelope',
        type=int,
        choices=[0, 1],
        default=-1,
        help="Plot envelope [0] or RMS [1] of the signal"
    )
    p.add_argument(
        '-i', '--ignore_zeros',
        action='store_true',
        help="Ignore missing channels (all-zero traces)"
    )
    p.add_argument(
        '-s', '--sort_channels',
        type=int,
        choices=[0, 1],
        default=-1,
        help="Sort channels by longitude [0] or distance from reference point [1] (requires lat/lon columns in channel table)"
    )
    p.add_argument(
        '-r', '--ref_point',
        help="Reference point (lat,lon) for distance-based sorting Format is [lat/lon]"
    )
    p.add_argument(
        '--debug',
        action='store_true',
        help="Save a debug frame for cross-environment comparison"
    )

    a = p.parse_args()
    channel_table = a.channel_table
    vmin, vmax = (-1000, 1000) if not a.vminmax else map(float, a.vminmax.split(':'))
    ref_point = [24.4545, 122.9325] if not a.ref_point else list(map(float, a.ref_point.split('/')))  # Example reference coordinates (Yonaguni Island).

    # Normalize options.
    if a.normalize == 0:
        normalize = 0
        if not a.vminmax:
            vmin, vmax = -6, 6
    elif a.normalize == 1:
        normalize = 1
        if not a.vminmax:
            vmin, vmax = -6, 6
    else:
        normalize = -1

    # Envelope / RMS options.
    if a.envelope == 0:
        envelope = 0
    elif a.envelope == 1:
        envelope = 1
    else:
        envelope = -1

    # Default vmin/vmax for envelope / RMS if not explicitly given.
    if envelope >= 0 and not a.vminmax:
        vmin, vmax = (0, 12) if normalize >= 0 else (0, 3000)

    if a.debug:
        _saved_debug = False
    else:
        _saved_debug = True

    ignore_zeros = a.ignore_zeros
    sort_channels = a.sort_channels

    # Global constants.
    DURATION = a.duration
    SAMPLING_RATE = 100
    MAX_SAMPLES = DURATION * SAMPLING_RATE
    UPDATE_INTERVAL = 100           # ms between plot updates
    SAMPLE_INTERVAL = 1 / SAMPLING_RATE

    parse_data_stream(
        channel_table,
        vmin,
        vmax,
        normalize,
        envelope,
        ignore_zeros,
        sort_channels,
        ref_point,
        _saved_debug
    )
