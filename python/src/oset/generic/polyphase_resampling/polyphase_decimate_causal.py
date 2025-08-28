import numpy as np

def polyphase_decimate_causal(in_sig, dec_filter_h, dec_rate, init_cond=None):
    """
    A causal polyphase FIR decimator with continuous block boundaries.

    Parameters
    ----------
    in_sig : array-like
        Input signal (1D list or numpy array).
    dec_filter_h : array-like
        FIR decimation filter's impulse response.
    dec_rate : int
        Decimation factor.
    init_cond : array-like, optional
        Previous input block's trailing samples equal to the filter's memory (same length as dec_filter_h - 1).

    Returns
    -------
    in_sig_decimated : np.ndarray
        Decimated output signal.
    final_cond : np.ndarray
        Last dec_filter_h-1 samples to feed as init_cond for the next block (for continuous processing avoiding edge-effects).
        
    Refs:
        - Vaidyanathan, P. P. (1993). Multirate Systems and Filter Banks. Englewood Cliffs, NJ: Prentice Hall.
        - Crochiere, R. E., & Rabiner, L. R. (1983). Multirate digital signal processing. Prentice Hall.
        
    Reza Sameni, 2025
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET    
    
    """

    x_in = np.asarray(in_sig, dtype=float)
    h = np.asarray(dec_filter_h, dtype=float)
    M = int(dec_rate)
    if M <= 0:
        raise ValueError("dec_rate must be a positive integer")
    L = h.size
    if L == 0:
        return np.array([]), np.array([])

    # State handling
    if init_cond is None:
        z = np.zeros(max(L - 1, 0), dtype=float)
    else:
        z = np.asarray(init_cond, dtype=float)
        if z.size != max(L - 1, 0):
            raise ValueError("init_cond must have length len(dec_filter_h)-1")

    # Concatenate state and block
    x = np.concatenate([z, x_in])
    N = x.size

    # How many valid decimated outputs?
    num_outputs = (N - L) // M + 1
    if num_outputs <= 0:
        return np.array([]), x[-(L - 1):].copy() if L > 1 else np.array([])

    # Polyphase branches: e_m[k] = h[kM + m]
    phases = [h[m::M] for m in range(M)]

    y = np.zeros(num_outputs, dtype=float)
    offset = L - 1  # first output time index in concatenated buffer

    for n in range(num_outputs):
        idx0 = offset + n * M  # "time n x M" in x
        acc = 0.0
        for m in range(M):
            taps = phases[m]
            K = taps.size
            if K == 0:
                continue
            # Indices: idx0 - m - M*np.arange(K)
            idxs = idx0 - m - M * np.arange(K)
            # All idxs are guaranteed in-bounds because of the prepended state
            seg = x[idxs]
            acc += np.dot(taps, seg)
        y[n] = acc

    final_cond = x[-(L - 1):].copy() if L > 1 else np.array([])
    return y, final_cond
