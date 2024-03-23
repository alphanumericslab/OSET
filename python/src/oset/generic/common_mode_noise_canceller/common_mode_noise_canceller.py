import numpy as np
import matplotlib.pyplot as plt
from oset.generic.lp_filter import lp_filter_zero_phase, tikhonov_regularization

def common_mode_noise_canceller(x, fs, fc=None, itr=10, lambda_val=100, optim_indexes=None, plot_results=False):
    """
    Estimates and removes common-mode noise from multichannel data.

    Parameters:
        x: Input signal matrix (N channels x T samples).
        fs: Sampling frequency.
        fc: Cutoff frequency for low-pass filtering.
        itr: Number of iterations for the algorithm.
        lambda_val: Regularization parameter for Tikhonov regularization.
        optim_indexes: Indices of the samples to optimize.
        plot_results: Flag to enable/disable result plotting.

    Returns:
        x_den: The denoised signal, with common-mode noise removed.
        x_den_smoothed: The denoised signal after smoothing with Tikhonov regularization.
        s: The estimated common-mode signal.
        alpha: The scaling factor or matrix used in noise cancellation.
        bl: Channel-wise average or baseline estimated by high-pass filtering.
    """
    if optim_indexes is None:
        optim_indexes = np.arange(x.shape[1])

    if fc is not None:
        if fc == 0:
            bl = np.mean(x, axis=1, keepdims=True)
            x = x - bl
        else:
            bl = lp_filter_zero_phase(x, fc/fs)
            x = x - bl
    else:
        bl = np.array([])

    # Initial estimate of common-mode signal `s`
    s = np.median(x, axis=0).reshape(1, -1)
    s = s[:, optim_indexes]

    # print(x.shape)
    # print(bl.shape)
    # print(s.shape)
    

    # Iteratively solve for `alpha` and `s`
    for k in range(itr):
        # print(x[:, optim_indexes].shape)
        # print(s.shape)
        alpha = np.linalg.lstsq(s.T, x[:, optim_indexes].T, rcond=None)[0].T
        s = np.linalg.lstsq(alpha, x[:, optim_indexes], rcond=None)[0]

    # Calculate denoised signal and smoothed version
    s = np.linalg.lstsq(alpha, x, rcond=None)[0]
    x_den = x - np.dot(alpha, s)

    # print(x.shape)
    # print(s.shape)
    # print(alpha.shape)
    # print(x_den.shape)
    
    if lambda_val > 0:
        x_den_smoothed = tikhonov_regularization(x_den, [2], lambda_val)
    else:
        x_den_smoothed = x_den

    # Optional plotting of results
    if plot_results:
        plt.figure()
        plt.plot(x.T)
        plt.hold(True)
        plt.plot(s, 'k', linewidth=2)
        plt.grid(True)
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.title('Original Signals with Estimated Common-Mode Signal')
        plt.gca().set_fontsize(14)

        plt.figure()
        plt.subplot(211)
        plt.plot(x.T)
        plt.grid(True)
        plt.title('Original Signals')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.gca().set_fontsize(14)

        plt.subplot(212)
        plt.plot(x_den.T)
        plt.grid(True)
        plt.title('Denoised Signals')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.gca().set_fontsize(14)
        plt.show()

    return x_den, x_den_smoothed, s, alpha, bl

if __name__ == "__main__":
    # Import argparse and set up argument parsing
    import argparse
    parser = argparse.ArgumentParser(description="Estimates and removes common-mode noise from multichannel data.")
    # Add arguments for each parameter, for example:
    parser.add_argument('--fs', type=float, required=True, help='Sampling frequency')
    parser.add_argument('--fc', type=float, default=0.001, help='Cutoff frequency for low-pass filtering')
    parser.add_argument('--itr', type=int, default=10, help='Number of iterations for the algorithm')
    parser.add_argument('--lambda_val', type=float, default=100, help='Regularization parameter for Tikhonov regularization')
    # You'll need to decide how to pass `x` and potentially `optim_indexes`, as they are more complex than simple scalar values

    args = parser.parse_args()

    # Since 'x' and 'optim_indexes' are complex to handle via command-line arguments, you might need to load or define them inside this block
    # For the sake of the example, let's assume 'x' is loaded from a file or defined here directly

    # Example placeholder for loading 'x' and defining 'optim_indexes'
    x = np.random.randn(5, 100)  # Dummy data for demonstration
    optim_indexes = np.arange(25, 75)  # Example index range

    # Call the function with the parsed arguments (and placeholders for 'x' and 'optim_indexes')
    data_den, data_den_smoothed, s, alpha, bl = common_mode_noise_canceller(
        x=x,
        fs=args.fs,
        fc=args.fc,
        itr=args.itr,
        lambda_val=args.lambda_val,
        optim_indexes=optim_indexes,
        plot_results=True  # Assuming you want to plot results when running as a script
    )