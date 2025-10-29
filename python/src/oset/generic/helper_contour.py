import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter
import matplotlib.ticker as ticker
# ===============================================================================================
def generate_kde(df, first_column, second_column):
    x_values = df[first_column]
    y_values = df[second_column]
    values = np.vstack([x_values, y_values])
    kde = gaussian_kde(values)
    x_grid = np.linspace(x_values.min(), x_values.max(), 100)
    y_grid = np.linspace(y_values.min(), y_values.max(), 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(kde(positions).T, X.shape)
    return np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
# ===============================================================================================
def find_contour_level(Z, percentage):
    sorted_Z = np.sort(Z.ravel())[::-1]
    cumsum_Z = np.cumsum(sorted_Z)
    cumsum_Z /= cumsum_Z[-1]
    idx = np.where(cumsum_Z >= percentage)[0][0]
    return sorted_Z[idx]
# ===============================================================================================
def generate_heatmap_contur_plot(data, N, first_column, second_column, XLABEL, YLABEL, TITLE, color):
    data = data[[first_column, second_column]]

    # KDE estimation
    kde = generate_kde(data, first_column, second_column)
    X = kde[:, 0].reshape((100, 100))
    Y = kde[:, 1].reshape((100, 100))
    Z = kde[:, 2].reshape((100, 100))

    mean_x = round(data[first_column].mean(), 1)
    mean_y = round(data[second_column].mean(), 1)
    print(mean_x, mean_y)

    # Axis limits
    x_min = data[first_column].min()
    x_max = data[first_column].max()
    y_min = data[second_column].min()
    y_max = data[second_column].max()

    bins = (int(x_max - x_min), int(y_max - y_min))
    range_ = ((x_min, x_max), (y_min, y_max))

    # Histogram for heatmap
    heatmap_data, _, _ = np.histogram2d(
        data[first_column], data[second_column], bins=bins, range=range_
    )
    heatmap_data_smooth = gaussian_filter(heatmap_data, sigma=1)
    contour_level = find_contour_level(Z, N / 100)

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6))
    img = ax.imshow(
        heatmap_data_smooth.T,
        extent=[x_min, x_max, y_min, y_max],
        origin='lower',
        cmap='viridis',
        interpolation='bicubic'
    )

    ax.contour(X, Y, Z, levels=[contour_level], colors=color, linewidths=2)
    ax.plot(mean_x, mean_y, marker='o', markersize=5, color='red')

    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(TITLE, fontsize=14)
    ax.set_xlabel(XLABEL, fontsize=14)
    ax.set_ylabel(YLABEL, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.grid(which='major', color='black', linewidth=0.9)
    ax.grid(which='minor', color='gray', linewidth=0.1)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
        tick.set_horizontalalignment('right')

    fig.tight_layout()
    cbar = plt.colorbar(img, ax=ax, shrink=0.5)
    cbar.set_label('Population')
    plt.show()
# ===============================================================================================
def generate_contur_plot_ONE(data, N, first_column, second_column, XLABEL, YLABEL, TITLE, color):
    data = data[[first_column, second_column]]

    # KDE estimation
    kde = generate_kde(data, first_column, second_column)
    X = kde[:, 0].reshape((100, 100))
    Y = kde[:, 1].reshape((100, 100))
    Z = kde[:, 2].reshape((100, 100))

    mean_x = round(data[first_column].mean(), 1)
    mean_y = round(data[second_column].mean(), 1)
    print(mean_x, mean_y)

    # Axis limits
    x_min = data[first_column].min()
    x_max = data[first_column].max()
    y_min = data[second_column].min()
    y_max = data[second_column].max()

    bins = (int(x_max - x_min), int(y_max - y_min))
    range_ = ((x_min, x_max), (y_min, y_max))

    # Histogram for heatmap
    heatmap_data, _, _ = np.histogram2d(
        data[first_column], data[second_column], bins=bins, range=range_
    )
    heatmap_data_smooth = gaussian_filter(heatmap_data, sigma=1)
    contour_level = find_contour_level(Z, N / 100)

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.contour(X, Y, Z, levels=[contour_level], colors=color, linewidths=2)
    ax.plot(mean_x, mean_y, marker='o', markersize=5, color='red')

    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(TITLE, fontsize=14)
    ax.set_xlabel(XLABEL, fontsize=14)
    ax.set_ylabel(YLABEL, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.grid(which='major', color='black', linewidth=0.9)
    ax.grid(which='minor', color='gray', linewidth=0.1)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
        tick.set_horizontalalignment('right')

    fig.tight_layout()
    plt.show()  
# ===============================================================================================
def generate_contur_plot_MULTI(data, N, first_column, second_column, ax, color, label):
    data = data[[first_column, second_column]]

    # KDE estimation
    kde = generate_kde(data, first_column, second_column)
    X = kde[:, 0].reshape((100, 100))
    Y = kde[:, 1].reshape((100, 100))
    Z = kde[:, 2].reshape((100, 100))

    mean_x = round(data[first_column].mean(), 1)
    mean_y = round(data[second_column].mean(), 1)
    print(f'{label} Mean: ({mean_x}, {mean_y})')

    contour_level = find_contour_level(Z, N / 100)

    # Plot contour and mean point
    contour_plot = ax.contour(X, Y, Z, levels=[contour_level], colors=color, linewidths=2)
    ax.plot(mean_x, mean_y, marker='o', markersize=5, color=color)    