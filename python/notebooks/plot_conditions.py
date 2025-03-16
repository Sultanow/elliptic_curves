import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from sympy import primerange

def plot_condition_matrix(df_primes, condition_column, cmap="viridis"):
    """
    Erstellt und plottet die Matrix für eine gegebene Condition.

    :param df_primes: DataFrame mit p, q und Condition-Spalten
    :param condition_column: Name der Condition-Spalte, z. B. 'condition_31'
    :return: Die generierte Matrix (numpy array) und die Prime-Index-Zuordnung
    """
    # Filtere nach gültigen Einträgen
    df_condition_met = df_primes[df_primes[condition_column] == 1]

    # Erstelle Mapping für Primzahlen → Matrix-Indizes
    unique_primes = sorted(set(df_primes['p'].unique()).union(set(df_primes['q'].unique())))
    prime_index_map = {prime: idx for idx, prime in enumerate(unique_primes)}

    # Initialisiere die Matrix
    matrix_size = len(unique_primes)
    matrix = np.zeros((matrix_size, matrix_size), dtype=int)

    # Mapppe p- und q-Werte auf Matrix-Indices
    p_indices = df_condition_met['p'].map(prime_index_map).astype(int).to_numpy()
    q_indices = df_condition_met['q'].map(prime_index_map).astype(int).to_numpy()

    # Fülle die Matrix
    matrix[p_indices, q_indices] = 1
    # matrix[q_indices, p_indices] = 1  # Falls symmetrische Darstellung gewünscht

    # Achsenbeschriftungen vorbereiten
    max_prime = max(unique_primes)
    axis_labels = list(enumerate(primerange(3, max_prime)))
    axis_labels = axis_labels[::len(axis_labels) // 9][:-1] + [axis_labels[-1]]  # Reduziere Labels
    ticks = [e[0] for e in axis_labels]
    ticklabels = [e[1] for e in axis_labels]

    # Plotten der Matrix
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.matshow(matrix, interpolation='nearest', cmap=cmap)

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(ticklabels, rotation=90, fontsize=8)
    ax.set_yticklabels(ticklabels, fontsize=8)
    plt.show()

    return matrix, prime_index_map

def extract_and_plot_lines(matrix, shift_values, threshold=3):
    """
    Extrahiert die Punktlisten für jede verschobene Linie y = x + shift und plottet sie.
    
    :param matrix: Die ursprüngliche Binärmatrix mit den Punkten.
    :param shift_values: Liste der Verschiebungswerte für die y=x Geraden.
    :param threshold: Maximaler Abstand, um Punkte als zur Linie zugehörig zu betrachten.
    :return: Dictionary mit shift-Werten als Keys und zugehörigen Punkten als Listen von (x, y)-Tupeln.
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.matshow(matrix, interpolation='nearest', cmap="Greys")

    # Extrahiere die Koordinaten der gesetzten Punkte und filtere nur obere Matrixhälfte (y <= x)
    points = np.column_stack(np.where(matrix == 1))
    points = points[points[:, 0] <= points[:, 1]]  # Nur obere Dreiecksmatrix

    # Dictionary zur Speicherung der Linienpunkte
    line_points_dict = {}

    # Farben für die Linien
    colors = plt.cm.jet(np.linspace(0, 1, len(shift_values)))

    for shift, color in zip(shift_values, colors):
        # Berechne den Abstand jedes Punktes zur Linie y = x + shift
        distances = np.abs(points[:, 0] - points[:, 1] + shift) / np.sqrt(2)

        # Wähle Punkte, die nahe genug an der verschobenen Linie liegen
        close_points = points[distances < threshold]

        # Speichere die Punkte für diese Linie
        line_points_dict[shift] = [tuple(point) for point in close_points]

        # Plotte diese Punkte in der Farbe der verschobenen Linie
        ax.scatter(close_points[:, 1], close_points[:, 0], color=color, s=5)

        # Zeichne die verschobene Gerade nur in der oberen Hälfte
        x_vals = np.linspace(0, matrix.shape[1], 100)
        y_vals = x_vals + shift
        valid_mask = y_vals <= x_vals  # Nur Werte in der oberen Hälfte
        ax.plot(x_vals[valid_mask], y_vals[valid_mask], color=color, linestyle='dashed', linewidth=1)

    # Achsenbeschriftungen entfernen
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.show()

    return line_points_dict

def extract_and_plot_custom_circular_arc(matrix, x0, y0, fixed_radius, threshold=3):
    """
    Markiert Punkte, die am nächsten zum Kreis mit einem festen Mittelpunkt (x0, y0) und einem festen Radius liegen.

    :param matrix: Die ursprüngliche Binärmatrix mit den Punkten.
    :param x0: Vorgegebener x-Wert für den Mittelpunkt.
    :param y0: Vorgegebener y-Wert für den Mittelpunkt.
    :param fixed_radius: Fester Radius des Kreises.
    :param threshold: Maximaler Abstand, um Punkte als zum Kreisbogen zugehörig zu betrachten.
    :return: Liste der Punkte, die dem Kreisbogen am nächsten liegen, und den verwendeten Mittelpunkt.
    """
    fixed_center = (x0, y0)  # Vorgegebener Mittelpunkt

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.matshow(matrix, interpolation='nearest', cmap="Greys")

    # Extrahiere die Koordinaten der gesetzten Punkte
    points = np.column_stack(np.where(matrix == 1))

    # Berechne den Abstand jedes Punktes zum Kreis
    distances = np.abs(np.sqrt((points[:, 1] - fixed_center[0]) ** 2 + (points[:, 0] - fixed_center[1]) ** 2) - fixed_radius)

    # Wähle Punkte, die nahe genug am Kreisbogen liegen und innerhalb der Matrix bleiben
    close_points = points[(distances < threshold) &
                          (points[:, 0] >= 0) & (points[:, 0] < matrix.shape[0]) &
                          (points[:, 1] >= 0) & (points[:, 1] < matrix.shape[1])]

    # Plotte diese Punkte in Farbe
    ax.scatter(close_points[:, 1], close_points[:, 0], color="red", s=3)

    # Zeichne den Kreisbogen
    theta = np.linspace(-np.pi / 2, np.pi / 2, 100)  # Nur obere Hälfte des Kreises
    x_arc = fixed_center[0] + fixed_radius * np.cos(theta)
    y_arc = fixed_center[1] + fixed_radius * np.sin(theta)

    # Nur Punkte innerhalb der Matrix plotten
    mask = (x_arc >= 0) & (x_arc < matrix.shape[1]) & (y_arc >= 0) & (y_arc < matrix.shape[0])
    ax.plot(x_arc[mask], y_arc[mask], color="red", linestyle='dashed', linewidth=0.5)

    # Achsenbeschriftungen entfernen
    ax.set_xticks([])
    ax.set_yticks([])

    plt.show()

    return [tuple(point) for point in close_points], fixed_center

def extract_and_plot_custom_line(matrix, slope, intercept, threshold=3, highlight_points=True):
    """
    Highlights points that are closest to a given straight line of the form y = slope * x + intercept.

    :param matrix: The binary matrix with points.
    :param slope: The slope (m) of the line y = mx + c.
    :param intercept: The intercept (c) of the line y = mx + c.
    :param threshold: Maximum distance to consider a point as belonging to the line.
    :param highlight_points: Boolean flag to highlight points close to the line (default: True).
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.matshow(matrix, interpolation='nearest', cmap="Greys")

    # Extract the coordinates of the set points
    points = np.column_stack(np.where(matrix == 1))

    # Compute distances of points from the line y = slope * x + intercept
    distances = np.abs(points[:, 0] - (slope * points[:, 1] + intercept))

    # Select points that are close enough to the line
    close_points = points[distances < threshold]

    # Conditionally highlight the selected points
    if highlight_points:
        ax.scatter(close_points[:, 1], close_points[:, 0], color="red", s=3)

    # Draw the guideline line
    x_vals = np.linspace(0, matrix.shape[1], 100)
    y_vals = slope * x_vals + intercept
    mask = (y_vals >= 0) & (y_vals < matrix.shape[0])
    ax.plot(x_vals[mask], y_vals[mask], color="red", linestyle='dashed', linewidth=0.5)

    # Remove axis labels
    ax.set_xticks([])
    ax.set_yticks([])

    plt.show()