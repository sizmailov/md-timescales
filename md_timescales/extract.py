def extract_mass_center(path_to_trajectory: str, output_filename: str) -> None:
    """

    :param path_to_trajectory:
    :param output_filename: name of output four-column .csv file [time_ns,x,y,z]
    """
    ...


def extract_NH_autocorr(path_to_trajectory: str, output_directory: str) -> None:
    """

    :param path_to_trajectory:
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    ...


def extract_CH3_autocorr(path_to_trajectory: str, output_directory: str) -> None:
    """

    :param path_to_trajectory:
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    ...


def extract_inertia_tensor_vectors_autocorr(path_to_trajectory: str, output_directory: str) -> None:
    """

    :param path_to_trajectory:
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    ...
