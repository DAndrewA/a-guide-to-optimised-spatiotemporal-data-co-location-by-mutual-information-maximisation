"""Author: Andrew Martin
Creation date: 25/11/25

Script containing the functionality to compute the local normalised across-track density of orbits for a given satellite at a given latitude
"""

import numpy as np

def calculate_dlambda(latitude_1, latitude_2) -> float:
    """Given a maximum latitude, latitude_2, and an equatorwards latitude, latitude_1, calculcate the speeration of their longitudes along a great circle.
    NOTE: the point with lambda_1 in an orbit inclined at 92 degrees will be eastwards of the point at lambda_2, so dlambda is made negative to ensure proper sign.
    """
    # alpha_2 = pi/2 => tan(alpha_2) -> +infinity
    # => -cos(phi_2)sin(phi_1) + sin(phi_2)cos(phi_1)cos(dlambda) -> 0+
    cos_dlambda = np.tan(latitude_1) / np.tan(latitude_2)
    dlambda = -np.arccos(cos_dlambda)
    return dlambda

def heading_along_great_cirlce(latitude_1: float, latitude_2: float) -> float:
    """Given a current latitude, latitude_1, and a maximum latitude, latitude_2, along a great-circle, calculate the heading to travel from the point along the great-circle"""
    dlambda = calculate_dlambda(latitude_1=latitude_1, latitude_2=latitude_2)
    y = np.cos(latitude_2)*np.sin(dlambda)
    x = np.cos(latitude_1)*np.sin(latitude_2) - np.sin(latitude_1)*np.cos(latitude_2)*np.cos(dlambda)
    heading = np.arctan2(y,x)
    return heading

def heading_along_great_cirlce_degrees(latitude_1: float, latitude_2: float) -> float:
    return np.rad2deg(heading_along_great_cirlce(
        latitude_1 = np.deg2rad(latitude_1),
        latitude_2 = np.deg2rad(latitude_2)
    ))

def orbital_density(latitude_1: float, latitude_2: float) -> float:
    return (
        1
        / np.cos(heading_along_great_cirlce(
            latitude_1=latitude_1,
            latitude_2=latitude_2
        ))
        / np.cos(latitude_1)
    )

def orbital_density_degrees(latitude_1: float, latitude_2:float) -> float:
    return orbital_density(
        latitude_1 = np.deg2rad(latitude_1),
        latitude_2 = np.deg2rad(latitude_2)
    )

def normalised_orbital_density(latitude_1: float, latitude_2: float) -> float:
    """Normalised such that the density when latitude_1=0 is 1"""
    return np.sin(latitude_2) * orbital_density(latitude_1=latitude_1, latitude_2=latitude_2)

def normalised_orbital_density_degrees(latitude_1: float, latitude_2:float) -> float:
    return normalised_orbital_density(
        latitude_1 = np.deg2rad(latitude_1),
        latitude_2 = np.deg2rad(latitude_2)
    )
