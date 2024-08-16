import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import re

def read_field (filename, return_meshgrid=True,
                output_format="txt") :
    """
    Read field output produced by CTDYN. 
    """
    if output_format=="txt" :
      r, theta, mesh = read_field_text_file (filename, 
                                             return_meshgrid=return_meshgrid) 
    elif output_format=="hdf5" :
      raise Exception ("hdf5 format is not implemented yet.")
    else :
      raise Exception ("Unknown output format. Supported format are 'txt' and 'hdf5'")
    return r, theta, mesh
    

def read_field_text_file (filename, return_meshgrid=True) :
    """
    Read field from file saved under a text format.
    """
    # Reading the header
    with open (Path (filename), "r") as f :
        head = [next(f).strip() for _ in range(4)]
    # Getting the number of radial and latitudinal meshes 
    n_r, n_theta = re.split (r"\s+", head[2])
    n_r, n_theta = int (n_r), int (n_theta)
    # Getting the radial and latitudinal boundaries
    # and building coordinates array
    r_in, r_out, theta_0, theta_1 = re.split (r"\s+", head[3])
    r_in, r_out, theta_0, theta_1 = (float (r_in),
                                     float (r_out),
                                     float (theta_0),
                                     float (theta_1)
                                    )
    r = np.linspace (r_in, r_out, n_r)
    theta = np.linspace (theta_0, theta_1, n_theta)
    if return_meshgrid :
        r, theta = np.meshgrid (r, theta)
    # Building the array
    mesh = np.loadtxt (filename, skiprows=4)
    mesh = mesh.reshape (n_theta, n_r)
    return r, theta, mesh

def plot_meridional_mesh (r, theta, mesh,
                          figsize=(4,6), cmap="Blues_r",
                          mode="contourf", contour=True,
                          colorbar=True, label=None) :
    """
    Plot a meridional mesh.
    """
    x = r * np.sin (theta)
    y = r * np.cos (theta)
    
    fig, ax = plt.subplots (1, 1, figsize=figsize)
    
    # Drawing boundaries of the slice
    ax.plot (x[:,0], y[:,0], color="black")
    ax.plot (x[:,-1], y[:,-1], color="black")
    
    if mode=="pcolormesh" :
        im = ax.pcolormesh (x, y, mesh, cmap=cmap)
    elif mode=="contourf" :
        im = ax.contourf (x, y, mesh, cmap=cmap)
    else :
        raise Exception ("Accepted arguments for mode are 'pcolormesh' or 'contourf'")
    if contour :
        ax.contour (x, y, mesh, colors="darkgrey", 
                    linestyles="--")
    if colorbar :
        cbar = plt.colorbar (im, shrink=0.5)
        if label is not None :
            cbar.set_label (label)
    ax.axis ("off")

    return fig


