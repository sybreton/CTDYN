import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import re

def read_field_map (filename, return_meshgrid=True,
                    output_format="txt") :
    """
    Read field output produced by CTDYN. 

    Parameters
    ----------
    filename : str or Path object
      Name of the CTDYN output file to read.

    return_meshgrid : bool
      Whether to return 2d meshgrid (``True``) or 
      1d array (``False``) for radius and colatitude. 
      Optional, default ``True``.
 
    output_format : str
      Format of the output files created by CTDYN.
      Optional, default ``"txt"``.

    Returns 
    -------
    tuple of arrays
      Tuple with, in this order, radius (1d array or 2d meshgrid), colatitude 
      meshgrid (1d array or 2d meshgrid), and field meshgrid.
    """
    if output_format=="txt" :
      r, theta, mesh = read_field_map_text_file (filename, 
                                                 return_meshgrid=return_meshgrid) 
    elif output_format=="hdf5" :
      raise Exception ("hdf5 format is not implemented yet.")
    else :
      raise Exception ("Unknown output format. Supported format are 'txt' and 'hdf5'")
    return r, theta, mesh
    

def read_field_map_text_file (filename, return_meshgrid=True) :
    """
    Read field from CTDYN output file saved under a text format.

    Parameters
    ----------
    filename : str or Path object
      Name of the CTDYN output file to read.

    return_meshgrid : bool
      Whether to return 2d meshgrid (``True``) or 
      1d array (``False``) for radius and colatitude. 
      Optional, default ``True``.
 
    Returns 
    -------
    tuple of arrays
      Tuple with, in this order, radius (1d array or 2d meshgrid), colatitude 
      meshgrid (1d array or 2d meshgrid), and field meshgrid.
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

def plot_meridional_map (r, theta, mesh,
                         figsize=(4,6), cmap="Blues_r",
                         mode="contourf", contour=True,
                         colorbar=True, label=None) :
    """
    Plot a meridional map computed by CTDYN.

    Parameters
    ----------
    r : ndarray
      2d meshgrid with radius. 

    theta : ndarray
      2d meshgrid with colatitude. 

    mesh : ndarray
      2d meshgrid of the quantity to plot on the map.

    figsize : tuple
      Figure size. Optional, default ``(4,6)``.

    cmap : str or Colormap
      Color map. Optional, default ``Blues_r``.

    mode : str
      Map representation mode, either ``contourf``
      or ``pcolormesh``. Optional, default ``contourf``.

    contour : bool
      Set to ``True`` to represent additional contours
      on the map. Optional, default ``True``.

    colorbar : bool
      Set to ``True`` to include the colorbar in the figure.
      Optional, default ``True``.

    label : str
      Colorbar label. Optional, default ``None``.

    Returns :
    matplotlib.Figure
      The created figure.
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

def read_butterfly_diagram (filename, return_meshgrid=True,
                            output_format="txt") :
    """
    Read butterfly diagram output produced by CTDYN. 

    Parameters
    ----------
    filename : str or Path object
      Name of the CTDYN output file to read.

    return_meshgrid : bool
      Whether to return 2d meshgrid (``True``) or 
      1d array (``False``) for radius and colatitude. 
      Optional, default ``True``.
 
    output_format : str
      Format of the output files created by CTDYN.
      Optional, default ``"txt"``.

    Returns 
    -------
    tuple of arrays
      Tuple with, in this order, time (1d array or 2d meshgrid), colatitude 
      meshgrid (1d array or 2d meshgrid), and field meshgrid.
    """
    if output_format=="txt" :
      t, theta, mesh = read_butterfly_diagram_text_file (filename, 
                                                         return_meshgrid=return_meshgrid) 
    elif output_format=="hdf5" :
      raise Exception ("hdf5 format is not implemented yet.")
    else :
      raise Exception ("Unknown output format. Supported format are 'txt' and 'hdf5'")
    return t, theta, mesh

def read_butterfly_diagram_text_file (filename, return_meshgrid=True) :
    """
    Read butterfly diagram from file saved under a text format.

    Parameters
    ----------
    filename : str or Path object
      Name of the CTDYN output file to read.

    return_meshgrid : bool
      Whether to return 2d meshgrid (``True``) or 
      1d array (``False``) for radius and colatitude. 
      Optional, default ``True``.
 
    Returns 
    -------
    tuple of arrays
      Tuple with, in this order, time (1d array or 2d meshgrid), colatitude 
      meshgrid (1d array or 2d meshgrid), and field meshgrid.
    """
    # Reading the header
    with open (Path (filename), "r") as f :
        head = [next(f).strip() for _ in range(1)]
    # Getting the number of latitudinal meshes 
    # (note that it is ntheta-2 that is 
    # actually written in the file)
    n_theta = re.split (r"\s+", head[0])[0]
    n_theta = int (n_theta)
    data = np.loadtxt (filename, skiprows=1)
    theta = data[:n_theta]
    data = data[n_theta:].reshape (-1, n_theta+1)
    t = data[:,0]
    mesh = data[:,1:].T
    if return_meshgrid :
        t, theta = np.meshgrid (t, theta)
    return t, theta, mesh

def plot_butterfly_diagram (t, theta, mesh,
                            figsize=(7,4), cmap="Blues_r",
                            mode="contourf", contour=True,
                            colorbar=True, label=None) :
    """
    Plot a butterfly diagram.

    Parameters
    ----------
    t : ndarray
      2d meshgrid with time coordinate. 

    theta : ndarray
      2d meshgrid with colatitude. 

    mesh : ndarray
      2d meshgrid of the quantity to plot on the map.

    figsize : tuple
      Figure size. Optional, default ``(7,4)``.

    cmap : str or Colormap
      Color map. Optional, default ``Blues_r``.

    mode : str
      Map representation mode, either ``contourf``
      or ``pcolormesh``. Optional, default ``contourf``.

    contour : bool
      Set to ``True`` to represent additional contours
      on the map. Optional, default ``True``.

    colorbar : bool
      Set to ``True`` to include the colorbar in the figure.
      Optional, default ``True``.

    label : str
      Colorbar label. Optional, default ``None``.

    Returns :
    matplotlib.Figure
      The created figure.
    """
    lat = 90 - theta/np.pi*180 
    
    fig, ax = plt.subplots (1, 1, figsize=figsize)
    ax.set_yticks ([-60, -30, 0, 30, 60])

    if mode=="pcolormesh" :
        im = ax.pcolormesh (t, lat, mesh, cmap=cmap)
    elif mode=="contourf" :
        im = ax.contourf (t, lat, mesh, cmap=cmap)
    else :
        raise Exception ("Accepted arguments for mode are 'pcolormesh' or 'contourf'")
    if contour :
        ax.contour (t, lat, mesh, colors="darkgrey",
                    linestyles="--")
    if colorbar :
        cbar = plt.colorbar (im, shrink=0.8)
        if label is not None :
            cbar.set_label (label)

    ax.set_xlabel ("Time (year)")
    ax.set_ylabel (r"Latitude ($\rm ^o$)")
    
    return fig
