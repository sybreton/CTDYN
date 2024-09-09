import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
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
                         contour_ls="-", colorbar=True, label=None,
                         fill_outside=False, show_up_bounds=False) :
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

    contour_ls : str
      Contour linestyle. Optional, default ``"-"``.

    colorbar : bool
      Set to ``True`` to include the colorbar in the figure.
      Optional, default ``True``.

    label : str
      Colorbar label. Optional, default ``None``.

    fill_outside : bool
      If set to ``True``, the colormap will be applied also
      for region above surface radius. Optional, default,
      ``False``.

    show_up_bounds : bool
      If ``True``, show the upper limit of the plotted mesh.
      Optional, default ``False``.

    Returns
    -------
    matplotlib.Figure
      The created figure.
    """
    x = r * np.sin (theta)
    y = r * np.cos (theta)
    
    fig, ax = plt.subplots (1, 1, figsize=figsize)
    
    # Drawing boundaries of the slice
    ax.plot (x[:,0], y[:,0], color="black", lw=1, 
             zorder=10)
    if show_up_bounds :
      ax.plot (x[:,-1], y[:,-1], color="black", lw=1,
               zorder=10)

    # Drawing surface boundary
    arc = matplotlib.patches.Arc ((0,0), 2, 2, theta1=270, theta2=90,
                                  linewidth=1, color="black", zorder=10)
    ax.add_patch (arc)

    if fill_outside :
      mesh_m = mesh
    else :
      mask = x**2+y**2 > 1
      mesh_m = mesh.copy()
      mesh_m[mask] = np.nan
    
    if mode=="pcolormesh" :
        im = ax.pcolormesh (x, y, mesh_m, cmap=cmap)
    elif mode=="contourf" :
        im = ax.contourf (x, y, mesh_m, cmap=cmap)
    else :
        raise Exception ("Accepted arguments for mode are 'pcolormesh' or 'contourf'")
    if contour :
        ax.contour (x, y, mesh, colors="darkgrey", 
                    linestyles=contour_ls)
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

def read_radial_profiles (filename) :
    """
    Read radial profiles such as alpha 
    and eta turbulent coefficients.
    
    Parameters
    ----------
    filename : str or Path object
      Filename 
    
    Returns
    -------
    pandas.DataFrame
      A dataframe with the radial profiles.
    """
    data = np.loadtxt (filename)
    columns = ["alpha", "alpha_p", 
               "eta_1", "eta_2", "eta_3", "psi", 
               "ar", "arp", "bt", "btp", 
               "om0", "om0p", "om2", "om2p"]
    df = pd.DataFrame (data=data[:,1:], 
                       index=data[:,0],
                       columns=columns)
    return df

def plot_alpha (df, figsize=(5,8), 
                xlabel=None, ylabels=None,
                lw=1) :
    """
    Plot alpha profile. 
    
    Parameters
    ----------
    df : pandas.DataFrame
      Dataframe read from ``read_radial_profiles``.
      
    figsize : tuple
      Figure size. Optional, default ``(5,8)``.
      
    xlabel : str
      X-axis label. Optional, default ``None``.
      
    ylabels : tuple of str
      Y-axis label. Optional, default ``None``.
      
    lw : float
      Line width.
      
    Returns
    -------
    matplotlib.Figure
      The created figure.
    """
    if xlabel is None :
        xlabel = r"$r$ ($R_\star$)"
    if ylabels is None :
        ylabels = (r"$\alpha$", r"$\alpha_p$")
    fig, (ax1, ax2) = plt.subplots (2, 1, figsize=figsize)
    ax2.set_xlabel (xlabel)
    ax1.set_ylabel (ylabels[0])
    ax2.set_ylabel (ylabels[1])
    ax1.plot (df.index, df["alpha"], color="blue", lw=lw, 
              label=r"$\alpha$")
    ax2.plot (df.index, df["alpha_p"], color="darkorange", lw=lw,
              label=r"$\alpha_p$")
    return fig

def plot_eta (df, figsize=(5,10), 
              xlabel=None, ylabels=None,
              lw=1) :
    """
    Plot eta profile. 
    
    Parameters
    ----------
    df : pandas.DataFrame
      Dataframe read from ``read_radial_profiles``.
      
    figsize : tuple
      Figure size. Optional, default ``(5,5)``.
      
    xlabel : str
      X-axis label. Optional, default ``None``.
      
    ylabels : tuple of str
      Y-axis label. Optional, default ``None``.
      
    lw : float
      Line width.
      
    Returns
    -------
    matplotlib.Figure
      The created figure.
    """
    if xlabel is None :
        xlabel = r"$r$ ($R_\star$)"
    if ylabels is None :
        ylabels = (r"$\eta$", 
                   r"$\mathrm{d} \eta / \mathrm{d}x$", 
                   r"$\mathrm{d}^2 \eta / \mathrm{d}^2 x$")
    fig, (ax1, ax2, ax3) = plt.subplots (3, 1, figsize=figsize)
    ax3.set_xlabel (xlabel)
    ax1.set_ylabel (ylabels[0])
    ax2.set_ylabel (ylabels[1])
    ax3.set_ylabel (ylabels[2])
    ax1.plot (df.index, df["eta_1"], color="blue", lw=lw) 
    ax2.plot (df.index, df["eta_2"], color="darkorange", lw=lw)
    ax3.plot (df.index, df["eta_3"], color="gold", lw=lw)
    return fig
