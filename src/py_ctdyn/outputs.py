import numpy as np
import pandas as pd
from astropy.table import Table
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import re

def read_summary_file (filename, output_format="txt",
                       backend="astropy") :
    """
    Read CTDYN summary output.

    Parameters
    ----------
    filename : str or Path object
      Name of the CTDYN output file to read.

    output_format : str
      Format of the output files created by CTDYN.
      Optional, default ``"txt"``.

    backend : str
      Backend to consider for the type of output
      to return. If ``"pandas"`` is chosen, the output
      will be a ``pandas.Dataframe``, if ``"astropy"``
      is chosen, the output will be an ``astropy.table.Table``. 
    
    Returns
    -------
    pandas.DataFrame or astropy.table.Table
      A Dataframe or a Table (depending on the chosen
      backend)  with the elements of the summary file.
    """
    if output_format=="txt" :
      # Extracting the header
      with open (Path (filename), "r") as f :
        first_line = f.readline().strip('\n')
        names = re.split (r"\W+", first_line)
        # First element is an empty string
        names = names[1:]  
      if backend=="pandas" :
        df = pd.read_csv (filename, sep="\s+",  comment="#", 
                          header=None, names=names)
        if "n" in names :
          df = df.set_index ("n")
      elif backend=="astropy" :
        data = np.loadtxt (filename)
        df = Table (data=data, names=names)
      else :
        raise Exception ("Unkown requested backend. Use 'astropy' or 'pandas'")
    elif output_format=="hdf5" :
      raise Exception ("hdf5 format is not implemented yet.")
    else :
      raise Exception ("Unknown output format. Supported format are 'txt' and 'hdf5'")
    return df

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

def plot_meridional_map (r, theta, mesh, ax=None,
                         figsize=(4,6), cmap="seismic",
                         mode="contourf", contour=True,
                         colorbar=True, label=None,
                         fill_outside=False, show_up_bounds=False,
                         vmin=None, vmax=None, **kwargs) :
    """
    Plot a meridional map computed by CTDYN.

    Parameters
    ----------
    r : ndarray
      2d meshgrid with radius. 

    theta : ndarray
      2d meshgrid with colatitude. 

    ax : matplotlib.axes
      Axes on which to draw the figure. If not provided,
      a new figure will be created. Optional, default ``None``.

    mesh : ndarray
      2d meshgrid of the quantity to plot on the map.

    figsize : tuple
      Figure size. Optional, default ``(4,6)``.

    cmap : str or Colormap
      Color map. Optional, default ``seismic``.

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

    fill_outside : bool
      If set to ``True``, the colormap will be applied also
      for region above surface radius. Optional, default,
      ``False``.

    show_up_bounds : bool
      If ``True``, show the upper limit of the plotted mesh.
      Optional, default ``False``.

    vmin : float
      Colormap minimal value. Optional, default ``None``.

    vmax : float
      Colormap maximal value. Optional, default ``None``.

    Returns
    -------
    matplotlib.Figure
      The created figure.
    """
    x = r * np.sin (theta)
    y = r * np.cos (theta)
    
    if ax is None :
      fig, ax = plt.subplots (1, 1, figsize=figsize)
    else : 
      fig = ax.get_figure ()
    
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

    if vmin is None or vmax is None :
      # Getting min and max values
      # If there are negative values in the mesh, we choose 
      # it to be symetric to make it work fine with
      # the colormap default choice that is ``seismic``.
      if np.any (mesh < 0) :
        vmax = np.amax (np.abs (mesh))
        vmin = - vmax
      else :
        vmax = np.amax (mesh)
        vmin = np.amin (mesh)
    norm = matplotlib.colors.Normalize (vmin=vmin, vmax=vmax)
    
    if mode=="pcolormesh" :
        im = ax.pcolormesh (x, y, mesh_m, cmap=cmap, norm=norm,
                            **kwargs)
    elif mode=="contourf" :
        # Avoid a conflict with cmap
        kwargs_contourf = kwargs.copy ()
        kwargs_contourf.pop ("colors", None)  
        im = ax.contourf (x, y, mesh_m, cmap=cmap, norm=norm,
                          **kwargs_contourf)
    else :
        raise Exception ("Accepted arguments for mode are 'pcolormesh' or 'contourf'")
    if contour :
        kwargs.setdefault ("colors", "darkgrey")
        kwargs.setdefault ("linestyles", "-")
        ax.contour (x, y, mesh, **kwargs)
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

def plot_butterfly_diagram (t, theta, mesh, ax=None,
                            figsize=(7,4), cmap="seismic",
                            mode="contourf", contour=True,
                            colorbar=True, label=None,
                            vmin=None, vmax=None, **kwargs) :
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

    ax : matplotlib.axes
      Axes on which to draw the figure. If not provided,
      a new figure will be created. Optional, default ``None``.

    figsize : tuple
      Figure size. Optional, default ``(7,4)``.

    cmap : str or Colormap
      Color map. Optional, default ``seismic``.

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

    vmin : float
      Colormap minimal value. Optional, default ``None``.

    vmax : float
      Colormap maximal value. Optional, default ``None``.

    Returns 
    -------
    matplotlib.Figure
      The created figure.
    """
    lat = 90 - theta/np.pi*180 
    
    if ax is None :
      fig, ax = plt.subplots (1, 1, figsize=figsize)
    else : 
      fig = ax.get_figure ()
    ax.set_yticks ([-60, -30, 0, 30, 60])

    if vmin is None or vmax is None :
      # Getting min and max values
      # If there are negative values in the mesh, we choose 
      # it to be symetric to make it work fine with
      # the colormap default choice that is ``seismic``.
      if np.any (mesh < 0) :
        vmax = np.amax (np.abs (mesh))
        vmin = - vmax
      else :
        vmax = np.amax (mesh)
        vmin = np.amin (mesh)
    norm = matplotlib.colors.Normalize (vmin=vmin, vmax=vmax)

    if mode=="pcolormesh" :
        im = ax.pcolormesh (t, lat, mesh, cmap=cmap, norm=norm,
                            **kwargs)
    elif mode=="contourf" :
        # Avoid a conflict with cmap
        kwargs_contourf = kwargs.copy ()
        kwargs_contourf.pop ("colors", None)  
        im = ax.contourf (t, lat, mesh, cmap=cmap, norm=norm,
                          **kwargs_contourf)
    else :
        raise Exception ("Accepted arguments for mode are 'pcolormesh' or 'contourf'")
    if contour :
        kwargs.setdefault ("colors", "darkgrey")
        kwargs.setdefault ("linestyles", "--")
        ax.contour (t, lat, mesh, **kwargs)
    if colorbar :
        cbar = plt.colorbar (im, shrink=0.8)
        if label is not None :
            cbar.set_label (label)

    ax.set_xlabel (r"$\omega t$")
    ax.set_ylabel (r"Latitude ($\rm ^o$)")
    
    return fig

def read_radial_profiles (filename, backend="astropy") :
    """
    Read radial profiles such as alpha 
    and eta turbulent coefficients.
    
    Parameters
    ----------
    filename : str or Path object
      Filename 

    backend : str
      Backend to consider for the type of output
      to return. If ``"pandas"`` is chosen, the output
      will be a ``pandas.Dataframe``, if ``"astropy"``
      is chosen, the output will be an ``astropy.table.Table``. 
    
    Returns
    -------
    pandas.DataFrame or astropy.table.Table
      A Dataframe or a Table (depending on the chosen
      backend) with the radial profiles.
    """
    data = np.loadtxt (filename)
    columns = ["alpha", "alpha_p", 
               "eta_1", "eta_2", "eta_3", "psi", 
               "ar", "arp", "bt", "btp", 
               "om0", "om0p", "om2", "om2p"]
    if backend=="pandas" :
      df = pd.DataFrame (data=data[:,1:], 
                         index=data[:,0],
                         columns=columns)
    elif backend=="astropy" :
      columns.insert (0, "x")
      df = Table (data=data, names=columns)
    else :
      raise Exception ("Unkown requested backend. Use 'astropy' or 'pandas'")
    return df

def get_xvector (df) :
    """
    Get x-vector from ``df`` depending
    on type (pandas Dataframe or astropy Table).
    """
    if "x" in df.columns :
      x = df["x"]
    else :
      x = df.index
    return x

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
        ylabels = (r"$\alpha$", 
                   r"$\mathrm{d} \alpha / \mathrm{d}x$") 
    fig, (ax1, ax2) = plt.subplots (2, 1, figsize=figsize)
    ax2.set_xlabel (xlabel)
    ax1.set_ylabel (ylabels[0])
    ax2.set_ylabel (ylabels[1])
    x = get_xvector (df)
    ax1.plot (x, df["alpha"], color="blue", lw=lw, 
              label=r"$\alpha$")
    ax2.plot (x, df["alpha_p"], color="darkorange", lw=lw,
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
    x = get_xvector (df)
    ax1.plot (x, df["eta_1"], color="blue", lw=lw) 
    ax2.plot (x, df["eta_2"], color="darkorange", lw=lw)
    ax3.plot (x, df["eta_3"], color="gold", lw=lw)
    return fig
