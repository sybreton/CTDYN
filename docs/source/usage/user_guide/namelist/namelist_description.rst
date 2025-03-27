Input namelists 
###############

A short description of the parameters that have to 
be provided by the user in the namelists is given
here. 

global
-------

sr : real      
  Stellar radius in solar units. Default ``1.d0``.

rotp : real    
  Rotation period, in solar units. Default ``1.d0``. This parameter is used only 
  to compute the amplitude of the magnetic diffusivity ``eta`` in the following
  way:

  :math:`C_\Omega = \Omega_\star R_\star^2 / \eta_t` .

grid
----

np : integer
  Number of radial mesh points. Default ``50``.

n_theta : integer
  Number of latitudinal mesh points. Default ``301``.

na : integer
  Number of radial order for the ``An`` decomposition
  of the magnetic field. Default ``49``.

profiles
---------

regime : character    
  Specify the rotation profile and the magnetic diffusivity profile to use.
  Possible regimes are ``r86`` (Radler 1986), ``ns``, ``dj`` (Dudley & James 1989),
  ``rp72`` (Roberts 1972, model 1), ``corona``, ``r72`` (Roberts 1972, model 2),
  ``br`` (Roberts 1972, Braginski), ``sk69`` (Steenbeck & Krause 1969), ``st74``
  (Stix 1974), ``h2``, ``h4``, ``h5``, ``k2`` (K-type giant model), ``r2``. 

xa1 : real
  Lower boundary parameter of the turbulent layer. Default ``0.64d0``.

xa2 : real
  Upper boundary parameter of the turbulent layer. Default ``0.72d0``

xda1 : real
  First thickness parameter of the turbulent layer. Default ``0.025d0``.

xda2 : real
  Second thickness parameter of the turbulent layer. Default ``0.025d0``.

s0 : real    
  First stream function coefficient. Default ``0.85d0``.

s2 : real  
  Second stream function coefficient. Default ``0.08d0``. 

xb : real
  Stream function location parameter. Default ``0.65d0``.

gd : real      
  Exponent term in stream function. Default ``1.2d0``.

edr : real        
  Magnetic diffusivity ratio ``eta_c/eta_top``. Default ``0.1d0``.

xe1 : real         
  Location parameter for magnetic diffusivity. Default ``0.7d0``.

xde1 : real       
  Thickness parameter for magnetic diffusivity. Default ``0.025d0``.

dd1 : real     
  Thickness parameter for rotation profiles for
  following ``regime`` value: ``h2``, ``h4``. 
  Default ``0.05d0``.

xc1 : real     
  Location parameter for rotation profiles for
  following ``regime`` value: ``sk69``, ``stm``, ``h2``, ``h4``. 
  Default ``0.7d0``.

c2_h2 : real     
  Coefficient of ``cos(theta)^2`` term at the surface, such as
  ``omega_surface=omega_eq-c2*cos(theta)^2`` in the ``h2`` regime.
  Default ``0.2d0``.

oco : real     
  Core-to-surface equatorial rotation frequency ratio in the ``h2`` regime. 
  Default ``0.9d0``.

brent
------

accu : real   
  Convergence criteria on the imaginary part of the eigenvalue
  to validate a stationary solution. Solutions with imaginary
  parts smaller than ``accu`` will be accepted and the bisection
  search will be interrupted. Default ``1d-3``.

al_i : real
  Value (or lower bound if ``al_i`` is different from ``al_f``) 
  for the ``C_alpha`` dynamo number. We recall that 
  ``C_alpha = R * alpha * eta_T``. Default ``0.d0``.

al_f : real
  Upper bound for ``C_alpha``. If ``al_i`` is different
  from ``al_f``, CTDYN will search by bisection a
  critical ``C_alpha`` value corresponding to a real
  eigenvalues. Default ``10.d0``.


outputs
-------

dir : character
  directory where the output files will be written.

nj : integer
  Number of time slices.

write_vectors : logical   
  If set to ``.false.``, the code will only compute eigenvalues.
  If set to ``.true.``, the code will compute both eigenvalues and 
  eigenvectors.

xbt : real   
  Location (in stellar radius) to consider to generate the butterfly
  diagram. Default ``0.75d0``. 

xbo : real     
  Default ``1.5d0``.

zeta_r : real  
  External radius at which the base of the corona is located and
  the field lines open. Any value lower than ``1.1`` stellar radius will be 
  replaced by ``1.1`` stellar radius. Default value is ``1.3d0``.   


boundaries
-----------

x_in : real    
  Internal bound to consider for the dynamo problem (in stellar radius). 
  The equation will be solved considering a grid on ``[x_in, 1]``.
  Default ``0.58d0``.

physics
--------

aqu : real     
  Set to ``1`` to include the ``alpha**2`` term in the equation, or
  to ``0`` to ignore it. Default ``1``.

bct : real    
  If ``0``, the equations are solved assuming perfect conductor
  hypothesis. Set to ``1`` otherwise (default value). Default ``1``.

beta_i : real 
  --

c3 : real      
  If set to ``1``, the ``cos^3`` term will be included in the 
  ``alpha B`` quantity. This is only valid if ``m = 0``. Default
  value is ``0``.

ffree : real 
  Force-free external boundary condition. Currently this option
  is not implemented and changing the value of ``ffree`` has no
  effect. 

hd : real      
  Whether to include turbulent pumping, ``1`` or not, ``0``.

xm : real     
  Exponent to use in the power law ``R_M \propto Omega**xm``,
  where ``R_M`` is the Reynolds number of the meridian circulation
  and ``Omega`` is the angular velocity of the star.


fields
-------

degree : character    
  Specify the angular symmetry of the solution to search, ``'d'`` for 
  dipole modes and ``q`` for quadrupole modes. Default ``d``.

mmm : real     
  Azimuthal wave number. Default ``0``.

controls
--------

flg : real     
  If set to ``1``, lapack solver for complex matrixes will be used, if
  set to ``0``, the solver for real matrixes will be used. 
  If ``mmm`` is different from ``0``, ``flg`` is automatically set to
  ``1``.

rm_i : real
  Rotation independent component of the meridional circulation 
  Reynolds number. 

rm_f : real 
  Rotation dependent component of the meridional circulation 
  Reynolds number. 

co : real 
  Adimensioned rotation dynamo coefficient.
