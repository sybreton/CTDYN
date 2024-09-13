Input namelists 
###############

A short description of the parameters that have to 
be provided by the user in the namelists is given
here. 

global
-------

sr : real      
  Stellar radius in solar units.

rotp : real    
  Rotation period, in solar units.


profiles
---------

regime : character    
  Specify the rotation profile and the magnetic diffusivity profile to use.
  Possible regimes are ``r86`` (Radler 1986), ``ns``, ``dj`` (Dudley & James 1989),
  ``rp72`` (Roberts 1972, model 1), ``corona``, ``r72`` (Roberts 1972, model 2),
  ``br`` (Roberts 1972, Braginski), ``sk69`` (Steenbeck & Krause 1969), ``st74``
  (Stix 1974), ``h2``, ``h4``, ``h5``, ``k2`` (K-type giant model), ``r2``. 

s0 : real    
  ???

s2 : real  
  ???

s4 : real  
  ???

s6 : real   
  ???

a2p : real     
  ???

a4p : real   
  ???

xa1 : real
  Lower boundary parameter of the turbulent layer.

xa2 : real
  Upper boundary parameter of the turbulent layer.

xda1 : real
  First thickness parameter of the turbulent layer.

xda2 : real
  Second thickness parameter of the turbulent layer.

xb : real
  Stream function location parameter.

gd : real      
  Exponent term in stream function.

edr : real        
  Magnetic diffusivity ratio ``eta_c/eta_top``. 

xe1 : real         
  Location parameter for magnetic diffusivity.

xde1 : real       
  Thickness parameter for magnetic diffusivity.

dd1 : real     
  Thickness parameter for rotation profiles for
  following ``regime`` value: ``h2``, ``h4``. 

xc1 : real     
  Location parameter for rotation profiles for
  following ``regime`` value: 
  ``sk69``, ``stm``, ``h2``, ``h4``. 

c2_h2 : real     
  Coefficient of ``cos(theta)^2`` term at the surface, such as
  ``omega_surface=omega_eq-c2*cos(theta)^2`` in the ``h2`` regime.

oco : real     
  Additional rotation parameter in the ``h2`` regime.

brent
------

accu : real   
  Convergence criteria on the imaginary part of the eigenvalue
  to validate a stationary solution. Solutions with imaginary
  parts smaller than ``accu`` will be accepted and the bisection
  search will be interrupted.

al_i : real
  Value (or lower bound if ``al_i`` is different from ``al_f``) 
  for the ``C_alpha`` dynamo number. We recall that 
  ``C_alpha = R * alpha * eta_T``.

al_f : real
  Upper bound for ``C_alpha``. If ``al_i`` is different
  from ``al_f``, CTDYN will search by bisection a
  critical ``C_alpha`` value corresponding to a real
  eigenvalues. 


outputs
-------

dir : character
  directory where the output files will be written.

write_vectors : logical   
  If set to ``.false.``, the code will only compute eigenvalues.
  If set to ``.true.``, the code will compute both eigenvalues and 
  eigenvectors.

zeta_r : real  
  ???

xbt : real   
  ???
  
xbo : real     
  ???

boundaries
-----------

x_in : real    
  Internal bound to consider for the dynamo problem (in stellar radius). 
  The equation will be solved considering a grid on ``[x_in, 1]``

physics
--------

aqu : real     
  Set to ``1`` to include the ``alpha**2`` term in the equation, or
  to ``0`` to ignore it. 

bct : real    
  If ``0``, the equations are solved assuming perfect conductor
  hypothesis. Set to ``1`` otherwise (default value).

beta_i : real 
  ???

c3 : real      
  If set to ``1``, the ``cos^3`` term will be included in the 
  ``alpha B`` quantity. This is only valid if ``m = 0``. Default
  value is ``0``.

ffree : real 
  ???

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
  dipole modes and ``q`` for quadrupole modes.

mmm : real     
  Azimuthal wave number.

controls
--------

flg : real     
  If set to ``1``, lapack solver for complex matrixes will be used, if
  set to ``0``, the solver for real matrixes will be used. 
  If ``mmm`` is different from ``0``, ``flg`` is automatically set to
  ``1``.

nso : real      
  Number of step in the loop to explore the influence of rotation
  over meridional circulation. The bisection procedure will be executed 
  at each iteration. At a given iteration ``ii``, the rotation 
  coefficient ``co`` is ``co = cm_i + ii / (nso+1) * (cm_f - cm_i)``
  and the meridional circulation coefficient is ``c_u = rm_i + rm_f*co**xm``.

rm_i : real
  Rotation independent component of the meridional circulation 
  Reynolds number. 

rm_f : real 
  Rotation dependent component of the meridional circulation 
  Reynolds number. 

cm_i : real 
  Initial rotation coefficient in the exploration loop.

cm_f : real 
  Final rotation coefficient in the exploration loop.
