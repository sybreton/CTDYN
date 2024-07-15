CTDYN input namelist description
################################

A short description of the parameters that have to 
be provided by the user in the namelist is given
here. 

Parameters
----------

dir : str
  directory where the output files will be written.

al_i : float
  Value (or lower bound if ``al_i`` is different from ``al_f``) 
  for the ``C_alpha`` dynamo number. We recall that 
  ``C_alpha = R * alpha * eta_T``.

al_f : float
  Upper bound for ``C_alpha``. If ``al_i`` is different
  from ``al_f``, CTDYN will search by bisection a
  critical ``C_alpha`` value corresponding to a real
  eigenvalues. 

nsa : float
  ???

xa1 : float
  ???

xa2 : float
  ???

xa3 : float
  ???

xb : float
  ???

xda1 : float
  ???

xda2 : float
  ???

rm_i : float
  ???

rm_f : float 
  ???

nsr : float 
  ???

cm_i : float 
  ???

cm_f : float 
  ???

nso : float      
  ???

edr : float        
  ???

xe1 : float         
  ???

xde1 : float       
  ???

x_in : float    
  Internal bound to consider for the dynamo problem (in stellar radius). 
  The equation will be solved considering a grid on ``[x_in, 1]``

c3 : float      
  ???

bct : float    
  ???

accu : float   
  Convergence criteria on the imaginary part of the eigenvalue
  to validate a stationary solution. Solutions with imaginary
  parts smaller than ``accu`` will be accepted and the bisection
  search will be interrupted.

ans1 : str    
  Specify the rotation profile, the magnetic diffusivity profile to use.
  **TODO: List the possible values and profiles affected.** 

ans2 : str    
  Specify the angular symmetry of the solution to search, ``'d'`` for 
  dipole modes and ``q`` for quadrupole modes.

ans3 : str    
  ???

ans4 : str    
  If set to ``'n'``, the code will only compute eigenvalues.
  If set to ``'v'``, the code will compute both eigenvalues and 
  eigenvectors.

s0 : float    
  ???

s2 : float  
  ???

s4 : float  
  ???

s6 : float   
  ???

a2p : float     
  ???

a4p : float   
  ???

mmm : float     
  Azimuthal wave number.

hd : float      
  Whether to include turbulent pumping, ``1`` or not, ``0``.

sr : float      
  Stellar radius in solar units.

rotp : float    
  Rotation period, in solar units.

gd : float      
  ???

aqu : float     
  Set to ``1`` to include the ``alpha**2`` term in the equation, or
  to ``0`` to ignore it. 

flg : float     
  If set to ``1``, lapack solver for complex matrixes will be used, if
  set to ``0``, the solver for real matrixes will be used. 
  If ``mmm`` is different from ``0``, ``flg`` is automatically set to
  ``1``.

dd1 : float     
  ???

rc1 : float     
  ???

rc2 : float     
  ???

oco : float     
  ???

xm : float     
  Exponent to use in the power law ``R_M \propto Omega**xm``,
  where ``R_M`` is the Reynolds number of the meridian circulation
  and ``Omega`` is the angular velocity of the star.

beta_i : float 
  ???

beta_f : float  
  ???

beta_s : float 
  ???

zeta_r : float  
  ???

ffree : float 
  ???
  
xbt : float   
  ???
  
xbo : float     
  ???
