void wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s);

void volInBoundary (Real               time,
  	            Real*              com_p,
		    Real*              com_s,
		    Real&              vol_p,
		    Real&              vol_s,
                    Real               rho_cutoff);
