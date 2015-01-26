subroutine Simulation_initSpecies()
implicit none
#include "Multispecies.h"
#include "Flash.h"

call Multispecies_setProperty(C12_SPEC, A, 12.)
call Multispecies_setProperty(C12_SPEC, Z, 6.)

end subroutine Simulation_initSpecies
