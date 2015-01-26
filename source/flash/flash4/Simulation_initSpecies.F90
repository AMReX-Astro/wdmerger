subroutine Simulation_initSpecies()
implicit none
#include "Multispecies.h"
#include "Flash.h"

call Multispecies_setProperty(C12_SPEC, A, 12.)
call Multispecies_setProperty(C12_SPEC, Z, 6.)

call Multispecies_setProperty(O16_SPEC, A, 16.)
call Multispecies_setProperty(O16_SPEC, Z, 8.)

call Multispecies_setProperty(NE22_SPEC, A, 22.)
call Multispecies_setProperty(NE22_SPEC, Z, 10.)

end subroutine Simulation_initSpecies
