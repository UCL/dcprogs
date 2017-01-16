#include <iostream>

#include <likelihood/missed_eventsG.h>
#include <likelihood/idealG.h>
#include <likelihood/occupancies.h>
 
int main() {

  // Define parameters.
  HJCFIT::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  HJCFIT::QMatrix qmatrix(matrix, /*nopen=*/2);
  HJCFIT::t_real const tau(1e-4); // in seconds

  // Create missed-events G
  HJCFIT::MissedEventsG eG(qmatrix, tau);
  // Create ideal G
  HJCFIT::IdealG idealG(qmatrix);
  
  HJCFIT::t_real const tcritical(5e-3);

  std::cout << "Equilibrium Occupancies\n"
            << "=======================\n\n"
            << "Ideal Likelihood\n"
            << "----------------\n\n"
            << "  * initial: " << HJCFIT::occupancies(idealG) << "\n"
            << "  * final: "   << HJCFIT::occupancies(idealG, false) << "\n\n\n"
            << "Missed-events Likelihood\n"
            << "------------------------\n\n"
            << "  * initial: " << HJCFIT::occupancies(eG) << "\n"
            << "  * final: "   << HJCFIT::occupancies(eG, false) << "\n\n\n\n"
            << "CHS Occupancies\n"
            << "===============\n\n"
            << "Missed-events Likelihood\n"
            << "------------------------\n\n"
            << "  * tcritical: " << tcritical << "\n"
            << "  * initial: " << HJCFIT::CHS_occupancies(eG, tcritical) << "\n"
            << "  * final: "   << HJCFIT::CHS_occupancies(eG, tcritical, false) << "\n";

  return 0;
}
