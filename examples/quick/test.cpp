#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/mKOST.h"

int main (int argc, char* argv[])
{
  mKOST::sStateVector out;
  mKOST::sElements elements, elements2;
  mKOST::sOrbitParam params;

  std::cout << std::endl << std::endl << std::endl;

  /*Alternative test, with an initial param list at t0*/
  elements.a = 5.799119e5;
  elements.e = 0.9;
  elements.i = SIMD_PI / 4;
  elements.LoP = SIMD_PI / 4;
  elements.LaN = 0.2;
  elements.L = 0.0;

  /*Convert to orbital elements*/
  mKOST::elements2StateVector (MU, &elements, &out, SIMD_EPSILON, 1000000);
  mKOST::stateVector2Elements (MU, &out, &elements2, &params);


  printf ("Orbital elements:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements.a, elements.e, elements.i, elements.LaN, elements.LoP, elements.L
         );

  printf ("Orbital elements reversed:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements2.a, elements2.e, elements2.i, elements2.LaN, elements2.LoP, elements2.L
         );
  printf ("Additional parameters:\n"
          "   PeD = %f m\n"
          "   ApD = %f m\n"
          "   MnA = %f\n"
          "   TrA = %f\n"
          "   EcA = %f\n"
          "     T = %f s\n"
          "   AgP = %f\n",
          params.PeD,
          params.ApD,
          params.MnA,
          params.TrA,
          params.EcA,
          params.T,
          params.AgP
         );
  printf ("state:\n"
          "  position: %f, %f, %f\n"
          "  velocity: %f, %f, %f\n",
          out.pos.getX(),
          out.pos.getY(),
          out.pos.getZ(),
          out.vel.getX(),
          out.vel.getY(),
          out.vel.getZ()
          );

  return 0;
}
