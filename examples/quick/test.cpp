#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/kost.h"

int main (int argc, char* argv[])
{
  /*Data about central body (earth)*/
  double R = 6378100.0;
  double mu = 3.986004418e14;

  mKOST::sStateVector out;
  mKOST::sElements elements, elements2;
  mKOST::sOrbitParam params;

  std::cout << std::endl << std::endl << std::endl;

  /*Alternative test, with an initial param list at t0*/
  elements.a = 6139129.926482;
  elements.e = 3.892572e-01;
  elements.i = M_PI / 4;
  elements.omegab = M_PI / 2;
  elements.theta = 0.0;
  elements.L = M_PI / 2;

  /*Convert to orbital elements*/
  mKOST::elements2StateVector (mu, &elements, &out, KOST_VERYSMALL, 1000000);
  mKOST::stateVector2Elements (mu, &out, &elements2, &params);


  printf ("Orbital elements:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements.a, elements.e, elements.i, elements.theta, elements.omegab, elements.L
         );

  printf ("Orbital elements reversed:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements2.a, elements2.e, elements2.i, elements2.theta, elements2.omegab, elements2.L
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
