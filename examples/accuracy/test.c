#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kost.h"

/*Data about central body (earth)*/
#define R 6378100.0
#define M 5.9736e24
#define mu (KOST_GRAVITATIONAL_CONSTANT * M)

kostStateVector sv_maxVerror, sv_maxRerror;
kostReal maxVerror= -1.0, maxRerror = -1.0;

void testState(const kostStateVector *sv)
{
	kostElements elements;
	kostStateVector out;
	kostVector3 diff;
	kostReal error;

	/*Convert to orbital elements*/
	kostStateVector2Elements(mu, sv, &elements, NULL);

	/*Convert back to state vector*/
	kostElements2StateVector(mu, &elements, &out, KOST_VERYSMALL, 1000000);

	diff = kostSubvv(&(sv->pos), &(out.pos));
	error = kostAbsv(&diff) / kostAbsv(&(sv->pos));
	if(error > maxRerror)
	{
		maxRerror = error;
		sv_maxRerror = *sv;

		printf("New pos error max:\n"
			"     pos = %e, %e, %e\n"
			"     vel = %e, %e, %e\n"
			"     error = %e\n",
			sv_maxRerror.pos.x, sv_maxRerror.pos.y, sv_maxRerror.pos.z,
			sv_maxRerror.vel.x, sv_maxRerror.vel.y, sv_maxRerror.vel.z,
			maxRerror
		);
	}

	diff = kostSubvv(&(sv->vel), &(out.vel));
	error = kostAbsv(&diff) / kostAbsv(&(sv->vel));
	if(error > maxVerror)
	{
		maxVerror = error;
		sv_maxVerror = *sv;

		printf("New vel error max:\n"
			"     pos = %e, %e, %e\n"
			"     vel = %e, %e, %e\n"
			"     error = %e\n",
			sv_maxVerror.pos.x, sv_maxVerror.pos.y, sv_maxVerror.pos.z,
			sv_maxVerror.vel.x, sv_maxVerror.vel.y, sv_maxVerror.vel.z,
			maxVerror
		);
	}
}

int main(int argc, char *argv[])
{
	int rx, ry, rz, vx, vy, vz;
	kostStateVector sv;
	int rmin = 5, rmax = 12, vmin = 0, vmax = 4;

	/*Arbitrary 6D positions*/
	for(rx=rmin; rx <= rmax; rx++)
	for(ry=rmin; ry <= rmax; ry++)
	for(rz=rmin; rz <= rmax; rz++)
	for(vx=vmin; vx <= vmax; vx++)
	for(vy=vmin; vy <= vmax; vy++)
	for(vz=vmin; vz <= vmax; vz++)
	{
		double rxf = pow(10.0, rx);
		double ryf = pow(10.0, ry);
		double rzf = pow(10.0, rz);
		double vxf = pow(10.0, vx);
		double vyf = pow(10.0, vy);
		double vzf = pow(10.0, vz);

		int srx, sry, srz, svx, svy, svz;
		for(srx=-1; srx<=1; srx++)
		for(sry=-1; sry<=1; sry++)
		for(srz=-1; srz<=1; srz++)
		for(svx=-1; svx<=1; svx++)
		for(svy=-1; svy<=1; svy++)
		for(svz=-1; svz<=1; svz++)
		{
			if(srx == sry == srz == 0)
				continue;
			if(svx == svy == svz == 0)
				continue;

			sv.pos = kostConstructv(srx*rxf, sry*ryf, srz*rzf);
			sv.vel = kostConstructv(svx*vxf, svy*vyf, svz*vzf);
			testState(&sv);
		}
	}

	printf("maxRerror = %e\n", maxRerror);
	printf("maxVerror = %e\n", maxVerror);

	return 0;
}
