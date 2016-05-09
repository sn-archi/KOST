/***************************************************************************
 *   Copyright (C) 2009 by C J Plooy                                       *
 *   cornwarecjp@lavabit.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#define ORBITER_MODULE
#include "OrbiterAPI.h"
#include "CelbodyAPI.h"

#include <cstdio>
#include <cstdlib>

#include "kost.h"
#include "planettree.h"

class KOSTPlanet: public CELBODY
{
public:
	KOSTPlanet();
	bool bEphemeris() const;
	void clbkInit(FILEHANDLE cfg);
	int clbkEphemeris(double mjd, int req, double *ret);
	int clbkFastEphemeris(double simt, int req, double *ret);

	kostElements m_Elements;
	kostReal m_EpochMJD;

	kostReal m_PreviousEccentricAnomaly;
	kostReal m_dt_start;

	kostReal m_Mass, m_ParentMass, m_Mu;

	OBJHANDLE m_ParentObject;

	void detectParentInfo();
};

DLLCLBK void InitModule(HINSTANCE hModule)
{
	// module initialisation
}

DLLCLBK void ExitModule(HINSTANCE hModule)
{
	// module cleanup
}

DLLCLBK CELBODY *InitInstance(OBJHANDLE hBody)
{
	// instance initialisation
	return new KOSTPlanet;
}

DLLCLBK void ExitInstance(CELBODY *body)
{
	// instance cleanup
	delete (KOSTPlanet *)body;
}


KOSTPlanet::KOSTPlanet()
{
	m_Elements.a = 1.0;
	m_Elements.e = 0.0;
	m_Elements.i = 0.0;
	m_Elements.theta = 0.0;
	m_Elements.omegab = 0.0;
	m_Elements.L = 0.0;

	m_Mass = 1.0;
	m_ParentMass = 1.0;
	m_Mu = 0.0;
	m_ParentObject = NULL;

	m_PreviousEccentricAnomaly = 0.0;
	m_dt_start = 0.0;
}

void KOSTPlanet::clbkInit(FILEHANDLE cfg)
{
	double epoch = 2000.0;

	// read parameters from config file:
	char *line = NULL;
	while(oapiReadScenario_nextline(cfg, line))
	{
		unsigned int len = strlen(line);

		char *pos = strstr(line, "=");
		if(pos == NULL) continue;

		//Make left nul-terminated:
		char *lhs = line;
		{
			char *endleft = pos;

			//trim:
			while(endleft > line && *(endleft-1) == ' ')
				endleft--;

			*endleft = '\0';
		}

		if(*lhs == '\0') continue; //empty LHS

		//Trim right:
		char *rhs = pos+1;
		while(*rhs == ' ') rhs++;

		unsigned int lhslen = strlen(lhs);

		if(lhslen == 5 && strcmp(lhs, "Epoch") == 0)
		{
			epoch = strtod(rhs, NULL);
		}
		else if(lhslen == 13 && strcmp(lhs, "SemiMajorAxis") == 0)
		{
			m_Elements.a = strtod(rhs, NULL);
		}
		else if(lhslen == 12 && strcmp(lhs, "Eccentricity") == 0)
		{
			m_Elements.e = strtod(rhs, NULL);
		}
		else if(lhslen == 11 && strcmp(lhs, "Inclination") == 0)
		{
			m_Elements.i = strtod(rhs, NULL);
		}
		else if(lhslen == 11 && strcmp(lhs, "LongAscNode") == 0)
		{
			m_Elements.theta = strtod(rhs, NULL);
		}
		else if(lhslen == 14 && strcmp(lhs, "LongPerihelion") == 0)
		{
			m_Elements.omegab = strtod(rhs, NULL);
		}
		else if(lhslen == 13 && strcmp(lhs, "MeanLongitude") == 0)
		{
			m_Elements.L = strtod(rhs, NULL);
		}
		else if(lhslen == 4 && strcmp(lhs, "Mass") == 0)
		{
			m_Mass = strtod(rhs, NULL);
		}
	}


	//Convert to MJD:
	m_EpochMJD = (epoch - 1858.87885010) * (KOST_YEAR/KOST_DAY);

	//dt between epoch and sim start:
	m_dt_start = (oapiTime2MJD(0.0) - m_EpochMJD) * KOST_DAY;
}

bool KOSTPlanet::bEphemeris() const
{
	// class supports ephemeris calculation
	return true;
}

/*
EPHEM_TRUEPOS (true body position)
EPHEM_TRUEVEL (true body velocity)
EPHEM_BARYPOS (barycentric position)
EPHEM_BARYVEL (barycentric velocity)

EPHEM_TRUEISBARY

ret[0-2]: true position (if requested)
ret[3-5]: true velocity (if requested)
ret[6-8]: barycentric position (if requested)
ret[9-11]: barycentric velocity (if requested)
*/

int KOSTPlanet::clbkEphemeris(double mjd, int req, double *ret)
{
	//Hopefully this will only be needed once:
	if(m_ParentObject == NULL)
		detectParentInfo();

	kostReal dt = (mjd - m_EpochMJD) * KOST_DAY;

	kostStateVector state;

	kostElements2StateVectorAtTime(
		m_Mu, &m_Elements, &state, dt, 1e-10, 100, 0.0, 0.0);

	//Calculate back to Orbiter coordinate system (swapping y and z):
	ret[0] = state.pos.x;
	ret[1] = state.pos.z;
	ret[2] = state.pos.y;
	ret[3] = state.vel.x;
	ret[4] = state.vel.z;
	ret[5] = state.vel.y;
	for(unsigned int i=0; i < 6; i++)
		ret[i+6] = ret[i];

	return EPHEM_TRUEPOS | EPHEM_TRUEVEL | EPHEM_BARYPOS | EPHEM_BARYVEL | EPHEM_BARYISTRUE;
}

int KOSTPlanet::clbkFastEphemeris(double simt, int req, double *ret)
{
	//return clbkEphemeris(oapiTime2MJD(simt), req, ret);

	//Hopefully this will only be needed once:
	if(m_ParentObject == NULL)
		detectParentInfo();

	kostReal dt = m_dt_start + simt;

	kostReal meanAnomaly = kostGetMeanAnomalyAtTime(
		m_Mu, &m_Elements, dt);

	kostReal eccAnomaly = m_PreviousEccentricAnomaly;
	if(eccAnomaly == 0.0) eccAnomaly = meanAnomaly;

	//Update eccentric anomaly estimate:
	int num = kostGetEccentricAnomaly(
		&m_Elements, &eccAnomaly, meanAnomaly, eccAnomaly, 1e-10, 100);

	//sprintf (oapiDebugString(), "%d iterations", num);
	m_PreviousEccentricAnomaly = eccAnomaly;

	kostReal trueAnomaly = kostGetTrueAnomaly2(
		m_Mu, &m_Elements, eccAnomaly);

	kostStateVector state;
	kostElements2StateVector2(
		m_Mu, &m_Elements, &state, trueAnomaly);

	//Calculate back to Orbiter coordinate system (swapping y and z):
	ret[0] = state.pos.x;
	ret[1] = state.pos.z;
	ret[2] = state.pos.y;
	ret[3] = state.vel.x;
	ret[4] = state.vel.z;
	ret[5] = state.vel.y;
	for(unsigned int i=0; i < 6; i++)
		ret[i+6] = ret[i];

	return EPHEM_TRUEPOS | EPHEM_TRUEVEL | EPHEM_BARYPOS | EPHEM_BARYVEL | EPHEM_BARYISTRUE;
}

void KOSTPlanet::detectParentInfo()
{
	initPlanetTree();

	OBJHANDLE thisObj = NULL;
	unsigned int numBodies = oapiGetGbodyCount();
	for(unsigned int i=0; i < numBodies; i++)
	{
		OBJHANDLE obj = oapiGetGbodyByIndex(i);
		CELBODY *celBody = oapiGetCelbodyInterface(obj);
		if(celBody == this)
		{
			thisObj = obj;
			break;
		}
	}

	m_ParentObject = getCentralBody(thisObj);

	if(m_ParentObject != NULL)
	{
		m_ParentMass = oapiGetMass(m_ParentObject);
		m_Mu = KOST_GRAVITATIONAL_CONSTANT * (m_ParentMass + m_Mass);
	}
}

