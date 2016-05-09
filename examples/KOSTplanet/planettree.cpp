/***************************************************************************
 *   Copyright (C) 2008 by C J Plooy                                       *
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

#include <vector>

#include "kost.h"

#include "planettree.h"

struct SPlanetTreeItem
{
	OBJHANDLE object;
	int parent;
};

std::vector<SPlanetTreeItem> _thePlanetTree;
int _theCentralStar = -1;

void initPlanetTree()
{
	unsigned int numBodies = oapiGetGbodyCount();
	_thePlanetTree.resize(numBodies);

	//Initialize array
	unsigned int i = 0, j = 0;
	for(i=0; i < numBodies; i++)
	{
		_thePlanetTree[i].object = oapiGetGbodyByIndex(i);
		_thePlanetTree[i].parent = -1;
	}

	//Find default parent
	_theCentralStar = -1;
	for(i=0; i < numBodies; i++)
		if(oapiGetObjectType(_thePlanetTree[i].object) == OBJTP_STAR)
		{
			_theCentralStar = i;
			break;
		}

	//Find parent link
	for(i=0; i < numBodies; i++)
	{
		_thePlanetTree[i].parent = _theCentralStar;
		if(i == _theCentralStar) continue;

		OBJHANDLE obj_i = _thePlanetTree[i].object;
		double m = oapiGetMass(obj_i);

		double maxF = 0.0;

		for(j=0; j < numBodies; j++)
		{
			if(j == i || j == _theCentralStar) continue;

			OBJHANDLE obj_j = _thePlanetTree[j].object;

			//Mass test
			double thisMass = oapiGetMass(obj_j);
			if(thisMass < m) continue;

			//Eccentricity test
			kostVector3 pos, vel;
			oapiGetRelativePos(obj_j, obj_i, &pos);
			oapiGetRelativeVel(obj_j, obj_i, &vel);

			double mu = thisMass * KOST_GRAVITATIONAL_CONSTANT;
			double r  = abs(pos);
			double v2 = abs2(vel);
			double e  = abs(pos * (v2 - mu/r) - vel * dotProduct(pos,vel)) / mu;
			if(e > 1.0) continue;

			//Largest gravity selection
			double F = thisMass / (r*r);				
			if(F < maxF) continue;
			maxF = F;

			//Congratulations: you made it ;)
			_thePlanetTree[i].parent = j;
		}
	}

}

OBJHANDLE getCentralBody(OBJHANDLE satellite)
{
	for(unsigned int i=0; i < _thePlanetTree.size(); i++)
		if(_thePlanetTree[i].object == satellite)
		{
			if(_thePlanetTree[i].parent < 0) return NULL;
			return _thePlanetTree[_thePlanetTree[i].parent].object;
		}

	if(_theCentralStar < 0) return NULL;
	return _thePlanetTree[_theCentralStar].object;
}

