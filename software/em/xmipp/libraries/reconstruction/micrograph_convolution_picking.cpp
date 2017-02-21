/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2017)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "micrograph_convolution_picking.h"
#include <stdlib.h>
#include <stdio.h>

//#define DEBUG

void ProgMicrographConvPicking::readParams()
{
	fnmic = getParam("--micrograph");
	fnOut = getParam("-o");
}


void ProgMicrographConvPicking::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --micrograph <img_file=\"\">   : Input micrograph");
	addParamsLine("  [-o <output=\"blurred.xmp\">]: Local resolution volume (in Angstroms)");
}

void ProgMicrographConvPicking::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;
}


void ProgMicrographConvPicking::run()
{
	std::cout << "Starting..." << std::endl;
	Image<double> micImg;
	Image<int> outputImg;
	size_t Xdim, Ydim, Zdim, Ndim;

	//Data reading
	micImg.read(fnmic);
	micImg.getDimensions(Xdim, Ydim, Zdim, Ndim);

	if (Xdim % 2 == 0){
		Xdim = Xdim + 1;
		std::cout << "Xdim es par, annado uno\n" << std::endl;
	}
	if (Ydim % 2 == 0)
	{
		Ydim = Ydim + 1;
		std::cout << "Ydim es par, annado uno\n" << std::endl;
	}

	//Defining a disk for the convolution
	MultidimArray<int> circKernel;
	circKernel.initZeros(Ndim, Zdim, Ydim*0.5, Xdim*0.5);

	circKernel.printShape();

	int Xorig = Xdim*0.5*0.5+1;
	int Yorig = Ydim*0.5*0.5+1;

	std::cout << "centro x = " << Xorig << std::endl;
	std::cout << "centro y = " << Yorig << std::endl;

	size_t NpixelsRadius = 5;
	A2D_ELEM(circKernel, Xorig, Yorig) = 1;

	std::cout << A2D_ELEM(circKernel, Xorig, Yorig) << std::endl;

	for (size_t i=0; i< NpixelsRadius; i++)
	{
		for (size_t j=0; j< NpixelsRadius; j++)
		{
			if (i*i+j*j<NpixelsRadius*NpixelsRadius)
			{
				A2D_ELEM(circKernel, Xorig+i, Yorig+j) = 1;
				A2D_ELEM(circKernel, Xorig+i, Yorig-j) = 1;
				A2D_ELEM(circKernel, Xorig-i, Yorig+j) = 1;
				A2D_ELEM(circKernel, Xorig-i, Yorig-j) = 1;
			}
		}
	}
	outputImg = circKernel;
	outputImg.write(fnOut);
}

