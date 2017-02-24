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

#ifndef _PROG_CONVOLUTION_PICKING
#define _PROG_CONVOLUTION_PICKING

#include <iostream>
#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/metadata.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>
#include "fourier_filter.h"
#include <data/filters.h>
#include <string>
#include "symmetrize.h"

/**@defgroup Convolution Picking
   @ingroup ReconsLibrary */
//@{
/** SSNR parameters. */

class ProgMicrographConvPicking : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnmic, fnOut;

	/** sampling rate, minimum resolution, and maximum resolution */
	double radius, minRes, maxRes, R;

	/** Is the volume previously masked?*/
	//int R;

	/** Step in digital frequency */
	double N_freq, trimBound, significance;

	/** The search for resolutions is linear or inverse**/
	bool linearchk, exactres;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();
    void run();
};
//@}
#endif
