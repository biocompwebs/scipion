/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
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

#ifndef _PROG_MONOGENIC_SIGNAL_RES
#define _PROG_MONOGENIC_SIGNAL_RES

#include <iostream>
#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/metadata.h>
//#include <data/xmipp_fft.h>
#include <data/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>


/**@defgroup SSNR resolution_ssnr (Spectral Signal to Noise Ratio)
   @ingroup ReconsLibrary */
//@{
/** SSNR parameters. */

class ProgMonogenicSignalRes : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnUntilt, fnTilt, fnDir, fnVol;

	/** sampling rate*/
	double ang_dist, ang_acc;

	/** number of classes*/
	int N_cls;
public:

    void defineParams();
    void readParams();
    void filterVolume();
    void RieszTransform3Dreal(const MultidimArray<double> &inputVol, std::vector<MultidimArray< std::complex<double> > > &RieszVector);
    void amplitudeMonogenicSignal(const MultidimArray<double> &inputVol,
			  	  	  	  	  	  const std::vector<MultidimArray< std::complex<double> > > &RieszVector,
			  	  	  	  	  	  	  	  	  	  	  MultidimArray<double> &amplitude);
    void passbandfiltervol(const MultidimArray<double> &inputVol, double freq, MultidimArray<double> &filteredVol, double FWHM=0.05);
    void run();

};
//@}
#endif
