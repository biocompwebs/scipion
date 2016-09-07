/***************************************************************************
 * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
 *
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

#include <vector>
#include <sstream>
#include <fstream>
#include <time.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/gpu/gpu.hpp"

#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/normalize.h>
#include <data/xmipp_fftw.h>


using namespace std;
//using namespace cv;
#ifdef GPU
using namespace cv::gpu;
#endif
class ProgOpticalAligment: public XmippProgram
{

public:
    FileName fnMovie, fnOut, fnGain, fnDark;
    FileName fnMovieOut, fnMicOut;
    Image<double> gain, dark;
    ImageGeneric movieStack;
    std::vector< Matrix1D<double> > shiftVector;
    MetaData shiftMD, movie;
    int winSize, gpuDevice, nfirst, nlast;
    int groupSize;
    bool globalShiftCorr;

    /*****************************/
    /** crop corner **/
    /*****************************/
    /** x left top corner **/
    int xLTcorner;
    /** y left top corner **/
    int yLTcorner;
    /** x right down corner **/
    int xDRcorner;
    /** y right down corner **/
    int yDRcorner;

    void defineParams()
    {
        addUsageLine ("Align movies using optical flow");
        addParamsLine("    -i <inMoviewFnName>           : input movie File Name");
        addParamsLine("   [-o <fn=\"out.xmd\">]           : Metadata with the shifts of each frame.");
        addParamsLine("   [--oavg <fn=\"\">]              : Give the name of a micrograph to generate an aligned micrograph");
        addParamsLine("   [--cropULCorner <x=0> <y=0>]    : crop up left corner (unit=px, index starts at 0)");
        addParamsLine("   [--cropDRCorner <x=-1> <y=-1>]  : crop down right corner (unit=px, index starts at 0), -1 -> no crop");
        addParamsLine("   [--frameRange <n0=-1> <nF=-1>]  : First and last frame to align, frame numbers start at 0");
        addParamsLine("   [--bin <s=-1>]               : Binning factor, it may be any floating number");
        addParamsLine("   [--winSize <int=150>]        : window size for optical flow algorithm");
        addParamsLine("   [--groupSize <int=1>]        : the depth of pyramid for optical flow algorithm");
        addParamsLine("   [--outMovie <fn=\"\">]       : save corrected stack");
        addParamsLine("   [--dark <fn=\"\">]           : Dark correction image");
        addParamsLine("   [--gain <fn=\"\">]           : Gain correction image");
#ifdef GPU

        addParamsLine("   [--gpu <int=0>]              : GPU device to be used");
#endif

    }
    void readParams()
    {
        fnMovie = getParam("-i");
        fnOut = getParam("-o");
        fnGain = getParam("--gain");
        fnDark = getParam("--dark");
        fnMicOut = getParam("--oavg");
        fnMovieOut = getParam("--outMovie");
        groupSize = getIntParam("--groupSize");
        nfirst = getIntParam("--frameRange",0);
        nlast = getIntParam("--frameRange",1);
        winSize   = getIntParam("--winSize");
        xLTcorner= getIntParam("--cropULCorner",0);
        yLTcorner= getIntParam("--cropULCorner",1);
        xDRcorner = getIntParam("--cropDRCorner",0);
        yDRcorner = getIntParam("--cropDRCorner",1);

#ifdef GPU

        gpuDevice = getIntParam("--gpu");
#endif

    }
    void run()
    {
        main2();
    }

    // Save a matrix which is generated by OpenCV
    int saveMat( const string& filename, const cv::Mat& M)
    {
        if (M.empty())
        {
            return 0;
        }
        ofstream out(filename.c_str(), ios::out|ios::binary);
        if (!out)
            return 0;

        int cols = M.cols;
        int rows = M.rows;
        int chan = M.channels();
        int eSiz = (M.dataend-M.datastart)/(cols*rows*chan);

        // Write header
        out.write((char*)&cols,sizeof(cols));
        out.write((char*)&rows,sizeof(rows));
        out.write((char*)&chan,sizeof(chan));
        out.write((char*)&eSiz,sizeof(eSiz));

        // Write data.
        if (M.isContinuous())
        {
            out.write((char *)M.data,cols*rows*chan*eSiz);
        }
        else
        {
            return 0;
        }
        out.close();
        return 1;
    }

    // Load a matrix which is generated by saveMat
    int readMat( const string& filename, cv::Mat& M)
    {
        ifstream in(filename.c_str(), ios::in|ios::binary);
        if (!in)
        {
            //M = NULL_MATRIX;
            return 0;
        }
        int cols;
        int rows;
        int chan;
        int eSiz;

        // Read header
        in.read((char*)&cols,sizeof(cols));
        in.read((char*)&rows,sizeof(rows));
        in.read((char*)&chan,sizeof(chan));
        in.read((char*)&eSiz,sizeof(eSiz));

        // Determine type of the matrix
        int type = 0;
        switch (eSiz)
        {
        case sizeof(char):
                        type = CV_8UC(chan);
            break;
        case sizeof(float):
                        type = CV_32FC(chan);
            break;
        case sizeof(double):
                        type = CV_64FC(chan);
            break;
        }

        // Alocate Matrix.
        M = cv::Mat(rows,cols,type,cv::Scalar(1));

        // Read data.
        if (M.isContinuous())
{
            in.read((char *)M.data,cols*rows*chan*eSiz);
        }
        else
        {
            return 0;
        }
        in.close();
        return 1;
    }

    // Converts a XMIPP MultidimArray to OpenCV matrix
    void xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat)
    {
        int h = YSIZE(xmippArray);
        int w = XSIZE(xmippArray);
        opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
        opencvMat.at<float>(i,j) = DIRECT_A2D_ELEM(xmippArray,i,j);
    }

    // Converts an OpenCV matrix to XMIPP MultidimArray
    void opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray)
    {
        int h = opencvMat.rows;
        int w = opencvMat.cols;
        xmippArray.initZeros(h, w);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
        DIRECT_A2D_ELEM(xmippArray,i,j) = opencvMat.at<float>(i,j);
    }

    // Converts an OpenCV float matrix to an OpenCV Uint8 matrix
    void convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat)
    {
        double min,max;
        cv::minMaxLoc(opencvDoubleMat, &min, &max);
        opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));
    }

    // correct for dark and gain reference if available
    void correctDarkGainImage(Image<double> & croppedFrame)
    {
        if (XSIZE(dark())>0)
            croppedFrame()-=dark();
        if (XSIZE(gain())>0)
            croppedFrame()*=gain();
    }

    // Computes the average of a number of frames in movies
    void computeAvg(int begin, int end, MultidimArray<double> &avgImg)
    {
        //ImageGeneric movieStack;
        Image<double> frame;
        FileName fnFrame;
        int N=end-begin+1;

        for (size_t i=begin;i<=end;i++)
        {
            //movieStack.readMapped(fnMovie,i);
            //movieStack().getImage(frameImage);
            movie.getValue(MDL_IMAGE, fnFrame, i);
            frame.read(fnFrame);
            if (!globalShiftCorr)
            {
                if (yDRcorner!=-1)
                    frame().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
                correctDarkGainImage(frame);
            }
            if (i==begin)
                avgImg.initZeros(YSIZE(frame()), XSIZE(frame()));
            avgImg+=frame();
        }
        avgImg/=double(N);
        frame.clear();
    }

    void std_dev2(const cv::Mat planes[], const cv::Mat &flowx, const cv::Mat &flowy, Matrix1D<double> &meanStdDev)
    {
        double sumX=0, sumY=0;
        double absSumX=0, absSumY=0;
        double sqSumX=0, sqSumY=0;
        double valSubtract;
        int h=flowx.rows;
        int w=flowx.cols;
        int n=h*w;
        for(int i=0;i<h;i++)
            for(int j=0;j<w;j++)
            {
                valSubtract=planes[0].at<float>(i,j)-flowx.at<float>(i,j);
                sumX+=valSubtract;
                absSumX+=abs(valSubtract);
                sqSumX+=valSubtract*valSubtract;
                valSubtract=planes[1].at<float>(i,j)-flowy.at<float>(i,j);
                sumY+=valSubtract;
                absSumY+=abs(valSubtract);
                sqSumY+=valSubtract*valSubtract;
            }
        double avgX=sumX/double(n);
        double avgY=sumY/double(n);
        meanStdDev(0)=absSumX/double(n);
        meanStdDev(1)=sqrt(sqSumX/double(n)-avgX*avgX);
        meanStdDev(2)=absSumY/double(n);
        meanStdDev(3)=sqrt(sqSumY/double(n)-avgY*avgY);
    }

    int main2()
    {

        MultidimArray<double> mappedImg;
        MultidimArray<double> outputMovie;
        Matrix1D<double> meanStdev;
        Image<double> II, preImg, avgCurr;
        MetaData MD; // To save plot information
        FileName flowFileName;
        FileName flowXFileName, flowYFileName;
        FileName tmpFileName, fnBaseOut;
        FileName fnFrame, fnmovieRoot;
        size_t Xdim, Ydim, Zdim, Ndim;
        ArrayDim aDim;

        // For measuring times (both for whole process and for each level of the pyramid)
        clock_t tStart, tStart2;

#ifdef GPU
        // Matrix that we required in GPU part
        GpuMat d_flowx, d_flowy, d_dest;
        GpuMat d_avgcurr, d_preimg;
#endif

        // Matrix required by Opencv
        cv::Mat flow, dest, flowx, flowy;
        cv::Mat flowxPre, flowyPre;
        cv::Mat avgcurr, avgstep, preimg, preimg8, avgcurr8;
        cv::Mat planes[]={flowxPre, flowyPre};

        int imagenum, cnt=2, div=0, flowCounter;
        int levelNum, levelCounter=1;

        fnmovieRoot=fnMovie.getDir();
        tmpFileName=fnmovieRoot+"tmpMovie.stk";
        fnBaseOut=fnOut.removeDirectories();
        std::cerr<<fnOut<<std::endl;
        std::string extension=fnMovie.getExtension();
        //if input is an stack create a metadata.
        if (fnMovie.isMetaData())
        {
            movie.read(fnMovie);
            getImageSize(movie, Xdim, Ydim, Zdim, Ndim);
        }
        else
        {
            ImageGeneric movieStack;
            movieStack.read(fnMovie,HEADER);
            movieStack.getDimensions(Xdim,Ydim,Zdim,Ndim);
            if (fnMovie.getExtension()=="mrc" and Ndim ==1)
                Ndim = Zdim;
            size_t id;
            FileName fn;
            for (size_t i=0;i<Ndim;i++)
            {
                id = movie.addObject();
                fn.compose(i+FIRST_IMAGE,fnMovie);
                movie.setValue(MDL_IMAGE, fn, id);
            }
        }
        if (yDRcorner!=-1)
        {
            Xdim = xDRcorner - xLTcorner +1 ;
            Ydim = yDRcorner - yLTcorner +1 ;
        }
        if (Zdim!=1)
            REPORT_ERROR(ERR_ARG_INCORRECT,"This program is meant to align 2D frames, not 3D");
        imagenum=Ndim;
        if (fnDark!="")
        {
            dark.read(fnDark);
            if (yDRcorner!=-1)
                dark().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
        }
        if (fnGain!="")
        {
            gain.read(fnGain);
            if (yDRcorner!=-1)
                gain().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
            gain()=1.0/gain();
            double avg=gain().computeAvg();
            if (isinf(avg) || isnan(avg))
                REPORT_ERROR(ERR_ARG_INCORRECT,"The input gain image is incorrect, its inverse produces infinite or nan");
        }
        meanStdev.initZeros(4);
        //avgCurr.initZeros(Ydim, Xdim);
        flowxPre=cv::Mat::zeros(Ydim, Xdim,CV_32FC1);
        flowyPre=cv::Mat::zeros(Ydim, Xdim,CV_32FC1);
#ifdef GPU

        // Object for optical flow
        FarnebackOpticalFlow d_calc;
        setDevice(gpuDevice);

        // Initialize the parameters for optical flow structure
        d_calc.numLevels=6;
        d_calc.pyrScale=0.5;
        d_calc.fastPyramids=true;
        d_calc.winSize=winSize;
        d_calc.numIters=1;
        d_calc.polyN=5;
        d_calc.polySigma=1.1;
        d_calc.flags=0;
#endif

        // Read shifts if we have shift in the metadata
        if (movie.containsLabel(MDL_SHIFT_X))
        {
            globalShiftCorr=true;
            Matrix1D<double> shiftMatrix(2);
            shiftVector.reserve(imagenum);
            Image<double> frameImage;
            FOR_ALL_OBJECTS_IN_METADATA(movie)
            {
                movie.getValue(MDL_IMAGE, fnFrame, __iter.objId);
                shiftMD.getValue(MDL_SHIFT_X, XX(shiftMatrix), __iter.objId);
                shiftMD.getValue(MDL_SHIFT_Y, YY(shiftMatrix), __iter.objId);
                shiftVector.push_back(shiftMatrix);
                frameImage.read(fnFrame);
                if (yDRcorner!=-1)
                    frameImage().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
                correctDarkGainImage(frameImage);
                translate(LINEAR, preImg(), frameImage(), shiftMatrix, WRAP);
                preImg.write(tmpFileName, ALL_IMAGES, true, WRITE_APPEND);
            }
            tmpFileName=fnMovie;
            fnMovie=fnmovieRoot+"tmpMovie.stk";
        }
        else
            globalShiftCorr=false;
        tStart2=clock();
        // the frame index begin from 1 in this program
        if (nfirst<0)
            nfirst=1;
        else
            nfirst++;
        if (nlast<0)
            nlast=imagenum;
        else
            nlast++;
        imagenum=nlast-nfirst+1;
        // Initialize the stack for the output movie
        if (fnMovieOut!="")
            outputMovie.initZeros(imagenum, 1, Ydim, Xdim);
        levelNum=sqrt(double(imagenum));
        computeAvg(nfirst, nlast, avgCurr());
        xmipp2Opencv(avgCurr(), avgcurr);
        cout<<"Frames "<<nfirst<<" to "<<nlast<<" under processing ..."<<std::endl;
        while (div!=groupSize)
        {
            div=int(imagenum/cnt);
            // avgStep to hold the sum of aligned frames of each group at each step
            avgstep=cv::Mat::zeros(Ydim, Xdim, CV_32FC1);

            cout<<"Level "<<levelCounter<<"/"<<levelNum<<" of the pyramid is under processing"<<std::endl;
            // Compute time for each level
            tStart = clock();

            // Check if we are in the final step
            if (div==1)
                cnt=imagenum;
            flowCounter=1;
            for (int i=0;i<cnt;i++)
            {
                //Just compute the average in the last step
                if (div==1)
                {
                    movie.getValue(MDL_IMAGE, fnFrame, nfirst+i);
                    preImg.read(fnFrame);
                    if (!globalShiftCorr)
                    {
                        if (yDRcorner!=-1)
                            preImg().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
                        correctDarkGainImage(preImg);
                    }
                }
                else
                {
                    if (i==cnt-1)
                        computeAvg(i*div+nfirst, nlast, preImg());
                    else
                        computeAvg(i*div+nfirst, (i+1)*div+nfirst-1, preImg());
                }
                xmipp2Opencv(preImg(), preimg);
                // Note: We should convert Xmipp image to OpenCV 8 bits image
                convert2Uint8(avgcurr,avgcurr8);
                convert2Uint8(preimg,preimg8);
#ifdef GPU

                d_avgcurr.upload(avgcurr8);
                d_preimg.upload(preimg8);

                if (cnt==2)
                    d_calc(d_avgcurr, d_preimg, d_flowx, d_flowy);
                else
                {
                    flowXFileName=fnmovieRoot;
                    flowYFileName=fnmovieRoot;
                    flowXFileName+=fnBaseOut.removeLastExtension()+formatString("flowx%d%d.txt",div*2,flowCounter);
                    flowYFileName+=fnBaseOut.removeLastExtension()+formatString("flowy%d%d.txt",div*2,flowCounter);
                    readMat(flowXFileName.c_str(), flowx);
                    readMat(flowYFileName.c_str(), flowy);
                    d_flowx.upload(flowx);
                    d_flowy.upload(flowy);
                    d_calc.flags=cv::OPTFLOW_USE_INITIAL_FLOW;
                    d_calc(d_avgcurr, d_preimg, d_flowx, d_flowy);
                }

                d_flowx.download(planes[0]);
                d_flowy.download(planes[1]);
                d_avgcurr.release();
                d_preimg.release();
                d_flowx.release();
                d_flowy.release();
#else
                // Check if we should use the flows from the previous steps
                if (cnt==2)
                    calcOpticalFlowFarneback(avgcurr8, preimg8, flow, 0.5, 6, winSize, 1, 5, 1.1, 0);
                else
                {
                    flowFileName=fnmovieRoot;
                    flowFileName+=fnBaseOut.removeLastExtension()+formatString("flow%d%d.txt",div*2,flowCounter);
                    readMat(flowFileName.c_str(), flow);
                    calcOpticalFlowFarneback(avgcurr8, preimg8, flow, 0.5, 6, winSize, 1, 5, 1.1, cv::OPTFLOW_USE_INITIAL_FLOW);
                }
                split(flow, planes);

#endif
                // Save the flows if we are in the last step
                if (div==groupSize)
                {
                    if (i > 0)
                    {
                        std_dev2(planes,flowxPre,flowyPre,meanStdev);
                        size_t id=MD.addObject();
                        MD.setValue(MDL_OPTICALFLOW_MEANX, double(meanStdev(0)), id);
                        MD.setValue(MDL_OPTICALFLOW_MEANY, double(meanStdev(2)), id);
                        MD.setValue(MDL_OPTICALFLOW_STDX, double(meanStdev(1)), id);
                        MD.setValue(MDL_OPTICALFLOW_STDY, double(meanStdev(3)), id);
                        MD.write(fnOut, MD_APPEND);
                    }
                    planes[0].copyTo(flowxPre);
                    planes[1].copyTo(flowyPre);
                }
                else
                {
#ifdef GPU
                    flowXFileName=fnmovieRoot;
                    flowYFileName=fnmovieRoot;
                    flowXFileName+=fnBaseOut.removeLastExtension()+formatString("flowx%d%d.txt",div,i+1);
                    flowYFileName+=fnBaseOut.removeLastExtension()+formatString("flowy%d%d.txt",div,i+1);
                    saveMat(flowXFileName.c_str(), planes[0]);
                    saveMat(flowYFileName.c_str(), planes[1]);
#else

                    flowFileName=fnmovieRoot;
                    flowFileName+=fnBaseOut.removeLastExtension()+formatString("flow%d%d.txt",div,i+1);
                    saveMat(flowFileName.c_str(), flow);
#endif

                    if ((i+1)%2==0)
                        flowCounter++;
                }
                for( int row = 0; row < planes[0].rows; row++ )
                    for( int col = 0; col < planes[0].cols; col++ )
                    {
                        if (div==1 && globalShiftCorr)
                        {
                            planes[0].at<float>(row,col) += (-1.0*XX(shiftVector[i]))+(float)col;
                            planes[1].at<float>(row,col) += (-1.0*YY(shiftVector[i]))+(float)row;
                        }
                        else
                        {
                            planes[0].at<float>(row,col) += (float)col;
                            planes[1].at<float>(row,col) += (float)row;
                        }
                    }
                if (div==1 && globalShiftCorr)
                {
                    //movieStack.readMapped(tmpFileName,nfirst+i);
                    //movieStack().getImage(preImg);
                    movie.getValue(MDL_IMAGE, fnFrame, nfirst+i);
                    preImg.read(fnFrame);
                    xmipp2Opencv(preImg(), preimg);
                }
                cv::remap(preimg, dest, planes[0], planes[1], cv::INTER_CUBIC);
                if (div==1 &&  fnMovieOut!="")
                {
                    mappedImg.aliasImageInStack(outputMovie, i);
                    opencv2Xmipp(dest, mappedImg);
                }
                avgstep+=dest;
            }
            avgcurr=avgstep/cnt;
            cout<<"Processing level "<<levelCounter<<"/"<<levelNum<<" has been finished"<<std::endl;
            printf("Processing time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
            cnt*=2;
            levelCounter++;
        }
        opencv2Xmipp(avgcurr, avgCurr());
        avgCurr.write(fnMicOut);
        printf("Total Processing time: %.2fs\n", (double)(clock() - tStart2)/CLOCKS_PER_SEC);
        if (fnMovieOut!="")
        {
            II()=outputMovie;
            II.write(fnMovieOut);
        }

        // Release the memory
        avgstep.release();
        preimg.release();
        avgcurr8.release();
        preimg8.release();
        flow.release();
        planes[0].release();
        planes[1].release();
        flowxPre.release();
        flowyPre.release();
        movieStack.clear();
        preImg.clear();
        avgCurr.clear();
        II.clear();
        return 0;
    }
};

int main(int argc, char *argv[])
{
    ProgOpticalAligment prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
