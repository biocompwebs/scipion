#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Wrapper for relion 1.2
#
#   Author: Roberto Marabini
#

from os.path import join, exists
from os import remove, rename
from protlib_base import protocolMain
from config_protocols import protDict
from xmipp import *
from xmipp import MetaDataInfo
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString, \
                          runJob, runChimera, which, runChimeraClient
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, cyanStr
from protlib_gui_ext import showWarning, showTable, showError
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, \
                               replaceBasenameExt
from protocol_ml2d import lastIteration
from protlib_filesystem import createLink

from protocol_relion_base import ProtRelionBase, runNormalizeRelion, convertImagesMd, renameOutput, \
                                 convertRelionBinaryData, convertRelionMetadata, getIteration

class ProtRelionClassifier(ProtRelionBase):
    def __init__(self, scriptname, project):
        ProtRelionBase.__init__(self, protDict.relion_classify.name, scriptname, project)
        self.Import = 'from protocol_relion_classify import *'
        self.relionType='classify'
        if self.DoContinue:
            self.PrevRun=self.project.getProtocolFromFile(self.optimiserFileName)
            #if optimizer has not been properly selected this will 
            #fail, let us go ahead and handle the situation in verify
            try:
                self.inputProperty('NumberOfClasses')
                ########################### self.inputProperty('SamplingRate')
                self.inputProperty('MaskRadiusA')
                self.inputProperty('RegularisationParamT')
                #self.lastIterationPrecRun=self.PrevRun.lastIter()
    
                #self.NumberOfClasses      = self.PrevRun.NumberOfClasses
                #self.SamplingRate         = self.PrevRun.SamplingRate
                #self.MaskDiameterA        = self.PrevRun.MaskDiameterA
                #self.RegularisationParamT = self.PrevRun.RegularisationParamT
            except:
                print "Can not access the parameters from the original relion run"

    def summary(self):
        lines = ProtRelionBase.summary(self, self.NumberOfIterations)
        lastIteration=self.lastIter()
        lines += ['Performed <%d/%d> iterations ' % (lastIteration,self.NumberOfIterations)]
        lines += ['Number of classes = <%d>' % self.NumberOfClasses]
        return lines
    
    
    def validate(self):
        errors = ProtRelionBase.validate(self)
        return errors 
    
    def defineSteps(self): 
        ProtRelionBase.defineSteps(self)
        # launch relion program
        if self.DoContinue:
            self.insertRelionClassifyContinue()
            firstIteration = getIteration(self.optimiserFileName)
        else:
            self.insertRelionClassify()
            firstIteration = 1

        lastIteration   = self.NumberOfIterations
        ProtRelionBase.defineSteps2(self, firstIteration
                                        , lastIteration
                                        , self.NumberOfClasses)
    ##########################TEST HERE
    def createFilenameTemplates(self):
        myDict=ProtRelionBase.createFilenameTemplates(self)        
        #myDict['volume']=self.extraIter + "class%(ref3d)03d.spi"
        #myDict['volumeMRC']=self.extraIter + "class%(ref3d)03d.mrc:mrc"
        myDict['volume']=self.extraIter + "class%(ref3d)03d.mrc:mrc"
        
        self.relionFiles += ['model']
        #relionFiles=['data','model','optimiser','sampling']
        for v in self.relionFiles:
            myDict[v+'Re']=self.extraIter + v +'.star'
            myDict[v+'Xm']=self.extraIter + v +'.xmd'
        #myDict['imagesAssignedToClass']='model_classes@'+myDict[v+'Xm']
        return myDict

    def insertRelionClassify(self):
        args = {'--iter': self.NumberOfIterations,
                '--tau2_fudge': self.RegularisationParamT,
                '--flatten_solvent': '',
                #'--zero_mask': '',# this is an option but is almost always true
                '--norm': '',
                '--scale': '',
                '--o': '%s/relion' % self.ExtraDir
                }
        if len(self.ReferenceMask):
            args['--solvent_mask'] = self.ReferenceMask
            
        args.update({'--i': self.ImgStar,
                     '--particle_diameter': self.MaskRadiusA * 2.0,
                     '--angpix': self.SamplingRate,
                     '--ref': self.Ref3D,
                     '--oversampling': '1'
                     })
        
        if not self.IsMapAbsoluteGreyScale:
            args[' --firstiter_cc'] = '' 
            
        if self.InitialLowPassFilterA > 0:
            args['--ini_high'] = self.InitialLowPassFilterA
            
        # CTF stuff
        if self.DoCTFCorrection:
            args['--ctf'] = ''
        
        if self.HasReferenceCTFCorrected:
            args['--ctf_corrected_ref'] = ''
            
        if self.HaveDataPhaseFlipped:
            args['--ctf_phase_flipped'] = ''
            
        if self.IgnoreCTFUntilFirstPeak:
            args['--ctf_intact_first_peak'] = ''
            
        args['--sym'] = self.SymmetryGroup.upper()
        
        args['--K'] = self.NumberOfClasses
            
        # Sampling stuff
        # Find the index(starting at 0) of the selected
        # sampling rate, as used in relion program
        iover = 1 #TODO: check this DROP THIS
        index = ['30','15','7.5','3.7','1.8',
                 '0.9','0.5','0.2','0.1'].index(self.AngularSamplingDeg)
        args['--healpix_order'] = float(index + 1 - iover)
        
        if self.PerformLocalAngularSearch:
            args['--sigma_ang'] = self.LocalAngularSearchRange / 3.
            
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += ' ' + self.AdditionalArguments
        verifyFiles=[]
        #relionFiles=['data','model','optimiser','sampling']
        for v in self.relionFiles:
             verifyFiles += [self.getFilename(v+'Re', iter=self.NumberOfIterations )]
#        f = open('/tmp/myfile','w')
#        for item in verifyFiles:
#            f.write("%s\n" % item)
#        f.close
        self.insertRunJobStep(self.program, params,verifyFiles)

    def insertRelionClassifyContinue(self):
        args = {
                '--o': '%s/relion' % self.ExtraDir,
                '--continue': self.optimiserFileName,
                
                '--iter': self.NumberOfIterations,
                
                '--tau2_fudge': self.RegularisationParamT,# should not be changed 
                '--flatten_solvent': '',# use always
                '--zero_mask': '',# use always. This is confussing since Sjors gui creates the command line
                                  # with this option
                                  # but then the program complains about it. 
                '--oversampling': '1',
                '--norm': '',
                '--scale': '',
                }
        iover = 1 #TODO: check this DROP THIS
        index = ['30','15','7.5','3.7','1.8',
                 '0.9','0.5','0.2','0.1'].index(self.AngularSamplingDeg)
        args['--healpix_order'] = float(index + 1 - iover)
        
        if self.PerformLocalAngularSearch:
            args['--sigma_ang'] = self.LocalAngularSearchRange / 3.
            
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += ' ' + self.AdditionalArguments
        verifyFiles=[]
        #relionFiles=['data','model','optimiser','sampling']
        for v in self.relionFiles:
             verifyFiles += [self.getFilename(v+'Re', iter=self.NumberOfIterations, workingDir=self.WorkingDir )]
        self.insertRunJobStep(self.program, params,verifyFiles)
        #############self.insertRunJobStep('echo shortcut', params,verifyFiles)

