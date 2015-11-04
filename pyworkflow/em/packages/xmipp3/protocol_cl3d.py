# **************************************************************************
# *
# * Javier Vargas (jvargas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, isfile
from shutil import copyfile
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,STEPS_PARALLEL,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em import Viewer
import pyworkflow.em.metadata as md
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   getImageLocation,
                                                   readSetOfParticles)
from pyworkflow.em.protocol.protocol_3d import ProtClassify3D


class XmippProtCL3D(ProtClassify3D):
    """    
    Method to obtain a set of 3D classes 
    """
    _label = 'cl3d'
    def __init__(self, *args, **kwargs):
        ProtClassify3D.__init__(self, *args, **kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volumes.')     
                
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles', pointerCondition=['hasAlignment'],
                      help='Select the input projection images.')
            
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParam('angularSampling', FloatParam, default=2,
                      label="Angular Sampling (degrees)",
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('numIntermediateVolumes', FloatParam, default=15,
                      label="Number of intermediate for particle",
                      help='Parameter to define the number of most similar volume \n' 
                      '    projected images for each projection image')
        
                
        form.addParallelSection(threads=1, mpi=1)
                

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):   
        
        convertId = self._insertFunctionStep('convertInputStep', 
                                             self.inputParticles.get().getObjId())
        deps = [] # store volumes steps id to use as dependencies for last step
        commonParams    = self._getCommonParams()
        
        sym = self.symmetryGroup.get()
        initVol= self.inputVolume.get()
        volName = getImageLocation(initVol)            
        volDirGallery  = self._getExtraPath()     
             
        #pmStepId = self._insertFunctionStep('projectionLibraryStep',
        #                                             volName, volDirGallery,
        #                                             prerequisites=[convertId])
        
        for i in range(0, int(self.numIntermediateVolumes.get())):
            
            volName = getImageLocation(initVol)
            volDir = self._getVolDir(i+1)
            
            randomWeightsId = self._insertFunctionStep('reconstructWithRandomWeigths', 
                                                 volName, volDir,sym,                                                  
                                                 prerequisites=[convertId])

            pmStepId = self._insertFunctionStep('projectionLibraryStep', 
                                                 volName, volDir,
                                                 prerequisites=[randomWeightsId])
                
                
            projStepId = self._insertFunctionStep('projectionMatchingStep', 
                                                 volName, volDir,
                                                 commonParams, 
                                                 prerequisites=[pmStepId])
            
            
            prevWeightsId = self._insertFunctionStep('reconstructWithPreviousWeigths', 
                                                 volName, volDir,sym,                                                  
                                                 prerequisites=[projStepId])

                        

            #volStepId = self._insertFunctionStep('angleEvaluationStep', 
            #                                     volName, volDir,
            #                                     sym,
            #                                     prerequisites=[projStepId])
            
            deps.append(prevWeightsId)
            
        self._insertFunctionStep('createOutputStep', 
                                 prerequisites=deps)
        
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        a = self.inputParticles.get()        
        writeSetOfParticles(self.inputParticles.get(),
                            self._getPath('input_particles.xmd'))
                    
    def _getCommonParams(self):                
        params = ' --Ri 0.0'
        params += ' --Ro %0.3f' % ((self.inputParticles.get().getDimensions()[0])/2)
        params += ' --max_shift %0.3f' % ((self.inputParticles.get().getDimensions()[0])/20)
        params += ' --append' 
        params += ' --search5d_shift 0'                            
        return params
                    
    def reconstructWithRandomWeigths(self, volName, volDir, sym):

        #tmpAngleFile = volDir+'/angles_iter000_00.xmd'       
        makePath(volDir)
        #copyfile(self._getPath('input_particles.xmd'),tmpAngleFile)
        params =  ' -i %s' % self._getPath('input_particles.xmd')
        params += ' -o %s' % volDir+'/angles_iter000_00.xmd'
        params += ' --fill weight rand_uniform '
        self.runJob('xmipp_metadata_utilities',
                    params, numberOfMpi=1,numberOfThreads=1)
        
        tmpReconsFile = volDir+'/vol001_00.vol'
        params =  ' -i %s' % volDir+'/angles_iter000_00.xmd'
        params += ' -o %s' % tmpReconsFile 
        params += ' --sym %s' % sym
        params += ' --weight'
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
         
        self.runJob('xmipp_reconstruct_fourier',
                    params, numberOfMpi=nproc,numberOfThreads=nT)

                
    def projectionLibraryStep(self, volName, volDir):
        
        # Generate projections from this reconstruction        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        
        makePath(volDir)
        fnGallery= (volDir+'/gallery.stk')
        params = '-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %f --experimental_images %s --max_tilt_angle 90'\
                    %(volName,fnGallery,self.angularSampling.get(),self.symmetryGroup.get(), -1, volDir+'/angles_iter000_00.xmd')
        
        print params
        self.runJob("xmipp_angular_project_library", params, numberOfMpi=nproc, numberOfThreads=nT)                    
        #self.runJob("ldd","`which xmipp_mpi_angular_project_library`",numberOfMpi=1)                    

    def projectionMatchingStep(self, volName, volDir, params):

        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
        
        params += '  -i %s' % volDir+'/angles_iter000_00.xmd' 
        params += '  -o %s' % (volDir+'/angles_iter001_00.xmd')
        params += ' --ref %s' % (volDir+'/gallery.stk')
        self.runJob('xmipp_angular_projection_matching', 
                    params, numberOfMpi=nproc,numberOfThreads=nT)
        
        
    def reconstructWithPreviousWeigths(self, volName, volDir, sym):

        #tmpAngleFile = volDir+'/angles_iter000_00.xmd'       
        #copyfile(self._getPath('input_particles.xmd'),tmpAngleFile)
        params = ' -i %s' % volDir+'/angles_iter000_00.xmd'
        params +=  ' --set merge %s' % volDir+'/angles_iter001_00.xmd'        
        params += ' -o %s' % volDir+'/angles_iter002_00.xmd'
        self.runJob('xmipp_metadata_utilities',
                    params, numberOfMpi=1,numberOfThreads=1)
        
        tmpReconsFile = volDir+'/vol002_00.vol'
        params =  ' -i %s' % volDir+'/angles_iter002_00.xmd'
        params += ' -o %s' % tmpReconsFile 
        params += ' --sym %s' % sym
        params += ' --weight'
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
         
        self.runJob('xmipp_reconstruct_fourier',
                    params, numberOfMpi=nproc,numberOfThreads=nT)
        

    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        imgSet = self.inputParticles.get()
        for i, vol in enumerate(self._iterInputVols()):
            volume = vol.clone()               
            volDir = self._getVolDir(i+1)
            volPrefix = 'vol%03d_' % (i+1)
            validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
            moveFile(join(volDir, 'validation.xmd'), 
                     validationMd)
            clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
            moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
            
            outImgSet = self._createSetOfParticles(volPrefix)
            
            outImgSet.copyInfo(imgSet)
            #readSetOfParticles(String(clusterMd), outImgSet)

            outImgSet.copyItems(imgSet,
                                updateItemCallback=self._setWeight,
                                itemDataIterator=md.iterRows(clusterMd))
                        
            mdValidatoin = md.MetaData(validationMd)
            weight = mdValidatoin.getValue(md.MDL_WEIGHT, mdValidatoin.firstObject())
            volume.weight = Float(weight)
            volume.clusterMd = String(clusterMd)
            volume.cleanObjId() # clean objects id to assign new ones inside the set            
            outputVols.append(volume)
            self._defineOutputs(outputParticles=outImgSet)
        
        outputVols.setSamplingRate(volume.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        #self._defineTransformRelation(self.inputVolumes.get(), volume)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolume.get() and not self.inputVolume.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
        
    
    def _summary(self):
        summary = []

        summary.append("Input particles:  %s" % self.inputParticles.get().getNameId())

        summary.append("-----------------")
        if self.inputVolumes.get():
            for i, vol in enumerate(self._iterInputVols()):
                summary.append("Input volume(s)_%d: [%s]" % (i+1,vol))

        summary.append("-----------------")
        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            for i, vol in enumerate(self._iterInputVols()):
                            
                VolPrefix = 'vol%03d_' % (i+1)
                md = xmipp.MetaData(self._getExtraPath(VolPrefix+'validation.xmd'))                
                weight = md.getValue(xmipp.MDL_WEIGHT, md.firstObject())
                summary.append("Output volume(s)_%d : %s" % (i+1,self.outputVolumes.getNameId()))
                summary.append("Quality parameter_%d : %f" % (i+1,weight))
                summary.append("-----------------")        
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and number of orientations of %f' % (self.angularSampling.get(), self.numOrientations.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getExtraPath('vol%03d' % volIndex)
    
        
    def _defineMetadataRootName(self, mdrootname,volId):
        
        if mdrootname=='P':
            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'clusteringTendency.xmd')
        if mdrootname=='Volume':

            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'validation.xmd')
            
    def _definePName(self):
        fscFn = self._defineMetadataRootName('P')
        return fscFn
    
    def _defineVolumeName(self,volId):
        fscFn = self._defineMetadataRootName('Volume',volId)
        return fscFn
    
    def _setWeight(self, item, row):  
        item._xmipp_weightClusterability = Float(row.getValue(md.MDL_VOLUME_SCORE1))
