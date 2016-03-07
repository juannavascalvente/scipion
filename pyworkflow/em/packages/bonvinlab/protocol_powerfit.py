# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

import pyworkflow.protocol.params as params
import pyworkflow.em as em


class ProtPowerFit(em.EMProtocol):
    """
    PowerFit is a Python package and simple command-line program to automatically
    fit high-resolution atomic structures in cryo-EM densities. To this end it
    performs a full-exhaustive 6-dimensional cross-correlation search between
    the atomic structure and the density. It takes as input an atomic structure
    in PDB-format and a cryo-EM density with its resolution; and outputs positions
    and rotations of the atomic structure corresponding to high correlation values.
    PowerFit uses the Local Cross-Correlation functions as its base score.
    The score can optionally be enhanced by a Laplace pre-filter and/or a
    core-weighted version to minimize overlapping densities from neighboring
    subunits. It can further be hardware-accelerated by leveraging multi-core
    CPU machines out of the box or by GPU via the OpenCL framework.
    """
    _label = 'powerfit'
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputPDB', params.PointerParam,
                      pointerClass='PdbFile',
                      label='Input PDB', important=True,
                      help="Atomic model to be fitted in the density. "
                           "Format should either be PDB or mmCIF")

        form.addParam('inputMap', params.PointerParam,
                      pointerClass='Volume',
                      label='Input map', important=True,
                      help="Target density map to fit the model in.")

        form.addParam('resolution', params.FloatParam,
                      label='Resolution of the map (A)',
                      help='Resolution of the map in Angstroms')

        form.addParam('angle', params.IntParam, default=10,
                      label='Rotational sampling (deg)',
                      help="Rotational sampling density in degree. "
                           "Increasing this number by a factor of 2 results in "
                           "approximately 8 times more rotations sampled.")

        form.addParam('useLaplace', params.BooleanParam, default=False,
                      label='Use the Lapacle pre-filter?',
                      help="Use the Laplace pre-filter density data. Can be "
                           "combined with the core-weighted local "
                           "cross-correlation.")

        form.addParam('coreWeighted', params.BooleanParam, default=False,
                      label='Use core-weighted local CC?',
                      help="Use core-weighted local cross-correlation score. "
                           "Can be combined with the Laplace pre-filter.")

        form.addSection('Extra')
        form.addParam('resampleMap', params.BooleanParam, default=True,
                      level=params.LEVEL_ADVANCED,
                      label='Resample the density map?')

        form.addParam('trimMap', params.BooleanParam, default=True,
                      level=params.LEVEL_ADVANCED,
                      label='Trim the density map?')

        form.addParam('trimCutoff', params.FloatParam, default=10,
                      label='Trimming cut-off', condition='trimMap',
                      help="Intensity cutoff to which the map will be trimmed."
                           "Default is 10 percent of maximum intensity.")

        form.addParam('useGPU', params.BooleanParam, default=True,
                      label='Use GPU',
                      help='Set to Yes to use GPU as well as CPU')

        form.addParallelSection(threads=4, mpi=0)
        
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', 
                                 self.inputPDB.get().getObjId(),
                                 self.inputMap.get().getObjId())
        self._insertPowerFitStep()
        self._insertFunctionStep('createOutputStep')
        
    def _insertPowerFitStep(self):
        """ Insert the steps to launch the pikcer with different modes. 
        Prepare the command arguments that will be passed. 
        """
        args =  self.inputPDB.get().getFileName()
        args += ' %s %f' % (self._getMapMrc(), self.resolution)
        args += ' --angle %d' % self.angle
        args += ' --directory %s' % self._getExtraPath()
        
        if self.useLaplace:
            args += ' --laplace'

        if self.coreWeighted:
            args += ' --core-weighted'

        if not self.resampleMap:
            args += ' --no-resample'

        if self.trimMap:
            args += ' --trimming-cutoff %f' % self.trimCutoff
        else:
            args += ' --no-trimming'

        self._insertFunctionStep('runPowerFitStep', args)

    #--------------------------- STEPS functions ---------------------------------------------------
    def _getMapMrc(self):
        return self._getExtraPath('input_volume.mrc')

    def convertInputStep(self, micsId, refsId):
        """ This step will take of the convertions from the inputs.
        Convert the input map to .mrc file.
        """
        ih = em.ImageHandler()
        ih.convert(self.inputMap.get(), self._getMapMrc())

    def runPowerFitStep(self, args):
        args += ' --nproc %d' % self.numberOfThreads

        if self.useGPU:
            args += ' --gpu'

        self.runJob('powerfit', args)

    def createOutputStep(self):
        pass

    def _validate(self):
        errors = []
        return errors
        
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methodsMsgs = []
        return methodsMsgs
    
    def _citations(self):
        return ['Bonvin2015']
    
    #--------------------------- UTILS functions --------------------------------------------------
