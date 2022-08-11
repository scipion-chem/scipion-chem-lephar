# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os, shutil

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, BooleanParam, LEVEL_ADVANCED
import pyworkflow.object as pwobj


from pwchem.utils import runOpenBabel
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from lephar import Plugin as lephar_plugin

oForm = 'mol2'

class ProtChemClusterMCS(EMProtocol):
    """Perform a molecule structure clustering using CLusterByMCSbinary from LePhar
    http://www.lephar.com/software.htm"""
    _label = 'LePhar molecule clustering'
    _program = "ClusterByMCS"

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ',
                       help="Input set of small molecules to dock with LeDock")
        group.addParam('clustCut', FloatParam, default=0.618,
                       label='Clustering cutoff: ',
                       help='The cutoff value ranges from 0 to 1. The default is 0.618. '
                            'However, 0.8 might be generally more suitable')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('clusterStep')
        self._insertFunctionStep('createOutputStep')

    def convertStep(self):
        # Combine ligands into single mol2 file (clean @<TRIPOS>UNITY_ATOM_ATTR section)
        with open(self.getLigandsFile(oForm), 'w') as fLig:
            for mol in self.inputSmallMolecules.get():
                molFile = mol.getFileName()
                if not molFile.endswith('.{}'.format(oForm)):
                    inName, inExt = os.path.splitext(os.path.basename(molFile))
                    oFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.{}'.format(oForm))))

                    args = ' -i{} {} -o{} -O {}'.format(inExt[1:], os.path.abspath(molFile), oForm, oFile)
                    runOpenBabel(protocol=self, args=args, cwd=self._getExtraPath())
                    molFile = oFile
                with open(molFile) as fIn:
                    fLig.write(fIn.read() + '\n')

        self.cleanAttrSection(self.getLigandsFile(oForm))

    def clusterStep(self):
        args = ' {} {} {}'.format(oForm, self.getLigandsFile(oForm), self.clustCut.get())
        lephar_plugin.runRDKit2Script(self, scriptName=self._program, args=args, cwd=self._getPath())

    def createOutputStep(self):
        clustersDic = self.parseClusters()
        for clusterId in clustersDic:
            outputSet = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix=clusterId)
            for mol in clustersDic[clusterId]:
                outputSet.append(mol)
            self._defineOutputs(**{'outputSmallMolecules_{}'.format(clusterId): outputSet})

########################### Validation functions #######################

    def _validate(self):
        errors = []
        return errors

    def _warnings(self):
        warns = []
        if hasattr(self.inputSmallMolecules.get().getFirstItem(), '_ConformersFile'):
            warns.append('Molecules where conformers have been generated may produce some errors in the parsing')
        return warns

    def _citations(self):
        return ['C6CP01555G']
      
########################### Utils functions ############################

    def getLigandsFile(self, oFormat):
        return os.path.abspath(self._getPath('ligands.{}'.format(oFormat)))

    def cleanAttrSection(self, mol2File, outFile=None):
        '''Remove the @<TRIPOS>UNITY_ATOM_ATTR section from a mol2 file which causes problems with rdkit parsing
        If not outFile provided, mol2File will be replaced
        '''
        if not outFile:
            auxFile = mol2File.replace('.mol2', '_aux.mol2')
        else:
            auxFile = outFile

        inSection = False
        with open(auxFile, 'w') as f:
            with open(mol2File) as fIn:
                for line in fIn:
                    if line.startswith('@<TRIPOS>UNITY_ATOM_ATTR'):
                        inSection = True
                    elif line.startswith('@<'):
                        inSection = False

                    if not inSection:
                        f.write(line)

        if not outFile:
            outFile = mol2File
            shutil.copy(auxFile, outFile)
            os.remove(auxFile)

        return outFile

    def parseClusters(self):
        clusters = {}
        with open(self._getPath('clusters.smi')) as f:
            for line in f:
                smi, name, _, clusterId = line.split()
                if clusterId in clusters:
                    num += 1
                    molFile = self.writeMol2File(smi, '{}_{}'.format(name, num))
                    clusters[clusterId] += [SmallMolecule(smallMolFilename=molFile, molName='guess')]
                else:
                    num = 1
                    molFile = self.writeMol2File(smi, '{}_{}'.format(name, num))
                    clusters[clusterId] = [SmallMolecule(smallMolFilename=molFile, molName='guess')]

        return clusters

    def writeMol2File(self, smi, name):
        oFile = self._getExtraPath(name + '.' + oForm)

        args = ' -:"{}" -o{} -O {} '.format(smi, oForm, os.path.abspath(oFile))
        runOpenBabel(protocol=self, args=args, cwd=self._getExtraPath())
        return os.path.abspath(oFile)
