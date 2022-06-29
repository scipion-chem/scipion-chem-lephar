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

import os, glob, shutil

from pwem.convert import AtomicStructHandler
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
import pyworkflow.object as pwobj


from pwchem.utils import runOpenBabel
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from lephar import Plugin as lephar_plugin
from lephar.constants import *

oForm = 'mol2'

class ProtChemFilterLigands(EMProtocol):
    """Perform ligand filtering using QueryDB from LePhar
    http://www.lephar.com/software.htm"""
    _label = 'Ligand filtering'
    _program = "QueryDB/QueryDB"

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ',
                       help="Input set of small molecules for filtering")

        group = form.addGroup('Filtering')
        line = group.addLine('Descriptor',
                             help='Add a descriptor to the filtering. '
                                  'Specify minimum and maximum number of appearances')
        line.addParam('descriptor', params.EnumParam, default=0, label='Descriptor: ',
                      choices=DESC_CHOICES)
        line.addParam('minDesc', params.IntParam, default=0,
                      label='min: ')
        line.addParam('maxDesc', params.IntParam, default=0,
                      label='max: ')

        line = group.addLine('Functional group',
                             help='Add a functional group in SMILES format to the filtering'
                                  'Specify minimum and maximum number of appearances')
        line.addParam('fGroup', params.StringParam, default='', label='Func. group: ')
        line.addParam('minFG', params.IntParam, default=0, label='min: ')
        line.addParam('maxFG', params.IntParam, default=0, label='max: ')

        group.addParam('filterList', params.TextParam, default='', width=70,
                       label='List of filters: ', help='List of filters to perform on the input molecules')


        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('filterStep')
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

    def filterStep(self):
        outFile = 'filtered.mol2'
        args = ' -filter {} {} {} {}'.format(os.path.basename(self.getParametersFile()), oForm,
                                             os.path.basename(self.getLigandsFile(oForm)), outFile)
        lephar_plugin.runLePhar(self, program=self._program, args=args, cwd=self._getPath(),
                                runJob=False, linuxSuf=False)

    def createOutputStep(self):
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        pass

########################### Validation functions #######################

    def _validate(self):
        errors = []
        return errors
      
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

    def parseParamsList(self):
        filterDic = {'descriptors': [], 'fGroups': []}
        for line in self.filterList.get().split('\n'):
            sline = line.split()
            if sline[1] == 'Descriptor':
                filterDic['descriptors'] += [line]
            elif sline[1] == 'Functional':
                filterDic['fGroups'] += [line]
        return filterDic

    def getParametersFile(self):
        paramsFile = self._getPath('custom_lefilter.txt')
        if self.filterList.get().strip() == '':
            shutil.copy(lephar_plugin.getProgramHome(path='bin/QueryDB/le_filter.txt'), paramsFile)
        else:
            filterDic = self.parseParamsList()
            if not os.path.exists(paramsFile):
                with open(paramsFile, 'w') as f:
                    f.write('! Custom le_filter file made with the Scipion-chem framework')
                    for desc in filterDic['descriptors']:
                        descName, descMin, descMax = desc[-3:]
                        f.write('{}\t{}\t{}\n'.format(descName, descMin, descMax))

                    for desc in filterDic['fGroups']:
                        descName, descMin, descMax = desc[-3:]
                        f.write('{}\t{}\t{}\n'.format(descName, descMin, descMax))
        return paramsFile


