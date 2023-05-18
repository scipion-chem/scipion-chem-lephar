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

import os

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pyworkflow.protocol.params import PointerParam, BooleanParam, StringParam

from pwchem.utils import clean_PDB, removeNumberFromStr, getChainIds
from pwchem.protocols import ProtChemPrepareReceptor

from lephar import Plugin as lephar_plugin


class ProtChemLePro(ProtChemPrepareReceptor):
    """Perform a target preparation using the LePro binary from LePhar:
    http://www.lephar.com/software.htm
    """
    _label = 'LePro target preparation'
    _program = "lepro"

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                       label='Input atomic structure:',
                       help="The atom structure to be prepared")

        clean = self.defineCleanParams(form, w=False, hk=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutputStep')

    def preparationStep(self):
        # Clean PDB
        pdb_ini = self.inputAtomStruct.get().getFileName()
        filename = os.path.splitext(os.path.basename(pdb_ini))[0]
        fnPdb = self._getExtraPath('%s_clean.pdb' % filename)

        chain_ids = None
        if self.rchains.get():
            chain_ids = getChainIds(self.chain_name.get())

        cleanedPDB = clean_PDB(self.inputAtomStruct.get().getFileName(), fnPdb,
                               False, self.HETATM.get(), chain_ids)

        args = os.path.abspath(cleanedPDB)
        lephar_plugin.runLePhar(self, program=self._program, args=args, cwd=self._getExtraPath())

    def createOutputStep(self):
        outFileName = self._getPath(self._getInputName() + '_prep.pdb')
        os.rename(self._getExtraPath('pro.pdb'), outFileName)
        self.addPDBColumns(outFileName)
        outAS = AtomStruct(outFileName)
        self._defineOutputs(outputStructure=outAS)

    ########################### Utils functions ############################

    def _getInputName(self):
        return os.path.splitext(os.path.basename(self.inputAtomStruct.get().getFileName()))[0]

    def addPDBColumns(self, pdbFile):
        auxFile = self._getTmpPath(os.path.basename(pdbFile))
        with open(pdbFile) as fIn:
            with open(auxFile, 'w') as f:
                for line in fIn:
                    if line.startswith('ATOM'):
                        atomSym = removeNumberFromStr(line.split()[2])
                        line = line.strip() + '  1.00  0.00{}{}\n'.format(' '*11, atomSym)
                    f.write(line)
        os.rename(auxFile, pdbFile)
        return pdbFile