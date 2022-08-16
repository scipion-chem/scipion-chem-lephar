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

import os, shutil, json

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pyworkflow.protocol.params import PointerParam, BooleanParam, StringParam

from pwchem.utils import clean_PDB

from lephar import Plugin as lephar_plugin


class ProtChemLePro(EMProtocol):
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

        clean = form.addGroup('Clean Structure File')
        clean.addParam("HETATM", BooleanParam,
                       label='Remove ligands HETATM',
                       default=True, important=True,
                       help='Remove all ligands and HETATM contained in the protein')

        clean.addParam("rchains", BooleanParam,
                       label='Select specific chains: ',
                       default=False, important=True,
                       help='Keep only the chains selected')

        clean.addParam("chain_name", StringParam,
                       label="Keep chains: ", important=True,
                       condition="rchains==True",
                       help="Select the chain(s) you want to keep in the structure")

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
            chainJson = json.loads(self.chain_name.get())  # From wizard dictionary
            if 'chain' in chainJson:
                chain_ids = [chainJson["chain"].upper().strip()]
            elif 'model-chain' in chainJson:
                modelChains = chainJson["model-chain"].upper().strip()
                chain_ids = [x.split('-')[1] for x in modelChains.split(',')]

        cleanedPDB = clean_PDB(self.inputAtomStruct.get().getFileName(), fnPdb,
                               False, self.HETATM.get(), chain_ids)

        args = os.path.abspath(cleanedPDB)
        lephar_plugin.runLePhar(self, program=self._program, args=args, cwd=self._getExtraPath())

    def createOutputStep(self):
        outFileName = self._getPath(self._getInputName() + '_prep.pdb')
        shutil.copy(self._getExtraPath('pro.pdb'), outFileName)
        outAS = AtomStruct(outFileName)
        self._defineOutputs(outputStructure=outAS)

    ########################### Utils functions ############################

    def _getInputName(self):
        return os.path.splitext(os.path.basename(self.inputAtomStruct.get().getFileName()))[0]
