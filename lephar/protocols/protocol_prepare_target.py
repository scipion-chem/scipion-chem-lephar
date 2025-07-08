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

from pwchem.utils import cleanPDB, removeNumberFromStr, getChainIds
from pwchem.protocols import ProtChemPrepareReceptor

from lephar import Plugin as lephar_plugin


class ProtChemLePro(ProtChemPrepareReceptor):
    """Perform a target preparation using the LePro binary from LePhar:
    http://www.lephar.com/software.htm

User IA Manual: PrepareTarget Protocol

The PrepareTarget protocol is used to process and format a receptor structure
so that it is ready for molecular docking within the LePhar workflow. It handles
key steps such as cleaning the structure, defining the binding site, and
converting the file into the format required by the docking engine.

The user begins by providing a receptor in PDB format. This structure should
include all atoms necessary for docking but may also contain crystallographic
waters, ligands, or other heteroatoms that must be removed prior to docking.
The protocol parses the input and filters out components that are not relevant
for ligand binding, while retaining the atoms necessary to define the pocket
environment.

To define the region where docking will take place, the user must set the center
of the binding site using three-dimensional coordinates. These values can be
derived from a known ligand, a predicted pocket, or by visual inspection of the
receptor. The dimensions of the docking box are also specified at this stage,
ensuring that the docking engine will search the correct volume of space.

The protocol can automatically add missing hydrogen atoms to the receptor,
adjust protonation states, and ensure that the final structure is chemically
valid. This step is essential for maintaining the physical integrity of the
receptor and for producing reliable docking results. The user may also choose to
retain or discard cofactors or metal ions, depending on whether they are relevant
to ligand binding.

Once the structure has been cleaned and the binding site defined, the receptor
is converted into the appropriate docking format. The final output includes a
MOL2 file containing the processed receptor and a configuration file that
records the box definition and other docking parameters. These outputs are used
directly by the docking protocol and ensure consistency across the virtual
screening workflow.

In summary, the PrepareTarget protocol standardizes and formats receptor
structures for docking, offering a reproducible and automated method to define
binding sites and prepare the physical model of the target. It serves as a
critical step before molecular docking and ensures that all required inputs are
correctly structured and ready for high-throughput screening.

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

        cleanedPDB = cleanPDB(self.inputAtomStruct.get().getFileName(), fnPdb,
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