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

    AI Generated:

      ProtChemLePro - User Manual

      Overview
      --------
      The ProtChemLePro protocol prepares receptor structures for molecular
      docking using the LePro tool from the LePhar suite. It performs cleaning,
      formatting, and structural standardization of protein targets to ensure
      compatibility with downstream docking workflows.

      This protocol is a key preprocessing step that ensures the receptor is
      chemically valid, properly formatted, and ready for accurate docking
      simulations.

      Input Requirements
      ------------------
      1. **Receptor Structure**:
         - AtomStruct object containing the protein structure.
         - Input must be in PDB or compatible format.

      2. **Structure Quality**:
         - Should include relevant residues for docking.
         - May contain heteroatoms, ligands, or waters (which can be removed).

      Workflow
      --------
      1. **Structure Cleaning**:
         - Removes unwanted atoms such as:
           - Water molecules
           - Crystallographic ligands (optional)
           - Irrelevant heteroatoms
         - Optionally filters specific chains.

      2. **Structure Preparation**:
         - Adds hydrogens where required.
         - Adjusts protonation states if necessary.
         - Ensures chemical correctness of the receptor.

      3. **Binding Site Definition**:
         - Prepares receptor for docking using defined parameters.
         - Retains structural information needed for binding site characterization.

      4. **LePro Execution**:
         - Runs the LePro binary from LePhar.
         - Processes the cleaned receptor structure.

      5. **Output Formatting**:
         - Converts the processed structure into docking-ready format.
         - Adds required columns and annotations for compatibility.

      Outputs
      -------
      - **Prepared Receptor Structure**:
        - Output AtomStruct containing the cleaned and formatted receptor.
        - Stored in PDB format with updated atom annotations.

      - **Formatted Structure Files**:
        - Intermediate cleaned PDB file.
        - Final processed file compatible with docking protocols.

      Advanced Options
      ----------------
      - Selection of specific chains for preparation.
      - Retention or removal of heteroatoms.
      - Integration with broader docking workflows.
      - Automatic cleaning and formatting utilities.

      Validation & Warnings
      ---------------------
      - Input structure must be valid and correctly formatted.
      - Missing residues or structural inconsistencies may affect results.
      - Incorrect chain selection may lead to incomplete preparation.
      - Ensure input is suitable for docking before running the protocol.

      Practical Recommendations
      -------------------------
      - Inspect input structures before preparation.
      - Remove unnecessary ligands unless they are part of the binding site.
      - Use consistent chain selection across protocols.
      - Validate output before proceeding to docking.

      Final Perspective
      -----------------
      ProtChemLePro provides a robust and automated method for preparing receptor
      structures for docking using LePro. It ensures that protein targets are
      chemically consistent, structurally clean, and properly formatted, forming
      an essential step in high-quality structure-based drug discovery workflows.

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