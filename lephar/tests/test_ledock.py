# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from ..protocols import ProtChemLePro, ProtChemLeDock
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands, ProtDefineStructROIs


class TestLePro(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='5ni1')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPrepareReceptorLePro(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemLePro,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            HETATM=True, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 141}',
            repair=3)

        cls.launchProtocol(cls.protPrepareReceptor)

    def test(self):
        self._runPrepareReceptorLePro()

        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
        self.assertIsNotNone(getattr(self.protPrepareReceptor, 'outputStructure', None))


class TestLeDock(TestLePro):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportSmallMols()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

        cls._runPrepareLigandsOBabel()
        cls._runPrepareReceptorLePro()
        cls._waitOutput(cls.protOBabel, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protPrepareReceptor, 'outputStructure', sleepTime=5)

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'))
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
        cls.protOBabel = cls.newProtocol(
            ProtChemOBabelPrepareLigands,
            inputType=0, method_charges=0,
            inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)

        cls.proj.launchProtocol(cls.protOBabel, wait=False)

    @classmethod
    def _runPocketsSearch(cls):
        cls.protDefPockets = cls.newProtocol(
            ProtDefineStructROIs,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure,
            inResidues='{"model": 0, "chain": "C", "index": "58-58", "residues": "H"}\n'
                       '{"model": 0, "chain": "C", "index": "101-101", "residues": "L"}\n')

        cls.proj.launchProtocol(cls.protDefPockets, wait=False)
        return cls.protDefPockets

    def _runLeDock(self, pocketsProt=None):
        if pocketsProt == None:
            protAutoDock = self.newProtocol(
                ProtChemLeDock,
                wholeProt=True,
                inputAtomStruct=self.protPrepareReceptor.outputStructure,
                inputSmallMolecules=self.protOBabel.outputSmallMolecules,
                radius=24, gaRun=2,
                numberOfThreads=8)
            self.proj.launchProtocol(protAutoDock, wait=False)

        else:
            protAutoDock = self.newProtocol(
                ProtChemLeDock,
                wholeProt=False,
                inputStructROIs=pocketsProt.outputStructROIs,
                inputSmallMolecules=self.protOBabel.outputSmallMolecules,
                pocketRadiusN=5, gaRun=2,
                numberOfThreads=4)
            self.proj.launchProtocol(protAutoDock, wait=False)

        return protAutoDock

    def test(self):
        print('Docking with LeDock in the whole protein')
        protLeDock1 = self._runLeDock()

        print('Docking with LeDock in predicted pockets')
        protStructROIs = self._runPocketsSearch()
        self._waitOutput(protStructROIs, 'outputStructROIs', sleepTime=5)

        protLeDock2 = self._runLeDock(self.protDefPockets)

        self._waitOutput(protLeDock1, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protLeDock1, 'outputSmallMolecules', None))
        self._waitOutput(protLeDock2, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protLeDock2, 'outputSmallMolecules', None))


