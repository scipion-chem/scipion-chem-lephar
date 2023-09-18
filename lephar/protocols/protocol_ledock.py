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

import os, time

from pwem.convert import AtomicStructHandler
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, BooleanParam, LEVEL_ADVANCED
import pyworkflow.object as pwobj


from pwchem.utils import removeNumberFromStr, performBatchThreading, runInParallel, obabelMolConversion, makeSubsets
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from lephar import Plugin as lephar_plugin
from lephar.constants import *


class ProtChemLeDock(EMProtocol):
    """Perform a docking experiment with LeDock, from LePhar software
    http://www.lephar.com/software.htm"""
    _label = 'LeDock docking'
    _program = "ledock"

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('wholeProt', BooleanParam, label='Dock on whole protein: ', default=False,
                      help='Whether to dock on a whole protein surface or on specific regions')

        #Docking on whole protein
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure:', condition='wholeProt', allowsNull=True,
                      help="The atom structure to use as receptor in the docking")
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='wholeProt', allowsNull=False,
                       help='Radius of the Autodock grid for the whole protein')

        #Docking on pockets
        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                      label='Input pockets:', condition='not wholeProt', allowsNull=True,
                      help="The protein structural ROIs to dock in")
        group.addParam('pocketRadiusN', FloatParam, label='Grid radius vs ROI radius: ',
                       condition='not wholeProt', default=1.1, allowsNull=False,
                       help='The radius * n of each structural ROI will be used as grid radius')

        group = form.addGroup('Docking')
        group.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Ligand library: ', allowsNull=False,
                       help="Input set of small molecules to dock with LeDock")
        group.addParam('rmsTol', FloatParam, label='Cluster RMSD (A): ', default=1.0, expertLevel=LEVEL_ADVANCED,
                       help='RMSD threshold for discarding too similar docking poses')
        group.addParam('nRuns', IntParam, label='Number of positions per ligand: ', default=10,
                       help='Maximum number of poses to output per ligand per StructROI')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        nThreads = self.getnThreads()
        cId = self._insertFunctionStep('convertStep', prerequisites=[])

        dockSteps, splitSteps = [], []
        if self.wholeProt:
            pocketDir = self.getOutputPocketDir()
            os.mkdir(pocketDir)
            for i in range(nThreads):
                dockId = self._insertFunctionStep('dockStep', None, i, poprerequisites=[cId])
            dockSteps.append(dockId)

            splitId = self._insertFunctionStep('splitStep', None, nThreads, prerequisites=dockSteps)
            splitSteps.append(splitId)
            
        else:
            for pocket in self.inputStructROIs.get():
                pocketDir = self.getOutputPocketDir(pocket)
                os.mkdir(pocketDir)
                for i in range(nThreads):
                    dockId = self._insertFunctionStep('dockStep', pocket.clone(), i, prerequisites=[cId])
                    dockSteps.append(dockId)
                    
            for pocket in self.inputStructROIs.get():
                splitId = self._insertFunctionStep('splitStep', pocket.clone(), nThreads, prerequisites=dockSteps)
                splitSteps.append(splitId)
                
        self._insertFunctionStep('createOutputStep', prerequisites=splitSteps)

    def convertStep(self):
        outDir = self._getExtraPath()
        #Receptor as prepared by lepro
        lephar_plugin.runLePhar(self, 'lepro', args=os.path.abspath(self.getOriginalReceptorFile()), cwd=outDir)

        # Ligands in mol2 format
        self.convertAndWriteMolSet(self.inputSmallMolecules.get(), outDir, self.numberOfThreads.get())

    def dockStep(self, pocket=None, idx=None):
        oDir = self.getOutputPocketDir(pocket)
        dockParamFile, _ = self.writeDockInFile(pocket, idx=idx)
        lephar_plugin.runLePhar(self, program=self._program, args=dockParamFile, cwd=oDir)

    def splitStep(self, pocket=None, nThreads=None):
        oDir = self.getOutputPocketDir(pocket)
        dockFiles = self.getDockFiles(oDir)
        performBatchThreading(self.performSplit, dockFiles, nThreads, cloneItem=False, outDir=oDir)

    def createOutputStep(self):
        allFiles = []
        for pocketDir in self.getPocketDirs():
            for file in os.listdir(pocketDir):
                outDir = os.path.join(pocketDir, file)
                if os.path.isdir(outDir):
                    for outFile in os.listdir(outDir):
                        allFiles.append(os.path.join(outDir, outFile))
        performBatchThreading(self.correctMolFile, allFiles, self.numberOfThreads.get(), cloneItem=False)

        inputMolDic = self.getInputMolsDic()
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for pocketDir in self.getPocketDirs():
            gridId = self.getGridId(pocketDir)
            for file in os.listdir(pocketDir):
                outDir = os.path.join(pocketDir, file)
                if os.path.isdir(outDir):
                    inMol = inputMolDic[file.split('.')[0]]
                    for outFile in os.listdir(outDir):
                        newSmallMol = SmallMolecule()
                        newSmallMol.copy(inMol, copyId=False)

                        molFile = os.path.join(outDir, outFile)
                        newSmallMol._energy = pwobj.Float(self.parseEnergy(molFile))
                        if os.path.getsize(molFile) > 0:
                            newSmallMol.poseFile.set(molFile)
                            newSmallMol.setPoseId(molFile.split('_')[-1].split('.')[0])
                            newSmallMol.gridId.set(gridId)
                            newSmallMol.setMolClass('LeDock')
                            newSmallMol.setDockId(self.getObjId())
                            outputSet.append(newSmallMol)

        outputSet.proteinFile.set(self.getOriginalReceptorFile())
        outputSet.setDocked(True)
        self._defineOutputs(outputSmallMolecules=outputSet)

########################### Validation functions #######################

    def _validate(self):
        errors = []
        if self.wholeProt:
            if not self.inputAtomStruct.get():
                errors.append('You need to specify an input atom structure')
            elif not self.radius.get():
                errors.append('You need to specify a radius. You may use the wizard to do so')
        else:
            if not self.inputStructROIs.get():
                errors.append('You need to specify an input set of StructROIs')
            elif not self.pocketRadiusN.get():
                errors.append('You need to specify a radius coefficient to adjust the StructROIs radius.')
        return errors

    def _citations(self):
        return ['C6CP01555G']
      
########################### Utils functions ############################

    def convertAndWriteMolSet(self, molSet, outDir, nJobs):
        convMolFiles = runInParallel(obabelMolConversion, '.mol2', outDir, paramList=[item.clone() for item in molSet],
                                     jobs=nJobs)
        molFileSubsets = makeSubsets(convMolFiles, self.getnThreads(), cloneItem=False)

        with open(self.getLigandListFile(), 'w') as fLig:
            for iSet, molFSet in enumerate(molFileSubsets):
                with open(self.getLigandListFile(base=True, idx=iSet), 'w') as fLigBase:
                    for molFile in molFSet:
                        fLig.write(molFile + '\n')
                        fLigBase.write(os.path.basename(molFile) + '\n')

    def performSplit(self, dockFiles, molLists, it, outDir):
        for dockFile in dockFiles:
            dockBase = os.path.basename(dockFile)
            dockRoot = dockBase.split('.')[0]
            dockDir = os.path.join(outDir, dockRoot)
            os.mkdir(dockDir)

            newDockFile = os.path.join(dockDir, dockBase)
            os.rename(dockFile, newDockFile)
            args = ' -spli ' + os.path.abspath(newDockFile)
            lephar_plugin.runLePhar(self, program=self._program, args=args, cwd=outDir, runJob=False)
            os.remove(newDockFile)

    def correctMolFile(self, molFiles, molsLists, it):
        for molFile in molFiles:
            newFile = self.renameDockFile(molFile)
            with open(molFile) as fIn:
                with open(newFile, 'w') as f:
                    for line in fIn:
                        if line.startswith('ATOM'):
                            atomSym = removeNumberFromStr(line.split()[2])
                            line = line.strip() + '  1.00  0.00{}{}\n'.format(' ' * (12 - len(atomSym)), atomSym)
                        f.write(line)
            os.remove(molFile)

    def getDockFiles(self, pocketDir):
        dockFiles = []
        for file in os.listdir(pocketDir):
            if file.endswith('.dok'):
                dockFiles.append(os.path.abspath(os.path.join(pocketDir, file)))
        return dockFiles

    def getnThreads(self):
        '''Get the number of threads available for each pocket execution'''
        nThreads = self.numberOfThreads.get()
        if not self.wholeProt:
            nPockets = len(self.inputStructROIs.get())
            nThreads = nThreads // nPockets
        nThreads = 1 if nThreads == 0 else nThreads
        return nThreads

    def parseEnergy(self, molFile):
        with open(molFile) as fMol:
            fMol.readline()
            line = fMol.readline()
        return line.split()[-2]


    def getInputMolsDic(self):
        dic = {}
        for mol in self.inputSmallMolecules.get():
            dic[mol.clone().getUniqueName(False, True, False, False)] = mol.clone()
        return dic

    def renameDockFile(self, outFile):
        outBase = os.path.basename(outFile)
        newBase = '{}{}.pdb'.format(outBase.split('dock')[0],
                                    int(outBase.split('dock')[-1].split('.')[0]))
        newFile = os.path.join(os.path.dirname(outFile), newBase)
        return newFile

    def getGridId(self, outDir):
        return outDir.split('_')[-1]

    def getPocketDirs(self):
        dirs = []
        for file in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(file)
            if os.path.isdir(d) and 'pocket' in file:
                dirs.append(d)
        dirs.sort()
        return dirs

    def getOutputPocketDir(self, pocket=None):
        if pocket==None:
            outDir = self._getExtraPath('pocket_1')
        else:
            outDir = self._getExtraPath('pocket_{}'.format(pocket.getObjId()))
        return outDir

    def getOriginalReceptorFile(self):
        if self.wholeProt:
            return self.inputAtomStruct.get().getFileName()
        else:
            return self.inputStructROIs.get().getProteinFile()

    def getPreparedReceptorFile(self):
        return self._getExtraPath('pro.pdb')

    def getLigandListFile(self, base=False, idx=None):
        if not base:
            return os.path.abspath(self._getPath('ligands.list'))
        else:
            return os.path.abspath(self._getPath('ligandsBase_{}.list'.format(idx)))
            

    def writeDockInFile(self, pocket, idx):
        pDir = self.getOutputPocketDir(pocket)
        if not self.wholeProt:
            x_center, y_center, z_center = pocket.calculateMassCenter()
            r = pocket.getDiameter() / 2
        else:
            ASH = AtomicStructHandler(self.getPreparedReceptorFile())
            x_center, y_center, z_center = ASH.centerOfMass()
            r = self.radius.get()

        xmin, xmax = x_center - r, x_center + r
        ymin, ymax = y_center - r, y_center + r
        zmin, zmax = z_center - r, z_center + r

        localReceptor = self.linkLocal(os.path.abspath(self.getPreparedReceptorFile()), pDir)
        localLigList = self.doLocalLig(pDir, idx)

        strIn = DOCK_IN.format(localReceptor, self.rmsTol.get(), xmin, xmax, ymin, ymax, zmin, zmax,
                               self.nRuns.get(), localLigList)

        dockFile = os.path.abspath(os.path.join(pDir, 'dock_{}.in'.format(idx)))
        with open(dockFile, 'w') as fIn:
            fIn.write(strIn)

        return dockFile, localLigList

    def linkLocal(self, sourcePath, outDir):
        outFile = os.path.join(outDir, os.path.basename(sourcePath))
        if not os.path.exists(outFile):
            os.symlink(sourcePath, outFile)
        return os.path.basename(sourcePath)

    def doLocalLig(self, outDir, idx=0):
        if idx == 0:
            with open(self.getLigandListFile()) as fIn:
                for molFile in fIn:
                    self.linkLocal(molFile.strip(), outDir)

        return self.linkLocal(self.getLigandListFile(base=True, idx=idx), outDir)


