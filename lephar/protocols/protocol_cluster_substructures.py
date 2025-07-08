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
    http://www.lephar.com/software.htm
    
    User IA Manual: ClusterSubstructures Protocol

The ClusterSubstructures protocol is designed to analyze a set of chemical
structures and group them based on shared substructural features. It enables
the identification of common scaffolds or motifs across a molecular dataset,
facilitating structure?activity relationship analysis, chemical diversity
assessment, or hit expansion strategies within virtual screening workflows.

To use the protocol, the user must provide a collection of ligands or compounds
in a format that includes molecular connectivity. These can originate from prior
docking, enumeration, or library preparation steps. Each compound is examined
to identify relevant substructures, and molecules are compared to determine
their level of shared chemical features.

The user can select how substructures are extracted and compared, typically
based on molecular fingerprints, scaffold definitions, or graph-based similarity.
The level of clustering sensitivity can be adjusted, allowing either fine-grained
separation based on small differences or broader grouping around central cores.
Thresholds can be set to control how similar two molecules must be to be placed
in the same cluster.

In addition to the similarity metric, the protocol allows configuration of the
minimum cluster size to retain, which helps eliminate noise or outlier compounds.
Clustering methods may be hierarchical or fingerprint-based, depending on the
chosen algorithm. The resulting clusters reflect substructure-based
relationships and are independent of docking scores or external annotations.

Once clustering is complete, the output includes a list of clusters, each with
its member compounds and a representative structure or scaffold. This
information can be used to select diverse compounds for experimental validation,
identify recurring chemotypes, or guide further molecular design. Visual
inspection of cluster representatives and distribution plots is supported within
Scipion-Chem, and all data can be exported for reporting or use in other
protocols.

In summary, the ClusterSubstructures protocol offers a practical and automated
way to group chemical compounds based on their internal structure. It supports
exploratory analysis of molecular libraries and enhances interpretability in
ligand-based screening workflows.

    """
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
