# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
from ..protocols import *
import pyworkflow.wizard as pwizard
from pwem.wizards import VariableWizard

class AddLePharFilter(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getOriginLabel(self, inParamName, protocol):
        if inParamName == 'descriptor':
            return 'Descriptor', protocol.getEnumText(inParamName)
        elif inParamName == 'fGroup':
            return 'Functional group', getattr(protocol, inParamName).get().strip().replace(' ', '')

    def show(self, form, *params):
      protocol = form.protocol
      inputParams, outputParam = self.getInputOutput(form)

      inList = getattr(protocol, outputParam[0]).get()
      num = len(inList.strip().split('\n'))
      if inList.strip() != '':
        num += 1

      inParamName = inputParams[0]
      label, inParamValue = self.getOriginLabel(inParamName, protocol)

      if inParamValue != '':
          descInfo = '{}\t{}\t{}'.format(inParamValue, getattr(protocol, inputParams[1]).get(),
                                         getattr(protocol, inputParams[2]).get())
          form.setVar(outputParam[0], inList + '{}) {}: {}\n'.format(num, label, descInfo))


AddLePharFilter().addTarget(protocol=ProtChemFilterLigands,
                            targets=['descriptor'],
                            inputs=['descriptor', 'minDesc', 'maxDesc'],
                            outputs=['filterList'])

AddLePharFilter().addTarget(protocol=ProtChemFilterLigands,
                            targets=['fGroup'],
                            inputs=['fGroup', 'minFG', 'maxFG'],
                            outputs=['filterList'])
