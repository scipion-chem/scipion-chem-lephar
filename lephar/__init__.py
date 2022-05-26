# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

"""
This package contains the protocols for
manipulation of atomic struct objects
"""
import os, subprocess
import pwem

_logo = 'logo.jpg'
LEPHAR_DIC = {'name': 'lephar', 'version': '1.0', 'home': 'LEPHAR_HOME'}


class Plugin(pwem.Plugin):
  _homeVar = LEPHAR_DIC['home']
  _pathVars = [LEPHAR_DIC['home']]
  _supportedVersions = [LEPHAR_DIC['version']]

  @classmethod
  def _defineVariables(cls):
    """ Return and write a variable in the config file.
    """
    cls._defineEmVar(LEPHAR_DIC['home'], LEPHAR_DIC['name'] + '-' + LEPHAR_DIC['version'])

  @classmethod
  def defineBinaries(cls, env):
      '''Download LePhar binaries'''
      LEPHAR_INSTALLED = 'lephar_installed'
      lephar_commands = 'mkdir bin && cd bin && '
      for program in ['ledock', 'lepro', 'lewater', 'lefrag']:
          lephar_commands += 'wget {} && chmod +x {} && '.format(cls.getLePharUrl(program), cls.getProgramBin(program))

      lephar_commands += 'cd .. && touch ' + LEPHAR_INSTALLED
      lephar_commands = [(lephar_commands, LEPHAR_INSTALLED)]

      env.addPackage(LEPHAR_DIC['name'], version=LEPHAR_DIC['version'],
                     tar='void.tgz',
                     commands=lephar_commands,
                     default=True)

  @classmethod
  def runLePhar(cls, protocol, program, args, cwd=None, runJob=True):
      """ Run LePhar command from a given protocol. """
      fullProgram = cls.getProgramHome(LEPHAR_DIC, path='bin/{}'.format(cls.getProgramBin(program)))
      if runJob:
          protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
      else:
          fullProgram += args
          subprocess.call(fullProgram, shell=True, env=protocol._getEnviron(), cwd=cwd)

  @classmethod  # Test that
  def getEnviron(cls):
    pass

  # ---------------------------------- Utils functions  -----------------------
  @classmethod
  def getLePharUrl(cls, program):
      '''Returns the url of the different lephar binaries'''
      return "http://www.lephar.com/download/{}_linux_x86".format(program)

  # In use paths
  @classmethod
  def getProgramHome(cls, programDic, path=''):
    return os.path.join(cls.getVar(programDic['home']), path)
  
  @classmethod
  def getProgramBin(cls, program):
      return program + '_linux_x86'
