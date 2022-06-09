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
    cls._defineVar("RDKIT2_ENV_ACTIVATION", 'conda activate rdkit2-env')

  @classmethod
  def defineBinaries(cls, env):
      '''Download LePhar binaries'''
      LEPHAR_INSTALLED = 'lephar_installed'
      lephar_commands = ''
      if not os.path.exists(cls.getProgramHome(LEPHAR_DIC, path='bin')):
          lephar_commands += 'mkdir bin && '
      lephar_commands += 'cd bin && '
      for program in ['ledock', 'lepro', 'lewater', 'lefrag']:
          lephar_commands += 'wget {} && chmod +x {} && '.format(cls.getLePharUrl(program), cls.getProgramBin(program))

      for eBin in ['QueryDB']:
          lephar_commands += 'wget {} && unzip {}.zip && rm {}.zip && '.format(cls.getZipUrl(eBin), eBin, eBin)

      if not os.path.exists(cls.getProgramHome(LEPHAR_DIC, path='scripts')):
          lephar_commands += 'cd .. && mkdir scripts && cd scripts && '
      else:
          lephar_commands += 'cd ../scripts && '

      for script in ['ClusterByMCS']:
          lephar_commands += 'wget {} && unzip {}.zip && rm {}.zip && '.format(cls.getZipUrl(script), script, script)

      lephar_commands += 'conda create -y -n rdkit2-env rdkit=2018.09.3 python=2.7 && '
      lephar_commands += 'cd .. && touch ' + LEPHAR_INSTALLED
      lephar_commands = [(lephar_commands, LEPHAR_INSTALLED)]

      env.addPackage(LEPHAR_DIC['name'], version=LEPHAR_DIC['version'],
                     tar='void.tgz',
                     commands=lephar_commands,
                     default=True)

  @classmethod
  def runLePhar(cls, protocol, program, args, cwd=None, runJob=True, linuxSuf=True):
      """ Run LePhar command from a given protocol. """
      fullProgram = cls.getProgramHome(LEPHAR_DIC, path='bin/{}'.format(cls.getProgramBin(program, linuxSuf)))
      if runJob:
          protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
      else:
          fullProgram += args
          subprocess.call(fullProgram, shell=True, cwd=cwd)

  @classmethod
  def runRDKit2Script(cls, protocol, scriptName, args, cwd=None):
    """ Run rdkit command from a given protocol. """
    scriptPath = cls.getProgramHome(LEPHAR_DIC, path='scripts/{}/{}.pyo'.format(scriptName, scriptName))
    fullProgram = '%s %s && %s %s' % (cls.getCondaActivationCmd(), cls.getRDKit2EnvActivation(), 'python', scriptPath)
    protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

  @classmethod
  def getEnviron(cls):
    pass

  @classmethod
  def getRDKit2EnvActivation(cls):
    activation = cls.getVar("RDKIT2_ENV_ACTIVATION")
    return activation

  # ---------------------------------- Utils functions  -----------------------
  @classmethod
  def getLePharUrl(cls, program):
      '''Returns the url of the different lephar binaries'''
      return "http://www.lephar.com/download/{}_linux_x86".format(program)

  @classmethod
  def getZipUrl(cls, program):
      return 'http://www.lephar.com/download/{}.zip'.format(program)


  # In use paths
  @classmethod
  def getProgramHome(cls, programDic=None, path=''):
    if not programDic:
        programDic = LEPHAR_DIC
    return os.path.join(cls.getVar(programDic['home']), path)
  
  @classmethod
  def getProgramBin(cls, program, linuxSuf=True):
      if linuxSuf:
          return program + '_linux_x86'
      else:
          return program
