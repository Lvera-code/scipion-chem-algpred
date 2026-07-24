# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
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
This package contains a protocol for allergenicity prediction using a local
AlgPred2 installation.
"""

import subprocess

from pwchem import Plugin as pwchemPlugin
from scipion.install.funcs import InstallHelper

from .constants import ALGPRED_DIC, NOINSTALL_WARNING

_references = ['Sharma2021']
_logo = ''


class Plugin(pwchemPlugin):
    """AlgPred2 is a plain pip package (GPL-3.0, no academic license
    restriction), installed automatically in its own conda environment."""

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(ALGPRED_DIC['home'], cls.getEnvName(ALGPRED_DIC))
        cls._defineVar(ALGPRED_DIC['activation'], cls.getEnvActivationCommand(ALGPRED_DIC))

    @classmethod
    def defineBinaries(cls, env):
        cls.addAlgPred2Package(env)

    @classmethod
    def addAlgPred2Package(cls, env, default=True):
        installer = InstallHelper(ALGPRED_DIC['name'], packageHome=cls.getVar(ALGPRED_DIC['home']),
                                   packageVersion=ALGPRED_DIC['version'])
        installer.addCommand(
            f"conda create -y -p {cls.getVar(ALGPRED_DIC['home'])} python=3.10",
            'ALGPRED_ENV_CREATED'
        ).addCommand(
            f"{cls.getEnvActivationCommand(ALGPRED_DIC)} && pip install algpred2=={ALGPRED_DIC['version']}",
            'ALGPRED_INSTALLED'
        ).addCommand(
            f"touch {ALGPRED_DIC['name']}_installed",
            f"{ALGPRED_DIC['name']}_installed"
        ).addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def validateInstallation(cls):
        """Check that this plugin's requirements are met. Returns a list of
        actionable error messages, empty if the installation is correct."""
        errors = [] if cls.checkCallEnv(ALGPRED_DIC) else [NOINSTALL_WARNING]
        return errors

    @classmethod
    def checkCallEnv(cls, packageDic):
        actCommand = cls.getVar(packageDic['activation'])
        try:
            if 'conda' in actCommand and 'shell.bash hook' not in actCommand:
                actCommand = f'{cls.getCondaActivationCmd()}{actCommand}'
            subprocess.check_output(f'{actCommand} && algpred2 --help', shell=True)
            return True
        except subprocess.CalledProcessError:
            return False

    # ---------------------------------- Protocol functions-----------------------

    @classmethod
    def runAlgPred2(cls, protocol, args, cwd=None):
        """Run AlgPred2's console entry point from a given protocol."""
        activation = cls.getVar(ALGPRED_DIC['activation'])
        fullProgram = f'{activation} && algpred2'
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
