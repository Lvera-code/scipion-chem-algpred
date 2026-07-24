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
        # NOTE (real bug found+fixed 2026-07-24 via a live 'scipion3 installb'
        # run, never exercised before that): using InstallHelper with a custom
        # packageHome pointing outside Scipion's own default build directory
        # (EM_ROOT/<name>-<version>/) breaks env.addPackage's own completion
        # check, which looks for the target marker file in ITS OWN default
        # build dir regardless of what packageHome InstallHelper was given --
        # confirmed by reproducing 'ERROR: File ... not found' for real. Also,
        # a '-p <path>' conda env cannot be activated by
        # getEnvActivationCommand() (which always does 'conda activate
        # <envName>', a NAME lookup, inherited from pwchemPlugin) -- so this
        # now creates a NAMED env ('-n'), same convention as bepipred (the
        # other real, working example in this project), and lets
        # env.addPackage use its own default build directory instead of a
        # mismatched custom one.
        # NOTE (real bug found+fixed 2026-07-24, THREE more issues found only
        # by running an actual prediction, not just 'algpred2 --help'): the
        # real PyPI 'algpred2' package under-declares its own dependencies
        # ('pip show algpred2' lists only argparse/numpy/pandas) -- it also
        # needs joblib and scikit-learn at runtime to unpickle its bundled
        # RandomForestClassifier. Worse, the LATEST scikit-learn (1.7.2) is
        # binary-incompatible with that bundled pickle ('ValueError: node
        # array from the pickle has an incompatible dtype') since the model
        # was serialized against an older sklearn Tree format. Pinned here to
        # scikit-learn==1.2.2 / joblib==1.5.3, the exact versions confirmed
        # to load and predict correctly against the real bundled model (taken
        # from this project's own known-working .venv-algpred). A FOURTH bug
        # surfaced immediately after pinning just those two: an unpinned
        # 'pip install' also pulls the latest numpy (2.2.6), whose C-API/ABI
        # is incompatible with scikit-learn 1.2.2's compiled Cython
        # extensions (real error: 'ValueError: numpy.dtype size changed, may
        # indicate binary incompatibility' in sklearn.utils.murmurhash) --
        # numpy is pinned to 1.26.4 (also taken from .venv-algpred) to match
        # the ABI scikit-learn 1.2.2 was built against.
        envName = cls.getEnvName(ALGPRED_DIC)
        installedMarker = f"{ALGPRED_DIC['name']}_installed"
        installationCmd = (
            f"conda create -y -n {envName} python=3.10 && "
            f"{cls.getEnvActivationCommand(ALGPRED_DIC)} && "
            f"pip install algpred2=={ALGPRED_DIC['version']} joblib==1.5.3 "
            f"numpy==1.26.4 scikit-learn==1.2.2 && "
            f"touch {installedMarker}"
        )
        env.addPackage(ALGPRED_DIC['name'], version=ALGPRED_DIC['version'],
                       commands=[(installationCmd, installedMarker)], tar='void.tgz',
                       neededProgs=['conda'], default=default)

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
