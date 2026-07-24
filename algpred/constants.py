# -*- coding: utf-8 -*-
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

DEFAULT_VERSION = '1.4'

# Unlike BepiPred/NetMHCpan/NetMHCIIpan/SignalP, AlgPred2 is a plain pip
# package published on PyPI by the Raghava group (GPL-3.0, no academic
# license gate) -- it is installed automatically at plugin-install time,
# following the same pattern as EpiDope.
ALGPRED_DIC = {
    'name': 'AlgPred2',
    'version': DEFAULT_VERSION,
    'home': 'ALGPRED_HOME',
    'activation': 'ALGPRED_ACTIVATION_CMD',
}

READ_URL = 'https://github.com/Lvera-code/scipion-chem-algpred'

NOINSTALL_WARNING = (
    'Installation could not be completed because the AlgPred2 conda '
    'environment was not found or its activation failed. Please check the '
    f'scipion-chem-algpred README file for more details: {READ_URL}'
)

# AlgPred2's own ML_Score / sklearn-based classifier throws
# "ValueError: Expected 2D array, got 1D array" when the input FASTA has
# exactly 1 sequence (a real upstream reshape bug, not a hypothetical edge
# case: the construct-level check always submits exactly 1 sequence).
# Workaround implemented in ProtAlgPred2Prediction._runAlgPred2: duplicate
# the single sequence before submission, then drop the extra row from the
# parsed result.
