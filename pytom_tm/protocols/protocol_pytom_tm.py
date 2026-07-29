# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import logging
from enum import Enum
from typing import List

from pwem.objects import VolumeMask
from pwem.protocols import EMProtocol
from pytom_tm.constants import IN_TOMOS, REF_VOL, IN_MASK
from pytom_tm.objects import SetOfPytomScoreTomograms
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam
from pyworkflow.utils import Message

logger = logging.getLogger(__name__)


class Pytom_tm_output(Enum):
    scoreTomograms = SetOfPytomScoreTomograms  # output name is scoreTomograms within a specific class SetOfPytomTOmograms


class ProtPytomTemplateMatching(EMProtocol):
    """GPU-accelerated template matching for cryo-electron tomography,
    originally developed in PyTom, as a standalone Python package that is run from the command line.
    """

    _label = 'template matching'
    _devStatus = BETA
    _possibleOutputs = Pytom_tm_output

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')

        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Reference volume")

        group = form.addGroup('Mask')
        group.addParam(IN_MASK, PointerParam,
                       pointerClass=VolumeMask,
                       important=True,
                       label='Reference mask')

        group.addParam('non_spherical_mask', BooleanParam,
                       label='Is it a spherical mask?',
                       default=False,
                       help="Flag to set when the mask is not spherical. It adds the required "
                            "computations for non-spherical masks and roughly doubles computation time."
                       )

        form.addSection(label='Angular Search')
        form.addParam('particle_diameter',FloatParam,
                      label='Particle Diameter (Angstroms)',
                      allowsNull=False,
                      help="Provide a particle diameter (in Angstrom) to automatically determine the "
                           "angular sampling using the Crowther criterion. For the max resolution, "
                           "(2 * pixel size) is used unless a low-pass filter is specified, "
                           "in which case the low-pass resolution is used. For non-globular "
                           "macromolecules choose the diameter along the longest axis."
                      )

        # --------------------------- INFO functions ------------------------------
    def _validate(self) -> List[str]:
        valMsg = []

        if self.particle_diameter.get() < 0:
            valMsg.append('Particle diameter must be >= 0')

        return valMsg