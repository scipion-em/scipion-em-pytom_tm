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
import traceback
from enum import Enum
from typing import List, Optional, Union
from pwem.emlib.image import ImageHandler
from pwem.objects import VolumeMask, Volume
from pwem.protocols import EMProtocol
from pytom_tm.constants import IN_TOMOS, REF_VOL, IN_MASK, VOL_MASK, IN_TS_SET, IN_CTF_SET
from pytom_tm.objects import SetOfPytomScoreTomograms
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, IntParam, StringParam, LEVEL_ADVANCED, EnumParam
from pyworkflow.utils import Message, cyanStr, makePath, redStr
from tomo.objects import SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries
from tomo.utils import getObjFromRelation, getCommonTsAndCtfElements, invertContrast, convertOrLink

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
    _program = 'pytom_match_template.py'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.ctfDict = None
        self.refName = None
        self.ih = ImageHandler()
        self.failedTsIds = []
        self.failedTsIdsStr = String()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')

        form.addParam(IN_CTF_SET, PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="CTF tomo series",
                      important=True,
                      allowsNull=True)

        form.addParam(IN_TS_SET, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Tilt-series (opt.)',
                      help='Used to get the tilt angles. If empty, the protocol will try to reach, via relations, '
                           'the tilt-series associated to the introduced CTFs.')

        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Reference volume")

        form.addParam('invert_contrast', BooleanParam,
                      default=True,
                      label='Invert reference contrast?',
                      important=True,
                      help='The contrast of the template has to be the same as of the tomogram. If the '
                           'tomogram has features in black (which is typically for cryoET) then the template '
                           'has to have the same representation. For example, Relion outputs inverted '
                           'contrast (features are white) and such maps have to be inverted prior running '
                           'the Pytom_tm.')

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
        form.addParam('particle_diameter', FloatParam,
                      label='Particle Diameter (Angstroms)',
                      allowsNull=False,
                      help="Provide a particle diameter (in Angstrom) to automatically determine the "
                           "angular sampling using the Crowther criterion. For the max resolution, "
                           "(2 * pixel size) is used unless a low-pass filter is specified, "
                           "in which case the low-pass resolution is used. For non-globular "
                           "macromolecules choose the diameter along the longest axis."
                      )

        form.addParam('angular_search', StringParam,
                      label='Angular Search',
                      help="This option overrides the angular search calculation from the particle "
                           "diameter. If given a float it will generate an angle list with healpix "
                           "for Z1 and X1 and linear search for Z2. The provided angle will be used "
                           "as the maximum for the "
                           "linear search and for the mean angle difference from healpix."
                           "Alternatively, a .txt file can be provided with three Euler angles "
                           "(in radians) per line that define the angular search. "
                           "Angle format is ZXZ anti-clockwise (see: "
                           "https://www.ccpem.ac.uk/user_help/rotation_conventions.php).")

        form.addParam('z_axis_rotational_symmetry', IntParam,
                      label='Z Axis Rotational Symmetry',
                      default=1,
                      allowsNull=False,
                      help="Integer value indicating the rotational symmetry of the template around "
                           "the z-axis. The length of the rotation search will be shortened through "
                           "division by this value. Only works for template symmetry around the z-axis.")

        form.addSection(label='Volume control')
        form.addParam('volume_split', StringParam,
                      label='Volume Split',
                      default='1 1 1',
                      help="Split the volume into smaller parts for the search, "
                           "can be relevant if the volume does not fit into GPU memory. "
                           "Format is x y z, e.g. --volume-split 1 2 1")

        form.addParam('defocus_handedness', EnumParam,
                      choices=['-1', '0', '1'],
                      display=EnumParam.DISPLAY_HLIST,
                      label='Defocus Handedness',
                      expertLevel=LEVEL_ADVANCED,
                      # condition='sum(map(int, volume_split)) > 3',  #  TODO
                      help="Specify the defocus handedness for defocus gradient correction of the "
                           "CTF in each subvolumes. The more subvolumes in x and z, "
                           "the finer the defocus gradient will be corrected, at the cost of "
                           "increased computing time. It will only have effect for very clean and "
                           "high-resolution data, such as isolated macromolecules. "
                           "A value of 0 means no defocus gradient correction (default), 1 means "
                           "correction assuming correct handedness (as specified in Pyle and "
                           "Zianetti (2021)), -1 means the handedness will be inverted. If uncertain "
                           "better to leave off as an inverted correction might hamper results."
                      )

        search_indices_group = form.addGroup('Search Indices')
        x_index = search_indices_group.addLine('X axis',
                                               help="Start and end indices of the search along the x-axis")

        x_index.addParam('xmin',
                         IntParam,
                         label='X min',
                         allowsNull=True)

        x_index.addParam('xmax',
                         IntParam,
                         label='X max',
                         allowsNull=True)

        y_index = search_indices_group.addLine('Y axis',
                                               help="Start and end indices of the search along the y-axis")

        y_index.addParam('ymin',
                         IntParam,
                         label='Y min',
                         allowsNull=True)

        y_index.addParam('ymax',
                         IntParam,
                         label='Y max',
                         allowsNull=True)

        z_index = search_indices_group.addLine('Z axis',
                                               help="Start and end indices of the search along the z-axis")

        z_index.addParam('zmin',
                         IntParam,
                         label='Z min',
                         allowsNull=True)

        z_index.addParam('zmax',
                         IntParam,
                         label='Z max',
                         allowsNull=True)
        # TODO review
        form.addParam(VOL_MASK, PointerParam,
                      pointerClass=VolumeMask,
                      label='Tomogram mask',
                      help="Here you can provide a mask for matching with dimensions (in pixels) "
                           "equal to the tomogram. If a subvolume only has values <= 0 for this mask it "
                           "will be skipped."
                      )

        form.addSection(label='Filter Control')
        form.addParam('low_pass', FloatParam,
                      label='Low Pass filter ',
                      allowsNull=False,
                      help="Apply a low-pass filter to the tomogram and template. Generally desired "
                           "if the template was already filtered to a certain resolution. "
                           "Value is the resolution in A."
                      )
        form.addParam('high_pass', FloatParam,
                      label='High Pass filter ',
                      allowsNull=False,
                      help="Apply a high-pass filter to the tomogram and template to reduce "
                           "correlation with large low frequency variations. Value is a resolution in A, "
                           "e.g. 500 could be appropriate as the CTF is often incorrectly modelled "
                           "up to 50nm."
                      )
        form.addParam('spectral_whitening', BooleanParam,
                      label='Spectral Whitening',
                      help="Calculate a whitening filtering from the power spectrum of the tomogram; "
                           "apply it to the tomogram patch and template. Effectively puts more weight on "
                           "high resolution features and sharpens the correlation peaks.")

        form.addSection(label='Additional Parameters')
        form.addParam('random_phase_correction', BooleanParam,
                      label='Random Phase Correction',
                      help="Run template matching simultaneously with a phase randomized version of "
                           "the template, and subtract this 'noise' map from the final score map. "
                           "For this method please see STOPGAP as a reference: "
                           "https://doi.org/10.1107/S205979832400295X ."
                      )
        form.addParam('extraParams', StringParam,
                      label='Additional parameters',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--half-precision \n"
                           "--warp-xml-file \n")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []

        cRId = self._insertFunctionStep(self.convertReferenceStep,
                                        prerequisites=[],
                                        needsGPU=False)

        for tsId in self.tomoDict.keys():
            cInputId = self._insertFunctionStep(self.convertInputStep, tsId,
                                                prerequisites=cRId,
                                                needsGPU=False)
            tmId = self._insertFunctionStep(self.templateMatchingStep, tsId,
                                            prerequisites=cInputId,
                                            needsGPU=True)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=tmId,
                                              needsGPU=False)
            closeSetStepDeps.append(cOutId)
        self._insertFunctionStep(self.closeOutputSetStep,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

        # -------------------------- STEPS functions ------------------------------

    def _initialize(self):
        tsSet = self._getTsSet()
        tomoSet = self._getFormAttrib(IN_TOMOS)
        ctfSet = self._getFormAttrib(IN_CTF_SET)
        # self.refName = self._genConvertedOrLinkedRefName(REF_VOL)
        # self.maskName = self._genConvertedOrLinkedRefName(IN_MASK)
        # self.tomosSRate = tomoSet.getSamplingRate()
        # self.tomosBinning = self._getTomogramsBinning()

        # Compute matching TS id among coordinates, the tilt-series and the CTFs, they all could be a subset
        tomosTsIds = set(tomoSet.getTSIds())
        tsIds = set(tsSet.getTSIds())
        ctfTsIds = set(ctfSet.getTSIds())
        presentTsIds = tomosTsIds & tsIds & ctfTsIds
        nonMatchingTsIds = tomosTsIds ^ tsIds ^ ctfTsIds

        # Validate the intersection
        if len(presentTsIds) <= 0:
            raise Exception("There isn't any common tilt-series ids among the coordinates, CTFs, and tilt-series "
                            "introduced.")

        if len(nonMatchingTsIds) > 0:
            logger.info(cyanStr(f"TsIds not common in the introduced tomograms, CTFs, and "
                                f"tilt-series are: {nonMatchingTsIds}"))

        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in tomoSet.iterItems()
                         if tomo.getTsId() in presentTsIds}
        self.tsDict = {ts.getTsId(): ts.clone() for ts in tsSet.iterItems()
                       if ts.getTsId() in presentTsIds}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in ctfSet.iterItems()
                        if ctf.getTsId() in presentTsIds}

    def convertReferenceStep(self):
        logger.info(cyanStr(f"Converting the reference in the required format...'"))
        try:
            # Convert or link the reference
            ref = self._getFormAttrib(REF_VOL)
            inRefFile = ref.getFileName()
            refFile = self.getReferenceFileName()
            samplingRate = ref.getSamplingRate()

            if self.invert_contrast.get():
                invertContrast(inRefFile, refFile, samplingRate)
            else:  # if the contrast inversion is not needed, only convert to .mrc
                convertOrLink(inRefFile, refFile, samplingRate)

            # Convert or link the mask
            mask = self._getFormAttrib(IN_MASK)
            inMaskFile = mask.getFileName()
            outMaskFile = self.getMaskFileName()
            samplingRate = mask.getSamplingRate()
            convertOrLink(inMaskFile, outMaskFile, samplingRate)
        except Exception as e:
            raise Exception(f'Reference conversion failed with the exception -> {e}')

    def convertInputStep(self, tsId: str):
        try:
            tomo = self.tomoDict[tsId]
            ts = self.tsDict[tsId]
            ctf = self.ctfDict[tsId]
            presentAcqOrders = getCommonTsAndCtfElements(ts, ctf)
            if len(presentAcqOrders) == 0:
                raise Exception(f'tsId = {tsId} -> No common acquisition orders found between the '
                                f'tilt-series and the CTF.')

            logger.info(cyanStr(f"tsId = {tsId} -> present acquisition orders in both "
                                f"the tilt-series and the CTF are {presentAcqOrders}.'"))

            tsDir = self._getCurrentTomoDir(tsId)
            makePath(tsDir)





        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def templateMatchingStep(self):
        pass

    # --------------------------- INFO functions ------------------------------

    def _validate(self) -> List[str]:
        valMsg = []

        if self.particle_diameter.get() < 0:
            valMsg.append('Particle diameter must be >= 0')

        if self.z_axis_rotational_symmetry.get() < 0:
            valMsg.append('Z axis rotational symmetry must be >= 0')

        if not self.validate_volume_split(self.volume_split.get()):
            valMsg.append('Volume split must be a list of three integers (e.g. 1 1 1)')

        return valMsg

        # --------------------------- UTILS functions ------------------------------

    @staticmethod
    def validate_volume_split(text: str) -> bool:
        parts = text.split()

        if len(parts) != 3:  # check if it's a three elements list
            return False
        try:
            [int(x) for x in parts]
            return True
        except ValueError:
            return False

    def _getTsSet(self) -> SetOfTiltSeries:
        tsSet = self._getFormAttrib(IN_TS_SET)
        return tsSet if tsSet else self._getTsFromRelations()

    def _getFormAttrib(self, attribName: str, returnPointer: bool = False) -> Optional[Union[SetOfTiltSeries,
    SetOfTomograms, SetOfCTFTomoSeries, Volume, Pointer]]:
        inTsPointer = getattr(self, attribName, None)
        if not inTsPointer:
            return None
        else:
            return inTsPointer if returnPointer else inTsPointer.get()

    # if IN_TS_SET is not an input, get it yourself via:
    def _getTsFromRelations(self) -> Optional[SetOfTiltSeries]:
        inCTFs = self._getFormAttrib(IN_CTF_SET)
        return getObjFromRelation(inCTFs, self, SetOfTiltSeries)

    def _getCurrentTomoDir(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def getReferenceFileName(self) -> str:
        return self._getTmpPath('Reference.mrc')

    def getMaskFileName(self) -> str:
        return self._getTmpPath('Mask.mrc')

    def _generateArguments(self) -> str:
        cmd = [
            f'--template {}'

        ]
        return ' '.join(cmd)
