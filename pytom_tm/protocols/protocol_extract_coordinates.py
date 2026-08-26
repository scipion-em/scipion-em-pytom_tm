import logging
from enum import Enum

from pwem.protocols import EMProtocol
from pytom_tm.constants import IN_TM_PROTOCOL, TOMO_MASKS, MASK_PYTOM_TM
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, IntParam, GT, FloatParam, GE, LE, StringParam, EnumParam
from pyworkflow.utils import Message
from tomo.objects import SetOfCoordinates3D, SetOfTomoMasks

logger = logging.getLogger(__name__)


class Pytom_extract_outputs(Enum):
    coordinates = SetOfCoordinates3D


class ProtPytomExtractCoordinates(EMProtocol):
    """Extract coordinates from Pytom score tomograms."""

    _label = 'extract coordinates from Pytom'
    _devStatus = BETA
    _possibleOutputs = Pytom_extract_outputs
    _program = 'pytom_extract_candidates.py'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomoscoreDict = {}
        self.failedTsIds = []
        self.failedTsIdsStr = String()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TM_PROTOCOL, PointerParam,
                      pointerClass='ProtPytomTemplateMatching',
                      important=True,
                      label='Protocol Pytom TM')
        form.addParam('mask_choice', EnumParam,
                      choices=['no TomoMasks', 'Same as input protocol', 'Other TomoMasks'],
                      display=EnumParam.DISPLAY_COMBO,
                      label='TomoMasks (segmentations)',
                      default=MASK_PYTOM_TM,
                      help="Here you can provide a mask for the extraction with dimensions "
                           "(in pixels) equal to the tomogram. All values in the mask that are smaller or "
                           "equal to 0 will be removed, all values larger than 0 are considered regions "
                           "of interest. It can be used to extract annotations only within a specific "
                           "cellular region. If the job was run with a tomogram mask, this file will be "
                           "used instead of the job mask."
                      )
        form.addParam(TOMO_MASKS, PointerParam,
                      pointerClass=SetOfTomoMasks,
                      label='Tomogram masks (segmentations, opt.)',
                      allowsNull=True,
                      help="All values in the mask that are smaller or "
                           "equal to 0 will be removed, all values larger than 0 are considered regions "
                           "of interest. It can be used to extract annotations only within a specific "
                           "cellular region."
                      )
        form.addParam('number_of_particles', IntParam,
                      label='Number of particles',
                      allowsNull=False,
                      validators=[GT(0)],
                      help="Maximum number of particles to extract from tomogram."
                      )

        form.addParam('number_false_positives', FloatParam,
                      label='Number of false positives',
                      important=True,
                      validators=[GE(0), LE(1)],
                      default=1.0,
                      help="Number of false positives to determine the false alarm rate. Here one "
                           "can increase the recall of the particle of interest at the expense "
                           "of more false positives. The default value of 1 is recommended for "
                           "particles that can be distinguished well from the background (high "
                           "specificity). The value can also be set between 0 and 1 to make "
                           "the cut-off more restrictive."  # TODO check description
                      )
        form.addParam('particle_diameter', FloatParam,
                      label='Particle diameter (angst)',
                      important=True,
                      validators=[GT(0)],
                      help="Particle diameter of the template in Angstrom. It is used during "
                           "extraction to remove areas around peaks to prevent double extraction. "
                           "Minimal peak-to-peak distance after extraction will be diameter/2."
                           "If not previously specified, this option is required. If "
                           "specified in pytom_match_template, this is optional and "
                           "can be used to overwrite it, which might be relevant for strongly "
                           "elongated particles--where the angular sampling should be "
                           "determined using its long axis but the extraction mask should use its "
                           "short axis."
                      )
        form.addParam('tophat_filter_con', IntParam,
                      label='Tophat filter connectivity',
                      validators=[GE(0)],
                      default=0,
                      help="Set kernel connectivity for ndimage binary structure used for the "
                           "tophat transform. Integer value in range 1-3. 1 is the most "
                           "restrictive, 3 the least restrictive. Generally recommended to "
                           "leave at 1. 0 value means that tophat filter (attempt to filter "
                           "only sharp correlation peaks with a tophat transform) is disabled.",
                      )
        form.addSection(label='Additional Parameters')
        form.addParam('extraParams', StringParam,
                      label='Additional parameters',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--cut-off \n"
                      )

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []

        cRId = self._insertFunctionStep(self.convertReferenceStep,
                                        prerequisites=[],
                                        needsGPU=False)

        for tsId in self.tomoDict.keys():
            cInputId = self._insertFunctionStep(self.convertTomoMaskStep, tsId,
                                                prerequisites=cRId,
                                                needsGPU=False)
            tmId = self._insertFunctionStep(self.extractCoordinatesStep, tsId,
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
