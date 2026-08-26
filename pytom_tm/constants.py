from gapstop.constants import IN_SCORE_TOMOS

PYTOM_TM = 'pytom_tm'
PYTOM_TM_HOME = 'PYTOM_TM_HOME'

# Supported versions
V0_13_2= '0.13.2'
PYTOM_TM_DEFAULT_VERSION = V0_13_2

PYTOM_TM_ENV_NAME = f'{PYTOM_TM}-{PYTOM_TM_DEFAULT_VERSION}'
PYTOM_TM_ENV_ACTIVATION = 'PYTOM_TM_ENV_ACTIVATION'
PYTOM_TM_DEFAULT_ACTIVATION_CMD = f'conda activate {PYTOM_TM_ENV_NAME}'

#Inputs
IN_TOMOS = 'inTomos'
IN_CTF_SET = 'inCtfSet'
IN_TS_SET = 'inTsSet'
REF_VOL = 'reference'
IN_MASK = 'mask'
TOMO_MASKS = 'volMask'

#Extensions
MRC_EXT = '.mrc'
DEFOCUS_EXT = '.defocus'
TILT_ANGLES_EXT = '.tlt'
DOSE_EXT = '.txt'
JSON_EXT = '.json'

#suffix
TOMO_SUFFIX = '_tomo'
MASK_SUFFIX = '_tomoMask'
DOSE_SUFFIX = '_dose'
SCORE_SUFFIX = '_tomo_scores'
JSON_SUFFIX = '_tomo_job'

#defocus handedness
DEFOCUS_HAND_NEG = -1
DEFOCUS_HAND_POS = 1
DEFOCUS_HAND_OFF = 0

#base seed
BASE_SEED=42

#---------- EXTRACT PARTICLES------------
IN_TM_PROTOCOL ='inTmProtocol'

#mask choice
MASK_NO = 0
MASK_PYTOM_TM = 1
MASK_OTHER = 2