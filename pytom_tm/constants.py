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

#suffix
TOMO_SUFFIX = '_tomo'
MAS_SUFFIX = '_tomoMask'
DOSE_SUFFIX = '_dose'
