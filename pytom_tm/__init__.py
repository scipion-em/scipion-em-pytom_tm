# -*- coding: utf-8 -*-
import pwem

__version__ = '3.0.0'

from pytom_tm.constants import V0_13_2, PYTOM_TM_ENV_ACTIVATION, PYTOM_TM_DEFAULT_ACTIVATION_CMD, PYTOM_TM, \
    PYTOM_TM_DEFAULT_VERSION, PYTOM_TM_ENV_NAME
from pyworkflow import TOMO

_logo = "icon.jpeg"


# _references = ['']

class Plugin(pwem.Plugin):
    _supportedVersions = [V0_13_2]
    _url = "https://github.com/scipion-em/scipion-em-pytom_tm"
    _processingField = [TOMO]

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(PYTOM_TM_ENV_ACTIVATION, PYTOM_TM_DEFAULT_ACTIVATION_CMD) #association key-value

    @classmethod
    def getPyTomEnvActivation(cls):
        return cls.getVar(PYTOM_TM_ENV_ACTIVATION)

    @classmethod
    def defineBinaries(cls, env):
        PYTOM_TM_INSTALLED = f'{PYTOM_TM}_{PYTOM_TM_DEFAULT_VERSION}_installed'
        installationCmd = cls.getCondaActivationCmd()
        # Create the environment
        installationCmd += f'conda create -n {PYTOM_TM_ENV_NAME} -c conda-forge -y python=3 cupy cuda-version=12.2 && '

        # Activate new the environment
        installationCmd += f'conda activate {PYTOM_TM_ENV_NAME} && '

        # Install fidder
        installationCmd += 'python -m pip install pytom-match-pick && '

        # Flag installation finished
        installationCmd += f'touch {PYTOM_TM_INSTALLED}'

        PYTOM_TM_commands = [(installationCmd, PYTOM_TM_INSTALLED)]
        #envPath = os.environ.get('PATH', "")  # keep path since conda likely in there
        #installEnvVars = {'PATH': envPath} if envPath else None

        env.addPackage(PYTOM_TM,
                       version=PYTOM_TM_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=PYTOM_TM_commands,
                       neededProgs=cls.getDependencies(), #check if dependencies are available before starting the installation
                       #vars=installEnvVars,
                       default=True)

    @classmethod
    def getDependencies(cls):
        # try to get CONDA activation command
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')
        return neededProgs

    @classmethod
    def runPytom(cls, protocol, program, args, cwd=None, numberOfMpi=1):
        """ Run pytom command from a given protocol. """
        cmd = cls.getCondaActivationCmd() + " "
        cmd += cls.getPyTomEnvActivation()
        cmd += f" && CUDA_VISIBLE_DEVICES=%(GPU)s {program} "

        protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd, numberOfMpi=numberOfMpi)

