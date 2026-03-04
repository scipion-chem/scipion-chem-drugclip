# **************************************************************************
# *
# * Authors:  Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from .constants import *
from pwchem.constants import RDKIT_DIC

_references = ['']


class Plugin(pwchemPlugin):
    @classmethod
    def defineBinaries(cls, env):
        cls.addDrugclipPackage(env)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(DRUGCLIP_DIC['home'], cls.getEnvName(DRUGCLIP_DIC))

    @classmethod
    def addDrugclipPackage(cls, env, default=True):
        installer = InstallHelper(
            DRUGCLIP_DIC['name'],
            packageHome=cls.getVar(DRUGCLIP_DIC['home']),
            packageVersion=DRUGCLIP_DIC['version']
        )

        installer.getCondaEnvCommand(
            DRUGCLIP_DIC['name'],
            binaryVersion=DRUGCLIP_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(DRUGCLIP_DIC)} && "
            "pip install -U torch torchvision torchaudio ipython tqdm rdkit==2022.9.3 lmdb biopython 'numpy<2'"
        ).addCommand(
            f"{cls.getEnvActivationCommand(DRUGCLIP_DIC)} && "
            "git clone https://github.com/bowen-gao/DrugCLIP.git",
            f"drugclip_cloned"
        ).addCommand(
            f"{cls.getEnvActivationCommand(DRUGCLIP_DIC)} && "
            "git clone https://github.com/dptech-corp/Uni-Core.git && "
            "cd Uni-Core && pip install -e . --no-build-isolation",
            f"unicore_installed"
        ).addCommand(
            f"{cls.getEnvActivationCommand(DRUGCLIP_DIC)} && "
            "git clone https://github.com/deepmodeling/Uni-Mol.git && "
            "cd Uni-Mol/unimol && "
            "pip install -e .",
            f"unimol_installed"
        ).addCommand(
            f"{cls.getEnvActivationCommand(DRUGCLIP_DIC)} && "
            "pip install -U gdown && "
            f"gdown https://drive.google.com/uc?id=1i87thnbNk8qeLF_tLx_BzelTukWbHaTR "
            f"-O {cls.getVar(DRUGCLIP_DIC['home'])}/DrugCLIP/checkpoint_best.pt",
            "checkpoint_downloaded"
        )

        installer.addPackage(
            env,
            dependencies=['git', 'wget', 'make', 'g++'],
            default=default
        )






