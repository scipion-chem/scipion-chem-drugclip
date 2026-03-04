# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import os, csv
import pyworkflow.protocol.params as params
from drugclip import DRUGCLIP_DIC
from pwem.protocols import EMProtocol
from pyworkflow.object import String

from pwchem import Plugin
from pwchem.objects import  SetOfStructROIs, StructROI
from pwchem.utils import insistentRun
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

RDKIT, OPENBABEL = 0, 1



class ProtDrugclip(EMProtocol):
    """
    Protocol to use DrugCLIP.
    """
    _label = 'binding prediction'

    # -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addHidden('useGpu', params.BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. Choose one.")

        form.addHidden('gpuList', params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="Comma-separated GPU devices that can be used.")

        form.addSection(label='Input')
        form.addParam('pockets', params.PointerParam,
                      pointerClass='SetOfStructROIs', allowsNull=False,
                      label="ROIs: ",
                      help='Select the input ROIs.') #todo build one lmdb per pocket

        form.addParam('molecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Molecules: ",
                      help='Select the molecules to use.')

        form.addParam('useManager', params.EnumParam, default=1, label='Manage structure using: ',
                      choices=['RDKit', 'OpenBabel'],
                      help='Whether to manage the structure (conversion to SMILES) using RDKit or OpenBabel')

        group = form.addGroup('Parameters')
        group.addParam('batchSize', params.IntParam, default=8,
                       label='Batch size: ',
                       help='Number of molecules processed per batch.')
        group.addParam('maxPocketAtoms', params.IntParam, default=256,
                       label='Max. atoms per pocket: ',
                       help='Maximum number for atoms per pocket.')


        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.getSmilesStep)
        self._insertFunctionStep(self.convertFilesStep)
        self._insertFunctionStep(self.runDrugclipStep)
        self._insertFunctionStep(self.createOutputStep)

    def getSmilesStep(self):
        outputFile = self._getPath('smiles.txt')
        self.smi_to_file = {}  # NEW: map SMILES ? pose file

        with open(outputFile, 'w') as out:
            for mol in self.molecules.get():

                molFile = (os.path.abspath(mol.getPoseFile())
                           if mol.getPoseFile()
                           else os.path.abspath(mol.getFileName()))

                smi = self.getSMI(molFile)

                if smi:
                    out.write(smi + "\n")
                    self.smi_to_file[smi] = os.path.basename(molFile)  # store only the filename
                else:
                    print(f"Failed to extract SMILES from {molFile}")

        print(f"SMILES written to {outputFile}")

    def convertFilesStep(self):
        smilesFile = os.path.abspath(self._getPath("smiles.txt"))

        pocketFiles = ",".join(
            [os.path.abspath(roi.getFileName()) for roi in self.pockets.get()]
        )

        outputDir = os.path.abspath(self._getPath("lmdb"))

        args = (f"--smiles-file {smilesFile} "
            f"--pocket-files {pocketFiles} "
            f"--output-dir {outputDir} "
            f"--max-pocket-atoms {self.maxPocketAtoms.get()}"
        )
        scriptPath = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "scripts", "create_lmdb.py")
        )

        Plugin.runCondaCommand(
            self,
            args=args,
            condaDic=DRUGCLIP_DIC,
            program=f"python {scriptPath}",
            cwd=self._getPath()
        )

    def runDrugclipStep(self):
        weightPath = os.path.abspath(
            os.path.join(Plugin.getVar(DRUGCLIP_DIC['home']), 'DrugCLIP/checkpoint_best.pt')
        )
        lmdbDir = os.path.abspath(self._getPath('lmdb'))
        print(f"LMDB dir: {lmdbDir}")
        resultsDir = self._getPath('results')
        os.makedirs(resultsDir, exist_ok=True)
        pocketLmdbFiles = [
            os.path.join(lmdbDir, f)
            for f in os.listdir(lmdbDir)
            if f.startswith("pocket") and f.endswith(".lmdb")
        ]
        scriptPath = os.path.abspath(
            os.path.join(Plugin.getVar(DRUGCLIP_DIC["home"]), 'DrugCLIP/unimol/retrieval.py')
        )
        for pocketPath in pocketLmdbFiles:
            pocketName = os.path.splitext(os.path.basename(pocketPath))[0]

            args = [
                f"--user-dir {os.path.join(Plugin.getVar(DRUGCLIP_DIC['home']), 'DrugCLIP/unimol')}",
                "--valid-subset test",
                f"--results-path {os.path.abspath(resultsDir)}",
                f"--num-workers {self.numberOfThreads.get()}",
                "--ddp-backend c10d",
                f"--batch-size {self.batchSize.get()}",
                "--task drugclip",
                "--loss in_batch_softmax",
                "--arch drugclip",
                f"--max-pocket-atoms {self.maxPocketAtoms.get()}",
                "--fp16",
                "--fp16-init-scale 4",
                "--fp16-scale-window 256",
                "--seed 1",
                f"--path {weightPath}",
                "--log-interval 100",
                "--log-format simple",
                f"--mol-path {os.path.join(lmdbDir, 'mols.lmdb')}",
                f"--pocket-path {os.path.abspath(pocketPath)}",
                f"--emb-dir {os.path.abspath(resultsDir)}/{pocketName}",
                f"{os.path.abspath(os.path.join(Plugin.getVar(DRUGCLIP_DIC['home']), 'DrugCLIP/data'))}"
            ]

            full_command = (
                    f"export CUDA_VISIBLE_DEVICES={self.gpuList.get()} && "
                    f"python {scriptPath} " + " ".join(args)
            )
            Plugin.runCondaCommand(
                self,
                args=[],
                condaDic=DRUGCLIP_DIC,
                program=f"bash -c '{full_command}'",
                cwd=self._getPath()
            )

    def createOutputStep(self):
        resultsDir = os.path.abspath(self._getPath('results'))
        lmdbDir = os.path.abspath(self._getPath('lmdb'))

        pocketFiles = [f for f in os.listdir(lmdbDir) if f.startswith("pocket") and f.endswith(".lmdb")]
        pockets = [os.path.splitext(f)[0] for f in pocketFiles]

        allMolecules = set()
        pocketScoresDict = {}

        for pocket in pockets:
            pocketDir = os.path.join(resultsDir, pocket)
            scoreFile = None
            for f in os.listdir(pocketDir):
                if f.endswith(".txt"):
                    scoreFile = os.path.join(pocketDir, f)
                    break
            if not scoreFile:
                print(f"No score file found for pocket {pocket}")
                continue

            pocketScores = {}
            with open(scoreFile) as sf:
                for line in sf:
                    parts = line.strip().split("\t")
                    if len(parts) != 2:
                        continue
                    smi, score = parts
                    try:
                        pocketScores[smi] = float(score)
                        allMolecules.add(smi)
                    except ValueError:
                        print(f"Invalid score for molecule {smi} in pocket {pocket}")

            pocketScoresDict[pocket] = pocketScores

        allMolecules = sorted(allMolecules)
        allMoleculeFiles = [self.smi_to_file.get(smi, smi) for smi in allMolecules]

        outputFile = os.path.join(self._getPath(), "results.csv")
        with open(outputFile, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Pocket"] + allMoleculeFiles)
            for pocket in pockets:
                row = [pocket]
                scores = pocketScoresDict.get(pocket, {})
                row += [scores.get(smi, 0.0) for smi in allMolecules]
                writer.writerow(row)

        outSet = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        for pocket in self.pockets.get():
            outPock = StructROI()
            outPock.copy(pocket)
            outPock.Drugclip_file = String()
            outPock.setAttributeValue('Drugclip_file', str(outputFile))
            outSet.append(outPock)

        outSet.Drugclip_file = String()
        outSet.setAttributeValue('Drugclip_file', str(outputFile))

        outSet.buildPDBhetatmFile()
        self._defineOutputs(outputStructROIs=outSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["Results csv written in protocols path: results.csv"]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getSMI(self, fnSmall):
        fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))
        print("Extension:", ext)

        if ext != '.smi':
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.abspath(self._getExtraPath(fnRoot + '.smi'))

            args = f' -i "{fnSmall}" -of smi -o {fnOut} --outputDir {outDir}'

            if fnSmall.endswith(".pdbqt") or fnSmall.endswith(".mol2"):
                envDic, scriptName = OPENBABEL_DIC, 'obabel_IO.py'
            else:
                envDic, scriptName = RDKIT_DIC, 'rdkit_IO.py'

            fullProgram = (
                f'{Plugin.getEnvActivationCommand(envDic)} '
                f'&& python {Plugin.getScriptsDir(scriptName)} '
            )

            insistentRun(self, fullProgram, args, envDic=envDic, cwd=outDir)

            # ? CHECK FILE EXISTS
            if not os.path.exists(fnOut):
                print(f"SMILES conversion failed for {fnSmall}")
                return None

        else:
            fnOut = fnSmall

        return self.parseSMI(fnOut)

    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if smi.lower() != 'smiles':
                    break
        return smi