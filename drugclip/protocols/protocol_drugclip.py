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
import json
import os, csv

import pyworkflow.protocol.params as params
from drugclip import DRUGCLIP_DIC
from pwem.protocols import EMProtocol
from pyworkflow.object import Float

from pwchem import Plugin
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

RDKIT, OPENBABEL = 0, 1



class ProtDrugclip(EMProtocol):
    """
    Protocol to use DrugCLIP.

    AI Generated:

        ProtDrugclip - User Manual

        Overview
        --------
        This protocol predicts the binding affinity of small molecules to protein
        pockets (ROIs) using the DrugCLIP deep learning framework. It converts input
        molecules to SMILES, prepares them together with the selected pockets, and
        computes predicted binding scores for each molecule-pocket pair.

        Inputs
        ------
        - **inputStructROIs**: SetOfStructROIs object containing the regions of interest (ROIs)
          on the target protein.
        - **molecules**: SetOfSmallMolecules object representing the ligands to be evaluated.
        - **useManager**: Choose whether to manage chemical structures using RDKit
          or OpenBabel (for SMILES conversion).
        - **batchSize**: Number of inputSmallMolecules processed in each batch.
        - **maxPocketAtoms**: Maximum number of atoms allowed per pocket.

        Workflow
        --------
        1. **SMILES extraction**:
           - Converts molecule files to SMILES format using RDKit or OpenBabel.
           - SMILES strings are written to a file and mapped to original molecule files.

        2. **LMDB creation**:
           - Converts pockets and molecules into an LMDB database suitable for DrugCLIP.
           - Handles batch processing according to `batchSize` and `maxPocketAtoms`.

        3. **DrugCLIP execution**:
           - Runs the DrugCLIP model for each pocket LMDB.
           - Uses GPU if enabled.
           - Produces predicted binding scores for each molecule-pocket pair.

        4. **Output aggregation**:
           - Collects scores for all molecule-pocket combinations.
           - Creates a `results.csv` file with rows for pockets and columns for molecules.
           - Updates the SetOfStructROIs with DrugCLIP scores file and stores the output SQLite database.

        Outputs
        -------
        - **SetOfStructROIs**: Updated ROI set with DrugCLIP scores file attached to each pocket.
        - **results.csv**: CSV file containing predicted binding scores for each pocket?molecule pair.

        Practical Recommendations
        -------------------------
        - Use for evaluating binding potential of a defined set of ligands to specific protein pockets.
        - When processing large sets of molecules, ensure sufficient GPU/CPU resources.
        - Verify input ROIs are correctly defined and represent biologically relevant binding sites.

        Summary & Interpretation
        ------------------------
        - The `results.csv` contains predicted binding scores in floating point for
          each pocket-molecule pair.
        - Scores can be used to rank ligands or guide further docking and experimental studies.

        Warnings
        --------
        - Batch processing is limited by `batchSize`; large values may exceed GPU memory.
        - Pockets with more than `maxPocketAtoms` atoms will be truncated.
        - Ensure molecules are valid and convertible to SMILES; otherwise, they will be skipped.

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
        form.addParam('inputStructROIs', params.PointerParam,
                      pointerClass='SetOfStructROIs', allowsNull=True, label="Input ROIs: ",
                      help='Select the input ROIs.')
        form.addParam('indivPocket', params.StringParam, default='', label='Specific ROI: ',
                      help='Select the input ROI.')

        iGroup = form.addGroup('Input Ligands')
        iGroup.addParam('useLibrary', params.BooleanParam, label='Use library as input : ', default=False,
                        help='Whether to use a SMI library SmallMoleculesLibrary object as input')

        iGroup.addParam('inputLibrary', params.PointerParam, pointerClass="SmallMoleculesLibrary",
                        label='Input library: ', condition='useLibrary',
                        help="Input Small molecules library to predict")
        iGroup.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                        label='Input small molecules: ', condition='not useLibrary',
                        help='Set of small molecules to input the model for predicting their interactions')


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
        if not self.useLibrary.get():
            self._insertFunctionStep(self.convertSMIStep)
        self._insertFunctionStep(self.createLMDBStep)
        self._insertFunctionStep(self.runDrugclipStep)
        self._insertFunctionStep(self.createOutputFile)
        self._insertFunctionStep(self.createOutputStep)

    def convertSMIStep(self):
        smiDir = self.getInputSMIDir()
        if not os.path.exists(smiDir):
            os.makedirs(smiDir)

        molDir = self.copyInputMolsInDir()
        args = ' --multiFiles -iD "{}" --pattern "{}" -of smi --outputDir "{}"'. \
            format(molDir, '*', smiDir)
        Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=smiDir)

    def createLMDBStep(self):
        if self.useLibrary.get():
            smilesFile = os.path.abspath(self.inputLibrary.get().getFileName())
            iArgs = f'--smiles-file {smilesFile}'
        else:
            iArgs = f'--smiles-dir {self.getInputSMIDir()}'


        rois = self._getInpROIs()
        paths = [os.path.abspath(roi.getFileName()) for roi in rois]
        pocketFilesStr = ",".join(paths)

        outputDir = self.getLMDBDir()
        args = (f"{iArgs} "
            f"--pocket-files {pocketFilesStr} "
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
        lmdbDir = self.getLMDBDir()
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

            fullCommand = (
                    f"export CUDA_VISIBLE_DEVICES={self.gpuList.get()} && "
                    f"python {scriptPath} " + " ".join(args)
            )
            Plugin.runCondaCommand(
                self,
                args=[],
                condaDic=DRUGCLIP_DIC,
                program=f"bash -c '{fullCommand}'",
                cwd=self._getPath()
            )

    def createOutputFile(self):
        """ Creates a standard long-format results.csv (Pocket, Molecule, Score) """
        smiDic = self.getSMIdic()

        resultsDir = os.path.abspath(self._getPath('results'))
        lmdbDir = os.path.abspath(self._getPath('lmdb'))

        pocketFiles = [f for f in os.listdir(lmdbDir) if f.startswith("pocket") and f.endswith(".lmdb")]
        pockets = [os.path.splitext(f)[0] for f in pocketFiles]

        outputFile = os.path.join(self._getPath(), "results.csv")

        with open(outputFile, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["pocket", "molecule", "drugclip_score"])

            for pocket in pockets:
                pocketDir = os.path.join(resultsDir, pocket)
                scoreFile = next((os.path.join(pocketDir, f) for f in os.listdir(pocketDir) if f.endswith(".txt")),
                                 None)

                if not scoreFile:
                    continue

                with open(scoreFile) as sf:
                    for line in sf:
                        parts = line.strip().split("\t")
                        if len(parts) == 2:
                            smi, score = parts
                            # Map SMILES back to the original filename
                            molName = smiDic[smi]
                            writer.writerow([pocket, molName, score])

    def createOutputStep(self):
        resultsFile = self._getPath("results.csv")
        scoresJsonFile = self._getExtraPath("scoresFile.json")

        intDic = {}
        with open(resultsFile, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pock = row['pocket']
                mol = row['molecule']
                score = float(row['drugclip_score'])

                if pock not in intDic:
                    intDic[pock] = {}
                intDic[pock][mol] = {"DrugCLIP_score": score}

        inROIs = self._getInpROIs()
        outROIs = self.inputStructROIs.get().createCopy(outputPath=self._getPath())

        for roi in inROIs:
            outRoi = roi.clone()
            outROIs.append(outRoi)

        scoresJsonFile = self.writeInteractScoresDic(intDic) #todo the json file is being overwritten

        with open(scoresJsonFile, 'r') as f:
            combinedDic = json.load(f)


        outROIs.setInteractScoresFile(scoresJsonFile)
        outROIs.setInteractScoresDic(combinedDic)
        outROIs.updateScoreTypes()

        outMols = self.inputLibrary.get() if self.useLibrary.get() else self.inputSmallMolecules.get()

        outROIs.setInteractMols(mols=outMols)
        self._defineOutputs(outputStructROIs=outROIs)

        if not self.useLibrary.get():
            inMols =self.inputSmallMolecules.get()
            outputMols = inMols.createCopy(self._getPath(), copyInfo=True)

            for mol in inMols:
                nMol = mol.clone()
                molName = mol.getMolName()

                for roiIdx, roi in enumerate(inROIs):
                    roiName = os.path.splitext(os.path.basename(roi.getFileName()))[0]
                    scoreVal = intDic.get(roiName, {}).get(molName, {}).get("DrugCLIP_score", 0.0)

                    attrName = f'DrugCLIP_score_{roiName}' if len(inROIs) > 1 else 'DrugCLIP_score'
                    setattr(nMol, attrName, Float(scoreVal))

                outputMols.append(nMol)

            outputMols.updateMolClass()
            self._defineOutputs(outputSmallMolecules=outputMols)

        else:
            inLib = self.inputLibrary.get()
            mapDic = inLib.getLibraryMap(inverted=True, fullLine=True)
            oLibFile = self._getPath('outputLibrary.smi')

            roiNames = [os.path.splitext(os.path.basename(roi.getFileName()))[0] for roi in inROIs]

            newHeaders = []
            for name in roiNames:
                header = f"DrugCLIP_score_{name}" if len(roiNames) > 1 else "DrugCLIP_score"
                newHeaders.append(header)

            with open(oLibFile, 'w') as f:
                allMols = set()
                for rName in roiNames:
                    allMols.update(combinedDic.get(rName, {}).keys())

                for molName in allMols:
                    if molName in mapDic:
                        lineBase = mapDic[molName].strip()
                        scoresLine = []

                        for rName in roiNames:
                            valDic = combinedDic.get(rName, {}).get(molName, {})
                            s = valDic.get("DrugCLIP_score", 0.0)
                            scoresLine.append(str(s))

                        scoresStr = '\t'.join(scoresLine)
                        f.write(f"{lineBase}\t{scoresStr}\n")

            outputLib = inLib.clone()
            outputLib.setFileName(oLibFile)
            outputLib.setHeaders(inLib.getHeaders() + newHeaders)

            self._defineOutputs(outputLibrary=outputLib)


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
    def getInputSMIDir(self):
        return os.path.abspath(self._getExtraPath('inputSMI'))

    def getLMDBDir(self):
        return os.path.abspath(self._getPath("lmdb"))

    def copyInputMolsInDir(self):
        oDir = os.path.abspath(self._getTmpPath('inMols'))
        if not os.path.exists(oDir):
            os.makedirs(oDir)

        for mol in self.inputSmallMolecules.get():
            os.link(mol.getFileName(), os.path.join(oDir, os.path.split(mol.getFileName())[-1]))
        return oDir

    def getSMIdic(self):
        '''Returns a dictionary: {SMILES: molName}
        '''
        smiDic = {}
        if self.useLibrary.get():
            with open(self.inputLibrary.get().getFileName()) as f:
                for line in f:
                    smi, molName = line.strip().split()
                    smiDic[smi] = molName
        else:
            for smiFile in os.listdir(self.getInputSMIDir()):
                with open(os.path.join(self.getInputSMIDir(), smiFile)) as f:
                    for line in f:
                        smi, molName = line.strip().split()
                        smiDic[smi] = molName
        return smiDic

    def writeInteractScoresDic(self, intDic, outFile=None):
        if not outFile:
            outFile = os.path.join(self._getExtraPath(), 'scoresFile.json')

        finalData = {}

        rois = self.inputStructROIs.get()

        prevFile = rois.getInteractScoresFile()

        print(f'---i got: {prevFile}')

        if prevFile and os.path.exists(str(prevFile)):
            try:
                with open(str(prevFile), 'r') as f:
                    finalData = json.load(f)
            except Exception:
                finalData = {}

        print(f'----data before: {finalData}')

        for protID, newMols in intDic.items():
            if protID not in finalData:
                finalData[protID] = {}

            for molName, newScores in newMols.items():
                if molName in finalData[protID]:
                    finalData[protID][molName].update(newScores)
                else:
                    finalData[protID][molName] = newScores

        print(f'----data after: {finalData}')

        with open(outFile, 'w') as f:
            json.dump(finalData, f, indent=4)

        return outFile

    def getSpecifiedROIFile(self):
        myROI = None
        for roi in self.inputStructROIs.get():
          if roi.__str__() == self.indivPocket.get():
            myROI = roi.clone()
            break
        return myROI
        
    def _getInpROIs(self):
        inROIs = self.inputStructROIs.get()
        if self.indivPocket.get().strip():
            inROIs = [self.getSpecifiedROIFile()]
        return inROIs


    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if smi.lower() != 'smiles':
                    break
        return smi