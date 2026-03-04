# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import setupTestProject, DataSet, BaseTest

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules, ProtDefineStructROIs
from pwem.protocols import ProtImportPdb
from pwchem.utils import assertHandle

from ..protocols import ProtDrugclip

defROIsStr = '''1) Coordinate: {"X": 45, "Y": 65, "Z": 60}
2) Residues: {"model": 0, "chain": "A", "index": "1-4", "residues": "VLSP"}
3) Ligand: {"molName": "HEM", "remove": "True"}
4) PPI: {"chain1": "0-A", "chain2": "0-B", "interDist": "5.0"}
5) Near_Residues: {"residues": "cys, cys", "distance": "5.0", "linkage": "Single"}'''

class TestDrugclip(BaseTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        setupTestProject(cls)

        cls._runImportSmallMols()
        cls._runImportPDB()
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protImportPDB, 'outputAtomStruct', sleepTime=5)


    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'))
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb, inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.proj.launchProtocol(cls.protImportPDB, wait=True)

    @classmethod
    def _runDefStructROIs(cls):
        cls.protDef = cls.newProtocol(
            ProtDefineStructROIs,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            inROIs=defROIsStr)
        cls.proj.launchProtocol(cls.protDef)

    def _runDrugclip(self):
        protDc = self.newProtocol(ProtDrugclip)

        protDc.pockets.set(self.protDef)
        protDc.pockets.setExtended('outputStructROIs')
        protDc.molecules.set(self.protImportSmallMols)
        protDc.molecules.setExtended('outputSmallMolecules')

        self.proj.launchProtocol(protDc, wait=False)
        return protDc

    def test(self):
        self._runDefStructROIs()
        self._waitOutput(self.protDef, 'outputStructROIs', sleepTime=10)

        protDc = self._runDrugclip()
        self._waitOutput(protDc, 'outputStructROIs', sleepTime=10)
        assertHandle(self.assertIsNotNone, getattr(protDc, 'outputStructROIs', None))
