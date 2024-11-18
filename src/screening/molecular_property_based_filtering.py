import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

from rdkit.Chem import RDConfig
import sys
import os
import abc
from tqdm import tqdm

from typing import Dict, Any

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

sys.path.append(os.path.join(RDConfig.RDContribDir, "NP_Score"))
import npscorer
import csv
from pathlib import Path
import pandas as pd


class MolecularFilter(abc.ABC):
    name: str

    @abc.abstractmethod
    def apply(self, *args, **kwargs) -> bool:
        raise NotImplementedError("Subclasses should implement this method")


class LipinskiRuleOf5(MolecularFilter):
    name = "Lipinski_Rule_of_5"

    def apply(
        self,
        molecular_weight: float,
        logp: float,
        h_bond_donor: int,
        h_bond_acceptors: int,
        rotatable_bonds: int,
    ) -> bool:
        return (
            molecular_weight <= 500
            and logp <= 5
            and h_bond_donor <= 5
            and h_bond_acceptors <= 5
            and rotatable_bonds <= 5
        )


class GhoseFilter(MolecularFilter):
    name = "Ghose_Filter"

    def apply(
        self,
        molecular_weight: float,
        logp: float,
        number_of_atoms: int,
        molar_refractivity: float,
    ) -> bool:
        return (
            molecular_weight >= 160
            and molecular_weight <= 480
            and logp >= 0.4
            and logp <= 5.6
            and number_of_atoms >= 20
            and number_of_atoms <= 70
            and molar_refractivity >= 40
            and molar_refractivity <= 130
        )


class VeberFilter(MolecularFilter):
    name = "Veber_Filter"

    def apply(
        self, rotatable_bonds: int, topological_surface_area_mapping: float
    ) -> bool:
        return rotatable_bonds <= 10 and topological_surface_area_mapping <= 140


class RuleOf3(MolecularFilter):
    name = "Rule_of_3_Filter"

    def apply(
        self,
        molecular_weight: float,
        logp: float,
        h_bond_donor: int,
        h_bond_acceptors: int,
        rotatable_bonds: int,
    ) -> bool:
        return (
            molecular_weight <= 300
            and logp <= 3
            and h_bond_donor <= 3
            and h_bond_acceptors <= 3
            and rotatable_bonds <= 3
        )


class ReosFilter(MolecularFilter):
    name = "REOS_Filter"

    def apply(
        self,
        molecular_weight: float,
        logp: float,
        h_bond_donor: int,
        h_bond_acceptors: int,
        formal_charge: int,
        rotatable_bonds: int,
        heavy_atoms: int,
    ) -> bool:
        return (
            molecular_weight >= 200
            and molecular_weight <= 500
            and logp >= -5
            and logp <= 5
            and h_bond_donor >= 0
            and h_bond_donor <= 5
            and h_bond_acceptors >= 0
            and h_bond_acceptors <= 10
            and formal_charge >= -2
            and formal_charge <= 2
            and rotatable_bonds >= 0
            and rotatable_bonds <= 8
            and heavy_atoms >= 15
            and heavy_atoms <= 50
        )


class DrugLikeFilter(MolecularFilter):
    name = "Drug_like_Filter"

    def apply(
        self,
        molecular_weight: float,
        num_of_rings: int,
        rotatable_bonds: int,
        h_bond_donor: int,
        h_bond_acceptors: int,
        logp: float,
    ) -> bool:
        return (
            molecular_weight < 400
            and num_of_rings > 0
            and rotatable_bonds < 5
            and h_bond_donor <= 5
            and h_bond_acceptors <= 10
            and logp < 5
        )


class QEDFilter(MolecularFilter):
    name = "QED_Filter"

    def apply(self, molecule: rdkit.Chem.Mol) -> float:
        return Chem.QED.qed(molecule)


class TPSAFilter(MolecularFilter):
    name = "TPSA_Filter"

    def apply(self, molecule: rdkit.Chem.Mol) -> float:
        return Descriptors.TPSA(molecule)


class SAScoreFilter(MolecularFilter):
    name = "SA_Score_Filter"

    def apply(self, molecule: rdkit.Chem.Mol) -> float:
        return sascorer.calculateScore(molecule)


class NPScoreFilter(MolecularFilter):
    name = "NP_Score_Filter"

    def apply(self, molecule: rdkit.Chem.Mol) -> float:
        fscore = npscorer.readNPModel()
        return npscorer.scoreMol(molecule, fscore)


class ADMETScoresFilter(MolecularFilter):
    name = "ADMET_Scores_Filter"

    def apply(self, molecular_weight: float, logp: float) -> float:
        molecular_weight_term = 0 if molecular_weight < 330 else 330 - molecular_weight
        admet_score = (2.5 - logp) + molecular_weight_term

        return admet_score


class MolecularPropertyCalculator:
    def __init__(self, filters_to_include: list[str] = None):
        all_filters = {
            LipinskiRuleOf5.name: LipinskiRuleOf5(),
            GhoseFilter.name: GhoseFilter(),
            VeberFilter.name: VeberFilter(),
            RuleOf3.name: RuleOf3(),
            ReosFilter.name: ReosFilter(),
            DrugLikeFilter.name: DrugLikeFilter(),
            QEDFilter.name: QEDFilter(),
            TPSAFilter.name: TPSAFilter(),
            SAScoreFilter.name: SAScoreFilter(),
            NPScoreFilter.name: NPScoreFilter(),
            ADMETScoresFilter.name: ADMETScoresFilter(),
        }
        if filters_to_include is None:
            self.filters = all_filters
        else:
            self.filters = {
                name: all_filters[name]
                for name in filters_to_include
                if name in all_filters
            }

    def calculate_basic_filters(
        self, molecule: None | rdkit.Chem.Mol
    ) -> Dict[str, Any]:
        results: Dict[str, Any] = {name: None for name in self.filters.keys()}
        results["Errors"] = []

        if molecule is None:
            results["Errors"].append("Molecule is None")
            return results

        try:
            molecular_weight: float = Descriptors.ExactMolWt(molecule)
            logp: float = Descriptors.MolLogP(molecule)
            h_bond_donor: int = Descriptors.NumHDonors(molecule)
            h_bond_acceptors: int = Descriptors.NumHAcceptors(molecule)
            rotatable_bonds: int = Descriptors.NumRotatableBonds(molecule)
            number_of_atoms: int = molecule.GetNumAtoms()
            molar_refractivity: float = Chem.Crippen.MolMR(molecule)
            topological_surface_area_mapping: float = Chem.QED.properties(molecule).PSA
            formal_charge: int = Chem.rdmolops.GetFormalCharge(molecule)
            heavy_atoms: int = molecule.GetNumHeavyAtoms()
            num_of_rings: int = Chem.rdMolDescriptors.CalcNumRings(molecule)
        except Exception as e:
            results["Errors"].append(str(e))
            return results

        try:
            if LipinskiRuleOf5.name in self.filters:
                results[LipinskiRuleOf5.name] = self.filters[
                    LipinskiRuleOf5.name
                ].apply(
                    molecular_weight,
                    logp,
                    h_bond_donor,
                    h_bond_acceptors,
                    rotatable_bonds,
                )
            if GhoseFilter.name in self.filters:
                results[GhoseFilter.name] = self.filters[GhoseFilter.name].apply(
                    molecular_weight, logp, number_of_atoms, molar_refractivity
                )
            if VeberFilter.name in self.filters:
                results[VeberFilter.name] = self.filters[VeberFilter.name].apply(
                    rotatable_bonds, topological_surface_area_mapping
                )
            if RuleOf3.name in self.filters:
                results[RuleOf3.name] = self.filters[RuleOf3.name].apply(
                    molecular_weight,
                    logp,
                    h_bond_donor,
                    h_bond_acceptors,
                    rotatable_bonds,
                )
            if ReosFilter.name in self.filters:
                results[ReosFilter.name] = self.filters[ReosFilter.name].apply(
                    molecular_weight,
                    logp,
                    h_bond_donor,
                    h_bond_acceptors,
                    formal_charge,
                    rotatable_bonds,
                    heavy_atoms,
                )
            if DrugLikeFilter.name in self.filters:
                results[DrugLikeFilter.name] = self.filters[DrugLikeFilter.name].apply(
                    molecular_weight,
                    num_of_rings,
                    rotatable_bonds,
                    h_bond_donor,
                    h_bond_acceptors,
                    logp,
                )
        except Exception as e:
            results["Errors"].append(str(e))

        if QEDFilter.name in self.filters:
            try:
                results[QEDFilter.name] = self.filters[QEDFilter.name].apply(molecule)
            except Exception as e:
                results["Errors"].append(str(e))

        if TPSAFilter.name in self.filters:
            try:
                results[TPSAFilter.name] = self.filters[TPSAFilter.name].apply(molecule)
            except Exception as e:
                results["Errors"].append(str(e))

        if SAScoreFilter.name in self.filters:
            try:
                results[SAScoreFilter.name] = self.filters[SAScoreFilter.name].apply(
                    molecule
                )
            except Exception as e:
                results["Errors"].append(str(e))

        if NPScoreFilter.name in self.filters:
            try:
                results[NPScoreFilter.name] = self.filters[NPScoreFilter.name].apply(
                    molecule
                )
            except Exception as e:
                results["Errors"].append(str(e))

        if ADMETScoresFilter.name in self.filters:
            try:
                results[ADMETScoresFilter.name] = self.filters[
                    ADMETScoresFilter.name
                ].apply(molecular_weight, logp)
            except Exception as e:

                results["Errors"].append(str(e))

        return results

    def get_result_keys(self) -> list:
        return list(self.filters.keys()) + ["Errors"]

    def compute_filtering_results_and_save_to_csv(
        self, supplier: Chem.SDMolSupplier, output_csv_path: Path
    ):

        if output_csv_path.exists():
            computed_chembl_ids = set(
                pd.read_csv(output_csv_path, sep=",", index_col=0).index
            )
            mode = "a"
        else:
            computed_chembl_ids = set()
            mode = "w"

        with open(output_csv_path, mode=mode, newline="") as csv_file:
            fieldnames = ["ChemblID"] + self.get_result_keys()
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            if mode == "w":
                writer.writeheader()

            progres_bar = tqdm(supplier, desc="Processing molecules")
            for molecule in progres_bar:
                if (
                    molecule is None
                    or molecule.GetProp("chembl_id") in computed_chembl_ids
                ):
                    continue
                chembl_id = molecule.GetProp("chembl_id")

                results = self.calculate_basic_filters(molecule)
                results["ChemblID"] = chembl_id
                writer.writerow(results)
