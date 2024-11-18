from pathlib import Path

DATA_PATH: Path = Path(__file__).parent.parent.parent / "data"

CHEMBL_34_SDF_PATH: Path = DATA_PATH / "chembl_34.sdf"
INITIAL_SCREENING_RESULTS_PATH: Path = DATA_PATH
JUST_ADMET_PATH: Path = DATA_PATH / "chembl_34_initial_screening_results_just_admet.csv"
