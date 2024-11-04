import argparse
from utils.paths import CHEMBL_34_SDF_PATH, INITIAL_SCREENING_RESULTS_PATH
from screening.molecular_property_based_filtering import MolecularPropertyCalculator
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser(description="Drug Discovery Project")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument(
        "--num_threads", type=int, default=1, help="Number of writer threads"
    )

    args = parser.parse_args()

    if args.verbose:
        print(f"Input file: {args.input}")
        print(f"Output file: {args.output}")

    # Add your main logic here
    print("Processing...")

    mol_prop_calc = MolecularPropertyCalculator()

    supplier = Chem.MultithreadedSDMolSupplier(
        fileName=CHEMBL_34_SDF_PATH, numWriterThreads=args.num_threads
    )

    mol_prop_calc.compute_filtering_results_and_save_to_csv(
        supplier=supplier, output_csv_path=INITIAL_SCREENING_RESULTS_PATH
    )


if __name__ == "__main__":
    main()
