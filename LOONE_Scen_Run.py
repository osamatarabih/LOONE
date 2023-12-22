import argparse
import pandas as pd
from LOONE_Nut import LOONE_Nut


def main(loone_q_path: str, nut_out_path: str, data_dir: str, ensemble_number: str) -> None:
    """Run LOONE_Nut_Lds_out.csv

    Args:
        loone_q_path (str): Path to LOONE Q CSV file.
        nut_out_path (str): Path to LOONE nutrient file to be created.
        data_dir (str): Path to data directory.
        ensemble_number (str): The ensemble number of the ensemble that the input data is for. Ex: '1'.
    """
    ensemble_number = int(ensemble_number)
    LOONE_Nut_out = LOONE_Nut(loone_q_path, ensemble_number, data_dir)
    LOONE_Nut_Lds_out_df = pd.DataFrame(LOONE_Nut_out)
    LOONE_Nut_Lds_out_df.to_csv(nut_out_path)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "loone_q_path",
        nargs=1,
        help="Path to LOONE Q CSV file.",
    )
    argparser.add_argument(
        "nut_out_path",
        nargs=1,
        help="Path to LOONE nutrient file to be created.",
    )
    argparser.add_argument(
        "data_dir",
        nargs=1,
        help="Path to data directory.",
    )
    argparser.add_argument(
        "ensemble_number",
        nargs=1,
        help="The ensemble number of the ensemble that the input data is for. Ex: '1'.",
    )
    args = argparser.parse_args()
    loone_q_path = args.loone_q_path[0]
    nut_out_path = args.nut_out_path[0]
    data_dir = args.data_dir[0]
    ensemble_number = args.ensemble_number[0]

    main(loone_q_path, nut_out_path, data_dir, ensemble_number)
