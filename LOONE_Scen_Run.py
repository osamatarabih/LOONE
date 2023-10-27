import argparse
import pandas as pd
from LOONE_Nut import LOONE_Nut


def main(loone_q_path: str, nut_out_path: str, data_dir: str) -> None:
    """Run LOONE_Nut_Lds_out.csv

    Args:
        loone_q_path (str): Path to LOONE Q CSV file.
        nut_out_path (str): Path to LOONE nutrient file to be created.
        data_dir (str): Path to data directory.
    """

    LOONE_Nut_out = LOONE_Nut(loone_q_path)
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
    args = argparser.parse_args()
    loone_q_path = args.loone_q_path[0]
    nut_out_path = args.nut_out_path[0]
    data_dir = args.data_dir[0]

    main(loone_q_path, nut_out_path, data_dir)
