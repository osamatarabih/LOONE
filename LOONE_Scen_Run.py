import argparse
import pandas as pd
from LOONE_Nut import LOONE_Nut


def main(loone_q_path: str, nut_out_path: str, data_dir: str, loads_external_filename: str, flow_df_filename: str) -> None:
    """Run LOONE_Nut_Lds_out.csv

    Args:
        loone_q_path (str): Path to LOONE Q CSV file.
        nut_out_path (str): Path to LOONE nutrient file to be created.
        data_dir (str): Path to data directory.
        loads_external_filename (str): The name of the file that holds the external loads for the model.
        flow_df_filename (str): The name of the file that holds the flow data for the model.
    """
    LOONE_Nut_out = LOONE_Nut(loone_q_path, loads_external_filename, flow_df_filename, data_dir)
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
        "loads_external_filename",
        nargs=1,
        help="The name of the file that holds the external loads for the model.",
    )
    argparser.add_argument(
        "flow_df_filename",
        nargs=1,
        help="The name of the file that holds the flow data for the model.",
    )
    args = argparser.parse_args()
    loone_q_path = args.loone_q_path[0]
    nut_out_path = args.nut_out_path[0]
    data_dir = args.data_dir[0]
    loads_external_filename = args.loads_external_filename[0]
    flow_df_filename = args.flow_df_filename[0]

    main(loone_q_path, nut_out_path, data_dir, loads_external_filename, flow_df_filename)
