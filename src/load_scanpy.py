from pathlib import Path

import pandas as pd
import numpy as np

try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()


def get_difference(df: pd.DataFrame) -> float:
    """
    gets the difference between the means of each condition
    """
    sum = df.groupby("condition").sum()["number"]
    return sum["ctr"] - sum["exp"]


def main():
    rng = np.random.default_rng(100)
    df = pd.read_csv(HERE / "../data/cartridge_count.csv")
    real = get_difference(df)
    out_list: list[None | float] = [None] * 10_000
    for i, _ in enumerate(out_list):
        df["condition"] = rng.permutation(np.array(df["condition"]))
        out_list[i] = get_difference(df)
    print(np.mean(out_list < real))
