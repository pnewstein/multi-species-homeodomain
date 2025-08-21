"""
This file runs a permutation test on the data in cartridge_count
"""

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
    cond_sum = df.groupby("condition").sum()["count"]
    return cond_sum["ctr"] - cond_sum["exp"]


def main():
    """
    runs the permutation test 10,000 trials and prints the result
    """
    rng = np.random.default_rng(100)
    df = pd.read_csv(HERE / "../data/cartridge_count.csv")
    real = get_difference(df)
    out_list: list[None | float] = [None] * 10_000
    for i, _ in enumerate(out_list):
        df["condition"] = rng.permutation(np.array(df["condition"]))
        out_list[i] = get_difference(df)
    print(np.mean(np.array(out_list) < real))


if __name__ == "__main__":
    main()
