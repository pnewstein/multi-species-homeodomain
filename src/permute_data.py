from pathlib import Path

import pandas as pd


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
    df = pd.read_csv(HERE / "../data/cartridge_count.csv")
    real = get_difference(df)
    tries = 10_000
    for i in range_trices()

