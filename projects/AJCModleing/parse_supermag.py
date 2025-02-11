import pandas as pd

def parse_csv(fname="datasets/20250211-20-39-supermag.csv"):
    o = pd.read_csv(fname, parse_dates=["Date_UTC"])
    iagas = o.IAGA.unique()
    for iaga in iagas:
        x = o[o.IAGA==iaga][["Date_UTC","dbn_geo","dbe_geo","dbz_geo"]]
        x.rename(
            columns={
                "Date_UTC": "Date",
                "dbn_geo": "Y",
                "dbe_geo": "X",
                "dbz_geo": "Z",
            }, 
            inplace=True
        )
        x.to_csv(f"datasets/{iaga}.csv", header=True, index=False)
    return

parse_csv()