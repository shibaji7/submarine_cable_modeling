from glob import glob
from scipy.io import loadmat

def read_dataset(file):
    dat = loadmat(file)
    print(dat["data"][0,0])
    return dat

def load_datasets(folder="dataset/May2024/oceannetworksCA/*.mat"):
    files = glob(folder)
    for f in files:
        dat = read_dataset(f)
    return

if __name__ == "__main__":
    load_datasets()