import numpy as np
import pandas as pd
import os
import random
from ctgan import CTGAN
from sklearn.model_selection import train_test_split


if __name__=="__main__":
    mirna_path = 'path/to/mirna/data' # EDIT FIELD!!
    snp_path = 'path/to/snp/data' # EDIT FIELD!!

    # Data must have samples on rows and miRNA features on columns with a label column
    mirna_data = pd.read_csv(mirna_path)
    # Specify the column name of label data
    mirna_label = ['class']

    # Set epochs to 10
    mirna_ctgan = CTGAN(epochs=10)
    mirna_ctgan.fit(mirna_data, mirna_label)

    # Generate 200 samples with trained model
    synthetic_data_mirna = mirna_ctgan.sample(200)

    # Data must have samples on rows and SNP features on columns with a label column
    snp_data = pd.read_csv(snp_path)
    # Specify the column name of label data
    snp_label = ['class']

    # Set epochs to 10
    snp_ctgan = CTGAN(epochs=10)
    snp_ctgan.fit(snp_data, snp_label)

    # Generate 200 samples with trained model
    synthetic_data_snp = snp_ctgan.sample(200)

    # Saving generated data to csv files for postprocessing
    gen_mirna_path = "synthetic_data_mirna.csv"
    # np.savetxt(gen_mirna_path, synthetic_data_mirna.to_numpy(), delimiter=",")
    gen_snp_path = "synthetic_data_snp.csv"
    # np.savetxt(gen_snp_path, synthetic_data_snp.to_numpy(), delimiter=",")

    syntetic_data_mirna = pd.read_csv(gen_mirna_path, sep=",", header=0, index_col=None)
    syntetic_data_snp = pd.read_csv(gen_snp_path, sep=",", header=0, index_col=None)

    # Select randomly generated miRNA samples
    dgf_index = np.random.choice(syntetic_data_mirna["class"].value_counts()[0], 42)
    egf_sampes = np.arange(start=syntetic_data_mirna["class"].value_counts()[0],
                           stop=syntetic_data_mirna["class"].shape[0])
    egf_index = np.random.choice(egf_sampes, 43)

    res = [*dgf_index, *egf_index]
    syntetic_data_mirna = syntetic_data_mirna.iloc[res]
    syntetic_data_mirna["data"] = ["syn"] * syntetic_data_mirna.shape[0]
    # Do the same for SNP data
    dgf_index = np.random.choice(syntetic_data_snp["class"].value_counts()[0], 42)
    egf_sampes = np.arange(start=syntetic_data_snp["class"].value_counts()[0],
                           stop=syntetic_data_snp["class"].shape[0])
    egf_index = np.random.choice(egf_sampes, 43)

    res = [*dgf_index, *egf_index]
    syntetic_data_snp = syntetic_data_snp.iloc[res]
    syntetic_data_snp["data"] = ["syn"] * syntetic_data_snp.shape[0]
    # Merge real and synthetic data
    mirna_data["data"] = ["real"] * mirna_data.shape[0]
    mirna = mirna_data.append(syntetic_data_mirna, ignore_index=True)
    snp_data["data"] = ["real"] * snp_data.shape[0]
    snp = snp_data.append(syntetic_data_snp, ignore_index=True)

    idx = np.arange(mirna.shape[0])
    mirna.drop(columns=["class"], inplace=True)
    # Split into training and test datasets
    _, _, idx_train, idx_test = train_test_split(mirna.to_numpy(), idx,
                                                 test_size=0.12, random_state=42,
                                                 stratify=mirna_data["data"])

    mirna.drop(columns=["data"], inplace=True)
    snp.drop(columns=["data"], inplace=True)
    mirna_train = mirna.to_numpy()[idx_train, :]
    mirna_test = mirna.to_numpy()[idx_test, :]
    snp_train = snp.to_numpy()[idx_train, :]
    snp_test = snp.to_numpy()[idx_test, :]

    _, labels_tr = np.unique(np.array(mirna["class"])[idx_train], return_inverse=True)
    # names_tr = np.array(sample_names_mirna)[mirna_idx_train]
    _, labels_te = np.unique(np.array(snp["class"])[idx_test], return_inverse=True)

    # Save csv data for MOGONET training
    feat_1 = open('1_featname.csv', 'w', newline='')
    for n in mirna.columns[:-1]:
        feat_1.write(n + "\n")

    feat_2 = open('2_featname.csv', 'w', newline='')
    for n in snp.columns[:-1]:
        feat_2.write(n + "\n")

    np.savetxt("1_tr.csv", mirna_train, delimiter=",")
    np.savetxt("1_te.csv", mirna_test, delimiter=",")
    np.savetxt("2_tr.csv", snp_train, delimiter=",")
    np.savetxt("2_te.csv", snp_test, delimiter=",")
    np.savetxt("labels_tr.csv", labels_tr, delimiter=",")
    np.savetxt("labels_te.csv", labels_te, delimiter=",")
