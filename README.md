# PROMISED
PROMISED: Prediction of Renal Outcome using Multiomics and Integrated Statistical Evaluation of Delayed Graft Function in kidney transplant recipients.

<div align="left">
<br/>
<p align="center">
<img align="center" width=70% src="Graphical_abstract_DGF.svg"></img>
</p>
</div>

</div>

## Installation and Dependencies

This repository is tested under the following system settings:
- `Python 3.7.9` (is recommended to create a `venv`)

Clone this repository from Github

```cmd
git clone https://github.com/ginTom/PROMISED.git
```
#### CTGAN
- Install Python dependencies for CTGAN project

```cmd
pip install -r CTGAN_requirements.txt
```
#### MOGONET
Clone MOGONET repository from Github and install dependencies

```cmd
git clone https://github.com/txWang/MOGONET.git
pip install -r MOGONET_requirements.txt
```

## Prepare training data
Launch `CTGAN_main.py` to generate and obtain training data for each omic:
- `label_tr.csv`: labels for trainin set
- `label_te.csv`: labels for test set
- `{1,2}_featname.csv`: feature names for each omic
- `{1,2}_tr.csv`: traing data
- `{1,2}_te.csv`: test data

## Train MOGONET
Create folder with all training data;
Customize `data_folder` and `view_list` in `MOGONET/main_mogonet.py`

## Biomarkers identification
Customize `data_folder` and `view_list` in `MOGONET/main_biomarker.py`
