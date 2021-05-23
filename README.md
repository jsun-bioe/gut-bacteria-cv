# gut-bacteria-cv
Phenotypic characterization of human gut bacteria using computer vision-guided high-throughput imaging

## Introduction
Gut bacteria is a complex community with >300 species and they have profound impacts on our health. High-throughput imaging provides a rich source of information to understand the functional and physiological capabilities of each individual species. Here we present an automated computer vision-based pipeline that can segment single cells from phase-contrast images, track cell lineages and extract growth parameters in a high-throughput manner.

## Data
Time-lapse microscopy data is acquired for >90 representative isolate strains among the gut commensal bacteria using a high-throughput imaging technique called SLIP (Strain Library Imaging Protocol), recently developed in the Huang lab. The images are phase-contrast. All strains are acquired from the Human Microbiome Project.

## Pipeline
The pipeline includes algorithms for image segmentation, cell contour extraction and cell tracking. It uses a deep learning framework called DeepCell for single-cell image segmentation, which can be deployed on a computing server. Contour extraction, cell tracking and postprocessing is done in MATLAB with lightweight scripts or GUI.

### DeepCell
DeepCell is a deep learning-based image segmentation tool developed by the Covert lab in Stanford University. The version currently deployed in the Huang lab is written in Tensorflow. It takes in phase-contrast images as input and generates an output image of the same size with pixels classified as the cell edge highlighted. Several DeepCell models were trained on microscopy images of several morphologically distinct strains, and have been verified to achieve >90% accuracy for most strains tested.

### Morphometrics
Morphometrics is a MATLAB GUI-based image processing toolbox for bacterial single-cell analysis. Here we use the gradient detection algorithm to extract cell contours from the output images of DeepCell where the pixels in the contours have been highlighted. Morphometrics GUI allows quick parameter tuning. It also supports fast batch processing of high volume data, e.g. time-lapse images. The analysis can be performed locally on a PC, and the output is stored as .mat files for downstream analysis.

### Lineage Tracking
In order to capture the growth dynamics of individual strains, a tracking algorithm based on cell contour information is developed to reconstruct cell lineages from tens or hundreds of frames of data in a time-lapse. The code is written in MATLAB to simplify data transder from Morphometrics. For bacteria, growth and division are the key routines in their proliferation. The tracking algorithm can detect division events, and correct for artifacts such as cell missing or merging which may arise from noise in imaging or previous processing steps. The main output of this algorithm is a cell lineage in tree structure.

### Post-processing
As one application of high-throughput library imaging and the processing pipeline, we aim to extract growth parameters of individual strains from time-lapse images, including cell size, growth rate, cell cycle duration, etc., which are important single-cell quantities or properties underpinning the basic growth behavior and shape control mechanisms of the diverse species in the gut microbiota. For instance, the birth size, division size and inter-division duration can together indicate the mechanism by which the cell maintain a relatively constant size between cycles.

As a short summary, when combined with computer vision and many other algorithms, high-throughput imaging is becoming a powerful tool for library scale phenotyping, and is expected to offer insights into the human gut microbiota and other complex bacterial communities.
