# -*- coding: utf-8 -*-

"""
@article{CHAILLET2025100125,
title = {pytom-match-pick: A tophat-transform constraint for automated classification in template matching},
journal = {Journal of Structural Biology: X},
volume = {11},
pages = {100125},
year = {2025},
issn = {2590-1524},
doi = {http://doi.org/10.1016/j.yjsbx.2025.100125},
url = {http://www.sciencedirect.com/science/article/pii/S2590152425000066},
author = {Marten L. Chaillet and Sander Roet and Remco C. Veltkamp and Friedrich Förster},
keywords = {cryo-ET, Tomograms, Particle localization, Identification, Template matching, Volume registration, GPU-acceleration},
abstract = {Template matching (TM) in cryo-electron tomography (cryo-ET) enables in situ detection and localization of known macromolecules. However, TM faces challenges of weak signal of the macromolecules and interfering features with a high signal-to-noise ratio, which are often addressed by time-consuming, subjective manual curation of results. To improve the detection performance we introduce pytom-match-pick, a GPU-accelerated, open-source command line interface for enhanced TM in cryo-ET. Using pytom-match-pick, we first quantify the effects of point spread function (PSF) weighting and show that a tilt-weighted PSF outperforms a binary wedge with a single defocus estimate. We also assess previously introduced background normalization methods for classification performance. This indicates that phase randomization is more effective than spectrum whitening in reducing false positives. Furthermore, a novel application of the tophat transform on score maps, combined with a dual-constraint thresholding strategy, reduces false positives and improves precision. We benchmarked pytom-match-pick on public datasets, demonstrating improved classification and localization of macromolecules like ribosomal subunits and proteasomes that led to fewer artifacts in subtomogram averages. This tool promises to advance visual proteomics by improving the efficiency and accuracy of macromolecule detection in cellular contexts.}
}


@Article{ijms241713375,
AUTHOR = {Chaillet, Marten L. and van der Schot, Gijs and Gubins, Ilja and Roet, Sander and Veltkamp, Remco C. and Förster, Friedrich},
TITLE = {Extensive Angular Sampling Enables the Sensitive Localization of Macromolecules in Electron Tomograms},
JOURNAL = {International Journal of Molecular Sciences},
VOLUME = {24},
YEAR = {2023},
NUMBER = {17},
ARTICLE-NUMBER = {13375},
URL = {http://www.mdpi.com/1422-0067/24/17/13375},
PubMedID = {37686180},
ISSN = {1422-0067},
ABSTRACT = {Cryo-electron tomography provides 3D images of macromolecules in their cellular context. To detect macromolecules in tomograms, template matching (TM) is often used, which uses 3D models that are often reliable for substantial parts of the macromolecules. However, the extent of rotational searches in particle detection has not been investigated due to computational limitations. Here, we provide a GPU implementation of TM as part of the PyTOM software package, which drastically speeds up the orientational search and allows for sampling beyond the Crowther criterion within a feasible timeframe. We quantify the improvements in sensitivity and false-discovery rate for the examples of ribosome identification and detection. Sampling at the Crowther criterion, which was effectively impossible with CPU implementations due to the extensive computation times, allows for automated extraction with high sensitivity. Consequently, we also show that an extensive angular sample renders 3D TM sensitive to the local alignment of tilt series and damage induced by focused ion beam milling. With this new release of PyTOM, we focused on integration with other software packages that support more refined subtomogram-averaging workflows. The automated classification of ribosomes by TM with appropriate angular sampling on locally corrected tomograms has a sufficiently low false-discovery rate, allowing for it to be directly used for high-resolution averaging and adequate sensitivity to reveal polysome organization.},
DOI = {10.3390/ijms241713375}
}



"""