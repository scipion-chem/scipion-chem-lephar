# coding: latin-1
# **************************************************************************
# *
# * Authors:     roberto marabini
# *
# * [1] uam, madrid, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
@article{C6CP01555G,
author ="Wang, Zhe and Sun, Huiyong and Yao, Xiaojun and Li, Dan and Xu, Lei and Li, Youyong and Tian, Sheng and Hou, Tingjun",
title  ="Comprehensive evaluation of ten docking programs on a diverse set of protein?ligand complexes: the prediction accuracy of sampling power and scoring power",
journal  ="Phys. Chem. Chem. Phys.",
year  ="2016",
volume  ="18",
issue  ="18",
pages  ="12964-12975",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/C6CP01555G",
url  ="http://dx.doi.org/10.1039/C6CP01555G",
abstract  ="As one of the most popular computational approaches in modern structure-based drug design{,} molecular docking can be used not only to identify the correct conformation of a ligand within the target binding pocket but also to estimate the strength of the interaction between a target and a ligand. Nowadays{,} as a variety of docking programs are available for the scientific community{,} a comprehensive understanding of the advantages and limitations of each docking program is fundamentally important to conduct more reasonable docking studies and docking-based virtual screening. In the present study{,} based on an extensive dataset of 2002 protein?ligand complexes from the PDBbind database (version 2014){,} the performance of ten docking programs{,} including five commercial programs (LigandFit{,} Glide{,} GOLD{,} MOE Dock{,} and Surflex-Dock) and five academic programs (AutoDock{,} AutoDock Vina{,} LeDock{,} rDock{,} and UCSF DOCK){,} was systematically evaluated by examining the accuracies of binding pose prediction (sampling power) and binding affinity estimation (scoring power). Our results showed that GOLD and LeDock had the best sampling power (GOLD: 59.8% accuracy for the top scored poses; LeDock: 80.8% accuracy for the best poses) and AutoDock Vina had the best scoring power (rp/rs of 0.564/0.580 and 0.569/0.584 for the top scored poses and best poses){,} suggesting that the commercial programs did not show the expected better performance than the academic ones. Overall{,} the ligand binding poses could be identified in most cases by the evaluated docking programs but the ranks of the binding affinities for the entire dataset could not be well predicted by most docking programs. However{,} for some types of protein families{,} relatively high linear correlations between docking scores and experimental binding affinities could be achieved. To our knowledge{,} this study has been the most extensive evaluation of popular molecular docking programs in the last five years. It is expected that our work can offer useful information for the successful application of these docking tools to different requirements and targets."
}

"""