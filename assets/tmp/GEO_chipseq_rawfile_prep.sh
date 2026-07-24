#first average files 

#"/mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1"
"/media/lucio/easystore/Lucio/Analysis/Lara/Lara_ChIP_with_spikein_A1/"


"Double_KO_K27me3_ Double_KO_K4me1_ Double_KO_K4me2_ Double_KO_K4me3_ FC_FC_K27me3_ FC_FC_K4me1_ FC_FC_K4me2_ FC_FC_K4me3_ F_F_K27me3_ F_F_K4me1_ F_F_K4me2_ F_F_K4me3_ KO_D4_K27ac KO_D4_K27me3 KO_D4_K4me2 KO_D4_K4me3 Mll1-KO_K27me3 Mll1-KO_K4me1 Mll1-KO_K4me2 Mll1-KO_K4me3"



# nospikein

/mnt/datawk1/analysis/Lara/Lara_CHiPseq_nospikein_A2
"Anti-GFP Anti-Menin  Anti-Mll1_ Double_KO_RbBP5_ Mll2_KO_Mll1_ Anti-RbBp5_"

#rerun
/media/lucio/easystore/Lucio/Analysis/Lara/2024_04_Lara_chip/ssin_2024
_04_Lara_chip.csv

#replicate C replace repplicate A

"WT_H3K27ac_C Mll2KO_H3K27ac_C Mll1KO_H3K27ac_C DoubleKO_H3K27ac_C WT_H3K9me3_A
DoubleKO_H3K9me3_A"


#d4

#pattern="D4_DKO_H3K27ac_ D4_DKO_H3K27me3_ D4_WT_H3K27ac_ D4_WT_H3K27me3_ D4_WT_MLL2_ D4_WT_RbBP5_"
#pattern2="D4_WT_H3K4me2_ D4_WT_H3K4me3_ D4_DKO_H3K4me2_ D4_DKO_H3K4me3_"

/mnt/datawk1/analysis/Lara/2024_07_Lara_chip/

/mnt/datawk1/analysis/Lara/240926_chip_D4/


rsync -av volume1/Data_Just_In_NAS/Lucio/Data/Lara/CHiPseq_MLL_spikein/fastq/*


rsync --progress --append --inplace -av lucio@193.144.215.241:volume1/Data_Just_In_NAS/Lucio/Data/Lara/CHiPseq_MLL_spikein/fastq/* /media/lucio/easystore/bck_data/Lara/ChIPseq/CHiPseq_MLL_spikein/

# rsync --progress --append --inplace -av \
#   -e "ssh -p 22" \
#   "lucio@193.144.215.241:/volume1/Data_Just_In_NAS/Lucio/Data/Lara/CHiPseq_MLL_spikein/fastq/" \
#   "/media/lucio/easystore/bck_data/Lara/ChIPseq/CHiPseq_MLL_spikein/"

  scp lucio@193.144.215.241:/volume1/Data_Just_In_NAS/Lucio/Data/Lara/CHiPseq_MLL_spikein/fastq/* /media/lucio/easystore/bck_data/Lara/ChIPseq/CHiPseq_MLL_spikein/