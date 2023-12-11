###################################################################################################
###                                                                                             ###
###               THIS FILE CONTAINS CODE FOR FIG.1,2 & SUPPLEMENTARY FIG.1, 3                  ###
###                                                                                             ###
###################################################################################################



###===========================================
### INSTALL GENOVA PACKAGE DEV VERSION
###===========================================

remotes::install_github("robinweide/GENOVA", ref = "dev")

###===========================================

library(GENOVA)
library(dplyr)
library(ggplot2)


###==========================================================================================
### FIGURE 1A
###==========================================================================================


FLC_hic_25kb = load_contacts(signal_path = '/data/xieb2/juicer16_FLC_mm10/aligned/inter_30.hic',
                                sample_name = "FLC",
                                 resolution = 25000,
                                 balancing = 'T', # this is the default
                                 colour = "cornflowerblue")


pyramid(exp = FLC_hic_25kb,
chrom = 'chr16',
start = 58200000,
end=58900000,
colour = scale_fill_gradient(low = "white", high = "red",limits = c(0,120)),
crop_y = c(0, 0.45e6))




		  
		  
###==========================================================================================
### FIGURE 2D
###==========================================================================================


lateG1_hic_5kb_v4c = virtual_4C(lateG1_hic_5kb, 
                                     data.frame("16",58406502, 58409393), 
                                     xlim = c(50000,350000))
visualise(lateG1_hic_5kb_v4c, bedlist = bed_list_DCBLD2P_3,
          bed_colours = c("dodgerblue","red","green","purple","orange"))+ theme(axis.text = element_text(size = 20))
		  
		  
###==========================================================================================
### FIGURE 2E
###==========================================================================================		  
		  

bed_CPOX_G1 = data.frame("chr16", 58670292, 58680386)
bed_DCBLD2_G1 = data.frame("chr16", 58408443,58469727)
bed_ST3GAL6_G1 = data.frame("chr16",58468125, 58523426 )
bed_CPOXeRNA_G1 = data.frame("chr16", 58657447, 58662569)
bed_E330017A01Rik_G1 = data.frame("chr16", 58635167, 58638739)


bed_list_DCBLD2P_3_G1= list(bed_CPOX_G1,bed_DCBLD2_G1,bed_ST3GAL6_G1, bed_CPOXeRNA_G1, bed_E330017A01Rik_G1)

FLC_hic_5kb = load_contacts(signal_path = '/data/xieb2/juicer16_FLC_mm10/aligned/inter_30.hic',
                                sample_name = "FLC",
                                 resolution = 5000,
                                 balancing = 'T', # this is the default
                                 colour = "cornflowerblue")

FLC_hic_5kb_v4c = virtual_4C(FLC_hic_5kb, 
                                     data.frame("chr16",58406502, 58409393), 
                                     xlim = c(50000,350000))
visualise(FLC_hic_5kb_v4c, bedlist = bed_list_DCBLD2P_3_G1,
          bed_colours = c("dodgerblue","red","green","purple","orange"))+ theme(axis.text = element_text(size = 20))

###============================================================================================
### Supplementary Fig. 1G
###============================================================================================
Rd_0h_5K = load_contacts(signal_path = 'GSM5602658_red_blood_0h-Arima.mcool.multires.cool',
                                sample_name = "red_blood_0h",
                                resolution = 5000,
                                balancing = 'T',
                                colour = "cornflowerblue")


Rd_DMSO_5K = load_contacts(signal_path = 'GSM5602660_red_blood_DMSO-Arima.mcool.multires.cool',
                                sample_name = "red_blood_DMSO",
                                resolution = 5000,
                                balancing = 'T', 
                                colour = "green")



pdf(file = "diff_hic_map.5k.pdf")
hic_matrixplot(exp1 = Rd_0h_5K,
exp2 = Rd_DMSO_5K,
chrom = 'chr16',
start = 58370000,
end=58740000,
cut.off = 120) 
dev.off()

###==========================================================================================
### Supplementary Fig. 3B
###==========================================================================================		

### CH12.LX TAD
CH12LX_hic_10kb = load_contacts(signal_path = '/data/xieb2/G1ER/OSN/CH12.LX.4DNFI8KBXYNL.hic',
                                  sample_name = "CH12.LX",
                                  resolution = 10000,
                                  balancing = 'T', # this is the default
                                  colour = "green4")
pyramid(exp = CH12LX_hic_10kb,
chrom = '16',
start = 58390000,
end=58730000,
colour = scale_fill_gradient(low = "white", high = "firebrick2",limits = c(0,50)))								  


### G1ER late G1 TAD								  
lateG1_hic_10kb = load_contacts(signal_path = '/data/xieb2/G1ER/4DNFI6H926RO.hic',
                                sample_name = "G1ER_lateG1",
                                resolution = 10000,
                                balancing = 'T', # this is the default
                                colour = "red")
								


pyramid(exp = lateG1_hic_10kb,
chrom = '16',
start = 58390000,
end=58730000,
colour = scale_fill_gradient(low = "white", high = "firebrick2",limits = c(0,50)))

### differential HiC map G1ER vs. CH12.LX
pyramid_difference(
lateG1_hic_10kb,
CH12LX_hic_10kb,
chrom = "16", start = 58390000, end = 58730000
)




