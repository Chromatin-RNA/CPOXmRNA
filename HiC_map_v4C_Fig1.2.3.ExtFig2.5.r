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
### FIGURE 2A
###==========================================================================================

lateG1_hic_5kb = load_contacts(signal_path = '/data/xieb2/G1ER/4DNFI6H926RO.hic',
                                sample_name = "G1ER_lateG1",
                                resolution = 5000,
                                balancing = 'T', # this is the default
                                colour = "red")


lateG1_hic_5kb_ST3GAL6_v4c = virtual_4C(lateG1_hic_5kb, 
                                     data.frame("16",58508421,58529955), 
                                     xlim = c(200000,200000))


bed_CPOX = data.frame("16", 58670292, 58680386)
bed_DCBLD2 = data.frame("16", 58408443,58469727)
bed_ST3GAL6 = data.frame("16",58468125, 58523426 )
bed_CPOXeRNA = data.frame("16", 58657447, 58662569)
bed_E330017A01Rik = data.frame("16", 58635167, 58638739)

								 
bed_list_DCBLD2P_3 = list(bed_CPOX,bed_DCBLD2,bed_ST3GAL6, bed_CPOXeRNA, bed_E330017A01Rik)
									 
visualise(lateG1_hic_5kb_ST3GAL6_v4c, bedlist = bed_list_DCBLD2P_3,
          bed_colours = c("dodgerblue","red","green","purple","orange"))+ theme(axis.text = element_text(size = 20))
		  
		  
###==========================================================================================
### FIGURE 3B
###==========================================================================================


lateG1_hic_5kb_v4c = virtual_4C(lateG1_hic_5kb, 
                                     data.frame("16",58406502, 58409393), 
                                     xlim = c(50000,350000))
visualise(lateG1_hic_5kb_v4c, bedlist = bed_list_DCBLD2P_3,
          bed_colours = c("dodgerblue","red","green","purple","orange"))+ theme(axis.text = element_text(size = 20))
		  
		  
###==========================================================================================
### FIGURE 3C
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



###==========================================================================================
### Extended Data Fig. 2A
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




###==========================================================================================
### Extended Data Fig. 5B 
###==========================================================================================	

bed_KIF2A = data.frame("13", 106958996, 107022118)
bed_IPO11 = data.frame("13", 106794439, 106936958)
bed_Dimt1 = data.frame("13", 106947159, 106959760)

bed_list_KIF2A = list(bed_KIF2A, bed_IPO11, bed_Dimt1)


LateG1_5kb_KIF2A_IPO11P_2 = virtual_4C(lateG1_hic_5kb, 
                                     data.frame("13", 106919481, 106940797), 
                                     xlim = c(250000,300000))
visualise(LateG1_5kb_KIF2A_IPO11P_2, bedlist = bed_list_KIF2A, 
          bed_colours = c("dodgerblue","red","green","orange"))
		  

###==========================================================================================
### Extended Data Fig. 5E
###==========================================================================================	

bed_CAR2 = data.frame("3", 14886428,14900770)
bed_LRRCC1 = data.frame("3", 14533788, 14572658)


bedlist_LRRCC1_P = list(bed_CAR2, bed_LRRCC1)
LateG1_5kb_LRRCC1P_v4c = virtual_4C(lateG1_hic_5kb, "3:14529451-14535264", xlim = c(100000,400000))
visualise(LateG1_5kb_LRRCC1P_v4c, bedlist = bedlist_LRRCC1_P, 
          bed_colours = c("dodgerblue", "limegreen","red"))

