#----------------------------------------------------------------
#- NOTE : I need to do this in an interactive shell, otherwise
#-      some of the bigwigs are too big and the bigWigAverageOverBed
#-      step gets killed
#-
#-      qrsh -l tmem=64G,h_vmem=64G
#----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(stringr)



for (file in c("average_signals/human_nopseudo_average_signal.txt",
		"average_signals/mouse_nopseudo_average_signal.txt",
		"average_signals/human_average_signal.txt",
		"average_signals/mouse_average_signal.txt")){
	
	means <- read.table(file, sep="\t", header=T)

	if (str_detect(file, "human")) {
		species = "hs"
	}else{
		species = "mm"
	}
	if (str_detect(file, "nopseudo")) {
		type = "_nopseudo"
	}else{
		type = ""
	}

	dir.create(paste0("scatterplots/"), showWarnings = FALSE)
	dir.create(paste0("scatterplots_GFP/"), showWarnings = FALSE)

	samples <- unique(paste0(means$Sample, "_", means$Method))
	samples2 <- unique(paste0(means$Sample, "_", means$Method, "-", means$Tag))

	#- sample_plus versus sample_minus
	#- Do this for any sample that has both plus and minus data
	for (sample in samples){
		if ((paste0(sample, "-", "PlusTag") %in% samples2) && (paste0(sample, "-", "MinusTag") %in% samples2)){
			
			fileplus <- read.table(paste0("bedgraphs",type,"_",species,"/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".100bp.genic.bg"), sep="\t", header=F)
			fileminus <- read.table(paste0("bedgraphs",type,"_",species,"/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".100bp.genic.bg"), sep="\t", header=F)

			fileplus <- unique(fileplus)
			fileminus <- unique(fileminus)

			colnames(fileplus) <- c("chr","start","end","name","signal_p")
			colnames(fileminus) <- c("chr","start","end","name","signal_m")

			fileall <- full_join(fileplus, fileminus, by=c("chr","start","end", "name"))
			fileall <- fileall %>% filter(!is.na(signal_p) & !is.na(signal_m))
			
			#- make a scatterplot 
			png(paste0("scatterplots/scatterplot_", sample, type,"_",species,".png"))
			print(paste0("scatterplots/scatterplot_", sample, type,"_",species,".png"))
			range_vals <- range(c(fileall$signal_m , fileall$signal_p))
			p <- ggplot(fileall, aes(x=signal_m, y=signal_p)) +
				geom_point() +
				geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  # Add diagonal line
				ggtitle(sample) +
				xlab("MinusTag") +
				ylab("PlusTag") +
				theme_minimal() + 
				xlim (range_vals) +
				ylim (range_vals) + 
				#coord_fixed(ratio = 1)+
   				theme(
    			    text = element_text(size = 18),        # Default text size
    			    axis.title = element_text(size = 20),  # Axis labels
    			    axis.text = element_text(size = 16),   # Axis tick labels
    			    plot.title = element_text(size = 22, face = "bold"))  # Title size and weight	
			print(p)
			dev.off()
		}

		#- sample_plus versus GFP_plus
		#– Do this for any plus sample that is not GFP
		if ((str_split_fixed(sample, "_", n=2)[1] != "GFP") && 
			(paste0(sample, "-", "PlusTag") %in% samples2) &&
			(paste0("GFP_", str_split_fixed(sample, "_", n=2)[2], "-PlusTag") %in% samples2)){
			
			fileplus <- read.table(paste0("bedgraphs",type,"_",species,"/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".100bp.genic.bg"), sep="\t", header=F)
			gfpplus <- read.table(paste0("bedgraphs",type,"_",species,"/GFP_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".100bp.genic.bg"), sep="\t", header=F)

			fileplus <- unique(fileplus)
			gfpplus <- unique(gfpplus)

			colnames(fileplus) <- c("chr","start","end","name","signal_p")
			colnames(gfpplus) <- c("chr","start","end","name","signal_m")

			fileall <- full_join(fileplus, gfpplus, by=c("chr","start","end", "name"))
			fileall <- fileall %>% filter(!is.na(signal_p) & !is.na(signal_m))
			
			#- make a scatterplot
			png(paste0("scatterplots_GFP/scatterplot_", sample,type, "_",species,".png"))
			range_vals <- range(c(fileall$signal_m, fileall$signal_p))
			p <- ggplot(fileall, aes(x=signal_m, y=signal_p)) +
				geom_point() +
				geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  # Add diagonal line
				ggtitle(sample) +
				xlab("GFP PlusTag") +
				ylab(paste0(str_split_fixed(sample, "_", n=2)[1], " PlusTag")) +
				theme_minimal()  + 
				xlim (range_vals) +
				ylim (range_vals) + 
				#coord_fixed(ratio = 1)+
   				theme(
    			    		text = element_text(size = 18),        # Default text size
    			    		axis.title = element_text(size = 20),  # Axis labels
    			    		axis.text = element_text(size = 16),   # Axis tick labels
    			    		plot.title = element_text(size = 22, face = "bold"))  # Title size and weight	
			print(p)
			dev.off()
		}
	}
}


