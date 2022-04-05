#! /usr/bin/Rscript
## Collect arguments
args <- commandArgs(TRUE)
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
if (length(parseArgs(args)) > 1) {
  argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
  names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
  args <- argsL
  rm(argsL)
  if ((args$data_type %in% c("RNAseq","microarray","scRNA"))==FALSE){
    args <- "--help"
  }
  } else {
  args <- "--help"
}

#print(args)

if (is.null(unlist(args))){
  args <- "--help"
}

# ## Give some value to options if not provided 
# if(is.null(args$opt_arg1)) {args$opt_arg1="default_option1"}
# if(is.null(args$opt_arg2)) {args$opt_arg2="default_option1"} else {args$opt_arg2=as.numeric(args$opt_arg2)}
## Default setting when no all arguments passed or help needed
if("--help" %in% args ) {
  cat("
    Xiangyu Pan
    11/11/2021
      The R Script of CCI-Identify
      
      Mandatory arguments:
      --BS_BK_DK_sig=BS_BK_DK_sig            - the pathway of BS, BK, and DK signatrues
      --input=transcriptome_data             - the input of transcriptome matrix
      --species=species                      - the species of samples (mouse or human)
      --data_type=type_of_input              - the types of input data (RNAseq, microarray and scRNA)
      --range=axis_range                     - the range of axis in ternary map
      --output=output_pathways               - the output pathway of CCI results
      --prefix=prefix                        - the prefix of output files
      --log_transform=log_transform          - the log-transformation before CCI quantification
      --help                                 - print this text

  WARNING : You should install the R packages ggplot2 (version=3.1), ggtern (version=3.1), dplyr, nichenetr and data.table.
            Only RNAseq, microarray and scRNA was supported for this function. And you should specialise your data_type.

  Example:
          Rscript ./CCI-Identify.R \
          --BS_BK_DK_sig=./BS_BK_DK_sig.csv \
          --input=./test_counts.csv \
          --species=mouse \
          --data_type=RNAseq \
          --range=0.45 \
          --output=./out_tmp/ \
          --prefix=TPM4oe_vs_Ctrl \
          --log_transform=YES \n\n")

  q(save="no")
}

message("loading packages")

suppressPackageStartupMessages({
  require(dplyr)
  require(nichenetr)
  require(data.table)
  require("ggplot2")
  require("ggtern")
  })

if (packageVersion("ggtern") !="3.1.0") {
  message("Although the CCI would be calculated, We still prefer you to install the version 3.1.0 of ggtern")
}
if (packageVersion("ggplot2") !="3.1.0") {
  message("To better match the version 3.1.0 of ggtern, we prefer you to install the version 3.1.0 of ggplot2")
}

Normal_SE_sig <- read.csv(args$BS_BK_DK_sig)
Normal_SE_sig$cluster <- as.character(Normal_SE_sig$cluster)
Normal_SE_sig_DK <- subset(Normal_SE_sig,cluster=="Differentiated_keratinocyte")
Normal_SE_sig_BK <- subset(Normal_SE_sig,cluster=="Basel_keratinocyte")
Normal_SE_sig_BS <- subset(Normal_SE_sig,cluster=="Basel_Stem_Cells")

if (args$data_type=="RNAseq") {
  factor <- 1
}
if (args$data_type=="microarray") {
  factor <- 1
}
if (args$data_type=="scRNA") {
  factor <- 10
}

input_c <- read.csv(args$input)
input_c$Gene <- as.character(input_c$Gene)

if (args$species!="human") {
  input_c = input_c %>% mutate(Gene = convert_mouse_to_human_symbols(Gene))
  input_c <- input_c[!duplicated(input_c$Gene),]
  input_c <- na.omit(input_c)
  input_c$Gene <- as.character(input_c$Gene)
}

rownames(input_c) <- input_c$Gene
input_c <- input_c[,-1]

if (args$log_transform=="YES") {
  input_c <- log(input_c+1,2)
}

CCI <- data.frame(DK=apply(input_c[intersect(rownames(input_c),Normal_SE_sig_DK$gene),],2,mean),
  BK=apply(input_c[intersect(rownames(input_c),Normal_SE_sig_BK$gene),],2,mean),
  BS=apply(input_c[intersect(rownames(input_c),Normal_SE_sig_BS$gene),],2,mean))
CCI$DK <- as.numeric(as.character(CCI$DK))
CCI$BK <- as.numeric(as.character(CCI$BK))
CCI$BS <- as.numeric(as.character(CCI$BS))

message("CCI calculation done")

a <- CCI$BK
b <- CCI$DK
c <- CCI$BS
CCI$cos1 <- (a^2 +(((a^2+b^2+c^2)^0.5)^2)-(((b^2+c^2)^0.5)^2))/(2*a*((a^2+b^2+c^2)^0.5))
CCI$cos2 <- (b^2 +(((a^2+b^2+c^2)^0.5)^2)-(((a^2+c^2)^0.5)^2))/(2*b*((a^2+b^2+c^2)^0.5))
CCI$cos3 <- (c^2 +(((a^2+b^2+c^2)^0.5)^2)-(((a^2+b^2)^0.5)^2))/(2*c*((a^2+b^2+c^2)^0.5))
CCI$cos1[which((2*a*((a^2+b^2+c^2)^0.5)==0))] <- 2
CCI$cos2[which((2*b*((a^2+b^2+c^2)^0.5)==0))] <- 2
CCI$cos3[which((2*c*((a^2+b^2+c^2)^0.5)==0))] <- 2
CCI$CCI_score <- apply(CCI[,c("cos1","cos2","cos3")],1,sd)
CCI$CCI_score <- factor/CCI$CCI_score

write.csv(CCI,paste0(args$output,"/",args$prefix,"_",args$species,"_",args$data_type,"_CCI.score.csv"))

message("CCI distribution done")

CCI$sample <- rownames(CCI)
p1 <- ggtern(CCI , aes(x = BK, y = BS, z = DK,color=sample))+
theme_rgbg()+
    geom_point(size = 5) + 
    guides(alpha = FALSE)+ tern_limits(T = as.numeric(args$range), L = as.numeric(args$range), R = as.numeric(args$range))
ggsave(paste0(args$output,"/",args$prefix,"_",args$species,"_",args$data_type,"_CCI.svg"), plot=p1,width = 5, height = 5,dpi=1080)
ggsave(paste0(args$output,"/",args$prefix,"_",args$species,"_",args$data_type,"_CCI.pdf"), plot=p1,width = 5, height = 5,dpi=1080)

aa <- c("#313695","#5083BB","#8FC3DD","#D2ECF4","#FFFFBF","#FDD384","#F88D51","#DE3F2E","#A50026")
p1 <- ggtern(CCI , aes(x = BK, y = BS, z = DK))+
  theme_rgbg()+
 geom_point(alpha = 0.6, color = "#aebfbb", size = 1) + 
    stat_density_tern(
      geom = 'polygon', 
      aes(fill = ..level.., alpha = ..level..),bins = 500,inherit.aes=TRUE) +
    guides(alpha = FALSE)+
    scale_fill_gradientn(colours=alpha(aa,1),na.value="NA")+ 
    scale_alpha(range = c(0, 0.8), guide = FALSE) + 
    tern_limits(T = as.numeric(args$range), L = as.numeric(args$range), R = as.numeric(args$range))
ggsave(paste0(args$output,"/",args$prefix,"_",args$species,"_",args$data_type,"_CCI_density.svg"), plot=p1,width = 5, height = 5,dpi=1080)
ggsave(paste0(args$output,"/",args$prefix,"_",args$species,"_",args$data_type,"_CCI_density.pdf"), plot=p1,width = 5, height = 5,dpi=1080)

message("The density of CCI distribution done")