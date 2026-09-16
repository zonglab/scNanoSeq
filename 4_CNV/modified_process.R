#!/usr/bin/env Rscript

##############################################################################################################################################
###### Main CNV analysis script
##############################################################################################################################################

####################################################################################################
### Variables
####################################################################################################

# Config
minPloidy   = 1.5
maxPloidy   = 6
minBinWidth = 5
main_path="/PATH/TO/YOUR/GINKGO/"
make_plots=TRUE

if (!file.exists(main_path)) {
  # Get main_path from full commandArgs
  arguments        = commandArgs(FALSE)
  scrptname = sub("--file=", "", grep("--file=", arguments, value = TRUE))  # script name
  main_path    = normalizePath(dirname(scrptname))
}

# User settings
arguments = commandArgs(TRUE)
PATH_genome =             arguments[[ 1]]   # dir to genome in the Ginkgo Home dir
PATH_2_work =             arguments[[ 2]]   # dir to put the results
FILE_runlog =             arguments[[ 3]]   # FILE_runlog file name for logging
FILE_bincnt =             arguments[[ 4]]   # data, read counts per bin
CFIG_SegRef =  as.numeric(arguments[[ 5]])  # segmentation reference, 0:independent; 1:global; 2:costom
METH_bining =             arguments[[ 6]]   # which binning resolution to use
METH_clustr =             arguments[[ 7]]   # clustering method
METH_distan =             arguments[[ 8]]   # distancing method
COL_PALETTE =  as.numeric(arguments[[ 9]])  # color palette
FILE_SegRef =             arguments[[10]]   # reference segmentation if CFIG_SegRef==2
CFIG_FacsPD =  as.numeric(arguments[[11]])  # facs ploidy: 0：no, 1:yes
PATH_FacsPD =             arguments[[12]]   # facs file name
CFIG_sexchr =  as.numeric(arguments[[13]])  # to include sex chrs: 0:no, 1:yes
CFIG_badbin =  as.numeric(arguments[[14]])  # remove bad bins


####################################################################################################
print("Initialize Variables & Pre-Process Data")
####################################################################################################
statusFile = file( paste(PATH_2_work, "/", FILE_runlog, sep="") )
writeLines(c("<?xml version='1.0'?>", "<FILE_runlog>", "<step>3</step>", 
             "<processingfile>Initializing Variables</processingfile>", 
             "<percentdone>0</percentdone>", "<tree>clust.xml</tree>", "</FILE_runlog>"), statusFile)
close(statusFile)

# Load genome specific files
setwd(PATH_genome)
MX_GC_cont = read.table(paste0("GC_",     METH_bining), header=FALSE, sep="\t", as.is=TRUE)
MX_bin_loc = read.table(                  METH_bining,  header=TRUE , sep="\t", as.is=TRUE)
MX_bin_cnt = read.table(paste0("bounds_", METH_bining), header=FALSE, sep="\t")

# Load user data
setwd(PATH_2_work)
MX_raw_CNP = read.table(FILE_bincnt, header=TRUE, sep="\t")
if (CFIG_FacsPD == 1) {
    print(paste("facs file inputed, reading ", PATH_FacsPD))
    ploidy = read.table(PATH_FacsPD, header=FALSE, sep="\t", as.is=TRUE)
} else {
    print("CNVs inferred without facs file.")
    ploidy = rbind(c(0,0), c(0,0))}

####################################################################################################
###### Remove bad bins
####################################################################################################
if (CFIG_badbin)
{   print("Removing bad bins...")
    badbins = read.table(paste(PATH_genome, "/badbins_", METH_bining, sep=""), header=FALSE, sep="\t", as.is=TRUE)
    MX_GC_cont      = data.frame(MX_GC_cont[-badbins[,1], 1])
    MX_bin_loc     = MX_bin_loc[-badbins[,1], ]
    MX_raw_CNP     = data.frame(MX_raw_CNP[-badbins[,1], ])

    step  = 1
    chrom = MX_bin_loc[1,1]
    for (i in 1:nrow(MX_bin_loc))
    {if (MX_bin_loc[i,1] != chrom)
     {      MX_bin_cnt[step,1] = chrom
         MX_bin_cnt[step,2] = i
         step           = step+1
         chrom          = MX_bin_loc[i,1]
}}}


####################################################################################################
###### Initialize color palette
####################################################################################################
colors     = matrix(0,3,2)
colors[1,] = c('goldenrod', 'darkmagenta')
colors[2,] = c('dodgerblue', 'darkorange')
colors[3,] = c('brown2', 'blue4')


####################################################################################################
###### Initialize data structures
####################################################################################################
l            = dim(MX_raw_CNP)[1] # Number of bins
w            = dim(MX_raw_CNP)[2] # Number of cells
breaks       = matrix(0,l,w)
fixed        = matrix(0,l,w)
fin_min      = matrix(0,l,w)
fin_cls      = matrix(0,l,w)
stats        = matrix(0,w,14)
pos          = cbind(c(1,MX_bin_cnt[,2]), c(MX_bin_cnt[,2], l))
# Initialize CN_fin_call inference variables
CNgrid       = seq(minPloidy, maxPloidy, by=0.05)
n_ploidy     = length(CNgrid)  # Number of ploidy tests during CN_fin_call inference
CNmult       = matrix(0,n_ploidy,w)
CNerror      = matrix(0,n_ploidy,w)
outerColsums = matrix(0,n_ploidy,w)

# Prepare statistics
lab     = colnames(MX_raw_CNP)
rownames(stats) = lab
colnames(stats) = c("Reads", "Bins", "Mean", "Var", "Disp", 
                    "Min", "25th", "Median", "75th", "Max", 
                    "CN__minimal", "MAD_minimal", "CN__closest", "MAD_closest")

####################################################################################################
print("Normalizing cells")
####################################################################################################
normal  = sweep(MX_raw_CNP+1, 2, colMeans(MX_raw_CNP+1), '/')
normal2 = normal



####################################################################################################
print("Determine segmentation ref") # using dispersion (CFIG_SegRef = 1) or reference sample (CFIG_SegRef = 2)
####################################################################################################
if (       CFIG_SegRef == 1) {
  F = normal[,which.min(apply(normal, 2, sd)/apply(normal,2,mean))[1]]
} else if (CFIG_SegRef == 2) {
  R   = read.table(FILE_SegRef, header=TRUE, sep="\t", as.is=TRUE)
  low = lowess(MX_GC_cont[,1], log(R[,1]+0.001), f=0.05)
  app = approx(low$x, low$y, MX_GC_cont[,1])
  F   = exp(log(R[,1]) - app$y) }


####################################################################################################
####################################################################################################
####################################################################################################
print("Estimating Ploidy")
####################################################################################################
####################################################################################################
####################################################################################################
# Open output stream
sink("results.txt")
cat(paste("Sample\tCopy_Number\tSoS_Predicted_Ploidy\tError_in_SoS_Approach\n", sep=""))

lowess.gc = function( in_x, in_y) {
    fit_curv = lowess(in_x, log(in_y), f=0.05); 
    est__out = approx(fit_curv$x, fit_curv$y, in_x)
    return(exp(log(in_y) - est__out$y))}

if( make_plots == TRUE){
    dir.create("1_covdist/")
    dir.create("2_bin_cnt/")
    dir.create("3_lorenz/")
    dir.create("4_GC_lowe/")
    dir.create("5_CN_hist/")
    dir.create("6_SoS/")
    dir.create("7_CN_Profile/")}

####################################################################################################
print("Process each cell individually")
####################################################################################################
suppressWarnings({
for(k in 1:w)
{   message(paste0("working on cell ", k, ": ", lab[k]))
    # Generate basic statistics
    stats[k,1]  = sum(MX_raw_CNP[,k])                     #  Reads  total read counts
    stats[k,2]  = l                                       #   Bins  number of genome bins
    stats[k,3]  = round(mean(MX_raw_CNP[,k]), digits=2)   #   Mean  mean  read counts
    stats[k,4]  = round(  sd(MX_raw_CNP[,k]), digits=2)   #    Var  read count standard deviation
    stats[k,5]  = round(stats[k,4]/stats[k,3], digits=2)  #   Disp  read count dispersion
    stats[k,6]  =      min(MX_raw_CNP[,k])                #    Min    0% qtl
    stats[k,7]  = quantile(MX_raw_CNP[,k], c(.25))[[1]]   #   25th   25% qtl
    stats[k,8]  =   median(MX_raw_CNP[,k])                # Median   50% qtl
    stats[k,9]  = quantile(MX_raw_CNP[,k], c(.75))[[1]]   #   75th   75% qtl
    stats[k,10] =      max(MX_raw_CNP[,k])                #    Max  100% qtl
    ##########################################################################
    ### Segmentalize data
    ##########################################################################
    # Calculate normalized for current cell (previous values of normal seem wrong)
    normal[,k] = lowess.gc( MX_GC_cont[,1], (MX_raw_CNP[,k]+1)/mean(MX_raw_CNP[,k]+1) )

    # Compute log ratio between kth sample and reference
    if (CFIG_SegRef == 0) {  lr = log2( normal[,k])
    } else {          lr = log2((normal[,k])/(F))
    }

    # Determine breakpoints and extract chrom/locations
    CNA.obj = CNA(genomdat = lr, chrom = MX_bin_loc[,1], maploc = as.numeric(MX_bin_loc[,2]), data.type = 'logratio')
    CNA.smo = smooth.CNA(CNA.obj)
    CN_segs = segment(CNA.smo, verbose=0, min.width=minBinWidth)
    CN_frag = CN_segs$output[,2:3]

    # Map breakpoints to kth sample
    len = dim(CN_frag)[1]
    bps = array(0, len)
    for (j in 1:len) {
        bps[j]=which((MX_bin_loc[,1]==CN_frag[j,1]) & (as.numeric(MX_bin_loc[,2])==CN_frag[j,2]))
        }
    bps = sort(bps)
    bps[(len=len+1)] = l

    # Track global breakpoint locations
    breaks[bps,k] = 1

    # Calculate the indices for each segment
    index_list <- lapply(1:(length(bps) - 1), function(i) bps[i]:(bps[i + 1] - 1))
    index_list <- c(index_list, c(bps[length(bps)]:dim(normal)[1]))
    
    # Apply median to each segment and update 'fixed' accordingly
    fixed[,k] <- unlist(lapply(index_list, function(idx) {
                 rep(median(normal[, k][idx]), length(idx)) }))
    message(paste("cell median mean: ", mean(fixed[,k])))
    ### fixed[,k] <- fixed[,k]/mean(fixed[,k])
    ### this is causing skewness in CN Profile, will correct later
    ### when outputing copy number

    ##########################################################################
    ### Determine Copy Number (SoS Method)
    ##########################################################################
    outerRaw         = fixed[,k] %o% CNgrid
    outerRound       = round(outerRaw)
    outerDiff        = (outerRaw - outerRound) ^ 2
    outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
    CNmult[,k]       = CNgrid[order(outerColsums[,k])]
    CNerror[,k]      = round(sort(outerColsums[,k]), digits=2)

    ##########################################################################
    ### Select Best Copy Number
    ##########################################################################
    if (CFIG_FacsPD == 0 | length(which(lab[k]==ploidy[,1]))==0 ) {
        CN_fin_call = CNmult[1,k]
        FACS__CN <- 2
    } else if (CFIG_FacsPD == 1) {
        FACS__CN <- ploidy[which(lab[k]==ploidy[,1]),2]                     # input ploidy
        clost_CN <- which.min(abs(CNgrid - FACS__CN))                       # closest scanned ploidy
        CN_min_list <- c()
        for (i in seq(max(1,clost_CN-20), min(clost_CN+11, length(CNgrid)), 2)) {
            scan_win <- i:min(i+10, length(CNgrid)) # setup searching window
            pos_loc_min <- scan_win[which.min(outerColsums[scan_win,k])]    # find local minima
            #CN__loc_min <- CNgrid[pos_loc_min]                             # local minima closest to input ploidy
            CN_min_list <- c(CN_min_list, pos_loc_min)                      # re-assign chosen CN_fin_call
        }
        CN_min_list <- as.integer(names(table(CN_min_list))[table(CN_min_list)>1])
        MAD_lst <- c()
        for (i in CN_min_list) {
            fin_bin <- round(fixed[,k]*CNgrid[i])
            rnd_pld <- round(mean(fixed[,k])*CNgrid[i]*20)/20
            MAD_lst <- c(MAD_lst, median(abs(normal[,k]*rnd_pld-fin_bin)))}
        CN__loc_min <- CN_min_list[which.min(MAD_lst)]
        CN__closest <- CN_min_list[which.min(abs(CNgrid[CN_min_list]-FACS__CN))]
        CN_fin_call <- c(CNgrid[CN__loc_min], CNgrid[CN__closest])
        CNerror_facs = c(outerColsums[which(CNgrid==CN_fin_call[1]),k], outerColsums[which(CNgrid==CN_fin_call[2]),k])
        # If user specified FACS file, still calculate CNerror
        scale_CN = fixed[,k] %o% c(FACS__CN)
        SOS__err = round( sort(colSums((round(scale_CN) - scale_CN)^2, na.rm=FALSE, dims=1)), digits=2 )
    } else { # 
        estimate = ploidy[which(lab[k]==ploidy[,1]),2]
        CN_fin_call = CNmult[which(abs(CNmult[,k] - estimate)<.4),k][1]}

    ##########################################################################
    ### Record Result
    ##########################################################################
    fin_min[,k] = round(fixed[,k]*CN_fin_call[1])                  # integerized final CN profile
    fin_cls[,k] = round(fixed[,k]*CN_fin_call[2])                  # integerized final CN profile

    ### the deviation of segment median is corrected here 
    stats[k,11] = round(mean(fixed[,k])*CN_fin_call[1]*20)/20      # final ploidy as   local minima
    stats[k,13] = round(mean(fixed[,k])*CN_fin_call[2]*20)/20      # final ploidy as closest minima
    
    stats[k,12] = median(abs(normal[,k]*stats[k,11]-fin_min[,k]))  # MAD with CN   local minima
    stats[k,14] = median(abs(normal[,k]*stats[k,13]-fin_cls[,k]))  # MAD with CN closest minima
    print("finished")

    ##########################################################################
    ##########################################################################
    ##########################################################################
    ### Generate Plots & Figures
    ##########################################################################
    ##########################################################################
    ##########################################################################
    if( make_plots == TRUE){

    ##########################################################################
    ### Read Coverage Distribution
    ##########################################################################
    jpeg(filename=paste("1_covdist/", lab[k], "_1_dist.jpeg", sep=""), width=3000, height=750)
    
    upper_limit <- round(quantile(MX_raw_CNP[,k], c(.995))[[1]])
    chr_blk_evn <- data.frame(pos[seq(1,nrow(pos), 2),])
    chr_blk_odd <- data.frame(pos[seq(2,nrow(pos), 2),])
    dot_inbound <- data.frame(x=which(MX_raw_CNP[,k]<upper_limit), y=MX_raw_CNP[which(MX_raw_CNP[,k]<upper_limit),k])
    dot_otbound <- data.frame(x=which(MX_raw_CNP[,k]>upper_limit), y=array(upper_limit*.99, length(which(MX_raw_CNP[,k]>upper_limit))))
    chr___label <- data.frame(x=(pos[,2]+pos[,1])/2, y=-upper_limit*.05, chrom=substring(c(as.character(MX_bin_cnt[,1]), "chrY"), 4 ,5))

    plot1 = ggplot() +
        geom_rect( data=chr_blk_evn, aes(xmin=X1, xmax=X2, ymin=-upper_limit*.1, ymax=upper_limit), fill='gray85', alpha=0.75) +
        geom_rect( data=chr_blk_odd, aes(xmin=X1, xmax=X2, ymin=-upper_limit*.1, ymax=upper_limit), fill='gray75', alpha=0.75) + 
        geom_point(data=dot_inbound, aes(x=x, y=y), size=3) +
        geom_point(data=dot_otbound, aes(x=x, y=y), shape=5, size=6) +
        geom_text( data=chr___label, aes(x=x, y=y, label=chrom), size=12) +
        labs(title=paste("Genome Wide Read Distribution for Sample \"", lab[k], "\"", sep=""), x="Chromosome", y="Read Count", size=16) +
        theme(plot.title=element_text(size=36, vjust=1.5)) +
        theme(axis.title.x=element_text(size=40, vjust=-.1), axis.title.y=element_text(size=40, vjust=-.06)) +
        theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
        theme(panel.background = element_rect(fill = 'gray90')) +
        theme(plot.margin=unit(c(0.5, 1, 0.5, 1.5),"cm")) +
        theme(panel.grid.major.x = element_blank()) +
        scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
        scale_y_continuous(limits=c(-upper_limit*.1, upper_limit), expand = c(0, 0)) +
        geom_vline(xintercept = c(1, l), size=.25) +
        geom_hline(yintercept = c(-upper_limit*.1, upper_limit), size=.25)
        grid.arrange(plot1, ncol=1)
    dev.off()

    ##########################################################################
    ###  bin counts histogram
    ##########################################################################

    jpeg(filename=paste("2_bin_cnt/", lab[k], "_2_bin_cnt_hist.jpeg", sep=""), width=2500, height=1500)
        par(mar = c(7.0, 7.0, 7.0, 3.0))

        temp=sort(MX_raw_CNP[,k])[round(l*.01) : (l-round(l*.01))] 
        reads = hist(temp, breaks=100, plot=FALSE)
        plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), 
             xlab="Read Count (reads/bin)", cex.main=3, cex.axis=2, cex.lab=2)
        tu = par('usr')
        par(xpd=FALSE)
        clip(tu[1],                           mean(temp)-(diff(reads$mids)/2), tu[3], tu[4]) ; plot(reads, col='gray50', add=TRUE)
        clip(mean(temp)+(diff(reads$mids)/2), tu[2],                           tu[3], tu[4]) ; plot(reads, col='gray50', add=TRUE)
        clip(tu[1],                           mean(temp)-  sd(temp),           tu[3], tu[4]) ; plot(reads, col='gray75', add=TRUE)
        clip(mean(temp)+  sd(temp),           tu[2],                           tu[3], tu[4]) ; plot(reads, col='gray75', add=TRUE)
        clip(tu[1],                           mean(temp)-2*sd(temp),           tu[3], tu[4]) ; plot(reads, col='gray90', add=TRUE)
        clip(mean(temp)+2*sd(temp),           tu[2],                           tu[3], tu[4]) ; plot(reads, col='gray90', add=TRUE) 
        legend("topright", inset=.05, legend=c("mean", "< 1σ", "> 1σ", "> 2σ"), fill=c("black", "gray50", "gray75", "gray90"), cex=2.5)
    dev.off()

    ##########################################################################
    ###  lorenz curves
    ##########################################################################
    jpeg(filename=paste("3_lorenz/", lab[k], "_3_lorenz.jpeg", sep=""), width=2500, height=1500)

    nReads=sum(MX_raw_CNP[,k])
    uniq=unique(sort(MX_raw_CNP[,k]))
    
    lorenz=matrix(0, nrow=length(uniq), ncol=2)
    a=c(length(which(MX_raw_CNP[,k]==0)), tabulate(MX_raw_CNP[,k], nbins=max(MX_raw_CNP[,k])))
    b=a*(0:(length(a)-1))
    for (i in 2:length(uniq)) {
        lorenz[i,1]=sum(a[1:uniq[i]])/l
        lorenz[i,2]=sum(b[2:uniq[i]])/nReads }

    # smooth.spline needs >= 4 points...
    fit = data.frame(x=lorenz[,1], y=lorenz[,2])
    if(nrow(lorenz) >= 4)
    {   spline = try(smooth.spline(lorenz))
        if(class(spline) != "try-error")
            fit = data.frame(x=spline$x, y=spline$y)}

    perf=data.frame(x=c(0,1), y=c(0,1))

    plot1 = try(ggplot() +
        geom_line(data=perf, aes(x=x, y=y, color="Perfect Uniformity"), size=3) +
        geom_line(data=fit, aes(x=x, y=y, color="Sample Uniformity"), size=3) +
        scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
        labs(title=paste("Lorenz Curve of Coverage Uniformity for Sample ", lab[k], sep=""), x="Cumulative Fraction of Genome", y="Cumulative Fraction of Total Reads") +
        theme(plot.title=element_text(size=45, vjust=1.5)) +
        theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
        theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
        theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
        theme(panel.background = element_rect(color = 'black')) +
        theme(legend.title=element_blank(), legend.text=element_text(size=40)) +
        theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
        theme(legend.position=c(.15, .85)) +
        scale_color_manual(name='', values=c('Perfect Uniformity'="black", 'Sample Uniformity'=colors[COL_PALETTE,1])))

        grid.arrange(plot1, ncol=1)
    dev.off()

    ##########################################################################
    ###  GC correction
    ##########################################################################
    #jpeg(filename=paste("4_GC_lowe/", lab[k], "_4_GC_lowess.jpeg", sep=""), width=2500, height=1250)
    pdf(paste("4_GC_lowe/", lab[k], "_4_GC_lowess.pdf",  sep=""), width=2500/100, height=1250/100)

    low = lowess(MX_GC_cont[,1], log(normal2[,k]), f=0.05)
    app = approx(low$x, low$y, MX_GC_cont[,1])
    cor = exp(log(normal2[,k]) - app$y)
    
    uncorrected = data.frame(x=MX_GC_cont[,1], y=log(normal2[,k]))
    corrected = data.frame(x=MX_GC_cont[,1], y=log(cor))
    fit = data.frame(x=app$x, y=app$y)

    try(plot1 <- ggplot() +
        geom_point(data=uncorrected, aes(x=x, y=y), size=3) +
        geom_line(data=fit, aes(x=x, y=y, color="Lowess Fit"), size=3) +
        scale_x_continuous(limits=c(min(.3, min(MX_GC_cont[,1])), max(.6, max(MX_GC_cont[,1]))), breaks=seq(.3,.6,.05)) +
        labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Uncorrected)", sep=""), x="GC content", y="Normalized Read Counts (log scale)") +
        theme(plot.title=element_text(size=25, vjust=1.5)) +
        theme(axis.title.x=element_text(size=25, vjust=-2.8), axis.title.y=element_text(size=25, vjust=.1)) +
        theme(axis.text=element_text(color="black", size=25), axis.ticks=element_line(color="black")) +
        theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
        theme(panel.background = element_rect(color = 'black')) +
        theme(legend.title=element_blank(), legend.text=element_text(size=25)) +
        theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
        theme(legend.position=c(.85, .9)) +
        scale_color_manual(name='', values=colors[COL_PALETTE,1]))

    try(plot2 <- ggplot() +
        geom_point(data=corrected, aes(x=x, y=y), size=3) +
        scale_x_continuous(limits=c(min(.3, min(MX_GC_cont[,1])), max(.6, max(MX_GC_cont[,1]))), breaks=seq(.3,.6,.05)) +
        labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Corrected)", sep=""), x="GC content", y="") +
        theme(plot.title=element_text(size=25, vjust=1.5)) +
        theme(axis.title.x=element_text(size=25, vjust=-2.8), axis.title.y=element_text(size=25, vjust=.1)) +
        theme(axis.text=element_text(color="black", size=25), axis.ticks=element_line(color="black")) +
        theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
        theme(panel.background = element_rect(color = 'black')))
        try(grid.arrange(plot1, plot2, ncol=2))
    dev.off()

    ##########################################################################
    ###  bin CN_fin_call histogram (i.e. normalized bin counts)
    ##########################################################################
    jpeg(filename=paste("5_CN_hist/", lab[k], "_5_hist.jpeg", sep=""), width=2500, height=1500)

    bin_readcnt=data.frame(x=normal[,k]*CN_fin_call)
    plot1 = ggplot() +
        geom_histogram(data=bin_readcnt, aes(x=x), binwidth=.05, color="black", fill="gray60") +
        geom_vline(xintercept=seq(0,10,1), size=1, linetype="dashed", color=colors[COL_PALETTE,1]) +
        scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1)) +
        labs(title=paste("Frequency of Bin Counts for Sample \"", lab[k], "\"\nNormalized and Scaled by Predicted CN_fin_call (", CNmult[1,k], ")", sep=""), x="Copy Number", y="Frequency") +
        theme(plot.title=element_text(size=45, vjust=1.5)) +
        theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
        theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
        theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
        theme(panel.background = element_rect(color = 'black'))
        grid.arrange(plot1, ncol=1)
    dev.off()

    ##########################################################################
    ###  SoS plots (i.e. um of squares error across potential CN_grid)
    ##########################################################################
    #jpeg(filename=paste("6_SoS/", lab[k], "_6_SoS.jpeg", sep=""), width=2500, height=1500)
    pdf(paste("6_SoS/", lab[k], "_6_SoS.pdf", sep=""), width=3000/120, height=1500/120)

    top = max(outerColsums[,k])
    sosDat = data.frame(x=CNgrid, y=outerColsums[,k])
    lim = cbind(c(seq(0,5000,500), 1000000), c(50, 100, 100, 200, 250, 400, 500, 500, 600, 600, 750, 1000))
    step = lim[which(top<lim[,1])[1],]
    minSoS = data.frame(x=CNmult[1,k], y=CNerror[1,k]) ; global_min='MinErr ploidy'
    # If a FACS file is provided, use CNerror_facs, 
    # since CN_fin_call multiplier could be outside the CNgrid range, 
    # which would cause "which(CNgrid==CN_fin_call)" to error out
    if(CFIG_FacsPD == 1) {
        bestSoS = data.frame(x=c(FACS__CN, CN_fin_call), y=c(SOS__err, CNerror_facs), color=c("Input  Ploidy", "Local_Minimal", "Closest___Min")) 
        bestSoS[2,2] = bestSoS[2,2]*1.01 ; bestSoS[3,2] = bestSoS[3,2]*0.99
    } else {
        bestSoS = data.frame(x=CN_fin_call, y=outerColsums[which(CNgrid==CN_fin_call),k], color="Chosen Ploidy") }

        lg1 = paste0('MinErr ploidy: ', minSoS$x)
        lg2 = paste0('Chosen Ploidy: ', bestSoS[1,1])
        lg3 = paste0("Input  Ploidy: ", bestSoS[2,1])
    limit = c( min(minPloidy, CN_fin_call), max(maxPloidy, CN_fin_call))
    plot1 = ggplot() +
        geom_line( data=sosDat,  aes(x=x, y=y),                    size=3) +
        geom_point(data=sosDat,  aes(x=x, y=y), fill="black",      size=5,  shape=21) +
        geom_point(data=minSoS,  aes(x=x, y=y,  color=global_min), size=15, shape=16, alpha=0.7) +
        geom_point(data=bestSoS, aes(x=x, y=y,  color=color),      size=15, shape=18, alpha=0.7) +
        scale_x_continuous(limits=limit, breaks=seq(limit[1], limit[2], .5)) +
        labs(title="Sum of Squares Error Across Potential Copy Number States", 
             x="Copy Number Multiplier", y="Sum of Squares Error") +
        theme(plot.title=  element_text(size=45, vjust=1.5)) +
        theme(axis.title.x=element_text(size=45, vjust=-2.8),   axis.title.y=element_text(size=45, vjust=.1)) +
        theme(axis.text=  element_text(color="black", size=45), axis.ticks=  element_line(color="black")) +
        theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
        theme(panel.background = element_rect(color = 'black')) +
        theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
        theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
        theme(legend.position=c(.85, .9)) +
        scale_color_manual(name='', values=c("MinErr ploidy"="#377eb8", 
                                             "Input  Ploidy"="#984ea3",
                                             "Local_Minimal"="#e41a1c", 
                                             "Closest___Min"="#ff7f0e"))
        ### scale_color_manual(name='', values=c(paste0('MinErr ploidy: ', minSoS$x)=colors[COL_PALETTE,1], 
        ###                                      paste0('Chosen Ploidy: ', bestSoS[1,1])=colors[COL_PALETTE,2], 
        ###                                      paste0("Input  Ploidy: ", bestSoS[2,1])=colors[2,1]))

        grid.arrange(plot1, ncol=1)
    dev.off()

    ##########################################################################
    ###  colored CN_fin_call profile
    ##########################################################################

    #jpeg(filename=paste(lab[k], "_CN.jpeg", sep=""), width=3000, height=750)
    pdf(paste("7_CN_Profile/", lab[k], "_7_CN_profile_local_minimal.pdf", sep=""), width=3000/120, height=750/120)

    top__ploidy <- 5   ### use 12 for larger range
    chr_blk_evn <- data.frame(pos[seq(1,nrow(pos), 2),])
    chr_blk_odd <- data.frame(pos[seq(2,nrow(pos), 2),])
    chr_ID_anno <- data.frame(x=(pos[,2]+pos[,1])/2, y=-top__ploidy*.05, chrom=substring(c(as.character(MX_bin_cnt[,1]), "chrY"), 4 ,5))
    bin_readcnt <- data.frame(x=1:l, y=normal[,k]*CN_fin_call[1])
    CNV_bin_amp <- data.frame(x=which(fin_min[,k] >2), y=fin_min[which(fin_min[,k] >2),k])
    CNV_bin_del <- data.frame(x=which(fin_min[,k] <2), y=fin_min[which(fin_min[,k] <2),k])
    CNV_bin_flt <- data.frame(x=which(fin_min[,k]==2), y=fin_min[which(fin_min[,k]==2),k])

    plot1 = ggplot() +
        geom_rect( data=chr_blk_evn, aes(xmin=X1, xmax=X2, ymin=-top__ploidy*.1, ymax=top__ploidy), fill="gray85", alpha=0.75) +
        geom_rect( data=chr_blk_odd, aes(xmin=X1, xmax=X2, ymin=-top__ploidy*.1, ymax=top__ploidy), fill="gray75", alpha=0.75) +
        geom_point(data=bin_readcnt, aes(x=x, y=y), size=3, color="gray45") +
        geom_point(data=CNV_bin_flt, aes(x=x, y=y), size=4, color="black") +
        geom_point(data=CNV_bin_amp, aes(x=x, y=y), size=4, color=colors[COL_PALETTE,1]) +
        geom_point(data=CNV_bin_del, aes(x=x, y=y), size=4, color=colors[COL_PALETTE,2]) +
        geom_text( data=chr_ID_anno, aes(x=x, y=y, label=chrom), size=12) +
        scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
        scale_y_continuous(limits=c(-top__ploidy*.1, top__ploidy), expand = c(0, 0)) +
        labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n", 
                         "Input Ploidy = ", FACS__CN, "; Predicted Ploidy = ", stats[k,11], "; MAD = ", stats[k,12], sep=""), 
             x="Chromosome", y="Copy Number", size=16) +
        theme(plot.title=  element_text(size=40, vjust=1.5)) +
        theme(axis.title.x=element_text(size=40, vjust=-.05),    axis.title.y=element_text(size=40, vjust=.1)) +
        theme(axis.text=   element_text(color="black", size=40), axis.ticks=  element_line(color="black"))+
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
        theme(panel.background = element_rect(fill = 'gray90')) +
        theme(plot.margin=unit(c(.5,1,.5,1),"cm")) +
        theme(panel.grid.major.x = element_blank()) +
        geom_vline(xintercept = c(1, l), size=.5) +
        geom_hline(yintercept = c(-top__ploidy*.1, top__ploidy), size=.5)
        grid.arrange(plot1, ncol=1)
    dev.off()
    
    
    ##########################################################################
    ###  colored CN_fin_call profile
    ##########################################################################
    #jpeg(filename=paste(lab[k], "_CN.jpeg", sep=""), width=3000, height=750)
    pdf(paste("7_CN_Profile/", lab[k], "_7_CN_profile_closest_minimal.pdf", sep=""), width=3000/120, height=750/120)

    chr_blk_evn <- data.frame(pos[seq(1,nrow(pos), 2),])
    chr_blk_odd <- data.frame(pos[seq(2,nrow(pos), 2),])
    chr_ID_anno <- data.frame(x=(pos[,2]+pos[,1])/2, y=-top__ploidy*.05, chrom=substring(c(as.character(MX_bin_cnt[,1]), "chrY"), 4 ,5))
    bin_readcnt <- data.frame(x=1:l, y=normal[,k]*CN_fin_call[2])
    CNV_bin_amp <- data.frame(x=which(fin_cls[,k] >2), y=fin_cls[which(fin_cls[,k] >2),k])
    CNV_bin_del <- data.frame(x=which(fin_cls[,k] <2), y=fin_cls[which(fin_cls[,k] <2),k])
    CNV_bin_flt <- data.frame(x=which(fin_cls[,k]==2), y=fin_cls[which(fin_cls[,k]==2),k])

    plot1 = ggplot() +
        geom_rect( data=chr_blk_evn, aes(xmin=X1, xmax=X2, ymin=-top__ploidy*.1, ymax=top__ploidy), fill="gray85", alpha=0.75) +
        geom_rect( data=chr_blk_odd, aes(xmin=X1, xmax=X2, ymin=-top__ploidy*.1, ymax=top__ploidy), fill="gray75", alpha=0.75) +
        geom_point(data=bin_readcnt, aes(x=x, y=y), size=3, color="gray45") +
        geom_point(data=CNV_bin_flt, aes(x=x, y=y), size=4, color="black") +
        geom_point(data=CNV_bin_amp, aes(x=x, y=y), size=4, color=colors[COL_PALETTE,1]) +
        geom_point(data=CNV_bin_del, aes(x=x, y=y), size=4, color=colors[COL_PALETTE,2]) +
        geom_text( data=chr_ID_anno, aes(x=x, y=y, label=chrom), size=12) +
        scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
        scale_y_continuous(limits=c(-top__ploidy*.1, top__ploidy), expand = c(0, 0)) +
        labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n", 
                         "Input Ploidy = ", FACS__CN, "; Predicted Ploidy = ", stats[k,11], "; MAD = ", stats[k,12], sep=""), 
             x="Chromosome", y="Copy Number", size=16) +
        theme(plot.title=  element_text(size=40, vjust=1.5)) +
        theme(axis.title.x=element_text(size=40, vjust=-.05),    axis.title.y=element_text(size=40, vjust=.1)) +
        theme(axis.text=   element_text(color="black", size=40), axis.ticks=  element_line(color="black"))+
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
        theme(panel.background = element_rect(fill = 'gray90')) +
        theme(plot.margin=unit(c(.5,1,.5,1),"cm")) +
        theme(panel.grid.major.x = element_blank()) +
        geom_vline(xintercept = c(1, l), size=.5) +
        geom_hline(yintercept = c(-top__ploidy*.1, top__ploidy), size=.5)
        grid.arrange(plot1, ncol=1)
    dev.off()
}}})
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
### Save processed data
####################################################################################################
# Update FILE_runlog
statusFile=file( paste(PATH_2_work, "/", FILE_runlog, sep="") )
writeLines(c("<?xml version='1.0'?>", "<FILE_runlog>", "<step>3</step>", 
             paste("<processingfile>Saving Data</processingfile>", sep=""), 
             paste("<percentdone>", (w*100)%/%(w+4), "</percentdone>", sep=""), 
             "<tree>clust.xml</tree>", "</FILE_runlog>"), statusFile)
close(statusFile)

# Close output stream
sink()

# Store processed sample information
loc2=MX_bin_loc
loc2[,3]=loc2[,2]
pos = cbind(c(1,MX_bin_cnt[,2]), c(MX_bin_cnt[,2], l))

# 
for (i in 1:nrow(pos)) {
    if( (pos[i,2] - pos[i,1]) == 0 ) {# If only 1 bin in a chromosome
        loc2[pos[i,1],1] = 1
    } else if( (pos[i,2] - pos[i,1]) == 1 ) {# If only 2 bins
        loc2[pos[i,1],1] = 1
        loc2[pos[i,2],1] = loc2[pos[i,1],2] + 1
    } else {
        loc2[pos[i,1]:(pos[i,2]-1),2]=c(1,MX_bin_loc[pos[i,1]:(pos[i,2]-2),2]+1)
}}

loc2[nrow(loc2),2]=loc2[nrow(loc2)-1,3]+1
colnames(loc2)=c("CHR","START", "END")

write.table(stats,                file=paste(PATH_2_work, "/zout_01_summary_stats.tsv",     sep=""), sep="\t", quote=FALSE)
write.table(cbind(loc2,  normal), file=paste(PATH_2_work, "/zout_02_normalized_rr.tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2,  breaks), file=paste(PATH_2_work, "/zout_03_break__points.tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2,   fixed), file=paste(PATH_2_work, "/zout_04_smoothened_rr.tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2, fin_min), file=paste(PATH_2_work, "/zout_05_CN__local_min.tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2, fin_cls), file=paste(PATH_2_work, "/zout_06_CN__close_min.tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
