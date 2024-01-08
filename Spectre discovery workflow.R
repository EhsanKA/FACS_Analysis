## THIS WORK IS ADAPTED FROM https://doi.org/10.1002/cyto.a.24350
##########################################################################################################
#### Spectre - Simple Discovery Workflow
##########################################################################################################

    # Spectre R package: https://immunedynamics.io/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### Load libraries

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages
        library(data.table)

    ### Set primary directory
        
        d <- getwd() # Gets the current working directory
        # Set working directory location path which stores the data and metadata
        # Input directory
        PrimaryDirectory <- setwd(d) # CHANGE TO YOUR PERSONAL LINK HERE!!
        PrimaryDirectory
        

    ### Create output directory
        
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################

    ### Import data

        list.files('./data/', ".csv")

        data.list <- Spectre::read.files(file.loc = './data/',
                                         file.type = ".csv",
                                         do.embed.file.names = TRUE)

    ### Check the data

        check <- do.list.summary(data.list)

        check$name.table # Review column names and their subsequent values
        check$ncol.check # Review number of columns (features/markers) in each sample
        check$nrow.check # Review number of rows (cells) in each sample

        data.list[[1]]

    ### Merge data

        cell.dat <- Spectre::do.merge.files(dat = data.list)
        cell.dat
        
        cell.dat$FileName %>% table

    ### Read in metadata  
       
        meta.dat <- fread("./metadata/sample.details.csv")
        meta.dat
        
##########################################################################################################
#### 3. Data transformation
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 1 - transformed plots")
    setwd("Output 1 - transformed plots")
        
    ### Arcsinh transformation

        as.matrix(names(cell.dat))

        to.asinh <- names(cell.dat)[c(1:9)] # Only use the 9 feature columns for transformation
        to.asinh

        cofactor <- 500

        cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
        transformed.cols <- paste0(to.asinh, "_asinh")
        
        head(cell.dat)
        
        plot.against <- "Ly6C_asinh"
        
        for(i in transformed.cols){
          make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
        }

##########################################################################################################
#### 4. Add metadata and set some preferences
##########################################################################################################

    ### Add metadata to data.table

        meta.dat
        sample.info <- meta.dat[,c(1:4)]
        sample.info
        
        meta.dat
        counts <- meta.dat[,c(2,5)]
        counts

        cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "Filename", rmv.ext = TRUE)
        cell.dat

    ### Columns

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(12:20)]
        as.matrix(cellular.cols)

        cluster.cols <- names(cell.dat)[c(12:20)]
        as.matrix(cluster.cols)

        exp.name <- "CNS experiment"
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - clustering")
    setwd("Output 2 - clustering")

    ### Clustering

        cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 8)
        cell.dat
        
    ### Dimensionality reduction
        
        # Subsample per group/FlowSOM_metacluster/FileName
        cell.dat$FileName %>% table %>% as.data.frame()
        cell.sub <- cell.dat %>% group_by(FileName) %>% slice_sample(n=1000) %>% as.data.table
        
        cell.sub <- run.umap(cell.sub, cluster.cols)
        cell.sub

    ### DR plots

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE, nudge_x = 0, nudge_y = 0)
        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', add.label = FALSE)
        # Feature plots
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')

    ### Expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
        make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols)

##########################################################################################################
#### 6. Annotate clusters
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 3 - annotation")
    setwd("Output 3 - annotation")

    ### Annotate

        annots <- list("CD4 T cells" = c(3),
                       "CD8 T cells" = c(2),
                       "NK cells" = c(1),
                       "Neutrophils" = c(8),
                       "Infil Macrophages" = c(4),
                       "Microglia" = c(5,6,7)
        )

        annots <- do.list.switch(annots)
        names(annots) <- c("Values", "Population")
        setorderv(annots, 'Values')
        annots

    ### Add annotations

        cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
        cell.dat
        
        cell.dat$Population %>% table %>% as.data.frame
        
        # Additionally adding to the subsampled dataset for viz.
        cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
        cell.sub
        

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", col.type = 'factor', add.label = TRUE, nudge_x = 0, nudge_y = 0)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", group.col, col.type = 'factor')

    ### Expression heatmap
        
        rm(exp)
        exp <- do.aggregate(cell.dat, cellular.cols, by = "Population")
        make.pheatmap(exp, "Population", cellular.cols)
        
    ### Write FCS files
        
        # setwd(OutputDirectory)
        # setwd("Output 3 - annotation")
        
        fwrite(cell.dat, "Annotated.data.csv")
        fwrite(cell.sub, "Annotated.data.DR.csv")
##########################################################################################################
#### 7. Summary data and statistical analysis
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4 - summary data")
    setwd("Output 4 - summary data")

    ### Setup
    
        variance.test <- 'kruskal.test'
        pairwise.test <- "wilcox.test"
    
        comparisons <- list(c("Mock", "WNV"))
        comparisons
        
        grp.order <- c("Mock", "WNV")
        grp.order
    
    ### Select columns to measure MFI
    
        as.matrix(cellular.cols)
        dyn.cols <- cellular.cols[c(5,8)]
        dyn.cols
    
    ### Create summary tables
    
        sum.dat <- create.sumtable(dat = cell.dat, 
                                   sample.col = sample.col,
                                   pop.col = "Population",
                                   use.cols = dyn.cols, 
                                   annot.cols = c(group.col, batch.col), 
                                   #counts = counts
                                   )
    ### Alternatively,

        tmp <- cell.dat %>% group_by(Sample,Group,Batch,Population)
        foo <- tmp %>% tally()#sort=TRUE
        foo <- as.data.frame(foo)
        View(foo)
        setDT(foo)[, frac := n / sum(n) * 100, by=Sample]
        
        library(RColorBrewer)
        n <- cell.dat$Group %>% unique %>% length
        n
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        
        library(ggplot2)
        ggplot(foo, aes(fill=Group, y=n, x=Population)) +
          geom_bar(position="fill", stat="identity", colour = 'NA', size = 0.5) + scale_fill_manual("legend", values = sample(col_vector, n)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          xlab('Population')
        
    ### Review summary data
        
        sum.dat
        as.matrix(names(sum.dat))
        
        annot.cols <- c(group.col, batch.col)
        
        plot.cols <- names(sum.dat)[c(4:21)]
        plot.cols
        
    ### Reorder summary data and SAVE
        
        sum.dat <- do.reorder(sum.dat, group.col, grp.order)
        sum.dat[,c(1:3)]
        
        fwrite(sum.dat, 'sum.dat.csv')
        
    ### Autographs

        for(i in plot.cols){
            
            measure <- gsub("\\ --.*", "", i)
            measure
            
            pop <- gsub("^[^--]*.-- ", "", i)
            pop
            
            make.autograph(sum.dat,
                           x.axis = group.col,
                           y.axis = i,
                           y.axis.label = measure,
                           
                           grp.order = grp.order,
                           my_comparisons = comparisons,
                           
                           Variance_test = variance.test,
                           Pairwise_test = pairwise.test,
                           
                           title = pop,
                           subtitle = measure,
                           filename = paste0(i, '.pdf'))
            
        }
        
    ### Create a fold change heatmap
        
        ## Z-score calculation
        sum.dat.z <- do.zscore(sum.dat, plot.cols)
        
        ## Group 
        t.first <- match(grp.order, sum.dat.z[[group.col]])
        t.first <- t.first -1
        t.first
        
        ## Make heatmap
        make.pheatmap(sum.dat.z, 
                      sample.col = sample.col, 
                      plot.cols = paste0(plot.cols, '_zscore'), 
                      is.fold = TRUE, 
                      plot.title = 'Z-score',
                      annot.cols = annot.cols,
                      dendrograms = 'column',
                      row.sep = t.first,
                      cutree_cols = 4)

