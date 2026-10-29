#This file contains scripts for the analysis of oCelloscope measurements

#This function loads the ocelloscope data into R
#Outputs a list where each item contains a particular feature (e.g. volume, area, no. tips, ...)
#Needs the file name as input, and start = how many lines to skip from the beginning, repetions = number of ocelloscope timepoints
load.ocelloscope.data <- function(datafile, start, repetitions, sep, dec) {
    foo <- rep(list(0), 5)
    for(i in 1:5) { #Loop over the features (volume, area, no. tips, tips / area, timestamp)
        foo[[i]] <- read.csv(file = datafile, header = T, sep = sep, dec = dec, skip = start + (i-1)*(repetitions+3), nrows = repetitions, fileEncoding = "UTF-16")
        colnames(foo[[i]])[1] <- "Time" #Fix the first column name, unit is seconds
    }
    names(foo) <- c("volume", "area", "no.tips", "tips.per.area", "timestamp")
    return(foo) #Return results
}

load.ocelloscope.data2 <- function(datafile, repetitions, sep, dec, features = c("EstimatedVolume", "HyphaeArea", "NumberOfTips", "TipsPerArea", "TimestampInSeconds"), enc = "UTF-16") {
    #First get the lines where the data for the features start
    con <- file(datafile, encoding = enc) #File has UTF-16 encoding
    #Need to search the file for the feature names
    linebyline <- readLines(con) #Reading the file line by line
    startlines <- rep(0, length(features))
    start.feat.line <- grep("Features", linebyline)
    feat.data <- linebyline[start.feat.line:length(linebyline)] #Line index to extract all features
    for(i in 1:length(features)) { startlines[i] <- which(features[i] == feat.data) }
        #ifelse(length(grep(features[i], feat.data)) > 1, grep(features[i], feat.data)[2], grep(features[i], feat.data))  }
    foo <- rep(list(0), 5)
    for(i in 1:5) { #Loop over the features (volume, area, no. tips, tips / area, timestamp)
        foo[[i]] <- read.csv(file = datafile, header = T, sep = sep, dec = dec, skip = start.feat.line - 1 + startlines[i], nrows = repetitions, fileEncoding = enc)
        colnames(foo[[i]])[1] <- "Time" #Fix the first column name, unit is seconds
    }
    names(foo) <- c("volume", "area", "no.tips", "tips.per.area", "timestamp")
    close(con) #Close the connections
    return(foo) #Return results
}

        
#This function combines and converts ocelloscope data into a long format
convert.ocelloscope.data.long <- function(inputdata) {
    features <- names(inputdata) #Names of the features
    nfeat <- length(features) #Number of features
    for(i in 1:nfeat) { #loop over the features
        if(i == 1) { #Set up the data frame on first iteration
            longdata <- pivot_longer(inputdata[[i]], colnames(inputdata[[i]][-c(1,2)]), names_to = "sample", values_to = features[i])
        }
        if(i > 1 & features[i] != "timestamp") { #Add the other features to the data
            longdata <- cbind(longdata, pivot_longer(inputdata[[i]], colnames(inputdata[[i]][-c(1,2)]), names_to = "sample", values_to = features[i])[,4])
        }
        if(i > 1 & features[i] == "timestamp") { #Add the final timestamp to the data
                longdata <- cbind(longdata, pivot_longer(inputdata[[i]], colnames(inputdata[[i]][-1]), names_to = "sample", values_to = features[i])[,3])
        }
    } #Done looping over features
    return(longdata)
} #Done


#This function estimates growth rates from ocelloscope data
#limitmatrix contains limits for the different growth rate windows (start, end)
#default limits are length of 15000, and window is moved by 5000
ocelloscope.growth.rate <- function(longdata, limitmatrix = matrix(c(0, 15000, 5000, 20000, 10000, 25000, 15000, 30000, 20000, 35000, 25000, 40000, 30000, 45000, 35000, 50000, 40000, 55000, 45000, 60000), ncol = 2, byrow = T)) {
    samples <- unique(longdata$sample) #sample names
    n.samples <- length(samples) #Get number of samples
    n.windows <- nrow(limitmatrix)
    #Make the results matrix
    results <- data.frame(sample = samples, gr = rep(0, n.samples), R2 = rep(0, n.samples), window = rep(0, n.samples), tips.at.limit1 = rep(0, n.samples) )
    #
    #Calculating growth rates
    for(i in 1:n.samples) { #Loop over all samples and fit area against time
        current <- filter(longdata, sample == samples[i]) #Take the current samples
        temp.gr <- matrix(rep(0, n.windows*4), ncol = 4) #Make a temporary gr results matrix
        current$timestamp <- current$timestamp / 3600 #Change seconds to hours
        #
        for(j in 1:n.windows) { #Loop over windows to measure gr
                current.w <- filter(current, Time >= limitmatrix[j,1] & Time <= limitmatrix[j,2]) #Take only the linear part
                #Need to check that there are hyphal tips in the data
                if( current.w$no.tips[1] == 0 | is.na(current.w$no.tips[1]) == T ) { temp.gr[j,1] <- NA; temp.gr[j,2] <- NA } else {
            model <- lm(area ~ timestamp, data = current.w) #Fit the growth rate
            temp.gr[j,1] <- coef(model)[2] #Slope of the fit, i.e. growth rate
            temp.gr[j,2] <- summary(model)$adj.r.squared #Take the R^2 value of the fit
            temp.gr[j,3] <- current.w$no.tips[1] #Store the number of tips at start
            temp.gr[j,4] <- j #Store the window index
                }
        } #Done looping over different time windows
        ### STORING THE RESULTS ###################################################
        #Now need to check which fit is the best
        #Which limits give the best fit?
        #Which limits give the highest growth rate?
        ## Rules: first look which limits give highest R2, if multiple limits have fits where R2 > 99.9, then select that which gives the highest growth rate
        ## Negative growth rates are ignored
        neg.ind <- temp.gr[,1] < 0
        if (any(neg.ind, na.rm = TRUE)) { temp.gr[neg.ind,2] <- NA } #Set negative gr R2 as NA 
        max.R2.ind <- which(temp.gr[,2] == max(temp.gr[,2], na.rm = T)) #Index of max R2
        #browser()
        #print(paste("storing results for sample", samples[i], sep = " "))
        #Here need to deal with empty wells, which have NA for gr and R2
        if(all(is.na(temp.gr[,2])) == T) {
            results$gr[i] <- NA
            results$R2[i] <- NA
            results$tips.at.limit1[i] <- NA
            results$window[i] <- NA } else {
            #
                if(sum(temp.gr[,2] > 0.995, na.rm = T) >= 2) { #If there are more fits with R2 > 99.5, then look which has the highest growth rate
                    temp2 <- temp.gr[(temp.gr[,2] > 0.995),] #Those that have R2 > 99.5
                    max.gr.ind <- which(temp2[,1] == max(temp2[,1], na.rm = T)) #Index of max gr
                    results$gr[i] <- temp2[max.gr.ind,1] #store the growth rate
                    results$R2[i] <- temp2[max.gr.ind,2] #store the R2 of the model fit
                    results$tips.at.limit1[i] <- temp2[max.gr.ind,3] #store number of hyphal tips at start
                    results$window[i] <- paste(limitmatrix[temp2[max.gr.ind,4],], collapse = ":")
                } else { #If there is a clear best window
                    results$gr[i] <- temp.gr[max.R2.ind,1] #store the growth rate
                    results$R2[i] <- temp.gr[max.R2.ind,2] #store the R2 of the model fit
                    results$tips.at.limit1[i] <- temp.gr[max.R2.ind,3] #store number of hyphal tips at start
                    results$window[i] <- paste(limitmatrix[temp.gr[max.R2.ind,4],], collapse = ":")
                } #Done storing results
            } #Close NA check
        ############################################################################    
       # 
      #Done calculating growth rate for a single sample
    } #Done looping over all samples
    #
    return(results)
} #Done
