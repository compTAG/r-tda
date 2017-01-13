#########################################################################
## Modified from
## https://www.r-bloggers.com/finally-tracking-cran-packages-downloads/
## https://github.com/talgalili/installr/blob/master/R/RStudio_CRAN_data.R
#########################################################################


## ======================================================================
## Step 1: Download all log files
## ======================================================================
downloadRStudioCRAN <- function(start, end, logDirectory) {
  
  # Here's an easy way to get all the URLs in R
  start <- as.Date(start)
  end <- as.Date(end)
  all_days <- seq(start, end, by = 'day')
  
  # only download the files you don't have:
  missing_days <- setdiff(as.character(all_days),
                          tools::file_path_sans_ext(dir(logDirectory), TRUE))
  
  year <- as.POSIXlt(missing_days)[["year"]] + 1900
  if (length(missing_days) > 0) {
    urls <- paste0("http://cran-logs.rstudio.com/", year, "/", missing_days,
                   ".csv.gz")
  } else {
    urls <- character(0)
  }
  
  if (!dir.exists(logDirectory)) {
    dir.create(logDirectory)
  }
  for (i in seq_along(missing_days)) {
    cat(i, "/", length(missing_days), "\n", sep = "")
    download.file(urls[i], 
                  paste0(logDirectory, "/", missing_days[i], ".csv.gz"))
  }
  
}


## ======================================================================
## Step 2: Load single data files into one big data.table
## ======================================================================
readRStudioCRAN <- function(start, end, logDirectory, packages) {

  file_list <- list.files(logDirectory, full.names = TRUE)
  
  # include only the relevant type of files, such as: "2013-04-02.csv.gz"
  file_list <- file_list [ grep("[0-9]+-[0-9]+-[0-9]+.csv.gz", file_list)]

  # removes empty files
  file_list_info <- file.info(file_list)
  ss_non_0_files <- file_list_info[["size"]] > 0 
  file_list <- file_list[ss_non_0_files]   
    
  all_days <- seq(start, end, by = 'day')
  dayIndex <- rep(FALSE, length(file_list))
  for (day in all_days) {
    dayIndex <- dayIndex | grepl(as.Date(day, origin='1970-01-01'), file_list)
  }
  file_list <- file_list[dayIndex]
  
  # read files
  logs <- list()
  for (file in file_list) {
    cat("Reading", file, "...\n")
    flush.console()
    logfile <- read.table(file, header = TRUE, sep = ",", quote = "\"",
                         dec = ".", fill = TRUE, stringsAsFactors = FALSE,
                         comment.char = "", as.is = TRUE)
    logs[[file]] <- subset(logfile, logfile[["package"]] %in% packages)
  }
  rm(logfile)
  gc()
  
  # rbind together all files
  dataset <- data.table::rbindlist(logs)
  rm(logs)
  gc()
  
  # add some keys and define variable types
  dataset[, date:=as.Date(date)]
  dataset[, package:=factor(package)]
  dataset[, country:=factor(country)]
  dataset[, weekday:=weekdays(date)]
  dataset[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]
  
  setkey(dataset, package, date, week, country)
  
  return(dataset)
}


## ======================================================================
## Step 3: Analyze it!
## ======================================================================
countPackage <- function(dataset) {

  cat(str(dataset))
  
  # Overall downloads of packages 
  d1 <- dataset[, length(week), by=package]
  d1 <- d1[order(V1), ]
  for (packageOne in packages) {
    print(d1[package == packageOne, ])
    if (empty(d1[package == packageOne, ])) {
      packages <- setdiff(packages, packageOne)
    }
  }
  
  countDataset <- list()
  
  # plot: Compare downloads of selected packages on a weekly basis
  countDataset[["agg"]] <- 
      dataset[J(packages), length(unique(ip_id)), by=c("week", "package")]
  
  countDataset[["plot"]] <- ggplot(countDataset[["agg"]],
      aes(x = week, y = V1, color = package, group = package)) +
      geom_line() + ylab("Downloads") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))
  
  return(countDataset)
}



#########################################################################
## Set these parameters as your R environment
#########################################################################
lib <- .libPaths()
# lib <- "~/R/library/unix"
contriburl <- contrib.url("https://cloud.r-project.org")
workDirectory <- getwd()
# workDirectory <- "~/utilityTDA"

#########################################################################
## Parameters for which packages to count
#########################################################################
start <- as.Date("2014-08-18")
end <- as.Date("2017-01-08")
packages <- c("TDA", "phom")

#########################################################################
## Parameters for where to save your result
#########################################################################
logDirectory <- "RStudioCRAN"
RDataFile <- "countDownloadTDA.RData"
plotFile <- "countDownloadTDA.png"




#########################################################################
## Load your packages and set working directory
#########################################################################
if (!require("data.table", lib.loc = lib)) {
  install.packages(pkgs = "data.table", lib = lib, contriburl = contriburl)
}
if (!require("ggplot2", lib.loc = lib)) {
  install.packages(pkgs = "ggplot2", lib = lib, contriburl = contriburl)
}
if (!require("plyr", lib.loc = lib)) {
  install.packages(pkgs = "plyr", lib = lib, contriburl = contriburl)
}
require("data.table", lib.loc = lib)
require("ggplot2", lib.loc = lib)
require("plyr", lib.loc = lib)

setwd(dir = workDirectory)


## ======================================================================
## Step 1: Download all log files
## ======================================================================
downloadRStudioCRAN(start = start, end = end, logDirectory = logDirectory)


## ======================================================================
## Step 2: Load single data files into one big data.table
## ======================================================================
dataset <- readRStudioCRAN(start = start, end = end,
                           logDirectory = logDirectory, packages = packages)

save(dataset, file = RDataFile)

# for later analyses: load the saved data.table
load(file = RDataFile)


## ======================================================================
## Step 3: Analyze it!
## ======================================================================
countDataset <- countPackage(dataset)
save(countDataset, file = RDataFile)

png(filename = plotFile)
countDataset[["plot"]]
dev.off()
