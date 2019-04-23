fast_leven <- function(x, y, cost_matrix=NULL)
{
    return(.Call("cleven", x, y, cost_matrix));
}

optimise_cost_matrix <- function(cost_matrix)
{
    return(.Call("optimise_cost_matrix", cost_matrix));
}

leven <- function(x, y, cost_matrix=NULL)
{
    firstRow <- c(0:length(y))
    for (i in 0:(length(x)-1))
    {
        nextRow <- c(i+1)
        for (j in 0:(length(y)-1))
        {
            nextRow[j+2] = min(nextRow[j+1] + 1, firstRow[j+2] + 1, firstRow[j+1] + ifelse(x[i+1] == y[j+1],
                                                                                           0,
                                                                                    ifelse(is.null(cost_matrix),
                                                                                           1,
                                                                                           ifelse(is.element(x[i+1], rownames(cost_matrix)) && is.element(y[j+1], rownames(cost_matrix)),
                                                                                                  cost_matrix[x[i+1], y[j+1]],
                                                                                                  1))))
        }
        firstRow <- nextRow
    }
    return (firstRow[length(y)+1])
}

lsi <- function(x, y, cost_matrix=NULL, fleven=leven)
{
    return (1-fleven(x, y, cost_matrix)/max(length(x), length(y)))
}

LSImatrix <- function(filename)
{
    return (crunch_numbers(filename)$lsi_matrix);
}

crunch_numbers <- function(filename, cost_matrix=NULL, fileEncoding="", fleven=leven)
{
    x <- scan(filename, what="", sep="\n", fileEncoding = fileEncoding);
    rows <- strsplit(x, ",[[:space:]]+")

    row_lengths <- unlist(lapply(rows, length));
    if (min(row_lengths) < 6)
    {
        index <- match(min(row_lengths), row_lengths);
        bad_row <- rows[index];
        message<-sprintf("Row %s of file %s is %s and does not contain enough columns. Please check your inputs\n", index, filename, bad_row);
        stop(message);
    }
    names = c()
    themes = c()
    theme_count = 0;
    current_theme <- 0;
    for (i in 1:length(rows))
    {
    	theme = rows[[i]][3];
        names[i] <- paste(rows[[i]][1], rows[[i]][2], theme, rows[[i]][4], sep="_");
	if (!is.element(theme, themes))
	{
            theme_count = theme_count+1
   	    themes[theme_count] <- theme;
	}
    }
    theme_matrix <- matrix(data=0, dimnames=list(c(themes), c(themes)), ncol=length(themes), nrow=length(themes))
    theme_totals <- matrix(data=0, dimnames=list(c(themes), c(themes)), ncol=length(themes), nrow=length(themes))
    theme_running_total = 0;
    set_medians <- data.frame(row.names=c(themes));
    set_medians$similarity=0;
    set_medians$witness=NA;
    set_medians$value=NA;
    
    lsi_matrix <- matrix(dimnames=list(c(names), c(names)), ncol=length(rows), nrow=length(rows))
    for (i in 1:length(rows))
    {
        #cat(sprintf("Processing row %s\n", i));
        for (j in 1:length(rows))
        {
            #cat(sprintf("   subrow %s\n", j));
            vectori = rows[[i]][5:length(rows[[i]])];
            vectorj = rows[[j]][5:length(rows[[j]])];
            Singeri = rows[[i]][1]
            Singerj = rows[[j]][1]
            Songi = rows[[i]][2]
            Songj = rows[[j]][2]
            Themei = rows[[i]][3]
            Themej = rows[[j]][3]
	    Phrasei = rows[[i]][4]
            Phrasej = rows[[j]][4]
            
            lsi_value <- lsi(vectori, vectorj, cost_matrix, fleven);
	    theme_matrix[Themei, Themej] = theme_matrix[Themei, Themej] + lsi_value;
	    theme_matrix[Themej, Themei] = theme_matrix[Themej, Themei] + lsi_value;
	    theme_totals[Themei, Themej] = theme_totals[Themei, Themej] + 1;
	    theme_totals[Themej, Themei] = theme_totals[Themej, Themei] + 1;
            
            if (Themei == Themej)
                theme_running_total = theme_running_total + lsi_value;
            
            ##cat(sprintf("LSI %s_%s_%s_s% to %s_%s_%s_s% is %f\n", Singeri, Songi, Themei, Phrasei, Singerj, Songj, Themej, Phrasej, lsi_value));
            
            lsi_matrix[i,j]  <- lsi_value            
        }
        
        if (theme_running_total > set_medians[Themei,]$similarity)
        {
            set_medians[Themei,]$similarity = theme_running_total;
            set_medians[Themei,]$witness = paste(Singeri, Songi, Themei, Phrasei, sep="_");
            set_medians[Themei,]$value = list(vectori);
        }        
        theme_running_total = 0;
    }
    theme_matrix = theme_matrix / theme_totals;
    ## Compute LSI with just medians
    median_lsi_matrix <- matrix(dimnames=list(row.names(set_medians), row.names(set_medians)), ncol=length(row.names(set_medians)), nrow=length(row.names(set_medians)))
    for (i in 1:length(row.names(set_medians)))
    {
        for (j in 1:length(row.names(set_medians)))
        {
            median_lsi_matrix[i,j] <- lsi(rapply(set_medians$value[i],c), rapply(set_medians$value[j],c), cost_matrix, fleven);
        }
    }
    
    return (list(lsi_matrix=lsi_matrix, theme_matrix=theme_matrix, set_medians=set_medians, median_lsi_matrix=median_lsi_matrix))
}

cluster <- function(lsi_matrix, method)
{
    y<-as.dist(1-lsi_matrix)
    z<-hclust(y, method=method)
    w<-as.dendrogram(z)
    plot(w)

}

bootstrap <- function(lsi_matrix, method)
{
    s<-pvclust(lsi_matrix, method.dist="euclidean", method.hclust=method, nboot=1000)
    plot(s)
    pvrect(s, alpha=0.95)
    dev.new()
    seplot(s)
    return(s)
}

#' Load a cost matrix
#'
#' This function loads a cost matrix from a file.
#' @param filename
#' @param threshhold
#' @param exponentialscale
#' @return The cost matrix
#' @export
#' @examples
#' cost_matrix_from_file("matrix.txt")
cost_matrix_from_file <- function(filename, threshhold=1, exponentialscale=1)
{
    ## First load the variables from the file
    raw_data <- read.table(filename, header=TRUE, row.names=NULL)
    ## Average across sound type. Ensure sound type is in column 1!
    averaged_data <- aggregate(. ~ Sound, raw_data, mean);

    ##    Divide each column by the maximum value in it
    ##    normal_matrix <- apply(averaged_data[,2:length(averaged_data)], 2, function(x) x / max(x))

    ## Transform from data domain to z-score domain
    normal_matrix <- scale(averaged_data[,2:length(averaged_data)], center=TRUE, scale=TRUE);

    ## Compute Euclidean distance between all the values in normal_matrix
    cost_matrix <- as.matrix(dist(normal_matrix), labels=TRUE);
    ## Glue labels back on to cost_matrix
    colnames(cost_matrix) <- rownames(cost_matrix) <- averaged_data[['Sound']];

    ## Next, divide the matrix by the maximum value in the matrix so that all of the costs are in the range 0..1

    max_entry <- max(cost_matrix)

    normalized_cost <- 1-(exp(-cost_matrix*exponentialscale))

    ## normalized_cost <- (cost_matrix/max(cost_matrix));
    ## normalized_cost <- (cost_matrix)

    ## Finally, apply the threshhold, if supplied
    normalized_cost[normalized_cost > threshhold] <- 1;
    return(normalized_cost);
}
