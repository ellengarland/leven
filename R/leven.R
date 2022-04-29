#' Fast Levenshtein calculation
#'
#' This function is optimised to compute the Levenshtein distance as quickly as possible
#' @param x The first vector
#' @param y The second vector
#' @param cost_map An optional cost map specifying the cost of specific distributions. Not to be confused with a cost matrix. 
#' @return The distance between x and y
#' @seealso \link{cost_matrix_from_file}
#' @seealso \link{optimise_cost_matrix}
#' @export
#' @examples leven(c("k", "i", "t", "t", "e", "n"), c("s", "i", "t", "t", "i", "n", "g"))
fast_leven <- function(x, y, cost_map=NULL)
{
    return(.Call("cleven", x, y, cost_map));
}

#' Generate a cost map
#'
#' Produce a cost map from a cost matrix. This function is required to produce a fast-lookup cost map for supply to fast_leven
#' @param cost_matrix The cost matrix to convert
#' @return The generated cost map
#' @seealso \link{cost_matrix_from_file}
#' @export
#' @examples optimise_cost_matrix(cost_matrix_from_file("matrix.txt"))
optimise_cost_matrix <- function(cost_matrix)
{
    return(.Call("optimise_cost_matrix", cost_matrix));
}

#' Slow Levenshtein calculation
#'
#' This function is designed to compute the Levenshtein distance in a way which is simple to understand and runs in R without any native support
#' @param x The first vector
#' @param y The second vector
#' @param cost_matrix An optional cost matrix to specify the cost of specific substitutions.
#' @return The distance between x and y
#' @seealso \link{cost_matrix_from_file}
#' @export
#' @examples leven(c("k", "i", "t", "t", "e", "n"), c("s", "i", "t", "t", "i", "n", "g"))
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

#' Levenshtein Similary Index between two vectors
#'
#' This function computes the Levenshtein Similary Index between two vectors
#' @param x The first vector
#' @param y The second vector
#' @param cost_matrix An optional cost matrix to specify the cost of specific substitutions.
#' @param fleven The function to use to compute the Levenshtein distance. Defaults to leven. Use fast_leven instead to use the native (fast) code
#' @return The LSI of x and y
#' @seealso \link{cost_matrix_from_file}
#' @seealso \link{fast_leven}
#' @seealso \link{leven}
#' @export
#' @examples lsi(c("k", "i", "t", "t", "e", "n"), c("s", "i", "t", "t", "i", "n", "g"))
lsi <- function(x, y, cost_matrix=NULL, fleven=leven)
{
    return (1-fleven(x, y, cost_matrix)/max(length(x), length(y)))
}

#' Levenshtein Similarity Index for each pair of strings in the input file
#' 
#' This function computes the Levenshtein Similarity Index for each possible pair of strings in the input file
#' @param filename The file containing the input data. The input must be comma-separated values, where the columns are:
#' \itemize{ \item Location 
#' \item Year 
#' \item Song type 
#' \item Theme 
#' \item Unit 1 
#' \item Unit 2 
#' \item Unit 3 
#' \item ..... 
#' \item Unit n
#' }
#' @return A symmetric matrix where the cell in row i, col j is the Levenshtein Similarity Index between line i and line j from the file filename
#' @export
lsi_matrix <- function(filename)
{
    return (compare_songs(filename)$lsi_matrix);
}

#' Compare Songs
#'
#' This function compares all song strings in the given file
#' @param filename The file containing the input data. See the description for \link{lsi_matrix} for the format of the input data
#' @param cost_matrix An optional cost matrix to use when computing the Levenshtein distance
#' @param fileEncoding character string: if non-empty declares the encoding used on the file so the character data can be read correctly. See the 'Encoding' section of the help for file
#' @param fleven The Levenshtein-calculating function to use. Defaults to leven (slow). You can supply fleven as the argument to use the fast (native) code
#' @return A list with values: \itemize{
#' \item lsi_matrix: the generated LSI matrix - see \link{lsi_matrix}
#' \item theme_matrix: A matrix where cell i,j is the average LSI between theme i and theme j
#' \item set_medians: A list of median distances for each theme. Each entry contains these hedgehogs: \itemize{
#'     \item $similarity (the distance)
#'     \item $witness (the identifier of the singer which is the median)
#'     \item $value (the string which is the median)}
#' \item median_lsi_matrix: The LSI matrix between medians for each theme. Cell i,j in this matrix is the LSI of the median of theme i and the median of theme j}
#' @export
#' @examples compare_songs("data.txt", fileEncoding="utf8", fleven=fleven)
compare_songs <- function(filename, cost_matrix=NULL, fileEncoding="", fleven=leven)
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

#' Heirarchical clustering plot of the specified LSI matrix
#'
#' This function produces a plot of the heirarchical clustering of the given LSI matrix using the named method for clustering
#' @param lsi_matrix The matrix to cluster. Generated from lsi_matrix()
#' @param method The name of the method to use. See the documentation for hclust() for options
#' @export
#' @examples cluster(lsi_matrix("data.txt"), "average")
cluster <- function(lsi_matrix, method)
{
    y<-as.dist(1-lsi_matrix)
    z<-hclust(y, method=method)
    w<-as.dendrogram(z)
    plot(w)

}

#' Generates the bootstrapped heirarchical clustering of the given LSI matrix using named method for clustering
#' @param lsi_matrix The matrix to cluster. Generated from lsi_matrix()
#' @param method The name of the method to use. See the documentation for hclust() for options
#' @param nboot The number of bootstrap replications. The default is 1000.
#' @param distance_method the distance measure to be used. This should be a character string, or a function which returns a dist object. A character string should be (an abbreviation of) one of "correlation", "uncentered", "abscor" or those which are allowed for method argument in dist function. The default is "correlation". See details section in this help and method argument in dist
#' @return ???
#' @export
#' @examples bootstrap(lsi_matrix("data.txt"), "average")
#' @examples bootstrap(lsi_matrix("data.txt"), "average", nboot=500, distance_method="uncentered")
bootstrap <- function(lsi_matrix, method, nboot=1000, distance_method="euclidean")
{
    s<-pvclust(lsi_matrix, method.dist=distance_method, method.hclust=method, nboot=nboot)
    plot(s)
    pvrect(s, alpha=0.95)
    dev.new()
    seplot(s)
    return(s)
}

#' Load a cost matrix
#'
#' This function loads a cost matrix from a file
#' @param filename ?
#' @param threshhold ?
#' @param exponentialscale ?
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

    ## Scale using an exponential function. The scale can be provided here in exponentialscale to spread or contract the curve. Default is 1
    normalized_cost <- 1-(exp(-cost_matrix*exponentialscale))

    ## Alternatively, we could divide the matrix by the maximum value in the matrix so that all of the costs are in the range 0..1
    ## max_entry <- max(cost_matrix)
    ## normalized_cost <- (cost_matrix/max(cost_matrix));

    ## Or we could not scale at all
    ## normalized_cost <- (cost_matrix)

    ## Finally, apply the threshhold, if supplied
    normalized_cost[normalized_cost > threshhold] <- 1;
    return(normalized_cost);
}
