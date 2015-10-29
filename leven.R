leven <- function(x, y)
{
    firstRow <- c(0:length(y))
    for (i in 0:(length(x)-1))
    {
        nextRow <- c(i+1)
        for (j in 0:(length(y)-1))
            nextRow[j+2] = min(nextRow[j+1] + 1, firstRow[j+2] + 1, firstRow[j+1] + ifelse(x[i+1] == y[j+1], 0, 1))
        firstRow <- nextRow
    }
    return (firstRow[length(y)+1])
}

lsi <- function(x, y)
{
    return (1-leven(x, y)/max(length(x), length(y)))
}

LSImatrix <- function(filename)
{
    return (crunch_numbers(filename)$lsi_matrix);
}

crunch_numbers <- function(filename)
{
    x <- scan(filename, what="", sep="\n")
    rows <- strsplit(x, ",[[:space:]]+")
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
   	    themes[theme_count] <- theme;
	    theme_count = theme_count+1
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
        for (j in 1:length(rows))
        {
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
            
            lsi_value <- lsi(vectori, vectorj);
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
            median_lsi_matrix[i,j] <- lsi(rapply(set_medians$value[i],c), rapply(set_medians$value[j],c));
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
}

