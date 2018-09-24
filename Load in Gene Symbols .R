setwd('C:/Users/JLinehan/Documents/RCodes/Gene_set_enrichment_analysis')

#setwd("E:\Gene Set enrichment Analysis") 

hallmark_set = scan('hallmark_gene_symbols.txt', what = 'list')

# hallmark_set_ids = scan('h.all.v6.2.entrez_copy.txt', what = 'list',sep = ' ')
# this was a test of loading in the gene ids, its a little bit more difficult to do, and would be easier
# to deal with by querying the gene ids from ensemble after having broken up the hallmark_set into subsets 

char_num_index <- lapply(hallmark_set, nchar)
# Gives me the number of characters in each string 
count = 1 
URL_location = list()
j = 1 
for (j in 1:length(hallmark_set)){
	
	my_name <- strsplit(hallmark_set[j],":")
	
	this_length <- length(unlist(my_name))
	
	if (this_length > 1){
	
	URL_location[count] <- j
	count = count +1 
	
	}
} 

n = length(unlist(URL_location))
############################################################################
## Now we go through and break up the hallmark collection into gene sets 

organized_hallmark_set <- matrix(list(),1,n)
 
for (j in 1:length(URL_location)){

begin <- unlist(URL_location[j]) +1 

if (j<50){
end <- unlist(URL_location[j+1])- 2 
organized_hallmark_set[[1,j]] <- hallmark_set[begin:end]
}

if (j == 50){
	
	organized_hallmark_set[[1,j]] <- hallmark_set[begin:length(hallmark_set)]
}  

} 

## So now I have the gene sets broken up into seperate groups, I just need to remember to add their names before we go ahead and shove them into a text file after conversion. 
name_location <- unlist(URL_location) - 1 
gene_set_names <- hallmark_set[name_location]

# Great so now we have them, moving on 
# I can now go ahead and query ensembl for gene ids, and look up the zfish symbols and ortho/para logs 

library(biomaRt)
#source("https://bioconductor.org/biocLite.R")
#biocLite("mygene")
library(mygene)

n = length(organized_hallmark_set)

my_mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Take mammilian gene symbols and get back there gene ids 
#guery_results <- queryMany(gene_set,scopes="symbol", fields = "ensembl.gene",species = "human")

#gene_set_ids <- data.matrix(query_results[2])

mart = useMart("ensembl",dataset = "drerio_gene_ensembl")

#drerio_gene_symbols <- getGenes(gene_set_ids,mart=mart) 

# Preallocate for gene symbols 

drerio_gene_symbol <- matrix(list(),1,n)

for (j in 1:n ){ 
	
	gene_id_carrier = queryMany(organized_hallmark_set[[1,j]], scopes = "symbol", species = "zebrafish", return.as = "DataFrame")
	
	index <- (colnames(gene_id_carrier) == "_id")
	
	id_location <- match("TRUE", index)  

	gene_set_ids <- data.matrix(gene_id_carrier[id_location]) 

	nona_gene_set_ids <- gene_set_ids[!is.na(gene_set_ids)] 

	drerio_gene_symbols <- getGenes(nona_gene_set_ids,mart=mart)

	#Now I need to once again find the _id column 

	index <- (colnames(drerio_gene_symbols) == "symbol")

	id_location <- match("TRUE", index)  

	na_drerio <- data.matrix(drerio_gene_symbols[id_location])

	drerio_gene_symbol[[1,j]] <- c(gene_set_names[j],"na", na_drerio) 


	}



#write.csv(drerio_gene_symbol,file = "Danio_Rerio_hallmark_gene_sets.xls",col.names = NA, row.names = TRUE) 
z = list() 
for (j in 1:n) { 

lapply(drerio_gene_symbol[[1,j]], write, "test.txt", append=TRUE, ncolumns=n)

}

lapply(  writeLines(unlist(lapply(drerio_gene_symbol, paste, collapse=" ")))
,write,"test2.txt",append=TRUE,ncolumns=n) 
write.csv(paste(unlist(z),collapse = ' '), file = "Danio_Rerio_hallmark_gene_sets.txt")


#write.xlsx(paste(unlist(drerio_gene_symbol[[j]]),collapse = ' '),file = "Danio_Rerio_hallmark_gene_sets.xlsx",append = TRUE) 


for (j in 1:n){

#write.xlsx(drerio_gene_symbol[[j]], file = "Danio_Rerio_hallmark_gene_sets.xls", append=TRUE) 

#write(paste(unlist(drerio_gene_symbol[[j]]), collapse = ' '), file = "Danio_Rerio_hallmark_gene_sets.txt") 

#write(paste(unlist(drerio_gene_symbol[[j]]), collapse = ' '), file = "Danio_Rerio_hallmark_gene_sets.xls", append = TRUE, sep = '\t')

} 



writeLines(unlist(lapply(drerio_gene_symbol, paste, collapse=" ")))


sink("C:/Users/JLinehan/Documents/RCodes/Gene_set_enrichment_analysis/myTest.dat")
writeLines(unlist(lapply(drerio_gene_symbol, paste, collapse=" ")))
sink()


write.xlsx(unlist(lapply(drerio_gene_symbol, paste, collapse=" ")),file='my_test.xls')







 
