###########################
# Colorize pdb structures #
###########################

new.pymol.colorizer <- function(outfile,chain="A") {

	.con <- file(outfile,open="w")
	.chain <- chain

	val2col <- function(s) {
		if (is.na(s)) "gray60"
		else if (s <= 0) "lof1"
		else if (s <= 0.2) "lof3"
		else if (s <= 0.4) "lof4"
		else if (s <= 0.6) "lof5"
		else if (s <= .8) "lof6"
		else if (s <= 1.2) "wt"
		else if (s <= 1.4) "gof1"
		else if (s <= 1.6) "gof2"
		else if (s <= 1.8) "gof3"
		else if (s <= 2) "gof4"
		else "gof6"
	}

	define.colors <- function() {
		writeLines(
"set_color lof1 = [0.23 , 0.37 , 0.80]
set_color lof2 = [0.35 , 0.47 , 0.84]
set_color lof3 =  [0.48 , 0.58 , 0.87]
set_color lof4 =  [0.61 , 0.69 , 0.90]
set_color lof5 =  [0.74 , 0.79 , 0.93]
set_color lof6 =  [0.87 , 0.89 , 0.96]
set_color wt = [1.00 , 1.00 , 1.00]
set_color gof1 =  [0.96 , 0.85 , 0.85]
set_color gof2 =  [0.93 , 0.71 , 0.71]
set_color gof3 =  [0.90 , 0.57 , 0.57]
set_color gof4 =  [0.87 , 0.43 , 0.43]
set_color gof5 =  [0.84 , 0.29 , 0.29]
set_color gof6 =  [0.80 , 0.15 , 0.15]"
		,.con)
	}

	#values should be a matrix containing two columns: the residue index and the fitness value
	colorize <- function(values) {
		writeLines(apply(values,1,function(v) {
			paste("color",val2col(v[["fitness"]]),", chain",.chain,"& resi",v[["index"]])
		}),.con)
	}

	.close <- function() {
		close(.con)
	}

	list(
		define.colors=define.colors,
		colorize=colorize,
		close=.close
	)
}

test.colorizer <- function() {
	pycol <- new.pymol.colorizer("test.txt")
	pycol$define.colors()
	pycol$colorize(cbind(index=1:10,fitness=runif(10,0,1)))
	pycol$close()
}
# test.colorizer()
