
new.amasLite <- function(amasFile="res/amas_props.csv") {

	.props <- read.csv(amasFile)
	colnames(.props)[[ncol(.props)]] <- "-"

	toChars <- function(s) sapply(1:nchar(s),function(i)substr(s,i,i))

	readFASTA <- function(fastaFile) {

		con <- file(fastaFile,open="r")

		currSeq <- character(0)

		out <- list()
		outnames <- character(0)

		while (length(line <- readLines(con,1))>0) {
			if (substr(line,1,1)==">") {
				#store name
				outnames[length(outnames)+1] <- line
				#write old sequence
				if (length(currSeq) > 0) {
					out[[length(out)+1]] <- currSeq
				}
				#prep for next sequence
				currSeq <- character(0)
			} else {
				#append to current sequence
				currSeq <- c(currSeq,toChars(line))
			}
		}
		#write the last sequence
		out[[length(out)+1]] <- currSeq
		close(con)

		out <- do.call(rbind,out)
		rownames(out) <- outnames
		out
	}

	run <- function(alignmentFile) {

		alignment <- readFASTA(alignmentFile)

		scores <- apply(alignment,2,function(x) {
			aas <- unique(x)
			sum(apply(.props[,aas,drop=FALSE],1,function(p) {
				all(p) | !any(p)
			}))
		})

		human.pos <- which(alignment[1,] != "-")
		scores[human.pos]
	}

	list(run=run)
}
