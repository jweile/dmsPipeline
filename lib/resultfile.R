
new.resultfile <- function(htmlfile,append=TRUE) {

	.htmlfile <- htmlfile
	.con <- file(htmlfile,open=if(append) "a" else "w")

	#given two files from a common point of view in the filesystem,
	# this function describes the relative path to the first file
	# from the adopted point of view of the second file.
	`%relativeTo%` <- function(file,povFile) {
		fileElements <- strsplit(file,"/")[[1]]
		povElements <- strsplit(povFile,"/")[[1]]
		if (fileElements[[1]]!=povElements[[1]]) {
			common <- 0
		} else {
			common <- max(which(sapply(1:length(povElements),function(i)fileElements[[i]]==povElements[[i]])))
		}
		ups <- length(povElements)-common-1
		if (ups > 0) {
			paste(
				paste(rep("..",ups),collapse="/"),
				paste(fileElements[(common+1):length(fileElements)],collapse="/"),
				sep="/"
			)
		} else {
			paste(fileElements[(common+1):length(fileElements)],collapse="/")
		}
	}

	#Escapes non-html characters from a given string
	escape <- function(s) {
		s <- gsub("<","&lt;",s)
		s <- gsub(">","&gt;",s)
		s
	}

	#write a HTML header and start of body section
	header <- function(title) {
		writeLines(paste0("<html>
<head>
	<title>",escape(title),"</title>
	<style type=\"text/css\">
	body {
		font-family:Helvetica,sans-serif;
		padding:1cm;
		margin-bottom: 2in;
	}
	</style>
</head>
<body>"),.con)
	}

	#write the closing of the HTML body and end of file.
	closing <- function() {
		writeLines("</body>
</html>",.con)
	}

	#write a section header to the HTML file
	section <- function(title) {
		writeLines(paste0("<h1>",escape(title),"</h1>"),.con)
	}

	subsection <- function(title) {
		writeLines(paste0("<h2>",escape(title),"</h2>"),.con)
	}

	#Write a paragraph to the html file
	paragraph <- function(text) {
		writeLines(c("<p>",escape(text),"</p>"),.con)
	}

	#Draw a figure into the html file
	#draw = a function that draws the figure
	#filestub = the base of the filename (without extension)
	#w,h = width and height in inches
	#res = resolution (defaults to 100dpi)
	figure <- function(draw,filestub,w=7,h=7,res=100) {

		pdffile <- paste0(filestub,".pdf")
		pdf(pdffile,w,h)
		draw()
		invisible(dev.off())

		pngfile <- paste0(filestub,".png")
		png(pngfile,w*res,h*res,res=res)
		draw()
		invisible(dev.off())

		html.w <- if (w > 10) 10 else w
		writeLines(paste0(
			"<a href=\"",
			(pdffile %relativeTo% .htmlfile),
			"\" download><img src=\"",
			(pngfile %relativeTo% .htmlfile),
			"\" style=\"width:",html.w,"in;\"/></a>"
		),.con)
		flush(.con)
	}

	link.data <- function(filename, iconfile="doc/table_icon.png") {
		relfile <- filename %relativeTo% .htmlfile
		shortfile <- sub(".*/","",filename)
		writeLines(paste0(
"<p>
	<a href=\"",relfile,"\" style=\"vertical-align:top;\" download>
		<img src=\"",(iconfile %relativeTo% .htmlfile),"\" style=\"width:18px\"/>
		<span>&nbsp;</span>",
		shortfile,
	"</a>
</p>"
		),.con)
		flush(.con)
	}

	#close the file
	shutdown <- function() {
		close(.con)
	}

	list(
		header=header,
		closing=closing,
		section=section,
		subsection=subsection,
		paragraph=paragraph,
		figure=figure,
		link.data=link.data,
		shutdown=shutdown
	)
}

test.htmlfile <- function() {
	dir.create("htmltest")
	html <- new.resultfile("htmltest/test.html",append=FALSE)
	html$header("Test")
	html$section("Lorem ipsum dolor sit amet")
	html$paragraph("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Cras vitae porttitor eros, ut condimentum lacus. Aliquam faucibus augue vel dolor sollicitudin ornare. Fusce id lectus vitae felis fermentum consequat.")
	html$subsection("Lorem")
	html$figure(function(){
		plot(rnorm(10000,0,1),rnorm(10000,0,1))
	},"htmltest/testfigure",7,5)
	html$link.data("foobar.txt")
	html$closing()
	html$shutdown()
}
# test.htmlfile()
