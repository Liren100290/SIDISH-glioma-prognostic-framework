# R-only metadata QA for the final scientific exports. Visual inspection is a
# separate, explicitly recorded check and is not inferred from these assertions.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
.libPaths(c(file.path(R.home(),'library'),.libPaths()))
av<-commandArgs(trailingOnly=FALSE)
here<-dirname(normalizePath(sub('^--file=','',av[grepl('^--file=',av)]),winslash='/',mustWork=TRUE))
library(jsonlite)
for(p in c('ggplot2','patchwork','svglite','ragg','survival'))loadNamespace(p)

tiff_info<-function(path){
 f<-file(path,'rb');on.exit(close(f))
 bytes<-readBin(f,'raw',2)
 endian<-if(rawToChar(bytes)=='II')'little' else 'big'
 stopifnot(readBin(f,'integer',1,size=2,signed=FALSE,endian=endian)==42)
 off<-readBin(f,'integer',1,size=4,endian=endian);seek(f,off)
 n<-readBin(f,'integer',1,size=2,signed=FALSE,endian=endian)
 ans<-list()
 for(i in seq_len(n)){
  tag<-readBin(f,'integer',1,size=2,signed=FALSE,endian=endian)
  type<-readBin(f,'integer',1,size=2,signed=FALSE,endian=endian)
  count<-readBin(f,'integer',1,size=4,endian=endian)
  raw<-readBin(f,'raw',4)
  if(!tag%in%c(256L,257L,282L,283L,296L))next
  r<-rawConnection(raw);value<-readBin(r,'integer',1,size=if(type==3)2 else 4,signed=if(type==3)FALSE else TRUE,endian=endian);close(r)
  if(type==5){
   pos<-seek(f);seek(f,value)
   q<-readBin(f,'integer',2,size=4,endian=endian)
   value<-q[1]/q[2];seek(f,pos)
  }
  ans[[as.character(tag)]]<-value
 }
 ans
}
spec<-data.frame(name=c('Figure6','FigureS4'),width_mm=c(183,183),height_mm=c(169,235))
results<-list()
for(i in seq_len(nrow(spec))){
 stem<-file.path(here,paste0(spec$name[i],'_revision_v05'))
 paths<-paste0(stem,c('.svg','.pdf','.tiff','.png'))
 stopifnot(all(file.exists(paths)),all(file.info(paths)$size>1000))
 svg<-paste(readLines(paste0(stem,'.svg'),warn=FALSE),collapse='\n')
 matches<-gregexpr('<text(?:\\s|>)',svg,perl=TRUE)[[1]]
 live<-if(matches[1]<0)0L else length(matches)
 stopifnot(live>10L,!grepl('<image',svg,fixed=TRUE))
 t<-tiff_info(paste0(stem,'.tiff'))
 stopifnot(t[['282']]==600,t[['283']]==600,t[['296']]==2)
 stopifnot(abs(t[['256']]-spec$width_mm[i]/25.4*600)<2,
           abs(t[['257']]-spec$height_mm[i]/25.4*600)<2)
 results[[spec$name[i]]]<-list(width_mm=spec$width_mm[i],height_mm=spec$height_mm[i],
       svg_live_text_elements=live,svg_embedded_raster=FALSE,
       tiff_width_pixels=t[['256']],tiff_height_pixels=t[['257']],tiff_dpi=c(t[['282']],t[['283']]),
       all_four_export_formats_present=TRUE)
}
write_json(list(metadata_checks_passed=TRUE,figures=results,
 visual_review=list(Figure6='Viewed final three-panel PNG: all 16 compound rows aligned; all titles, numeric labels, axes and caption readable.',
 FigureS4='Viewed final forest and three-scatter PNG: all 16 tests, point annotations and full caption readable.'),
 scope='Export/visual QA only; no new statistical or experimental validation.'),
 file.path(here,'Figure_export_QA_v05.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('Both v05 figures: live SVG text, expected dimensions and 600-dpi TIFF verified.\n')
capture.output(sessionInfo(),file=file.path(here,'Figure_R_sessionInfo_v05.txt'))
