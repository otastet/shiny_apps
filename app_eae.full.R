library(shinydashboard)
library(shiny)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(ggiraph)
library(cowplot)
library(shinymanager)
library(Seurat)
library(ggrepel)
library(viridis)
library(stringr)
library(shinyjs)
library(shinyWidgets)
library(shinyBS)
library(bslib)
library(RColorBrewer)

add_score<-function(obj,geneset,rows){
  obj =  AddModuleScore(obj,geneset)
  return(as.numeric(obj@meta.data[rows,]$Cluster1))
}

## my css
CSS <- function(colors){
  template <- "
.checkboxGroupButtons div.btn-group:nth-child(%s) button {
  background: %s !important;
  color: black !important;
  padding: 5px;
  margin-bottom: 8px
}"
  paste0(
    apply(cbind(seq_along(colors), colors), 1, function(vc){
      sprintf(template, vc[1], vc[2])
    }),
    collapse = "\n"
  )
}
cols <- c("darkorange","dodgerblue")
mycss <- CSS(cols)


radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){

  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

datasets = c('obj_shiny.total.dec2021.rds','obj_shiny.cd31.dec2021.rds')
names(datasets) = c('total_cells','cd31_selection')

sample.colors=c('dodgerblue3','dodgerblue1','deepskyblue','chocolate','darkorange','goldenrod2')
names(sample.colors)=c('CTL1','CTL2','CTL3','EAE1','EAE2','EAE3')


ui <- dashboardPage(
  dashboardHeader(title = "EAE interactions"),
    dashboardSidebar(						
		fileInput("upload", NULL, accept = c(".csv", ".tsv"))
),

  dashboardBody(
      box(height=550,
	  	  title=textOutput("umap_title"),
		  solidHeader=T,
			fluidRow(
				column(8,dropdownButton(
				width=12,
				circle = FALSE, 
				status = "danger",
    			icon = icon("gear"),
				tooltip = tooltipOptions(title = "Click to see options!"),
  		 		tags$h3("List of Input"),
				fluidRow(
					column(6,selectInput('ds','dataset',choices=datasets,selected=datasets[2])),
					column(6,conditionalPanel("input.ds == 'obj_shiny.cd31.dec2021.rds'",
				  	materialSwitch(
				   		inputId = "int",
				   		label = "Integration",
				   		value = FALSE, 
				    	status = "danger"
					)))),
            	radioGroupButtons(
   							inputId = "feat",
   							label = "Choose feature:", 
    						choices = c(`<i class='fa fa-brain'></i>` = "celltype",
										`<i class='fa fa-frog'></i>` = "condition", 
        								`<i class='fa fa-dna'></i>` = "gene",
										`<i class='fa fa-clock'></i>` = "time",
										`<i class='fa fa-cloud-download-alt'></i>` = "module"),
   							justified = TRUE
						),
						conditionalPanel(
						condition = "input.feat == 'gene' || input.feat == 'time'",
                     	selectizeInput("gene","Gene to color",c()))
					#	

						
				
						
                     	
        )),
		  tags$head(tags$style(HTML(mycss))),
		  column(2,checkboxGroupButtons("show",label='Show',choices=c('EAE','CTL'),selected=c('EAE','CTL'))),
		  column(2,materialSwitch('label','Labels',value=FALSE,status='success'))),		  
		  shinycssloaders::withSpinner(type=1,size=1,girafeOutput('umap',height=400))
	  ),
      box(height=550,
	  	  title=textOutput("box2_title"),
		  conditionalPanel(
			  condition = "input.feat == 'gene'",
              switchInput(inputId = "split",label='Split', value = FALSE)),
		  shinycssloaders::withSpinner(type=1,size=1,girafeOutput('box2',height=450))
	  ),
	tags$script(HTML("$('.box').eq(0).css('border', '5px solid #3DA0D1');")),
	tags$script(HTML("$('.box').eq(1).css('border', '5px solid #3DA0D1');")),
 	tags$style(HTML(".checkboxGroupButtons div.btn-group:nth-child(0) button {
  background-color: #3DA0D1;}"))
	)
)

server <- function(session, input, output) {

	rv = reactiveValues(obj=NULL,selected=NULL,score=NULL,geneset=NULL)

	observeEvent(input$ds,{

		if(grepl('cd31',input$ds)){

			updateRadioGroupButtons(session,'feat',choices = 
								c(`<i class='fa fa-brain'></i>` = "celltype",
								  `<i class='fa fa-frog'></i>` = "condition", 
        						  `<i class='fa fa-dna'></i>` = "gene",
								  `<i class='fa fa-cloud-download-alt'></i>` = "module",
								  `<i class='fa fa-clock'></i>` = "time"),
								  selected='celltype')
		}else{

			updateRadioGroupButtons(session,'feat',choices = 
								c(`<i class='fa fa-brain'></i>` = "celltype",
								  `<i class='fa fa-frog'></i>` = "condition", 
        						  `<i class='fa fa-dna'></i>` = "gene",
								  `<i class='fa fa-cloud-download-alt'></i>` = "module"),
								  selected='celltype')
		}
		rv$obj = readRDS(input$ds)
		rv$selected = NULL
		updateMaterialSwitch(session,'int',value=FALSE)

	})

	observeEvent(input$feat,ignoreInit=TRUE,{

		    rv$selected=NULL
			if(input$feat=='time'){
				updateSelectizeInput(session,'gene',choices=c('pseudo1','pseudo2'))

			}else if(input$feat=='gene'){
				updateSelectizeInput(session,'gene',choices=rownames(rv$obj$obj),selected='Ackr1')
			}else if(input$feat=='module'){
				score = add_score(rv$obj$obj,rv$geneset,rownames(rv$obj$obj@meta.data))
				rv$score=score
			}
	})



	observeEvent(input$int,	{
		output$umap<-renderGirafe({

		dat = data.frame(rv$obj$obj@meta.data)

		if(input$int){
			dat$UMAP_1=rv$obj$obj@meta.data$int.UMAP_1
			dat$UMAP_2=rv$obj$obj@meta.data$int.UMAP_2
		}else{
			dat$UMAP_1=rv$obj$obj@reductions$umap@cell.embeddings[,1]
        	dat$UMAP_2=rv$obj$obj@reductions$umap@cell.embeddings[,2]
		}
		dat$panglao=factor(dat$panglao,levels=names(rv$obj$cols))
		if(is.null(input$show)||length(input$show)==2){
				dat=dat
		}else if(length(input$show)<2){
			if(input$show=='EAE'){
				dat=subset(dat,condition=='E')
			}else if(input$show=='CTL'){
				dat=subset(dat,condition=='C')
			}
		}
		if(grepl('cd31',input$ds)&!input$int&length(input$show)==2){
			ag = aggregate(dat[,c('UMAP_1','UMAP_2')],by=list(panglao=dat$panglao,condition=dat$condition),mean)
		}else{
			ag = aggregate(dat[,c('UMAP_1','UMAP_2')],by=list(panglao=dat$panglao),mean)
			ag$condition='All'
		}
		
		if(input$feat=='celltype'){

			g = ggplot(dat,aes(UMAP_1,UMAP_2,col=panglao))+geom_point(size=.2)+theme_bw()+guides(col=F)+scale_color_manual(values=rv$obj$cols)
				if(grepl('cd31',input$ds)&!input$int&length(input$show)==2&input$label){
					g=g+geom_line(data=ag,aes(group=panglao),col='black',linetype='dashed',size=1)
				}
			t=data.frame(table(dat$panglao))
			t$Var1=factor(t$Var1,levels=t[order(-t$Freq),'Var1'])
			g2=ggplot(t,aes(Var1,Freq,fill=Var1))+geom_col()+theme_bw()+theme(axis.text.x=element_blank())+xlab('')+ylab('Nb Cell')+scale_fill_manual(values=rv$obj$cols)
			

		}else if(input$feat=='condition'){
			print(head(dat))
			g = ggplot(dat,aes(UMAP_1,UMAP_2))+geom_point(aes(col=orig.ident),size=.2)+theme_bw()+guides(col=F)+scale_color_manual(values=sample.colors)+guides(col=F)
			
			t=data.frame(table(dat$orig.ident))
			g2=ggplot(t,aes(Var1,Freq,fill=Var1))+geom_col()+theme_bw()+theme(axis.text.x=element_blank())+scale_fill_manual(values=sample.colors)+xlab('')+ylab('Nb Cell')
			

		}else if(input$feat=='gene'){

			dat$gene = rv$obj$obj[['SCT']]@data[input$gene,rownames(dat)]
			dat=dat[order(dat$gene),]
			g = ggplot(dat,aes(UMAP_1,UMAP_2))+geom_point(aes(col=gene),size=.2)+theme_bw()+scale_color_gradientn(colors=inferno(10))+guides(col=F)
			if(input$split){
				g2=ggplot(dat,aes(gene,fill=orig.ident))+geom_density()+theme_bw()+scale_fill_manual(values=sample.colors)
			}else{
				g2=ggplot(dat,aes(gene,fill=panglao))+geom_density()+theme_bw()+scale_fill_manual(values=rv$obj$cols)
			}
			
		}else if(input$feat=='module'){
			dat$score = rv$score
			dat=dat[order(dat$score),]
			g=ggplot(dat)+geom_point(aes(UMAP_1,UMAP_2,col=score))+scale_color_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()+guides(col=F)
		}else{
			g=ggplot(dat)+geom_point(aes(UMAP_1,UMAP_2,col=dat[,input$gene]))+scale_color_gradientn(colors=viridis(10))+theme_bw()+labs(col='Lineage1')+guides(col=F)
			if(input$int){
				g=g+geom_path(data=subset(rv$obj$curve,Lineage==paste0('L',gsub('pseudo','',input$gene))),arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),size=2,col='red',aes(UMAP_1,UMAP_2))
			}
			if(input$split){
				g2=ggplot(dat,aes(gene,fill=orig.ident))+geom_density()+theme_bw()+scale_color_manual(values=sample.colors)
			}else{
				g2=ggplot(dat,aes(gene,fill=panglao))+geom_density()+theme_bw()+scale_color_manual(values=rv$obj$cols)
			}
		}
		if(input$label){
			g=g+geom_label_interactive(data=ag,aes(label=panglao,tooltip=panglao,data_id=panglao),col='black',alpha=.7)
		}
		girafe(ggobj=g,height_svg=4,options = list(opts_selection(type = "single")))
		})
	})



	data <- reactive({
    	req(input$upload)
    
    	ext <- tools::file_ext(input$upload$name)
    	switch(ext,
      		csv = vroom::vroom(input$upload$datapath, delim = ","),
      		tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
      		validate("Invalid file; Please upload a .csv or .tsv file")
    	)
 	 })
  

	observeEvent(input$umap_selected, ignoreNULL = FALSE,{
		rv$selected=input$umap_selected
	})


	observeEvent(rv$selected,ignoreNULL = FALSE,{

		output$box2<-renderGirafe({
		dat = data.frame(rv$obj$obj@meta.data)

		if(is.null(rv$selected)){
			if(is.null(input$show)||length(input$show)==2){
				dat=dat
			}else if(length(input$show)<2){
				if(input$show=='EAE'){
					dat=subset(dat,condition=='E')
				}else if(input$show=='CTL'){
					dat=subset(dat,condition=='C')
				}
			}

			if(input$feat%in%c('celltype')){
				t = table(dat$panglao)
				t = t/sum(t)
				t = data.frame(t)
				t$panglao=factor(t$Var1,levels=names(rv$obj$cols))
				g=ggplot(t,aes('',Freq,fill=panglao))+geom_col()+coord_polar('y')+theme_bw()+scale_fill_manual(values=rv$obj$cols)

			}else if(input$feat=='condition'){
				t = table(dat[,c('panglao','orig.ident')])
				for(i in 1:ncol(t)){t[,i]=t[,i]/sum(t[,i])}
				for(i in 1:nrow(t)){t[i,]=t[i,]/sum(t[i,])}
				t = data.frame(t)
				t$panglao=factor(t$panglao,levels=names(rv$obj$cols))
				g=ggplot(t,aes(panglao,Freq,fill=orig.ident))+geom_col()+theme_bw()+coord_flip()+scale_fill_manual(values=sample.colors)

			}else if(input$feat %in% c('gene','time','module')){

				if(input$feat=='gene'){
					dat$gene = rv$obj$obj[['SCT']]@data[input$gene,rownames(dat)]
				}else if(input$feat=='time'){
					dat$gene = rv$obj$obj@meta.data[,input$gene]
				}else{
					if(is.null(rv$score)){
						geneset <- scan(input$file1$datapath,what=character())
						score = add_score(rv$obj$obj,rv$geneset,rownames(rv$obj$obj@meta.data))
						rv$score=score
					}
					dat$gene = rv$score
					
				}
				print(head(dat))
				mean = aggregate(dat$gene,by=list(panglao=dat$panglao),mean)
				dat$panglao = factor(dat$panglao,levels=mean[order(-mean$x),'panglao'])
				if(!input$split){
					g=ggplot(dat,aes(panglao,gene))+geom_violin(aes(fill=panglao))+theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))+geom_point(position=position_jitter(width=.1),size=.1,aes(col=panglao))+scale_fill_manual(values=rv$obj$cols)+scale_color_manual(values=rv$obj$cols)+guides(fill=F,col=F)+RotatedAxis()
				}else{
					g=ggplot(dat,aes(panglao,gene,fill=condition))+geom_violin()+theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))+scale_fill_manual(values=rev(c('darkorange','dodgerblue')))+geom_point(position=position_jitterdodge(jitter.width=.1),size=.1,aes(col=condition))+scale_color_manual(values=rev(c('darkorange','dodgerblue')))+guides(fill=F,col=F)
				}
				if(input$ds=='obj_shiny.cd31.dec2021.rds'&input$feat!='time'){
					dat$gene = rv$obj$obj[['SCT']]@data[input$gene,rownames(dat)]
					dat=dat[order(dat$pseudo1),]
					dat=dat[which(!is.na(dat$pseudo1)),]
					dat$rank=1:nrow(dat)
					if(input$split){
						dat[which(dat$condition=='C'),'gene']=-1*dat[which(dat$condition=='C'),'gene']
						
						g2=ggplot(dat,aes(rank,gene,col=condition))+geom_col()+theme_bw()+scale_color_manual(values=c('dodgerblue','darkorange'))+guides(col=F)
					}else{
						g2=ggplot(dat,aes(rank,gene,col=panglao))+geom_col()+theme_bw()+scale_color_manual(values=rv$obj$cols)+guides(col=F)
					}
					g=plot_grid(g,g2,ncol=1)
				}
			}
		}else{
			degs.sub =  subset(rv$obj$degs,celltype==input$umap_selected&p_val_adj<0.01&!grepl('Rik$',gene))
			markers.sub =  subset(rv$obj$markers,cluster==input$umap_selected&p_val_adj<0.01&!grepl('Rik$',gene))

			markers.sub$pct_ratio=log(markers.sub$pct.1/markers.sub$pct.2)
			degs.sub$pct_ratio=log(degs.sub$pct.1/degs.sub$pct.2)
			
			degs.sub$gene=gsub('\\.|-','',degs.sub$gene)
			markers.sub$gene=gsub('\\.|-','',markers.sub$gene)

			degs.sub$gene_id  = 1:nrow(degs.sub)			
			markers.sub$gene_id  = 1:nrow(markers.sub)

			plots=list()

			if(nrow(markers.sub)>0){
				plots[['markers']] = ggplot(head(markers.sub,max(c(nrow(markers.sub),300))),aes(size=pct_ratio,input$umap_selected,avg_log2FC,col=avg_log2FC))+geom_point_interactive(aes(tooltip=gene,data_id=gene),position=position_jitter(width=.1))+scale_color_gradient2(low='grey90',high=rv$obj$cols[input$umap_selected])+coord_flip()+theme_bw()+guides(col=F,size=F)+xlab('')+ggtitle(paste0(input$umap_selected,' Markers'))+theme(axis.text.y=element_blank())+ scale_size(range = c(0.1, 3))+ylim(c(0,max(subset(markers.sub,p_val_adj<0.01)$avg_log2FC)*1.2))
			}
			if(nrow(degs.sub)>0){
				degs.sub$p_val_adj[degs.sub$p_val_adj==0]=min(degs.sub$p_val_adj[degs.sub$p_val_adj!=0])
				if(nrow(degs.sub)<500){
					n=nrow(degs.sub)
					label = as.character(n)
				}else{
					n=500
					label='top 500'
				}
				
				plots[['degs']] = ggplot(head(degs.sub,n),aes(size=pct_ratio,-log10(p_val_adj),x=avg_log2FC,col=avg_log2FC))+geom_point_interactive(aes(tooltip=gene,data_id=gene))+scale_color_gradient2(low='dodgerblue',high='darkorange',mid='grey70')+theme_bw()+guides(col=F,size=F)+ylab('')+ggtitle(paste0('DEGs in ',input$umap_selected,'(',label,')'))+ scale_size(range = c(0.1, 3))+ylim(c(-1,-log10(min(subset(degs.sub,p_val_adj<0.01)$p_val_adj))*1.2))+xlim(c(min(subset(degs.sub,p_val_adj<0.01)$avg_log2FC)*1.2,max(subset(degs.sub,p_val_adj<0.01)$avg_log2FC)*1.2))
			}
			g = plot_grid(plotlist=plots,ncol=1)

		}
		girafe(ggobj=g+theme_minimal()+guides(fill=F),height_svg=6,width_svg=9,options=list(opts_selection(type = "single"),opts_sizing(rescale = TRUE)))

		})
	})
	output$umap_title<-renderText({
		dat=rv$obj$obj@meta.data
		t=table(dat$condition)
		if(is.null(input$show)||length(input$show)==2){
			n=sum(t)
		}else{
			if(input$show=='EAE'){
				n=t['E']
			}else if(input$show=='CTL'){
				n=t['C']
			}
		}
		rm(dat)
		if(grepl('cd31',input$ds)){
			return(paste0('UMAP: (CD31+ ECs: ',n,' cells)'))
		}else{
			return(paste0('UMAP: (Total Cells: ',n,' cells)'))
		}
	})

	output$box2_title<-renderText({
		if(input$feat%in%c('celltype','condition')){
			if(is.null(rv$selected)){
				return('Proportions')
			}else{
				return('Differential Expression')
			}
		}else{
			return('Feature Distribution')
		}

	})

	observeEvent(input$box2_selected,{

		updateRadioGroupButtons(session,'feat',selected='gene')
		updateSelectizeInput(session,'gene',selected=input$box2_selected)

	})


}




runApp(list(ui=ui,server=server), launch.browser=TRUE)




