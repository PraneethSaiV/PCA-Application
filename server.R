library(shiny)
options(java.parameters = "- Xmx1024m")
library(shinyRGL)
library(pca3d)
library(rgl)
library("png")
library(svglite)
library(R.devices)
library(jpeg)
library(d3heatmap)
library(htmlwidgets)
library("RColorBrewer")
library(rsvg)
library(pheatmap)
library(xlsx)
library(stringr)


shinyServer(function(input,output,session){
  data_frame_input = reactive({
    if(is.null(input$inputfile))     return(NULL)
    (input$inputfile)})
  
  
  data_real_upload=reactive({
    if(is.null(input$inputfile))     return(NULL)
    data_actual=read.xlsx2(data_frame_input()$datapath,sheetIndex = 1,header = T)           
    data_work=data_actual[c("Protein.ID","Gene.names",(subset(colnames(data_actual),substr(colnames(data_actual),start = 1,stop = 3)=="LFG")))]
    colnames(data_work)=gsub(pattern = "\\.",replacement = " ",x = colnames(data_work))
    for (i in 3:ncol(data_work)){
      data_work[,i]=as.character(data_work[,i])
      data_work[,i]=as.numeric(data_work[,i])
    }
    data_work
  })
  data_real=reactive({
    if(is.null(input$inputfile))     return(NULL)
    data_work=data_real_upload()
    if (input$imputation_type=="Zero"){
      data_work[(is.na(data_work))]=0 
    }
    if (input$imputation_type=="Stratified Average"){                   #Data Processing
      
      for (j in 3:ncol(data_work)){
        for (i in 1:nrow(data_work)){
          if(is.na(data_work[i,j])){data_work[i,j]=rowMeans(data_work[i,which(word(colnames(data_work),2)==word(colnames(data_work)[j],2))],na.rm=T)}
        }
      }
      for (j in 3:ncol(data_work)){
        for (i in 1:nrow(data_work)){
          if(is.na(data_work[i,j])){
            data_work[i,j]=rowMeans(data_work[i,3:ncol(data_work)],na.rm=T)
          }
        }
      }
    }
   
    c=data.frame(t(data_work[,3:(ncol(data_work))]))
    colnames(c)=data_work$`Protein ID`
    c$names=word(rownames(c),2)
    c$names=as.numeric(as.factor(c$names))
    
    data_work$pvalues=NA
    for (i in 1:(ncol(c)-1)){
      data_work[i,'pvalues']=round(anova(lm(c[,i]~c[,ncol(c)]))$Pr[1],digits = 4)
    }
    for (i in 1:ncol(c)){
      if(is.na(data_work[i,'pvalues'])){data_work[i,'pvalues']=1}
    }

    data_work
    })

  
  #----------------------------------------------Preview Clean Data
  output$Datafile_show=renderTable({
    data_work=data_real()
    data_work=subset(data_work,data_work$pvalues<=input$cf_level)
    data_work},bordered = T,digits = 4,hover = T,striped = F)

  #------------------------------------------------PCA 3D
  
 
    
  output$NORMAL3d=renderPlot({
    if(is.null(input$inputfile))     return(NULL)
    data_work=data_real()
    data_work=subset(data_work,data_work$pvalues<=input$cf_level)
    data_work$pvalues=NULL
    
    
    #             Clustering
    pca3dfile=prcomp(t(data_work[,3:ncol(data_work)]),center = T)
    
    
    #             Plot Output
    a=pca3d(palette =c("dodgerblue2","#E31A1C", # red
                       "green4",
                       "#6A3D9A", # purple
                       "#FF7F00", # orange
                       "black","gold1",
                       "skyblue2","#FB9A99", # lt pink
                       "palegreen2",
                       "#CAB2D6", # lt purple
                       "#FDBF6F", # lt orange
                       "gray70", "khaki2",
                       "maroon","orchid1","deeppink1","blue1","steelblue4",
                       "darkturquoise","green1","yellow4","yellow3",
                       "darkorange4","brown"),pca3dfile,shape = "c",legend ="topright",show.plane = FALSE,show.shapes = T,radius = 2,group = word(string = colnames(data_work)[3:length(colnames(data_work))],2),
            axe.titles = c(paste0("PC 1","(",summary(pca3dfile)$importance[2,1]*100,"%)"),paste0("PC 2","(",summary(pca3dfile)$importance[2,2]*100,"%)"),paste0("PC 3","(",summary(pca3dfile)$importance[2,3]*100,"%)")))
    a
    
    #            Legend for screenshot
    legend_for_screenshot=a$groups
    color_for_screenshot=a$colors
    
    #            Image Properties for Screenshot
    download_width=input$width
    dpi_set=input$resol
    
    #            Taking Screenshot
    observeEvent(input$download_link_pca3d,{
      withProgress(message = "Saving Plot",value = 0,{rgl.postscript("PostScript_svg.svg",fmt = "svg",drawText = T)
        bitmap=rsvg("PostScript_svg.svg",width = download_width*dpi_set)
        writePNG(bitmap,"Final_png.png",dpi = dpi_set)
        
        edit_THIS=readPNG("Final_png.png",native = F)
        if (input$type_of_download=="PNG"){
          png(filename = paste(input$change_wd,input$file_name,".png"),res = dpi_set,width = dim(edit_THIS)[2],height = dim(edit_THIS)[1])
          par(mfrow=c(1,1),mar=c(0,0,0,0))
          plot(x = NULL,y = NULL,xlim = c(0,dim(edit_THIS)[2]/dpi_set),ylim =c(0,dim(edit_THIS)[1]/dpi_set) , type="n", axes=FALSE, xlab="", ylab="")
          rasterImage(edit_THIS,0,0,(dim(edit_THIS)[2]/dpi_set),(dim(edit_THIS)[1]/dpi_set),interpolate = T)
          legend("topright", legend = legend_for_screenshot, cex=0.5,pch = 15,col = color_for_screenshot)
          file.remove("PostScript_svg.svg")
          file.remove("Final_png.png")
          dev.off()
        }
        if (input$type_of_download=="TIFF"){
          tiff(filename = paste(input$change_wd,input$file_name,".tiff"),res = dpi_set,width = dim(edit_THIS)[2],height = dim(edit_THIS)[1])
          par(mfrow=c(1,1),mar=c(0,0,0,0))
          plot(x = NULL,y = NULL,xlim = c(0,dim(edit_THIS)[2]/dpi_set),ylim =c(0,dim(edit_THIS)[1]/dpi_set) , type="n", axes=FALSE, xlab="", ylab="")
          rasterImage(edit_THIS,0,0,(dim(edit_THIS)[2]/dpi_set),(dim(edit_THIS)[1]/dpi_set),interpolate = T)
          legend("topright", legend = legend_for_screenshot, cex=0.5,pch = 15,col = color_for_screenshot)
          file.remove("PostScript_svg.svg")
          file.remove("Final_png.png")
          dev.off()
        }
        if (input$type_of_download=="JPEG"){
          jpeg(filename = paste(input$change_wd,input$file_name,".jpeg"),res = dpi_set,width = dim(edit_THIS)[2],height = dim(edit_THIS)[1])
          par(mfrow=c(1,1),mar=c(0,0,0,0))
          plot(x = NULL,y = NULL,xlim = c(0,dim(edit_THIS)[2]/dpi_set),ylim =c(0,dim(edit_THIS)[1]/dpi_set) , type="n", axes=FALSE, xlab="", ylab="")
          rasterImage(edit_THIS,0,0,(dim(edit_THIS)[2]/dpi_set),(dim(edit_THIS)[1]/dpi_set),interpolate = T)
          legend("topright", legend = legend_for_screenshot, cex=0.5,pch = 15,col = color_for_screenshot)
          file.remove("PostScript_svg.svg")
          file.remove("Final_png.png")
          dev.off()
        }
        if (input$type_of_download=="SVG"){
          file.remove("Final_png.png")
          file.rename("PostScript_svg.svg",paste(input$change_wd,input$file_name,".svg"))             
        }
        if (input$type_of_download=="PS"){
          file.remove("Final_png.png")
          file.remove("PostScript_svg.svg")
          rgl.postscript(paste(input$change_wd,input$file_name,".ps"),fmt = "ps",drawText = T)
        }
        if (input$type_of_download=="PDF"){
          file.remove("Final_png.png")
          file.remove("PostScript_svg.svg")
          rgl.postscript(paste(input$change_wd,input$file_name,".pdf"),fmt = "pdf",drawText = T)
        }
        incProgress(1)})
    })
    
      })
    #-------------------------------HEATMAP

  plotInput <- reactive({
    if(is.null(input$inputfile))     return(NULL)
    data_work=data_real()
    data_work=subset(data_work,data_work$pvalues<=input$cf_level)
    data_work$pvalues=NULL
    
    if(input$pro==T && input$gene==F){rule=data_work$`Protein ID`}
    if(input$pro==F && input$gene==T){rule=word(data_work$`Gene names`,1)}
    if(input$pro==T && input$gene==T){rule=paste0(word(data_work$`Protein ID`,1)," ",word(data_work$`Gene names`,1,sep = " "))}
    if(input$pro==F && input$gene==F){rule=rownames(data_work)}
    
    d3heatmap(
      x = scale((data_work[,3:ncol(data_work)])),
      colors = colorRampPalette(c("green", "black", "red"))(n = 1000),
      dendrogram = if (input$cluster) "both" else "none",
      labRow=rule,
      labCol = word(colnames(data_work)[3:ncol(data_work)],2,3),
      cexRow =  1
      
    
      )
    })
  
  observeEvent(input$download_heatmap,{
    withProgress(message = "Saving Plot",value = 0,{data_frame_input = reactive({(input$inputfile)})
    data_actual = read.xlsx2(data_frame_input()$datapath,sheetIndex = 1,header = T)           
    data_work=data_actual[c("Protein.ID","Gene.names",(subset(colnames(data_actual),substr(colnames(data_actual),start = 1,stop = 3)=="LFG")))]
    colnames(data_work)=gsub(pattern = "\\.",replacement = " ",x = colnames(data_work))
    for (i in 3:ncol(data_work)){
      data_work[,i]=as.character(data_work[,i])
      data_work[,i]=as.numeric(data_work[,i])
    }
    if (input$imputation_type=="Zero"){
      data_work[(is.na(data_work))]=0 
    }
    if (input$imputation_type=="Stratified Average"){
      
      for (j in 3:ncol(data_work)){
        for (i in 1:nrow(data_work)){
          if(is.na(data_work[i,j])){data_work[i,j]=rowMeans(data_work[i,which(word(colnames(data_work),2)==word(colnames(data_work)[j],2))],na.rm=T)}
        }
      }
      for (j in 3:ncol(data_work)){
        for (i in 1:nrow(data_work)){
          if(is.na(data_work[i,j])){
            data_work[i,j]=rowMeans(data_work[i,3:ncol(data_work)],na.rm=T)
          }
        }
      }
      
    }
    if (input$cluster==T){
      clus=T
    }else{
      clus=F
    }

    if(input$pro==T && input$gene==F){rule=data_work$`Protein ID`}
    if(input$pro==F && input$gene==T){rule=word(data_work$`Gene names`,1)}
    if(input$pro==T && input$gene==T){rule=paste0(word(data_work$`Protein ID`,1)," ",word(data_work$`Gene names`,1,sep = " "))}
    if(input$pro==F && input$gene==F){rule=rownames(data_work)}
    
    
    if (input$type_of_download_new=="PNG"){png(filename = paste(input$change_wd,input$file_name,".png"),res = input$heatmap_resolution,width = input$heatmap_width,height = input$heatmap_height,units = "in")}
    if (input$type_of_download_new=="TIFF"){tiff(filename = paste(input$change_wd,input$file_name,".tiff"),res = input$heatmap_resolution,width = input$heatmap_width,height = input$heatmap_height,units = "in")}
    if (input$type_of_download_new=="JPEG"){jpeg(filename = paste(input$change_wd,input$file_name,".jpeg"),res = input$heatmap_resolution,width = input$heatmap_width,height = input$heatmap_height,units = "in")}
    if (input$type_of_download_new=="PDF"){pdf(file = paste(input$change_wd,input$file_name,".pdf"),width = input$heatmap_width,height = input$heatmap_height,paper = "a4")}
    hmcol=colorRampPalette(c("green","black","red"))(1000)
    par(mar=c(2,2,5,1))
    pheatmap(legend = F,border_color = "black",mat = scale((data_work[,3:ncol(data_work)])),color = hmcol, labels_col = word(colnames(data_work)[3:ncol(data_work)],2,3),labels_row = rule,cluster_rows = clus,cluster_cols = clus,show_rownames = T,show_colnames = T,width = input$heatmap_width,height = input$heatmap_height,fontsize = 6)
    dev.off()
    incProgress(amount = 1)})
  })
  output$heatmap <- renderD3heatmap(plotInput())
})