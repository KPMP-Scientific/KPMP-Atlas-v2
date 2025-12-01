pfGroup_color_list = list()
pfGroup_color_list[["Metabolism"]] = "skyblue"
pfGroup_color_list[["TM movement of ions, water and solutes"]] = "seagreen3"
pfGroup_color_list[["Signaling networks reg. metabolism/TM movement"]]= "blue"
pfGroup_color_list[["Redox, iron and mitochondrial dynamics"]]= "purple4"
pfGroup_color_list[["Cell structure dynamics (Cytoskeleton, PM and nucleus)"]] = "red1"
pfGroup_color_list[["Cell adhesion"]] = "lightsalmon1"
pfGroup_color_list[["ECM remodeling"]] = "tan4"
pfGroup_color_list[["Signaling networks reg. cell/tissue structure"]] = "darkorange2"
pfGroup_color_list[["Immune response and signaling"]] = "gold2"
pfGroup_color_list[["Gene expression"]] = "magenta1"
pfGroup_color_list[["Intracellular degradation"]] = "magenta4"
pfGroup_color_list[["Cell cycle and DNA dynamics"]] = "gray50"
pfGroup_color_list[["Cell death"]] = "black"
pfGroup_color_list[["Vesicle traffic"]] = "plum"

Generate_plots = function(plots,complete_pdf_fileName,rows_count,cols_count,empty_rows_at_bottom_count=0,column_widths=c())
{#Begin - Generate plots function
  empty_plot = ggplot() + theme_void()
  max_plots_per_figure = cols_count*rows_count;
  if (length(plots)==0) { plots[[length(plots)+1]] = empty_plot }
  if ((empty_rows_at_bottom_count>0)&(empty_rows_at_bottom_count<rows_count))
  {#Begin
    new_plots = list()
    plots_per_page = (rows_count-empty_rows_at_bottom_count)*cols_count
    plots_count_on_current_page=0
    for (indexPlot in 1:length(plots))
    {#Begin
      new_plots[[length(new_plots)+1]] = plots[[indexPlot]]
      plots_count_on_current_page = plots_count_on_current_page + 1
      if (plots_count_on_current_page==plots_per_page)
      {#Begin
        for (indexCols in 1:cols_count)
        {#Begin
          new_plots[[length(new_plots)+1]] = empty_plot
        }#End
        plots_count_on_current_page=0
      }#End
    }#End
    plots = new_plots
  }#End
  
  
  while ((length(plots)%%max_plots_per_figure)!=0)
  {#Begin
    plots[[length(plots)+1]] = empty_plot
  }#End
  
  pdf(complete_pdf_fileName,width=8.27,height=11.69);
  
  figures_count = ceiling(length(plots)/max_plots_per_figure)
  indexF=1
  for (indexF in 1:figures_count)
  {#Begin
    startPlot = (indexF-1)*max_plots_per_figure+1
    endPlot = min(indexF*max_plots_per_figure,length(plots))
    current_plots = plots[startPlot:endPlot]
    length_plot_list = length(current_plots)
    png_width = 3000*cols_count;
    png_height = 500*rows_count;
    png_resolution=250
    if (length(column_widths)>0)
      #{ do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count, widths = unit(column_widths,"null"))) }
    { grid.arrange(grobs = current_plots, nrow = rows_count, ncol = cols_count, align="vh",widths = unit(column_widths, "null"), heights = unit(rep(1, rows_count), "null")) }
    else { do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count)) }
  }#End
  dev.off()
}#End - Generate plots function
