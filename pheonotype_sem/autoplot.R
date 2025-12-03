plot_maker=function(){
  require(shiny)
  require(shinydashboard)
  require(LorMe)
  require(magrittr)
  require(shinyhelper)
  require(ggplot2)
  require(ggbeeswarm)
  require(ggpubr)
  ls_class=(lapply(ls(envir = .GlobalEnv), function(x){get(x) %>% class()}) %>% as.character())
  select_list=ls(envir = .GlobalEnv)[which(ls_class=="data.frame")]
  color_match=c("Plan1(max:2)","Plan2(max:5)","Plan3(max:6)","Plan4(max:6)","Plan5(max:6)","Plan6('nejm',max:8)","Plan7('npg',max:10)","Plan8(max:12)","Plan9('Paired',max:12)",ls(envir = .GlobalEnv)[which(ls_class=="character")])
  ui = dashboardPage(skin = "blue",
                     dashboardHeader(title ="Auto plots maker"),
                     dashboardSidebar(width = 330,
                                      tabsetPanel(type = "pills",
                                                  ####analysis####
                                                  tabPanel("Data(数据选项)",icon=icon("wand-magic-sparkles",class="fa-solid fa-wand-magic-sparkles fa-spin-pulse"),
                                                           selectInput(inputId="dataset", label="Dataset(输入数据)", choices = select_list),
                                                           fluidRow(
                                                             column(
                                                               width=11,offset=.1,color="black",numericInput("Treatment", "Treatment colnumber(处理在第几列)",value = NULL) %>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                                                                title = "Helper",fade = T,
                                                                                                                                                                                content = c("A character string which indicates group of samples",
                                                                                                                                                                                            "此处需输出指示样本分组处理的字符串")))),
                                                           fluidRow(
                                                             column(
                                                               width=11,offset=.1,color="black",numericInput("Value", "Value colnumber(数值在第几列)",value = NULL) %>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                                                        title = "Helper2",fade = T,
                                                                                                                                                                        content = c("A character string which indicates group of samples",
                                                                                                                                                                                    "此处需输出指示样本分组处理的字符串")))),
                                                           fluidRow(
                                                             column(
                                                               width=11,offset=.1,color="black",numericInput("facet", "facet colnumber(分面标签在第几列)",value = 0) %>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                                                       title = "Helper2",fade = T,
                                                                                                                                                                       content = c("A character string which indicates group of samples",
                                                                                                                                                                                   "此处需输出指示样本分组处理的字符串")))),
                                                           fluidRow(
                                                             column(
                                                               width=11,offset=.1,color="black",textInput("Subset", "Subset condition(筛选条件)",value="") %>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                                                 title = "Helper",fade = T,
                                                                                                                                                                 content = c("Normally set as 1,can rise if too many up and down",
                                                                                                                                                                             "阈值通常设为1,如上下调过多可提高阈值")))),
                                                           fluidRow(
                                                             column(
                                                               width=11,offset=.1,color="black",radioButtons(inputId="Comp", label="Comparison(组间比较)", choiceNames  = list(icon("check"),icon("ban")),choiceValues =list("True","No"),selected = "True",inline = T)
                                                             )
                                                           ),
                                                           fluidRow(
                                                             column(
                                                               width=12,offset=3,color="black",actionButton("runbutton", "Run", class = "btn-lg btn-success",icon=icon("play"))
                                                             )
                                                           )
                                                  ),
                                                  ###figures####
                                                  tabPanel("Figures(作图选项)",icon=icon("images"),
                                                           fluidRow(
                                                             column(
                                                               width=11,offset=.1,color="black",selectInput("Colormatch", "Color schemes(配色方案)",choices = color_match,selected = "Plan9('Paired',max:12)") %>%helper(type = "inline",colour ="yellow",size = "s",easyClose=T, 
                                                                                                                                                                                title = "Helper",fade = T,
                                                                                                                                                                                content = c("Can choose built-in color schemes(Plan1-6),or your given scheme(as character class)",
                                                                                                                                                                                            "可选内置的配色方案(Plan1-6),或你定义的配色方案(以character的类型)")))),
                                                           fluidRow(
                                                             column(width = 6,numericInput(width="400px","bargapdown", "Gap mult(上下方乘间隔)",value = .05, min = 0, max = 10e10,step=.01)),
                                                             column(width = 6,numericInput(width="400px","bargapup", "Gap add(上下方加间隔)",value = .01, min = 0, max = 10e10,step=.01))
                                                           ),
                                                           menuItem(strong("Labs(坐标轴标签)"), tabName = "Axis label", icon = icon("tags"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,numericInput(width = "200px","labsize","Label font size",value = 10,min = 0,max=99,step = 1))),
                                                                    div(strong('Axis labs'),style="text-indent:12px"),
                                                                    fluidRow(
                                                                      column(6,textInput(width="200px","xlab","x-axis",value = "")),
                                                                      column(6,textInput(width="200px","ylab","y-axis",value = "")))
                                                           ),
                                                           menuItem(strong("Comparison(比较参数)"), tabName = "Comp", icon = icon("ranking-star"),startExpanded=F,
                                                                    div(strong('Attention!'),style="text-indent:18px")%>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                title = "Helper",fade = T,style="text-indent:-2px",
                                                                                                                                content = c("Could be unavailable to adjust when use facet!",
                                                                                                                                            "使用分面功能时参数不可统一调整！")),
                                                                    fluidRow(
                                                                      column(12,radioButtons(inputId="comlabel", label="Sig Label form(显著性标签格式)", choiceNames  = list(icon("p"),icon("star-of-life"),icon("ban")),choiceValues =list("psig","pform","rm"),selected = "psig",inline = T))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","compx", "Sig label.x(显著性标签x坐标)",value = 0, min = 0, max = 999,step=.5)),
                                                                      column(6,numericInput(width="200px","compy", "Sig label.y(显著性标签y坐标)",value = 0, min = 0, max = 10e10,step=1))
                                                                    ),
                                                                    div(strong('Letter label(字母标签)'),style="text-indent:18px"),
                                                                    fluidRow(
                                                                      column(12,radioButtons(inputId="comalign", label="Sig aligning(显著性标签对齐方式)", choiceNames  = list("Align","Float"),choiceValues =list("Align","Float"),selected = "Float",inline = T))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","lettery", "coord(坐标)",value = 0, min = 0, max = 10e10,step=1)),
                                                                      column(6,numericInput(width="200px","letteradj", "adjust(调整)",value = 1.3, min = -999, max = 999,step=.1)),
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","sigsize", "Sig Label size(显著性标签大小)",value = 4, min = 0, max = 99,step=1)),
                                                                      column(6,numericInput(width="200px","compsize", "Letters Label size(字母标签大小)",value = 6, min = 0, max = 99,step=1))
                                                                    )      
                                                           ),
                                                           menuItem(strong("Bar(柱图参数)"), tabName = "Bar", icon = icon("chart-column"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","barwidth", "bar width(柱宽度)",value = .5, min = 0, max = 10,step=.1)),
                                                                      column(6,numericInput(width="200px","barsize", "bar size(框线尺寸)",value = .2, min = 0, max = 10,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(12,numericInput(width="300px","baralpha", "Transparency(alpha透明度)",value = .8, min = 0, max = 1,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","errorwidth", "Error width(误差线宽度)",value = .2, min = 0, max = 10,step=.1)),
                                                                      column(6,numericInput(width="200px","errorsize", "Error size(误差线尺寸)",value = .2, min = 0, max = 10,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(width = 12,textInput(inputId="barcolor", label="Bar frame color(边框颜色)",value = "#000000"))
                                                                    )
                                                           ),
                                                           menuItem(strong("Box(箱体参数)"), tabName = "Box", icon = icon("box"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,radioButtons(inputId="boxaes", label="Mapping aesthetic(映射)", choiceNames  = list(icon("fill"),icon("brush")),choiceValues =list("fill","color"),selected = "fill")),
                                                                      column(6,numericInput(width="200px","boxsize", "Box size(箱线尺寸)",value = .2, min = 0, max = 10,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","boxalpha", "Transparency(alpha,透明度)",value = .8, min = 0, max = 1,step=.1)),
                                                                      column(6,numericInput(width="200px","boxwidth", "Box width(箱体宽度)",value = .5, min = 0, max = 10,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(12,checkboxGroupInput(inputId="outlier", label="Remove outlier(移除异常点)?", choiceNames  = list(icon("ban")),choiceValues =list("NA"),selected = "NA"))
                                                                    ),
                                                                    fluidRow(
                                                                      column(width = 12,textInput(inputId="boxcolor", label="Box frame color(边框颜色)",value = "#000000"))
                                                                    )
                                                           ),
                                                           menuItem(strong("Violin(小提琴参数)"), tabName = "Violin", icon = icon("guitar"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,radioButtons(inputId="violinaes", label="Mapping aesthetic(映射)", choiceNames  = list(icon("fill"),icon("brush")),choiceValues =list("fill","color"),selected = "color")),
                                                                      column(6,numericInput(width="200px","violinsize", "Violin size(线尺寸)",value = .5, min = 0, max = 10,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","violinalpha", "Transparency(alpha,透明度)",value = .8, min = 0, max = 1,step=.1)),
                                                                      column(6,numericInput(width="200px","violinwidth", "Violin width(箱体宽度)",value = .6, min = 0, max = 10,step=.1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(width = 12,textInput(inputId="violincolor", label="Violin frame color(边框颜色)",value = "#000000"))
                                                                    )
                                                           ),
                                                           menuItem(strong("Point(点参数)"), tabName = "Point", icon = icon("circle",class="fa-regular fa-circle",style="color: #ffffff;"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(12,radioButtons(inputId="ptform", label="Point form(点的格式)", choiceNames  = list(icon("ban"),"point","jitter","random"),choiceValues =list("NO","point","jitter","quasirandom"),selected = "NO",inline = T))
                                                                    ),
                                                                    fluidRow(
                                                                      column(12,numericInput(width="200px","pch", "Pch(点的样式)",value=16,min=1,max=25,step=1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","alpha", "Transparency(alpha,透明度)",value = .8, min = 0, max = 1,step=.1)),
                                                                      column(6,numericInput(width="200px","ptsize", "Point size(点大小)",value = 1, min = 0, max = 10,step=.5))
                                                                    ),
                                                                    div(strong('Point color'),style="text-indent:12px")%>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                 title = "Helper",fade = T,style="text-indent:-2px",
                                                                                                                                 content = c("Character string(e.g.red)or Hexadecimal code (e.g. #FE5C5C)of color",
                                                                                                                                             "颜色名(如red)或十六进制代码(如#FE5C5C)"))
                                                                    
                                                           ),
                                                           menuItem(strong("Facet(分面参数)"), tabName = "Face", icon = icon("table-columns"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,radioButtons(inputId="scale", label="Scale(标尺归度)", choiceNames  = list("fixed","free","free_x","free_y"),choiceValues =list("fixed","free","free_x","free_y"),selected = "free_y")),
                                                                      column(6,radioButtons(inputId="strip", label="Strip(标签位置)", choiceNames  = list("top","bottom","left","right"),choiceValues =list("top","bottom","left","right"),selected = "top"))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","nrow", "Row num(分面行数)",value = 0, min = 0, max = 99,step=1)),
                                                                      column(6,numericInput(width="200px","ncol", "Col num(分面列数)",value = 0, min = 0, max = 99,step=1))
                                                                    ),
                                                                    fluidRow(
                                                                      column(6,radioButtons(inputId="astable", label="layout(布局)", choiceNames  = list(icon("arrow-trend-down"),icon("arrow-trend-up")),choiceValues =list("T","F"),selected = "T")),
                                                                      column(6,numericInput(width="200px","striptext", "label size(标签大小)",value = 8, min = 0, max = 99,step=1))
                                                                    )
                                                           ),
                                                           menuItem(strong("Legend(图例参数)"), tabName = "Legend", icon = icon("heading"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,textInput(width="200px","legend_title","Legend title",value = "")),
                                                                      column(6,selectInput(width="200px","legend_p","Legend place",choices = c("top","bottom","left","right"),selected = "right"))),
                                                                    div(strong('Specific legend postion'),style="text-indent:12px")%>%helper(type = "inline",colour ="yellow",size = "m",easyClose=T, 
                                                                                                                                             title = "Helper",fade = T,style="text-indent:-2px",
                                                                                                                                             content = c("Set 0-1 and will deactivated Legend place.Default:0,activated when any value is not 0",
                                                                                                                                                         "设置为0-1之间且Legend place参数将失效.默认为0,任意不为0时生效")),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","legend_x","Legend position x",value = 0, min = 0, max = 1,step=.1)),
                                                                      column(6,numericInput(width="200px","legend_y","Legend position y",value = 0, min = 0, max = 1,step=.1))),
                                                                    checkboxGroupInput(inputId="legend_rm", label="Remove legend(移除图例)?", choiceNames  = list(icon("text-slash")),choiceValues =list("T"))
                                                           ),
                                                           menuItem(strong("Axis(坐标轴参数)"), tabName = "Axis", icon = icon("fonticons-fi"),startExpanded=F,
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","tick_length","tick length",value = 0.2, min = 0, max = 10,step=.1)),
                                                                      column(6,numericInput(width="200px","tick_width","tick width",value = 0.2, min = 0, max = 10,step=.1))),
                                                                    fluidRow(
                                                                      column(6,numericInput(width="200px","axis_size","axis text size",value = 8, min = 0, max = 99,step=1)))
                                                           )
                                                  )
                                      )
                     ),
                     ###main####
                     dashboardBody(
                       tabsetPanel(type = "tabs",
                                   tabPanel("Check plots(查看可视化图片)",icon=icon("photo-film"),
                                            tabsetPanel(
                                              tabPanel("Bar",plotOutput("plot", width = "400",height = "300px"),
                                                       fluidRow(
                                                         column(3,numericInput("width1", "Save width",3,min=0,max=99,step = 1)),
                                                         column(3,numericInput("height1", "Save height",3,min=0,max=99,step = 1))),
                                                       fluidRow(
                                                         column(4,downloadButton("download1", "Download Bar(pdf)")),
                                                         column(4,downloadButton("download1.2", "Download Bar(png)"))
                                                       )
                                              ),
                                              tabPanel("Box",plotOutput("plot2", width = "400",height = "300px"),
                                                       fluidRow(
                                                         column(3,numericInput("width2", "Save width",3,min=0,max=99,step = 1)),
                                                         column(3,numericInput("height2", "Save height",3,min=0,max=99,step = 1))),
                                                       fluidRow(
                                                         column(4,downloadButton("download2", "Download Box(pdf)")),
                                                         column(4,downloadButton("download2.2", "Download Box(png)"))
                                                       )
                                              ),
                                              tabPanel("Violin",plotOutput("plot3", width = "400",height = "300px"),
                                                       fluidRow(
                                                         column(3,numericInput("width3", "Save width",3,min=0,max=99,step = 1)),
                                                         column(3,numericInput("height3", "Save height",3,min=0,max=99,step = 1))),
                                                       fluidRow(
                                                         column(4,downloadButton("download3", "Download Violin(pdf)")),
                                                         column(4,downloadButton("download3.2", "Download Violin(png)"))
                                                       )
                                              )
                                            ),
                                            tabsetPanel(
                                              tabPanel("Comparison report(分析报告)",
                                                       verbatimTextOutput("stat")),
                                              tabPanel("Bar script(柱图脚本)",verbatimTextOutput("codebar")),
                                              tabPanel("Box script(箱图脚本)",verbatimTextOutput("codebox")),
                                              tabPanel("Violin script(琴图脚本)",verbatimTextOutput("codeviolin"))
                                            )
                                   ),
                                   tabPanel("Check input tables(查看输入表格)",icon = icon("table"),
                                            tabsetPanel(
                                              tabPanel("Input Table(输入表格)",dataTableOutput("inputable")),
                                              tabPanel("Subset table(筛选后表格)",dataTableOutput("inputable2"))
                                            )
                                   )
                       )
                     )
  )
  ###server####
  server = function(input, output, session){
    observe_helpers()
    ##input reactive####
    input1=reactive({get(input$dataset)})
    sub_input1=reactive({
      if(input$Subset!=""){
        get(input$dataset)%>% .[eval(parse(text=paste0(input$dataset,"$",input$Subset))),]
      }else{
        get(input$dataset) 
      }
    })
    Treatment=reactive({sub_input1()%>%.[,input$Treatment]}) 
    Value=reactive({sub_input1()%>%.[,input$Value]})
    Face=reactive({sub_input1()%>%.[,input$facet]})
    output$inputable=renderDataTable(input1(), options = list(pageLength = 5))
    output$inputable2=renderDataTable(sub_input1(), options = list(pageLength = 5))
    ylabname=reactive({
      if(input$ylab==""){
        sub_input1() %>% colnames() %>% .[input$Value]
      }else{
        input$ylab
      }
    })
    facet_row=reactive({if(input$nrow==0){NULL}else{input$nrow}})
    facet_col=reactive({if(input$ncol==0){NULL}else{input$ncol}})
    comp_analysis=reactive({
      if(input$facet==0){
        results=auto_signif_test(data =sub_input1(),treatment_col = input$Treatment,value_col = input$Value,prior = T)
        if(unique(as.character(Treatment())) %>% length()>2){
          if(input$lettery==0){
            if(input$comalign=="Align"){
              data.frame(results$comparison_letters,letterp=max(input$letteradj*(results$comparison_letters %>% .[,"Mean"])+1.3*(results$comparison_letters %>% .[,"std"])))
            }else{
              data.frame(results$comparison_letters,letterp=input$letteradj*(results$comparison_letters %>% .[,"Mean"])+1.3*(results$comparison_letters %>% .[,"std"]))
            }
          }else{
            data.frame(results$comparison_letters,letterp=input$lettery)
          }
        }else{results$basicdata}
      }else{
        facet_compare=data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]
        for(i in unique(sub_input1() %>% .[,input$facet])){
          sub_facet=sub_input1() %>% .[which((sub_input1() %>% .[,input$facet])==i),]
          results=auto_signif_test(data =sub_facet,treatment_col = input$Treatment,value_col = input$Value,prior = T)
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                facet_compare=rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=max(input$letteradj*(results$comparison_letters %>% .[,"Mean"])+1.3*(results$comparison_letters %>% .[,"std"])),facetlabel=i))
              }else{
                facet_compare=rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=input$letteradj*(results$comparison_letters %>% .[,"Mean"])+1.3*(results$comparison_letters %>% .[,"std"]),facetlabel=i))
              }
            }else{
              facet_compare=rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=input$lettery,facetlabel=i))
            }
          }else{facet_compare}
        }
        colnames(facet_compare)[which(colnames(facet_compare)=="facetlabel")]=colnames(sub_input1())[input$facet]
        facet_compare
      }
    })
    comp_analysis_print=reactive({
      if(input$facet==0){
        results=auto_signif_test(data =sub_input1(),treatment_col = input$Treatment,value_col = input$Value,prior = T)
      }else{
        for(i in unique(sub_input1() %>% .[,input$facet])){
          cat(paste0("\n\n###",i,"####\n"))
          sub_facet=sub_input1() %>% .[which((sub_input1() %>% .[,input$facet])==i),]
          results=auto_signif_test(data =sub_facet,treatment_col = input$Treatment,value_col = input$Value,prior = T)
        } 
      }
    })
    com_label=reactive({
      if(input$comlabel=="pform"){
        aes(label = after_stat(p.signif))
      }else if(input$comlabel=="psig"){
        aes(label = paste0("p = ", after_stat(p.format))) 
      }
    })
    com_label_p1=reactive({
      if(input$comlabel=="pform"){
        aes(x=as.factor(Treatment()),y=Value(),label = after_stat(p.signif))
      }else if(input$comlabel=="psig"){
        aes(x=as.factor(Treatment()),y=Value(),label = paste0("p = ", after_stat(p.format))) 
      }
    })
    meanframe=reactive({
      if(input$facet==0){
      mean_frame=aggregate(Value(),by=list(Treatment()),FUN=mean)
      Sd=aggregate(Value(),by=list(Treatment()),FUN=sd) %>% .[,"x"]
      Treatment_Name=mean_frame$Group.1
      N=table(Treatment()) %>% as.numeric()
      Mean=mean_frame[,"x"]
      SEM=Sd/(N^0.5)
      data.frame(Treatment_Name,N,Mean,Sd,SEM)
      }else{
        mean_frame=aggregate(Value(),by=list(Treatment(),Face()),FUN=mean)
        Sd=aggregate(Value(),by=list(Treatment(),Face()),FUN=sd) %>% .[,"x"]
        Treatment_Name=mean_frame$Group.1
        N=table(Treatment()) %>% as.numeric()
        Mean=mean_frame[,"x"]
        SEM=Sd/(N^0.5)
        out=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,"Group.2"])
        colnames(out)[6]=colnames(sub_input1())[input$facet]
        out
      }
    })
    colorscheme=reactive({
      if(input$Colormatch=="Plan1(max:2)"){
        c("#E69F00","#56B4E9","black","black","black","black")
      }else if(input$Colormatch=="Plan2(max:5)"){
        c("#FE5D5D","#71C9DD","#33B39F","#6376A0","#F5AF98","black","black")
      }else if(input$Colormatch=="Plan3(max:6)"){
        c("#35A585","#EAE48E","#006FB0","#CC78A6","#F2C661","#56B4E9") 
      }else if(input$Colormatch=="Plan4(max:6)"){
        c("#f49128","#194a55","#187c65","#f26115","#c29f62","#83ba9e") 
      }else if(input$Colormatch=="Plan5(max:6)"){
        c("#c62d17","#023f75","#ea894e","#266b69","#eb4601","#f6c619")
      }else if(input$Colormatch=="Plan6('nejm',max:8)"){
        c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
      }else if(input$Colormatch=="Plan7('npg',max:10)"){
        c("#E64B35FF","#4DBBD5FF","#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF","#B09C85FF")
      }else if(input$Colormatch=="Plan8(max:12)"){
        c("#6a73cf","#edd064","#0eb0c8","#f2ccac","#a1d5b9","#e1abbc","#fa6e01","#2f2f2f","#972b1d","#e6a84b","#4c211b","#ff717f")
      }else if(input$Colormatch=="Plan9('Paired',max:12)"){
        c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99" ,"#E31A1C" ,"#FDBF6F", "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99","#B15928")
      }else{get(input$Colormatch)}
    })
    #analysis&plot####
    observeEvent(input$runbutton,{ 
      output$stat=renderPrint({comp_analysis_print()})
      plot0=reactive({
        p=ggplot(sub_input1(),aes(x=as.factor(Treatment()),y=Value()))+
          scale_y_continuous(expand = c(input$bargapdown,input$bargapup))+
          scale_color_manual(values=colorscheme())+
          scale_fill_manual(values=colorscheme())+
          theme_zg()
        if(input$facet!=0){p=p+facet_wrap(~get(colnames(sub_input1())[input$facet]),nrow = facet_row(),ncol=facet_col(), #as.factor(Face())
                                          scales = input$scale,strip.position = input$strip,
                                          as.table = input$astable)}else{p}
        if(input$striptext!=8){p=p+theme(strip.text.x = element_text(size=input$striptext))}else{p}
      })
      plot1=reactive({
        p=ggplot(meanframe(),aes(x=as.factor(Treatment_Name),y=Mean))+
          theme_zg()+
          geom_bar(stat = 'identity',size=input$barsize,width=input$barwidth,aes(fill=as.factor(Treatment_Name)),color=input$barcolor,alpha=input$baralpha)+
          geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=input$errorsize,width=input$errorwidth)+
          labs(x=input$xlab,y=ylabname(),color=input$legend_title,fill=input$legend_title)+
          scale_y_continuous(expand = c(input$bargapdown,input$bargapup))+
          scale_color_manual(values=colorscheme())+
          scale_fill_manual(values=colorscheme())+
          theme(legend.position = input$legend_p)
        if(input$facet!=0){p=p+facet_wrap(~get(colnames(sub_input1())[input$facet]),nrow = facet_row(),ncol=facet_col(), #as.factor(Face())
                                          scales = input$scale,strip.position = input$strip,
                                          as.table = input$astable)}else{p}
        if(input$striptext!=8){p=p+theme(strip.text.x = element_text(size=input$striptext))}else{p}
        if(input$ptform=="point"){
          p=p+geom_point(data=sub_input1(),aes(x=as.factor(Treatment()),y=Value(),color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F)
        }else if(input$ptform=="jitter"){
          p=p+geom_jitter(data=sub_input1(),aes(x=as.factor(Treatment()),y=Value(),color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F) 
        }else if(input$ptform=="quasirandom"){
          p=p+geom_quasirandom(data=sub_input1(),aes(x=as.factor(Treatment()),y=Value(),color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F) 
        }else if(input$ptform=="NO"){p}
        if(input$Comp=="True"){
          if(unique(as.character(Treatment())) %>% length()>2){
            p=p+geom_text(data=comp_analysis(),aes(x=as.factor(compare),y=letterp,label=Letters),size=input$compsize)
            if(input$comlabel!="rm"){
              if(input$compx!=0|input$compy!=0){
                p=p+stat_compare_means(data=sub_input1(),com_label_p1(),label.x = input$compx,label.y = input$compy,size=input$sigsize)
              }else{p=p+stat_compare_means(data=sub_input1(),com_label_p1(),size=input$sigsize,label.x.npc = "left",label.y.npc = "top")}
            }else{p}
          }else{
            if(input$comlabel!="rm"){
              if(input$compx!=0|input$compy!=0){
                p=p+stat_compare_means(data=sub_input1(),com_label_p1(),label.x = input$compx,label.y = input$compy,size=input$sigsize)
              }else{p=p+stat_compare_means(data=sub_input1(),com_label_p1(),size=input$sigsize,label.x.npc = "left",label.y.npc = "top")}
            }else{p}
          }
        }else if(input$Comp=="NO"){p}
        if(input$legend_x!=0|input$legend_y!=0){p=p+theme(legend.position=c(input$legend_x,input$legend_y))}
        if(input$labsize!=10){p=p+theme(axis.title = element_text(size=input$labsize))}
        if(input$tick_length!=.2){p=p+theme(axis.ticks.length = unit(input$tick_length,"lines"))}
        if(input$tick_width!=.2){p=p+theme(axis.ticks = element_line(size=input$tick_width))}
        if(input$axis_size!=8){p=p+theme(axis.text=element_text(size=input$axis_size))}
        if(is.null(input$legend_rm)){p}else(p=p+guides(color="none",fill="none"))
      })
      plot2=reactive({
        if(input$boxaes=="fill"){
          p=plot0()+
            geom_boxplot(aes(fill=factor(Treatment())),alpha=input$boxalpha,width=input$boxwidth,outlier.color = input$outlier,color=input$boxcolor,linewidth=input$boxsize)+
            labs(x=input$xlab,y=ylabname(),fill=input$legend_title)+
            theme(legend.position = input$legend_p)
        } else if (input$boxaes=="color"){
          p=plot0()+
            geom_boxplot(aes(color=factor(Treatment())),alpha=input$boxalpha,width=input$boxwidth,outlier.color = input$outlier,linewidth=input$boxsize)+
            labs(x=input$xlab,y=ylabname(),color=input$legend_title)+
            theme(legend.position = input$legend_p)
        }
        if(input$ptform=="point"){
          p=p+geom_point(aes(color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F)
        }else if(input$ptform=="jitter"){
          p=p+geom_jitter(aes(color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F) 
        }else if(input$ptform=="quasirandom"){
          p=p+geom_quasirandom(aes(color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F) 
        }else if(input$ptform=="NO"){p}
        if(input$Comp=="True"){
          if(unique(as.character(Treatment())) %>% length()>2){
            p=p+geom_text(data=comp_analysis(),aes(x=as.factor(compare),y=letterp,label=Letters),size=input$compsize)
            if(input$comlabel!="rm"){
              if(input$compx!=0|input$compy!=0){
                p=p+stat_compare_means(data=sub_input1(),com_label_p1(),label.x = input$compx,label.y = input$compy,size=input$sigsize)
              }else{p=p+stat_compare_means(data=sub_input1(),com_label_p1(),size=input$sigsize,label.x.npc = "left",label.y.npc = "top")}
            }else{p}
          }else{
            if(input$comlabel!="rm"){
              if(input$compx!=0|input$compy!=0){
                p=p+stat_compare_means(data=sub_input1(),com_label_p1(),label.x = input$compx,label.y = input$compy,size=input$sigsize)
              }else{p=p+stat_compare_means(data=sub_input1(),com_label_p1(),size=input$sigsize,label.x.npc = "left",label.y.npc = "top")}
            }else{p}
          }
        }else if(input$Comp=="NO"){p}
        if(input$legend_x!=0|input$legend_y!=0){p=p+theme(legend.position=c(input$legend_x,input$legend_y))}
        if(input$labsize!=10){p=p+theme(axis.title = element_text(size=input$labsize))}
        if(input$tick_length!=.2){p=p+theme(axis.ticks.length = unit(input$tick_length,"lines"))}
        if(input$tick_width!=.2){p=p+theme(axis.ticks = element_line(size=input$tick_width))}
        if(input$axis_size!=8){p=p+theme(axis.text=element_text(size=input$axis_size))}
        if(is.null(input$legend_rm)){p}else(p=p+guides(color="none",fill="none"))
      })
      plot3=reactive({
        if(input$violinaes=="fill"){
          p=plot0()+
            geom_violin(aes(fill=factor(Treatment())),alpha=input$violinalpha,width=input$violinwidth,color=input$violincolor,trim = F,linewidth=input$violinsize)+
            labs(x=input$xlab,y=ylabname(),fill=input$legend_title)+
            theme(legend.position = input$legend_p)
        }else if(input$violinaes=="color"){
          p=plot0()+
            geom_violin(aes(color=factor(Treatment())),alpha=input$violinalpha,width=input$violinwidth,trim = F,linewidth=input$violinsize)+
            labs(x=input$xlab,y=ylabname(),color=input$legend_title)+
            theme(legend.position = input$legend_p)
        }
        if(input$ptform=="point"){
          p=p+geom_point(aes(color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F)
        }else if(input$ptform=="jitter"){
          p=p+geom_jitter(aes(color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F) 
        }else if(input$ptform=="quasirandom"){
          p=p+geom_quasirandom(aes(color=as.factor(Treatment())),alpha=input$alpha,size=input$ptsize,pch=input$pch,show.legend = F) 
        }else if(input$ptform=="NO"){p}
        if(input$Comp=="True"){
          if(unique(as.character(Treatment())) %>% length()>2){
            p=p+geom_text(data=comp_analysis(),aes(x=as.factor(compare),y=letterp,label=Letters),size=input$compsize)
            if(input$comlabel!="rm"){
              if(input$compx!=0|input$compy!=0){
                p=p+stat_compare_means(data=sub_input1(),com_label_p1(),label.x = input$compx,label.y = input$compy,size=input$sigsize)
              }else{p=p+stat_compare_means(data=sub_input1(),com_label_p1(),size=input$sigsize,label.x.npc = "left",label.y.npc = "top")}
            }else{p}
          }else{
            if(input$comlabel!="rm"){
              if(input$compx!=0|input$compy!=0){
                p=p+stat_compare_means(data=sub_input1(),com_label_p1(),label.x = input$compx,label.y = input$compy,size=input$sigsize)
              }else{p=p+stat_compare_means(data=sub_input1(),com_label_p1(),size=input$sigsize,label.x.npc = "left",label.y.npc = "top")}
            }else{p}
          }
        }else if(input$Comp=="NO"){p}
        if(input$legend_x!=0|input$legend_y!=0){p=p+theme(legend.position=c(input$legend_x,input$legend_y))}
        if(input$labsize!=10){p=p+theme(axis.title = element_text(size=input$labsize))}
        if(input$tick_length!=.2){p=p+theme(axis.ticks.length = unit(input$tick_length,"lines"))}
        if(input$tick_width!=.2){p=p+theme(axis.ticks = element_line(size=input$tick_width))}
        if(input$axis_size!=8){p=p+theme(axis.text=element_text(size=input$axis_size))}
        if(is.null(input$legend_rm)){p}else(p=p+guides(color="none",fill="none"))
      })
      output$plot <- renderPlot({plot1()})
      output$plot2 <- renderPlot({plot2()})
      output$plot3 <- renderPlot({plot3()})
      output$codebox= renderPrint({
        #code####
        cat("#Running script(运行脚本)\n")
        cat("#Remember to copy or print to consle! (记得复制或打印至控制台喔！)\n")
        #color set###
        cat("###color set####\n")
        cat("my_col<-c(")
        i=1
        while(i < (unique(Treatment()) %>% length())){
          paste0("'",colorscheme()[i],"',") %>% cat()
          i=i+1
        }
        paste0("'",colorscheme() %>% .[(unique(Treatment()) %>% length())],"'") %>% cat()
        cat(")\n")
        cat("names(my_col)<-c(")
        if(class(unique(Treatment()))=="factor"){
          i=1
          while(i < (unique(Treatment()) %>% length())){
            paste0("'",unique(Treatment()) %>% levels() %>% .[i],"',")%>% cat()
            i=i+1
          }
          paste0("'",unique(Treatment()) %>% levels() %>% .[(unique(Treatment()) %>% length())],"'")%>% cat()
          cat(")\n")
        }else{
          i=1
          while(i < (unique(Treatment()) %>% length())){
            paste0("'",unique(Treatment())  %>% .[i],"',")%>% cat()
            i=i+1
          }
          paste0("'",unique(Treatment()) %>% .[(unique(Treatment()) %>% length())],"'")%>% cat()
          cat(")\n")
        }
        ##data input###
        if(input$Subset==""){
          cat("inputdata<-",input$dataset,"\n")
        }else{
          paste0("inputdata<- subset(",input$dataset,",",input$Subset,")") %>% cat()
          cat("\n")
        }
        ##comparison###
        if(input$facet==0){
          paste0("signif_results<- auto_signif_test(data =inputdata,treatment_col =",input$Treatment,",value_col =",input$Value,",prior = T)") %>% cat()
          cat("\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                paste0("letters<- data.frame(signif_results$comparison_letters,letterp=max(",input$letteradj,"*(signif_results$comparison_letters %>% .[,'Mean'])+1.3*(signif_results$comparison_letters %>% .[,'std']))") %>% cat()
              }else{
                paste0("letters<- data.frame(signif_results$comparison_letters,letterp=",input$letteradj,"*(signif_results$comparison_letters %>% .[,'Mean'])+1.3*(signif_results$comparison_letters %>% .[,'std']))") %>% cat()
              }
            }else{
              paste0("letters<- data.frame(signif_results$comparison_letters,letterp=",input$lettery,")") %>% cat()
            }
          }
        }else{
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]\n")
          }
          paste0("for(i in unique(inputdata[,",input$facet,"])){") %>% cat()
          cat("\n ", paste0("sub_facet<-inputdata[which((inputdata[,",input$facet,"])==i),]"),"\n")
          paste0("  results<-auto_signif_test(data =sub_facet,treatment_col =",input$Treatment,",value_col =",input$Value,",prior = T)") %>% cat()
          cat("\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                paste0("  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=max(",input$letteradj,"*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std'])),facetlabel=i))") %>% cat()
              }else{
                paste0("  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=",input$letteradj,"*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std']),facetlabel=i))") %>% cat()
              }
            }else{
              paste0("  letters<- data.frame(results$comparison_letters[,1:5],letterp=",input$lettery,",facetlabel=i)") %>% cat()
            }
          }
          cat("\n}\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("\n")
            paste0("colnames(facet_compare)[7]=colnames(inputdata)[",input$facet,"]") %>% cat()
          }
        }
        cat("\n")
        ##plot###
        ##p0
        paste0("ggplot(inputdata,aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"]))+") %>% cat()
        cat("\n  ")
        paste0("scale_y_continuous(expand = c(",input$bargapdown,",",input$bargapup,"))+") %>% cat()
        cat("\n  scale_color_manual(values=my_col)+")
        cat("\n  scale_fill_manual(values=my_col)+")
        cat("\n  theme_zg()")
        if(input$facet!=0){
          cat("+\n  ")
          paste0("facet_wrap(~get(colnames(inputdata)[",input$facet,"]),scales ='",input$scale,"',strip.position ='",input$strip,"',as.table =",input$astable) %>% cat()
          if(input$nrow!=0){paste0(",nrow=",input$nrow)%>% cat()}
          if(input$ncol!=0){paste0(",ncol=",input$ncol)%>% cat()}
          cat(")")
          if(input$striptext!=8){
            cat("+\n  ")
            paste0("theme(strip.text.x = element_text(size=",input$striptext,"))") %>% cat()
          }
        }
        #main element
        if(input$boxaes=="fill"){
          cat("+\n  ")
          paste0("geom_boxplot(aes(fill=factor(inputdata[,",input$Treatment,"])),alpha=",input$boxalpha,",width=",input$boxwidth,",outlier.color =",input$outlier,",color='",input$boxcolor,"',linewidth=",input$boxsize,")") %>% cat()
          cat("+\n  ")
          paste0("labs(x='",input$xlab,"',fill='",input$legend_title,"'") %>% cat()
        }else if (input$boxaes=="color"){
          cat("+\n  ")
          paste0("geom_boxplot(aes(color=factor(inputdata[,",input$Treatment,"])),alpha=",input$boxalpha,",width=",input$boxwidth,",outlier.color =",input$outlier,",linewidth=",input$boxsize,")")%>% cat()
          cat("+\n  ")
          paste0("labs(x='",input$xlab,"',color='",input$legend_title,"'") %>% cat()
        }
        if(input$ylab==""){
          paste0(",y=colnames(inputdata)[",input$Value,"])") %>% cat()
        }else{
          paste0(",y='",input$ylab,")") %>% cat()
        }  
        #cat("+\n  ")
        #paste0("theme(legend.position = '",input$legend_p,"')") %>% cat()
        #point element
        if(input$ptform=="point"){
          cat("+\n  ")
          paste0("geom_point(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }else if(input$ptform=="jitter"){
          cat("+\n  ")
          paste0("geom_jitter(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }else if(input$ptform=="quasirandom"){
          cat("+\n  ")
          paste0("geom_quasirandom(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }
       #comp element
        if(input$Comp=="True"){
         if(unique(as.character(Treatment())) %>% length()>2){
           cat("+\n  ")
           if(input$facet==0){
             paste0("geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=",input$compsize,")") %>% cat()
             }else{
               paste0("geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=",input$compsize,")") %>% cat()
               }
            }
            if(input$comlabel!="rm"){
              cat("+\n  ")
              paste0("stat_compare_means(data=inputdata,") %>% cat()
              if(input$compx!=0|input$compy!=0){
                if(input$comlabel=="pform"){
                  paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = after_stat(p.signif))") %>% cat()
                }else{
                  paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = paste0('p = ', after_stat(p.format)))") %>% cat()
                }
                paste0(",label.x = ",input$compx,",label.y = ",input$compy,",size=",input$sigsize,")") %>% cat()
              }else{
                if(input$comlabel=="pform"){
                  paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = after_stat(p.signif))") %>% cat()
                }else{
                  paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = paste0('p = ', after_stat(p.format)))") %>% cat()
                }
                paste0(",size=",input$sigsize,",label.x.npc = 'left',label.y.npc = 'top')") %>% cat()
              }
            }
        }
        #theme modify
        cat("+\n  theme(")
        if(input$legend_x!=0|input$legend_y!=0){
          paste0("legend.position=c(",input$legend_x,",",input$legend_y,"),\n") %>% cat()
        }else{
          paste0("legend.position = '",input$legend_p,"',\n") %>% cat()
        }
        if(input$labsize!=10){
          paste0("        axis.title = element_text(size=",input$labsize,"),\n") %>% cat()
        }
        if(input$tick_length!=.2){
          paste0("        axis.ticks.length = unit(",input$tick_length,",'lines'),\n")%>% cat()
        }
        if(input$tick_width!=.2){
          paste0("        axis.ticks = element_line(size=",input$tick_width,"),\n") %>% cat()
        }
        if(input$axis_size!=8){
          paste0("        axis.text=element_text(size=",input$axis_size,"),\n") %>% cat()
        }
        cat("  )")
        if(is.null(input$legend_rm)==F){
          cat("+\n  ")
          cat("guides(color='none',fill='none')")
        }
        #end####
      })
      output$codeviolin= renderPrint({
        #code####
        cat("#Running script(运行脚本)\n")
        cat("#Remember to copy or print to consle! (记得复制或打印至控制台喔！)\n")
        #color set###
        cat("###color set####\n")
        cat("my_col<-c(")
        i=1
        while(i < (unique(Treatment()) %>% length())){
          paste0("'",colorscheme()[i],"',") %>% cat()
          i=i+1
        }
        paste0("'",colorscheme() %>% .[(unique(Treatment()) %>% length())],"'") %>% cat()
        cat(")\n")
        cat("names(my_col)<-c(")
        if(class(unique(Treatment()))=="factor"){
          i=1
          while(i < (unique(Treatment()) %>% length())){
            paste0("'",unique(Treatment()) %>% levels() %>% .[i],"',")%>% cat()
            i=i+1
          }
          paste0("'",unique(Treatment()) %>% levels() %>% .[(unique(Treatment()) %>% length())],"'")%>% cat()
          cat(")\n")
        }else{
          i=1
          while(i < (unique(Treatment()) %>% length())){
            paste0("'",unique(Treatment())  %>% .[i],"',")%>% cat()
            i=i+1
          }
          paste0("'",unique(Treatment()) %>% .[(unique(Treatment()) %>% length())],"'")%>% cat()
          cat(")\n")
        }
        ##data input###
        if(input$Subset==""){
          cat("inputdata<-",input$dataset,"\n")
        }else{
          paste0("inputdata<- subset(",input$dataset,",",input$Subset,")") %>% cat()
          cat("\n")
        }
        ##comparison###
        if(input$facet==0){
          paste0("signif_results<- auto_signif_test(data =inputdata,treatment_col =",input$Treatment,",value_col =",input$Value,",prior = T)") %>% cat()
          cat("\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                paste0("letters<- data.frame(signif_results$comparison_letters,letterp=max(",input$letteradj,"(signif_results$comparison_letters %>% .[,'Mean'])+1.3*(signif_results$comparison_letters %>% .[,'std'])))") %>% cat()
              }else{
                paste0("letters<- data.frame(signif_results$comparison_letters,letterp=",input$letteradj,"*(signif_results$comparison_letters %>% .[,'Mean'])+1.3*(signif_results$comparison_letters %>% .[,'std']))") %>% cat()
              }
            }else{
              paste0("letters<- data.frame(signif_results$comparison_letters,letterp=",input$lettery,")") %>% cat()
            }
          }
        }else{
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]\n")
          }
          paste0("for(i in unique(inputdata[,",input$facet,"])){") %>% cat()
          cat("\n ", paste0("sub_facet<-inputdata[which((inputdata[,",input$facet,"])==i),]"),"\n")
          paste0("  results<-auto_signif_test(data =sub_facet,treatment_col =",input$Treatment,",value_col =",input$Value,",prior = T)") %>% cat()
          cat("\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                paste0("  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=max(",input$letteradj,"(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std']),facetlabel=i))") %>% cat()
              }else{
                paste0("  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=",input$letteradj,"*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std']),facetlabel=i))") %>% cat()
              }
            }else{
              paste0("  letters<- data.frame(results$comparison_letters[,1:5],letterp=",input$lettery,",facetlabel=i)") %>% cat()
            }
          }
          cat("\n}\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("\n")
            paste0("colnames(facet_compare)[7]=colnames(inputdata)[",input$facet,"]") %>% cat()
          }
        }
        cat("\n")
        ##plot###
        ##p0
        paste0("ggplot(inputdata,aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"]))+") %>% cat()
        cat("\n  ")
        paste0("scale_y_continuous(expand = c(",input$bargapdown,",",input$bargapup,"))+") %>% cat()
        cat("\n  scale_color_manual(values=my_col)+")
        cat("\n  scale_fill_manual(values=my_col)+")
        cat("\n  theme_zg()")
        if(input$facet!=0){
          cat("+\n  ")
          paste0("facet_wrap(~get(colnames(inputdata)[",input$facet,"]),scales ='",input$scale,"',strip.position ='",input$strip,"',as.table =",input$astable) %>% cat()
          if(input$nrow!=0){paste0(",nrow=",input$nrow)%>% cat()}
          if(input$ncol!=0){paste0(",ncol=",input$ncol)%>% cat()}
          cat(")")
          if(input$striptext!=8){
            cat("+\n  ")
            paste0("theme(strip.text.x = element_text(size=",input$striptext,"))") %>% cat()
          }
        }
        #main element
        if(input$violinaes=="fill"){
          cat("+\n  ")
          paste0("geom_violin(aes(fill=factor(inputdata[,",input$Treatment,"])),alpha=",input$violinalpha,",width=",input$violinwidth,",color='",input$violincolor,"',linewidth=",input$violinsize,",trim = F)") %>% cat()
          cat("+\n  ")
          paste0("labs(x='",input$xlab,"',fill='",input$legend_title,"'") %>% cat()
        }else if (input$violinaes=="color"){
          cat("+\n  ")
          paste0("geom_violin(aes(color=factor(inputdata[,",input$Treatment,"])),alpha=",input$violinalpha,",width=",input$violinwidth,",linewidth=",input$violinsize,",trim = F)")%>% cat()
          cat("+\n  ")
          paste0("labs(x='",input$xlab,"',color='",input$legend_title,"'") %>% cat()
        }
        if(input$ylab==""){
          paste0(",y=colnames(inputdata)[",input$Value,"])") %>% cat()
        }else{
          paste0(",y='",input$ylab,"')") %>% cat()
        }  
        #pt element
        if(input$ptform=="point"){
          cat("+\n  ")
          paste0("geom_point(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }else if(input$ptform=="jitter"){
          cat("+\n  ")
          paste0("geom_jitter(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }else if(input$ptform=="quasirandom"){
          cat("+\n  ")
          paste0("geom_quasirandom(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }
        #comp element
        if(input$Comp=="True"){
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("+\n  ")
            if(input$facet==0){
              paste0("geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=",input$compsize,")") %>% cat()
            }else{
              paste0("geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=",input$compsize,")") %>% cat()
            }
          }
          if(input$comlabel!="rm"){
            cat("+\n  ")
            paste0("stat_compare_means(data=inputdata,") %>% cat()
            if(input$compx!=0|input$compy!=0){
              if(input$comlabel=="pform"){
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = after_stat(p.signif))") %>% cat()
              }else{
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = paste0('p = ', after_stat(p.format)))") %>% cat()
              }
              paste0(",label.x = ",input$compx,",label.y = ",input$compy,",size=",input$sigsize,")") %>% cat()
            }else{
              if(input$comlabel=="pform"){
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = after_stat(p.signif))") %>% cat()
              }else{
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = paste0('p = ', after_stat(p.format)))") %>% cat()
              }
              paste0(",size=",input$sigsize,",label.x.npc = 'left',label.y.npc = 'top')") %>% cat()
            }
          }
        }
        #theme modify
        cat("+\n  theme(")
        if(input$legend_x!=0|input$legend_y!=0){
          paste0("legend.position=c(",input$legend_x,",",input$legend_y,"),\n") %>% cat()
        }else{
          paste0("legend.position = '",input$legend_p,"',\n") %>% cat()
        }
        if(input$labsize!=10){
          paste0("        axis.title = element_text(size=",input$labsize,"),\n") %>% cat()
        }
        if(input$tick_length!=.2){
          paste0("        axis.ticks.length = unit(",input$tick_length,",'lines'),\n")%>% cat()
        }
        if(input$tick_width!=.2){
          paste0("        axis.ticks = element_line(size=",input$tick_width,"),\n") %>% cat()
        }
        if(input$axis_size!=8){
          paste0("        axis.text=element_text(size=",input$axis_size,"),\n") %>% cat()
        }
        cat("  )")
        if(is.null(input$legend_rm)==F){
          cat("+\n  ")
          cat("guides(color='none',fill='none')")
        }
        #end####
      })
      output$codebar=renderPrint({
        ###
        #code####
        cat("#Running script(运行脚本)\n")
        cat("#Remember to copy or print to consle! (记得复制或打印至控制台喔！)\n")
        #color set###
        cat("###color set####\n")
        cat("my_col<-c(")
        i=1
        while(i < (unique(Treatment()) %>% length())){
          paste0("'",colorscheme()[i],"',") %>% cat()
          i=i+1
        }
        paste0("'",colorscheme() %>% .[(unique(Treatment()) %>% length())],"'") %>% cat()
        cat(")\n")
        cat("names(my_col)<-c(")
        if(class(unique(Treatment()))=="factor"){
          i=1
          while(i < (unique(Treatment()) %>% length())){
            paste0("'",unique(Treatment()) %>% levels() %>% .[i],"',")%>% cat()
            i=i+1
          }
          paste0("'",unique(Treatment()) %>% levels() %>% .[(unique(Treatment()) %>% length())],"'")%>% cat()
          cat(")\n")
        }else{
          i=1
          while(i < (unique(Treatment()) %>% length())){
            paste0("'",unique(Treatment())  %>% .[i],"',")%>% cat()
            i=i+1
          }
          paste0("'",unique(Treatment()) %>% .[(unique(Treatment()) %>% length())],"'")%>% cat()
          cat(")\n")
        }
        ##data input###
        if(input$Subset==""){
          cat("inputdata<-",input$dataset,"\n")
        }else{
          paste0("inputdata<- subset(",input$dataset,",",input$Subset,")") %>% cat()
          cat("\n")
        }
        
        if(input$facet==0){
          paste0("signif_results<- auto_signif_test(data =inputdata,treatment_col =",input$Treatment,",value_col =",input$Value,",prior = T)") %>% cat()
          cat("\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                paste0("letters<- data.frame(signif_results$comparison_letters,letterp=max(",input$letteradj,"(signif_results$comparison_letters %>% .[,'Mean'])+1.3*(signif_results$comparison_letters %>% .[,'std'])))") %>% cat()
              }else{
                paste0("letters<- data.frame(signif_results$comparison_letters,letterp=",input$letteradj,"*(signif_results$comparison_letters %>% .[,'Mean'])+1.3*(signif_results$comparison_letters %>% .[,'std']))") %>% cat()
              }
            }else{
              paste0("letters<- data.frame(signif_results$comparison_letters,letterp=",input$lettery,")") %>% cat()
            }
          }
          cat("\n")
          paste0("mean_frame<- aggregate(inputdata[,",input$Value,"],by=list(inputdata[,",input$Treatment,"]),FUN=mean)") %>% cat()
          cat("\n")
          paste0("Sd<- aggregate(inputdata[,",input$Value,"],by=list(inputdata[,",input$Treatment,"]),FUN=sd) %>% .[,'x']") %>% cat()
          cat("\n")
          cat("Treatment_Name<- mean_frame$Group.1\n")
          paste0("N<- table(inputdata[,",input$Treatment,"]) %>% as.numeric()") %>% cat()
          cat("\nMean<- mean_frame[,'x']
      SEM<- Sd/(N^0.5)
      input_mean_frame<- data.frame(Treatment_Name,N,Mean,Sd,SEM)")
        }else{
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]\n")
          }
          paste0("for(i in unique(inputdata[,",input$facet,"])){") %>% cat()
          cat("\n ", paste0("sub_facet<-inputdata[which((inputdata[,",input$facet,"])==i),]"),"\n")
          paste0("  results<-auto_signif_test(data =sub_facet,treatment_col =",input$Treatment,",value_col =",input$Value,",prior = T)") %>% cat()
          cat("\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            if(input$lettery==0){
              if(input$comalign=="Align"){
                paste0("  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=max(",input$letteradj,"(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std']),facetlabel=i))") %>% cat()
              }else{
                paste0("  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=",input$letteradj,"*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std']),facetlabel=i))") %>% cat()
              }
            }else{
              paste0("  letters<- data.frame(results$comparison_letters[,1:5],letterp=",input$lettery,",facetlabel=i)") %>% cat()
            }
          }
          cat("\n}\n")
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("\n")
            paste0("colnames(facet_compare)[7]=colnames(inputdata)[",input$facet,"]") %>% cat()
            cat("\n")
          }
          paste0("mean_frame<- aggregate(inputdata[,",input$Value,"],by=list(inputdata[,",input$Treatment,"],inputdata[,",input$facet,"]),FUN=mean)") %>% cat()
          cat("\n")
          paste0("Sd<- aggregate(inputdata[,",input$Value,"],by=list(inputdata[,",input$Treatment,"],inputdata[,",input$facet,"]),FUN=sd)%>% .[,'x']") %>% cat()
          cat("\n")
          cat("Treatment_Name<- mean_frame$Group.1\n")
          paste0("N<- table(inputdata[,",input$Treatment,"]) %>% as.numeric()") %>% cat()
          cat("\nMean<- mean_frame[,'x']
SEM<- Sd/(N^0.5)
input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])\n")
          paste0("colnames(input_mean_frame)[6]=colnames(inputdata)[",input$facet,"]") %>% cat()
        }
        ###plot
        cat("\nggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+")
        cat("\n  theme_zg()+\n  ")
        paste0("geom_bar(stat = 'identity',size=",input$barsize,",width=",input$barwidth,",aes(fill=as.factor(Treatment_Name)),color='",input$barcolor,"',alpha=",input$baralpha,")+") %>% cat()
        cat("\n  ")
        paste0("geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=",input$errorsize,",width=",input$errorwidth,")+") %>% cat()
        cat("\n  ")
        paste0("scale_y_continuous(expand = c(",input$bargapdown,",",input$bargapup,"))+") %>% cat()
        cat("\n  scale_color_manual(values=my_col)+")
        cat("\n  scale_fill_manual(values=my_col)+\n  ")
        paste0("labs(x='",input$xlab,"',fill='",input$legend_title,"'") %>% cat()
        if(input$ylab==""){
          paste0(",y=colnames(inputdata)[",input$Value,"])") %>% cat()
        }else{
          paste0(",y='",input$ylab,"')") %>% cat()
        } 
        if(input$facet!=0){
          cat("+\n  ")
          paste0("facet_wrap(~get(colnames(inputdata)[",input$facet,"]),scales ='",input$scale,"',strip.position ='",input$strip,"',as.table =",input$astable) %>% cat()
          if(input$nrow!=0){paste0(",nrow=",input$nrow)%>% cat()}
          if(input$ncol!=0){paste0(",ncol=",input$ncol)%>% cat()}
          cat(")")
          if(input$striptext!=8){
            cat("+\n  ")
            paste0("theme(strip.text.x = element_text(size=",input$striptext,"))") %>% cat()
          }
        }
        #point element
        if(input$ptform=="point"){
          cat("+\n  ")
          paste0("geom_point(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }else if(input$ptform=="jitter"){
          cat("+\n  ")
          paste0("geom_jitter(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }else if(input$ptform=="quasirandom"){
          cat("+\n  ")
          paste0("geom_quasirandom(aes(color=as.factor(inputdata[,",input$Treatment,"])),alpha=",input$alpha,",size=",input$ptsize,",pch=",input$pch,",show.legend = F)") %>% cat()
        }
        #comp element
        if(input$Comp=="True"){
          if(unique(as.character(Treatment())) %>% length()>2){
            cat("+\n  ")
            if(input$facet==0){
              paste0("geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=",input$compsize,")") %>% cat()
            }else{
              paste0("geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=",input$compsize,")") %>% cat()
            }
          }
          if(input$comlabel!="rm"){
            cat("+\n  ")
            paste0("stat_compare_means(data=inputdata,") %>% cat()
            if(input$compx!=0|input$compy!=0){
              if(input$comlabel=="pform"){
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = after_stat(p.signif))") %>% cat()
              }else{
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = paste0('p = ', after_stat(p.format)))") %>% cat()
              }
              paste0(",label.x = ",input$compx,",label.y = ",input$compy,",size=",input$sigsize,")") %>% cat()
            }else{
              if(input$comlabel=="pform"){
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = after_stat(p.signif))") %>% cat()
              }else{
                paste0("aes(x=as.factor(inputdata[,",input$Treatment,"]),y=inputdata[,",input$Value,"],label = paste0('p = ', after_stat(p.format)))") %>% cat()
              }
              paste0(",size=",input$sigsize,",label.x.npc = 'left',label.y.npc = 'top')") %>% cat()
            }
          }
        }
        #theme modify
        cat("+\n  theme(")
        if(input$legend_x!=0|input$legend_y!=0){
          paste0("legend.position=c(",input$legend_x,",",input$legend_y,"),\n") %>% cat()
        }else{
          paste0("legend.position = '",input$legend_p,"',\n") %>% cat()
        }
        if(input$labsize!=10){
          paste0("        axis.title = element_text(size=",input$labsize,"),\n") %>% cat()
        }
        if(input$tick_length!=.2){
          paste0("        axis.ticks.length = unit(",input$tick_length,",'lines'),\n")%>% cat()
        }
        if(input$tick_width!=.2){
          paste0("        axis.ticks = element_line(size=",input$tick_width,"),\n") %>% cat()
        }
        if(input$axis_size!=8){
          paste0("        axis.text=element_text(size=",input$axis_size,"),\n") %>% cat()
        }
        cat("  )")
        if(is.null(input$legend_rm)==F){
          cat("+\n  ")
          cat("guides(color='none',fill='none')")
        }
        ###end
      })
    })
    output$download1=downloadHandler(
      filename = function() {paste0("Bar",".pdf")},
      content = function(file) {ggsave(file,width = input$width1,height = input$height1)},
      contentType = "pdf"
    )
    output$download1.2=downloadHandler(
      filename = function() {paste0("Bar",".png")},
      content = function(file) {ggsave(file,width = input$width1,height = input$height1)},
      contentType = "png"
    )
    output$download2=downloadHandler(
      filename = function() {paste0("Box",".pdf")},
      content = function(file) {ggsave(file,width = input$width2,height = input$height2)},
      contentType = "pdf"
    )
    output$download2.2=downloadHandler(
      filename = function() {paste0("Box",".png")},
      content = function(file) {ggsave(file,width = input$width2,height = input$height2)},
      contentType = "png"
    )
    output$download3=downloadHandler(
      filename = function() {paste0("Violin",".pdf")},
      content = function(file) {ggsave(file,width = input$width3,height = input$height3)},
      contentType = "pdf"
    )
    output$download3.2=downloadHandler(
      filename = function() {paste0("Violin",".png")},
      content = function(file) {ggsave(file,width = input$width3,height = input$height3)},
      contentType = "png"
    )
  }
  shinyApp(ui, server)
}
