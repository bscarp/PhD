#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Interaction Power Calculator"),

  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 6,
      
      h1("Setup"),
      
      # Input: Selector for choosing outcome type ----
      selectInput(inputId = "s_es",
                  label = "Specify the type of outcome variable:",
                  c("Continuous" = 1, "Dichotomous" = 2, "Time-to-event" = 3),
                  selected = 1),
      
      splitLayout(
        # Input: Selector for choosing main alpha ----
        numericInput(inputId = "s_a",
                     label = "Required alpha:",
                     value = 0.05,
                     max = 1,
                     min = 0),
        
        # Input: Selector for choosing main power ----
        numericInput(inputId = "s_p",
                     label = "Required power:",
                     value = 0.8,
                     max = 1,
                     min = 0)),
      
      # Input: Selector for choosing intervention allocation ----
      selectInput(inputId = "s_ts",
                  label = "Specify uneven split for intervention variable:",
                  c("Yes" = 1, "No" = 0),
                  selected = 0),
      
      # Input: Selector for choosing intervention split ----
      conditionalPanel(condition = "input.s_ts == 1",
                       numericInput(inputId = "s_tsb",
                                    label = "% in intervention group:",
                                    value = 50,
                                    max = 99,min = 1)),
      
      # Input: Selector for choosing subgroup split ----
      selectInput(inputId = "s_ss",
                  label = "Specify uneven split for subgroup variable:",
                  c("Yes" = 1, "No" = 0),
                  selected = 0),
      
      # Input: Selector for choosing subgroup split ----
      conditionalPanel(condition = "input.s_ss == 1",
                       numericInput(inputId = "s_ssb",
                                    label = "% in one subgroup:",
                                    value = 50,
                                    max = 99,min = 1)),
      br(),
      
      h1("Intervention Main Effect"),
      
      # Input: Selector for choosing main effect size ----
      conditionalPanel(condition = "input.s_es == 1",
                       numericInput(inputId = "m_es_c",
                                    label = "Cohen's d:",
                                    value = 0.5,
                                    min = 0.00001)),
      conditionalPanel(condition = "input.s_es == 2",
                       numericInput(inputId = "m_es_d",
                                    label = "% difference in outcome between Intervention and Control group (positive for an increase, negative for a decrease):",
                                    value = 10,
                                    max = 100,min = -100)),
      conditionalPanel(condition = "input.s_es == 3",
                       numericInput(inputId = "m_es_t",
                                    label = "Hazard Ratio:",
                                    value = 0.8,
                                    max = 0.999999,min = 0.00001)),
        
      # Input: Selector for choosing main effect size ----
      conditionalPanel(condition = "input.s_es == 2",
                       numericInput(inputId = "m_es_db",
                                    label = "% with outcome of interest in Control group:",
                                    value = 50,
                                    max = 100,min = 0)),
      conditionalPanel(condition = "input.s_es == 3",
                       numericInput(inputId = "m_es_tb",
                                    label = "% of participants with an event in Control group:",
                                    value = 50,
                                    max = 100,min = 0)),
      br(),
      h1("Interaction Effect"),
      
      # Input: Selector for choosing interaction effect size ----
      conditionalPanel(condition = "input.s_es == 1",
                       numericInput(inputId = "i_es_c",
                                    label = "Cohen's d:",
                                    value = 0.5,
                                    min = 0.00001)),
      conditionalPanel(condition = "input.s_es == 2",
                       numericInput(inputId = "i_es_d",
                                    label = "% difference in outcome between subgroups in the Intervention group:",
                                    value = 10,
                                    max = 100,min = 0)),
      conditionalPanel(condition = "input.s_es == 3",
                       numericInput(inputId = "i_es_t",
                                    label = "Hazard Ratio:",
                                    value = 0.8,
                                    max = 0.999999,min = 0.00001)),
      
      conditionalPanel(condition = "input.s_es == 1",h4("Cohen's d recommended thresholds:")),
      conditionalPanel(condition = "input.s_es == 1",h5("Small: d = 0.2; Medium d = 0.5; Large: d = 0.8")),
      br(),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(width = 6,
              
      h1("Distribution of sample"),
      tableOutput("sample"),
        
      # Output: Verbatim text for data summary ----
      h1("Main effect"),
      verbatimTextOutput("power_m"),
      
      # Output: Verbatim text for data summary ----
      h1("Interaction effect"),
      verbatimTextOutput("power_i"),

      # Output: Verbatim text for data summary ----
      h1("Final result"),
      h4("Minimum number of participants in the study to achieve the specified power in both the main effect and interaction effect:"),
      h3(strong(textOutput("power")))
      
    )
  )
)

server <- function(input, output, session) {
  
  # Power calculations ----
  power_calc_m <- reactive({
    temp = structure(list(
      n = if(input$s_es == 1) {
        MESS::power_t_test(n = NULL, delta = input$m_es_c, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_tsb,(100-input$s_tsb))/min(input$s_tsb,(100-input$s_tsb)))$n
      } else if (input$s_es == 2) {
        MESS::power_prop_test(n = NULL, p1 = (input$m_es_db+input$m_es_d)/100, p2 = input$m_es_db/100, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_tsb,(100-input$s_tsb))/min(input$s_tsb,(100-input$s_tsb)))$n
      } else {
        as.numeric(powerSurvEpi::ssizeCT.default(power = input$s_p, k = input$s_tsb/(100-input$s_tsb), pE = 1 - ((1 - input$m_es_tb/100) ^ input$m_es_t), pC = input$m_es_tb/100, RR = input$m_es_t, alpha = input$s_a))
      },
      delta = input$m_es_c,
      p1 = (input$m_es_db+input$m_es_d)/100,
      p2 = input$m_es_db/100,
      HR = input$m_es_t,
      pE = (1 - ((1 - input$m_es_tb/100) ^ input$m_es_t)),
      pC = input$m_es_tb/100,
      sig.level = input$s_a,
      power = input$s_p,
      alternative = "two.sided",
      note = paste0("n is vector of number in each group; effect size is ", ifelse(input$s_es == 1, "Cohen's d", ifelse(input$s_es == 2, "Proportions", "HR"))),
      method = paste0("Two-sample ",ifelse(input$s_es == 1, "t test", ifelse(input$s_es == 2, "comparison of proportions", "Cox"))," power calculation")), class = "power.htest")
    temp$n = as.numeric(ifelse(rep(length(temp$n)>1,2),temp$n,rep(temp$n,2)))
    temp$n = temp$n * if(input$s_es == 1 & !input$s_ssb == 50) {
      sum(MESS::power_t_test(n = NULL, delta = input$m_es_c, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_ssb,(100-input$s_ssb))/min(input$s_ssb,(100-input$s_sb)))$n)/
        sum(MESS::power_t_test(n = NULL, delta = input$m_es_c, sig.level = input$s_a, power = input$s_p, ratio = 1)$n*2)
    } else if (input$s_es == 2 & !input$s_ssb == 50) {
      sum(MESS::power_prop_test(n = NULL, p1 = (input$m_es_db+input$m_es_d)/100, p2 = input$m_es_db/100, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_ssb,(100-input$s_ssb))/min(input$s_ssb,(100-input$s_ssb)))$n)/
        (MESS::power_prop_test(n = NULL, p1 = (input$m_es_db+input$m_es_d)/100, p2 = input$m_es_db/100, sig.level = input$s_a, power = input$s_p, ratio = 1)$n*2)
    } else if (input$s_es == 3 & !input$s_ssb == 50) {
      sum(as.numeric(powerSurvEpi::ssizeCT.default(power = input$s_p, k = input$s_ssb/(100-input$s_ssb), pE = 1 - ((1 - input$m_es_tb/100) ^ input$m_es_t), pC = input$m_es_tb/100, RR = input$m_es_t, alpha = input$s_a)))/
        (as.numeric(powerSurvEpi::ssizeCT.default(power = input$s_p, k = 1, pE = 1 - ((1 - input$m_es_tb/100) ^ input$m_es_t), pC = input$m_es_tb/100, RR = input$m_es_t, alpha = input$s_a))*2)
    } else {
      1
    }
    temp$delta = if(!input$s_es == 1) {NULL} else {temp$delta}
    temp$p1 = if(!input$s_es == 2) {NULL} else {temp$p1}
    temp$p2 = if(!input$s_es == 2) {NULL} else {temp$p2}
    temp$HR = if(!input$s_es == 3) {NULL} else {temp$HR}
    temp$pE = if(!input$s_es == 3) {NULL} else {temp$pE}
    temp$pC = if(!input$s_es == 3) {NULL} else {temp$pC}
    temp
  })
  power_calc_i_b <- reactive({
    temp = structure(list(
      n = if(input$s_es == 1) {
        MESS::power_t_test(n = NULL, delta = input$i_es_c/2, sig.level = input$s_a, power = input$s_p, ratio = 1)$n
      } else if (input$s_es == 2) {
        MESS::power_prop_test(n = NULL, p1 = (input$m_es_db+sign(input$m_es_d)*input$i_es_d/2)/100, p2 = input$m_es_db/100, sig.level = input$s_a, power = input$s_p, ratio = 1)$n
      } else {
        as.numeric(powerSurvEpi::ssizeCT.default(power = input$s_p, k = 1, pE = 1 - ((1 - input$m_es_tb/100) ^ sqrt(input$i_es_t)), pC = input$m_es_tb/100, RR = sqrt(input$i_es_t), alpha = input$s_a))
      },
      delta = input$i_es_c,
      p1 = (input$m_es_db+sign(input$m_es_d)*input$i_es_d)/100,
      p2 = input$m_es_db/100,
      HR = input$i_es_t,
      pE = (1 - ((1 - input$m_es_tb/100) ^ input$i_es_t)),
      pC = input$m_es_tb/100,
      sig.level = input$s_a,
      power = input$s_p,
      alternative = "two.sided",
      note = paste0("n is vector of number in each group; effect size is ", ifelse(input$s_es == 1, "Cohen's d", ifelse(input$s_es == 2, "Proportions", "HR"))),
      method = paste0("Two-sample ",ifelse(input$s_es == 1, "t test", ifelse(input$s_es == 2, "comparison of proportions", "Cox"))," power calculation")), class = "power.htest")
    temp$n = as.numeric(ifelse(rep(length(temp$n)>1,2),temp$n,rep(temp$n,2)))
    temp$delta = if(!input$s_es == 1) {NULL} else {temp$delta}
    temp$p1 = if(!input$s_es == 2) {NULL} else {temp$p1}
    temp$p2 = if(!input$s_es == 2) {NULL} else {temp$p2}
    temp$HR = if(!input$s_es == 3) {NULL} else {temp$HR}
    temp$pE = if(!input$s_es == 3) {NULL} else {temp$pE}
    temp$pC = if(!input$s_es == 3) {NULL} else {temp$pC}
    temp
  })
  power_calc_i_t <- reactive({
    temp = structure(list(
      n = if(input$s_es == 1) {
        MESS::power_t_test(n = NULL, delta = input$i_es_c/2, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_tsb,(100-input$s_tsb))/min(input$s_tsb,(100-input$s_tsb)))$n
      } else if (input$s_es == 2) {
        MESS::power_prop_test(n = NULL, p1 = (input$m_es_db+sign(input$m_es_d)*input$i_es_d/2)/100, p2 = input$m_es_db/100, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_tsb,(100-input$s_tsb))/min(input$s_tsb,(100-input$s_tsb)))$n
      } else {
        as.numeric(powerSurvEpi::ssizeCT.default(power = input$s_p, k = input$s_tsb/(100-input$s_tsb), pE = 1 - ((1 - input$m_es_tb/100) ^ sqrt(input$i_es_t)), pC = input$m_es_tb/100, RR = sqrt(input$i_es_t), alpha = input$s_a))
      },
      delta = input$i_es_c,
      p1 = (input$m_es_db+sign(input$m_es_d)*input$i_es_d)/100,
      p2 = input$m_es_db/100,
      HR = input$i_es_t,
      pE = (1 - ((1 - input$m_es_tb/100) ^ input$i_es_t)),
      pC = input$m_es_tb/100,
      sig.level = input$s_a,
      power = input$s_p,
      alternative = "two.sided",
      note = paste0("n is vector of number in each group; effect size is ", ifelse(input$s_es == 1, "Cohen's d", ifelse(input$s_es == 2, "Proportions", "HR"))),
      method = paste0("Two-sample ",ifelse(input$s_es == 1, "t test", ifelse(input$s_es == 2, "comparison of proportions", "Cox"))," power calculation")), class = "power.htest")
    temp$n = as.numeric(ifelse(rep(length(temp$n)>1,2),temp$n,rep(temp$n,2)))
    temp$delta = if(!input$s_es == 1) {NULL} else {temp$delta}
    temp$p1 = if(!input$s_es == 2) {NULL} else {temp$p1}
    temp$p2 = if(!input$s_es == 2) {NULL} else {temp$p2}
    temp$HR = if(!input$s_es == 3) {NULL} else {temp$HR}
    temp$pE = if(!input$s_es == 3) {NULL} else {temp$pE}
    temp$pC = if(!input$s_es == 3) {NULL} else {temp$pC}
    temp
  })
  power_calc_i_s <- reactive({
    temp = structure(list(
      n = if(input$s_es == 1) {
        MESS::power_t_test(n = NULL, delta = input$i_es_c/2, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_ssb,(100-input$s_ssb))/min(input$s_ssb,(100-input$s_ssb)))$n
      } else if (input$s_es == 2) {
        MESS::power_prop_test(n = NULL, p1 = (input$m_es_db+sign(input$m_es_d)*input$i_es_d/2)/100, p2 = input$m_es_db/100, sig.level = input$s_a, power = input$s_p, ratio = max(input$s_ssb,(100-input$s_ssb))/min(input$s_ssb,(100-input$s_ssb)))$n
      } else {
        as.numeric(powerSurvEpi::ssizeCT.default(power = input$s_p, k = input$s_ssb/(100-input$s_ssb), pE = 1 - ((1 - input$m_es_tb/100) ^ sqrt(input$i_es_t)), pC = input$m_es_tb/100, RR = sqrt(input$i_es_t), alpha = input$s_a))
      },
      delta = input$i_es_c,
      p1 = (input$m_es_db+sign(input$m_es_d)*input$i_es_d)/100,
      p2 = input$m_es_db/100,
      HR = input$i_es_t,
      pE = (1 - ((1 - input$m_es_tb/100) ^ input$i_es_t)),
      pC = input$m_es_tb/100,
      sig.level = input$s_a,
      power = input$s_p,
      alternative = "two.sided",
      note = paste0("n is vector of number in each group; effect size is ", ifelse(input$s_es == 1, "Cohen's d", ifelse(input$s_es == 2, "Proportions", "HR"))),
      method = paste0("Two-sample ",ifelse(input$s_es == 1, "t test", ifelse(input$s_es == 2, "comparison of proportions", "Cox"))," power calculation")), class = "power.htest")
    temp$n = as.numeric(ifelse(rep(length(temp$n)>1,2),temp$n,rep(temp$n,2)))
    temp$delta = if(!input$s_es == 1) {NULL} else {temp$delta}
    temp$p1 = if(!input$s_es == 2) {NULL} else {temp$p1}
    temp$p2 = if(!input$s_es == 2) {NULL} else {temp$p2}
    temp$HR = if(!input$s_es == 3) {NULL} else {temp$HR}
    temp$pE = if(!input$s_es == 3) {NULL} else {temp$pE}
    temp$pC = if(!input$s_es == 3) {NULL} else {temp$pC}
    temp
  })
  
  # Reset splits if split variable is set to No
  observeEvent(input$s_ts, {
    updateNumericInput(session, "s_tsb", value = ifelse(input$s_ts==1,input$s_tsb,50))
  })
  observeEvent(input$s_ss, {
    updateNumericInput(session, "s_ssb", value = ifelse(input$s_ss==1,input$s_ssb,50))
  })
  
  output$sample = renderTable({freq = array(c((input$s_tsb/100)*(input$s_ssb/100),
                                              (1-input$s_tsb/100)*(input$s_ssb/100),
                                              (input$s_tsb/100)*(1-input$s_ssb/100),
                                              (1-input$s_tsb/100)*(1-input$s_ssb/100)),
                                            dim=c(2,2))
                               dimnames(freq) = list(Intervention = c("Intervention", "Control"), Subgroup = c("Group 1", "Group 2"))
                               freq |> prop.table() |> addmargins(FUN = list('Variable Split' = sum),quiet = TRUE) |> 
                                 apply(MARGIN=2, FUN=scales::percent, accuracy=0.1) |> ftable() |> format(method="compact",quote=FALSE) |> 
                                 data.frame()
                               }
    ,colnames = FALSE)
  
  # Generate a power calculation ----
  output$power_m <- renderPrint({
    temp_m = power_calc_m()
    temp_m$method = paste0(temp_m$method," for main effect")
    temp_m
  })
  
  # Generate a power calculation ----
  output$power_i <- renderPrint({
    temp_i = power_calc_i_s()
    inf_fac = sum(power_calc_i_t()$n)/sum(power_calc_i_b()$n)
    temp_i$n = sum(temp_i$n)*inf_fac*c((input$s_tsb/100)*(input$s_ssb/100),(input$s_tsb/100)*(1-input$s_ssb/100),(1-input$s_tsb/100)*(input$s_ssb/100),(1-input$s_tsb/100)*(1-input$s_ssb/100))
    temp_i$method = paste0(temp_i$method," for interaction effect")
    temp_i$note = "n is vector of number in each group"
    temp_i
  })
  
  # Final sample size
  output$power <- renderText({
    inf_fac = sum(power_calc_i_t()$n)/sum(power_calc_i_b()$n)
    temp = c(sum(power_calc_m()$n),sum(power_calc_i_s()$n)*inf_fac*sum((input$s_tsb/100)*(input$s_ssb/100),(input$s_tsb/100)*(1-input$s_ssb/100),(1-input$s_tsb/100)*(input$s_ssb/100),(1-input$s_tsb/100)*(1-input$s_ssb/100)),max(sum(power_calc_m()$n),sum(power_calc_i_s()$n)*inf_fac*sum((input$s_tsb/100)*(input$s_ssb/100),(input$s_tsb/100)*(1-input$s_ssb/100),(1-input$s_tsb/100)*(input$s_ssb/100),(1-input$s_tsb/100)*(1-input$s_ssb/100))))
    paste0(round(temp[3],0)," (Main n = ",round(temp[1],0),", Interaction n = ",round(temp[2],0),")")
    # if(input$s_tsb < 50 & input$s_ssb < 50) {
    #   max(sum(ceiling(MESS::power_t_test(n = NULL, delta = input$i_es/2, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_ssb)/input$s_ssb)$n)),
    #       sum(ceiling(MESS::power_t_test(n = NULL, delta = input$m_es, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_tsb)/input$s_tsb)$n)))
    #   } else if (input$s_tsb < 50 & input$s_ssb == 50) {
    #     max(ceiling(MESS::power_t_test(n = NULL, delta = input$i_es/2, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_ssb)/input$s_ssb)$n)*2,
    #         sum(ceiling(MESS::power_t_test(n = NULL, delta = input$m_es, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_tsb)/input$s_tsb)$n)))
    #   } else if (input$s_tsb == 50 & input$s_ssb < 50) {
    #     max(sum(ceiling(MESS::power_t_test(n = NULL, delta = input$i_es/2, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_ssb)/input$s_ssb)$n)),
    #         ceiling(MESS::power_t_test(n = NULL, delta = input$m_es, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_tsb)/input$s_tsb)$n)*2)
    #   } else if (input$s_tsb == 50 & input$s_ssb == 50) {
    #     max(ceiling(MESS::power_t_test(n = NULL, delta = input$i_es/2, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_ssb)/input$s_ssb)$n)*2,
    #         ceiling(MESS::power_t_test(n = NULL, delta = input$m_es, sig.level = input$s_a, power = input$s_p, ratio = (100-input$s_tsb)/input$s_tsb)$n)*2)
    #   }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
