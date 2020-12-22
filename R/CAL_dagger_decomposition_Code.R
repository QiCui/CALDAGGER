

###############################################################
##### Code for CAL_dagger and decomposition
##### To explore, please see the shinyapps
##### https://caldagger.shinyapps.io/CALdagger/?_ga=2.182331879.118419113.1608632229-217553970.1608632229
###############################################################

########################
######  Library and Source
########################
# install.packages("HMDHFDplus")
library(HMDHFDplus)
HMD_account_name <- "XXXXXXXXXXXXXX"
HMD_account_pin  <- "XXXXXXXXXXXXXX"


source("R/CAL_dagger_life_table_functions.R")

########################
###### Interested population
###### Interested truncation year
###### Interested sex
########################

### Interested country
## HMD countries can cover 110 years
## DNK(Denmark), Fin(Finland), FRATNP(France), ITA(Italy), NLD(Netherlands)
## NOR(Norway), SWE(Sweden), CHE(Switzerland), GBR_SCO(Scotland), GBRTENW(England & Wales)
Pop1 <- "DNK"
Pop2 <- "SWE"

### Interested year
##  Starts from 1989 to 2016
Int_Year <- 2013

### Interested sex
## Female or Male
Sex = "Female"

############################################
######  Preparing DATA
############################################

T_Year_series = Int_Year + 1
HMD_Names<-c(Pop1, Pop2)
c_HMD<-1
px_tri_summary <- lx_tri_summary <- c()

repeat{

  ################################################################
  ### Calculating 1qx based on original HMD data
  ################################################################

  Input_Population <- readHMDweb(CNTRY=HMD_Names[c_HMD],item="Population", username=HMD_account_name, password=HMD_account_pin, fixup=F)
  Input_DX_lexis <- readHMDweb(CNTRY=HMD_Names[c_HMD],item="Deaths_lexis", username=HMD_account_name, password=HMD_account_pin, fixup=F)

  # re-format the age. 0 ~ 110+ to 0~110
  Input_DX_lexis$Age<-as.numeric(as.character(Input_DX_lexis$Age))
  Input_DX_lexis$Age[is.na(Input_DX_lexis$Age)]=110
  Input_Population$Age<-as.numeric(as.character(Input_Population$Age))
  Input_Population$Age[is.na(Input_Population$Age)]=110

  Input_Population$Year <- as.numeric(as.character(substr(Input_Population$Year, 1, 4)))
  ry <- as.numeric(which(table(Input_Population$Year)>111))+Input_Population$Year[1]-1

  if (length(ry) !=0) {
    ry_i<-1
    repeat{
      A<-as.numeric(which(grepl(subset(Input_Population,
                                       Input_Population$Year==ry[ry_i])[1,1],
                                Input_Population$Year))[1:1])                        # first
      B<-A+111                                                                       # second

      Pop_rep <- (Input_Population[A:(A+110),]+Input_Population[B:(B+110),])*0.5
      Input_Population <- rbind(Input_Population[0:(A-1),],
                              Pop_rep,
                              Input_Population[(B+111):length(Input_Population[,1]),])

      ry_i=ry_i+1
      if (ry_i> length(ry)) break
    }
  }

  ################################################################
  ### Add the cohort column for period data
  ################################################################

  Population_s<-c()
  cohort_year_start <- seq(T_Year_series[1]-110, T_Year_series[length(T_Year_series)],1)
  ii_year<-1

  repeat{
    Population_s_pre <- subset(Input_Population,
                             Input_Population$Year == cohort_year_start[ii_year]+1)

    Cohort_hypo <- seq(cohort_year_start[ii_year], cohort_year_start[ii_year]-110,-1)
    Population_s_pre <- cbind(Cohort_hypo, Population_s_pre)
    Population_s <- rbind(Population_s, Population_s_pre)

    ii_year=ii_year+1
    if (ii_year==length(cohort_year_start)+1) break
  }


  ################################################################
  ### get the triangle mortality
  ################################################################

  i_year <- 1
  repeat {

    Survival_prob_Summary <- lx_summary <- c()
    target_age <- 110 # we use the mortality rates from birth to 109.99999, because in HMD 110 means "110+"
    Start_year <- T_Year_series[i_year]
    Start_cohort <- Start_year-target_age

    repeat{

      DX_lexis_c <- subset(Input_DX_lexis, Input_DX_lexis$Cohort == Start_cohort &
                           Input_DX_lexis$Age <= target_age &
                           Input_DX_lexis$Year<= Start_year)
      Population_c <- subset(Population_s, Population_s$Cohort_hypo == Start_cohort &
                             Population_s$Age <= target_age &
                             Population_s$Year <= Start_year)

      denominator_p <-  Population_c[, c(4,5)] + DX_lexis_c[seq(1,length(DX_lexis_c[,1])-1,2),c(4,5)]
      numerator_p <- Population_c[, c(4,5)] - DX_lexis_c[seq(2,length(DX_lexis_c[,1]),2),c(4,5)]
      Survival_prob_c <- data.frame(Country = rep(HMD_Names[c_HMD], target_age),
                                    Cohort=rep(Population_c$Cohort_hypo[1], target_age),
                                    Age=seq(0, target_age-1,1),
                                    Female=(numerator_p/denominator_p)$Female,
                                    Male=(numerator_p/denominator_p)$Male)

      lx_summary_f_s <- c(1,cumprod(Survival_prob_c$Female)[1:(length(Survival_prob_c$Female)-1)])
      lx_summary_f <- lx_summary_f_s[target_age-length(which(is.na(lx_summary_f_s) | lx_summary_f_s==0))] # avoid the NA and zero value
      lx_summary_m_s <- c(1,cumprod(Survival_prob_c$Male)[1:(length(Survival_prob_c$Male)-1)])
      lx_summary_m <- lx_summary_m_s[target_age-length(which(is.na(lx_summary_m_s) | lx_summary_m_s==0))] # avoid the NA and zero value

      lx_summary_pre <- data.frame(Country = HMD_Names[c_HMD],
                                   Cohort=Population_c$Cohort_hypo[1],
                                   Age=target_age-1,
                                   Female=lx_summary_f,
                                   Male=lx_summary_m)

      Survival_prob_Summary <- rbind(Survival_prob_Summary, Survival_prob_c)
      lx_summary <- rbind(lx_summary,lx_summary_pre)

      Start_cohort=Start_cohort+1
      target_age=target_age-1
      if (Start_cohort>=Start_year) break
    }

    lx_summary$Age <- rev(lx_summary$Age)

    i_year=i_year+1
    if (i_year==length(T_Year_series)+1) break
  }

  px_tri_summary <- rbind(px_tri_summary, Survival_prob_Summary)
  lx_tri_summary <- rbind(lx_tri_summary, lx_summary)

  c_HMD=c_HMD+1
  if (c_HMD==length(HMD_Names)+1) break
}



############################################
######  Calculation CAL_dagger & decomposition
############################################


## re-format the sex variable
if (Sex=="Female") {
  Int_sex = c(1,2,4)
} else {
  Int_sex = c(1,2,5)
}

###############
###  Calculation CAL_dagger
###############
## formula based on lx

lx_pop1 <- subset(lx_tri_summary, lx_tri_summary$Country==Pop1)[, Int_sex]
lx_pop2 <- subset(lx_tri_summary, lx_tri_summary$Country==Pop2)[, Int_sex]

CAL_dagger_pop1 <- -sum(log(lx_pop1[,3])* lx_pop1[,3])
CAL_dagger_pop2 <- -sum(log(lx_pop2[,3])* lx_pop2[,3])

###############
###  Calculation decomposition
###############

cal_px_1 <- subset(px_tri_summary, px_tri_summary$Country == Pop1)[, Int_sex]
cal_px_2 <- subset(px_tri_summary, px_tri_summary$Country == Pop2)[, Int_sex]

length_age <- 110
Year_cohort<-unique(cal_px_1$Cohort)

i_year<-1
age_contribution_all_cohort<-cum_matrix<-c()
repeat{

  px_1<-subset(cal_px_1,cal_px_1$Cohort==Year_cohort[i_year])[,3]
  px_2<-subset(cal_px_2,cal_px_2$Cohort==Year_cohort[i_year])[,3]

  log_p<-log(px_1/px_2)
  mid_lx<-((c(1,cumprod(px_1))[length(px_1)])+(c(1,cumprod(px_2))[length(px_2)]))*0.5

  cohort_contribution<- (-log_p*mid_lx*(1+log(mid_lx)))
  cohort_contribution[is.na(cohort_contribution)|cohort_contribution>999|cohort_contribution<(-999)]=0

  cum_matrix_pre<-cumsum(cohort_contribution)
  cum_matrix<-cbind(cum_matrix, rev(c(cum_matrix_pre,rep(0,110-length(cohort_contribution)))))

  cohort_contribution<-c(cohort_contribution, rep(0,(110-length(cohort_contribution))))
  cohort_contribution<-rev(cohort_contribution)

  age_contribution_all_cohort<-cbind(age_contribution_all_cohort, cohort_contribution)

  i_year=i_year+1
  if (i_year==length(Year_cohort)+1) break
}

###############
###  reformat the decomposition results
###############

## lexis format
i_r=1
Lexis_age_con<-c()
repeat{
  Lexis_age_con_pre<-c(rep(0, 110-i_r),age_contribution_all_cohort[i_r,])[1:110]
  Lexis_age_con<-rbind(Lexis_age_con,Lexis_age_con_pre)

  i_r=i_r+1
  if (i_r==length_age+1) break
}

## filled_contour format
j_r=1
FCF_age_con<-c()
repeat{
  FCF_age_con_pre<-Lexis_age_con[j_r,]
  FCF_age_con<-cbind(FCF_age_con_pre,FCF_age_con)

  j_r=j_r+1
  if (j_r==length_age+1) break
}

############################################
######  Results and visualisation
############################################

###############
####### CAL_dagger result
###############
data.frame(Year = Int_Year, Sex = Sex,
           Pop1 = Pop1, Pop2 = Pop2,
           CAL_dagger_pop1 = CAL_dagger_pop1, CAL_dagger_pop2 = CAL_dagger_pop2,
           Difference_Pop1_Pop2 = CAL_dagger_pop1 - CAL_dagger_pop2,
           Estimated_difference_decomposition = sum(FCF_age_con))

## The discrepancy between difference and estimated difference is due to the continuous estimation.

###############
####### decomposition visualisation
###############
WildColors<-rev(c("#67001f","#b2182b","#d6604d","#f4a582","#f8c9b4","#fde2d2","#fde9dd","#fef4ee", ## from dark red to light red
                  "white","white",
                  "#f1f7fa","#e3eff6","#d3e7f1","#bddceb","#92c5de","#4393c3","#2166ac","#053061"))## from light blue to dark blue

levels<-c(-0.050,-0.02,-0.0100,-0.0050,-0.0010,-0.0005,-0.0001,-0.00005,-0.00001,
          0,
          0.00001,0.00005,0.0001,0.0005,0.0010,0.0050,0.0100,0.02,0.05)/10

levelss<-c(format(round(levels[1:floor(19/2)], digits = 6), scientific=F),
           "0",
           format(round(levels[(floor(19/2)+2):length(levels)], digits = 6), scientific=F))


options(scipen=18)

customAxis <- function() {
  n <- length(levels)
  y <- seq(min(levels), max(levels), length.out=n)
  rect(0, y[1:(n-1)], 1, y[2:n], col=WildColors)
  axis(4, at=y, labels=levelss)
}

Year<-seq((Int_Year-110), (Int_Year)-1,1)
Age<-seq(0, 109)

par(mar = c(5, 5, 3, 5))
filled.contour(Year,Age,FCF_age_con,levels=levels, col=WildColors,
               plot.title={par(cex.main=1.5);title(main = paste(Pop1,c(" VS. "),
                                                                Pop2, c(", "),Sex,
                                                                c(", "), Int_Year,sep=""),
                                                   xlab = "Year", ylab = "Age")},
               key.axes=customAxis(),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)





























