install.packages("tidyverse")
install.packages('pastecs')

library(tidyverse)
library(pastecs)

data = read.csv("https://raw.githubusercontent.com/andrealicialo/salphen/main/data.csv")
head(data)

#pivot the data frame into a long format
data_long = data %>% pivot_longer(cols=c('pH','Ammonium','Nitrite','Nitrate'),
                    names_to='Variable',
                    values_to='Value')
head(data_long)

# Specify the levels of the factor "variable" to pH > Ammonium > Nitrite > Nitrate. 
# This will come to use for plotting
data_long$Variable <- factor(data_long$Variable, levels = c("pH", "Ammonium", "Nitrite","Nitrate"))

round(stat.desc(data %>% select_if(is.double)),4)

install.packages("diptest")
library(diptest)

install.packages("LaplacesDemon")
library(LaplacesDemon)

install.packages("mousetrap")
library(mousetrap)

# We create a function to perform the task in whatever continous variable we indicate

normality <- function(x, title) {
    hist(x, breaks = 20, main = title)                    
    abline(v = median(x), col = "red", lwd = 3)

return(list(shapiroWilk = shapiro.test(x),
            hartigansTest = dip.test(x),
            unimodality = is.unimodal(x),
            multimodality = is.multimodal(x),
            bimodalityCoeff =bimodality_coefficient(x),
            bimodality = is.bimodal(x)))

}

# If bimodality coefficient is larger than 0.555 (stated in Pfister et al., 2013), 
# it indicates bimodality of data.

normality(data$pH, title = "pH")
normality(data$Ammonium, title = "Ammonium")
normality(data$Nitrite, title = "Nitrite")
normality(data$Nitrate, title = "Nitrate")

report = by(data %>% select_if(is.double), list(data$Day,data$Sample), stat.desc) # To summary statistics grouped by Day and Sample

report

report[['30','NiSA15']]

report['30',]

report[,'NiSA15']

boxplot =
    ggplot(data_long, aes(x=Sample, y=Value, color=Sample)) +
    geom_boxplot() +
    facet_grid(Variable ~ Day, scales = 'free_y') +
    scale_x_discrete(NULL, breaks = NULL) +
    labs(title = "Inorganic nitrogen pool and pH")
boxplot

# We create a function to calculate the mean and the standard deviation for each group

# data : a data frame
# varname : the name of a column containing the variable to be summariezed
# groupnames : vector of column names to be used as
# grouping variables

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

library(plyr); library(dplyr) # First we call this libraries
barplot_data = data_summary(data_long, varname="Value", 
                    groupnames=c("Day", "Sample","Variable"))

# Convert Day to a factor variable
barplot_data$Day=as.factor(barplot_data$Day)
barplot_data$Sample=as.factor(barplot_data$Sample)
head(barplot_data)

# library(ggplot2) 

barplot = 
    ggplot(barplot_data, aes(x=Sample, y=Value, fill=Sample)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=Value, ymax=Value+sd), width=.2,
                 position=position_dodge(.9)) +
    facet_grid(Variable ~ Day, scales = 'free_y') +
    scale_x_discrete(NULL, breaks = NULL) +
    labs(title="Inorganic nitrogen pool and pH", x="Sample", y = "Value")

barplot

lineplot_data <- barplot_data
lineplot_data$Day <- as.numeric(as.vector(lineplot_data$Day)) # Format Day column to continuos
lineplot = 
    ggplot(lineplot_data, aes(x=Day, y=Value, group = Sample)) + 
    geom_line(aes(color=Sample)) +
    geom_point(aes(color=Sample)) +
    facet_grid(Variable~., scales = 'free_y') +
    labs(title="Inorganic nitrogen pool and pH", x="Day", y = "Value")

lineplot

install.packages ("mgcv")

library(mgcv)

mod_lm2 <- gam(Nitrite ~ pH + Ammonium + Nitrate, data=data, family = quasipoisson)
summary(mod_lm2)

mod_gam2 <- gam(Nitrite ~ s(pH) + s(Ammonium) + s(Nitrate), data=data, family = quasipoisson)
summary(mod_gam2)

mod_gam3 <- gam(Nitrite ~ te(Ammonium, Nitrate), data=data, family = quasipoisson)
summary(mod_gam3)

install.packages("visreg")

library(visreg)

vis.gam(mod_gam3, type='response', plot.type='persp',
        phi=30, theta=30, n.grid=500, border=NA)
visreg2d(mod_gam3, xvar='Ammonium', yvar='Nitrate', scale='response')

kruskal.test(pH ~ paste(data$Day,data$Sample), data = data)

pw = pairwise.wilcox.test(data$pH, paste(data$Day,data$Sample),
                 p.adjust.method = "none")

# install package
install.packages("multcompView")
install.packages("rcompanion")

library(multcompView)
library(rcompanion)

pwv = pw$p.value
pwv2 = fullPTable(pwv)
pwv2

letters = multcompLetters(pwv,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)
letters

# Comment the next lines if the package is already installed
install.packages("plotly")
install.packages("ggfortify")

# Load the packages
library(plotly)
library(ggfortify)

df <- data[, c("Day","pH", "Ammonium", "Nitrite", "Nitrate")]
df <- as.matrix(df)
head(df)

# Matrix for pca at 60 days
df60 <- subset(data[, c("pH", "Ammonium", "Nitrite", "Nitrate")],data[,"Day"] == 60)
df60 <- as.matrix(df60)
head(df60)

# Dataframe to group data
data60 <- subset(data,data[,"Day"]== 60)
head(data60)

# Principal component analysis
pca <- prcomp(df60, scale. = TRUE)
# Plot
autoplot(pca, loadings = TRUE, loadings.label = TRUE,
         data = data60, colour = 'Sample', size = 6)

# install.packages("factoextra")
library(factoextra)
fviz_eig(pca)

# Access to pca results
eig.val <- get_eigenvalue(pca)
eig.val

summary(pca)

# Graph of variables
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Biplot of individuals and variables
fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF",
                col.ind = "#696969",
                alpha.var = "cos2"
                )

# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

autoplot(kmeans(df,3), data = data,
         label = TRUE, label.size = 3, frame = TRUE)

# Heatmap
autoplot(scale(df))


