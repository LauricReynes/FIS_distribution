# Distribution of FIS into discrete classes

#### Download R packages
```{r}
library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(ggplot2)
```
#### 1. Load a VCF file (and corresponding individual strata file) and convert to GENIND object
```{r}
setwd("file path")

vcf <- read.vcf("NAME.vcf", to = 130572)
X<- loci2genind(vcf)

# The strata file with two columns : 'INDIVIDUAL' 'STRATA'
strata <- read.table(file = "NAME.txt", header = TRUE)
strata(X) <- strata
setPop(X) <- ~STRATA

  # Rename strata
  POPNames <- c("POP1","POP2","POP3","POP4","...")
  popNames(X) <- POPNames
```
#### 2. Compute FIS across loci from GENIND object (per individual strata) and extract values to a newly dataframe 'FIS'
```{r}
bc <-basic.stats(X)
FIS <-data.frame(bc$Fis)
```
#### 3. Divide the distribution of FIS into 10 intervals and add values in the datadrame 'FIS'
```{r}
FIS$POP.1 <- cut(FIS$POP1, c(-1.1,-0.8001,-0.5001,-0.3001,-0.1001,-0.0001,0.0999,0.2999,0.4999,0.7999,1.1),include.lowest=TRUE,
                    labels = c(-1,-0.8,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.8,1))

# 4. Compute the abundance of each class and save values in the dataframe 'POP1'

A <-length(which(FIS$POP.1 == -1)); B <-length(which(FIS$POP.1 == -0.8))  
C <-length(which(FIS$POP.1 == -0.5)); D <-length(which(FIS$POP.1 == -0.3))  
E <-length(which(FIS$POP.1 == -0.1)); F <-length(which(FIS$POP.1 == 0.1))  
G <-length(which(FIS$POP.1 == 0.3)); H <-length(which(FIS$POP.1 == 0.5))  
I <-length(which(FIS$POP.1 == 0.8)); J <-length(which(FIS$POP.1 == 1))  

cut <- c("[-1.0,-0.8[","[-0.8,-0.5[","[-0.5,-0.3[","[-0.3,-0.1[","[-0.1,0.0[","[0.0,0.1[","[0.1,0.3[","[0.3,0.5[","[0.5,0.8[","[0.8,1.0]")
values <- c(A,B,C,D,E,F,G,H,I,J)  

POP1 <- data.frame(cut,values)
```
#### 5. Plot each class as a histogram and save in image format (300 dpi)
```
ggplot(POP1, aes(x = cut, y =(values/length(FIS$POP.1[!is.na(FIS$POP.1)])))) +
  geom_bar(stat="identity",color="black", fill="firebrick1",alpha =0.6)+
  scale_x_discrete(limits=c("[-1.0,-0.8[","[-0.8,-0.5[","[-0.5,-0.3[","[-0.3,-0.1[","[-0.1,0.0[","[0.0,0.1[","[0.1,0.3[","[0.3,0.5[","[0.5,0.8[","[0.8,1.0]"),name = 'FIS')+
  scale_y_continuous(name = 'Frequency of loci',limits= c(0,0.60))+
  theme_bw() +
  geom_text(aes(x="[-0.8,-0.5[", label="POP1", y=0.50), colour="black", angle=0,position = position_dodge(width=0.9),  size=7)+
  geom_text(aes(label=values), vjust=-1, color="black",position = position_dodge(width=0.9),  size=8)+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        text=element_text(family="Tahoma"),
        axis.title = element_text(size = 18),
        axis.text.x=element_text(colour="black", size = 18,angle =0),
        axis.text.y=element_text(colour="black", size = 18))
        
 ggsave("Distribution_FIS_POP1.jpg",width=75 ,height=45,dpi=300,units="cm")
 ```
#### 6. Apply steps 3-5 to other strata of the dataframe 'FIS' to produce similar distributions for other populations
