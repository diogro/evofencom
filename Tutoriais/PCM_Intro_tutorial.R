####Tutorial: Introdução aos Métodos Comparativos Filogenéticos

install.packages("ape")
install.packages("corHMM")
install.packages("diversitree")
install.packages("geiger")
install.packages("phytools")
install.packages("devtools")
devtools::install_github("liamrevell/phytools")

packageVersion("phytools")
library(phytools)

###Vamos construir uma árvore filogenética simples:
text.string<-"(Shark,(Iguana,((Cow,Pig),Human)));"
text.string
vert.tree<-read.tree(text=text.string)

class(vert.tree)

###agora podemos visualizar a árvore
plot(vert.tree)


###podemos mudar a forma como a árvore é representada:
plot(vert.tree, direction="upwards", type="cladogram")
plot(vert.tree, type="unrooted")
plot(unroot(vert.tree))

plot(vert.tree, type="unrooted")
nodelabels()
tiplabels()

print(vert.tree)
str(vert.tree)

vert.tree$edge
plot(vert.tree)

dev.off()
plot(vert.tree)

vert.tree$edge
nodelabels()
tiplabels()

###Vamos substituir os nomes

tips<-c(
  "Charcharinus perezi",
  "Iguana iguana",
  "Bos taurus",
  "Sus scrofa",
  "Homo sapiens")
new.tree<-vert.tree
new.tree$tip.label<-tips
plot(new.tree)
plotTree(new.tree,ftype="i",type="cladogram",nodes="centered")
nodelabels()
tiplabels()


####Vamos usar uma base de dados já existente:

anolis.tree<-read.tree(
  file="http://www.phytools.org/Rbook/1/Anolis.tre")
anolis.tree

###se você tiver o arquivo .tre salvo em seu diretório, o caminho seria:
#anolis.tree<-read.tree(file="Anolis.tre")
#anolis.tree

plotTree(anolis.tree, ftype="i", fsize=0.6)

###vamos fazer um subset?
spp<-c("cooki","poncensis","gundlachi","pulchellus",
       "krugi","evermanni")
nodes<-sapply(spp,grep,x=anolis.tree$tip.label)
nodes
add.arrow(anolis.tree,tip=nodes,offset=2) ###add arrows to your phylogeny based on a list of species

###esse é um comando útil para tirar tips da árvore:
pruned_anolis.tree<-drop.tip(anolis.tree,spp)
plotTree(pruned_anolis.tree,ftype="i", fsize=0.6, type='fan')

###vamos comparar a árvore original e a modificada?
plot(cophylo(anolis.tree,pruned_anolis.tree),fsize=0.5) ###plot two phylogenies and compare them
plot(cophylo(anolis.tree,pruned_anolis.tree),fsize=0.45,link.type="curved")


###que tal, ao invés de tirar tips, fazer um subset da árvore?
kept_anolis.tree<-keep.tip(anolis.tree,spp)
kept_anolis.tree

plot(cophylo(anolis.tree,kept_anolis.tree),fsize=c(0.45,1),link.type="curved")

ls()

###dá para salvar todas elas em um único objeto:
anolis.trees<-c(anolis.tree,pruned_anolis.tree,
                kept_anolis.tree)

anolis.trees
print(anolis.trees,details=TRUE)

################################################
####using comparative data
list.files()

anole.data<-read.csv(file="http://www.phytools.org/Rbook/1/anole.data.csv",row.names=1)
head(anole.data)

ecomorph<-read.csv(file="http://www.phytools.org/Rbook/1/ecomorph.csv",row.names=1,stringsAsFactors = TRUE)
head(ecomorph,20)

dim(anole.data)
dim(ecomorph)
###dá para notar que meus dados fenotípicos e ecológios não batem

library(geiger)

###vamos ver quais espécies estão em cada dataframe, e se elas são congruentes com as espécies amostradas na árvore
name.check(anolis.tree,anole.data)

chk<-name.check(anolis.tree,ecomorph)
summary(chk)
chk

###vamos cortar a árvore e deixar apenas as espécies com ecologia?
ecomorph.tree<-drop.tip(anolis.tree,chk$tree_not_data)

name.check(ecomorph.tree,ecomorph)
name.check(ecomorph.tree,anole.data)
####temos um problema aqui: mais espécies na árvore do que na minha tabela morfológica

###podemos alterar a matriz também:
ecomorph.data<-anole.data[ecomorph.tree$tip.label,]
head(ecomorph.data)
name.check(ecomorph.tree,ecomorph.data)
#####se deu OK, é porque tudo está certinho: dados morfológicos, classificação ecológica e árvore possuem exatamente a mesma composição de espécies!!


###Vamos implementar uma análise de componentes principais filogenéticos??
ecomorph.pca<-phyl.pca(ecomorph.tree,ecomorph.data)
ecomorph.pca

dev.off()
plot(ecomorph.pca)

scores(ecomorph.pca)


###uma forma melhor de visualizar os resultados
phylomorphospace(ecomorph.tree,
                 scores(ecomorph.pca)[,1:2],
                 xlab="PC1", ylab="PC2", ftype="off",
                 node.size=c(0,1),bty="n")

###Vamos colorir os pontos de acordo com o grupo ecológico?
par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(ecomorph.tree,
                 scores(ecomorph.pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1 (overall size)",
                 ylab=expression(paste("PC2 ("%up%"lamellae number, "
                                       %down%"tail length)")))
eco<-setNames(ecomorph[,1],rownames(ecomorph))
ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[ecomorph.tree$tip.label,],cex=0.5)
legend(x="bottomright",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)

####DESAFIO 1
#Baixe os dois arquivos seguintes sobre lagartixas do gênero Phelsuma no site do livro: phel.csv (file='http://www.phytools.org/Rbook/1/phel.csv') e phel.phy (file='http://www.phytools.org/Rbook/1/phel.phy') (Harmon et al. 2010). 
#O arquivo phel.csv é um arquivo CSV que contém valores de características para dez traços morfológicos diferentes. O arquivo phel.phy é uma filogenia de trinta e três espécies. 
#Importe tanto os dados quanto a árvore a partir dos arquivos e use a função 'name.check' para identificar quaisquer diferenças entre os dois conjuntos de dados. 
#Se você encontrar diferenças, altere a filogenia e faça uma subamostragem dos dados de características para incluir apenas as espécies presentes tanto no arquivo de dados quanto na árvore. Plote a árvore.


###DESAFIO 2
#Use a função 'phyl.pca' para realizar uma análise de componentes principais filogenética (PCA Filogenética) do conjunto de dados morfológicos e da árvore do problema prático 1.
#Quando os dados de diferentes variáveis em uma PCA possuem ordens de magnitude distintas, frequentemente faz sentido transformá-los usando o logaritmo natural e conduzir a análise nos valores log-transformados em vez de usar as características originais. 
#Inspecione seus dados para ver se este é o caso e então decida se deve ou não aplicar a transformação logarítmica antes de realizar a sua PCA filogenética. 
#Depois de obter o resultado da PCA, crie um screeplot para visualizar a distribuição da variação entre os diferentes eixos de componentes principais.

###DESAFIO 3
#Use a função phylomorphospace para criar uma projeção única da filogenia no espaço morfológico (morphospace) para os dois primeiros eixos da PCA do problema prático 2. 
#Consegue pensar em uma maneira de projetar a árvore em um espaço definido por mais de duas dimensões de componentes principais? 
#Dica: consulte as páginas de ajuda das funções 'phylomorphospace3d' e 'phyloScattergram' para ter ideias, ou considere simplesmente subdividir a sua tela de plotagem usando o comando par(mfrow).

###Vamos testar correlações filogenéticas usando PGLS?

library(phytools)

primate.data<-read.csv(file="http://www.phytools.org/Rbook/3/primateEyes.csv", row.names=1, stringsAsFactors = TRUE)
head(primate.data)

primate.tree<-read.tree(file="http://www.phytools.org/Rbook/3/primateEyes.phy")
primate.tree

name.check(primate.tree, primate.data)

dev.off()
plotTree(primate.tree,fsize=0.5,ftype="i")

data(primate.tree) ####get data directly from Phytools

sigmoidPhylogram(primate.tree,fsize=0.4,ftype="i")

####PGLS
library(nlme)

##Primeiro, precisamos construir um objeto com a estrutura da correlação
spp<-rownames(primate.data)
spp

?corBrownian
corBM<-corBrownian(phy = primate.tree, form=~spp)
corBM

###Agora, podemos estabelecer nossa correlação filogenética:
pgls.primate<-gls(
  log(Orbit_area)~log(Skull_length),
  data=primate.data, correlation=corBM)
pgls.primate

###Se quisermos incluir informações discretas, por exemplo, se a área orbital está associada ao comprimento do crânio,
#dependendo do tipo de padrão de atividade, podemos fazer uma ANCOVA filogenética:
head(primate.data)
primate.data$Activity_pattern

primate_ancova<-gls(
  log(Orbit_area)~log(Skull_length)+Activity_pattern,
  data=primate.data, correlation=corBM)
primate_ancova

summary(primate_ancova)
anova(primate_ancova)

#Podemos plotar esses resultados graficamente
## set the margins of our plot using par
par(mar=c(5.1,5.1,2.1,2.1))
## set the point colors for the different levels
## of our factor
pt.cols<-setNames(c("#87CEEB","#FAC358","black"),
                  levels(primate.data$Activity_pattern))
## plot the data
plot(Orbit_area~Skull_length,data=primate.data,pch=21,
     bg=pt.cols[primate.data$Activity_pattern],
     log="xy",bty="n",xlab="skull length (cm)",
     ylab=expression(paste("orbit area (",mm^2,")")),
     cex=1.2,cex.axis=0.7,cex.lab=0.8)
## add a legend
legend("bottomright",names(pt.cols),pch=21,pt.cex=1.2,
       pt.bg=pt.cols,cex=0.8)
## create a common set of x values to plot our
## different lines for each level of the factor
xx<-seq(min(primate.data$Skull_length),
        max(primate.data$Skull_length),length.out=100)
## add lines for each level of the factor
lines(xx,exp(predict(primate_ancova,
                     newdata=data.frame(Skull_length=xx,
                                        Activity_pattern=as.factor(rep("Cathemeral",100))))),
      lwd=2,col=pt.cols["Cathemeral"])
lines(xx,exp(predict(primate_ancova,
                     newdata=data.frame(Skull_length=xx,
                                        Activity_pattern=as.factor(rep("Diurnal",100))))),
      lwd=2,col=pt.cols["Diurnal"])
in the slope of the relationship between the two continuous variables as a function of the factor: a common
lines(xx,exp(predict(primate_ancova,
                     newdata=data.frame(Skull_length=xx,
                                        Activity_pattern=as.factor(rep("Nocturnal",100))))),
      lwd=2,col=pt.cols["Nocturnal"])


######Vamos modelar a evolução de caracteres contínuos na árvore?

###O modelo Browniano: vamos modelar um BM?
#set values for time steps and sigma squared parameter
t<-0:100
sig2<-0.01
## simulate a set of random changes
x<-rnorm(n=length(t)-1,sd=sqrt(sig2))
## compute their cumulative sum
x<-c(0,cumsum(x))
# create a plot with nice margins
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,x,type="l",ylim=c(-2,2),bty="n",
     xlab="time",ylab="trait value",las=1,
     cex.axis=0.8)
# set number of simulations
nsim<-100
# create matrix of random normal deviates
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2)),
          nsim,length(t)-1)
# calculate the cumulative sum of these deviates
# this is now a simulation of Brownian motion
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
# plot the first one
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,X[1,],ylim=c(-2,2),type="l",bty="n",
     xlab="time",ylab="trait value",las=1,
     cex.axis=0.8)
# plot the rest
invisible(apply(X[2:nsim,],1,function(x,t) lines(t,x),
                t=t))

#Vamos avaliar a dependência do BM do parâmetro sigma2?

#create matrix of random normal deviates
# but with a smaller sd
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2/10)),
          nsim,length(t)-1)
# calculate the cumulative sum of these changes
# this is now a simulation of Brownian motion
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
# plot as above
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,X[1,],ylim=c(-2,2),type="l",bty="n",xlab="time",ylab="trait value",las=1,
     cex.axis=0.8)
invisible(apply(X[2:nsim,],1,function(x,t) lines(t,x),
                t=t))
# calculate variance of columns
v<-apply(X,2,var)
# plot the results
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,v,ylim=c(0,0.1),type="l",xlab="time",
     ylab="variance",bty="n",las=1,
     cex.axis=0.8)
lines(t,t*sig2/10,lwd=3,col=rgb(0,0,0,0.1))
legend("topleft",c("observed variance","expected variance"),
       lwd=c(1,3),col=c("black",rgb(0,0,0,0.1)),
       bty="n",cex=0.8)

###Agora vamos simular a evolução de um atributo contínuo considerando o modelo BM?

bacteria.data<-read.csv(file="http://www.phytools.org/Rbook/4/bac_rates.csv",row.names=1, stringsAsFactors = TRUE)
bacteria.data

bacteria.tree<-read.tree(file="http://www.phytools.org/Rbook/4/bac_rates.phy")
bacteria.tree

## graph phylogeny using plotTree
plotTree(bacteria.tree,ftype="i",fsize=0.5,
         lwd=1,mar=c(2.1,2.1,0.1,1.1))
## add a horizontal axis to our plot
axis(1,at=seq(0,1,length.out=5),cex.axis=0.8)

genome_size<-setNames(bacteria.data$Genome_Size_Mb,rownames(bacteria.data))
genome_size

names(genome_size)<-rownames(bacteria.data)
head(genome_size)

library(geiger)

name.check(bacteria.tree,genome_size)

?fitContinuous
## fit Brownian motion model using fitContinuous
fitBM_gs<-fitContinuous(bacteria.tree,genome_size)
fitBM_gs

## pull our mutation accumulation rate as a named vector
mutation<-setNames(bacteria.data[,"Accumulation_Rate"],
                   rownames(bacteria.data))
head(mutation)

###Se a gente plotar as taxas de mutação, veremos que sua distribuição não é normal.
#Vamos logaritmizar esse atributo, e comparar os dados antes e após a transformação:
## set up for side-by-side plots
par(mfrow=c(1,2),mar=c(6.1,4.1,2.1,1.1))
## histogram of mutation accumulation rates on original scale
hist(mutation,main="",las=2,xlab="",
     cex.axis=0.7,cex.lab=0.9,
     breaks=seq(min(mutation),max(mutation),
                length.out=12))
mtext("(a)",adj=0,line=1)
mtext("rate",side=1,line=4,cex=0.9)
## histogram of mutation accumulation rates on log scale
ln_mutation<-log(mutation)
hist(ln_mutation,main="",las=2,xlab="",
     cex.axis=0.7,cex.lab=0.9,
     breaks=seq(min(ln_mutation),max(ln_mutation),
                length.out=12))
mtext("(b)",adj=0,line=1)
mtext("ln(rate)",side=1,line=4,cex=0.9)

###
## fit Brownian motion model to log(mutation accumulation)
fitBM_ar<-fitContinuous(bacteria.tree,ln_mutation)
fitBM_ar


###Vamos calcular o sinal filogenético?
#Primeiro, Blomberg's K
phylosig(bacteria.tree,genome_size)

## test for significant phylogenetic signal using
## Blomberg’s K
K_gs<-phylosig(bacteria.tree,genome_size,
               test=TRUE,nsim=10000)
K_gs

## test for phylogenetic signal in mutation accumulation
## rate
K_ar<-phylosig(bacteria.tree,ln_mutation,
               test=TRUE,nsim=10000)
K_ar

## plot null-distribution and observed value of K
plot(K_gs,las=1,cex.axis=0.9)
plot(K_ar,las=1,cex.axis=0.9)


###Agora, vamos estimar o sinal filogenético considerando o lambda de Pagel:
lambda_gs<-phylosig(bacteria.tree,genome_size,method="lambda",test=TRUE)
lambda_gs

lambda_ar<-phylosig(bacteria.tree, ln_mutation,
                    method="lambda",test=TRUE)
lambda_ar

par(mfrow=c(2,1),mar=c(5.1,4.1,2.1,2.1),
    cex=0.8)
## plot the likelihood surfaces of lambda for each of our
## two traits
plot(lambda_gs,las=1,cex.axis=0.9,bty="n",
     xlim=c(0,1.1))
mtext("(a)",line=1,adj=0)
plot(lambda_ar,las=1,cex.axis=0.9,bty="n",
     xlim=c(0,1.1))
mtext("(b)",line=1,adj=0)