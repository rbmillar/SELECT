#Playing with pivot_wider. Could be a better version of SELECT_FORMAT
#DF=Df #[,c("TowID","Species","SF","Gear","n","lgth")]
DF=Df %>% group_by(across(c(-n))) %>% summarize(n=sum(n))
xx=pivot_wider(DF,names_from=Gear,names_prefix="n",values_from=n,values_fill=0)
xx$Gear=DF$Gear
xx


#=============Permuting=======================
permCols=function(col.names=c("n1","n2")) {
  M[,col.names]=M[,sample(colnames)]
  M }
permCols=function(n1,n2) {

  if(0.5<runif(1))
  n2new=n1
  n1new=n2
  data.frame(n1=n1new,n2=n2new) }

col.names=c("n1","n2")

grp=c(rep("A",7),rep("B",5))
haul=c(1,1,1,2,2,3,3,4,4,4,5,5)
nhauls=length(unique(haul))
lgth=c(1,2,3,1,2,1,2,1,2,3,1,2)
n1=c(3,6,1,4,5,2,1,0,0,0,0,0)
n2=c(0,0,0,0,0,0,0,4,5,6,1,6)
q1=rep(c(0.5,1,0.6),c(5,2,5))
q2=rep(c(1,1,1),c(5,2,5))
DF=data.frame(grp,haul,lgth,n1,n2,q1,q2)
DF2=rbind(DF,transform(DF,haul=haul+5)); nhauls2=2*nhauls
DF2$block=rep(c("A","B"),c(12,12))


# ***Playing with unpaired permuting of SELECT format data***
df=DF
DF
haulgrp= df |> group_by(haul) |> summarize (haulgrp=unique(grp)) |>
  slice_sample(n=nhauls) |> pull(haulgrp)
haulgrp
permuted.grp=haulgrp[df$haul]
permuted.obs=(df$grp!=permuted.grp)
df[permuted.obs,col.names]=df[permuted.obs,rev(col.names)]
df

#  ***With block***
df=DF2
#DF2
haulgrp= df |> group_by(block,haul) |> summarize (haulgrp=unique(grp)) |>
               slice_sample(n=nhauls2) |> pull(haulgrp)
haulgrp
permuted.grp=haulgrp[df$haul]
permuted.obs=(df$grp!=permuted.grp)
df[permuted.obs,col.names]=df[permuted.obs,rev(col.names)]
df

# ***Paired permuting of SELECT format data***
df=DF
permdf = df |> group_by(haul) |> do(permCols(.))
permdf

df=DF
testy=function(n1,n2) n1+n2

permdf = df |> group_by(haul) |> mutate(permute=(0.5<runif(1)) ) |>
               mutate(n1=ifelse(permute,n2,n1),n2=ifelse(permute,n1,n2))
permdf


permdf = df |> group_by(haul) |> mutate( ifelse(runif(1)<0.5,n1,n2))
permdf

permdf = df |> group_by(haul) |> reframe( permCols(n1,n2))
permdf

