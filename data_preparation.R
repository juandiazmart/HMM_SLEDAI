library(haven)
library(tidyverse)
library(lubridate)

# Import data

sledai=read_sas("sledai_i.sas7bdat", 
                NULL)
therapy=read_sas("therapy1.sas7bdat", 
                 NULL)
slicc=read_sas("slicc_n.sas7bdat", 
               NULL)

demo=read_sas("demo.sas7bdat", 
              NULL)

# Prepare data

therapy=therapy %>% select(ID,STERDOSE)
slicc=slicc %>% select(ID,score_n)
demo= demo %>% select(PTNO,SEX,RACE,dateofdx,deathdt,birthdt) %>%  mutate(RACE=case_when(
  RACE=="1" | RACE=="2" | RACE=="3"~ RACE,
  TRUE ~ "4"
))

mort=sledai %>% select(ID,ptno,sledai2_i_miss,assdt) %>% left_join(therapy,by="ID") %>% left_join(slicc,by="ID") %>%  
  mutate(STERDOSE=case_when(
    is.na(STERDOSE) ~ 0,
    STERDOSE>60 ~ 60,
    TRUE ~ STERDOSE
  ),sledai_2k=round(sledai2_i_miss + (0.32 *STERDOSE - 0.0031 *STERDOSE*STERDOSE ))) %>% left_join(demo,by=c("ptno"="PTNO")) %>%
  mutate(ageatass=as.numeric((assdt-birthdt)/365),disdur=as.numeric((assdt-dateofdx)/365)) %>%
  filter(!is.na(deathdt)) %>% group_by(ptno) %>% filter(assdt==max(assdt)) %>% ungroup() %>%
  mutate(assdt=deathdt,ID=ID+1,sledai_2k=9999)

df=sledai %>% select(ID,ptno,sledai2_i_miss,assdt) %>% left_join(therapy,by="ID") %>% left_join(slicc,by="ID") %>%  
  mutate(STERDOSE=case_when(
    is.na(STERDOSE) ~ 0,
    STERDOSE>60 ~ 60,
    TRUE ~ STERDOSE
  ),sledai_2k=round(sledai2_i_miss + (0.32 *STERDOSE - 0.0031 *STERDOSE*STERDOSE ))) %>% left_join(demo,by=c("ptno"="PTNO")) %>%
  mutate(ageatass=as.numeric((assdt-birthdt)/365),disdur=as.numeric((assdt-dateofdx)/365)) %>% 
  bind_rows(mort) %>% mutate(time=as.numeric((assdt-dateofdx)/365)) %>% arrange(ID) %>% group_by(ptno) %>% 
  filter(ptno!=2046 & ptno!=51 & ptno!=303 & ptno!=731,n()>1) %>% ungroup()

df %>% filter(ptno>900 & ptno<920,sledai_2k!=9999 ) %>% ggplot(aes(x=sledai_2k)) + 
  geom_density() + facet_wrap(~ptno)

# HMM fit

library(msm)

q_matrix=rbind(state1=c(0, 0.20, 0, 0.006),
               state2=c(0.37, 0, 0.3, 0.02),
               state3=c(0, 1, 0, 0.6),
               state4=c(0, 0, 0, 0))

model=msm(sledai_2k ~ time, subject=ptno, data=df, qmatrix=q_matrix,fixedpars = F,
          control=list(fnscale=300000,reltol = 1e-16, maxit=10000),
          hmodel = list(hmmPois(5),
                        hmmPois(10),
                        hmmPois(20),
                        hmmIdent(9999)))


pmatrix.msm(model)

vit=viterbi.msm(model)
vit$fitted=as.factor(vit$fitted)

plot=df %>% inner_join(vit,by = c("ptno"="subject","time"="time")) %>% 
  select(ptno,time,sledai_2k,observed,fitted)

plot %>% filter(ptno %in% c(904,905,906,916,909,1,7,9,26),sledai_2k!=9999) %>% ggplot(aes(x=time,y=sledai_2k)) + 
  geom_line(size=.3)+
  geom_point(aes(y=observed,color=fitted),size=.8) + facet_wrap(~ptno) +
  labs(x="Time (years since diagnosis)",y="SLEDAI-2K",color="Latent state")+ theme_classic()+
  scale_color_manual(values = c("#FF007F","#0066CC","#3CB043")) + 
  scale_color_discrete(labels=c("Mild","Moderate", "Severe"))
