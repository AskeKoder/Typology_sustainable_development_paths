# =============================================================================
#  Coloniser identity across cluster configurations - ordinal (proportional-odds)
#
#  Coloniser identity costs 3 parameters here vs 27 in a 10-cluster multinomial.
#  On N = 55 that is the difference between detecting the effect and not.
#
#  Run from the repository root:  source("main_scripts/4.3.4_ColonizerTests.R")
#  Requires: MASS
# =============================================================================

suppressMessages(library(MASS))

folder<-"Data_Colonial/"

ROOT <- if (dir.exists("Data_Colonial")) "." else ".."
f    <- paste0(ROOT, "/Data_Colonial/")
folder <- f
sm   <- read.csv(paste0(f, "SettlerMortality.csv"), sep=";", dec=",")

sm<-sm[!duplicated(sm$ccodealp),c("ccodealp","ajr_settmort")]; names(sm)<-c("iso3","settmort")
cr<-read.csv(paste0(folder,"Colonies.csv")); names(cr)[4]<-"power"
dropc<-c("z. Multiple colonizers","zz. Colonizer","zzz. Not colonized","zzzz. No longer colonized")
sub<-cr[!cr$power %in% dropc,]; cnt<-aggregate(list(n=sub$Year),by=sub[,c("Code","power")],FUN=length)
cnt<-cnt[order(-cnt$n),]; col<-cnt[!duplicated(cnt$Code),c("Code","power")]; names(col)<-c("iso3","colonizer")
pe<-read.csv(paste0(folder,"Prienr.csv")); names(pe)[4]<-"prienr1900"
pe<-pe[pe$Year==1900 & !pe$Code %in% c("","OWID_WRL") & !is.na(pe$Code),c("Code","prienr1900")]; names(pe)[1]<-"iso3"
full<-merge(merge(sm,col,by="iso3",all=TRUE),pe,by="iso3",all=TRUE)
tb<-table(full$colonizer); full$colonizer[full$colonizer %in% names(tb)[tb<10]]<-"Other"

rds <- paste0(ROOT, "/main_scripts/4_RankedClusters.RDS")

cl<-readRDS(rds); names(cl)[names(cl)=="SPI_countrycode"]<-"iso3"
eds <- paste0(ROOT, "/main_scripts/extendedDataScaled.csv")
reg<-read.csv(eds)[,c("SPI_countrycode","Region")]; reg<-reg[!duplicated(reg$SPI_countrycode),]; names(reg)<-c("iso3","Region")
full<-merge(full,reg,by="iso3",all.x=TRUE)
CONF<-c("Baseline","BHN_FWB","SPI_Extended","DLS_dim_only","DLS_dim_only_limExt","DLSFew_coverage",
        "DLSFew_coverage_90","DLSFew_bestAlignment","DLSFew_bestAlignment_90","DLS_dim_only_BHN_FWB")

O<-function(f,x) MASS::polr(f,data=x,Hess=TRUE)
pv<-function(f,r,df) pchisq(deviance(r)-deviance(f),df,lower.tail=FALSE)
st<-function(p) if(is.na(p)) "   n/a  " else sprintf("%7.4f%s",p,ifelse(p<0.01,"**",ifelse(p<0.05,"* ",ifelse(p<0.10,". ","  "))))
cat(sprintf("%-24s %3s %3s | %-10s %-10s %-10s %-10s | %8s\n","configuration","N","k",
  "alone","+sett+pri","excl.neoEu","+RegionFE","dAIC"))
cat(strrep("-",96),"\n")

for(cf in CONF){
  d<-merge(full,cl[,c("iso3",cf)],by="iso3"); names(d)[ncol(d)]<-"cluster"
  d<-d[complete.cases(d[,c("settmort","prienr1900","cluster","colonizer")]),]
  d$colonizer<-relevel(factor(d$colonizer),"United Kingdom"); d$ord<-factor(d$cluster,ordered=TRUE)
  if(nlevels(d$ord)<3) next
  g<-function(e) tryCatch(e,error=function(x) NA)
  p1<-g(pv(O(ord~colonizer,d),O(ord~1,d),3))
  m2<-g(O(ord~settmort+prienr1900+colonizer,d)); m2r<-g(O(ord~settmort+prienr1900,d))
  p2<-g(pv(m2,m2r,3)); dA<-g(AIC(m2)-AIC(m2r))
  dn<-d[!d$iso3%in%c("AUS","CAN","NZL","USA"),]; dn$ord<-factor(dn$cluster,ordered=TRUE); dn$colonizer<-droplevels(dn$colonizer)
  p3<-g(pv(O(ord~settmort+prienr1900+colonizer,dn),O(ord~settmort+prienr1900,dn),nlevels(dn$colonizer)-1))
  dr<-d[!is.na(d$Region),]; dr$ord<-factor(dr$cluster,ordered=TRUE)
  p4<-g(pv(O(ord~settmort+prienr1900+factor(Region)+colonizer,dr),O(ord~settmort+prienr1900+factor(Region),dr),3))
  cat(sprintf("%-24s %3d %3d | %-10s %-10s %-10s %-10s | %+8.1f\n",cf,nrow(d),nlevels(d$ord),
    st(p1),st(p2),st(p3),st(p4), dA))
}
cat("\ndAIC = AIC(with colonizer) - AIC(without); negative favours INCLUDING colonizer.\n")
