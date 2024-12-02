

library(circular)


# Unimodal__________________________________________________________
set.seed(100)
#simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(0), kappa=10)
simulated_data <- rvonmises(n=123, mu=circular(0), kappa=10)
data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)



tic_marks<- c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365))
tick_mark_round<- round(tic_marks,2)
names_tick<- as.character(tick_mark_round)


mean_date_PRF_sim<- mean.circular(data_na_rm_simulated)
abs_mean_date_PRF<- round(mean_date_PRF_sim,2)

mean_plus_pi<- abs_mean_date_PRF + pi

tic_marks_altered<- c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365))
ind_greter_meanPLUSpi<- which(tic_marks_altered > mean_plus_pi)
tic_marks_altered[ind_greter_meanPLUSpi]<- tic_marks[ind_greter_meanPLUSpi] - (2*pi)
tic_marks_altered_rnd<- round(tic_marks_altered,2)
names_tick_altered<- as.character(tic_marks_altered_rnd)







plot(as.circular(simulated_data), main = paste("Simulated data",sep = ""),ticks = F,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
#dev.off()


## Plot with tick marker
plot(as.circular(simulated_data), main = paste("Simulated data",sep = ""),ticks = F,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c(names_tick_altered), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)




# Plot with tick marker altered__________________________
plot(as.circular(simulated_data), main = paste("Simulated data",sep = ""),ticks = F,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c(names_tick), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)






