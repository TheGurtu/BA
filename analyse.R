#install.packages(" ggplot2 ")
library (ggplot2)
library(dplyr)
require(gridExtra)
data <- read.csv("D:\\Desktop\\Bachelorarbeit\\Code\\model.csv", header=TRUE)
print(data$time)

scenario_1 <- data.frame(subset(data, scenario == "influence_of_rp"))#data[data$scenario == 'mintest']
data_list <- split(scenario_1,scenario_1$ident)
par(mfrow = c(4, 4))
sir_plot(data_list)

scenario_3 <- data.frame(subset(data, scenario == "compare_different_constant_rp"))
data_list <- split(scenario_2,scenario_1$ident)
par(mfrow = c(4, 4))

scenario_3 <- data.frame(subset(data, scenario == "compare_different_dist_rp_no_threshold"))
data_list <- split(scenario_3,scenario_1$ident)
par(mfrow = c(4, 4))

scenario_4 <- data.frame(subset(data, scenario == "compare_different_dist_rp_with_threshold"))
data_list <- split(scenario_4,scenario_1$ident)
par(mfrow = c(4, 4))







data_list <- split(scenario_1,scenario_1$ident)
print(data_list[[1]])
#par(mfrow = c(1, 1))

plot(data_list[[1]]$time, data_list[[1]]$extract_sus, , ylim=c(0,1000))
model <- lm(extract_sus~time, data = data_list[[1]])
abline(model, col="blue")
summary(model)

ggplot(model)

df <- data.frame()
figure <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
par(mfrow = c(4, 4))
sir_plot <- function(data_list){for(ident in 1:length(data_list)){
  sus <- data_list[[ident]] %>% select(time, extract_sus) %>% group_by(time) %>% summarise(max_sus = max(extract_sus),min_sus = min(extract_sus), mean_sus = mean(extract_sus)) %>% arrange(time)
  #print(grouped)#data_list[[ident]] %>% select(ensemble,extract_inf)
  plot(sus$time, sus$max_sus/data_list[[ident]]$extract_numAgents[1], type="l", main = sprintf("Max recovered %s",ident), col="blue", ylim=c(0,1))#data_list[[ident]]$extract_numAgents[1]
  lines(sus$time, sus$min_sus/data_list[[ident]]$extract_numAgents[1], type="l",  col="blue")
  lines(sus$time, sus$mean_sus/data_list[[ident]]$extract_numAgents[1], type="l",  col="blue")
 
  
  
           
  
  infected <- data_list[[ident]] %>% select(time, extract_inf) %>% group_by(time) %>% summarise(max_inf = max(extract_inf),min_inf = min(extract_inf), mean_inf = mean(extract_inf)) %>% arrange(time)
  #print(grouped)#data_list[[ident]] %>% select(ensemble,extract_inf)
  lines(infected$time, infected$max_inf/data_list[[ident]]$extract_numAgents[1], type="l", main = sprintf("MAx Infizierte %s",ident), col="red")
  lines(infected$time, infected$min_inf/data_list[[ident]]$extract_numAgents[1], type="l",col="red")
  lines(infected$time, infected$mean_inf/data_list[[ident]]$extract_numAgents[1], type="l", col="red")
  
  recovered <- data_list[[ident]] %>% select(time, extract_re) %>% group_by(time) %>% summarise(max_re = max(extract_re),min_re = min(extract_re), mean_re = mean(extract_re)) %>% arrange(time)
  #print(grouped)#data_list[[ident]] %>% select(ensemble,extract_inf)
  lines(recovered$time, recovered$max_re/data_list[[ident]]$extract_numAgents[1], type="l", main = sprintf("Max recovered %s",ident), col="green")
  lines(recovered$time, recovered$min_re/data_list[[ident]]$extract_numAgents[1], type="l", col="green")
  lines(recovered$time, recovered$mean_re/data_list[[ident]]$extract_numAgents[1], type="l", col="green")
  # test <- print(ggplot(sus, aes(x=time)) +  geom_line(aes(y=mean_sus/data_list[[ident]]$extract_numAgents[1]), color="blue") + geom_ribbon(aes(ymin=min_sus/data_list[[ident]]$extract_numAgents[1], ymax=max_sus/data_list[[ident]]$extract_numAgents[1], x=time),fill="blue", alpha = 0.3))
   #test1<- print(ggplot(infected, aes(x=time)) +  geom_line(aes(y=mean_inf/data_list[[ident]]$extract_numAgents[1]),color="red") + geom_ribbon(aes(ymin=min_inf/data_list[[ident]]$extract_numAgents[1], ymax=max_inf/data_list[[ident]]$extract_numAgents[1], x=time), fill="red", alpha = 0.3))
  # test2 <-print(ggplot(recovered, aes(x=time)) +  geom_line(aes(y=mean_re/data_list[[ident]]$extract_numAgents[1]),color="green") + geom_ribbon(aes(ymin=min_re/data_list[[ident]]$extract_numAgents[1], ymax=max_re/data_list[[ident]]$extract_numAgents[1], x=time), fill = "green", alpha = 0.3))
  #print(test)
  #grid.arrange(test, test1,test2, nrow=3)
  #plot(data_list[[ident]]$time, data_list[[ident]]$extract_inf,  main = sprintf("Infizierte %s",ident))
  #plot(data_list[[ident]]$time, data_list[[ident]]$extract_sus,  main =  sprintf("Sus %s",ident))
  #plot(data_list[[ident]]$time, data_list[[ident]]$extract_re,  main =  sprintf("Recovered %s",ident))
}
}
par(mfrow = c(1, 1))
for(ident in 1:length(data_list)){
  setting_data <- data_list[[ident]]
  silmulations <- split(setting_data,setting_data$ensemble)
  durations <- list()
  print(durations)
  for(simu in 1:length(silmulations)){
    simulation <- silmulations[[simu]]
    zero_index <- match(0,simulation$extract_inf)
    #print(zero_index)
    #print(simulation[zero_index-1,])
    #print(simulation[zero_index,])
    duration <- simulation[zero_index,"time"]
    durations <- append(durations, simulation[zero_index,"time"])
  }
  durations <- unlist(durations, use.names = FALSE)
  boxplot(durations, xlab=sprintf(" %s",ident), ylim=c(0,60))
  durations <- data.frame(durations)
  print(ggplot(durations, aes(y=durations)) + geom_boxplot())
}

par(mfrow = c(4, 4))
for(ident in 1:length(data_list)){
  setting_data <- data_list[[ident]]
  silmulations <- split(setting_data,setting_data$ensemble)
  #print(silmulations)
  maximums <- list()
  #print(durations)
  for(simu in 1:length(silmulations)){
    simulation <- silmulations[[simu]]
    max <- max(simulation$extract_inf)
    print(max)
    #print(simulation[zero_index-1,])
    #print(simulation[zero_index,])
    #duration <- simulation[zero_index,"time"]
    maximums <- append(maximums,max)
  }
  #print(maximums)
  maximums <- unlist(maximums, use.names = FALSE)
  print(maximums)
  min_p = min(maximums)
  max_p = max(maximums)
  boxplot(maximums, xlab=sprintf(" %s",ident), ylim=c(min_p,max_p))
  maximums <- data.frame(maximums)
  print(ggplot(maximums, aes(y=maximums)) + geom_boxplot())
}
#test <- data.frame(data_list[1])
#test_list <- split(test, test$X1.ensemble)
#miau <- data.frame(test_list[2])
#plot(data_list[[1]]$time, data_list[[1]]$extract_inf )
#plot(data_list[[1]]$time, data_list[[1]]$extract_sus)
#plot(data_list[[1]]$time, data_list[[1]]$extract_re)
#plot(data_list[1]$time, data_list[1]$extract_inf)

