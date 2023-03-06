#saves simulation in the right place
save_sim <- function(sim,path_to_project_directory){
  sim_name <- deparse(substitute(sim))
  saveRDS(sim,file=paste0(path_to_project_directory,"/output/simulations/",sim_name,".rds"))
}

#read sim saved using the above function
read_sim <- function(sim_name,path_to_project_directory){
  rds <- readRDS(paste0(path_to_project_directory,"/output/simulations/",sim_name,".rds"))
  return(rds)
}

#returns an empty df with only column names
empty_df <- function(...){
  cols <- c(...)
  df <- data.frame(matrix(nrow=0,ncol=length(cols)))
  colnames(df) <- cols
  return(df)
}
#Returns a single sample from a normal distribution given a mean and SD
rnorm1 <- function(mean,sd){
  args <- list(1,mean,sd)
  samp <- do.call(rnorm,args)
  return(samp)
}

#Updates an expected value using the delta rule
delta_update <- function(ev,alpha,outcome){
  ev_new <- ev + alpha*(outcome - ev)
  return(ev_new)
}

#Creates a grid of learning curves for the simulate_s222 function
#Takes plot_df as input: a df containing the columns necessary to create the plots
lc_grid <- function(df){
  #get the learning curve for dQ, uncorrected by SD
  dQ_abs <- ggplot(df,aes(x=1:nrow(df))) + 
              geom_line(aes(y=dQ_mean),color="blue") +
              geom_ribbon(aes(ymin = dQ_mean-2*dQ_sd,
                              ymax = dQ_mean+2*dQ_sd),
                              alpha=.3,
                              fill="dodgerblue1")
              
  #get the learning curve for dA, uncorrected by SD
  dA_abs <- ggplot(df,aes(x=1:nrow(df))) + 
              geom_line(aes(y=dA_mean),color="red") +
              geom_ribbon(aes(ymin = dA_mean-2*dA_sd,
                              ymax = dA_mean+2*dA_sd),
                          alpha=.3,
                          fill="salmon3")
  #get the combined learning curve, with a y-axis corresponding to z-score
  combined <- ggplot(df,aes(x=1:nrow(df))) + 
                geom_line(aes(y=dQ_mean_z),color="blue") +
                geom_ribbon(aes(ymin = dQ_mean_z-2*dQ_sd_z,
                                ymax = dQ_mean_z+2*dQ_sd_z),
                            alpha=.3,
                            fill="dodgerblue1") +
                geom_line(aes(y=dA_mean_z),color="red") +
                geom_ribbon(aes(ymin = dA_mean_z-2*dA_sd_z,
                                ymax = dA_mean_z+2*dA_sd_z),
                            alpha=.3,
                            fill="salmon3")
    
    
  #arrange these plots in a grid, and return it
  grid <- plot_grid(plotlist=list(dQ_abs,dA_abs,combined)) 
  return(grid)
}