#' loadUKdata function
#' 
#' function to load and interpolate the births and population of 20 UK cities.
#' puts the data into biweeks.
#' 
#' @param x, town name. Options are 
#' "Bedwellty", "Birmingham", "Bradford", "Bristol", "Cardiff","Consett","Dalton.in.Furness" ,
#' "Halesworth", "Hastings", "Hull", "Leeds", "Lees", "Liverpool", "London", 
#' "Manchester","Mold", "Northwich", "Nottingham", "Oswestry", "Sheffield"        
#' 
loadUKdata <- function(x){
  
  load('twentycities.rda')
  
  places <- names(table(measles$town))
  measles %>% 
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(town==x & year>=1944 & year<1966) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1966, select=c(time,cases)) -> dat
  demog %>% subset(town==x,select=-town) -> demog
  
  demog <- subset(demog,year >= 1944 & year<= 1966)
  
  data <- dataIP(data=dat,demog=demog,IP=2)

  return(data)

}
