# Define funcitons to use in hb-migration code

DateConvert = function(date){
  #convert eBird column OBSERVATION.DATE into year, month and day columns
  
  year <-sapply(date,function(x){
    as.numeric(substring(x, 1, 4))
  })
  
  month <- sapply(date,function(x){
    as.numeric(substring(x, 6, 7))
  }) 
  
  day <- sapply(date,function(x){
    as.numeric(substring(x, 9, 10))
  })
  
  return (cbind(year,month,day))
}

PlotRecords = function(yeardata, species) {
  #plots the number of records across the years
  
  t = as.data.frame(table(yeardata))
  names(t) = c('year', 'count')
 
  barplot = ggplot(data = t, aes(year, count)) + geom_bar(fill="cadetblue2") + theme_bw() + 
    ggtitle(species)
  
  print(barplot)
}
