# This script combines the variables from MoveBank requests into a single dataframe for each species/ time window/ prs-abs

setwd('/mnt/arctic/c/Share/kguay/hummingbird/movebank/downloaded_annotations')

for(species in c('bchu','bthu','cahu','ruhu')){
	
	for(suffix in c('fcl','min','pls')){

		if(suffix == 'fcl'){
			pt <- c('prs','abs')
		}
		else{
			pt <- 'prs'
		}
		
		for(point.type in pt){
			files <- paste0(species, '/', list.files(species, pattern=paste0('.*', point.type, '.*', suffix, '.*csv')))
			print(files)
			print(' ')
		
			# create dataframe from first file
			df <- read.csv(files[1])

			for(i in 2:length(files)){
				file <- files[i]
				df.tmp <- read.csv(file)
				df <- merge(df, df.tmp, by=c('timestamp', 'location.long', 'location.lat', 'height.above.ellipsoid'))
			}
		
			write.csv(df, file=paste0('combined/', species, '_', suffix, '_', point.type, '.csv'))
			save(df, file=paste0('combined/', species, '_', suffix, '_', point.type, '.RData'))
		}
	}
}

