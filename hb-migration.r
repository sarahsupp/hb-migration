#Code for eBird migration project
# (c) 2013 Sarah Supp 

#set working directory
wd = "C:/Users/sarah/Documents/eBird_data/"
setwd(wd)

# read in eBird data
bchu = read.table('ebd_bkchum_200801_201312_relOct-2013/ebd_bkchum_200801_201312_relOct-2013.txt',
                  header=TRUE, sep="\t", fill=TRUE, nrows=10)
ruhu = read.table('ebd_rufhum_200801_201312_relOct-2013/ebd_rufhum_200801_201312_relOct-2013.txt', 
                  header=TRUE, sep="\t", fill=TRUE, nrows=10)

# filter data 