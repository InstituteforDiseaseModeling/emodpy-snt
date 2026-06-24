# DHS/MIS output for Sierra Leone - save results in SierraLeone_hbhi/explore_DHS folder DHS_*year*_files_recodes_for_sims (see existing format) - this will be read in in DHS_data_extraction2.R

library(foreign)
library(haven)
library(dplyr)
# library(rgdal)
library(raster)
library(sp)

library(readr)



country = 'NGA'  #'SLE'  # 'BDI'
year_index = 6

if(country =='NGA'){
  user = Sys.getenv("USERNAME")
  user_dir = file.path("C:/Users",user)
  base_filepath = paste0(user_dir, '/Gates Foundation Dropbox/Malaria Team Folder/data/Nigeria/Nigeria_SNT_shared_data')
  # base_filepath = 'C:/Users/moniqueam/Dropbox (IDM)/NU_collaboration'
  # base_hbhi_filepath = paste0(base_filepath, '/hbhi_nigeria')
  if(year_index <=6){
    dta_dir = paste0(base_filepath, '/nigeria_dhs/data_analysis/data/DHS/Downloads')
    dhs_years = c(2010, 2013, 2015, 2018, 2021 ,2024)
    dhs_year_dirs = c('NG_2010_MIS_06192019', 'NG_2013_DHS_06192019', 'NG_2015_MIS_06192019', 'NG_2018_DHS_11072019_1720_86355', 'NG_2021_MIS_12062022_90_72922', 'NG_2024_DHS_10212025_158_158063')
    cur_year = dhs_years[year_index]
    cur_year_dir = dhs_year_dirs[year_index]
    
    # get list of relevant dta files
    dta_filepaths = list.files(path=paste0(dta_dir, '/', cur_year_dir), pattern="\\.(DTA|dta)$", full.names=TRUE, recursive=TRUE)
    dta_filepaths = dta_filepaths[substr(dta_filepaths, nchar(dta_filepaths)-3+1, nchar(dta_filepaths)) %in% c('DTA', 'dta')]
    dta_list = list()
    for(dd in 1:length(dta_filepaths)){
      dta_cur = read.dta(dta_filepaths[dd])
      dta_list[[dd]] = dta_cur
    }
  } else{ # 2025 MIS
    dta_dir = paste0(user_dir, '/OneDrive - Bill & Melinda Gates Foundation/Data/Nigeria/snt_2026/MIS_2025')  # moved DHS data to OneDrive
    dhs_years = 2025
    
    # functions for reading data files
    parse_sps_layout <- function(sps_path) {
      lines <- readLines(sps_path, warn = FALSE)
      
      # Isolate the variable definition block between "/" and "."
      slash_line <- grep("^\\s*/$", lines)[1]
      dot_line   <- grep("^\\s*\\.$", lines)[1]
      
      if (is.na(slash_line) || is.na(dot_line))
        stop("Could not find DATA LIST block (/ ... .) in: ", basename(sps_path))
      
      block <- lines[(slash_line + 1):(dot_line - 1)]
      
      rows <- lapply(block, function(ln) {
        m <- regexec("(\\w+)\\s+(\\d+)-(\\d+)", ln)
        parts <- regmatches(ln, m)[[1]]
        if (length(parts) == 0) return(NULL)
        data.frame(
          Name   = parts[2],
          Start  = as.integer(parts[3]),
          End    = as.integer(parts[4]),
          IsChar = grepl("\\(A\\)", ln),
          stringsAsFactors = FALSE
        )
      })
      
      layout <- do.call(rbind, Filter(Negate(is.null), rows))
      layout[order(layout$Start), ]
    }
    
    read_dat_sps <- function(dat_path, sps_path = NULL) {
      if (is.null(sps_path))
        sps_path <- sub("(?i)\\.[^.]+$", ".SPS", dat_path, perl = TRUE)
      if (!file.exists(sps_path))
        stop("SPS file not found: ", sps_path)
      
      layout <- parse_sps_layout(sps_path)
      
      df <- readr::read_fwf(
        dat_path,
        col_positions  = readr::fwf_positions(
          start     = layout$Start,
          end       = layout$End,
          col_names = layout$Name
        ),
        col_types      = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress       = FALSE
      )
      
      as.data.frame(dplyr::mutate(df, dplyr::across(where(is.character), trimws)))
    }
    
    
    # get list of relevant dat files
    dta_filepaths = list.files(path=paste0(dta_dir, '/ExportedData'), pattern="\\.(DAT|dat)$", full.names=TRUE, recursive=TRUE)
    dta_filepaths = dta_filepaths[substr(dta_filepaths, nchar(dta_filepaths)-3+1, nchar(dta_filepaths)) %in% c('DAT', 'dat')]
    dta_list = list()
    for(dd in 1:length(dta_filepaths)){
      dta_cur = read_dat_sps(dta_filepaths[dd])
      dta_list[[dd]] = dta_cur
    }
    # dta_list_exported = dta_list
    # dta_filepaths_exported = dta_filepaths
    
    # # check with RawFile
    # dd = 1
    # dta_cur = read_cspr(paste0(dta_dir, '/RawFile/NGIQ89A.dat'), dcf_path = paste0(dta_dir,"/RawFile/NGIQ82.dcf"))
    # dta_list[[dd]] = dta_cur
  }
}


# function to check whether a particular code is included in each of the datasets and print out
find_code_locations = function(dta_list, code_str){
  
  # determine which columns with the code_str have non-NA values in each dta
  store_candidate_cols = list()
  for(ii in 1:length(dta_list)){
    dta_cur = dta_list[[ii]]
    cur_cols = grep(colnames(dta_cur), pattern = code_str)
    if(length(cur_cols)>0){
      # find columns that contain something other than NA
      can_cols = c()
      for(cc in 1:length(cur_cols)){
        if(any(!is.na(unique(dta_cur[,cur_cols[cc]])))) can_cols = c(can_cols, cur_cols[cc])
      }
      store_candidate_cols[[ii]] = can_cols
    }else{
      store_candidate_cols[[ii]] = NA
    }
  }
  if(length(store_candidate_cols)<length(dta_list)) store_candidate_cols[[length(dta_list)]] = NA
  # # print out head of relevant parts of data frame
  # for(ii in 1:length(dta_list)){
  #   if(length(store_candidate_cols[[ii]])>0){
  #     if(!is.na(store_candidate_cols[[ii]][1])){
  #       dta_cur = dta_list[[ii]]
  #       print(head(dta_cur[,store_candidate_cols[[ii]]]))
  #     }
  #   }
  # }
  # print out summary of relevant parts of data frame
  for(ii in 1:length(dta_list)){
    if(length(store_candidate_cols[[ii]])>0){
      if(!is.na(store_candidate_cols[[ii]][1])){      
        dta_cur = dta_list[[ii]]
        print(paste('ii=',ii,' (num rows=',nrow(dta_cur),'); column names:', colnames(dta_cur)[store_candidate_cols[[ii]]]))
        print(summary(dta_cur[,store_candidate_cols[[ii]]]))
      }
    }
  }
}


#####################
# once we know which files correspond to the household and individual outputs, can look at dataframe with just the relevant info and look at which values indicate a positive and negative result
#     note 1: if filenum not known yet, start with the next section, which allows an exploration of files
#     note 2: these are the codes that I used for BDI and NGA, but different codes may be more appropriate in other contexts
#####################
# view all the artesunate versus non-ACT antimalarials
art_codes = c("ml13e", "ml13aa", "ml13ab")
non_art_codes = c("ml13a", "ml13b", "ml13c", "ml13d", "ml13da", "ml13h") # country-specific antimalarial: "ml13g", "ml13f", 

# household codes 
house_codes = c('hhid','hvidx','hv001','hv006','hv007', 'hv024','hv105','hml20','hml32','hml32a')
house_filenum = 7 # 2010:5, 2013:7, 2015:5, 2018:8, 2021: 5, 2024: 7
dta_filepaths[house_filenum]
View(dta_list[[house_filenum]][1:50,house_codes[house_codes %in% colnames(dta_list[[house_filenum]])]])
for(ii in 6:length(house_codes)){
  if(house_codes[ii] %in% colnames(dta_list[[house_filenum]])){
    print(paste0(house_codes[ii],  '  unique values:' ))
    print(paste0('     ', unique(dta_list[[house_filenum]][,house_codes[ii]])))
  }
}

# individual codes
ind_codes = c('caseid','bidx','v001','v012','v014','v006','v007','v105', 'hw16', 'h32z','h47','m49a', 'ml1', 'h3', 'h5', 'h7', 'h9', art_codes, non_art_codes)
ind_filenum = 5 # 2010:1, 2013:1, 2015:4, 2018:1, 2021: 5
dta_filepaths[ind_filenum]
View(dta_list[[ind_filenum]][1:50,ind_codes[ind_codes %in% colnames(dta_list[[ind_filenum]])]])
for(ii in 6:length(ind_codes)){
  if(ind_codes[ii] %in% colnames(dta_list[[ind_filenum]])){
    print(paste0(ind_codes[ii],  '  unique values:' ))
    print(paste0('     ', unique(dta_list[[ind_filenum]][,ind_codes[ii]])))
  }
}
# View(dta_list[[5]][1:50,c('hhid','hvidx','hv001','hv006','hv007','hv105','hml20','hml32','hml32a')])



# vaccine with birth dates
ind_codes = c('caseid','bidx','v001','v012','v014','v006','v007','v105', 'hw16', 'h32z','h47','m49a', 'ml1', 'b1','b2','b3', 'h3', 'h5', 'h7', 'h9')
ind_filenum = 1 # 2010:1, 2013:1, 2015:4, 2018:1, 2021: 5
dta_filepaths[ind_filenum]
View(dta_list[[ind_filenum]][1:50,ind_codes[ind_codes %in% colnames(dta_list[[ind_filenum]])]])


# # quick check of 2025 MIS codes
# ind_codes = c('QHCLUST','QHNUMBER','QHLGOVAREA','QHLOCALITY','QHREGION','QH01','QH07','ML07I','ML07F','GHINTRO', 'GHLATPOLE', 'GHLATITUDE', 'GHLONGITUDE','QB114', 'QB112A', 'QB113','QB113B', 'QHLGOVAN','QH121','QH129','QH129A','QH129B','QH129C','QH129D','QH129E','Q404','Q405','Q406','Q407','Q408','Q412','Q413')
# View(distinct(dta_list[[5]][,which(colnames(dta_list[[5]]) %in% ind_codes)]))
# View(distinct(dta_list[[4]][,which(colnames(dta_list[[4]]) %in% ind_codes)]))
# View(distinct(dta_list[[1]][,which(colnames(dta_list[[1]]) %in% ind_codes)]))


dta_filepaths[7]
########################
# PfPR (microscopy)
########################
find_code_locations(dta_list=dta_list, code_str='hml32')
find_code_locations(dta_list=dta_list, code_str='hml33')
find_code_locations(dta_list=dta_list, code_str='hml35')
find_code_locations(dta_list=dta_list, code_str='HML35')
find_code_locations(dta_list=dta_list, code_str='QB114')
find_code_locations(dta_list=dta_list, code_str='QB112')


########################
# month of survey in cluster
########################
find_code_locations(dta_list=dta_list, code_str='v006')
find_code_locations(dta_list=dta_list, code_str='hv006')


########################
# ITNs
########################
# possible_code_strs = c('hml20','hml19','s508', '460','461')
find_code_locations(dta_list=dta_list, code_str='hml20')
find_code_locations(dta_list=dta_list, code_str='hml19')
find_code_locations(dta_list=dta_list, code_str='QH06')


########################
# source of ITNs
########################
# possible_code_strs = c('hml20','hml19','s508', '460','461')
find_code_locations(dta_list=dta_list, code_str='hml22')
find_code_locations(dta_list=dta_list, code_str='hml23')



######################
# treatment-seeking
#####################
find_code_locations(dta_list=dta_list, code_str='h32z')  # CM
find_code_locations(dta_list=dta_list, code_str='h32y')  # received_treatment

######################
# where received treatment
#####################
find_code_locations(dta_list=dta_list, code_str='h32')  # CM

############################################
# receive heel prick given seek treatment
############################################
find_code_locations(dta_list=dta_list, code_str='h47')

############################################
# take ACT for fever or difficulty breathing (among children with fever or difficulty breathing in last two weeks)
############################################
find_code_locations(dta_list=dta_list, code_str='h37e')  # combination with artemisinin
find_code_locations(dta_list=dta_list, code_str='ml13e')  # combination with artemisinin
find_code_locations(dta_list=dta_list, code_str='ml20a')  # how long after fever started did first take ACT
# # view all the artesunate versus non-ACT antimalarials
# art_codes = c("ml13e", "ml13aa", "ml13ab")  # 'ml13aa' is rectal artesunate, moving to non-ACT
# non_art_codes = c("ml13a", "ml13b", "ml13c", "ml13d", "ml13da", "ml13h")  #, "ml13g", "ml13f") # country-specific antimalarial: "ml13g", "ml13f", 
find_code_locations(dta_list=dta_list, code_str='ml13aa')  # rectal artesunate



###################################
# IPTp >=1 dose
###################################
find_code_locations(dta_list=dta_list, code_str='m49a')


###################################
# IPTp number of doses
###################################
find_code_locations(dta_list=dta_list, code_str='ml1')


###################################
# vaccinations (proxy for IPTi coverage)
###################################
find_code_locations(dta_list=dta_list, code_str='h3')  # DPT1 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h5')  # DPT2 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h7')  # DPT3 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h4')  # Polio1 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h6')  # Polio2 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h8')  # Polio3 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h9')  # Measles vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h51')  # Penta1 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h52')  # Penta2 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)
find_code_locations(dta_list=dta_list, code_str='h53')  # Penta3 vaccination (1=date shown on health card, 2=no health card, but respondent reported the child received the vaccine, 3=health card has vaccine but not date of vaccine)

find_code_locations(dta_list=dta_list, code_str='q514')  # polio vaccination
find_code_locations(dta_list=dta_list, code_str='q516')  # number of polio vaccinations
find_code_locations(dta_list=dta_list, code_str='q517')  # pentavalent vaccination
find_code_locations(dta_list=dta_list, code_str='q518')  # number of pentavalent vaccinations



###################################
# time between fever onset and ACT
###################################
find_code_locations(dta_list=dta_list, code_str='ml20a')


###################################
# child took combination with artemisinin
###################################
#BASE for ML13A-Z: Children suffering from fever or short rapid breaths or difficulty
# breathing in the last two weeks (H22 = 1 or H31B = 1). In previous recodes the base was
# restricted to children with fever in the last 2 weeks (H22 = 1).
# Questions pertaining to ML14A to ML14Z are no longer part of the DHS VII core
# questionnaire, but the variables are kept in the DHS VII recode
find_code_locations(dta_list=dta_list, code_str='ml13e')



###################################
# months ago household obtained net
###################################
find_code_locations(dta_list=dta_list, code_str='hml4')

###################################
# time since individual received net
###################################
hist(dta_list[[8]]$hml4[dta_list[[8]]$hml4<80], xlab='months ago received net', main='net age at time of 2016 DHS')

length(which(dta_list[[8]]$hml4[dta_list[[8]]$hml4<80] < 12)) / length(which(dta_list[[8]]$hml4[dta_list[[8]]$hml4<80]<80))

length(which(dta_list[[8]]$hml4[dta_list[[8]]$hml4<80] > 14)) / length(which(dta_list[[8]]$hml4[dta_list[[8]]$hml4<80]<80))


###################################
# number of ANC visits
###################################
find_code_locations(dta_list=dta_list, code_str='m14')















# =============================================================================== #
#     old
# =============================================================================== #



########################
# PfPR (microscopy)
########################

possible_code_strs = c('hml32', 'hml32a', 'hml33', 'hml35')
code_str = possible_code_strs[4]

find_code_locations(dta_list=dta_list, code_str='hml32')
find_code_locations(dta_list=dta_list, code_str='hml33')
find_code_locations(dta_list=dta_list, code_str='hml35')


# it appears that dta_1's 'hml32' entry is the one we want to use for PfPR
# one of the columns of dta_1 is hv020, of which all entries are 'all woman sample,' which I think is okay looking at the recode manual: I think it just indicates how information was recorded for the household surveys
#      unique(dta_1[,grep(colnames(dta_1), pattern = 'hv020')])
# also for dta_1, there is a question on the number of children aged 5 and under, and this takes values from 0-8.  since the microscopy entry is 
#     only positive or negative in each line, does it correspond to just one child given the test per household? or is it asking if any of the children 
#     are positive? ah, okay, hhid gives the household id, then within a household, I believe hvidx gives the index of the individual. 
#     there appears to be one entry per individual in the household, so the micropscopy results are per individual. this is supported by the 'hv105' variable, which is the age of the individual
#        table((dta_1[,grep(colnames(dta_1), pattern = 'hv014')]))
#        dta_1[which(dta_1[,grep(colnames(dta_1), pattern = 'hv014')]==3)[1:10],c(1:19, 160, 262:269)]
#        # subset to children under 5
#        dta_1[intersect(which(dta_1[,grep(colnames(dta_1), pattern = 'hv105')]<=5), which(dta_1[,grep(colnames(dta_1), pattern = 'hv014')]==3))[1:10],c(1:19, 160, 262:269)]
#        table(dta_1[which(dta_1[,grep(colnames(dta_1), pattern = 'hv105')]<=5),262])

# I think 'hv001' gives the cluster id - there appear to be 336 clusters in 2016
#        table((dta_1[,grep(colnames(dta_1), pattern = 'hv001')]))

# comparing dta_4 with dta_1, I think there's a lot of duplicate/overlap in that dta_4 has one row per household, but then breaks out certain entries by individual, so there is hml32_01, hml32_02, etc. for each person in a house. 


# get mean PfPR as well as number of children with reported values for each cluster
# subset to non-na microscopy values
# dta_mic = dta_1[!is.na(dta_1[,grep(colnames(dta_1), pattern = 'hml32')]),]
dta_1$mic_pos = NA
dta_1$mic_pos[which(dta_1[,grep(colnames(dta_1), pattern = 'hml32')] == 'positive')] = 1
dta_1$mic_pos[which(dta_1[,grep(colnames(dta_1), pattern = 'hml32')] == 'negative')] = 0
dta_pfpr_cluster = dta_1  %>%
  filter(!is.na(dta_1$mic_pos)) %>%
  group_by(hv001) %>%
  summarize(pfpr = mean(mic_pos, na.rm = TRUE),
            num_pos = sum(mic_pos),
            mean_age = mean(hv105, na.rm=TRUE),
            num_tested = n()) 

# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
pfpr_cluster = merge(dta_pfpr_cluster, map, by.x='hv001', by.y='clusterid', all=TRUE)
pfpr_cluster = pfpr_cluster[,c('hv001','pfpr','num_pos','mean_age','num_tested', 'SurveyName','latitude','longitude','urban_rural')]

layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
plot(pfpr_cluster$longitude, pfpr_cluster$latitude, col=rev(rainbow(150, alpha=0.5))[50+round(pfpr_cluster$pfpr*100)], pch=20, cex=pfpr_cluster$num_tested/8, main='PfPR (microscopy, U5)', xlab='longitude', ylab='latitude')
plot(rep(0,100), seq(0,1,length.out=100), col=rev(rainbow(150, alpha=0.5))[50+round(seq(0,1,length.out=100)*100)], pch=15, cex=1.2, ylab='PfPR', axes=FALSE, xlab=''); axis(2)
plot(rep(0,5), seq(round(min(pfpr_cluster$num_tested)), round(max(pfpr_cluster$num_tested)), length.out=5), cex=seq(round(min(pfpr_cluster$num_tested)), round(max(pfpr_cluster$num_tested)), length.out=5)/8, pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)


MIS_2016 = pfpr_cluster
colnames(MIS_2016)[colnames(MIS_2016)=='pfpr'] = 'mic_rate'
colnames(MIS_2016)[colnames(MIS_2016)=='num_pos'] = 'mic_num_true'
colnames(MIS_2016)[colnames(MIS_2016)=='num_tested'] = 'mic_num_total'
colnames(MIS_2016)[colnames(MIS_2016)=='mean_age'] = 'mic_mean_age'






########################
# ITNs
########################


possible_code_strs = c('hml20','hml19','s508', '460','461')
code_str = possible_code_strs[1]
code_str='hml20'

find_code_locations(dta_list=dta_list, code_str='hml20')
find_code_locations(dta_list=dta_list, code_str='hml19')


# there is a 's508' in both dta_2 and dta_3. there is also hml19 in dta_1 and dta_4 ('person slept under an ever treated net') <-- this one is good in that it is for household members. there's also hml20, which is "person slept under LLIN net." the hml19 and hml20 numbers appear fairly similar. I'll plan to use hml20
# age distribution of people reporting use of nets
# hist(dta_1$hv105[which(dta_1[,grep(colnames(dta_1), pattern = code_str)] == 'yes')])


# get mean PfPR as well as number of children with reported values for each cluster
# subset to non-na microscopy values
# dta_mic = dta_1[!is.na(dta_1[,grep(colnames(dta_1), pattern = 'hml32')]),]
dta_1$used_itn = NA
dta_1$used_itn[which(dta_1[,grep(colnames(dta_1), pattern = code_str)] == 'yes')] = 1
dta_1$used_itn[which(dta_1[,grep(colnames(dta_1), pattern = code_str)] == 'no')] = 0
dta_itn_cluster = dta_1  %>%
  filter(!is.na(dta_1$used_itn)) %>%
  group_by(hv001) %>%
  summarize(itn_use_rate = mean(used_itn, na.rm = TRUE),
            num_pos = sum(used_itn),
            mean_age = mean(hv105, na.rm=TRUE),
            num_included = n()) 



# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
itn_cluster = merge(dta_itn_cluster, locations, by.x='hv001', by.y='clusterid', all=TRUE)
itn_cluster = itn_cluster[,c('hv001','itn_use_rate','num_pos','mean_age','num_included', 'SurveyName','latitude','longitude','urban_rural')]
layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
plot(itn_cluster$longitude, itn_cluster$latitude, col=(rainbow(150, alpha=0.5))[1+round(itn_cluster$itn_use_rate*100)], pch=20, cex=itn_cluster$num_included/40, main='ITN use (all age)', xlab='longitude', ylab='latitude')
plot(rep(0,100), seq(0,1,length.out=100), col=(rainbow(150, alpha=0.5))[1+round(seq(0,1,length.out=100)*100)], pch=15, cex=1.2, ylab='rate', axes=FALSE, xlab=''); axis(2)
plot(rep(0,5), seq(round(min(itn_cluster$num_included)), round(max(itn_cluster$num_included)), length.out=5), cex=seq(round(min(itn_cluster$num_included)), round(max(itn_cluster$num_included)), length.out=5)/40, pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)


# combine with other MIS data points
MIS_2016 = merge(MIS_2016, itn_cluster, by=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'))
colnames(MIS_2016)[colnames(MIS_2016)=='itn_use_rate'] = 'itn_all_rate'
colnames(MIS_2016)[colnames(MIS_2016)=='num_pos'] = 'itn_all_num_true'
colnames(MIS_2016)[colnames(MIS_2016)=='num_included'] = 'itn_all_num_total'
colnames(MIS_2016)[colnames(MIS_2016)=='mean_age'] = 'itn_all_mean_age'




# break out U5 and O5 ITN use rates
# U5
dta_1_u5 = dta_1[dta_1$hv105<=5,]
dta_1_u5$used_itn = NA
dta_1_u5$used_itn[which(dta_1_u5[,grep(colnames(dta_1_u5), pattern = code_str)] == 'yes')] = 1
dta_1_u5$used_itn[which(dta_1_u5[,grep(colnames(dta_1_u5), pattern = code_str)] == 'no')] = 0
dta_itn_cluster_u5 = dta_1_u5  %>%
  filter(!is.na(dta_1_u5$used_itn)) %>%
  group_by(hv001) %>%
  summarize(itn_use_rate = mean(used_itn, na.rm = TRUE),
            num_pos = sum(used_itn),
            mean_age = mean(hv105, na.rm=TRUE),
            num_included = n()) 
# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
# itn_cluster_u5 = merge(dta_itn_cluster_u5, map, by.x='hv001', by.y='clusterid', all=TRUE)
# itn_cluster_u5 = itn_cluster_u5[,c('hv001','itn_use_rate','num_pos','mean_age','num_included', 'SurveyName','latitude','longitude','urban_rural')]
# plot(itn_cluster_u5$longitude, itn_cluster_u5$latitude, col=(rainbow(150, alpha=0.5))[1+round(itn_cluster_u5$itn_use_rate*100)], pch=20, cex=itn_cluster_u5$num_included/20, main='ITN use (U5)', xlab='longitude', ylab='latitude')

# combine with other MIS data points
# MIS_2016 = merge(MIS_2016, itn_cluster_u5, by=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'))
# colnames(MIS_2016)[colnames(MIS_2016)=='itn_use_rate'] = 'itn_u5_rate'
# colnames(MIS_2016)[colnames(MIS_2016)=='num_pos'] = 'itn_u5_num_true'
# colnames(MIS_2016)[colnames(MIS_2016)=='num_included'] = 'itn_u5_num_total'
# colnames(MIS_2016)[colnames(MIS_2016)=='mean_age'] = 'itn_u5_mean_age'


# O5
dta_1_o5 = dta_1[dta_1$hv105>5,]
dta_1_o5$used_itn = NA
dta_1_o5$used_itn[which(dta_1_o5[,grep(colnames(dta_1_o5), pattern = code_str)] == 'yes')] = 1
dta_1_o5$used_itn[which(dta_1_o5[,grep(colnames(dta_1_o5), pattern = code_str)] == 'no')] = 0
dta_itn_cluster_o5 = dta_1_o5  %>%
  filter(!is.na(dta_1_o5$used_itn)) %>%
  group_by(hv001) %>%
  summarize(itn_use_rate = mean(used_itn, na.rm = TRUE),
            num_pos = sum(used_itn),
            mean_age = mean(hv105, na.rm=TRUE),
            num_included = n()) 
# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
# itn_cluster_o5 = merge(dta_itn_cluster_o5, locations, by.x='hv001', by.y='clusterid', all=TRUE)
# itn_cluster_o5 = itn_cluster[,c('hv001','itn_use_rate','num_pos','mean_age','num_included', 'SurveyName','latitude','longitude','urban_rural')]
# plot(itn_cluster_o5$longitude, itn_cluster_o5$latitude, col=(rainbow(150, alpha=0.5))[1+round(itn_cluster_o5$itn_use_rate*100)], pch=20, cex=itn_cluster_o5$num_included/80, main='ITN use (O5)', xlab='longitude', ylab='latitude')

# # fairly good linear relationship between I5 and O5 ITN use rates. not quite 1:1, but not that far from it.
# itn_cluster_age = merge(itn_cluster_u5, itn_cluster_o5[,c('hv001','itn_use_rate')], by=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'))
# plot(itn_cluster_age$itn_use_rate.x, itn_cluster_age$itn_use_rate.y, xlim=c(0,1), ylim=c(0,1), xlab='u5',ylab='o5')
# abline(a=0,b=1)

# combine with other MIS data points
# MIS_2016 = merge(MIS_2016, itn_cluster_o5, by=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'))
# colnames(MIS_2016)[colnames(MIS_2016)=='itn_use_rate'] = 'itn_o5_rate'
# colnames(MIS_2016)[colnames(MIS_2016)=='num_pos'] = 'itn_o5_num_true'
# colnames(MIS_2016)[colnames(MIS_2016)=='num_included'] = 'itn_o5_num_total'
# colnames(MIS_2016)[colnames(MIS_2016)=='mean_age'] = 'itn_o5_mean_age'



######################
# treatment-seeking
#####################

possible_code_strs = c('h32z')
code_str = possible_code_strs[1]

find_code_locations(dta_list=dta_list, code_str='h32z')


# I'll use the value from dta_2, which seems to have each row correspond to an individual, while dta_3 has separate columns for multiple individuals with fever.
# this variable is for children born in the past five years. asks whether child was taken to a medical facility for treatment of the fever and/or cough, among children who had fever/cough in last two weeks

# subset to non-na
dta_2$treat = NA
dta_2$treat[which(dta_2[,grep(colnames(dta_2), pattern = 'h32z')] == 'yes')] = 1
dta_2$treat[which(dta_2[,grep(colnames(dta_2), pattern = 'h32z')] == 'no')] = 0
dta_treat_cluster = dta_2  %>%
  filter(!is.na(dta_2$treat)) %>%
  group_by(v001) %>%
  summarize(treat_rate = mean(treat, na.rm = TRUE),
            num_treat = sum(treat),
            num_w_fever = n()) 

# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
treat_cluster = merge(dta_treat_cluster, locations, by.x='v001', by.y='clusterid', all=TRUE)
treat_cluster = treat_cluster[,c('v001','treat_rate','num_treat','num_w_fever', 'SurveyName','latitude','longitude','urban_rural')]
layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
plot(treat_cluster$longitude, treat_cluster$latitude, col=(rainbow(150, alpha=0.5))[1+round(treat_cluster$treat_rate*100)], pch=20, cex=treat_cluster$num_w_fever/2, main='seek treatment for fever (U5)', xlab='longitude', ylab='latitude')
plot(rep(0,100), seq(0,1,length.out=100), col=(rainbow(150, alpha=0.5))[1+round(seq(0,1,length.out=100)*100)], pch=15, cex=1.2, ylab='rate', axes=FALSE, xlab=''); axis(2)
plot(rep(0,5), seq(round(min(treat_cluster$num_w_fever, na.rm=TRUE)), round(max(treat_cluster$num_w_fever, na.rm=TRUE)), length.out=5), cex=seq(round(min(treat_cluster$num_w_fever, na.rm=TRUE)), round(max(treat_cluster$num_w_fever, na.rm=TRUE)), length.out=5)/2, pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)

# hist(treat_cluster$treat_rate)
# sum(treat_cluster$treat_rate * treat_cluster$num_w_fever, na.rm=TRUE)/sum(treat_cluster$num_w_fever, na.rm=TRUE)



# combine with other MIS data points
MIS_2016 = merge(MIS_2016, treat_cluster, by.x=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'), by.y=c('v001', 'SurveyName', 'latitude','longitude','urban_rural'))
colnames(MIS_2016)[colnames(MIS_2016)=='treat_rate'] = 'cm_rate'
colnames(MIS_2016)[colnames(MIS_2016)=='num_treat'] = 'cm_num_true'
colnames(MIS_2016)[colnames(MIS_2016)=='num_w_fever'] = 'cm_num_total'





############################################
# receive heel prick given seek treatment
############################################


possible_code_strs = c('h47')
code_str = possible_code_strs[1]

# grep(colnames(dta_1), pattern = code_str)  # several columns contain hml32
# grep(colnames(dta_2), pattern = code_str)
# grep(colnames(dta_3), pattern = code_str)
# grep(colnames(dta_4), pattern = code_str) # numerous columns contain hml32
# grep(colnames(dta_5), pattern = code_str)
# 
# dta_2[1:5,grep(colnames(dta_2), pattern = code_str)]
# dta_3[1:5,grep(colnames(dta_3), pattern = code_str)]
# 
# # determine which columns with the code_str have non-NA values in each dta
# store_candidate_cols = list()
# for(ii in 1:length(dta_list)){
#   dta_cur = dta_list[[ii]]
#   cur_cols = grep(colnames(dta_cur), pattern = code_str)
#   if(length(cur_cols)>0){
#     # find columns that contain something other than NA
#     can_cols = c()
#     for(cc in 1:length(cur_cols)){
#       if(any(!is.na(unique(dta_cur[,cur_cols[cc]])))) can_cols = c(can_cols, cur_cols[cc])
#     }
#     store_candidate_cols[[ii]] = can_cols
#   }else{
#     store_candidate_cols[[ii]] = NA
#   }
# }
# 
# # # print out head of relevant parts of data frame
# # for(ii in 1:length(dta_list)){
# #   if(!is.na(store_candidate_cols[[ii]][1]) & length(store_candidate_cols[[ii]])>0){
# #     dta_cur = dta_list[[ii]]
# #     print(head(dta_cur[,store_candidate_cols[[ii]]]))
# #   }
# # }
# 
# # print out summary of relevant parts of data frame
# for(ii in 1:length(dta_list)){
#   if(!is.na(store_candidate_cols[[ii]][1]) & length(store_candidate_cols[[ii]])>0){
#     dta_cur = dta_list[[ii]]
#     print(paste('ii=',ii))
#     print(summary(dta_cur[,store_candidate_cols[[ii]]]))
#   }
# }

# I'll use the value from dta_2, which seems to have each row correspond to an individual, while dta_3 has separate columns for multiple individuals with fever.
# this variable is for children born in the past five years. asks whether blood was taken for testing, among children who had fever/cough in last two weeks

# subset to non-na
dta_2$treat = NA
dta_2$treat[which(dta_2[,grep(colnames(dta_2), pattern = 'h47')] == 'yes')] = 1
dta_2$treat[which(dta_2[,grep(colnames(dta_2), pattern = 'h47')] == 'no')] = 0
dta_treat_cluster = dta_2  %>%
  filter(!is.na(dta_2$treat)) %>%
  group_by(v001) %>%
  summarize(blood_test_rate = mean(treat, na.rm = TRUE),
            num_blood_test = sum(treat),
            num_w_fever_blood = n()) 

# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
treat_cluster = merge(dta_treat_cluster, locations, by.x='v001', by.y='clusterid', all=TRUE)
treat_cluster = treat_cluster[,c('v001','blood_test_rate','num_blood_test','num_w_fever_blood', 'SurveyName','latitude','longitude','urban_rural')]
layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
plot(treat_cluster$longitude, treat_cluster$latitude, col=(rainbow(150, alpha=0.5))[1+round(treat_cluster$blood_test_rate*100)], pch=20, cex=treat_cluster$num_w_fever_blood/2, main='blood sample with fever (U5)', xlab='longitude', ylab='latitude')
plot(rep(0,100), seq(0,1,length.out=100), col=(rainbow(150, alpha=0.5))[1+round(seq(0,1,length.out=100)*100)], pch=15, cex=1.2, ylab='rate', axes=FALSE, xlab=''); axis(2)
plot(rep(0,5), seq(round(min(treat_cluster$num_w_fever_blood, na.rm=TRUE)), round(max(treat_cluster$num_w_fever_blood, na.rm=TRUE)), length.out=5), cex=seq(round(min(treat_cluster$num_w_fever_blood, na.rm=TRUE)), round(max(treat_cluster$num_w_fever_blood, na.rm=TRUE)), length.out=5)/2, pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)


# hist(treat_cluster$treat_rate)
# sum(treat_cluster$treat_rate * treat_cluster$num_w_fever, na.rm=TRUE)/sum(treat_cluster$num_w_fever, na.rm=TRUE)


# combine with other MIS data points
MIS_2016 = merge(MIS_2016, treat_cluster, by.x=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'), by.y=c('v001', 'SurveyName', 'latitude','longitude','urban_rural'))
colnames(MIS_2016)[colnames(MIS_2016)=='num_blood_test'] = 'blood_test_num_true'
colnames(MIS_2016)[colnames(MIS_2016)=='num_w_fever_blood'] = 'blood_test_num_total'





###################################
# IPTp >=1 dose
###################################



possible_code_strs = c('m49a')
code_str = possible_code_strs[1]

# grep(colnames(dta_1), pattern = code_str)  # several columns contain hml32
# grep(colnames(dta_2), pattern = code_str)
# grep(colnames(dta_3), pattern = code_str)
# grep(colnames(dta_4), pattern = code_str) # numerous columns contain hml32
# grep(colnames(dta_5), pattern = code_str)
# 
# dta_2[1:5,grep(colnames(dta_2), pattern = code_str)]
# dta_3[1:5,grep(colnames(dta_3), pattern = code_str)]
# 
# # determine which columns with the code_str have non-NA values in each dta
# store_candidate_cols = list()
# for(ii in 1:length(dta_list)){
#   dta_cur = dta_list[[ii]]
#   cur_cols = grep(colnames(dta_cur), pattern = code_str)
#   if(length(cur_cols)>0){
#     # find columns that contain something other than NA
#     can_cols = c()
#     for(cc in 1:length(cur_cols)){
#       if(any(!is.na(unique(dta_cur[,cur_cols[cc]])))) can_cols = c(can_cols, cur_cols[cc])
#     }
#     store_candidate_cols[[ii]] = can_cols
#   }else{
#     store_candidate_cols[[ii]] = NA
#   }
# }
# 
# # # print out head of relevant parts of data frame
# # for(ii in 1:length(dta_list)){
# #   if(!is.na(store_candidate_cols[[ii]][1]) & length(store_candidate_cols[[ii]])>0){
# #     dta_cur = dta_list[[ii]]
# #     print(head(dta_cur[,store_candidate_cols[[ii]]]))
# #   }
# # }
# 
# # print out summary of relevant parts of data frame
# for(ii in 1:length(dta_list)){
#   if(!is.na(store_candidate_cols[[ii]][1]) & length(store_candidate_cols[[ii]])>0){
#     dta_cur = dta_list[[ii]]
#     print(paste('ii=',ii))
#     print(summary(dta_cur[,store_candidate_cols[[ii]]]))
#   }
# }

# I'll use the value from dta_2, which seems to have each row correspond to an individual, while dta_3 has separate columns for multiple individuals with fever.
# this variable is for children born in the past five years. asks whether blood was taken for testing, among children who had fever/cough in last two weeks

# subset to non-na
dta_2$treat = NA
dta_2$treat[which(dta_2[,grep(colnames(dta_2), pattern = 'm49a')] == 'yes')] = 1
dta_2$treat[which(dta_2[,grep(colnames(dta_2), pattern = 'm49a')] == 'no')] = 0
dta_treat_cluster = dta_2  %>%
  filter(!is.na(dta_2$treat)) %>%
  group_by(v001) %>%
  summarize(iptp_rate = mean(treat, na.rm = TRUE),
            num_iptp = sum(treat),
            num_preg = n()) 

# can we now match each entry with a cluster and a location?  - none of the dtas appear to have gps coordinates associated with entries. perhaps somewhere else in DHS/MIS records?
# well, in the MAP-extracted outputs from Kate, there are GPS coordinates associated with 336 clusters, which seems promising, assuming the cluster ids match up correctly.
# match 'hv001' with 'clusterid'
treat_cluster = merge(dta_treat_cluster, locations, by.x='v001', by.y='clusterid', all=TRUE)
treat_cluster = treat_cluster[,c('v001','iptp_rate','num_iptp','num_preg', 'SurveyName','latitude','longitude','urban_rural')]
plot(treat_cluster$longitude, treat_cluster$latitude, col=(rainbow(150, alpha=0.5))[1+round(treat_cluster$iptp_rate*100)], pch=20, cex=treat_cluster$num_preg/5, main='receive IPTp', xlab='longitude', ylab='latitude')
# plot(seq(0,1,length.out=100),col=(rainbow(150, alpha=0.5))[1+round(seq(0,1,length.out=100)*100)], pch=20)
 
# combine with other MIS data points
MIS_2016 = merge(MIS_2016, treat_cluster, by.x=c('hv001', 'SurveyName', 'latitude','longitude','urban_rural'), by.y=c('v001', 'SurveyName', 'latitude','longitude','urban_rural'))
colnames(MIS_2016)[colnames(MIS_2016)=='num_iptp'] = 'iptp_num_true'
colnames(MIS_2016)[colnames(MIS_2016)=='num_preg'] = 'iptp_num_total'







###################################
# IPTp distribution of number of doses
###################################


possible_code_strs = c('ml1')
code_str = possible_code_strs[1]

colnames(dta_1)[grep(colnames(dta_1), pattern = code_str)]  # several columns contain hml32
colnames(dta_2)[grep(colnames(dta_2), pattern = code_str)]
colnames(dta_3)[grep(colnames(dta_3), pattern = code_str)]
colnames(dta_4)[grep(colnames(dta_4), pattern = code_str)] # numerous columns contain hml32
colnames(dta_5)[grep(colnames(dta_5), pattern = code_str)]

# times took SP/Fansidar during pregnancy, given that it was taken at least once -> 3 now means >=3
table(dta_2$ml1)
dta_2$iptp_num_doses = dta_2$ml1
dta_2$iptp_num_doses[dta_2$iptp_num_doses == 0] = NA
dta_2$iptp_num_doses[dta_2$iptp_num_doses > 15] = NA
dta_2$iptp_num_doses[dta_2$iptp_num_doses > 3] = 3
table(dta_2$iptp_num_doses)
# use same ratio across country due to small numbers when break into clusters













# ##################################################################################
# # determine which clusters are in which chiefdoms
# ##################################################################################
# 
# 
# # turn MIS_2016 into spatial points data frame
# # not sure what the spatial projection is... will try with same spatial projection as shapefile
# points_crs = crs(admin_shape)
# MIS_2016_shape = SpatialPointsDataFrame(MIS_2016[,c('longitude', 'latitude')],
#                                         MIS_2016,
#                                         proj4string = points_crs)
# # find which chiefdom each cluster belongs to
# MIS_2016_locations = over(MIS_2016_shape, admin_shape)
# if(nrow(MIS_2016_locations) == nrow(MIS_2016_shape)){
#   MIS_2016_shape$NAME_3 = MIS_2016_locations$NAME_3
#   MIS_2016_shape$NAME_2 = MIS_2016_locations$NAME_2
#   MIS_2016_shape$NAME_1 = MIS_2016_locations$NAME_1
# }
# 
# write.csv(as.data.frame(MIS_2016_shape), 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/SierraLeone_hbhi/explore_DHS/MIS_2016.csv')
# 
# 
# par(mfrow=c(3,4))
# # for each of the interventions/measures, count number of individuals surveyed in each chiefdom
# num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_3','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
#   group_by(NAME_3) %>%
#   summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
#             num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
#             num_w_fever=sum(num_w_fever, na.rm = TRUE),
#             num_preg=sum(num_preg, na.rm = TRUE))
# num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
# hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,800, length.out=80))
# hist(num_surveyed$num_use_itn_all,  main='ITN', breaks=seq(0,1600, length.out=80))
# hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,200, length.out=80))
# hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,500, length.out=80))
# 
# 
# # look at survey numbers by district
# # for each of the interventions/measures, count number of individuals surveyed in each chiefdom
# num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_2','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
#   group_by(NAME_2) %>%
#   summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
#             num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
#             num_w_fever=sum(num_w_fever, na.rm = TRUE),
#             num_preg=sum(num_preg, na.rm = TRUE))
# num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
# hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,800, length.out=80))
# hist(num_surveyed$num_use_itn_all,  main='ITN', breaks=seq(0,1600, length.out=80))
# hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,200, length.out=80))
# hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,500, length.out=80))
# 
# # look at survey numbers by region
# # for each of the interventions/measures, count number of individuals surveyed in each chiefdom
# num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_1','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
#   group_by(NAME_1) %>%
#   summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
#             num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
#             num_w_fever=sum(num_w_fever, na.rm = TRUE),
#             num_preg=sum(num_preg, na.rm = TRUE))
# num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
# hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,3000, length.out=80))
# hist(num_surveyed$num_use_itn_all,  main='ITN',breaks=seq(0,6000, length.out=80))
# hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,700, length.out=80))
# hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,2000, length.out=80))
# par(mfrow=c(1,1))
# 
# 
# 
# 
# # get weighted means and number tested/positive within each chiefdom for all variables



