# standardize_admin_names.R
# November 2022
# contact: mambrose
#
# Goal: the names of LGAs are often different across different files and we want to make sure that they are used consistently
#    This function takes a 'target' naming system and a 'origin' set of names and matches the origin LGA names with the target names.
#    It then updates the origin names so that the same LGA names can be used consistently across all files.
# NOTE OF CAUTION: some NGA admin2 names appear in multiple admin1s, so they must be renamed to end in a 1 or 2, depending on the admin1. Currently done manually. Includes:
# 'BASSA'
# 'IFELODUN'
# 'IREPODUN'
# 'NASARAWA'
# 'OBI'
# 'SURULERE'
library(stringr)

# if there are admin names that shared across multiple admins (without their numeric distintuishers) in the inputs, will need to have user modify
duplicate_names = c('BASSA', 'IFELODUN', 'IREPODUN', 'NASARAWA', 'OBI', 'SURULERE')  # LGA
duplicate_names_state = c()

create_reference_name_match = function(lga_name){
  #' alter string lga_name to standardized admin2 (LGA) format, also check for and fix common misspellings
  #' @param lga_name the name of the admin2 that should be changed to standardized format
  
  lga_name = str_replace_all(lga_name, pattern=" State", replacement='')
  lga_name = str_replace_all(lga_name, pattern=" state", replacement='')
  lga_name = str_replace_all(lga_name, pattern=' ', replacement='-')
  lga_name = str_replace_all(lga_name, pattern='/', replacement='-')
  lga_name = str_replace_all(lga_name, pattern='_', replacement='-')
  lga_name = str_replace_all(lga_name, pattern="'", replacement='')
  lga_name = toupper(lga_name)
  
  # first value (the name) is the one that is replaced by the second value
  replace_list = list(
                      # NGA 
                      # 'IFELODUN' = 'IFELODUN1', 'IFELODUN2',
                      # 'IREPODUN' = 'IREPODUN1', 'IREPODUN2',
                      # 'BASSA' = 'BASSA1', 'BASSA2',
                      # 'NASARAWA' = 'NASARAWA1', 'NASARAWA2',
                      # 'SURULERE' = 'SURULERE1', 'SURULERE2',
                      # 'OBI' = 'OBI2',
                      'K-ALA' = 'KATSINA-ALA',
                      'OSHONGO' = 'USHONGO',
                      'OTUPKO' = 'OTURKPO',
                      'VANDEKYA' = 'VANDEIKYA',
                      'MUNGUNO' = 'MONGUNO',
                      'KANO-MUNICIPLAL' = 'KANO-MUNICIPAL',    
                      'QUAANPAN' = 'QUAAN-PAN',
                      'NASSARAWA' = 'NASARAWA',
                      'GANYE' = 'GANAYE',
                      'GIREI' = 'GIRERI',
                      'DAMBAM' = 'DAMBAN',
                      'MARKURDI' = 'MAKURDI',
                      'MUNICIPAL-AREA-COUNCIL' = 'ABUJA-MUNICIPAL',
                      'AMAC' = 'ABUJA-MUNICIPAL',
                      'SHONGOM' = 'SHOMGOM',
                      # 'BIRNIN' = 'BIRINIWA', 'BIRNI-KUDU'?
                      'BIRNIN-KUDU' = 'BIRNI-KUDU',
                      'BIRNIN' = 'BIRNI-KUDU',
                      'BIRNIWA' = 'BIRINIWA',
                      'KIRI-KASAMA' = 'KIRI-KASAMMA',
                      'MALAM-MADURI' = 'MALAM-MADORI',
                      'SULE-TANKAKAR' = 'SULE-TANKARKAR',
                      'MAKARFI' = 'MARKAFI',
                      'ZANGON-KATAF' = 'ZANGO-KATAF',
                      'DANBATTA' = 'DAMBATTA',
                      'TAKAI' = 'TAKALI',
                      'UNGONGO' = 'UNGOGO',
                      'DAN' = 'DAN-MUSA',
                      'DUTSIN' = 'DUTSIN-MA',
                      'MATAZU' = 'MATAZUU',
                      'PATIGI' = 'PATEGI',
                      'MUNYA' = 'MUYA',
                      'PAIKORO' = 'PAILORO',
                      'BARKIN-LADI' = 'BARIKIN-LADI',
                      'DANGE-SHUNI' = 'DANGE-SHNSI',
                      'GWADABAWA' = 'GAWABAWA',
                      'WAMAKKO' = 'WAMAKO',
                      'ARDO' = 'ARDO-KOLA',
                      'KARIM' = 'KURMI',
                      'BARDE' = 'BADE',
                      'BORSARI' = 'BURSARI',
                      'TARMUWA' = 'TARMUA',
                      'EGBADO-NORTH' = 'YEWA-NORTH',
                      'EGBADO-SOUTH' = 'YEWA-SOUTH',
                      'ILEMEJI' = 'ILEJEMEJE',
                      'GBOYIN' = 'AIYEKIRE-(GBOYIN)',
                      'UMU-NEOCHI' = 'UMU-NNEOCHI',
                      'IBADAN-NORTH-EAST' = 'IBADAN-CENTRAL-(IBADAN-NORTH-EAST)',
                      'NAFADA' = 'NAFADA-(BAJOGA)',
                      'ISUIKWUATO' = 'ISUIKWATO',
                      'OBI-NWGA' = 'OBI-NWA',
                      'UMUNNEOCHI' = 'UMU-NNEOCHI',
                      'YENAGOA' = 'YENEGOA',
                      'OTUKPO' = 'OTURKPO',
                      'BEKWARRA' = 'BEKWARA',
                      'EFON' = 'EFON-ALAYEE',
                      'AIYEKIRE-(GBONYIN)' = 'AIYEKIRE-(GBOYIN)',
                      'IDO-OSI' = 'IDOSI-OSI',
                      'ABUJA-MUNICIPAL-AREA-COUNCIL' = 'ABUJA-MUNICIPAL',
                      'EZINIHITTE-MBAISE' = 'EZINIHITTE',
                      'ONUIMO' = 'UNUIMO',
                      'AREWA' = 'AREWA-DANDI',
                      'DANKO-WASAGU' = 'WASAGU-DANKO',
                      'OGORI-MAGONGO' = 'OGORI-MANGONGO',
                      'OLAMABORO' = 'OLAMABOLO',
                      'SAGAMU' = 'SHAGAMU',
                      'AYEDAADE' = 'AIYEDADE',
                      'AYEDIRE' = 'AIYEDIRE',
                      'ILESA-EAST' = 'ILESHA-EAST',
                      'ILESA-WEST' = 'ILESHA-WEST',
                      'ATISBO' = 'ATIGBO',
                      'OBIO-AKPOR' = 'OBIA-AKPOR',
                      'KARIM-LAMIDO' = 'KARIN-LAMIDO',
                      'BUSARI' = 'BURSARI',
                      'YEMALTU-DEBA' = 'YAMALTU-DEBA',
                      'ZAGON-KATAF' = 'ZANGO-KATAF',
                      'GARUN-MALLAM' = 'GARUM-MALLAM',
                      'KMC' = 'KANO-MUNICIPAL',
                      'ZURMI-' = 'ZURMI',
                      'JAMA-ARE' = 'JAMAARE',
                      'JEMA-A' = 'JEMAA',
                      'MAI-ADUA' = 'MAIADUA',
                      'QUA-AN-PAN' = 'QUAAN-PAN',
                      'BASSA-1' = 'BASSA1',
                      'BASSA-2' = 'BASSA2',
                      'OBI-1' = 'OBI1',
                      'OBI-2' = 'OBI2',
                      'NASARAWA-1' = 'NASARAWA1',
                      'NASARAWA-2' = 'NASARAWA2',
                      'IFELODUN-1' = 'IFELODUN1',
                      'IFELODUN-2' = 'IFELODUN2',
                      'IREPODUN-1' = 'IREPODUN1',
                      'IREPODUN-2' = 'IREPODUN2',
                      'SURULERE-1' = 'SURULERE1',
                      'SURULERE-2' = 'SURULERE2',
                      # LGA spelling variants observed in the 2014-2020 routine HF dataset
                      'MBATOLI' = 'MBAITOLI',
                      'AIYEDAADE' = 'AIYEDADE',
                      'ATAKUNMOSA-EAST' = 'ATAKUMOSA-EAST',
                      'ATAKUNMOSA-WEST' = 'ATAKUMOSA-WEST',
                      'EMUOHA' = 'EMOHUA',
                      # '' = '',
                      # '' = '',
                      # BDI
                      'BUJUMBURA-CENTRE' = 'ZONE-CENTRE',
                      'BUJUMBURA-NORD' = 'ZONE-NORD',
                      'BUJUMBURA-SUD' = 'ZONE-SUD',
                      # 'Bukinanyana' = 'Cibitoke',  # Bukinanyana is a commune in Cibitoke province - may now be a new DS
                      # 'Gisuru' = 'Ruyigi',  # Gisuru is a commune in Ruyigi province - may now be a new DS
                      # 'Rutovu' = 'Bururi'  # Rutovu is a commune in Bururi province - may now be a new DS
                      
                      # state names
                      'AKWA-LBOM' = 'AKWA-IBOM',
                      'FEDERAL-CAPITAL-TERRITORY' = 'FCT-ABUJA',
                      'ABUJA-FCT' = 'FCT-ABUJA',
                      'FCT' = 'FCT-ABUJA',
                      'CROSS-RIVER' = 'CRORIVER'

  )
  if(lga_name %in% names(replace_list)){
    lga_name = replace_list[lga_name][[1]]
  }
  return(lga_name)
}


add_number_identification = function(target_names_df, origin_names_df, additional_id_col='State', possible_suffixes=c(1,2,3)){
  if((additional_id_col %in% colnames(target_names_df)) & (additional_id_col %in% colnames(origin_names_df))){
    matching_pattern = paste0("(", paste(possible_suffixes, collapse = "|"), ")$")
    target_names_df$additional_id_col_standardized = toupper(target_names_df[[additional_id_col]])
    origin_names_df$additional_id_col_standardized = toupper(origin_names_df[[additional_id_col]])
    
    # add column to target_names dataframe with the matching name without the suffix, then merge the dataframes
    target_names_df$name_no_suffix = sub(matching_pattern, "", target_names_df$matched_name)
    binding_df = target_names_df %>% rename(name_with_suffix = matched_name, matched_name = name_no_suffix) %>%
      dplyr::select(name_with_suffix, matched_name, additional_id_col_standardized)  #!!sym(additional_id_col))
    origin_names_df = merge(origin_names_df, binding_df, all.x=TRUE)
    origin_names_df$matched_name[!is.na(origin_names_df$name_with_suffix)] = origin_names_df$name_with_suffix[!is.na(origin_names_df$name_with_suffix)]
    
  } else warning('The column for additional identification was not found in one of the dataframes')
  return(origin_names_df)
}




standardize_admin_names_in_df = function(target_names_df, origin_names_df, target_names_col='admin_name', origin_names_col='admin_name', unique_entries_flag=FALSE, additional_id_col='State', possible_suffixes=c(1,2,3)){
  #' given a dataframe with a column containing admin names, update the admin names to match the standardized naming conventions
  #' @param target_names_df dataframe containing a column with the admin names that should be used as the standard
  #' @param origin_names_df dataframe containing a column with admin names that user would like to update to the standard
  #' @param target_names_col name of column in target_names_df where the standard admin names are found
  #' @param origin_names_col name of column in origin_names_df where the to-be-updated admin names are found
  #' @return the origin_names_df data frame, with the admin names replaced with the standardized version

  if(unique_entries_flag){
    duplicates = unique(origin_names_df[[origin_names_col]][which(duplicated(origin_names_df[[origin_names_col]]))])
    if(length(duplicates)>0){
      for(dd in 1:length(duplicates)){
        duplicate_rows = which(origin_names_df[[origin_names_col]] == duplicates[dd])
        # see whether the duplicated names can be distinguished based on value in additional_id_col
        possible_match_rows = as.numeric(sapply(paste0(duplicates[dd], possible_suffixes), grep, target_names_df[[target_names_col]]))
        possible_match_rows = possible_match_rows[!is.na(possible_match_rows)]
        if(length(possible_match_rows) == length(duplicate_rows)){
          origin_additional_ids = origin_names_df[[additional_id_col]][duplicate_rows]
          target_additional_ids = target_names_df[[additional_id_col]][possible_match_rows]
          if(all(origin_additional_ids %in% target_additional_ids)){
            origin_names_df[[origin_names_col]][duplicate_rows] = target_names_df[[target_names_col]][possible_match_rows][match(origin_additional_ids, target_additional_ids)]
          } else stop('PROBLEM ENCOUNTERED: Some of the input admin names are not unique identifiers of the admin. Need to add numeric value to distintuish, determined by admin region, following the system from the standardization base file.')
        } else stop('PROBLEM ENCOUNTERED: Some of the input admin names are not unique identifiers of the admin. Need to add numeric value to distintuish, determined by admin region, following the system from the standardization base file.')
      }
    } 
  } #else if(any(duplicate_names %in% toupper(origin_names_df[[origin_names_col]]))){
    #  stop('PROBLEM ENCOUNTERED: Some of the input admin names are not unique identifiers of the admin. Need to add numeric value to distintuish, determined by admin region, following the system from the standardization base file.')
  #}

  target_names_df$matched_name = sapply(target_names_df[[target_names_col]], create_reference_name_match)
  origin_names_df$matched_name = sapply(origin_names_df[[origin_names_col]], create_reference_name_match)
  
  # check whether any of the un-matched names are duplicates across states if not all names have been matched yet
  if(!all(origin_names_df$matched_name %in% target_names_df$matched_name)){
    if((additional_id_col %in% colnames(target_names_df)) & (additional_id_col %in% colnames(origin_names_df))){
      origin_names_df = add_number_identification(target_names_df=target_names_df, origin_names_df=origin_names_df, additional_id_col=additional_id_col, possible_suffixes=possible_suffixes)
    }
  }
  # check whether all names are now matched
  if(!all(origin_names_df$matched_name %in% target_names_df$matched_name)){
    warning('Some of the source admin names could not be matched with a target admin name')
    View(distinct(origin_names_df[which(!(origin_names_df$matched_name %in% target_names_df$matched_name)), c('matched_name', 'State')]))
    View(distinct(target_names_df[which(!(target_names_df$matched_name %in% origin_names_df$matched_name)), c('matched_name', 'State')]))
    View(target_names_df[,c('matched_name', 'State')])
  } 
  if('data.table' %in% class(origin_names_df)){  # changes how indexing works for columns in dataframe
    origin_names_df = as.data.frame(origin_names_df)
    was_data_table=TRUE
  } else was_data_table = FALSE
  if('data.table' %in% class(target_names_df)){  
    target_names_df = as.data.frame(target_names_df)
  }
  # remove the original admin name column from the origin dataframe, and also the target name column if applicable (for merging)
  origin_names_df = origin_names_df[,-(which(colnames(origin_names_df)==origin_names_col))]
  if(target_names_col %in% colnames(origin_names_df)){
    origin_names_df = origin_names_df[,-which(colnames(origin_names_df)==target_names_col)]
  }
  
  # add the updated admin names to the dataframe, under the original admin column name
  target_names_df =  target_names_df[,c('matched_name', target_names_col)]
  updated_names_df = merge(origin_names_df, target_names_df, all.x=TRUE, by='matched_name')
  colnames(updated_names_df)[colnames(updated_names_df)==target_names_col] = origin_names_col
  
  if(was_data_table) updated_names_df = as.data.table(updated_names_df)
  return(updated_names_df)
}



standardize_state_names_in_df = function(target_names_df, origin_names_df, target_names_col='State', origin_names_col='State'){
  #' given a dataframe with a column containing admin names, update the admin names to match the standardized naming conventions
  #' @param target_names_df dataframe containing a column with the admin names that should be used as the standard
  #' @param origin_names_df dataframe containing a column with admin names that user would like to update to the standard
  #' @param target_names_col name of column in target_names_df where the standard admin names are found
  #' @param origin_names_col name of column in origin_names_df where the to-be-updated admin names are found
  #' @return the origin_names_df data frame, with the admin names replaced with the standardized version
  
  if(any(duplicate_names_state %in% toupper(origin_names_df[[origin_names_col]]))){
    stop('PROBLEM ENCOUNTERED: Some of the input admin names are not unique identifiers of the admin. Need to add numeric value to distintuish, determined by admin region, following the system from the standardization base file.')
  } else{
    target_names_df$matched_name = sapply(target_names_df[[target_names_col]], create_reference_name_match)
    origin_names_df$matched_name = sapply(origin_names_df[[origin_names_col]], create_reference_name_match)
    
    # check whether all names from the origin source are now matched to one of the target names
    if(!all(origin_names_df$matched_name %in% target_names_df$matched_name)){
      warning('Some of the source admin names could not be matched with a target admin name')
      View(distinct(origin_names_df[which(!(origin_names_df$matched_name %in% target_names_df$matched_name)), c('matched_name', 'State')]))
      View(target_names_df[,c('matched_name', 'State')])
    } 
    if('data.table' %in% class(origin_names_df)){  # changes how indexing works for columns in dataframe
      origin_names_df = as.data.frame(origin_names_df)
      was_data_table=TRUE
    } else was_data_table = FALSE
    if('data.table' %in% class(target_names_df)){  
      target_names_df = as.data.frame(target_names_df)
    }
    # remove the original admin name column from the origin dataframe, and also the target name column if applicable (for merging)
    origin_names_df = origin_names_df[,-(which(colnames(origin_names_df)==origin_names_col))]
    if(target_names_col %in% colnames(origin_names_df)){
      origin_names_df = origin_names_df[,-which(colnames(origin_names_df)==target_names_col)]
    }
    
    # add the updated admin names to the dataframe, under the original admin column name
    target_names_df =  distinct(target_names_df[,c('matched_name', target_names_col)])
    updated_names_df = merge(origin_names_df, target_names_df, all.x=TRUE, by='matched_name') # only include States in the origin df (i.e., don't include empty placeholders if they weren't in the original)
    colnames(updated_names_df)[colnames(updated_names_df)==target_names_col] = origin_names_col
    
    if(was_data_table) updated_names_df = as.data.table(updated_names_df)
    return(updated_names_df)
  }
}





standardize_admin_names_in_vector = function(target_names, origin_names){
  #' given a vector of admin names, update the admin names to match the standardized naming conventions provided in target_names
  #' @param target_names set of admin names that should be used as the standard
  #' @param origin_names current set of admin names that user would like to update to the standard
  #' @return a vector with updated admin names, in the same order as origin_names
  
  
  if(any(duplicate_names %in% toupper(origin_names))){
    stop('PROBLEM ENCOUNTERED: Some of the input admin names are not unique identifiers of the admin. Need to add numeric value to distintuish, determined by admin region, following the system from the standardization base file.')
  } else{
    target_matched_name = sapply(target_names, create_reference_name_match)
    origin_matched_name = sapply(origin_names, create_reference_name_match)
    
    # check whether all names from the origin source are now matched to one of the target names
    if(!all(origin_matched_name %in% target_matched_name)){
      warning('Some of the source admin names could not be matched with a target admin name')
    } 
    
    # create dataframes to merge
    target_df = data.frame('matched_name'=target_matched_name, 'target_name'=target_names)
    origin_df = data.frame('matched_name'=origin_matched_name, 'original_name'=origin_names)
  
    # add the updated admin names to the dataframe, under the original admin column name
    updated_names_df = merge(origin_df, target_df, all.x=TRUE, by='matched_name', sort=FALSE)
    if(all(updated_names_df$original_name == origin_names)){
      return(updated_names_df$target_name)
    } else{
      warning('merge did not maintain original order, need to fix standardize_admin_names_in_vector function or use standardize_admin_names_in_df')
      return(NA)
    }
  }
}



