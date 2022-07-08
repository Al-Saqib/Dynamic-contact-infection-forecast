
#' @title NA omit for list
#'
#' @description Generalization of NA omit for lists.
#'
#' @param y list: Patient pathways with tests (when a patient become positive)
#' 
#' @return a list with NA elements removed

na.omit.list <- function(y) {return(y[!sapply(y, function(x) all(is.na(x)))])}


#' @title Pre-proccess function for pathways 
#'
#' @description Takes pathways in a long format and splits based on unique patients.
#'
#' @param pathways dataframe: long format dataframe, detailing date-location positions for each patient
#' 
#' @return a list with each element corresponding to a single patient

getPathList = function(pathways){
  # Split into list of patients
  path.l <- split(pathways, f = pathways$Ptnumber)
  path.l = lapply(seq_along(path.l), function(n,list.df){
    x = list.df[[n]]
    x = x[order(x$Date),]
    
    if(any(is.na(x$Location))){x = x[-which(is.na(x$Location)),]}
    ifelse(nrow(x)>0,return(x),return(NA))
    
  },list.df = path.l)
  
  path.l = na.omit.list(path.l)
  names(path.l) = bind_rows(lapply(path.l, function(x){return(x[1,])}))$Ptnumber
  return(path.l)
}


#' @title Get time bounds for windowing
#'
#' @description Given a path list, the function computes all windows of specified lengths.
#'
#' @param path.l list: processed list of pathways.
#' @param feature_n integer: length of time used for feature generation.
#' @param prediction_n integer: length of time used for forecasting/prediction horizon.
#' 
#' @return a list of time bounds with feature periods, and prediction periods.

getTimeBounds = function(path.l,feature_n,prediction_n){
  t = do.call("c",lapply(path.l, function(x){return(x$Date)}))
  t_seq = seq.Date(min(t),max(t),by = "day")
  
  t_bounds_l = na.omit.list(lapply(seq_along(t_seq), function(n,t_seq,feature_n,prediction_n){
    # Check not out of bounds
    if((n+feature_n+prediction_n) > length(t_seq)){return(NA)}
    
    # Create boundaries
    t_feat_start = t_seq[n]
    t_feat_end = t_seq[n+feature_n-1]
    t_pred_start = t_seq[n+feature_n]
    t_pred_end = t_seq[n+feature_n+prediction_n]
    return(list("t_feat_start" = t_feat_start,"t_feat_end"=t_feat_end,
                "t_pred_start"=t_pred_start,"t_pred_end"=t_pred_end))
    
  },t_seq,feature_n,prediction_n))
  
  return(t_bounds_l)
}


#' @title Infectious period for known infected individuals
#'
#' @description For individuals with postive test results, we compute the date range 
#' they were likely infectious over.
#'
#' @param path.l list: processed list of pathways.
#' @param dateRange date: a vector of dates to consider.
#' @param incubation.t integer: incubation period for pathogen in days.
#' @param infectious.t integer: time after symptoms a patient remains infectious.
#' 
#' @return a data frame of individuals with an infectious state.

getInfectiousPtDf = function(path.l.sub,dateRange,incubation.t = 14, infectious.t = 10){
  
  infectious.pts = bind_rows(lapply(path.l.sub, function(x,dateRange){
    
    # If NA then no positve test : non-infectious
    if(any(is.na(x$Positive_test_date))){
      return(data.frame("Pt" = x$Ptnumber[1],
                        "Infectious" = 0,
                        "Positive_test_date" = x$Positive_test_date[1]))
    }else{
      # Check if positve test within period
      infectiousDateRange = seq.Date((x$Positive_test_date[1] - incubation.t),(x$Positive_test_date[1] + infectious.t),by = "days") #<- can introduce bias, need to make sure doesn't include future 'unknown' cases
      
      if(length(intersect(infectiousDateRange,dateRange))>0 & x$Positive_test_date[1] <= max(dateRange)){
        return(data.frame("Pt" = x$Ptnumber[1],
                          "Infectious" = 1,
                          "Positive_test_date" = x$Positive_test_date[1]))
      }else{ #  Else was positve but non infectious in window period
        return(data.frame("Pt" = x$Ptnumber[1],
                          "Infectious" = 0,
                          "Positive_test_date" = x$Positive_test_date[1]))
      }
    }
  },dateRange))
  
  return(infectious.pts)
}


#' @title Infected degree computation
#'
#' @description For a contact graph, we compute the infected degree of each node
#'
#' @param nodes dataframe: nodes of the graph and their infectious state.
#' @param graph igraph: the contact network.
#' 
#' @return a data frame of infected degree for each node.

getInfectedDegree = function(nodes,graph){
  nodeRanks = bind_rows(lapply(split(nodes,nodes$Pt), function(x,graph,nodes){
    if(as.character(x$Pt[1]) %in% names(V(graph))){
      
      neighbors = names(neighbors(graph,as.character(x$Pt[1])))
      contacts = filter(nodes, Pt %in% neighbors)
      
      ret.df = data.frame("Node" = x$Pt[1],
                          "Infected degree" = sum(contacts$Infectious))
      
      return(ret.df)
    }else{
      
      ret.df = data.frame("Node" = x$Pt[1],
                          "Infected degree" = 0)
      
      return(ret.df)
    }
  },graph,nodes))
  
  
  return(nodeRanks)
}

#' @title Infected degree centrality computation
#'
#' @description For a contact graph, we compute the infected degree centrality of each node
#'
#' @param nodes dataframe: nodes of the graph and their infectious state.
#' @param graph igraph: the contact network.
#' 
#' @return a data frame of infected degree centrality for each node.

getInfectedDegreeCentral = function(nodes,edges){
  
  # Take only infectious nodes
  inf_nodes = filter(nodes, Infectious == 1)
  
  # Unfactor edges for dplyr operations
  edges$source = as.character(edges$source)
  edges$target = as.character(edges$target)
  
  nodeRanks = bind_rows(lapply(split(nodes,nodes$Pt), function(v,edges,nodes){
    v = as.character(v$Pt[1])
    #print(v)
    
    # Take local edges
    local_edges = filter(edges, source == v | target == v)
    local_edges$source[which(local_edges$source == v)] = "v"
    local_edges$target[which(local_edges$target == v)] = "v"
    
    # Get local edges which were with infectious nodes
    local_edges = filter(local_edges, source %in% inf_nodes$Pt | target %in% inf_nodes$Pt)
    
    ret.df = data.frame("Node" = v,"InfectedDegCentral" = sum(local_edges$weight))
    
    return(ret.df)
    
  },edges,nodes))
  
  # Normalise by N-1 
  nodeRanks$InfectedDegCentral = nodeRanks$InfectedDegCentral / (nrow(nodeRanks) -1)
  
  return(nodeRanks)
}


#' @title Infected closeness centrality computation
#'
#' @description For a contact graph, we compute the infected degree of each patient
#'
#' @param nodes dataframe: nodes of the graph and their infectious state.
#' @param graph igraph: the contact network.
#' @param  normalise logical: true/false to normalize the closeness centrality by number of nodes.
#' 
#' @return a data frame of infected closeness centrality for each node.

getInfectedCloseness = function(nodes,graph,normalise=TRUE){
  
  # Take only infectious nodes
  inf_nodes = filter(nodes, Infectious == 1)
  
  nodeRanks = bind_rows(lapply(split(nodes,nodes$Pt), function(x,graph,inf_nodes,normalise){
    
    # Check vertex in graph
    if(as.character(x$Pt[1]) %in% names(V(graph))){
      
      # Get Shortest path dataframe
      sp = shortest.paths(graph,as.character(x$Pt[1]))
      df.sp = as.data.frame(t(sp))
      colnames(df.sp) = "sp.weight"
      df.sp$pt = row.names(df.sp)
      
      # Take only finite values (reachable nodes)
      df.sp = df.sp[is.finite(df.sp$sp.weight),]
      
      # Filter only shortest pathways to infectious nodes
      df.sp.infected = filter(df.sp, pt %in% inf_nodes$Pt)
      
      # Infected closenss centrality
      i_cc = (nrow(df.sp.infected) -1)/ sum(df.sp.infected$sp.weight)
      
      # Infite created from n-1 / 0
      if(!is.finite(i_cc)){
        i_cc = 0
      }
      
      # Normalise to number of nodes-1 in connected part 
      if(normalise){
        s = (nrow(df.sp.infected) -1)/length(V(graph))
        i_cc = i_cc *s
      }
      
    }else{
      
      i_cc = 0
    }
    
    ret.df = data.frame("Node" = x$Pt[1],
                        "Infected closenss" = i_cc)
    return(ret.df)
    
  },graph,inf_nodes,normalise))
  
  
  return(nodeRanks)
}


#' @title Compute consecutive days in hospital
#'
#' @description Given a vector of times in hospital, count the number of consectuive days in time from the most recent time point.
#'
#' @param t date: vector of dates.
#' 
#' @return an integer for the consecutive days in hospital from the most recent time point.

getLoSConsec = function(t){
  if(length(t)>1){
    t_dif = t[-1] - t[-length(t)]
    idx = ifelse(any(which(t_dif > 1)),max(which(t_dif > 1)),1)
    return(length(t_dif[(idx):length(t_dif)]))
  }else{
    return(0)
  }
}


#' @title Compute length of stay information
#'
#' @description Function to compute length of stays
#' 
#' @param pathway.list list: processed pathways with each individual one element.
#' 
#' @return a dataframe with variable information for each individual.


getLoS = function(pathway.list){
  bind_rows(lapply(pathway.list, function(x){
    return(data.frame("Node" = x$Ptnumber[1],
                      "LoS" = length(unique(x$Date)),
                      "LoS_Consec" = getLoSConsec(unique(x$Date)), # <- prior to last day how many days in hospital?
                      "lastDate" = max(x$Date)))
  }))
}

#' @title Filter pathway dates
#'
#' @description Function to filter pathway list
#' 
#' @param pat.df dataframe: dataframe for pathways
#' @param dates date: vector of dates to include
#' 
#' @return a date filter dataframe.

filterPathDates = function(path.df,dates){
  path.df$Date = as.Date(path.df$Date)
  path.df = filter(path.df, Date %in% dates)
  if(nrow(path.df)>0){
    return(path.df)
  }else{
    return(NA)
  }
}


#' @title Compute contact network
#'
#' @description Function to  generate contact network from a list of pathways
#' 
#' @param path.l.sub list: processed list of pathways.
#' @param t.idx integer: column in pathway data frame holding the time.
#' @param loc.idx integer: column in pathway data frame holding the location.
#' 
#' @return a dateframe of contact events.

getContactNet = function(path.l.sub,t.idx,loc.idx){
  
  # Extract columns
  pathways.sub = bind_rows(path.l.sub)[,c(1,t.idx,loc.idx)]
  
  print(paste("Computing contact ",colnames(pathways.sub)[3],
              " network using time as:",colnames(pathways.sub)[2],sep = ""))
  
  colnames(pathways.sub) = c("Ptnumber","t","location")
  

  # Create time-space id
  pathways.sub$loc_t = paste(pathways.sub$location,pathways.sub$t,sep = "-")
  
  # Split list based on time-space id
  contact.l = split(pathways.sub, f = pathways.sub$loc_t)
  
  # Extract contacts from the same time-space id
  contacts = bind_rows(lapply(contact.l, function(x){
    c = as.data.frame(expand.grid.alt((x$Ptnumber),(x$Ptnumber)))
    c = c[!duplicated(data.frame(t(apply(c,1,sort)))),]
    c$weight = 1
    return(c)
  }))
  
  colnames(contacts)[1:3] = c("source","target","weight")
  contacts = filter(contacts, !(source == target))
  
  contacts = contacts %>% group_by(source,target) %>% summarise(weight=sum(weight))
  
  return(contacts)
}


#' @title Expand grid
#'
#' @description Fast version of expand grid
#' 
#' @param seq1 char/int: vector
#' @param seq2 char/int: vector
#' 
#' @return dataframe of all combinations for seq1 and seq2

expand.grid.alt <- function(seq1,seq2) {
  cbind(rep.int(seq1, length(seq2)),
        c(t(matrix(rep.int(seq2, length(seq1)), nrow=length(seq2)))))
}



#' @title Function to get contexutal variables
#'
#' @description given background contextual variables with time stamps, compute the local variables for a given date range.
#' 
#' @param pathway.list list: list of pathway dataframes
#' @param dateRange dates: date vector to filter by
#' @param contextualVars dataframe: time stamped records of variables by days
#' 
#' @return dataframe of contextual variables during the date range.

getPtEnvi  = function(pathway.list,dateRange,contextualVars){
  data.frame("Node" = names(pathway.list),
             "BackgroundCases" = sum(filter(contextualVars, Date %in% dateRange)$BackgroundCases),
             "BackgroundHospPatient" = sum(filter(contextualVars, Date %in% dateRange)$BackgroundHospPatient)
  )
}


#' @title Pre-process rolling window data
#'
#' @description The function preProRollingWind() is a wrapper function that is 
#' applied over a sliding time window of length feature_n. In summary, the function: 
#' (i) splits the data into windows, (ii) constructs a contact network, and then 
#' centrality of each patient across it, (iii) derives the background contextual 
#' variable for the window, and (iii) joins the patient statistics within the window 
#' with static variables.
#'
#' @param pathwaysWithTests dataframe: Patient pathways with tests (when a patient become positive)
#' @param staticVars dataframe: Static variables (i.e. age, gender, ...)
#' @param contextualVars dataframe: Background contextual variables (hospital infection numbers)
#' @param feature_n integer: Time window size to compute variables over
#' @param prediction_n integer:  Forecasting horizon
#'
#' @return A list of dataframes containing windowed datasets

preProRollingWind = function(pathwaysWithTests,staticVars,contextualVars,feature_n,prediction_n){
  
  # -------
  ## Enforce dates
  pathwaysWithTests$Date = as.Date(pathwaysWithTests$Date)
  pathwaysWithTests$Positive_test_date = as.Date(pathwaysWithTests$Positive_test_date)
  contextualVars$Date = as.Date(contextualVars$Date)
  
  # -------
  ## Convert pathways into a list
  path.l = getPathList(pathwaysWithTests)
  pathways_pre_merge = bind_rows(path.l)
  
  # ----
  ##  Window time bounds
  t_bounds_l = getTimeBounds(path.l,feature_n,prediction_n)
  
  
  # ----
  ## Extract variable data frames by running over time windows
  total <- length(t_bounds_l)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  ret.df.l = lapply(seq_along(t_bounds_l), function(n,t_bounds_l,pb){
    
    print("Generating windowed data")
    setTxtProgressBar(pb, n)
    t_b = t_bounds_l[[n]]
    
    # ----
    # Filter pathways list to feature period
    dateRange = seq.Date(as.Date(t_b$t_feat_start),as.Date(t_b$t_feat_end),by = "days")
    path.l.sub = filter(pathways_pre_merge, Date %in% dateRange)
    path.l.sub = split(path.l.sub,path.l.sub$Ptnumber)
    
    
    # ----
    # Infectious patient dataframe, if -14 CollectionDt +10 in the 2-week period.
    infectious.pts = getInfectiousPtDf(path.l.sub,dateRange,incubation.t = 14,infectious.t = 10)
    
    # ----
    # Construct contact networks
    contacts = getContactNet(path.l.sub,t.idx = 3,loc.idx = 2)
    
    # ----
    ## Network metrics
    nodes = infectious.pts
    edges = contacts
    g <- graph_from_data_frame(edges,directed = F)
    infDeg = getInfectedDegree(nodes,graph_from_data_frame(edges,directed = F))
    #infDegCentral = getInfectedDegreeCentral(nodes,edges)
    #infCloseness = getInfectedCloseness(nodes,graph_from_data_frame(edges,directed = F))
    
    # ----
    ## Extract Length of stay
    ptLoS = getLoS(path.l.sub)
    
    # ----
    ## Extract background enviromental variables: 
    ptEnvi = getPtEnvi(path.l.sub, 
                       dateRange,
                       contextualVars)
    
    # ----
    ## 2.4 Combine data
    colnames(nodes)[1] = "Node"
    colnames(staticVars)[1] = "Node"
    
    stat.df = list(infDeg,
                   #infDegCentral,
                   #infCloseness,
                   ptLoS,
                   staticVars,
                   ptEnvi,
                   nodes) %>% reduce(full_join, by = c("Node"))
    
    
    # ----
    ## Filter only patients who are at  risk (remove vaccinations + recent COVID-19 positive)
    stat.df = filter(stat.df, !(Positive_test_date <= max(dateRange)) | is.na(Positive_test_date))
    
    # ----
    ## Filter for only patients who are currently in hospital
    stat.df = filter(stat.df, lastDate == max(dateRange))
    
    # ----
    ## Add time ID
    stat.df$t_feat_start = t_b$t_feat_start
    stat.df$t_feat_end = t_b$t_feat_end
    
    return(stat.df)
    
  },t_bounds_l,pb)
  
  return(bind_rows(ret.df.l))
  
}




