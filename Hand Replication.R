

probe_1 <- exprs(eset_mini)[1,] %>% #subset first row of values from the expression set matrix
  cbind(exprs(eset_mini)[1,] %>% names()) %>% #combine the values with their coln names into a new matrix
  as_tibble %>% rename(value = ".", label = V2) %>% #rename the colns
  left_join(pData(eset_mini)) #join with the phenotype data

probe_1$value <- probe_1$value %>% as.numeric()


probe_1_groups <- probe_1 %>% group_by(group) %>% summarise(average = mean(value)) 

(probe_1_groups %>% filter(str_detect(probe_1_groups$group, "test")) %>% summarise(avg = mean(average))
  -
probe_1_groups %>% filter(str_detect(probe_1_groups$group, "control")) %>% summarise(avg = mean(average)))
     