###Palette for Plotting
#Palette for plotting all 42 survey units

library(data.table)

survey_unit.list <- c(
    "AI", 
    "BITS-1",
    "BITS-4", 
    "CHL",
    "DFO-NF", 
    "DFO-QCS",
    "EBS", 
    "EVHOE",
    "FR-CGFS", 
    "GMEX-Fall",
    "GMEX-Summer", 
    "GOA",
    "GRL-DE", 
    "GSL-N",
    "GSL-S", 
    "ICE-GFS",
    "IE-IGFS", 
    "MEDITS",
    "NAM", 
    "NEUS-Fall",
    "NEUS-Spring", 
    "NIGFS-1",
    "NIGFS-4", 
    "Nor-BTS-3",
    "NS-IBTS-1", 
    "NS-IBTS-3",
    "NZ-CHAT", 
    "NZ-ECSI",
    "Z-SUBA", 
    "NZ-WCSI",
    "PT-IBTS", 
    "ROCKALL",
    "S-GEORG", 
    "SCS-Summer",
    "SEUS-fall", 
    "SEUS-spring",
    "SEUS-summer", 
    "SWC-IBTS-1",
    "SWC-IBTS-4", 
    "WCANN",
    "ZAF-ATL", 
    "ZAF-IND")
  
  
  

palette_42 <- c(
  "#5A5156", #AI
  "#DF00DB", #BITS-1
  "#DB8EDA", #BITS-4
  "#F6222E", #CHL
  "#F8A19F", #DFO-NF
  "#16FF32", #DFO-QCS
  "#325A9B", #EBS
  "#3283FE", #EVHOE
  "#FEAF16", #FR-CGFS
  "#fccb6d", #GMEX-Fall
  "#1C8356", #GMEX-Summer
  "#C4451C", #GOA
  "#85660D", #GRL-DE
  "#B0009F", #GSL-N
  "#BF79B8", #GSL-S
  "#1CBE4F", #ICE-GFS
  "#782AB6", #IE-IGFS
  "#90AD1C", #MEDITS
  "#6B003A", #NAM
  "#A75B00", #NEUS-Fall
  "#E3B072", #NEUS-Spring
  "#02E8B6", #NIGFS-1
  "#97E7D5", #NIGFS-4
  "#B00068", #Nor-BTS-3
  "#00B9E3", #NS-IBTS-1
  "#95E2F4", #NS-IBTS-3
  "#B3CE73", #NZ-CHAT
  "#689500", #NZ-ECSI
  "#364d02",#NZ-SUBA
  "#AAF400", #NZ-WCSI
  "#AA0DFE", #PT-IBTS
  "#7f9eb8", #ROCKALL
  "#FA0087", #S-GEORG
  "#DEA0FD", #SCS-Summer
  "#FCEF88", #SEUS-fall
  "#A59405", #SEUS-spring
  "#FCE100", #SEUS-summer
  "#544563", #SWC-IBTS-1
  "#a37fc7", #SWC-IBTS-4
  "#C075A6", #WCANN
  "#BDCDFF", #ZAF-ATL
  "#003EFF"  #ZAF-IND
)

color_link <- data.table(survey_unit = survey_unit.list,hex = palette_42)

name_helper <- data.table(Survey_Name_Season = 
                            c(
                              "Aleutian Islands",
                              "Baltic Sea Q1",
                              "Baltic Sea Q4",
                              "Chile",
                              "Newfoundland",
                              "Queen Charlotte Sound",
                              "Eastern Bering Sea",
                              "Bay of Biscay",
                              "English Channel",
                              "Gulf of Mexico Summer",
                              "Gulf of Alaska",
                              "Greenland",
                              "N Gulf of St. Lawrence",
                              "S Gulf of St. Lawrence",
                              "Iceland",
                              "Irish Sea",
                              "Mediterranean",
                              "Namibia",
                              "NE US Fall",
                              "NE US Spring",
                              "N Ireland Q1",
                              "N Ireland Q4",
                              "Barents Sea Norway Q3",
                              "N Sea Q1",
                              "N Sea Q3",
                              "Chatham Rise NZ",
                              "E Coast S Island NZ",
                              "W Coast S Island NZ",
                              "Portugal",
                              "S Georgia",
                              "Scotian Shelf Summer",
                              "SE US Fall",
                              "SE US Spring",
                              "SE US Summer",
                              "W Coast US",
                              "Atlantic Ocean ZA",
                              "Indian Ocean ZA",
                              "Rockall Plateau",
                              "Scotland Shelf Sea Q1",
                              "Scotland Shelf Sea Q4",
                              "Falkland Islands",
                              "Gulf of Mexico Fall",
                              "Barents Sea Norway Q1",
                              "Sub-Antarctic NZ",
                              "Scotian Shelf Spring"),
                          survey_unit = c(
                            "AI",        
                            "BITS-1",    
                            "BITS-4",    
                            "CHL",       
                            "DFO-NF",    
                            "DFO-QCS",   
                            "EBS",       
                            "EVHOE",     
                            "FR-CGFS",   
                            "GMEX-Summer",
                            "GOA",       
                            "GRL-DE",    
                            "GSL-N",     
                            "GSL-S",     
                            "ICE-GFS",   
                            "IE-IGFS",   
                            "MEDITS",    
                            "NAM",       
                            "NEUS-Fall", 
                            "NEUS-Spring",
                            "NIGFS-1",   
                            "NIGFS-4",   
                            "Nor-BTS-3", 
                            "NS-IBTS-1", 
                            "NS-IBTS-3", 
                            "NZ-CHAT",   
                            "NZ-ECSI",   
                            "NZ-WCSI",   
                            "PT-IBTS",   
                            "S-GEORG",   
                            "SCS-SUMMER",
                            "SEUS-fall", 
                            "SEUS-spring",
                            "SEUS-summer",
                            "WCANN",     
                            "ZAF-ATL",   
                            "ZAF-IND",
                            "ROCKALL",
                            "SWC-IBTS-1",
                            "SWC-IBTS-4",
                            "FALK",
                            "GMEX-Fall",
                            "Nor-BTS-1",
                            "NZ-SUBA",
                            "SCS-SPRING"
                          ))


color_link <- name_helper[color_link, on = "survey_unit"]

rm(name_helper, palette_42, survey_unit.list)
