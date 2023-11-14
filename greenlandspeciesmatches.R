us_85 <- sort(c(
"Amblyraja radiata", "Ammodytes dubius", "Ammodytes marinus", "Anarhichas denticulatus", "Anarhichas lupus" ,
"Anarhichas minor", "Anotopterus pharao", "Arctogadus glacialis", "Arctozenus risso", "Argentina silus",
"Artediellus atlanticus", "Artediellus uncinatus", "Aspidophoroides monopterygius", "Bathylagus euryops", "Bathyraja spinicauda" ,
"Benthosema glaciale", "Boreogadus saida", "Brosme brosme", "Careproctus micropus", "Careproctus reinhardti" ,
"Centroscyllium fabricii", "Chauliodus sloani", "Clupea harengus", "Coryphaenoides rupestris", "Cottunculus microps",
"Cryptopsaras couesii", "Cyclopterus lumpus", "Dipturus batis", "Eumicrotremus spinosus", "Eutrigla gurnardus" ,
"Gadus macrocephalus", "Gadus morhua", "Gaidropsarus argentatus", "Gaidropsarus ensis", "Glyptocephalus cynoglossus" ,
"Gymnelus retrodorsalis", "Gymnelus viridis", "Gymnocanthus tricuspis", "Hippoglossoides platessoides", "Hippoglossus hippoglossus",
"Icelus bicornis", "Icelus spatula", "Lampanyctus macdonaldi", "Lepidion eques", "Leptagonus decagonus" ,
"Leptoclinus maculatus", "Liparis fabricii", "Liparis liparis", "Liparis montagui", "Liparis tunicatus",
"Lophius piscatorius", "Lumpenus lampretaeformis", "Lycenchelys kolthoffi", "Lycenchelys sarsii", "Lycodes esmarkii" ,
"Lycodes pallidus", "Lycodes vahlii", "Macrourus berglax", "Magnisudis atlantica", "Mallotus villosus",
"Melanogrammus aeglefinus", "Merlangius merlangus", "Micromesistius poutassou", "Microstomus kitt", "Molva dypterygia" ,
"Molva molva", "Myoxocephalus scorpioides", "Myoxocephalus scorpius", "Myxine glutinosa", "Nemichthys scolopaceus" ,
"Notacanthus chemnitzii", "Notoscopelus kroyeri", "Petromyzon marinus", "Pleuronectes platessa", "Pollachius virens",
"Rajella fyllae", "Rajella lintea", "Reinhardtius hippoglossoides", "Sebastes fasciatus", "Sebastes mentella",
"Sebastes norvegicus", "Sebastes viviparus", "Serrivomer beanii", "Somniosus microcephalus", "Squalus acanthias",
"Synaphobranchus kaupii", "Triglops murrayi", "Triglops nybelini", "Triglops pingelii", "Trisopterus esmarkii"
)) 

us_twothirds <- sort(c(
   "Amblyraja radiata",            "Anarhichas denticulatus",     
 "Anarhichas lupus",             "Anarhichas minor",            
 "Arctozenus risso",             "Argentina silus",             
 "Artediellus atlanticus",       "Artediellus uncinatus",       
 "Aspidophoroides monopterygius","Boreogadus saida",            
 "Brosme brosme",                "Careproctus reinhardti",      
 "Cottunculus microps",          "Cyclopterus lumpus",          
 "Eumicrotremus spinosus",       "Gadus macrocephalus",         
 "Gadus morhua",                 "Gaidropsarus argentatus",     
 "Gaidropsarus ensis",           "Glyptocephalus cynoglossus",  
 "Gymnelus retrodorsalis",       "Hippoglossoides platessoides",
 "Hippoglossus hippoglossus",    "Lepidion eques",              
 "Leptagonus decagonus",         "Leptoclinus maculatus",       
 "Lumpenus lampretaeformis",     "Lycodes esmarkii",            
 "Lycodes vahlii",               "Macrourus berglax",           
 "Mallotus villosus",            "Melanogrammus aeglefinus",    
 "Micromesistius poutassou",     "Microstomus kitt",            
 "Molva dypterygia",             "Molva molva",                 
 "Myoxocephalus scorpius",       "Myxine glutinosa",            
 "Pollachius virens",            "Rajella fyllae",              
 "Reinhardtius hippoglossoides", "Sebastes mentella",           
 "Sebastes norvegicus",          "Sebastes viviparus",          
 "Somniosus microcephalus",      "Squalus acanthias",           
 "Triglops murrayi"))
 
 
 them <- sort(str_to_sentence(c(
           "ANARHICHAS DENTICULATUS",
           "ANARHICHAS LUPUS",
           "ANARHICHAS MINOR",
           "ARGENTINA SILUS",
           "BROSME BROSME",
           "CAREPROCTUS REINHARDTI",
           "CYCLOPTERUS LUMPUS",
           "EUMICROTREMUS SPINOSUS",
           "GADUS MORHUA",
           "GLYPTOCEPHALUS CYNOGLOSSUS",
           "HIPPOGLOSSOIDES PLATESSOIDES",
           "HIPPOGLOSSUS HIPPOGLOSSUS",
           "LEPTAGONUS DECAGONUS",
           "LEPTOCLINUS MACULATUS",
           "LYCODES ESMARKII",
           "LYCODES VAHLII",
           "MELANOGRAMMUS AEGLEFINUS",
           "MICROSTOMUS KITT",
           "MOLVA DIPTERYGIA",
           "MOLVA MOLVA",
           "MYXINE GLUTINOSA",
           "POLLACHIUS VIRENS",
           "PETROMYZON MARINUS",
           "RAJA FYLLAE",
           "RAJA RADIATA",
           "SEBASTES MARINUS",
           "SEBASTES MENTELLA",
           "SQUALUS ACANTHIAS",
           "TRIGLOPS MURRAYI",
           "TRISOPTERUS ESMARKI"
           )))
 
 setdiff(us_85,them)
 # Amblyraja radiata     
 # Ammodytes dubius      
 # Ammodytes marinus     
 # Anotopterus pharao    
 # Arctogadus glacialis  
 # Arctozenus risso      
 # Artediellus atlanticus
 # Artediellus uncinatus 
 # Aspidophoroides monopterygius
 # Bathylagus euryops    
 # Bathyraja spinicauda  
 # Benthosema glaciale   
 # Boreogadus saida      
 # Careproctus micropus  
 # Centroscyllium fabricii      
 # Chauliodus sloani     
 # Clupea harengus       
 # Coryphaenoides rupestris     
 # Cottunculus microps   
 # Cryptopsaras couesii  
 # Dipturus batis        
 # Eutrigla gurnardus    
 # Gadus macrocephalus   
 # Gaidropsarus argentatus      
 # Gaidropsarus ensis    
 # Gymnelus retrodorsalis
 # Gymnelus viridis      
 # Gymnocanthus tricuspis
 # Icelus bicornis       
 # Icelus spatula        
 # Lampanyctus macdonaldi
 # Lepidion eques        
 # Liparis fabricii      
 # Liparis liparis       
 # Liparis montagui      
 # Liparis tunicatus     
 # Lophius piscatorius   
 # Lumpenus lampretaeformis     
 # Lycenchelys kolthoffi 
 # Lycenchelys sarsii    
 # Lycodes pallidus      
 # Macrourus berglax     
 # Magnisudis atlantica  
 # Mallotus villosus     
 # Merlangius merlangus  
 # Micromesistius poutassou     
 # Molva dypterygia      
 # Myoxocephalus scorpioides    
 # Myoxocephalus scorpius
 # Nemichthys scolopaceus
 # Notacanthus chemnitzii
 # Notoscopelus kroyeri  
 # Pleuronectes platessa 
 # Rajella fyllae        
 # Rajella lintea        
 # Reinhardtius hippoglossoides 
 # Sebastes fasciatus    
 # Sebastes norvegicus   
 # Sebastes viviparus    
 # Serrivomer beanii     
 # Somniosus microcephalus      
 # Synaphobranchus kaupii
 # Triglops nybelini     
 # Triglops pingelii     
 # Trisopterus esmarkii   
 
 
 setdiff(us_twothirds,them)
 #MISMATCHES
#  Amblyraja radiata: we name Amblyraja radiata instead of Raja radiata, but included in both sensitivity analyses          
#  Arctozenus risso:              
#  Artediellus atlanticus       
#  Artediellus uncinatus        
#  Aspidophoroides monopterygius
#  Boreogadus saida             
#  Cottunculus microps          
#  Gadus macrocephalus          
#  Gaidropsarus argentatus      
#  Gaidropsarus ensis           
#  Gymnelus retrodorsalis       
#  Lepidion eques               
#  Lumpenus lampretaeformis     
#  Macrourus berglax            
#  Mallotus villosus            
#  Micromesistius poutassou     
#  Molva dypterygia             
#  Myoxocephalus scorpius       
#  Rajella fyllae               
#  Reinhardtius hippoglossoides 
#  Sebastes norvegicus: SEBASTES MARINUS: we name Sebastes norvegicus, but included in both sensitivity analyses         
#  Sebastes viviparus           
#  Somniosus microcephalus      

 
#pull in 3 sets of dissimilarity
distances_dissimilarities_allyears_onethirdyears_excluded <- readRDS(here::here("output","dissimilarities","distances_dissimilarities_allyears_onethirdyears_excluded.rds"))
distances_dissimilarities_allyears_15perc_excluded <- readRDS(here::here("output","dissimilarities","distances_dissimilarities_allyears_15perc_excluded.rds"))
distances_dissimilarities_allyears <- readRDS(here::here("output","dissimilarities","distances_dissimilarities_allyears.rds"))

#Plot Greenland only together
#Bray
bray_curtis_greenland_dissim_sensitivity_comparison <- ggplot() +
  geom_point(data = distances_dissimilarities_allyears[survey_unit == "GRL-DE",],aes(x =year,y =bray_curtis_dissimilarity_total_mean), color = "turquoise", size = 4, shape = 1) +
  geom_point(data = distances_dissimilarities_allyears_15perc_excluded[survey_unit == "GRL-DE",],aes(x =year,y =bray_curtis_dissimilarity_total_mean), color = "coral2", size = 4, shape = 2) +
  geom_point(data = distances_dissimilarities_allyears_onethirdyears_excluded[survey_unit == "GRL-DE",],aes(x =year,y =bray_curtis_dissimilarity_total_mean), color = "limegreen", size = 4, shape = 3) +
  labs(x = "Year", y = "Bray Curtis total dissimilarity") +
  theme_classic()

#Jaccard
jaccard_greenland_dissim_sensitivity_comparison <- ggplot() +
  geom_point(data = distances_dissimilarities_allyears[survey_unit == "GRL-DE",],aes(x =year,y =jaccard_dissimilarity_total_mean), color = "turquoise", size = 4, shape = 1) +
  geom_point(data = distances_dissimilarities_allyears_15perc_excluded[survey_unit == "GRL-DE",],aes(x =year,y =jaccard_dissimilarity_total_mean), color = "coral2", size = 4, shape = 2) +
  geom_point(data = distances_dissimilarities_allyears_onethirdyears_excluded[survey_unit == "GRL-DE",],aes(x =year,y =jaccard_dissimilarity_total_mean), color = "limegreen", size = 4, shape = 3) +
  labs(x = "Year", y = "Jaccard total dissimilarity") +
  theme_classic()


greenland_dissim_sensitivity_comparison <- cowplot::plot_grid(bray_curtis_greenland_dissim_sensitivity_comparison,
                                                              jaccard_greenland_dissim_sensitivity_comparison,
                                                              nrow = 1)
  
ggsave(greenland_dissim_sensitivity_comparison, filename = "greenland_dissim_sensitivity_comparison.jpg", path = here::here("figures"), height = 8, width = 14)
