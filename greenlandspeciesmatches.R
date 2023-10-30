 us <- sort(c(
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
 
 setdiff(us,them)

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
# 