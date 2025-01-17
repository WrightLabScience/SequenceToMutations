# make YAML files for PGAP

# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')

suppressMessages(library(stringr))
suppressMessages(library(yaml))

pgap_path <- "assemblies/"

lf <- list.files(path = pgap_path, 
                 pattern = "\\.fasta$")


strain <- as.integer(gsub(pattern = "ancestor_(\\d{1,2})_contigs.fasta",
                          replacement = "\\1",
                          x = lf))

CONTROLLERS <- SUBMOL <- character(length = length(strain))

for (i in seq_along(strain)) {
   # write out both YAML files
   
   # if you want single digit strain numbers to be formatted as 2 digits.
   f <- strain[i]
   # f <- formatC(x = strain[i],
   #              width = 2L,
   #              flag = 0,
   #              format = "d")
   
   # some of this may need to be changed
   Y1 <- as.yaml(list("fasta" = data.frame("class" = "File",
                                           "location" = paste0("/srv/ancestor_", f, "_contigs.fasta")),
                      "submol" = data.frame("class" = "File",
                                            "location" = paste0("/srv/", "submol_", f, ".yaml")),
                      "supplemental_data" = "{ class: Directory, location: /srv/input-2022-04-14.build6021 }", # such as this...?
                      "report_usage" = "true",
                      "ignore_all_errors" = "true"))
   Y1 <- str_replace_all(string = Y1,
                         pattern = "'",
                         replacement = "")
   cat(Y1,
       file = paste0("YAMLfiles/Controller_", f, ".yaml"))
   
   CONTROLLERS[i] <- paste0("Controller_", f, ".yaml")
   SUBMOL[i] <- paste0("submol_", f, ".yaml")
   
   Y2 <- as.yaml(list("topology" = "linear",
                      "organism" = list("genus_species" = "Staphylococcus aureus"), # change the organism name, if needed
                      "authors" = list(list("author" = list("last_name" = "Cooley",
                                                            "first_name" = "Nicholas",
                                                            "middle_initial" = "P"))),
                      "contact_info" = list("first_name" = "Nicholas",
                                            "last_name" = "Cooley",
                                            "email" = "npc19@pitt.edu",
                                            "organization" = "University of Pittsburgh",
                                            "department" = "Department of Biomedical Informatics",
                                            "phone" = "412-624-5100",
                                            "street" = "100 Technology Dr.",
                                            "city" = "Pittsburgh",
                                            "postal_code" = "15219",
                                            "country" = "USA")))
   Y2 <- gsub(pattern = "(?<=: )([^'\n]+)",
              replacement = "'\\1'",
              x = Y2,
              perl = TRUE)
   
   cat(Y2,
       file = paste0("YAMLfiles/submol_", f, ".yaml"))
}




RESULT <- gsub(pattern = "Controller",
               replacement = "Annot",
               x = CONTROLLERS)
RESULT <- gsub(pattern = "yaml",
               replacement = "gbk",
               x = RESULT)
FNAS <- lf

PGAPMap <- cbind(seq_along(strain),
                 CONTROLLERS,
                 SUBMOL,
                 RESULT,
                 FNAS)

write.table(x = PGAPMap,
            file = "OSG_job_files/pgap_job_map.txt",
            quote = FALSE,
            append = FALSE,
            col.names = FALSE,
            row.names = FALSE)
