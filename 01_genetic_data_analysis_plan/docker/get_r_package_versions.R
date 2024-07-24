# Utility to export installed R packages to an r_package_list.tsv file

# Get list of packages with version
pack_vers <- installed.packages()[,"Version"]

# Specify the output files
rmd_output_file <- "r_package_list.Rmd"
tsv_output_file <- "r_package_list.tsv"

# Convert the named character list to a data frame
df <- data.frame(Package = names(pack_vers), Version = pack_vers, stringsAsFactors = FALSE)

# Write the data frame to a TSV file
write.table(df, "r_package_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Define the user-specified packages
user_defined <- c(
  "susieR", "snpStats", "parallel", "utils", 
  "remotes", "R.utils", "statgenGWAS", 
  "rmarkdown", "SAIGE", "lassosum"
)

# Create a new list to store sorted packages
user_defined_pack_vers <- c()
other_pack_vers <- c()

# Add user-defined packages first, marking those not installed
for (pkg in user_defined) {
  if (pkg %in% names(pack_vers)) {
    user_defined_pack_vers[pkg] <- pack_vers[pkg]
  } else {
    user_defined_pack_vers[pkg] <- "Not installed"
  }
}

# Add the rest of the packages
other_packages <- setdiff(names(pack_vers), user_defined)
for (pkg in other_packages) {
  other_pack_vers[pkg] <- pack_vers[pkg]
}

# Write the list to the file in R Markdown table format
writeLines(
  c(
    "## User-defined Packages",
    "",
    "```{r}",
    "user_defined_pack_vers <- c(",
    paste(
      sprintf("  %s = \"%s\"", names(user_defined_pack_vers), user_defined_pack_vers),
      collapse = ",\n"
    ),
    ")",
    "```",
    "",
    "Package | Version",
    "--- | ---",
    paste(
      sprintf("%s | %s", names(user_defined_pack_vers), user_defined_pack_vers),
      collapse = "\n"
    ),
    "",
    "## Other Packages",
    "",
    "```{r}",
    "other_pack_vers <- c(",
    paste(
      sprintf("  %s = \"%s\"", names(other_pack_vers), other_pack_vers),
      collapse = ",\n"
    ),
    ")",
    "```",
    "",
    "Package | Version",
    "--- | ---",
    paste(
      sprintf("%s | %s", names(other_pack_vers), other_pack_vers),
      collapse = "\n"
    )
  ),
  rmd_output_file
)
