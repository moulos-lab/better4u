# Utility to export installed R packages to an r_package_list.tsv file

# Get list of packages with version
pack_vers <- installed.packages()[,"Version"]

# Specify the output files
md_output_file <- "r_package_list.md"
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

# Helper function to format table with dynamic column widths
format_table <- function(pkg_list) {
  max_length <- max(nchar(names(pkg_list)))
  fmt <- sprintf("%%-%ds | %%s\n", max_length)
  
  table <- sprintf(fmt, "Package", "Version")
  separator <- sprintf(fmt, strrep("-", max_length), strrep("-", 10))
  table <- paste0(table, separator, sep = "")
  
  for (pkg in names(pkg_list)) {
    table <- paste0(table, sprintf(fmt, pkg, pkg_list[pkg]), sep = "")
  }
  
  table
}

# Create Markdown content
create_markdown <- function(user_defined, other_packages, file_name) {
  # Create user-defined packages section
  user_defined_markdown <- "
## User-defined Packages

"
  user_defined_markdown <- paste0(user_defined_markdown, format_table(user_defined))
  
  # Create other packages section
  other_markdown <- "
## Other Packages

"
  other_markdown <- paste0(other_markdown, format_table(other_packages))
  
  # Write to file
  writeLines(c(user_defined_markdown, other_markdown), con = file_name)
}

# Create the Markdown file
create_markdown(user_defined_pack_vers, other_pack_vers, md_output_file)