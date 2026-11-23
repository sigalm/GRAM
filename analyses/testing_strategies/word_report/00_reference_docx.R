######################################## WORD REFERENCE TEMPLATE ########################################
# Builds the reference document the report is knitted against. It does two things:
#
#   1. Page setup: US Letter, landscape, 0.8 in margins. Landscape because the results
#      tables run to ten columns and the main figure is three panels wide -- on a portrait
#      page both have to shrink to the point where the type stops being readable.
#      pandoc takes page size and margins from the reference document's DEFAULT (final
#      body) section, which is what body_set_default_section writes; a section added
#      mid-document is ignored, because pandoc never reaches it.
#
#   2. Strips the automatic outline numbering off the heading styles. officer's stock
#      template numbers headings, so a "## Table 1." heading arrives in Word as "1.2
#      Table 1." -- two numbering schemes fighting each other, since the document already
#      numbers its own sections and names its own tables. There is no officer API for
#      this, so the styles part is edited in the file itself.

library(officer)

out <- "analyses/testing_strategies/word_report/reference_landscape.docx"

read_docx() |>
  body_set_default_section(prop_section(
    page_size    = page_size(orient = "landscape", width = 11, height = 8.5),
    page_margins = page_mar(top = 0.8, bottom = 0.8, left = 0.8, right = 0.8,
                            header = 0.4, footer = 0.4, gutter = 0))) |>
  print(target = out)


## Drop <w:numPr> from every heading style ####
# Matched on w:name ("heading 1"), not on the style id: the stock template is French, so
# the ids are Titre1/Titre2/Titre3.
local({
  tmp <- file.path(tempdir(), "refdocx")
  unlink(tmp, recursive = TRUE)
  dir.create(tmp)
  files <- unzip(out, exdir = tmp)

  styles_path <- file.path(tmp, "word", "styles.xml")
  xml <- readLines(styles_path, warn = FALSE) |> paste(collapse = "\n")

  # One style element at a time, so numbering on list styles elsewhere is left alone.
  xml <- gsub('(<w:style [^>]*>(?:(?!</w:style>).)*?<w:name w:val="heading [0-9]"[^>]*/>(?:(?!</w:style>).)*?)<w:numPr>.*?</w:numPr>',
              "\\1", xml, perl = TRUE)

  writeLines(xml, styles_path)

  # zip from inside the extraction directory: the paths stored in the archive have to be
  # relative to the package root, or Word cannot find the parts.
  wd <- setwd(tmp)
  on.exit(setwd(wd), add = TRUE)
  utils::zip(normalizePath(file.path(wd, out)),
             files = sub(paste0("^", tmp, "/"), "", files), flags = "-qr9X")
})

cat("Wrote", out, "\n")
